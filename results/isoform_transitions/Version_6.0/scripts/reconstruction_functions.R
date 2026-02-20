#!/usr/bin/env Rscript
#
# Reconstruction Functions v2 - Union Exon Based
#
# Refactored approach: Since event coordinates exactly match atomic union exon
# boundaries, we reconstruct by finding and applying the specific union exons
# that represent each event.
#
# Key principle:
#   Event coordinates → Union exons → Add/remove from comparator structure

library(tidyverse)

# ==============================================================================
# Core Helper: Find Union Exons for Event
# ==============================================================================

#' Parse comma-separated exon ranges and find matching union exons
#'
#' @param ranges_str Comma-separated ranges like "72000-72100" or "100-200, 300-400"
#' @param union_exons All union exons for the gene
#' @param gene_id Gene identifier
#' @return Tibble of matching union exons
parse_and_find_union_exons <- function(ranges_str, union_exons, gene_id) {
  if (is.na(ranges_str) || ranges_str == "") {
    return(tibble())
  }

  # Parse comma-separated ranges
  ranges_list <- strsplit(ranges_str, ",")[[1]]
  all_matches <- list()

  for (range_str in ranges_list) {
    coords <- as.integer(strsplit(trimws(range_str), "-")[[1]])
    range_start <- coords[1]
    range_end <- coords[2]

    # Find union exons within this range
    matches <- union_exons %>%
      filter(
        gene_id == !!gene_id,
        start >= range_start,
        end <= range_end
      ) %>%
      arrange(start)

    if (nrow(matches) > 0) {
      all_matches[[length(all_matches) + 1]] <- matches
    }
  }

  if (length(all_matches) > 0) {
    return(bind_rows(all_matches) %>% distinct())
  } else {
    return(tibble())
  }
}

#' Find union exons that match an event's coordinates
#'
#' @param event Event record
#' @param union_exons All union exons for the gene
#' @return Tibble of matching union exons
find_event_union_exons <- function(event, union_exons) {
  # For Alt_TSS/Alt_TES: use missing_terminal_exons
  if (event$event_type %in% c("Alt_TSS", "Alt_TES")) {
    if (is.na(event$missing_terminal_exons) || event$missing_terminal_exons == "") {
      return(tibble())
    }

    # Parse comma-separated ranges
    ranges_str <- strsplit(event$missing_terminal_exons, ",")[[1]]
    all_matches <- list()

    for (range_str in ranges_str) {
      coords <- as.integer(strsplit(trimws(range_str), "-")[[1]])
      range_start <- coords[1]
      range_end <- coords[2]

      # Find union exons that overlap this range, then clip to range boundaries.
      # Strict containment (start >= range_start) misses union exons that straddle
      # the event boundary — clipping ensures we capture the correct portion.
      matches <- union_exons %>%
        filter(
          gene_id == event$gene_id,
          start <= range_end,
          end >= range_start
        ) %>%
        mutate(
          start = pmax(start, range_start),
          end = pmin(end, range_end)
        ) %>%
        arrange(start)

      if (nrow(matches) > 0) {
        all_matches[[length(all_matches) + 1]] <- matches
      }
    }

    if (length(all_matches) > 0) {
      return(bind_rows(all_matches) %>% distinct())
    } else {
      return(tibble())
    }
  }

  # For other events: use five_prime/three_prime
  event_start <- min(event$five_prime, event$three_prime)
  event_end <- max(event$five_prime, event$three_prime)

  # Try strict containment first (preserves existing behavior)
  matches <- union_exons %>%
    filter(
      gene_id == event$gene_id,
      start >= event_start,
      end <= event_end
    ) %>%
    arrange(start)

  # Fallback: if strict containment finds nothing, try overlap + clip.
  # This handles union exons that straddle the event boundary.
  if (nrow(matches) == 0) {
    matches <- union_exons %>%
      filter(
        gene_id == event$gene_id,
        start <= event_end,
        end >= event_start
      ) %>%
      mutate(
        start = pmax(start, event_start),
        end = pmin(end, event_end)
      ) %>%
      arrange(start)
  }

  return(matches)
}

# ==============================================================================
# Union Exon Addition/Removal
# ==============================================================================

#' Add union exons to comparator structure
#'
#' @param exons Current exon structure
#' @param union_exons_to_add Union exons to insert
#' Merge adjacent or overlapping exons into continuous segments
#'
#' @param exons Exon dataframe with exon_start and exon_end
#' @param strand Strand (not used currently, but kept for consistency)
#' @return Merged exon dataframe
merge_adjacent_exons <- function(exons, strand) {
  if (nrow(exons) == 0) return(exons)

  # Sort by position
  exons <- exons %>% arrange(exon_start)

  merged <- list()
  current <- exons[1, ]

  for (i in seq_len(nrow(exons))[-1]) {
    next_exon <- exons[i, ]

    # Check if adjacent or overlapping (end of current >= start-1 of next)
    if (current$exon_end >= next_exon$exon_start - 1) {
      # Merge: extend current to cover next
      current$exon_end <- max(current$exon_end, next_exon$exon_end)
    } else {
      # Gap exists: save current and start new
      merged[[length(merged) + 1]] <- current
      current <- next_exon
    }
  }

  # Add final exon
  merged[[length(merged) + 1]] <- current

  bind_rows(merged)
}

#' @param strand Strand
#' @return Modified exon structure
add_union_exons <- function(exons, union_exons_to_add, strand) {
  if (nrow(union_exons_to_add) == 0) {
    return(exons)
  }

  # Convert union exons to exon format
  new_exons <- union_exons_to_add %>%
    mutate(
      exon_start = start,
      exon_end = end,
      transcript_id = exons$transcript_id[1]
    ) %>%
    select(chr, exon_start, exon_end, strand, gene_id, transcript_id)

  # Combine and sort
  combined <- bind_rows(exons, new_exons) %>%
    arrange(exon_start, exon_end) %>%
    distinct()

  return(combined)
}

#' Remove union exons from comparator structure
#'
#' @param exons Current exon structure
#' @param union_exons_to_remove Union exons to remove
#' @return Modified exon structure
remove_union_exons <- function(exons, union_exons_to_remove) {
  if (nrow(union_exons_to_remove) == 0) {
    return(exons)
  }

  # Remove exons that exactly match the union exon coordinates
  for (i in seq_len(nrow(union_exons_to_remove))) {
    ue <- union_exons_to_remove[i, ]
    exons <- exons %>%
      filter(!(exon_start == ue$start & exon_end == ue$end))
  }

  return(exons)
}

#' Modify exon boundary using event coordinates directly
#'
#' For splice site events (A3SS, A5SS, Partial_IR), extend or trim an existing
#' exon boundary. Target boundary is derived algebraically from the event's
#' five_prime/three_prime coordinates — no union exon lookup needed.
#'
#' @param exons Current exon structure
#' @param event Event record (uses five_prime, three_prime, event_type, direction, strand)
#' @return Modified exon structure
modify_exon_boundary <- function(exons, event) {
  strand <- event$strand
  event_type <- event$event_type
  direction <- event$direction

  fp <- event$five_prime
  tp <- event$three_prime

  # Determine which boundary to modify based on event type
  if (event_type %in% c("A5SS", "Partial_IR_5")) {
    # Donor differs (5' splice site)
    if (strand == "+") {
      # Plus: donor = exon end
      if (direction == "LOSS") {
        # Extend exon end to dom_end = max(fp, tp)
        boundary_coord <- min(fp, tp) - 1
        exon_idx <- which(abs(exons$exon_end - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_end[exon_idx] <- max(fp, tp)
        }
      } else {
        # Trim exon end to dom_end = min(fp, tp) - 1
        boundary_coord <- max(fp, tp)
        exon_idx <- which(abs(exons$exon_end - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_end[exon_idx] <- min(fp, tp) - 1
        }
      }
    } else {
      # Minus: donor = exon start
      if (direction == "LOSS") {
        # Extend exon start to dom_start = min(fp, tp) + 1
        boundary_coord <- max(fp, tp) + 1
        exon_idx <- which(abs(exons$exon_start - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_start[exon_idx] <- min(fp, tp) + 1
        }
      } else {
        # Trim exon start to dom_start = max(fp, tp) + 1
        boundary_coord <- min(fp, tp)
        exon_idx <- which(abs(exons$exon_start - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_start[exon_idx] <- max(fp, tp) + 1
        }
      }
    }
  } else if (event_type %in% c("A3SS", "Partial_IR_3")) {
    # Acceptor differs (3' splice site)
    if (strand == "+") {
      # Plus: acceptor = exon start
      if (direction == "LOSS") {
        # Extend exon start to dom_start = min(fp, tp)
        boundary_coord <- max(fp, tp) + 1
        exon_idx <- which(abs(exons$exon_start - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_start[exon_idx] <- min(fp, tp)
        }
      } else {
        # Trim exon start to dom_start = max(fp, tp) + 1
        boundary_coord <- min(fp, tp)
        exon_idx <- which(abs(exons$exon_start - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_start[exon_idx] <- max(fp, tp) + 1
        }
      }
    } else {
      # Minus: acceptor = exon end
      if (direction == "LOSS") {
        # Extend exon end to dom_end = max(fp, tp)
        boundary_coord <- min(fp, tp) - 1
        exon_idx <- which(abs(exons$exon_end - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_end[exon_idx] <- max(fp, tp)
        }
      } else {
        # Trim exon end to dom_end = min(fp, tp) - 1
        boundary_coord <- max(fp, tp)
        exon_idx <- which(abs(exons$exon_end - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_end[exon_idx] <- min(fp, tp) - 1
        }
      }
    }
  }

  return(exons)
}

#' Modify terminal exon for Alt_TSS/Alt_TES events
#'
#' LOSS: adds missing exon ranges from missing_terminal_exons
#' GAIN: removes gained sequence using five_prime (dominant's terminal boundary)
#'       to truncate/remove exons beyond the dominant's extent
#' @param exons Current exon structure
#' @param event Event record (uses five_prime for GAIN, missing_terminal_exons for LOSS)
#' @return Modified exon structure
modify_terminal_exon <- function(exons, event) {
  five_prime <- event$five_prime

  if (event$direction == "LOSS") {
    # ADD: insert each missing range as a new exon
    if (!is.na(event$missing_terminal_exons) && event$missing_terminal_exons != "") {
      ranges_str <- strsplit(event$missing_terminal_exons, ",")[[1]]
      for (rs in ranges_str) {
        coords <- as.integer(strsplit(trimws(rs), "-")[[1]])
        new_exon <- tibble(
          chr = exons$chr[1], exon_start = coords[1], exon_end = coords[2],
          strand = exons$strand[1], gene_id = exons$gene_id[1],
          transcript_id = exons$transcript_id[1]
        )
        exons <- bind_rows(exons, new_exon)
      }
      exons <- exons %>% arrange(exon_start) %>%
        distinct(exon_start, exon_end, .keep_all = TRUE)
    }
  } else {
    # GAIN: remove gained terminal sequence using five_prime (dominant's boundary)
    # Determine which side to cut based on event type and strand
    remove_high <- (event$event_type == "Alt_TES" && event$strand == "+") ||
                   (event$event_type == "Alt_TSS" && event$strand == "-")

    if (remove_high) {
      # Remove exons entirely beyond five_prime (high side)
      exons <- exons %>% filter(exon_start <= five_prime)
      # Truncate overlapping exon's end to dominant's boundary
      exons <- exons %>% mutate(exon_end = pmin(exon_end, five_prime))
    } else {
      # Remove exons entirely beyond five_prime (low side)
      exons <- exons %>% filter(exon_end >= five_prime)
      # Truncate overlapping exon's start to dominant's boundary
      exons <- exons %>% mutate(exon_start = pmax(exon_start, five_prime))
    }

    # Remove any degenerate exons (safety)
    exons <- exons %>% filter(exon_start <= exon_end)
  }

  return(exons)
}

# ==============================================================================
# Event Application Functions
# ==============================================================================

#' Apply a single event to reconstruct the dominant isoform
#'
#' RECONSTRUCTION PRINCIPLE: Reconstruction is performed using ONLY isoform pair
#' assignments and event descriptions (event type + associated coordinates/junctions).
#' Union exon lookups are used ONLY for IR events, where the intronic vs exonic split
#' structure genuinely requires atomic union exon granularity. All other event types
#' derive their target coordinates directly from the event's five_prime/three_prime
#' and missing_terminal_exons fields. This eliminates an entire class of bugs caused
#' by union exon boundary straddling.
#'
#' @param exons Current exon structure
#' @param event Event record
#' @param union_exons All union exons for the gene (used only by IR events)
#' @return Modified exon structure
apply_event_union_based <- function(exons, event, union_exons) {
  event_type <- event$event_type
  direction  <- event$direction

  # -------------------------------------------------------------------
  # SE / Missing_Internal: coordinate-based (no union exon lookup)
  # -------------------------------------------------------------------
  if (event_type %in% c("SE", "Missing_Internal")) {
    event_start <- min(event$five_prime, event$three_prime)
    event_end   <- max(event$five_prime, event$three_prime)

    if (direction == "GAIN") {
      # Comparator has extra exon → remove it by coordinate range
      exons <- exons %>% filter(!(exon_start >= event_start & exon_end <= event_end))
    } else {
      # Dominant has exon → add it directly from event coordinates
      new_exon <- tibble(
        chr = exons$chr[1], exon_start = event_start, exon_end = event_end,
        strand = exons$strand[1], gene_id = exons$gene_id[1],
        transcript_id = exons$transcript_id[1]
      )
      exons <- bind_rows(exons, new_exon) %>% arrange(exon_start)
    }
    return(exons)
  }

  # -------------------------------------------------------------------
  # A5SS / A3SS / Partial_IR: direct boundary modification (no UE lookup)
  # -------------------------------------------------------------------
  if (event_type %in% c("A5SS", "A3SS", "Partial_IR_5", "Partial_IR_3")) {
    exons <- modify_exon_boundary(exons, event)
    return(exons)
  }

  # -------------------------------------------------------------------
  # Alt_TSS / Alt_TES: coordinate-based terminal modification (no UE lookup)
  # -------------------------------------------------------------------
  if (event_type %in% c("Alt_TSS", "Alt_TES")) {
    exons <- modify_terminal_exon(exons, event)

    # Remove comp orphan terminal exons for both GAIN and LOSS.
    if (!is.na(event$orphan_terminal_exons) && event$orphan_terminal_exons != "") {
      ranges_list <- strsplit(event$orphan_terminal_exons, ",")[[1]]
      for (range_str in ranges_list) {
        coords <- as.integer(strsplit(trimws(range_str), "-")[[1]])
        exons <- exons %>% filter(!(exon_start >= coords[1] & exon_end <= coords[2]))
      }
    }
    return(exons)
  }

  # -------------------------------------------------------------------
  # IR events: direct coordinate reconstruction (no union exon lookup)
  #
  # Event coordinates contain everything needed:
  #   five_prime/three_prime = IR region extent (= retaining exon boundaries)
  #   ir_split_exons = the other isoform's exon coordinates in this region
  #     GAIN: ir_split = dominant's split exon coordinates
  #     LOSS: ir_split = comparator's split exon coordinates
  # -------------------------------------------------------------------
  if (event_type %in% c("IR", "IR_diff_5", "IR_diff_3", "IR_diff_5_3")) {
    ir_start <- min(event$five_prime, event$three_prime)
    ir_end <- max(event$five_prime, event$three_prime)

    # Get metadata for new exon creation
    meta_chr <- if (nrow(exons) > 0) exons$chr[1] else event$chr
    meta_strand <- if (nrow(exons) > 0) exons$strand[1] else event$strand
    meta_gene <- if (nrow(exons) > 0) exons$gene_id[1] else event$gene_id
    meta_tx <- if (nrow(exons) > 0) exons$transcript_id[1] else NA_character_

    if (direction == "GAIN") {
      # Comparator retained intron → split into dominant's exons
      # 1. Remove the comp exon spanning the IR region
      for (i in seq_len(nrow(exons))) {
        if (exons$exon_start[i] <= ir_start && exons$exon_end[i] >= ir_end) {
          exons <- exons[-i, ]
          break
        }
      }

      # 2. Add dominant's split exons directly from ir_split_exons coordinates
      if (!is.na(event$ir_split_exons) && event$ir_split_exons != "") {
        split_ranges <- strsplit(event$ir_split_exons, ",")[[1]]
        for (range_str in split_ranges) {
          coords <- as.integer(strsplit(trimws(range_str), "-")[[1]])
          new_exon <- tibble(
            chr = meta_chr, exon_start = coords[1], exon_end = coords[2],
            strand = meta_strand, gene_id = meta_gene, transcript_id = meta_tx
          )
          exons <- bind_rows(exons, new_exon)
        }
        exons <- exons %>% arrange(exon_start)
      } else {
        warning("IR GAIN event missing ir_split_exons field - cannot split properly")
      }

    } else {
      # IR LOSS: dominant retained intron → merge comp's split exons
      # 1. Remove comp exons overlapping the IR region
      exons <- exons %>% filter(!(exon_start <= ir_end & exon_end >= ir_start))

      # 2. Add one merged exon = the dominant's retained region
      merged_exon <- tibble(
        chr = meta_chr, exon_start = ir_start, exon_end = ir_end,
        strand = meta_strand, gene_id = meta_gene, transcript_id = meta_tx
      )
      exons <- bind_rows(exons, merged_exon) %>% arrange(exon_start)
    }
    return(exons)
  }

  return(exons)
}

# ==============================================================================
# Main Reconstruction Function
# ==============================================================================

#' Reconstruct dominant isoform from comparator + events
#'
#' PRINCIPLE: Reconstruction is performed using ONLY isoform pair assignments
#' and event descriptions (event type + associated coordinates). Union exon
#' lookups are reserved for IR events only.
#'
#' @param comparator_exons Exons from comparator isoform
#' @param events Events for this isoform pair
#' @param union_exons Union exon structure for the gene (used only by IR events)
#' @return Reconstructed dominant exons
reconstruct_dominant_v2 <- function(comparator_exons, events, union_exons) {
  # Start with comparator structure
  reconstructed <- comparator_exons

  # Two-phase reconstruction:
  # Phase 1: Internal events (IR, SE, splice sites) - merge after to consolidate structure
  # Phase 2: Terminal events (Alt_TSS, Alt_TES) - NO merge to keep new exons separate

  # Separate events into internal vs terminal
  internal_events_raw <- events %>%
    filter(event_type %in% c("IR", "IR_diff_5", "IR_diff_3", "IR_diff_5_3",
                              "SE", "Missing_Internal", "A5SS", "A3SS", "Partial_IR_5", "Partial_IR_3"))

  # Filter out Partial_IR events that overlap with IR events
  # (IR events handle the full retention, Partial_IR becomes redundant)
  ir_events <- internal_events_raw %>% filter(event_type %in% c("IR", "IR_diff_5", "IR_diff_3", "IR_diff_5_3"))
  partial_ir_events <- internal_events_raw %>% filter(event_type %in% c("Partial_IR_5", "Partial_IR_3"))
  other_events <- internal_events_raw %>% filter(!event_type %in% c("IR", "Partial_IR_5", "Partial_IR_3"))

  # Check each Partial_IR against IR events for overlap
  if (nrow(ir_events) > 0 && nrow(partial_ir_events) > 0) {
    keep_partial <- logical(nrow(partial_ir_events))
    for (i in seq_len(nrow(partial_ir_events))) {
      pe <- partial_ir_events[i, ]
      pe_start <- min(pe$five_prime, pe$three_prime)
      pe_end <- max(pe$five_prime, pe$three_prime)

      # Check if any IR event overlaps this Partial_IR
      overlaps <- FALSE
      for (j in seq_len(nrow(ir_events))) {
        ie <- ir_events[j, ]
        ie_start <- min(ie$five_prime, ie$three_prime)
        ie_end <- max(ie$five_prime, ie$three_prime)

        if (ie_start <= pe_end && ie_end >= pe_start) {
          overlaps <- TRUE
          cat(sprintf("  [DEBUG] Filtering %s [%d-%d] (overlaps with IR [%d-%d])\n",
                     pe$event_type, pe_start, pe_end, ie_start, ie_end))
          break
        }
      }
      keep_partial[i] <- !overlaps
    }
    partial_ir_events <- partial_ir_events[keep_partial, ]
  }

  # Combine filtered events
  internal_events <- bind_rows(ir_events, partial_ir_events, other_events) %>%
    arrange(
      case_when(
        event_type %in% c("IR", "IR_diff_5", "IR_diff_3", "IR_diff_5_3") ~ 1,
        event_type == "SE" ~ 2,
        TRUE ~ 3
      )
    )

  terminal_events <- events %>%
    filter(event_type %in% c("Alt_TSS", "Alt_TES")) %>%
    arrange(direction != "LOSS", event_type)  # LOSS before GAIN, then by event type

  # PHASE 1: Apply internal events
  for (i in seq_len(nrow(internal_events))) {
    event <- internal_events[i, ]
    reconstructed <- tryCatch({
      apply_event_union_based(reconstructed, event, union_exons)
    }, error = function(e) {
      warning(sprintf("Error applying %s event: %s", event$event_type, e$message))
      reconstructed
    })
  }

  # Pre-merge orphan removal: terminal events may mark comp exons as orphans
  # that need to be removed. If these orphans are adjacent to internally-modified
  # exons, the merge step would erroneously consolidate them. Remove first.
  if (nrow(terminal_events) > 0) {
    for (i in seq_len(nrow(terminal_events))) {
      te <- terminal_events[i, ]
      if (!is.na(te$orphan_terminal_exons) && te$orphan_terminal_exons != "") {
        ranges_list <- strsplit(te$orphan_terminal_exons, ",")[[1]]
        for (range_str in ranges_list) {
          coords <- as.integer(strsplit(trimws(range_str), "-")[[1]])
          reconstructed <- reconstructed %>%
            filter(!(exon_start >= coords[1] & exon_end <= coords[2]))
        }
      }
    }
  }

  # Merge after internal events to consolidate structure
  reconstructed <- reconstructed %>%
    arrange(exon_start, exon_end) %>%
    distinct(exon_start, exon_end, .keep_all = TRUE)

  if (nrow(reconstructed) > 0) {
    strand <- reconstructed$strand[1]
    reconstructed <- merge_adjacent_exons(reconstructed, strand)
  }

  # PHASE 2: Apply terminal events WITHOUT merging
  for (i in seq_len(nrow(terminal_events))) {
    event <- terminal_events[i, ]
    reconstructed <- tryCatch({
      apply_event_union_based(reconstructed, event, union_exons)
    }, error = function(e) {
      warning(sprintf("Error applying %s event: %s", event$event_type, e$message))
      reconstructed
    })
  }

  # Final cleanup: sort, deduplicate, and merge adjacent exons.
  # Terminal event reconstruction adds atomic union exons (which are split at every
  # boundary). Adjacent atomic UEs representing the same logical exon must be merged.
  reconstructed <- reconstructed %>%
    arrange(exon_start, exon_end) %>%
    distinct(exon_start, exon_end, .keep_all = TRUE)

  if (nrow(reconstructed) > 0) {
    strand <- reconstructed$strand[1]
    reconstructed <- merge_adjacent_exons(reconstructed, strand)
  }

  return(reconstructed)
}
