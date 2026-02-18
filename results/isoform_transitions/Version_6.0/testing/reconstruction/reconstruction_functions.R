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

      # Find union exons that exactly match this range
      matches <- union_exons %>%
        filter(
          gene_id == event$gene_id,
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

  # For other events: use five_prime/three_prime
  event_start <- min(event$five_prime, event$three_prime)
  event_end <- max(event$five_prime, event$three_prime)

  # Find union exons within the event range
  matches <- union_exons %>%
    filter(
      gene_id == event$gene_id,
      start >= event_start,
      end <= event_end
    ) %>%
    arrange(start)

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

#' Modify exon boundary using union exon
#'
#' For splice site events (A3SS, A5SS, Partial_IR), we extend or trim
#' an existing exon boundary rather than adding/removing whole exons
#'
#' @param exons Current exon structure
#' @param event_ues Union exons representing the event
#' @param event Event record
#' @return Modified exon structure
modify_exon_boundary <- function(exons, event_ues, event) {
  if (nrow(event_ues) == 0) {
    return(exons)
  }

  strand <- event$strand
  event_type <- event$event_type
  direction <- event$direction

  # Determine which boundary to modify based on event type
  if (event_type %in% c("A5SS", "Partial_IR_5")) {
    # Donor differs (5' splice site)
    if (strand == "+") {
      # Plus: donor = exon end
      if (direction == "LOSS") {
        # Extend exon end
        boundary_coord <- min(event$five_prime, event$three_prime) - 1
        exon_idx <- which(abs(exons$exon_end - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_end[exon_idx] <- max(event_ues$end)
        }
      } else {
        # Trim exon end
        boundary_coord <- max(event$five_prime, event$three_prime)
        exon_idx <- which(abs(exons$exon_end - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_end[exon_idx] <- min(event_ues$start) - 1
        }
      }
    } else {
      # Minus: donor = exon start
      if (direction == "LOSS") {
        # Extend exon start
        boundary_coord <- max(event$five_prime, event$three_prime) + 1
        exon_idx <- which(abs(exons$exon_start - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_start[exon_idx] <- min(event_ues$start)
        }
      } else {
        # Trim exon start
        boundary_coord <- min(event$five_prime, event$three_prime)
        exon_idx <- which(abs(exons$exon_start - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_start[exon_idx] <- max(event_ues$end) + 1
        }
      }
    }
  } else if (event_type %in% c("A3SS", "Partial_IR_3")) {
    # Acceptor differs (3' splice site)
    if (strand == "+") {
      # Plus: acceptor = exon start
      if (direction == "LOSS") {
        # Extend exon start
        boundary_coord <- max(event$five_prime, event$three_prime) + 1
        exon_idx <- which(abs(exons$exon_start - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_start[exon_idx] <- min(event_ues$start)
        }
      } else {
        # Trim exon start
        boundary_coord <- min(event$five_prime, event$three_prime)
        exon_idx <- which(abs(exons$exon_start - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_start[exon_idx] <- max(event_ues$end) + 1
        }
      }
    } else {
      # Minus: acceptor = exon end
      if (direction == "LOSS") {
        # Extend exon end
        boundary_coord <- min(event$five_prime, event$three_prime) - 1
        exon_idx <- which(abs(exons$exon_end - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_end[exon_idx] <- max(event_ues$end)
        }
      } else {
        # Trim exon end
        boundary_coord <- max(event$five_prime, event$three_prime)
        exon_idx <- which(abs(exons$exon_end - boundary_coord) <= 1)[1]
        if (!is.na(exon_idx)) {
          exons$exon_end[exon_idx] <- min(event_ues$start) - 1
        }
      }
    }
  }

  return(exons)
}

#' Modify terminal exon using union exon coordinates
#'
#' For Alt_TSS/Alt_TES with multiple regions: process each separately
#' - Adjacent regions: merge with existing exon
#' - Gapped regions: add as separate exons
#'
#' @param exons Current exon structure
#' @param event_ues Union exons representing the terminal difference
#' @param event Event record
#' @param union_exons All union exons (for re-parsing regions)
#' @return Modified exon structure
modify_terminal_exon <- function(exons, event_ues, event, union_exons = NULL) {
  if (nrow(exons) == 0) {
    return(exons)
  }

  # For LOSS events with multiple terminal regions, handle each separately
  if (!is.null(union_exons) && event$direction == "LOSS" &&
      !is.na(event$missing_terminal_exons) && event$missing_terminal_exons != "") {

    # Parse individual regions
    ranges_str <- strsplit(event$missing_terminal_exons, ",")[[1]]

    if (length(ranges_str) > 1) {
      # Parse and sort regions based on strand and event type
      regions <- lapply(ranges_str, function(rs) {
        coords <- as.integer(strsplit(trimws(rs), "-")[[1]])
        list(start = coords[1], end = coords[2], str = trimws(rs))
      })

      # Sort regions to process adjacent regions first (closest to comparator terminal exon)
      # For LOSS: we're adding regions that extend beyond comparator's terminal exon
      if (event$event_type == "Alt_TSS") {
        if (event$strand == "+") {
          # Plus TSS (low coord): dominant extends to LOWER coords than comparator
          # Process from high to low (closest to comparator first)
          regions <- regions[order(sapply(regions, function(r) r$start), decreasing = TRUE)]
        } else {
          # Minus TSS (high coord): dominant extends to HIGHER coords than comparator
          # Process from low to high (closest to comparator first)
          regions <- regions[order(sapply(regions, function(r) r$start))]
        }
      } else {  # Alt_TES
        if (event$strand == "+") {
          # Plus TES (high coord): dominant extends to HIGHER coords than comparator
          # Process from low to high (closest to comparator first)
          regions <- regions[order(sapply(regions, function(r) r$start))]
        } else {
          # Minus TES (low coord): dominant extends to LOWER coords than comparator
          # Process from high to low (closest to comparator first)
          regions <- regions[order(sapply(regions, function(r) r$start), decreasing = TRUE)]
        }
      }

      # Multiple regions - process each individually in sorted order
      for (region in regions) {
        region_start <- region$start
        region_end <- region$end

        # Find union exons for this specific region
        region_ues <- union_exons %>%
          filter(
            gene_id == event$gene_id,
            start >= region_start,
            end <= region_end
          )

        if (nrow(region_ues) == 0) next

        # Check if this region is adjacent to the terminal exon
        if (event$event_type == "Alt_TSS") {
          if (event$strand == "+") {
            # Plus: TSS = first exon
            terminal_exon <- exons[1, ]
            is_adjacent <- (region_end + 1 >= terminal_exon$exon_start - 1) &&
                          (region_start <= terminal_exon$exon_start + 1)

            if (is_adjacent) {
              # Merge: extend terminal exon
              exons$exon_start[1] <- min(region_start, terminal_exon$exon_start)
              exons$exon_end[1] <- max(region_end, terminal_exon$exon_end)
            } else {
              # Gap: add as separate exon
              new_exon <- tibble(
                chr = exons$chr[1],
                exon_start = region_start,
                exon_end = region_end,
                strand = exons$strand[1],
                gene_id = exons$gene_id[1],
                transcript_id = exons$transcript_id[1]
              )
              exons <- bind_rows(new_exon, exons)
            }
          } else {
            # Minus: TSS = last exon
            terminal_exon <- exons[nrow(exons), ]
            is_adjacent <- (region_start <= terminal_exon$exon_end + 1) &&
                          (region_end + 1 >= terminal_exon$exon_end - 1)

            if (is_adjacent) {
              # Merge: extend terminal exon
              exons$exon_start[nrow(exons)] <- min(region_start, terminal_exon$exon_start)
              exons$exon_end[nrow(exons)] <- max(region_end, terminal_exon$exon_end)
            } else {
              # Gap: add as separate exon
              new_exon <- tibble(
                chr = exons$chr[1],
                exon_start = region_start,
                exon_end = region_end,
                strand = exons$strand[1],
                gene_id = exons$gene_id[1],
                transcript_id = exons$transcript_id[1]
              )
              exons <- bind_rows(exons, new_exon)
            }
          }
        } else if (event$event_type == "Alt_TES") {
          if (event$strand == "+") {
            # Plus: TES = last exon
            terminal_exon <- exons[nrow(exons), ]
            is_adjacent <- (region_start <= terminal_exon$exon_end + 1) &&
                          (region_end + 1 >= terminal_exon$exon_end - 1)

            if (is_adjacent) {
              # Merge: extend terminal exon
              exons$exon_start[nrow(exons)] <- min(region_start, terminal_exon$exon_start)
              exons$exon_end[nrow(exons)] <- max(region_end, terminal_exon$exon_end)
            } else {
              # Gap: add as separate exon
              new_exon <- tibble(
                chr = exons$chr[1],
                exon_start = region_start,
                exon_end = region_end,
                strand = exons$strand[1],
                gene_id = exons$gene_id[1],
                transcript_id = exons$transcript_id[1]
              )
              exons <- bind_rows(exons, new_exon)
            }
          } else {
            # Minus: TES = first exon
            terminal_exon <- exons[1, ]
            is_adjacent <- (region_end + 1 >= terminal_exon$exon_start - 1) &&
                          (region_start <= terminal_exon$exon_start + 1)

            if (is_adjacent) {
              # Merge: extend terminal exon
              exons$exon_start[1] <- min(region_start, terminal_exon$exon_start)
              exons$exon_end[1] <- max(region_end, terminal_exon$exon_end)
            } else {
              # Gap: add as separate exon
              new_exon <- tibble(
                chr = exons$chr[1],
                exon_start = region_start,
                exon_end = region_end,
                strand = exons$strand[1],
                gene_id = exons$gene_id[1],
                transcript_id = exons$transcript_id[1]
              )
              exons <- bind_rows(new_exon, exons)
            }
          }
        }
      }

      # Sort exons after processing all regions
      exons <- exons %>% arrange(exon_start) %>% distinct(exon_start, exon_end, .keep_all = TRUE)
      return(exons)
    }
  }

  # For GAIN events with multiple terminal regions, remove/trim all affected exons
  if (!is.null(union_exons) && event$direction == "GAIN" &&
      !is.na(event$missing_terminal_exons) && event$missing_terminal_exons != "") {

    # Parse individual regions
    ranges_str <- strsplit(event$missing_terminal_exons, ",")[[1]]

    if (length(ranges_str) > 1) {
      # Parse all regions that need to be removed
      regions <- lapply(ranges_str, function(rs) {
        coords <- as.integer(strsplit(trimws(rs), "-")[[1]])
        list(start = coords[1], end = coords[2])
      })

      # Process each exon and determine if it should be removed or trimmed
      exons_to_keep <- list()

      for (i in seq_len(nrow(exons))) {
        exon <- exons[i, ]
        keep_exon <- TRUE
        modified_exon <- exon

        # Check against all regions
        for (region in regions) {
          # Check if exon is completely within region
          if (exon$exon_start >= region$start && exon$exon_end <= region$end) {
            # Completely covered: remove this exon
            keep_exon <- FALSE
            break
          }

          # Check if exon partially overlaps region
          if (exon$exon_start < region$end && exon$exon_end > region$start) {
            # Partial overlap: trim the exon
            if (event$event_type == "Alt_TSS") {
              if (event$strand == "+") {
                # Plus TSS: trim start if it overlaps
                if (exon$exon_start < region$end && exon$exon_start >= region$start) {
                  modified_exon$exon_start <- region$end + 1
                }
              } else {
                # Minus TSS: trim end if it overlaps
                if (exon$exon_end > region$start && exon$exon_end <= region$end) {
                  modified_exon$exon_end <- region$start - 1
                }
              }
            } else {  # Alt_TES
              if (event$strand == "+") {
                # Plus TES: trim end if it overlaps
                if (exon$exon_end > region$start && exon$exon_end <= region$end) {
                  modified_exon$exon_end <- region$start - 1
                }
              } else {
                # Minus TES: trim start if it overlaps
                if (exon$exon_start < region$end && exon$exon_start >= region$start) {
                  modified_exon$exon_start <- region$end + 1
                }
              }
            }
          }
        }

        # Validate and keep exon if it's still valid
        if (keep_exon && modified_exon$exon_start <= modified_exon$exon_end) {
          exons_to_keep[[length(exons_to_keep) + 1]] <- modified_exon
        }
      }

      if (length(exons_to_keep) > 0) {
        exons <- bind_rows(exons_to_keep) %>%
          arrange(exon_start) %>%
          distinct(exon_start, exon_end, .keep_all = TRUE)
      } else {
        exons <- exons[0, ]  # Empty tibble with same structure
      }

      return(exons)
    }
  }

  # Fall through to original logic for single-region or GAIN events
  if (nrow(event_ues) == 0) {
    return(exons)
  }

  strand <- event$strand
  event_type <- event$event_type
  direction <- event$direction

  # Helper: check if union exons overlap or are adjacent to terminal exon
  # Allow 1bp gap for adjacent exons (they should be extended, not replaced)
  check_overlap_or_adjacent <- function(terminal_exon, union_exons) {
    any(
      (union_exons$start <= terminal_exon$exon_end + 1) &
      (union_exons$end + 1 >= terminal_exon$exon_start)
    )
  }

  if (event_type == "Alt_TES") {
    # TES (Transcription End Site)
    # Plus strand: TES = last exon end
    # Minus strand: TES = first exon start

    if (strand == "+") {
      # Terminal exon is last exon
      terminal_exon <- exons[nrow(exons), ]

      # Check if union exon exactly matches the terminal exon
      exact_match <- (nrow(event_ues) == 1 &&
                     event_ues$start[1] == terminal_exon$exon_start &&
                     event_ues$end[1] == terminal_exon$exon_end)

      if (exact_match && direction == "GAIN") {
        # Union exon IS the terminal exon → remove it entirely
        exons <- exons[-nrow(exons), ]
      } else if (check_overlap_or_adjacent(terminal_exon, event_ues)) {
        # Overlapping: extend or trim
        if (direction == "LOSS") {
          exons$exon_end[nrow(exons)] <- max(event_ues$end)
        } else {
          exons$exon_end[nrow(exons)] <- min(event_ues$start) - 1
        }
      } else {
        # Non-overlapping: remove and replace
        exons <- exons[-nrow(exons), ]  # Remove last exon
        if (direction == "LOSS") {
          # Add dominant's terminal exon(s)
          exons <- add_union_exons(exons, event_ues, strand)
        }
        # For GAIN, just removing is sufficient (dominant doesn't have this region)
      }

    } else {
      # Minus strand: terminal exon is first exon
      terminal_exon <- exons[1, ]

      # Check if union exon exactly matches the terminal exon
      exact_match <- (nrow(event_ues) == 1 &&
                     event_ues$start[1] == terminal_exon$exon_start &&
                     event_ues$end[1] == terminal_exon$exon_end)

      if (exact_match && direction == "GAIN") {
        # Union exon IS the terminal exon → remove it entirely
        exons <- exons[-1, ]
      } else if (check_overlap_or_adjacent(terminal_exon, event_ues)) {
        # Overlapping: extend or trim
        if (direction == "LOSS") {
          exons$exon_start[1] <- min(event_ues$start)
        } else {
          exons$exon_start[1] <- max(event_ues$end) + 1
        }
      } else {
        # Non-overlapping: remove and replace
        exons <- exons[-1, ]  # Remove first exon
        if (direction == "LOSS") {
          # Add dominant's terminal exon(s)
          exons <- add_union_exons(exons, event_ues, strand)
        }
        # For GAIN, just removing is sufficient (dominant doesn't have this region)
      }
    }

  } else if (event_type == "Alt_TSS") {
    # TSS (Transcription Start Site)
    # Plus strand: TSS = first exon start
    # Minus strand: TSS = last exon end

    if (strand == "+") {
      # Terminal exon is first exon
      terminal_exon <- exons[1, ]

      # Check if union exon exactly matches the terminal exon
      exact_match <- (nrow(event_ues) == 1 &&
                     event_ues$start[1] == terminal_exon$exon_start &&
                     event_ues$end[1] == terminal_exon$exon_end)

      if (exact_match && direction == "GAIN") {
        # Union exon IS the terminal exon → remove it entirely
        exons <- exons[-1, ]
      } else if (check_overlap_or_adjacent(terminal_exon, event_ues)) {
        # Overlapping: extend or trim
        if (direction == "LOSS") {
          exons$exon_start[1] <- min(event_ues$start)
        } else {
          exons$exon_start[1] <- max(event_ues$end) + 1
        }
      } else {
        # Non-overlapping: remove and replace
        exons <- exons[-1, ]  # Remove first exon
        if (direction == "LOSS") {
          # Add dominant's terminal exon(s)
          exons <- add_union_exons(exons, event_ues, strand)
        }
        # For GAIN, just removing is sufficient (dominant doesn't have this region)
      }

    } else {
      # Minus strand: terminal exon is last exon
      terminal_exon <- exons[nrow(exons), ]

      # Check if union exon exactly matches the terminal exon (for GAIN, should remove entirely)
      exact_match <- (nrow(event_ues) == 1 &&
                     event_ues$start[1] == terminal_exon$exon_start &&
                     event_ues$end[1] == terminal_exon$exon_end)

      if (exact_match && direction == "GAIN") {
        # Union exon IS the terminal exon → remove it entirely
        exons <- exons[-nrow(exons), ]
      } else if (check_overlap_or_adjacent(terminal_exon, event_ues)) {
        # Overlapping: extend or trim
        if (direction == "LOSS") {
          exons$exon_end[nrow(exons)] <- max(event_ues$end)
        } else {
          exons$exon_end[nrow(exons)] <- min(event_ues$start) - 1
        }
      } else {
        # Non-overlapping: remove and replace
        exons <- exons[-nrow(exons), ]  # Remove last exon
        if (direction == "LOSS") {
          # Add dominant's terminal exon(s)
          exons <- add_union_exons(exons, event_ues, strand)
        }
        # For GAIN, just removing is sufficient (dominant doesn't have this region)
      }
    }
  }

  return(exons)
}

# ==============================================================================
# Event Application Functions
# ==============================================================================

#' Apply any event using union exons
#'
#' @param exons Current exon structure
#' @param event Event record
#' @param union_exons All union exons for the gene
#' @return Modified exon structure
apply_event_union_based <- function(exons, event, union_exons) {
  # Find union exons for this event
  event_ues <- find_event_union_exons(event, union_exons)

  if (nrow(event_ues) == 0) {
    warning(sprintf("No union exons found for %s event", event$event_type))
    return(exons)
  }

  event_type <- event$event_type
  direction <- event$direction

  # Apply based on event type and direction
  if (event_type %in% c("Alt_TSS", "Alt_TES")) {
    # Terminal events: extend or trim the terminal exon
    exons <- modify_terminal_exon(exons, event_ues, event, union_exons)

    # Add missing internal exons if present
    if (!is.na(event$missing_internal_exons) && event$missing_internal_exons != "") {
      missing_internal_ues <- parse_and_find_union_exons(event$missing_internal_exons, union_exons, event$gene_id)
      if (nrow(missing_internal_ues) > 0) {
        exons <- add_union_exons(exons, missing_internal_ues, event$strand)
      }
    }

  } else if (event_type == "SE") {
    # Skipped exon: add or remove entire exon(s)
    if (direction == "LOSS") {
      # Dominant has exon, comparator doesn't → add it
      exons <- add_union_exons(exons, event_ues, event$strand)
    } else {
      # Comparator has exon, dominant doesn't → remove it
      exons <- remove_union_exons(exons, event_ues)
    }

  } else if (event_type == "IR") {
    # Intron retention
    ir_start <- min(event$five_prime, event$three_prime)
    ir_end <- max(event$five_prime, event$three_prime)

    if (direction == "GAIN") {
      # Comparator retains intron → split into multiple exons (dominant)
      # Find and remove comparator's merged exon that spans the IR region
      for (i in seq_len(nrow(exons))) {
        if (exons$exon_start[i] <= ir_start && exons$exon_end[i] >= ir_end) {
          # This exon spans the entire IR region - remove it
          exons <- exons[-i, ]
          break
        }
      }

      # Use ir_split_exons field to find the correct exonic union exons
      if (!is.na(event$ir_split_exons) && event$ir_split_exons != "") {
        split_exons_ues <- parse_and_find_union_exons(event$ir_split_exons, union_exons, event$gene_id)
        if (nrow(split_exons_ues) > 0) {
          exons <- add_union_exons(exons, split_exons_ues, event$strand)
        } else {
          warning("ir_split_exons specified but no matching union exons found")
        }
      } else {
        warning("IR GAIN event missing ir_split_exons field - cannot reconstruct properly")
      }

    } else {
      # IR LOSS: Dominant has retention, comparator has split exons
      # Strategy: Only FILL IN the INTRONIC gaps between comparator's split exons
      # Don't touch exonic boundaries - let Alt_TSS/Alt_TES handle those

      # Get all union exons in the IR region
      all_ues_in_region <- event_ues

      # Get the comparator's exonic regions (from ir_split_exons)
      if (!is.na(event$ir_split_exons) && event$ir_split_exons != "") {
        exonic_ues <- parse_and_find_union_exons(event$ir_split_exons, union_exons, event$gene_id)

        # Find intronic union exons: those in the IR region but NOT in exonic regions
        if (nrow(exonic_ues) > 0 && nrow(all_ues_in_region) > 0) {
          # Identify intronic UEs by removing exonic UEs from all UEs in region
          intronic_ues <- all_ues_in_region %>%
            anti_join(exonic_ues, by = c("start", "end", "chr", "strand"))

          # Add only the intronic union exons to fill the gaps
          if (nrow(intronic_ues) > 0) {
            exons <- add_union_exons(exons, intronic_ues, event$strand)
          }
        }
      } else {
        # No split exons info - fill entire region (fallback to old behavior)
        merged_ue <- all_ues_in_region %>%
          summarize(
            chr = first(chr),
            start = min(start),
            end = max(end),
            strand = first(strand),
            gene_id = first(gene_id)
          )
        exons <- add_union_exons(exons, merged_ue, event$strand)
      }
    }

  } else if (event_type %in% c("A5SS", "A3SS", "Partial_IR_5", "Partial_IR_3")) {
    # Splice site variations: modify exon boundaries
    exons <- modify_exon_boundary(exons, event_ues, event)
  }

  return(exons)
}

# ==============================================================================
# Main Reconstruction Function
# ==============================================================================

#' Reconstruct dominant isoform from comparator + events (union exon based)
#'
#' @param comparator_exons Exons from comparator isoform
#' @param events Events for this isoform pair
#' @param union_exons Union exon structure for the gene
#' @return Reconstructed dominant exons
reconstruct_dominant_v2 <- function(comparator_exons, events, union_exons) {
  # Start with comparator structure
  reconstructed <- comparator_exons

  # Two-phase reconstruction:
  # Phase 1: Internal events (IR, SE, splice sites) - merge after to consolidate structure
  # Phase 2: Terminal events (Alt_TSS, Alt_TES) - NO merge to keep new exons separate

  # Separate events into internal vs terminal
  internal_events_raw <- events %>%
    filter(event_type != "Dual_boundary",
           event_type %in% c("IR", "SE", "A5SS", "A3SS", "Partial_IR_5", "Partial_IR_3"))

  # Filter out Partial_IR events that overlap with IR events
  # (IR events handle the full retention, Partial_IR becomes redundant)
  ir_events <- internal_events_raw %>% filter(event_type == "IR")
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
        event_type == "IR" ~ 1,
        event_type == "SE" ~ 2,
        TRUE ~ 3
      )
    )

  terminal_events <- events %>%
    filter(event_type %in% c("Alt_TSS", "Alt_TES")) %>%
    arrange(event_type)  # Alt_TSS before Alt_TES

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

  # Final cleanup: sort and deduplicate, but DON'T merge
  reconstructed <- reconstructed %>%
    arrange(exon_start, exon_end) %>%
    distinct(exon_start, exon_end, .keep_all = TRUE)

  return(reconstructed)
}
