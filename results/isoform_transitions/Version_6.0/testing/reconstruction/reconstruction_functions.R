#!/usr/bin/env Rscript
#
# Reconstruction Functions
#
# Functions to reconstruct dominant isoforms from comparator isoforms
# by applying inverse operations for each event type.
#
# Each function takes:
#   - exons: Current exon structure (tibble)
#   - event: Event record (single row from events file)
#   - Returns: Modified exon structure

library(tidyverse)

# ==============================================================================
# Helper Functions
# ==============================================================================

#' Find exon that contains a specific coordinate
#'
#' @param exons Exon tibble
#' @param coord Coordinate to search for
#' @param boundary Which boundary ("start" or "end")
#' @return Row number of matching exon, or NA
find_exon_with_boundary <- function(exons, coord, boundary = "either") {
  if (boundary == "start") {
    matches <- which(exons$exon_start == coord)
  } else if (boundary == "end") {
    matches <- which(exons$exon_end == coord)
  } else {  # either
    matches <- which(exons$exon_start == coord | exons$exon_end == coord)
  }

  if (length(matches) > 0) {
    return(matches[1])
  } else {
    return(NA)
  }
}

#' Find exon that overlaps a coordinate
#'
#' @param exons Exon tibble
#' @param coord Coordinate
#' @return Row number of overlapping exon, or NA
find_exon_containing <- function(exons, coord) {
  matches <- which(exons$exon_start <= coord & exons$exon_end >= coord)

  if (length(matches) > 0) {
    return(matches[1])
  } else {
    return(NA)
  }
}

#' Parse coordinate range string
#'
#' @param range_str String like "1000-1200,1500-1700"
#' @return Tibble with start and end columns
parse_coordinate_ranges <- function(range_str) {
  if (is.na(range_str) || range_str == "") {
    return(tibble(start = integer(), end = integer()))
  }

  ranges <- str_split(range_str, ",")[[1]]

  results <- map_dfr(ranges, function(r) {
    parts <- str_split(r, "-")[[1]]
    tibble(
      start = as.integer(parts[1]),
      end = as.integer(parts[2])
    )
  })

  return(results)
}

# ==============================================================================
# A5SS Reconstruction (Donor Differs)
# ==============================================================================

#' Apply A5SS event to reconstruct dominant
#'
#' A5SS: Donor (5' splice site) differs
#' Event coordinates (five_prime, three_prime) represent the exact region that differs
#'
#' @param exons Comparator exons
#' @param event Event record
#' @return Reconstructed exons
apply_a5ss <- function(exons, event) {
  strand <- event$strand
  direction <- event$direction

  if (strand == "+") {
    # Plus: donor = exon end
    if (direction == "LOSS") {
      # Comparator lost sequence: comparator ends at five_prime-1, extend to three_prime
      exon_idx <- which(abs(exons$exon_end - (event$five_prime - 1)) <= 1)[1]
      if (!is.na(exon_idx)) {
        exons$exon_end[exon_idx] <- event$three_prime
      }
    } else if (direction == "GAIN") {
      # Comparator gained sequence: comparator ends at three_prime, trim to five_prime-1
      exon_idx <- which(abs(exons$exon_end - event$three_prime) <= 1)[1]
      if (!is.na(exon_idx)) {
        exons$exon_end[exon_idx] <- event$five_prime - 1
      }
    }

  } else {  # Minus strand
    # Minus: donor = exon start
    if (direction == "LOSS") {
      # Comparator lost sequence: comparator starts at three_prime+1, extend to five_prime
      exon_idx <- which(abs(exons$exon_start - (event$three_prime + 1)) <= 1)[1]
      if (!is.na(exon_idx)) {
        exons$exon_start[exon_idx] <- event$five_prime
      }
    } else if (direction == "GAIN") {
      # Comparator gained sequence: comparator starts at five_prime, trim to three_prime+1
      exon_idx <- which(abs(exons$exon_start - event$five_prime) <= 1)[1]
      if (!is.na(exon_idx)) {
        exons$exon_start[exon_idx] <- event$three_prime + 1
      }
    }
  }

  return(exons)
}

# ==============================================================================
# A3SS Reconstruction (Acceptor Differs)
# ==============================================================================

#' Apply A3SS event to reconstruct dominant
#'
#' A3SS: Acceptor (3' splice site) differs
#' Event coordinates (five_prime, three_prime) represent the exact region that differs
#'
#' @param exons Comparator exons
#' @param event Event record
#' @return Reconstructed exons
apply_a3ss <- function(exons, event) {
  strand <- event$strand
  direction <- event$direction

  if (strand == "+") {
    # Plus: acceptor = exon start
    if (direction == "LOSS") {
      # Comparator lost sequence: comparator starts at three_prime+1, extend to five_prime
      exon_idx <- which(abs(exons$exon_start - (event$three_prime + 1)) <= 1)[1]
      if (!is.na(exon_idx)) {
        exons$exon_start[exon_idx] <- event$five_prime
      }
    } else if (direction == "GAIN") {
      # Comparator gained sequence: comparator starts at five_prime, trim to three_prime+1
      exon_idx <- which(abs(exons$exon_start - event$five_prime) <= 1)[1]
      if (!is.na(exon_idx)) {
        exons$exon_start[exon_idx] <- event$three_prime + 1
      }
    }

  } else {  # Minus strand
    # Minus: acceptor = exon end
    if (direction == "LOSS") {
      # Comparator lost sequence: comparator ends at three_prime-1, extend to five_prime
      exon_idx <- which(abs(exons$exon_end - (event$three_prime - 1)) <= 1)[1]
      if (!is.na(exon_idx)) {
        exons$exon_end[exon_idx] <- event$five_prime
      }
    } else if (direction == "GAIN") {
      # Comparator gained sequence: comparator ends at five_prime, trim to three_prime-1
      exon_idx <- which(abs(exons$exon_end - event$five_prime) <= 1)[1]
      if (!is.na(exon_idx)) {
        exons$exon_end[exon_idx] <- event$three_prime - 1
      }
    }
  }

  return(exons)
}

# ==============================================================================
# Partial IR Reconstruction (Same as A5SS/A3SS but larger bp_diff)
# ==============================================================================

#' Apply Partial_IR_5 event (partial retention at 5' splice site / donor)
#'
#' @param exons Comparator exons
#' @param event Event record
#' @return Reconstructed exons
apply_partial_ir_5 <- function(exons, event) {
  # Partial_IR_5 affects the donor (same as A5SS but ≥100bp)
  return(apply_a5ss(exons, event))
}

#' Apply Partial_IR_3 event (partial retention at 3' splice site / acceptor)
#'
#' @param exons Comparator exons
#' @param event Event record
#' @return Reconstructed exons
apply_partial_ir_3 <- function(exons, event) {
  # Partial_IR_3 affects the acceptor (same as A3SS but ≥100bp)
  return(apply_a3ss(exons, event))
}

# ==============================================================================
# IR Reconstruction (Intron Retention)
# ==============================================================================

#' Apply IR event to reconstruct dominant
#'
#' IR: One isoform has monoexonic region spanning multiple exons in the other
#'
#' Direction = GAIN: Comparator has retention (monoexonic)
#'   → Split comparator's long exon into multiple exons
#'
#' Direction = LOSS: Dominant has retention (monoexonic)
#'   → Merge comparator's multiple exons into one
#'
#' @param exons Comparator exons
#' @param event Event record
#' @param union_exons Union exon structure (to identify split points)
#' @return Reconstructed exons
apply_ir <- function(exons, event, union_exons = NULL) {
  direction <- event$direction
  strand <- event$strand

  # Get IR boundaries
  if (strand == "+") {
    ir_start <- event$five_prime
    ir_end <- event$three_prime
  } else {
    ir_start <- event$three_prime
    ir_end <- event$five_prime
  }

  if (direction == "GAIN") {
    # Comparator has intron retention → split into multiple exons
    # Find the long exon that spans this region
    exon_idx <- which(
      exons$exon_start <= ir_start &
      exons$exon_end >= ir_end
    )

    if (length(exon_idx) > 0) {
      exon_idx <- exon_idx[1]
      long_exon <- exons[exon_idx, ]

      # Use union exons to determine split points
      if (!is.null(union_exons)) {
        # Find union exons within this region
        ues_in_region <- union_exons %>%
          filter(
            chr == long_exon$chr,
            strand == long_exon$strand,
            start >= long_exon$exon_start,
            end <= long_exon$exon_end
          ) %>%
          arrange(start)

        if (nrow(ues_in_region) > 1) {
          # Split the long exon according to union exon boundaries
          new_exons <- list()

          for (i in seq_len(nrow(ues_in_region))) {
            new_exons[[i]] <- tibble(
              chr = long_exon$chr,
              exon_start = ues_in_region$start[i],
              exon_end = ues_in_region$end[i],
              strand = long_exon$strand,
              gene_id = long_exon$gene_id,
              transcript_id = long_exon$transcript_id
            )
          }

          # Replace long exon with split exons
          before_exons <- if (exon_idx > 1) exons[1:(exon_idx-1), ] else tibble()
          after_exons <- if (exon_idx < nrow(exons)) exons[(exon_idx+1):nrow(exons), ] else tibble()

          exons <- bind_rows(before_exons, bind_rows(new_exons), after_exons)
        }
      }
    }

  } else if (direction == "LOSS") {
    # Dominant has intron retention → merge comparator's multiple exons
    # Find all exons within the IR region
    exons_in_region <- which(
      exons$exon_start >= ir_start &
      exons$exon_end <= ir_end
    )

    if (length(exons_in_region) > 1) {
      # Merge these exons into one
      first_idx <- min(exons_in_region)
      last_idx <- max(exons_in_region)

      merged_exon <- exons[first_idx, ]
      merged_exon$exon_start <- min(exons$exon_start[exons_in_region])
      merged_exon$exon_end <- max(exons$exon_end[exons_in_region])

      # Replace multiple exons with merged exon
      before_exons <- if (first_idx > 1) exons[1:(first_idx-1), ] else tibble()
      after_exons <- if (last_idx < nrow(exons)) exons[(last_idx+1):nrow(exons), ] else tibble()

      exons <- bind_rows(before_exons, merged_exon, after_exons)
    }
  }

  return(exons)
}

# ==============================================================================
# SE Reconstruction (Skipped Exon)
# ==============================================================================

#' Apply SE event to reconstruct dominant
#'
#' SE: Exon present in one isoform but absent in the other
#'
#' Typically, dominant has the exon and comparator lacks it.
#' Reconstruction: Insert the skipped exon between flanking exons.
#'
#' @param exons Comparator exons
#' @param event Event record
#' @return Reconstructed exons
apply_se <- function(exons, event) {
  strand <- event$strand

  # SE coordinates store the missing exon boundaries
  if (strand == "+") {
    se_start <- event$five_prime
    se_end <- event$three_prime
  } else {
    se_start <- event$three_prime
    se_end <- event$five_prime
  }

  # Find insertion point (between flanking exons)
  # Insert after the exon that ends before se_start
  insert_after <- max(which(exons$exon_end < se_start), 0)

  # Create new exon
  new_exon <- tibble(
    chr = exons$chr[1],
    exon_start = se_start,
    exon_end = se_end,
    strand = exons$strand[1],
    gene_id = exons$gene_id[1],
    transcript_id = exons$transcript_id[1]
  )

  # Insert the exon
  if (insert_after == 0) {
    # Insert at beginning
    exons <- bind_rows(new_exon, exons)
  } else if (insert_after >= nrow(exons)) {
    # Insert at end
    exons <- bind_rows(exons, new_exon)
  } else {
    # Insert in middle
    before_exons <- exons[1:insert_after, ]
    after_exons <- exons[(insert_after+1):nrow(exons), ]
    exons <- bind_rows(before_exons, new_exon, after_exons)
  }

  return(exons)
}

# ==============================================================================
# Alt_TSS Reconstruction (Alternative Transcription Start Site)
# ==============================================================================

#' Apply Alt_TSS event to reconstruct dominant
#'
#' Alt_TSS: Transcription start site differs
#' Event coordinates (five_prime, three_prime) represent the exact region that differs
#'
#' @param exons Comparator exons (ordered TSS → TES)
#' @param event Event record
#' @return Reconstructed exons
apply_alt_tss <- function(exons, event) {
  direction <- event$direction
  strand <- event$strand

  if (strand == "+") {
    # Plus strand: TSS = exon start (lower coordinate)
    if (direction == "LOSS") {
      # Comparator lost sequence: extend first exon to include lost region
      # Comparator starts at three_prime+1, extend to five_prime
      if (nrow(exons) > 0 && abs(exons$exon_start[1] - (event$three_prime + 1)) <= 1) {
        exons$exon_start[1] <- event$five_prime
      }
    } else if (direction == "GAIN") {
      # Comparator gained sequence: trim first exon to exclude gained region
      # Comparator starts at five_prime, trim to three_prime+1
      if (nrow(exons) > 0 && abs(exons$exon_start[1] - event$five_prime) <= 1) {
        exons$exon_start[1] <- event$three_prime + 1
      }
    }

  } else {
    # Minus strand: TSS = exon end (higher coordinate)
    if (direction == "LOSS") {
      # Comparator lost sequence: extend first exon to include lost region
      # Comparator ends at three_prime-1, extend to five_prime
      if (nrow(exons) > 0 && abs(exons$exon_end[1] - (event$three_prime - 1)) <= 1) {
        exons$exon_end[1] <- event$five_prime
      }
    } else if (direction == "GAIN") {
      # Comparator gained sequence: trim first exon to exclude gained region
      # Comparator ends at five_prime, trim to three_prime-1
      if (nrow(exons) > 0 && abs(exons$exon_end[1] - event$five_prime) <= 1) {
        exons$exon_end[1] <- event$three_prime - 1
      }
    }
  }

  return(exons)
}

# ==============================================================================
# Alt_TES Reconstruction (Alternative Transcription End Site)
# ==============================================================================

#' Apply Alt_TES event to reconstruct dominant
#'
#' Alt_TES: Transcription end site differs
#' Event coordinates (five_prime, three_prime) represent the exact region that differs
#'
#' @param exons Comparator exons (ordered TSS → TES)
#' @param event Event record
#' @return Reconstructed exons
apply_alt_tes <- function(exons, event) {
  direction <- event$direction
  strand <- event$strand

  if (strand == "+") {
    # Plus strand: TES = exon end (higher coordinate)
    if (direction == "LOSS") {
      # Comparator lost sequence: extend last exon to include lost region
      # Comparator ends at five_prime-1, extend to three_prime
      last_idx <- nrow(exons)
      if (last_idx > 0 && abs(exons$exon_end[last_idx] - (event$five_prime - 1)) <= 1) {
        exons$exon_end[last_idx] <- event$three_prime
      }
    } else if (direction == "GAIN") {
      # Comparator gained sequence: trim last exon to exclude gained region
      # Comparator ends at three_prime, trim to five_prime-1
      last_idx <- nrow(exons)
      if (last_idx > 0 && abs(exons$exon_end[last_idx] - event$three_prime) <= 1) {
        exons$exon_end[last_idx] <- event$five_prime - 1
      }
    }

  } else {
    # Minus strand: TES = exon start (lower coordinate)
    if (direction == "LOSS") {
      # Comparator lost sequence: extend last exon to include lost region
      # Comparator starts at three_prime+1, extend to five_prime
      last_idx <- nrow(exons)
      if (last_idx > 0 && abs(exons$exon_start[last_idx] - (event$three_prime + 1)) <= 1) {
        exons$exon_start[last_idx] <- event$five_prime
      }
    } else if (direction == "GAIN") {
      # Comparator gained sequence: trim last exon to exclude gained region
      # Comparator starts at five_prime, trim to three_prime+1
      last_idx <- nrow(exons)
      if (last_idx > 0 && abs(exons$exon_start[last_idx] - event$five_prime) <= 1) {
        exons$exon_start[last_idx] <- event$three_prime + 1
      }
    }
  }

  return(exons)
}

# ==============================================================================
# Dual Boundary Event (SKIP for reconstruction)
# ==============================================================================

#' Handle Dual_boundary events (both boundaries differ)
#'
#' These are complex and typically handled by component events
#' (e.g., combination of A5SS and A3SS)
#'
#' @param exons Comparator exons
#' @param event Event record
#' @return Exons unchanged (skip)
apply_dual_boundary <- function(exons, event) {
  # Dual boundary events are typically decomposed into other events
  # Skip for now
  return(exons)
}

# ==============================================================================
# Main Reconstruction Function
# ==============================================================================

#' Reconstruct dominant isoform from comparator + events
#'
#' @param comparator_exons Exons from comparator isoform
#' @param events Events for this isoform pair (filtered to one transcript pair)
#' @param union_exons Union exon structure (optional, for IR reconstruction)
#' @return Reconstructed dominant exons
reconstruct_dominant <- function(comparator_exons, events, union_exons = NULL) {

  # Start with comparator structure
  reconstructed <- comparator_exons

  # Sort events by priority:
  # 1. Terminal events first (Alt_TSS, Alt_TES) - these define boundaries
  # 2. Internal events (A5SS, A3SS, Partial_IR, IR, SE)
  # 3. Skip Dual_boundary (handled by component events)
  events <- events %>%
    filter(event_type != "Dual_boundary") %>%
    arrange(
      case_when(
        event_type == "Alt_TSS" ~ 1,
        event_type == "Alt_TES" ~ 2,
        event_type == "IR" ~ 3,
        event_type == "SE" ~ 4,
        TRUE ~ 5
      )
    )

  # Apply each event
  for (i in seq_len(nrow(events))) {
    event <- events[i, ]

    reconstructed <- tryCatch({
      if (event$event_type == "A5SS") {
        apply_a5ss(reconstructed, event)

      } else if (event$event_type == "A3SS") {
        apply_a3ss(reconstructed, event)

      } else if (event$event_type == "Partial_IR_5") {
        apply_partial_ir_5(reconstructed, event)

      } else if (event$event_type == "Partial_IR_3") {
        apply_partial_ir_3(reconstructed, event)

      } else if (event$event_type == "IR") {
        apply_ir(reconstructed, event, union_exons)

      } else if (event$event_type == "SE") {
        apply_se(reconstructed, event)

      } else if (event$event_type == "Alt_TSS") {
        apply_alt_tss(reconstructed, event)

      } else if (event$event_type == "Alt_TES") {
        apply_alt_tes(reconstructed, event)

      } else {
        # Unknown or Dual_boundary - skip
        reconstructed
      }
    }, error = function(e) {
      warning(sprintf("Error applying %s event: %s", event$event_type, e$message))
      reconstructed  # Return unchanged on error
    })
  }

  # Sort exons by genomic position
  reconstructed <- reconstructed %>%
    arrange(exon_start, exon_end)

  return(reconstructed)
}
