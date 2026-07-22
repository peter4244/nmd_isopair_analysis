# Configuration
output_dir <- "results/isoform_transitions/v4.0_reference_based"
BATCH_SIZE <- 1000  # Save every 1000 genes
TSS_TES_TOLERANCE <- 20  # Minimum bp difference to call Alt_TSS/Alt_TES

# ============================================================================
# Helper Functions
# ============================================================================

#' Extract isoform's exons from union model
get_isoform_exons <- function(union_exons, isoform_id) {
  exon_records <- map_dfr(union_exons, function(union_exon) {
    variants_in_iso <- union_exon$variants %>%
      filter(isoform_id == !!isoform_id)

    if (nrow(variants_in_iso) == 0) {
      return(NULL)
    }

    variant <- variants_in_iso[1, ]

    tibble(
      exon_number = union_exon$exon_number,
      exon_type = union_exon$exon_type,
      start = variant$start,
      end = variant$end,
      is_first = variant$is_first,
      is_last = variant$is_last
    )
  })

  exon_records
}

#' Compare two isoforms using union model to detect events
detect_events_from_union <- function(union_exons, isoform_A_id, isoform_B_id, gene_strand) {
  exons_A <- get_isoform_exons(union_exons, isoform_A_id)
  exons_B <- get_isoform_exons(union_exons, isoform_B_id)

  if (nrow(exons_A) == 0 || nrow(exons_B) == 0) {
    return(tibble(event_type = character(),
                  exon_number = integer(),
                  direction = character(),
                  detail = character()))
  }

  events <- list()

  # Get first and last exons from each isoform
  first_A <- exons_A %>% filter(row_number() == 1)
  first_B <- exons_B %>% filter(row_number() == 1)
  last_A <- exons_A %>% filter(row_number() == n())
  last_B <- exons_B %>% filter(row_number() == n())

  # SKIP SINGLE-EXON GENE COMPARISONS
  # If both isoforms are single-exon (first == last), there's no splicing to analyze
  is_single_exon_A <- (first_A$exon_number == last_A$exon_number)
  is_single_exon_B <- (first_B$exon_number == last_B$exon_number)

  if (is_single_exon_A && is_single_exon_B) {
    # Both isoforms are single-exon - no splicing events possible
    return(tibble(event_type = character(),
                  exon_number = integer(),
                  direction = character(),
                  detail = character()))
  }

  # EXPLICIT Alt_TSS/Alt_TES DETECTION (outside union exon iteration)
  # This is necessary because first/last exons that differ by >20bp won't be
  # in the same union exon group, so we need to compare them explicitly

  # Check Alt_TSS (compare first exons)
  if (gene_strand == "+") {
    # Plus strand: TSS = start coordinate
    tss_diff <- abs(first_A$start - first_B$start)
    if (tss_diff > TSS_TES_TOLERANCE) {
      events[[length(events) + 1]] <- tibble(
        event_type = "Alt_TSS",
        exon_number = NA_integer_,  # First exons may be at different union positions
        direction = "none",  # No gain/loss for Alt_TSS - just different positions
        detail = sprintf("TSS diff: %dbp (A: %d, B: %d)",
                        tss_diff, first_A$start, first_B$start)
      )
    }
  } else if (gene_strand == "-") {
    # Minus strand: TSS = end coordinate
    tss_diff <- abs(first_A$end - first_B$end)
    if (tss_diff > TSS_TES_TOLERANCE) {
      events[[length(events) + 1]] <- tibble(
        event_type = "Alt_TSS",
        exon_number = NA_integer_,  # First exons may be at different union positions
        direction = "none",  # No gain/loss for Alt_TSS - just different positions
        detail = sprintf("TSS diff: %dbp (A: %d, B: %d)",
                        tss_diff, first_A$end, first_B$end)
      )
    }
  }

  # Check Alt_TES (compare last exons)
  if (gene_strand == "+") {
    # Plus strand: TES = end coordinate
    tes_diff <- abs(last_A$end - last_B$end)
    if (tes_diff > TSS_TES_TOLERANCE) {
      events[[length(events) + 1]] <- tibble(
        event_type = "Alt_TES",
        exon_number = NA_integer_,  # Last exons may be at different union positions
        direction = "none",  # No gain/loss for Alt_TES - just different positions
        detail = sprintf("TES diff: %dbp (A: %d, B: %d)",
                        tes_diff, last_A$end, last_B$end)
      )
    }
  } else if (gene_strand == "-") {
    # Minus strand: TES = start coordinate
    tes_diff <- abs(last_A$start - last_B$start)
    if (tes_diff > TSS_TES_TOLERANCE) {
      events[[length(events) + 1]] <- tibble(
        event_type = "Alt_TES",
        exon_number = NA_integer_,  # Last exons may be at different union positions
        direction = "none",  # No gain/loss for Alt_TES - just different positions
        detail = sprintf("TES diff: %dbp (A: %d, B: %d)",
                        tes_diff, last_A$start, last_B$start)
      )
    }
  }

  for (union_exon in union_exons) {
    exon_num <- union_exon$exon_number
    exon_type <- union_exon$exon_type

    variants_A <- union_exon$variants %>% filter(isoform_id == isoform_A_id)
    variants_B <- union_exon$variants %>% filter(isoform_id == isoform_B_id)

    has_A <- nrow(variants_A) > 0
    has_B <- nrow(variants_B) > 0

    # Case 1: Both isoforms have this exon
    if (has_A && has_B) {
      var_A <- variants_A[1, ]
      var_B <- variants_B[1, ]

      same_start <- var_A$start == var_B$start
      same_end <- var_A$end == var_B$end

      if (same_start && same_end) {
        # Constitutive exon
        events[[length(events) + 1]] <- tibble(
          event_type = "CONST",
          exon_number = exon_num,
          direction = "none",
          detail = sprintf("Exon %d identical", exon_num)
        )

      } else {
        # NOTE: Alt_TSS and Alt_TES are now handled explicitly above,
        # outside the union exon iteration, to avoid the grouping issue

        # CRITICAL: Skip terminal boundary differences for first/last exons
        # to prevent double-counting with explicit Alt_TSS/Alt_TES detection
        is_first_exon <- var_A$is_first || var_B$is_first
        is_last_exon <- var_A$is_last || var_B$is_last

        # For first exons: only process internal boundary (not TSS)
        # For last exons: only process internal boundary (not TES)
        # For single-exon genes (both first AND last): skip entirely
        if (is_first_exon && is_last_exon) {
          # Single-exon gene - all boundaries handled by Alt_TSS/Alt_TES
          next
        }

        # Adjust which boundaries to process for terminal exons
        # We skip TSS and TES boundaries (handled by explicit Alt_TSS/Alt_TES detection)
        # But we DO process internal splice sites (donors and acceptors)
        if (gene_strand == "+") {
          # Plus strand: TSS at start, TES at end
          process_start <- !is_first_exon  # Skip TSS (first exon start)
          process_end <- !is_last_exon     # Skip TES (last exon end)
        } else {
          # Minus strand: coordinates reversed - TES at start, TSS at end
          process_start <- !is_last_exon   # Skip TES (last exon start on minus)
          process_end <- !is_first_exon    # Skip TSS (first exon end on minus)
        }

        # Apply boundary processing restrictions
        check_start <- !same_start && process_start
        check_end <- !same_end && process_end

        # Internal exon with different coordinates - STRAND-AWARE
        # Key: On minus strand, event types are SWAPPED
        # Plus: different start = A3SS, different end = A5SS
        # Minus: different start = A5SS, different end = A3SS

        if (check_start && check_end) {
          # Both boundaries differ - always A5SS + A3SS
          if (gene_strand == "+") {
            direction_5 <- if_else(var_A$end > var_B$end, "gain", "loss")
            direction_3 <- if_else(var_A$start < var_B$start, "gain", "loss")
          } else {
            direction_5 <- if_else(var_A$start < var_B$start, "gain", "loss")
            direction_3 <- if_else(var_A$end > var_B$end, "gain", "loss")
          }

          events[[length(events) + 1]] <- tibble(
            event_type = "A5SS",
            exon_number = exon_num,
            direction = direction_5,
            detail = sprintf("5'SS diff")
          )

          events[[length(events) + 1]] <- tibble(
            event_type = "A3SS",
            exon_number = exon_num,
            direction = direction_3,
            detail = sprintf("3'SS diff")
          )

        } else if (check_start) {
          # Different start only
          if (gene_strand == "+") {
            # Plus: start = 3'SS (acceptor)
            direction <- if_else(var_A$start < var_B$start, "gain", "loss")
            event_name <- "A3SS"
          } else {
            # Minus: start = 5'SS (donor) - SWAPPED!
            direction <- if_else(var_A$start < var_B$start, "gain", "loss")
            event_name <- "A5SS"
          }

          events[[length(events) + 1]] <- tibble(
            event_type = event_name,
            exon_number = exon_num,
            direction = direction,
            detail = sprintf("%s: %d vs %d", event_name, var_A$start, var_B$start)
          )

        } else if (check_end) {
          # Different end only
          if (gene_strand == "+") {
            # Plus: end = 5'SS (donor)
            direction <- if_else(var_A$end > var_B$end, "gain", "loss")
            event_name <- "A5SS"
          } else {
            # Minus: end = 3'SS (acceptor) - SWAPPED!
            direction <- if_else(var_A$end > var_B$end, "gain", "loss")
            event_name <- "A3SS"
          }

          events[[length(events) + 1]] <- tibble(
            event_type = event_name,
            exon_number = exon_num,
            direction = direction,
            detail = sprintf("%s: %d vs %d", event_name, var_A$end, var_B$end)
          )
        }
      }

    } else if (has_A && !has_B) {
      # Case 2: Only A has this exon
      # Check if it's a first/last exon (already handled as Alt_TSS/TES)
      var_A <- variants_A[1, ]
      if (!var_A$is_first && !var_A$is_last) {
        # Only call SE for internal exons
        events[[length(events) + 1]] <- tibble(
          event_type = "SE",
          exon_number = exon_num,
          direction = "gain",
          detail = sprintf("Exon %d gained", exon_num)
        )
      }

    } else if (!has_A && has_B) {
      # Case 3: Only B has this exon
      # Check if it's a first/last exon (already handled as Alt_TSS/TES)
      var_B <- variants_B[1, ]
      if (!var_B$is_first && !var_B$is_last) {
        # Only call SE for internal exons
        events[[length(events) + 1]] <- tibble(
          event_type = "SE",
          exon_number = exon_num,
          direction = "loss",
          detail = sprintf("Exon %d lost", exon_num)
        )
      }
    }
  }

  if (length(events) == 0) {
    return(tibble(event_type = character(),
                  exon_number = integer(),
                  direction = character(),
                  detail = character()))
  }

  bind_rows(events)
}

cat("Helper functions loaded\n\n")

# ============================================================================
# Load Data
