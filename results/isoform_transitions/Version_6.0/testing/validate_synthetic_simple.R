#!/usr/bin/env Rscript
#
# Synthetic Gold Standard Validation (Simplified)
# Avoids Bioconductor namespace conflicts by using manual GTF parsing
#

library(tidyverse)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   GOLD STANDARD VALIDATION - Synthetic Test Data             ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# ==============================================================================
# Helper: Parse GTF manually
# ==============================================================================

parse_gtf <- function(file) {
  lines <- readLines(file)
  lines <- lines[!startsWith(lines, '#') & lines != '']

  results <- list()
  for (line in lines) {
    parts <- str_split(line, '\t')[[1]]

    # Extract attributes
    attrs_str <- parts[9]
    gene_id <- str_match(attrs_str, 'gene_id "([^"]+)"')[,2]
    transcript_id <- str_match(attrs_str, 'transcript_id "([^"]+)"')[,2]
    exon_number <- str_match(attrs_str, 'exon_number "([^"]+)"')[,2]

    results[[length(results) + 1]] <- tibble(
      seqnames = parts[1],
      type = parts[3],
      start = as.integer(parts[4]),
      end = as.integer(parts[5]),
      strand = parts[7],
      gene_id = gene_id,
      transcript_id = transcript_id,
      exon_number = as.integer(exon_number)
    )
  }

  bind_rows(results)
}

# ==============================================================================
# Source event detection functions
# ==============================================================================

cat("Loading event detection functions...\n")

# Source shared event detection library (used by Script 07 and validation)
source("../scripts/event_detection_functions.R")  # v1

# Legacy wrapper functions for backward compatibility with validation script
detect_a5ss <- function(exon_dom, exon_non_dom, strand) {
  result <- detect_shared_boundary_event(exon_dom, exon_non_dom, strand)
  return(list(detected = (result$event_type == "A5SS"), bp_diff = result$bp_diff))
}

detect_a3ss <- function(exon_dom, exon_non_dom, strand) {
  result <- detect_shared_boundary_event(exon_dom, exon_non_dom, strand)
  return(list(detected = (result$event_type == "A3SS"), bp_diff = result$bp_diff))
}

detect_partial_ir <- function(exon_dom, exon_non_dom, strand, is_first = FALSE, is_last = FALSE, has_overlap = FALSE) {
  result <- detect_shared_boundary_event(exon_dom, exon_non_dom, strand,
                                         is_first_exon = is_first,
                                         is_last_exon = is_last,
                                         terminal_has_overlap = has_overlap)
  return(list(detected = (result$event_type == "Partial_IR"), bp_diff = result$bp_diff))
}

# Note: detect_se for validation uses different logic than Script 07
# (checks both directions A vs B and B vs A)
detect_se <- function(struct_a, struct_b) {
  # SE: exon present in one isoform but not the other, flanked by comparable exons
  # Comparable = exact match OR overlap
  # Check both directions (A vs B and B vs A)

  # Helper: check if two exons overlap
  exons_overlap <- function(ex1_start, ex1_end, ex2_start, ex2_end) {
    return(ex1_start <= ex2_end && ex2_start <= ex1_end)
  }

  # Helper: check if exon exists in other isoform (exact or overlap)
  exon_comparable <- function(exon, other_struct) {
    for (j in seq_len(nrow(other_struct))) {
      other_exon <- other_struct[j, ]
      # Exact match OR overlap
      if ((exon$exon_start == other_exon$exon_start && exon$exon_end == other_exon$exon_end) ||
          exons_overlap(exon$exon_start, exon$exon_end, other_exon$exon_start, other_exon$exon_end)) {
        return(TRUE)
      }
    }
    return(FALSE)
  }

  n_se <- 0

  # Check exons in A that are missing in B
  for (i in seq_len(nrow(struct_a))) {
    exon_a <- struct_a[i, ]

    # Is this exon in B?
    in_b <- exon_comparable(exon_a, struct_b)

    if (!in_b && i > 1 && i < nrow(struct_a)) {
      # This is an internal exon not in B - check flanks are comparable
      prev_exon <- struct_a[i-1, ]
      next_exon <- struct_a[i+1, ]

      prev_comparable <- exon_comparable(prev_exon, struct_b)
      next_comparable <- exon_comparable(next_exon, struct_b)

      if (prev_comparable && next_comparable) n_se <- n_se + 1
    }
  }

  # Check exons in B that are missing in A
  for (i in seq_len(nrow(struct_b))) {
    exon_b <- struct_b[i, ]

    # Is this exon in A?
    in_a <- exon_comparable(exon_b, struct_a)

    if (!in_a && i > 1 && i < nrow(struct_b)) {
      # This is an internal exon not in A - check flanks are comparable
      prev_exon <- struct_b[i-1, ]
      next_exon <- struct_b[i+1, ]

      prev_comparable <- exon_comparable(prev_exon, struct_a)
      next_comparable <- exon_comparable(next_exon, struct_a)

      if (prev_comparable && next_comparable) n_se <- n_se + 1
    }
  }

  return(n_se)
}

cat("  ✓ Functions loaded\n\n")

# ==============================================================================
# Load test data
# ==============================================================================

cat("Loading synthetic test data...\n")

gtf_df <- parse_gtf("synthetic/TestData/exons/base_events.gtf")
expected <- read_tsv("synthetic/TestData/annotations/base_events.tsv",
                     comment = "#", show_col_types = FALSE)

cat(sprintf("  Parsed %d exon features\n", nrow(gtf_df)))
cat(sprintf("  Test genes: %d\n", length(unique(gtf_df$gene_id))))
cat(sprintf("  Expected events: %d\n\n", nrow(expected)))

# ==============================================================================
# Build structures
# ==============================================================================

cat("Building isoform structures...\n")

isoform_structures <- gtf_df %>%
  group_by(gene_id, transcript_id) %>%
  summarise(
    seqnames = first(seqnames),
    strand = first(strand),
    exons = list({
      exon_data <- tibble(
        exon_number = exon_number,
        exon_start = start,
        exon_end = end
      )
      # Use canonical biological ordering based on strand and coordinates
      order_exons_biological(exon_data, first(strand))
    }),
    .groups = "drop"
  )

cat(sprintf("  Created %d isoform structures\n\n", nrow(isoform_structures)))

# ==============================================================================
# Run validation on each test case
# ==============================================================================

cat("Running event detection on test cases...\n")
cat("═══════════════════════════════════════════════════════════\n\n")

results <- list()

for (i in seq_len(nrow(expected))) {
  test_case <- expected[i, ]

  cat(sprintf("[%d/%d] %s\n", i, nrow(expected), test_case$gene_id))
  cat(sprintf("      Expected: %s\n", test_case$event_type))

  iso_a <- test_case$isoform_A
  iso_b <- test_case$isoform_B

  # Get structures
  struct_a_data <- isoform_structures %>%
    filter(transcript_id == iso_a) %>%
    pull(exons) %>%
    .[[1]]

  struct_b_data <- isoform_structures %>%
    filter(transcript_id == iso_b) %>%
    pull(exons) %>%
    .[[1]]

  strand <- isoform_structures %>%
    filter(transcript_id == iso_a) %>%
    pull(strand)

  if (nrow(struct_a_data) == 0 || nrow(struct_b_data) == 0) {
    cat("      ⚠ Missing data\n\n")
    next
  }

  # TSS/TES (independent of splice site detection)
  # IMPORTANT: Structures are arranged by exon_number (biological order), not coordinate
  # So row 1 = biological first exon (TSS), last row = biological last exon (TES)
  # This is true for BOTH strands!
  tss_changed <- detect_tss_change(
    struct_a_data[1, ],
    struct_b_data[1, ],
    strand
  )
  tes_changed <- detect_tes_change(
    struct_a_data[nrow(struct_a_data), ],
    struct_b_data[nrow(struct_b_data), ],
    strand
  )

  # Check if first/last exons overlap (required for terminal splice site detection)
  first_exons_overlap <- calculate_overlap(
    struct_a_data[1, ],
    struct_b_data[1, ]
  )

  last_exons_overlap <- calculate_overlap(
    struct_a_data[nrow(struct_a_data), ],
    struct_b_data[nrow(struct_b_data), ]
  )

  # Count events and track directions
  n_a5ss <- 0
  n_a3ss <- 0
  n_partial_ir_5 <- 0
  n_partial_ir_3 <- 0
  n_ir <- 0
  n_se <- 0
  n_dual_boundary <- 0

  # Track event details with directions
  event_details <- list()

  # UNIFIED EVENT DETECTION (matching Script 07)
  # For each pair of exons, use detect_shared_boundary_event()
  # Exons are now in biological order (TSS → TES) via order_exons_biological()
  for (j in seq_len(nrow(struct_b_data))) {
    non_dom_exon <- struct_b_data[j, ]
    # Use biological_exon_number to determine position (1 = TSS, max = TES)
    is_first_non_dom <- (non_dom_exon$biological_exon_number == 1)
    is_last_non_dom <- (non_dom_exon$biological_exon_number == max(struct_b_data$biological_exon_number))

    for (k in seq_len(nrow(struct_a_data))) {
      dom_exon <- struct_a_data[k, ]
      is_first_dom <- (dom_exon$biological_exon_number == 1)
      is_last_dom <- (dom_exon$biological_exon_number == max(struct_a_data$biological_exon_number))

      # Determine if both exons are terminal and if they overlap
      is_both_first <- is_first_dom && is_first_non_dom
      is_both_last <- is_last_dom && is_last_non_dom
      terminal_overlap <- if (is_both_first) {
        first_exons_overlap
      } else if (is_both_last) {
        last_exons_overlap
      } else {
        FALSE
      }

      # Get flanking exons for IR boundary check
      flanking_dom <- if (nrow(struct_a_data) > 1) {
        struct_a_data %>% filter(biological_exon_number != dom_exon$biological_exon_number)
      } else {
        struct_a_data[0, , drop = FALSE]
      }

      flanking_non_dom <- if (nrow(struct_b_data) > 1) {
        struct_b_data %>% filter(biological_exon_number != non_dom_exon$biological_exon_number)
      } else {
        struct_b_data[0, , drop = FALSE]
      }

      # UNIFIED EVENT DETECTION
      # No strand conversion needed - exons are already in biological order
      event_result <- detect_shared_boundary_event(
        dom_exon, non_dom_exon, strand,
        is_first_exon = is_both_first,
        is_last_exon = is_both_last,
        is_first_exon_comp = is_first_non_dom,
        is_last_exon_comp = is_last_non_dom,
        terminal_has_overlap = terminal_overlap,
        flanking_exons_dom = flanking_dom,
        flanking_exons_non_dom = flanking_non_dom
      )

      # Apply terminal exon rules and count events
      if (event_result$event_type == "A5SS") {
        # A5SS: donor differs, acceptor shared
        # RULE: Cannot occur on last exons (TES is not a donor)
        if (is_both_last) {
          # Skip - TES is not a splice site
        } else if (is_both_first) {
          # First exons: only count if overlap ≥50%
          if (first_exons_overlap) {
            n_a5ss <- n_a5ss + 1
            event_details[[length(event_details) + 1]] <- list(
              type = "A5SS",
              direction = event_result$direction,
              bp_diff = event_result$bp_diff
            )
          }
        } else {
          # Internal exons: always count
          n_a5ss <- n_a5ss + 1
          event_details[[length(event_details) + 1]] <- list(
            type = "A5SS",
            direction = event_result$direction,
            bp_diff = event_result$bp_diff
          )
        }

      } else if (event_result$event_type == "A3SS") {
        # A3SS: acceptor differs, donor shared
        # RULE: Cannot occur on first exons (TSS is not an acceptor)
        if (is_both_first) {
          # Skip - TSS is not a splice site
        } else if (is_both_last) {
          # Last exons: only count if overlap ≥50%
          if (last_exons_overlap) {
            n_a3ss <- n_a3ss + 1
            event_details[[length(event_details) + 1]] <- list(
              type = "A3SS",
              direction = event_result$direction,
              bp_diff = event_result$bp_diff
            )
          }
        } else {
          # Internal exons: always count
          n_a3ss <- n_a3ss + 1
          event_details[[length(event_details) + 1]] <- list(
            type = "A3SS",
            direction = event_result$direction,
            bp_diff = event_result$bp_diff
          )
        }

      } else if (event_result$event_type == "Partial_IR_5") {
        # Partial_IR_5: ≥100bp difference at 5' splice site (donor)
        n_partial_ir_5 <- n_partial_ir_5 + 1
        event_details[[length(event_details) + 1]] <- list(
          type = "Partial_IR_5",
          direction = event_result$direction,
          bp_diff = event_result$bp_diff
        )

      } else if (event_result$event_type == "Partial_IR_3") {
        # Partial_IR_3: ≥100bp difference at 3' splice site (acceptor)
        n_partial_ir_3 <- n_partial_ir_3 + 1
        event_details[[length(event_details) + 1]] <- list(
          type = "Partial_IR_3",
          direction = event_result$direction,
          bp_diff = event_result$bp_diff
        )

      } else if (event_result$event_type == "IR") {
        # IR: detected via overlap-based logic (comparison exon spans into flanking)
        n_ir <- n_ir + 1
        event_details[[length(event_details) + 1]] <- list(
          type = "IR",
          direction = event_result$direction,
          bp_diff = event_result$bp_diff
        )

      } else if (event_result$event_type == "Dual_boundary") {
        # Dual_boundary: both boundaries differ
        n_dual_boundary <- n_dual_boundary + 1
        event_details[[length(event_details) + 1]] <- list(
          type = "Dual_boundary",
          direction = event_result$direction,
          bp_diff = event_result$bp_diff
        )
      }
    }
  }

  # IR (monoexonic vs multi-exonic spanning)
  for (j in seq_len(nrow(struct_b_data))) {
    if (detect_ir_simple(struct_b_data[j, ], struct_a_data)) {
      n_ir <- n_ir + 1
      # Comparison isoform (B) has the retained intron = GAIN
      event_details[[length(event_details) + 1]] <- list(
        type = "IR",
        direction = "GAIN",
        bp_diff = NA
      )
    }
  }
  for (k in seq_len(nrow(struct_a_data))) {
    if (detect_ir_simple(struct_a_data[k, ], struct_b_data)) {
      n_ir <- n_ir + 1
      # Dominant isoform (A) has the retained intron = comparison LOSS
      event_details[[length(event_details) + 1]] <- list(
        type = "IR",
        direction = "LOSS",
        bp_diff = NA
      )
    }
  }

  # SE (Skipped Exon)
  n_se <- detect_se(struct_a_data, struct_b_data)

  # Determine detected event
  detected_events <- character()
  if (tss_changed) detected_events <- c(detected_events, "Alt_TSS")
  if (tes_changed) detected_events <- c(detected_events, "Alt_TES")
  if (n_a5ss > 0) detected_events <- c(detected_events, "A5SS")
  if (n_a3ss > 0) detected_events <- c(detected_events, "A3SS")
  if (n_partial_ir_5 > 0) detected_events <- c(detected_events, "Partial_IR_5")
  if (n_partial_ir_3 > 0) detected_events <- c(detected_events, "Partial_IR_3")
  if (n_ir > 0) detected_events <- c(detected_events, "IR")
  if (n_se > 0) detected_events <- c(detected_events, "SE")

  detected_str <- if (length(detected_events) == 0) "none" else paste(detected_events, collapse = ",")

  # Check match (order-independent for multi-event cases)
  if (test_case$event_type == "none") {
    match <- detected_str == "none"
  } else if (test_case$event_type == "IR_monoexonic") {
    match <- n_ir > 0
  } else if (test_case$event_type == "Mixed") {
    match <- length(detected_events) > 0
  } else {
    expected_events <- str_split(test_case$event_type, ",")[[1]]
    match <- setequal(detected_events, expected_events)
  }

  # Format event details with directions
  event_details_str <- ""
  if (length(event_details) > 0) {
    details_parts <- sapply(event_details, function(e) {
      if (!is.null(e$direction) && e$direction != "") {
        sprintf("%s:%s", e$type, e$direction)
      } else {
        e$type
      }
    })
    event_details_str <- paste(details_parts, collapse = ",")
  }

  cat(sprintf("      Detected: %s\n", detected_str))
  if (event_details_str != "") {
    cat(sprintf("      Details: %s\n", event_details_str))
  }
  cat(sprintf("      Status: %s\n\n", ifelse(match, "✓ PASS", "✗ FAIL")))

  results[[i]] <- tibble(
    gene_id = test_case$gene_id,
    expected = test_case$event_type,
    detected = detected_str,
    event_details = event_details_str,
    match = match
  )
}

results_df <- bind_rows(results)

# ==============================================================================
# Summary
# ==============================================================================

cat("═══════════════════════════════════════════════════════════════════\n")
cat("VALIDATION SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

n_total <- nrow(results_df)
n_pass <- sum(results_df$match, na.rm = TRUE)
n_fail <- n_total - n_pass

cat(sprintf("Total test cases: %d\n", n_total))
cat(sprintf("✓ PASS: %d (%.1f%%)\n", n_pass, 100 * n_pass / n_total))
cat(sprintf("✗ FAIL: %d (%.1f%%)\n\n", n_fail, 100 * n_fail / n_total))

if (n_fail > 0) {
  cat("FAILED TEST CASES:\n")
  cat("─────────────────────────────────────────────────────────\n")
  failures <- results_df %>% filter(!match)
  print(failures, n = Inf)
  cat("\n")
}

cat("═══════════════════════════════════════════════════════════════════\n")

if (n_pass == n_total) {
  cat("✓✓✓ VALIDATION PASSED - All events correctly detected! ✓✓✓\n")
} else {
  cat("⚠ VALIDATION INCOMPLETE - Review failures above\n")
}

cat("═══════════════════════════════════════════════════════════════════\n\n")

# Save results
saveRDS(results_df, "synthetic_validation_results.rds")
write_csv(results_df, "synthetic_validation_results.csv")
cat("✓ Results saved to:\n")
cat("  - synthetic_validation_results.rds\n")
cat("  - synthetic_validation_results.csv\n\n")
