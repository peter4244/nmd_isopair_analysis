#!/usr/bin/env Rscript
#
# Detect Events and Save with Transcript IDs
#
# This script detects all splicing events between isoform pairs and saves them
# in a format suitable for reconstruction validation.
#
# Inputs:
#   - GTF file with isoform pairs
#   - Test annotations with expected comparisons
#
# Outputs:
#   - Events TSV file with all required columns for reconstruction

library(tidyverse)

# Source event detection functions
# Path relative to testing/ directory when run from there
source("../scripts/event_detection_functions.R")

# ==============================================================================
# Parse GTF
# ==============================================================================

#' Parse GTF file to extract exon information
#'
#' @param gtf_file Path to GTF file
#' @return Tibble with exon information
parse_gtf <- function(gtf_file) {
  cat("Reading GTF file:", gtf_file, "\n")

  lines <- readLines(gtf_file)
  lines <- lines[!startsWith(lines, '#') & lines != '']

  results <- list()
  for (line in lines) {
    parts <- strsplit(line, '\t')[[1]]

    # Only parse exon features
    if (parts[3] != "exon") {
      next
    }

    attrs_str <- parts[9]
    gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', attrs_str)
    transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', attrs_str)

    # Extract exon_number if present
    exon_number_str <- sub('.*exon_number "([^"]+)".*', '\\1', attrs_str)
    exon_number <- if (grepl('exon_number', attrs_str)) as.integer(exon_number_str) else NA_integer_

    results[[length(results) + 1]] <- tibble(
      chr = parts[1],
      exon_start = as.integer(parts[4]),  # 1-based
      exon_end = as.integer(parts[5]),    # 1-based
      strand = parts[7],
      gene_id = gene_id,
      transcript_id = transcript_id,
      exon_number = exon_number
    )
  }

  gtf_data <- bind_rows(results)

  cat(sprintf("  Parsed %d exon records\n", nrow(gtf_data)))
  cat(sprintf("  %d genes, %d transcripts\n",
              n_distinct(gtf_data$gene_id),
              n_distinct(gtf_data$transcript_id)))

  return(gtf_data)
}

# ==============================================================================
# Helper Functions
# ==============================================================================

#' Check if two exons have terminal overlap (≥50%)
#'
#' @param exon1 First exon
#' @param exon2 Second exon
#' @return Logical
has_terminal_overlap <- function(exon1, exon2) {
  # Calculate overlap
  overlap_start <- max(exon1$exon_start, exon2$exon_start)
  overlap_end <- min(exon1$exon_end, exon2$exon_end)

  if (overlap_start > overlap_end) return(FALSE)

  overlap_length <- overlap_end - overlap_start + 1
  len1 <- exon1$exon_end - exon1$exon_start + 1
  len2 <- exon2$exon_end - exon2$exon_start + 1

  overlap_pct <- overlap_length / min(len1, len2)

  return(overlap_pct >= 0.5)
}

# ==============================================================================
# Event Detection for Transcript Pair
# ==============================================================================

#' Detect all events between two transcripts
#'
#' @param dominant_exons Exons from dominant transcript
#' @param comparator_exons Exons from comparator transcript
#' @param gene_id Gene identifier
#' @param dominant_tid Dominant transcript ID
#' @param comparator_tid Comparator transcript ID
#' @param strand Gene strand
#' @return Tibble with detected events
detect_events_for_pair <- function(dominant_exons, comparator_exons,
                                   gene_id, dominant_tid, comparator_tid, strand) {

  # Order exons biologically (TSS → TES)
  dom_ordered <- order_exons_biological(dominant_exons, strand)
  comp_ordered <- order_exons_biological(comparator_exons, strand)

  # Initialize events list
  all_events <- list()

  # Check terminal boundaries (TSS and TES)
  first_dom <- dom_ordered[1, ]
  first_comp <- comp_ordered[1, ]
  last_dom <- dom_ordered[nrow(dom_ordered), ]
  last_comp <- comp_ordered[nrow(comp_ordered), ]

  # Detect Alt_TSS
  if (detect_tss_change(first_dom, first_comp, strand)) {
    # Determine direction and missing exons
    if (strand == "+") {
      dom_tss <- first_dom$exon_start
      comp_tss <- first_comp$exon_start
      five_prime <- min(dom_tss, comp_tss)
      three_prime <- max(dom_tss, comp_tss)
    } else {
      dom_tss <- first_dom$exon_end
      comp_tss <- first_comp$exon_end
      five_prime <- max(dom_tss, comp_tss)
      three_prime <- min(dom_tss, comp_tss)
    }

    # Direction: LOSS if dominant extends beyond comparator
    direction <- if (
      (strand == "+" && dom_tss < comp_tss) ||
      (strand == "-" && dom_tss > comp_tss)
    ) "LOSS" else "GAIN"

    # Missing terminal exons (only for LOSS)
    missing_exons <- ""
    if (direction == "LOSS") {
      # Use core function to properly walk exons and record only exonic regions
      missing_exons <- compute_missing_terminal_exons_tss(dom_ordered, comp_ordered, strand)
    }

    all_events[[length(all_events) + 1]] <- tibble(
      gene_id = gene_id,
      dominant_transcript_id = dominant_tid,
      comparator_transcript_id = comparator_tid,
      event_type = "Alt_TSS",
      direction = direction,
      chr = first_dom$chr,
      five_prime = five_prime,
      three_prime = three_prime,
      strand = strand,
      bp_diff = abs(dom_tss - comp_tss),
      missing_terminal_exons = missing_exons
    )
  }

  # Detect Alt_TES
  if (detect_tes_change(last_dom, last_comp, strand)) {
    # Similar logic for TES
    if (strand == "+") {
      dom_tes <- last_dom$exon_end
      comp_tes <- last_comp$exon_end
      five_prime <- min(dom_tes, comp_tes)
      three_prime <- max(dom_tes, comp_tes)
    } else {
      dom_tes <- last_dom$exon_start
      comp_tes <- last_comp$exon_start
      five_prime <- max(dom_tes, comp_tes)
      three_prime <- min(dom_tes, comp_tes)
    }

    direction <- if (
      (strand == "+" && dom_tes > comp_tes) ||
      (strand == "-" && dom_tes < comp_tes)
    ) "LOSS" else "GAIN"

    missing_exons <- ""
    if (direction == "LOSS") {
      # Use core function to properly walk exons and record only exonic regions
      missing_exons <- compute_missing_terminal_exons_tes(dom_ordered, comp_ordered, strand)
    }

    all_events[[length(all_events) + 1]] <- tibble(
      gene_id = gene_id,
      dominant_transcript_id = dominant_tid,
      comparator_transcript_id = comparator_tid,
      event_type = "Alt_TES",
      direction = direction,
      chr = last_dom$chr,
      five_prime = five_prime,
      three_prime = three_prime,
      strand = strand,
      bp_diff = abs(dom_tes - comp_tes),
      missing_terminal_exons = missing_exons
    )
  }

  # Detect internal events (A5SS, A3SS, Partial_IR, etc.)
  # Compare all overlapping exon pairs
  for (i in seq_len(nrow(comp_ordered))) {
    comp_exon <- comp_ordered[i, ]
    is_first_comp <- (i == 1)
    is_last_comp <- (i == nrow(comp_ordered))

    for (j in seq_len(nrow(dom_ordered))) {
      dom_exon <- dom_ordered[j, ]
      is_first_dom <- (j == 1)
      is_last_dom <- (j == nrow(dom_ordered))

      # Check if exons overlap
      overlaps <- (dom_exon$exon_start <= comp_exon$exon_end) &&
                  (dom_exon$exon_end >= comp_exon$exon_start)

      if (!overlaps) next

      # Check for terminal overlap
      terminal_overlap <- has_terminal_overlap(dom_exon, comp_exon)

      # Use unified event detection
      is_both_first <- is_first_dom && is_first_comp
      is_both_last <- is_last_dom && is_last_comp

      event_result <- detect_shared_boundary_event(
        dom_exon, comp_exon, strand,
        is_first_exon = is_both_first,
        is_last_exon = is_both_last,
        is_first_exon_comp = is_first_comp,
        is_last_exon_comp = is_last_comp,
        terminal_has_overlap = terminal_overlap,
        flanking_exons_dom = NULL,  # Simplified for now
        flanking_exons_non_dom = NULL
      )

      if (event_result$event_type != "none") {
        # Compute precise event coordinates: the actual exonic region that differs
        # five_prime and three_prime represent the gained/lost region

        event_type <- event_result$event_type
        direction <- if (!is.null(event_result$direction)) event_result$direction else "-"

        # Compute coordinates based on event type and direction
        if (event_type %in% c("A3SS", "Partial_IR_3")) {
          # Acceptor differs (3' splice site)
          if (strand == "+") {
            # Plus: acceptor = exon start
            if (direction == "LOSS") {
              # Comparator lost sequence: region from dom_start to comp_start-1
              five_prime <- dom_exon$exon_start
              three_prime <- comp_exon$exon_start - 1
            } else {
              # Comparator gained sequence: region from comp_start to dom_start-1
              five_prime <- comp_exon$exon_start
              three_prime <- dom_exon$exon_start - 1
            }
          } else {
            # Minus: acceptor = exon end
            if (direction == "LOSS") {
              # Comparator lost sequence: region from comp_end+1 to dom_end
              five_prime <- dom_exon$exon_end
              three_prime <- comp_exon$exon_end + 1
            } else {
              # Comparator gained sequence: region from dom_end+1 to comp_end
              five_prime <- comp_exon$exon_end
              three_prime <- dom_exon$exon_end + 1
            }
          }
        } else if (event_type %in% c("A5SS", "Partial_IR_5")) {
          # Donor differs (5' splice site)
          if (strand == "+") {
            # Plus: donor = exon end
            if (direction == "LOSS") {
              # Comparator lost sequence: region from comp_end+1 to dom_end
              five_prime <- comp_exon$exon_end + 1
              three_prime <- dom_exon$exon_end
            } else {
              # Comparator gained sequence: region from dom_end+1 to comp_end
              five_prime <- dom_exon$exon_end + 1
              three_prime <- comp_exon$exon_end
            }
          } else {
            # Minus: donor = exon start
            if (direction == "LOSS") {
              # Comparator lost sequence: region from dom_start to comp_start-1
              five_prime <- comp_exon$exon_start
              three_prime <- dom_exon$exon_start
            } else {
              # Comparator gained sequence: region from comp_start to dom_start-1
              five_prime <- dom_exon$exon_start
              three_prime <- comp_exon$exon_start
            }
          }
        } else {
          # For other event types (Dual_boundary, etc.), use union as before
          if (strand == "+") {
            five_prime <- min(dom_exon$exon_start, comp_exon$exon_start)
            three_prime <- max(dom_exon$exon_end, comp_exon$exon_end)
          } else {
            five_prime <- max(dom_exon$exon_end, comp_exon$exon_end)
            three_prime <- min(dom_exon$exon_start, comp_exon$exon_start)
          }
        }

        all_events[[length(all_events) + 1]] <- tibble(
          gene_id = gene_id,
          dominant_transcript_id = dominant_tid,
          comparator_transcript_id = comparator_tid,
          event_type = event_type,
          direction = direction,
          chr = dom_exon$chr,
          five_prime = five_prime,
          three_prime = three_prime,
          strand = strand,
          bp_diff = event_result$bp_diff,
          missing_terminal_exons = ""
        )
      }
    }
  }

  # Detect IR (monoexonic spanning multi-exonic)
  # Check if comparator exons span multiple dominant exons
  for (i in seq_len(nrow(comp_ordered))) {
    if (detect_ir_simple(comp_ordered[i, ], dom_ordered)) {
      comp_exon <- comp_ordered[i, ]
      all_events[[length(all_events) + 1]] <- tibble(
        gene_id = gene_id,
        dominant_transcript_id = dominant_tid,
        comparator_transcript_id = comparator_tid,
        event_type = "IR",
        direction = "GAIN",  # Comparator has retention
        chr = comp_exon$chr,
        five_prime = if (strand == "+") comp_exon$exon_start else comp_exon$exon_end,
        three_prime = if (strand == "+") comp_exon$exon_end else comp_exon$exon_start,
        strand = strand,
        bp_diff = NA_integer_,
        missing_terminal_exons = ""
      )
    }
  }

  # Check if dominant exons span multiple comparator exons
  for (i in seq_len(nrow(dom_ordered))) {
    if (detect_ir_simple(dom_ordered[i, ], comp_ordered)) {
      dom_exon <- dom_ordered[i, ]
      all_events[[length(all_events) + 1]] <- tibble(
        gene_id = gene_id,
        dominant_transcript_id = dominant_tid,
        comparator_transcript_id = comparator_tid,
        event_type = "IR",
        direction = "LOSS",  # Dominant has retention
        chr = dom_exon$chr,
        five_prime = if (strand == "+") dom_exon$exon_start else dom_exon$exon_end,
        three_prime = if (strand == "+") dom_exon$exon_end else dom_exon$exon_start,
        strand = strand,
        bp_diff = NA_integer_,
        missing_terminal_exons = ""
      )
    }
  }

  # Combine all events
  if (length(all_events) > 0) {
    return(bind_rows(all_events))
  } else {
    return(tibble())
  }
}

# ==============================================================================
# Main Pipeline
# ==============================================================================

main <- function(gtf_file, test_annotations_file, output_file) {
  cat("\n╔════════════════════════════════════════════════════════════════╗\n")
  cat("║   Event Detection and Export Pipeline                         ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")

  # Parse GTF
  gtf_data <- parse_gtf(gtf_file)

  # Load test annotations to get transcript pairs
  test_annot <- read_tsv(test_annotations_file, show_col_types = FALSE)

  cat(sprintf("\nProcessing %d test cases...\n", nrow(test_annot)))

  all_events <- list()

  for (i in seq_len(nrow(test_annot))) {
    test <- test_annot[i, ]

    cat(sprintf("\n[%d/%d] %s\n", i, nrow(test_annot), test$gene_id))

    # Extract exons for both transcripts
    dom_exons <- gtf_data %>%
      filter(gene_id == test$gene_id, transcript_id == test$isoform_A)

    comp_exons <- gtf_data %>%
      filter(gene_id == test$gene_id, transcript_id == test$isoform_B)

    if (nrow(dom_exons) == 0 || nrow(comp_exons) == 0) {
      cat("  Warning: Missing exons for one or both transcripts\n")
      next
    }

    # Get strand from GTF data
    gene_strand <- unique(dom_exons$strand)[1]

    # Detect events
    events <- detect_events_for_pair(
      dom_exons, comp_exons,
      test$gene_id, test$isoform_A, test$isoform_B, gene_strand
    )

    if (nrow(events) > 0) {
      cat(sprintf("  Detected %d events: %s\n",
                  nrow(events),
                  paste(unique(events$event_type), collapse = ", ")))
      all_events[[length(all_events) + 1]] <- events
    } else {
      cat("  No events detected\n")
    }
  }

  # Combine all events
  if (length(all_events) > 0) {
    final_events <- bind_rows(all_events)

    # Write events file
    cat(sprintf("\n\nWriting events to: %s\n", output_file))
    write_tsv(final_events, output_file)

    cat(sprintf("  Total events: %d\n", nrow(final_events)))
    cat(sprintf("  Event types: %s\n",
                paste(names(table(final_events$event_type)), collapse = ", ")))
  } else {
    cat("\nNo events detected!\n")
  }

  cat("\n═══════════════════════════════════════════════════════════════════\n")
  cat("EVENT DETECTION COMPLETE\n")
  cat("═══════════════════════════════════════════════════════════════════\n")
}

# ==============================================================================
# Command Line Interface
# ==============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 3) {
    cat("Usage: detect_and_save_events.R <gtf_file> <test_annotations> <output_events_file>\n")
    quit(status = 1)
  }

  main(args[1], args[2], args[3])
} else {
  # Interactive testing
  main(
    gtf_file = "../synthetic/TestData/exons/base_events.gtf",
    test_annotations_file = "../synthetic/TestData/annotations/base_events.tsv",
    output_file = "../synthetic/base_events_events.tsv"
  )
}
