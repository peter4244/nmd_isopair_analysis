#!/usr/bin/env Rscript
#
# Self-Contained Test Runner for Reconstruction Pipeline
#
# Takes a GTF + pairs file, detects events, builds union exons, reconstructs,
# and verifies PASS/FAIL per pair. Sources the same production functions that
# Scripts 07-08 use.
#
# Usage:
#   Rscript scripts/tests/run_tests.R                          # default test data
#   Rscript scripts/tests/run_tests.R <gtf_file> <pairs_file>  # custom test data
#
# Pairs file format (tab-delimited):
#   gene_id  isoform_A  isoform_B  chr  strand
#   isoform_A = dominant, isoform_B = comparator
#
# Exit code:
#   0 = all pairs pass
#   1 = any failure or error
#

library(tidyverse)

# ==============================================================================
# Resolve paths — works from any working directory
# ==============================================================================

# Find project root (Version_6.0) relative to this script's location
script_dir <- if (interactive()) {

  "scripts/tests"
} else {
  # When run via Rscript, find the script path
  args_all <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args_all[grep("--file=", args_all)])
  if (length(script_path) > 0) {
    dirname(script_path)
  } else {
    "scripts/tests"
  }
}

# Navigate to project root (Version_6.0)
project_root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)

# Source production function libraries
source(file.path(project_root, "scripts/core/event_detection_functions.R"))
source(file.path(project_root, "scripts/core/reconstruction_functions.R"))

# Default test data location
default_gtf <- file.path(project_root, "scripts/tests/test_data/test_cases.gtf")
default_pairs <- file.path(project_root, "scripts/tests/test_data/pairs.tsv")

# ==============================================================================
# Parse GTF (same logic as dev harness detect_and_save_events.R)
# ==============================================================================

parse_gtf <- function(gtf_file) {
  lines <- readLines(gtf_file)
  lines <- lines[!startsWith(lines, '#') & lines != '']

  results <- list()
  for (line in lines) {
    parts <- strsplit(line, '\t')[[1]]
    if (parts[3] != "exon") next

    attrs_str <- parts[9]
    gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', attrs_str)
    transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', attrs_str)

    exon_number_str <- sub('.*exon_number "([^"]+)".*', '\\1', attrs_str)
    exon_number <- if (grepl('exon_number', attrs_str)) as.integer(exon_number_str) else NA_integer_

    results[[length(results) + 1]] <- tibble(
      chr = parts[1],
      exon_start = as.integer(parts[4]),
      exon_end = as.integer(parts[5]),
      strand = parts[7],
      gene_id = gene_id,
      transcript_id = transcript_id,
      exon_number = exon_number
    )
  }

  bind_rows(results)
}

# ==============================================================================
# Build Atomic Union Exons (same algorithm as Script 03 / dev harness)
# ==============================================================================

build_union_exons_from_gtf_data <- function(gtf_data) {
  genes <- unique(gtf_data$gene_id)
  all_union_exons <- list()
  ue_counter <- 1

  for (gene in genes) {
    gene_exons <- gtf_data %>% filter(gene_id == !!gene)
    chr <- unique(gene_exons$chr)[1]
    strand <- unique(gene_exons$strand)[1]

    all_boundaries <- sort(unique(c(gene_exons$exon_start, gene_exons$exon_end)))
    if (length(all_boundaries) < 2) next

    boundary_types <- sapply(all_boundaries, function(b) {
      is_start <- b %in% gene_exons$exon_start
      is_end <- b %in% gene_exons$exon_end
      if (is_start && is_end) return("BOTH")
      if (is_start) return("START")
      return("END")
    })

    n <- length(all_boundaries)
    seg_starts <- c()
    seg_ends <- c()
    current_start <- all_boundaries[1]

    for (i in seq_len(n - 2) + 1) {
      b <- all_boundaries[i]
      if (boundary_types[i] == "START") {
        seg_starts <- c(seg_starts, current_start)
        seg_ends <- c(seg_ends, b - 1)
        current_start <- b
      } else if (boundary_types[i] == "END") {
        seg_starts <- c(seg_starts, current_start)
        seg_ends <- c(seg_ends, b)
        current_start <- b + 1
      } else {
        seg_starts <- c(seg_starts, current_start)
        seg_ends <- c(seg_ends, b)
        current_start <- b + 1
      }
    }

    seg_starts <- c(seg_starts, current_start)
    seg_ends <- c(seg_ends, all_boundaries[n])

    segments <- tibble(seg_start = seg_starts, seg_end = seg_ends) %>%
      filter(seg_start <= seg_end)

    # Filter to segments covered by at least one exon
    covered <- sapply(seq_len(nrow(segments)), function(j) {
      any(gene_exons$exon_start <= segments$seg_start[j] &
          gene_exons$exon_end >= segments$seg_end[j])
    })
    covered_segments <- segments[covered, ]

    for (i in seq_len(nrow(covered_segments))) {
      all_union_exons[[ue_counter]] <- tibble(
        chr = chr,
        start = covered_segments$seg_start[i],
        end = covered_segments$seg_end[i],
        union_exon_id = sprintf("ue_%d", ue_counter),
        strand = strand,
        gene_id = gene
      )
      ue_counter <- ue_counter + 1
    }
  }

  bind_rows(all_union_exons)
}

# ==============================================================================
# Verify Transcript (same logic as Script 08)
# ==============================================================================

verify_transcript <- function(original_exons, reconstructed_exons, strand) {
  if (nrow(original_exons) != nrow(reconstructed_exons)) {
    return(list(
      pass = FALSE,
      reason = sprintf("Exon count mismatch: %d original vs %d reconstructed",
                       nrow(original_exons), nrow(reconstructed_exons)),
      n_exons_original = nrow(original_exons),
      n_exons_reconstructed = nrow(reconstructed_exons)
    ))
  }

  original_sorted <- original_exons %>% arrange(exon_start, exon_end)
  reconstructed_sorted <- reconstructed_exons %>% arrange(exon_start, exon_end)
  n_exons <- nrow(original_sorted)

  for (i in seq_len(n_exons)) {
    orig  <- original_sorted[i, ]
    recon <- reconstructed_sorted[i, ]

    is_first_genomic <- (i == 1)
    is_last_genomic  <- (i == n_exons)

    if (strand == "+") {
      start_tol <- if (is_first_genomic) TSS_TOLERANCE else 0L
      end_tol   <- if (is_last_genomic)  TES_TOLERANCE else 0L
    } else {
      start_tol <- if (is_first_genomic) TES_TOLERANCE else 0L
      end_tol   <- if (is_last_genomic)  TSS_TOLERANCE else 0L
    }

    start_diff <- abs(orig$exon_start - recon$exon_start)
    end_diff   <- abs(orig$exon_end   - recon$exon_end)

    if (start_diff > start_tol || end_diff > end_tol) {
      return(list(
        pass = FALSE,
        reason = sprintf("Exon %d mismatch: orig [%d-%d] vs recon [%d-%d]",
                         i, orig$exon_start, orig$exon_end,
                         recon$exon_start, recon$exon_end),
        n_exons_original = nrow(original_exons),
        n_exons_reconstructed = nrow(reconstructed_exons)
      ))
    }
  }

  return(list(
    pass = TRUE,
    reason = "Match",
    n_exons_original = nrow(original_exons),
    n_exons_reconstructed = nrow(reconstructed_exons)
  ))
}

# ==============================================================================
# Main Test Pipeline
# ==============================================================================

#' @param ue_rds_file Optional path to pre-built union exons RDS file.
#'   When provided, loads production UEs instead of building from GTF.
#'   Use this to exactly reproduce production pipeline behavior.
run_tests <- function(gtf_file, pairs_file, ue_rds_file = NULL) {

  cat("\n")
  cat("======================================================================\n")
  cat("  RECONSTRUCTION TEST RUNNER\n")
  cat("======================================================================\n")
  cat(sprintf("  GTF:   %s\n", gtf_file))
  cat(sprintf("  Pairs: %s\n", pairs_file))
  if (!is.null(ue_rds_file)) {
    cat(sprintf("  UEs:   %s (production)\n", ue_rds_file))
  }
  cat("======================================================================\n\n")

  # 1. Parse GTF
  cat("Parsing GTF...\n")
  gtf_data <- parse_gtf(gtf_file)
  cat(sprintf("  %d exon records, %d genes, %d transcripts\n",
              nrow(gtf_data), n_distinct(gtf_data$gene_id),
              n_distinct(gtf_data$transcript_id)))

  # 2. Build or load union exons
  #    When production UEs are provided, use them for matching genes and
  #    fall back to GTF-built UEs for genes not in production (e.g., synthetic tests).
  cat("\nBuilding atomic union exons from GTF...\n")
  gtf_union_exons <- build_union_exons_from_gtf_data(gtf_data)
  cat(sprintf("  %d union exons from %d genes (from GTF)\n",
              nrow(gtf_union_exons), n_distinct(gtf_union_exons$gene_id)))

  if (!is.null(ue_rds_file)) {
    cat("Loading production union exons...\n")
    ue_raw <- readRDS(ue_rds_file)
    # Adapt column names from RDS format to reconstruction format
    if ("union_exon_start" %in% names(ue_raw)) {
      prod_ue <- ue_raw %>%
        rename(start = union_exon_start, end = union_exon_end)
      if ("seqnames" %in% names(prod_ue)) {
        prod_ue <- prod_ue %>% rename(chr = seqnames)
      }
    } else {
      prod_ue <- ue_raw
    }
    # Merge: production UEs for genes that exist there, GTF UEs for the rest
    prod_genes <- unique(prod_ue$gene_id)
    gtf_only_genes <- setdiff(unique(gtf_union_exons$gene_id), prod_genes)
    union_exons <- bind_rows(
      prod_ue,
      gtf_union_exons %>% filter(gene_id %in% gtf_only_genes)
    )
    cat(sprintf("  Merged: %d production genes + %d GTF-only genes = %d total UEs\n",
                length(intersect(unique(gtf_data$gene_id), prod_genes)),
                length(gtf_only_genes),
                nrow(union_exons)))
  } else {
    union_exons <- gtf_union_exons
  }

  # 3. Read pairs
  pairs <- read_tsv(pairs_file, show_col_types = FALSE)
  cat(sprintf("\nProcessing %d pairs...\n", nrow(pairs)))

  # 4. Process each pair: detect -> reconstruct -> verify
  results <- list()
  pass_count <- 0L
  fail_count <- 0L
  error_count <- 0L

  for (i in seq_len(nrow(pairs))) {
    pair <- pairs[i, ]

    dom_exons <- gtf_data %>%
      filter(gene_id == pair$gene_id, transcript_id == pair$isoform_A)
    comp_exons <- gtf_data %>%
      filter(gene_id == pair$gene_id, transcript_id == pair$isoform_B)

    result <- tryCatch({

      if (nrow(dom_exons) == 0 || nrow(comp_exons) == 0) {
        list(status = "ERROR", reason = "Missing exons for one or both transcripts",
             events = character(0), n_exons_original = nrow(dom_exons),
             n_exons_reconstructed = NA_integer_)
      } else {

        gene_strand <- unique(dom_exons$strand)[1]

        # a. Detect events
        events <- detect_events_for_pair(
          dom_exons, comp_exons,
          pair$gene_id, pair$isoform_A, pair$isoform_B, gene_strand
        )

        # b. Handle zero-event case: no structural difference detected,
        #    so comparator IS the dominant (within tolerance).
        #    Verify directly without reconstruction.
        if (nrow(events) == 0) {
          vr <- verify_transcript(dom_exons, comp_exons, gene_strand)
          if (vr$pass) {
            list(status = "PASS", reason = "No events (identical within tolerance)",
                 events = "",
                 n_exons_original = vr$n_exons_original,
                 n_exons_reconstructed = vr$n_exons_reconstructed)
          } else {
            list(status = "FAIL", reason = paste0("No events but mismatch: ", vr$reason),
                 events = "",
                 n_exons_original = vr$n_exons_original,
                 n_exons_reconstructed = vr$n_exons_reconstructed)
          }
        } else {

          # c. Get gene union exons
          gene_ue <- union_exons %>% filter(gene_id == pair$gene_id)

          # d. Reconstruct dominant from comparator + events
          comp_for_recon <- comp_exons %>%
            select(chr, exon_start, exon_end, strand, gene_id, transcript_id)

          reconstructed <- reconstruct_dominant_v2(comp_for_recon, events, gene_ue)

          # e. Verify
          if (nrow(reconstructed) == 0) {
            list(status = "ERROR", reason = "Reconstruction produced 0 exons",
                 events = paste(events$event_type, collapse = ", "),
                 n_exons_original = nrow(dom_exons),
                 n_exons_reconstructed = 0L)
          } else {
            vr <- verify_transcript(dom_exons, reconstructed, gene_strand)
            if (vr$pass) {
              list(status = "PASS", reason = vr$reason,
                   events = paste(events$event_type, collapse = ", "),
                   n_exons_original = vr$n_exons_original,
                   n_exons_reconstructed = vr$n_exons_reconstructed)
            } else {
              list(status = "FAIL", reason = vr$reason,
                   events = paste(events$event_type, collapse = ", "),
                   n_exons_original = vr$n_exons_original,
                   n_exons_reconstructed = vr$n_exons_reconstructed)
            }
          }
        }
      }
    }, error = function(e) {
      list(status = "ERROR", reason = paste0("Exception: ", e$message),
           events = character(0), n_exons_original = NA_integer_,
           n_exons_reconstructed = NA_integer_)
    })

    if (result$status == "PASS") pass_count <- pass_count + 1L
    else if (result$status == "FAIL") fail_count <- fail_count + 1L
    else error_count <- error_count + 1L

    results[[i]] <- tibble(
      gene_id = pair$gene_id,
      dominant_isoform_id = pair$isoform_A,
      comparator_isoform_id = pair$isoform_B,
      status = result$status,
      reason = result$reason,
      events = if (length(result$events) == 0) "" else result$events,
      n_exons_original = result$n_exons_original,
      n_exons_reconstructed = result$n_exons_reconstructed
    )
  }

  results_df <- bind_rows(results)

  # 5. Print summary
  total <- nrow(results_df)
  cat("\n")
  cat("======================================================================\n")
  cat("  RESULTS\n")
  cat("======================================================================\n")
  cat(sprintf("  Total:  %d\n", total))
  cat(sprintf("  PASS:   %d/%d (%.1f%%)\n", pass_count, total, 100 * pass_count / total))
  cat(sprintf("  FAIL:   %d/%d (%.1f%%)\n", fail_count, total, 100 * fail_count / total))
  cat(sprintf("  ERROR:  %d/%d (%.1f%%)\n", error_count, total, 100 * error_count / total))
  cat("======================================================================\n")

  # 6. Print failure details
  if (fail_count > 0) {
    cat("\nFAILURES:\n")
    cat("----------------------------------------------------------------------\n")
    failures <- results_df %>% filter(status == "FAIL")
    for (j in seq_len(nrow(failures))) {
      f <- failures[j, ]
      cat(sprintf("  [FAIL] %s: %s vs %s\n", f$gene_id, f$dominant_isoform_id, f$comparator_isoform_id))
      cat(sprintf("         Events: %s\n", f$events))
      cat(sprintf("         Reason: %s\n", f$reason))
    }
  }

  if (error_count > 0) {
    cat("\nERRORS:\n")
    cat("----------------------------------------------------------------------\n")
    errors <- results_df %>% filter(status == "ERROR")
    for (j in seq_len(nrow(errors))) {
      e <- errors[j, ]
      cat(sprintf("  [ERROR] %s: %s vs %s\n", e$gene_id, e$dominant_isoform_id, e$comparator_isoform_id))
      cat(sprintf("          Reason: %s\n", e$reason))
    }
  }

  if (pass_count == total) {
    cat("\n======================================================================\n")
    cat("  ALL TESTS PASSED\n")
    cat("======================================================================\n")
  }

  cat("\n")

  return(results_df)
}

# ==============================================================================
# CLI
# ==============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  # Parse --ue flag
  ue_rds_file <- NULL
  if ("--ue" %in% args) {
    ue_idx <- which(args == "--ue")
    if (ue_idx < length(args)) {
      ue_rds_file <- args[ue_idx + 1]
      args <- args[-c(ue_idx, ue_idx + 1)]
    } else {
      cat("Error: --ue requires a path argument\n")
      quit(status = 1)
    }
  }

  if (length(args) == 0) {
    gtf_file <- default_gtf
    pairs_file <- default_pairs
  } else if (length(args) == 2) {
    gtf_file <- args[1]
    pairs_file <- args[2]
  } else {
    cat("Usage: Rscript scripts/tests/run_tests.R [<gtf_file> <pairs_file>] [--ue <union_exons.rds>]\n")
    cat("  No arguments: uses default test data in scripts/tests/test_data/\n")
    cat("  --ue: load production union exons from RDS instead of building from GTF\n")
    quit(status = 1)
  }

  if (!file.exists(gtf_file)) {
    cat(sprintf("Error: GTF file not found: %s\n", gtf_file))
    quit(status = 1)
  }
  if (!file.exists(pairs_file)) {
    cat(sprintf("Error: Pairs file not found: %s\n", pairs_file))
    quit(status = 1)
  }
  if (!is.null(ue_rds_file) && !file.exists(ue_rds_file)) {
    cat(sprintf("Error: Union exons RDS not found: %s\n", ue_rds_file))
    quit(status = 1)
  }

  results <- run_tests(gtf_file, pairs_file, ue_rds_file)

  # Exit code: 0 if all pass, 1 otherwise
  if (all(results$status == "PASS")) {
    quit(status = 0)
  } else {
    quit(status = 1)
  }
}
