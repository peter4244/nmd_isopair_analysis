#!/usr/bin/env Rscript
#
# Verify Reconstruction
#
# This script verifies that reconstructed dominant isoforms exactly match
# the original dominant isoforms.
#
# Inputs:
#   - Original GTF (all isoforms)
#   - Reconstructed GTF (reconstructed dominants)
#   - Events file (to identify dominant transcript IDs)
#
# Outputs:
#   - Verification report (PASS/FAIL per transcript)
#   - Summary statistics

library(tidyverse)

# Tolerance for terminal exon boundaries — must match event_detection_functions.R
TSS_TOLERANCE <- 20
TES_TOLERANCE <- 20

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
      exon_start = as.integer(parts[4]),
      exon_end = as.integer(parts[5]),
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
# Verification Functions
# ==============================================================================

#' Verify that two exon structures match within tolerance
#'
#' Internal exon boundaries are checked for exact match. Terminal exon boundaries
#' (TSS and TES) allow up to TSS_TOLERANCE / TES_TOLERANCE bp difference, matching
#' the tolerance used in event detection.
#'
#' @param original_exons Original dominant exons
#' @param reconstructed_exons Reconstructed dominant exons
#' @param strand Gene strand ("+" or "-")
#' @param gene_id Gene ID (for reporting)
#' @param original_tid Original transcript ID
#' @param reconstructed_tid Reconstructed transcript ID
#' @return List with pass status and details
verify_transcript <- function(original_exons, reconstructed_exons, strand,
                              gene_id, original_tid, reconstructed_tid) {

  # Check exon count
  if (nrow(original_exons) != nrow(reconstructed_exons)) {
    return(list(
      pass = FALSE,
      reason = sprintf("Exon count mismatch: %d original vs %d reconstructed",
                      nrow(original_exons), nrow(reconstructed_exons)),
      n_exons_original = nrow(original_exons),
      n_exons_reconstructed = nrow(reconstructed_exons),
      first_mismatch_pos = NA,
      first_mismatch_detail = ""
    ))
  }

  # Sort both by genomic position
  original_sorted <- original_exons %>% arrange(exon_start, exon_end)
  reconstructed_sorted <- reconstructed_exons %>% arrange(exon_start, exon_end)

  n_exons <- nrow(original_sorted)

  # Compare each exon coordinate with strand-aware terminal tolerance
  for (i in seq_len(n_exons)) {
    orig  <- original_sorted[i, ]
    recon <- reconstructed_sorted[i, ]

    is_first_genomic <- (i == 1)
    is_last_genomic  <- (i == n_exons)

    # Assign tolerance per boundary based on strand
    # + strand: TSS = exon_start of first exon, TES = exon_end of last exon
    # - strand: TES = exon_start of first exon, TSS = exon_end of last exon
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
        reason = sprintf("Exon %d coordinate mismatch", i),
        n_exons_original = nrow(original_exons),
        n_exons_reconstructed = nrow(reconstructed_exons),
        first_mismatch_pos = i,
        first_mismatch_detail = sprintf(
          "Original: [%d-%d], Reconstructed: [%d-%d]",
          orig$exon_start, orig$exon_end,
          recon$exon_start, recon$exon_end
        )
      ))
    }
  }

  # All checks passed
  tol_note <- if (n_exons > 0 &&
    any(abs(original_sorted$exon_start - reconstructed_sorted$exon_start) > 0 |
        abs(original_sorted$exon_end   - reconstructed_sorted$exon_end)   > 0))
    " (within terminal tolerance)" else ""

  return(list(
    pass = TRUE,
    reason = paste0("Perfect match", tol_note),
    n_exons_original = nrow(original_exons),
    n_exons_reconstructed = nrow(reconstructed_exons),
    first_mismatch_pos = NA,
    first_mismatch_detail = ""
  ))
}

# ==============================================================================
# Main Verification
# ==============================================================================

#' Verify all reconstructed transcripts
#'
#' @param original_gtf Original GTF data
#' @param reconstructed_gtf Reconstructed GTF data
#' @param events Events data (to get transcript pairs)
#' @return Verification results
verify_all <- function(original_gtf, reconstructed_gtf, events) {

  # Get unique transcript pairs from events
  transcript_pairs <- events %>%
    select(gene_id, dominant_transcript_id, comparator_transcript_id) %>%
    distinct()

  cat(sprintf("\nVerifying %d transcript pairs...\n", nrow(transcript_pairs)))

  verification_results <- list()

  for (i in seq_len(nrow(transcript_pairs))) {
    pair <- transcript_pairs[i, ]

    cat(sprintf("\n[%d/%d] %s: %s\n",
                i, nrow(transcript_pairs),
                pair$gene_id,
                pair$dominant_transcript_id))

    # Extract original dominant exons
    original_exons <- original_gtf %>%
      filter(
        gene_id == pair$gene_id,
        transcript_id == pair$dominant_transcript_id
      ) %>%
      arrange(exon_start)

    # Extract reconstructed exons
    # Reconstructed transcript IDs have "_reconstructed" suffix
    reconstructed_tid <- paste0(pair$dominant_transcript_id, "_reconstructed")

    reconstructed_exons <- reconstructed_gtf %>%
      filter(
        gene_id == pair$gene_id,
        transcript_id == reconstructed_tid
      ) %>%
      arrange(exon_start)

    if (nrow(original_exons) == 0) {
      cat("  ERROR: Original dominant not found\n")
      verification_results[[i]] <- tibble(
        gene_id = pair$gene_id,
        dominant_transcript_id = pair$dominant_transcript_id,
        comparator_transcript_id = pair$comparator_transcript_id,
        status = "ERROR",
        reason = "Original dominant not found",
        n_exons_original = 0,
        n_exons_reconstructed = nrow(reconstructed_exons),
        first_mismatch_pos = NA_integer_,
        first_mismatch_detail = ""
      )
      next
    }

    if (nrow(reconstructed_exons) == 0) {
      cat("  ERROR: Reconstructed transcript not found\n")
      verification_results[[i]] <- tibble(
        gene_id = pair$gene_id,
        dominant_transcript_id = pair$dominant_transcript_id,
        comparator_transcript_id = pair$comparator_transcript_id,
        status = "ERROR",
        reason = "Reconstructed transcript not found",
        n_exons_original = nrow(original_exons),
        n_exons_reconstructed = 0,
        first_mismatch_pos = NA_integer_,
        first_mismatch_detail = ""
      )
      next
    }

    # Verify
    result <- verify_transcript(
      original_exons, reconstructed_exons,
      strand = original_exons$strand[1],
      pair$gene_id, pair$dominant_transcript_id, reconstructed_tid
    )

    status <- if (result$pass) "PASS" else "FAIL"
    cat(sprintf("  %s: %s\n", status, result$reason))

    if (!result$pass) {
      cat(sprintf("    %s\n", result$first_mismatch_detail))
    }

    verification_results[[i]] <- tibble(
      gene_id = pair$gene_id,
      dominant_transcript_id = pair$dominant_transcript_id,
      comparator_transcript_id = pair$comparator_transcript_id,
      status = status,
      reason = result$reason,
      n_exons_original = result$n_exons_original,
      n_exons_reconstructed = result$n_exons_reconstructed,
      first_mismatch_pos = result$first_mismatch_pos,
      first_mismatch_detail = result$first_mismatch_detail
    )
  }

  return(bind_rows(verification_results))
}

# ==============================================================================
# Main Pipeline
# ==============================================================================

main <- function(original_gtf_file, reconstructed_gtf_file, events_file, output_report) {
  cat("\n╔════════════════════════════════════════════════════════════════╗\n")
  cat("║   Reconstruction Verification Pipeline                        ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")

  # Load inputs
  original_gtf <- parse_gtf(original_gtf_file)
  reconstructed_gtf <- parse_gtf(reconstructed_gtf_file)

  cat("\nReading events file:", events_file, "\n")
  events <- read_tsv(events_file, show_col_types = FALSE)
  cat(sprintf("  %d events loaded\n", nrow(events)))

  # Verify
  results <- verify_all(original_gtf, reconstructed_gtf, events)

  # Write report
  cat("\n\nWriting verification report:", output_report, "\n")
  write_tsv(results, output_report)

  # Summary statistics
  cat("\n═══════════════════════════════════════════════════════════════════\n")
  cat("VERIFICATION SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════════\n\n")

  pass_count <- sum(results$status == "PASS")
  fail_count <- sum(results$status == "FAIL")
  error_count <- sum(results$status == "ERROR")
  total_count <- nrow(results)

  cat(sprintf("Total transcripts: %d\n", total_count))
  cat(sprintf("✓ PASS:  %d/%d (%.1f%%)\n",
              pass_count, total_count,
              100 * pass_count / total_count))
  cat(sprintf("✗ FAIL:  %d/%d (%.1f%%)\n",
              fail_count, total_count,
              100 * fail_count / total_count))
  cat(sprintf("⚠ ERROR: %d/%d (%.1f%%)\n\n",
              error_count, total_count,
              100 * error_count / total_count))

  if (pass_count == total_count) {
    cat("═══════════════════════════════════════════════════════════════════\n")
    cat("✓✓✓ PERFECT RECONSTRUCTION - All transcripts match exactly! ✓✓✓\n")
    cat("═══════════════════════════════════════════════════════════════════\n")
  } else {
    cat("Failures detected:\n\n")
    failures <- results %>% filter(status != "PASS")
    for (i in seq_len(nrow(failures))) {
      f <- failures[i, ]
      cat(sprintf("  [%s] %s: %s\n",
                  f$status,
                  f$dominant_transcript_id,
                  f$reason))
      if (f$first_mismatch_detail != "") {
        cat(sprintf("    %s\n", f$first_mismatch_detail))
      }
    }
  }

  cat(sprintf("\nReport saved to: %s\n", output_report))

  # Return exit code based on results
  if (pass_count == total_count) {
    quit(status = 0)
  } else {
    quit(status = 1)
  }
}

# ==============================================================================
# Command Line Interface
# ==============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 4) {
    cat("Usage: verify_reconstruction.R <original_gtf> <reconstructed_gtf> <events_file> <output_report>\n")
    quit(status = 1)
  }

  main(args[1], args[2], args[3], args[4])
} else {
  # Interactive testing
  main(
    original_gtf_file = "../synthetic/TestData/exons/base_events.gtf",
    reconstructed_gtf_file = "../synthetic/base_events_reconstructed.gtf",
    events_file = "../synthetic/base_events_events.tsv",
    output_report = "../synthetic/verification_report.tsv"
  )
}
