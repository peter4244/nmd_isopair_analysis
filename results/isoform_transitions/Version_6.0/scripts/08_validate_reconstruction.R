#!/usr/bin/env Rscript

################################################################################
# Script 08: Validate Reconstruction
################################################################################
#
# Purpose:
#   Validate the event detection + reconstruction pipeline by reconstructing
#   dominant isoforms from comparator exons + detected events, then comparing
#   against original dominant exon structures. This runs entirely in-memory
#   using RDS data (no GTF files).
#
# Input:
#   - data/splicing_choice_profiles.rds  (nested detailed_events per profile)
#   - data/isoform_structures_filtered.rds  (compact exon structures)
#   - data/union_exons_filtered.rds  (union exon coordinates)
#
# Output:
#   - data/reconstruction_verification.rds  (per-pair PASS/FAIL/ERROR results)
#   - logs/reconstruction_verification.log  (summary statistics)
#
################################################################################

library(tidyverse)

# ==============================================================================
# Input Validation Helpers
# ==============================================================================

validate_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required input file not found: %s\nHave you run the prerequisite scripts?", path))
  }
}

validate_columns <- function(df, required_cols, name) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing required columns: %s", name, paste(missing, collapse = ", ")))
  }
}

validate_overlap <- function(ids_a, ids_b, name_a, name_b, min_pct = 10) {
  overlap <- length(intersect(ids_a, ids_b))
  pct <- 100 * overlap / length(ids_a)
  cat(sprintf("  Cross-check: %d/%d %s found in %s (%.1f%%)\n",
              overlap, length(ids_a), name_a, name_b, pct))
  if (pct < min_pct) {
    stop(sprintf("FATAL: Only %.1f%% of %s found in %s. Files may be from different pipeline runs.",
                 pct, name_a, name_b))
  }
}

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 8: Validate Reconstruction                            ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Paths
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"
setwd(base_dir)

# ==============================================================================
# 0. Validate Inputs and Source Reconstruction Functions
# ==============================================================================

cat("Validating inputs...\n")
validate_file("data/splicing_choice_profiles.rds")
validate_file("data/isoform_structures_filtered.rds")
validate_file("data/union_exons_filtered.rds")
validate_file("scripts/reconstruction_functions.R")
cat("  All required input files found.\n\n")

cat("Loading reconstruction functions...\n")
source("scripts/reconstruction_functions.R")
cat("  reconstruct_dominant_v2() loaded\n\n")

# Verification tolerance — must match event_detection_functions.R
TSS_TOLERANCE <- 20
TES_TOLERANCE <- 20

# ==============================================================================
# 1. Load Data
# ==============================================================================

cat("Loading data...\n")

profiles <- readRDS("data/splicing_choice_profiles.rds")
isoform_structures_compact <- readRDS("data/isoform_structures_filtered.rds")
union_exons_raw <- readRDS("data/union_exons_filtered.rds")

cat(sprintf("  Profiles: %d\n", nrow(profiles)))
cat(sprintf("  Isoform structures (compact): %d\n", nrow(isoform_structures_compact)))
cat(sprintf("  Union exons: %d\n", nrow(union_exons_raw)))

# Column validation
validate_columns(profiles,
                 c("gene_id", "dominant_isoform_id", "non_dominant_isoform_id",
                   "n_events", "detailed_events"),
                 "splicing_choice_profiles")
validate_columns(isoform_structures_compact,
                 c("isoform_id", "gene_id", "seqnames", "strand", "exon_starts", "exon_ends"),
                 "isoform_structures_filtered")
validate_columns(union_exons_raw,
                 c("gene_id", "union_exon_id", "union_exon_start", "union_exon_end", "seqnames", "strand"),
                 "union_exons_filtered")

# Cross-file consistency checks
cat("\nCross-file consistency checks:\n")
validate_overlap(unique(profiles$gene_id),
                 unique(union_exons_raw$gene_id),
                 "profiles gene_ids", "union_exons")
profile_isoform_ids <- unique(c(profiles$dominant_isoform_id, profiles$non_dominant_isoform_id))
validate_overlap(profile_isoform_ids,
                 unique(isoform_structures_compact$isoform_id),
                 "profile isoform_ids", "isoform_structures")
cat("")

# ==============================================================================
# 2. Prepare Data
# ==============================================================================

cat("\nPreparing data for reconstruction...\n")

# --- 2a. Extract unique pairs and their events ---
cat("  Extracting events from profiles...\n")

# Each profile has a nested detailed_events tibble
# Extract all events into a single flat tibble
all_events <- profiles %>%
  filter(n_events > 0) %>%
  select(dominant_isoform_id, non_dominant_isoform_id, detailed_events) %>%
  unnest(detailed_events)

cat(sprintf("  Total events extracted: %d\n", nrow(all_events)))

# Get unique pairs (including those with zero events)
unique_pairs <- profiles %>%
  select(gene_id, dominant_isoform_id, non_dominant_isoform_id) %>%
  distinct()

cat(sprintf("  Unique pairs to validate: %d\n", nrow(unique_pairs)))

# --- 2b. Expand isoform structures to one-row-per-exon ---
cat("  Expanding isoform structures...\n")

isoform_structures <- isoform_structures_compact %>%
  rowwise() %>%
  mutate(
    exon_data = list(tibble(
      exon_number = seq_along(exon_starts),
      exon_start = exon_starts,
      exon_end = exon_ends
    ))
  ) %>%
  unnest(exon_data) %>%
  select(isoform_id, gene_id, seqnames, strand, exon_number, exon_start, exon_end) %>%
  ungroup()

cat(sprintf("  Expanded to %d exon rows\n", nrow(isoform_structures)))

# --- 2c. Rename union exons for reconstruction compatibility ---
# RDS format: union_exon_id, union_exon_start, union_exon_end, gene_id, seqnames, strand
# Reconstruction expects: union_exon_id, start, end, gene_id, chr, strand
cat("  Adapting union exon column names...\n")

union_exons <- union_exons_raw %>%
  rename(
    start = union_exon_start,
    end = union_exon_end,
    chr = seqnames
  )

cat(sprintf("  Union exons ready: %d rows\n", nrow(union_exons)))

# ==============================================================================
# 3. Verification Function (inlined from verify_reconstruction.R)
# ==============================================================================

#' Verify that two exon structures match within tolerance
#'
#' Internal exon boundaries require exact match. Terminal exon boundaries
#' (TSS and TES) allow up to TSS_TOLERANCE / TES_TOLERANCE bp difference.
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
# 4. Reconstruct + Verify (per pair)
# ==============================================================================

cat("\nReconstructing and verifying...\n")

batch_size <- 500
n_pairs <- nrow(unique_pairs)
n_batches <- ceiling(n_pairs / batch_size)

verification_results <- list()
pass_count <- 0L
fail_count <- 0L
error_count <- 0L

# Suppress reconstruction debug messages
suppressWarnings({

for (batch_idx in 1:n_batches) {
  start_idx <- (batch_idx - 1) * batch_size + 1
  end_idx <- min(batch_idx * batch_size, n_pairs)
  batch_pairs <- unique_pairs[start_idx:end_idx, ]

  for (row_i in seq_len(nrow(batch_pairs))) {
    pair <- batch_pairs[row_i, ]
    dom_id <- pair$dominant_isoform_id
    comp_id <- pair$non_dominant_isoform_id
    gene <- pair$gene_id

    result <- tryCatch({

      # --- a. Extract comparator exons ---
      comp_exons <- isoform_structures %>%
        filter(isoform_id == comp_id) %>%
        rename(chr = seqnames, transcript_id = isoform_id) %>%
        select(chr, exon_start, exon_end, strand, gene_id, transcript_id)

      if (nrow(comp_exons) == 0) {
        list(status = "ERROR", reason = "No comparator exons",
             n_exons_original = NA_integer_, n_exons_reconstructed = NA_integer_)
      } else {

        # --- b. Extract dominant exons (for verification) ---
        dom_exons <- isoform_structures %>%
          filter(isoform_id == dom_id) %>%
          arrange(exon_start)

        if (nrow(dom_exons) == 0) {
          list(status = "ERROR", reason = "No dominant exons",
               n_exons_original = NA_integer_, n_exons_reconstructed = NA_integer_)
        } else {

          strand <- dom_exons$strand[1]

          # --- c. Filter events for this pair ---
          pair_events <- all_events %>%
            filter(dominant_transcript_id == dom_id,
                   comparator_transcript_id == comp_id)

          # --- d. Filter union exons for this gene ---
          gene_union_exons <- union_exons %>%
            filter(gene_id == gene)

          # --- e. Reconstruct ---
          reconstructed <- reconstruct_dominant_v2(
            comp_exons, pair_events, gene_union_exons
          )

          # --- f. Verify ---
          if (nrow(reconstructed) == 0) {
            list(status = "ERROR", reason = "Reconstruction produced 0 exons",
                 n_exons_original = nrow(dom_exons), n_exons_reconstructed = 0L)
          } else {
            vr <- verify_transcript(dom_exons, reconstructed, strand)
            if (vr$pass) {
              list(status = "PASS", reason = vr$reason,
                   n_exons_original = vr$n_exons_original,
                   n_exons_reconstructed = vr$n_exons_reconstructed)
            } else {
              list(status = "FAIL", reason = vr$reason,
                   n_exons_original = vr$n_exons_original,
                   n_exons_reconstructed = vr$n_exons_reconstructed)
            }
          }
        }
      }
    }, error = function(e) {
      list(status = "ERROR", reason = paste0("Exception: ", e$message),
           n_exons_original = NA_integer_, n_exons_reconstructed = NA_integer_)
    })

    # Track counts
    if (result$status == "PASS") pass_count <- pass_count + 1L
    else if (result$status == "FAIL") fail_count <- fail_count + 1L
    else error_count <- error_count + 1L

    verification_results[[length(verification_results) + 1]] <- tibble(
      gene_id = gene,
      dominant_isoform_id = dom_id,
      comparator_isoform_id = comp_id,
      status = result$status,
      reason = result$reason,
      n_exons_original = result$n_exons_original,
      n_exons_reconstructed = result$n_exons_reconstructed
    )
  }

  # Progress report
  processed <- end_idx
  cat(sprintf("  Batch %d/%d: %d pairs processed (%.1f%%) | PASS: %d  FAIL: %d  ERROR: %d\n",
              batch_idx, n_batches, processed,
              100 * processed / n_pairs,
              pass_count, fail_count, error_count))
}

})  # end suppressWarnings

# Combine results
verification_df <- bind_rows(verification_results)

# ==============================================================================
# 5. Summary Statistics
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("VERIFICATION SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

total <- nrow(verification_df)
cat(sprintf("Total pairs verified: %d\n\n", total))
cat(sprintf("  ✓ PASS:  %d/%d (%.1f%%)\n", pass_count, total, 100 * pass_count / total))
cat(sprintf("  ✗ FAIL:  %d/%d (%.1f%%)\n", fail_count, total, 100 * fail_count / total))
cat(sprintf("  ⚠ ERROR: %d/%d (%.1f%%)\n", error_count, total, 100 * error_count / total))

if (fail_count > 0) {
  cat("\nFailed pairs:\n")
  failures <- verification_df %>% filter(status == "FAIL")
  for (i in seq_len(min(nrow(failures), 20))) {
    f <- failures[i, ]
    cat(sprintf("  [FAIL] %s: %s vs %s — %s\n",
                f$gene_id, f$dominant_isoform_id, f$comparator_isoform_id, f$reason))
  }
  if (nrow(failures) > 20) {
    cat(sprintf("  ... and %d more\n", nrow(failures) - 20))
  }
}

if (error_count > 0) {
  cat("\nError summary:\n")
  error_reasons <- verification_df %>%
    filter(status == "ERROR") %>%
    count(reason, sort = TRUE)
  for (i in seq_len(nrow(error_reasons))) {
    cat(sprintf("  %s: %d\n", error_reasons$reason[i], error_reasons$n[i]))
  }
}

if (pass_count == total) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════════\n")
  cat("  PERFECT RECONSTRUCTION — All pairs match!\n")
  cat("═══════════════════════════════════════════════════════════════════\n")
}

# ==============================================================================
# 6. Save Results
# ==============================================================================

cat("\nSaving results...\n")

saveRDS(verification_df, "data/reconstruction_verification.rds")
cat("  ✓ data/reconstruction_verification.rds\n")

# Save log
dir.create("logs", showWarnings = FALSE)
log_lines <- c(
  sprintf("Reconstruction Verification — %s", Sys.time()),
  sprintf("Total pairs: %d", total),
  sprintf("PASS:  %d (%.1f%%)", pass_count, 100 * pass_count / total),
  sprintf("FAIL:  %d (%.1f%%)", fail_count, 100 * fail_count / total),
  sprintf("ERROR: %d (%.1f%%)", error_count, 100 * error_count / total)
)
writeLines(log_lines, "logs/reconstruction_verification.log")
cat("  ✓ logs/reconstruction_verification.log\n")

cat("\n✓ Step 8 complete\n\n")
