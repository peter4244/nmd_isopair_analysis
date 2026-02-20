#!/usr/bin/env Rscript
################################################################################
# Test Script: Run 08_validate_reconstruction.R logic on GENCODE test data
################################################################################
#
# Adapts the GENCODE real_data GTF + events into the same code path used by
# Script 08, to verify that reconstruct_dominant_v2() + verify_transcript()
# produce the expected ~99.6% pass rate.
#
# Inputs (from development/reconstruction/real_data/):
#   - gencode_test.gtf        (all isoform exon structures)
#   - events_v4.tsv           (events detected by event_detection_functions.R)
#   - union_exons.tsv.gz      (atomic union exons)
#
################################################################################

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   TEST: Script 08 Logic on GENCODE Real Data                 ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"
setwd(base_dir)

# ==============================================================================
# 0. Source reconstruction functions (same as Script 08)
# ==============================================================================

source("scripts/reconstruction_functions.R")
cat("  Loaded reconstruction_functions.R\n")

# Verification tolerance — same as Script 08
TSS_TOLERANCE <- 20
TES_TOLERANCE <- 20

# ==============================================================================
# 1. Parse GENCODE GTF → isoform exon structures
# ==============================================================================

cat("\nParsing GENCODE GTF...\n")
gtf_file <- "development/reconstruction/real_data/gencode_test.gtf"

lines <- readLines(gtf_file)
lines <- lines[!startsWith(lines, '#') & lines != '']

gtf_results <- list()
for (line in lines) {
  parts <- strsplit(line, '\t')[[1]]
  if (parts[3] != "exon") next

  attrs_str <- parts[9]
  gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', attrs_str)
  transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', attrs_str)

  gtf_results[[length(gtf_results) + 1]] <- tibble(
    chr = parts[1],
    exon_start = as.integer(parts[4]),
    exon_end = as.integer(parts[5]),
    strand = parts[7],
    gene_id = gene_id,
    transcript_id = transcript_id
  )
}

isoform_structures <- bind_rows(gtf_results)
cat(sprintf("  Parsed %d exons, %d transcripts, %d genes\n",
            nrow(isoform_structures),
            n_distinct(isoform_structures$transcript_id),
            n_distinct(isoform_structures$gene_id)))

# ==============================================================================
# 2. Load events + union exons
# ==============================================================================

cat("\nLoading events...\n")
all_events <- read_tsv("development/reconstruction/real_data/events_v4.tsv",
                       show_col_types = FALSE)
cat(sprintf("  Events: %d\n", nrow(all_events)))

cat("Loading union exons...\n")
union_exons <- read_tsv("development/reconstruction/real_data/union_exons.tsv.gz",
                        show_col_types = FALSE,
                        col_names = c("chr", "start", "end", "union_exon_id", "strand", "gene_id"),
                        comment = "#")
cat(sprintf("  Union exons: %d\n", nrow(union_exons)))

# ==============================================================================
# 3. Verify_transcript function (identical to Script 08)
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
# 4. Reconstruct + Verify (same logic as Script 08)
# ==============================================================================

cat("\nReconstructing and verifying...\n")

# Get unique pairs from events
unique_pairs <- all_events %>%
  select(gene_id, dominant_transcript_id, comparator_transcript_id) %>%
  distinct()

cat(sprintf("  Unique pairs: %d\n", nrow(unique_pairs)))

verification_results <- list()
pass_count <- 0L
fail_count <- 0L
error_count <- 0L

suppressWarnings({

for (row_i in seq_len(nrow(unique_pairs))) {
  pair <- unique_pairs[row_i, ]
  dom_id <- pair$dominant_transcript_id
  comp_id <- pair$comparator_transcript_id
  gene <- pair$gene_id

  result <- tryCatch({

    # Extract comparator exons
    comp_exons <- isoform_structures %>%
      filter(transcript_id == comp_id) %>%
      select(chr, exon_start, exon_end, strand, gene_id, transcript_id)

    if (nrow(comp_exons) == 0) {
      list(status = "ERROR", reason = "No comparator exons",
           n_exons_original = NA_integer_, n_exons_reconstructed = NA_integer_)
    } else {

      # Extract dominant exons (for verification)
      dom_exons <- isoform_structures %>%
        filter(transcript_id == dom_id) %>%
        arrange(exon_start)

      if (nrow(dom_exons) == 0) {
        list(status = "ERROR", reason = "No dominant exons",
             n_exons_original = NA_integer_, n_exons_reconstructed = NA_integer_)
      } else {

        strand <- dom_exons$strand[1]

        # Filter events for this pair
        pair_events <- all_events %>%
          filter(dominant_transcript_id == dom_id,
                 comparator_transcript_id == comp_id)

        # Filter union exons for this gene
        gene_union_exons <- union_exons %>%
          filter(gene_id == gene)

        # Reconstruct
        reconstructed <- reconstruct_dominant_v2(
          comp_exons, pair_events, gene_union_exons
        )

        # Verify
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

  if (result$status == "PASS") pass_count <- pass_count + 1L
  else if (result$status == "FAIL") fail_count <- fail_count + 1L
  else error_count <- error_count + 1L

  verification_results[[length(verification_results) + 1]] <- tibble(
    gene_id = gene,
    dominant_transcript_id = dom_id,
    comparator_transcript_id = comp_id,
    status = result$status,
    reason = result$reason,
    n_exons_original = result$n_exons_original,
    n_exons_reconstructed = result$n_exons_reconstructed
  )

  # Progress every 500
  if (row_i %% 500 == 0 || row_i == nrow(unique_pairs)) {
    cat(sprintf("  %d/%d (%.1f%%) | PASS: %d  FAIL: %d  ERROR: %d\n",
                row_i, nrow(unique_pairs),
                100 * row_i / nrow(unique_pairs),
                pass_count, fail_count, error_count))
  }
}

})

verification_df <- bind_rows(verification_results)

# ==============================================================================
# 5. Summary
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("GENCODE TEST — VERIFICATION SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

total <- nrow(verification_df)
cat(sprintf("Total pairs: %d\n\n", total))
cat(sprintf("  PASS:  %d/%d (%.1f%%)\n", pass_count, total, 100 * pass_count / total))
cat(sprintf("  FAIL:  %d/%d (%.1f%%)\n", fail_count, total, 100 * fail_count / total))
cat(sprintf("  ERROR: %d/%d (%.1f%%)\n", error_count, total, 100 * error_count / total))

cat("\nExpected: ~4258/4274 PASS (99.6%), 0 FAIL, 16 ERROR\n")

if (fail_count > 0) {
  cat("\nFailed pairs:\n")
  failures <- verification_df %>% filter(status == "FAIL")
  for (i in seq_len(min(nrow(failures), 30))) {
    f <- failures[i, ]
    cat(sprintf("  [FAIL] %s: %s vs %s — %s\n",
                f$gene_id, f$dominant_transcript_id, f$comparator_transcript_id, f$reason))
  }
  if (nrow(failures) > 30) {
    cat(sprintf("  ... and %d more\n", nrow(failures) - 30))
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

cat("\n═══════════════════════════════════════════════════════════════════\n")

# Save for comparison
saveRDS(verification_df, "data/gencode_test_verification.rds")
cat("Saved: data/gencode_test_verification.rds\n")
