#!/usr/bin/env Rscript
library(tidyverse)

source("../scripts/event_detection_functions_v2.R")
source("visualization_functions.R")

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   V2 VALIDATION - All 33 Tests (Fixed Minus Strand)         ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Load data
gtf_df <- parse_gtf("synthetic/TestData/exons/base_events.gtf")
expected <- read_tsv("synthetic/TestData/annotations/base_events.tsv",
                     comment = "#", show_col_types = FALSE)

cat(sprintf("Loaded %d exons, %d test cases\n\n", nrow(gtf_df), nrow(expected)))

# Build isoform structures
isoform_structures <- gtf_df %>%
  group_by(gene_id, transcript_id) %>%
  arrange(exon_number) %>%
  summarise(
    seqnames = first(seqnames),
    strand = first(strand),
    exons = list(tibble(
      exon_number = exon_number,
      exon_start = start,
      exon_end = end
    )),
    .groups = "drop"
  )

# Run validation
results <- list()
for (i in seq_len(nrow(expected))) {
  test_case <- expected[i, ]

  cat(sprintf("[%d/%d] %s\n", i, nrow(expected), test_case$gene_id))

  iso_a <- test_case$isoform_A
  iso_b <- test_case$isoform_B

  struct_a <- isoform_structures %>% filter(transcript_id == iso_a) %>% pull(exons) %>% .[[1]]
  struct_b <- isoform_structures %>% filter(transcript_id == iso_b) %>% pull(exons) %>% .[[1]]
  strand <- isoform_structures %>% filter(transcript_id == iso_a) %>% pull(strand)

  if (nrow(struct_a) == 0 || nrow(struct_b) == 0) {
    cat("      ⚠ Missing data\n\n")
    next
  }

  events <- detect_splicing_events_v2(struct_a, struct_b, strand)
  detected <- summarize_events_v2(events)

  match <- (detected == test_case$event_type)

  results[[length(results) + 1]] <- tibble(
    gene_id = test_case$gene_id,
    expected = test_case$event_type,
    detected = detected,
    match = match
  )

  cat(sprintf("      Expected: %s\n", test_case$event_type))
  cat(sprintf("      Detected: %s\n", detected))
  cat(sprintf("      Status: %s\n\n", if(match) "✓ PASS" else "✗ FAIL"))
}

results_df <- bind_rows(results)

# Summary
cat("═══════════════════════════════════════════════════════════════════\n")
cat("VALIDATION SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

n_total <- nrow(results_df)
n_pass <- sum(results_df$match)
pct_pass <- round(100 * n_pass / n_total, 1)

cat(sprintf("Total: %d | ✓ PASS: %d (%.1f%%) | ✗ FAIL: %d\n\n",
            n_total, n_pass, pct_pass, n_total - n_pass))

if (n_pass < n_total) {
  cat("FAILED CASES:\n")
  print(results_df %>% filter(!match) %>% select(gene_id, expected, detected))
}

cat("\n")
if (n_pass == n_total) {
  cat("✓ ALL TESTS PASSED!\n")
} else {
  cat(sprintf("⚠ %d/%d tests need review\n", n_total - n_pass, n_total))
}
cat("═══════════════════════════════════════════════════════════════════\n")
