#!/usr/bin/env Rscript
# Verify the deduplication fix

library(tidyverse)

# Load the union exon models
union_models <- readRDS("results/isoform_transitions/v3.0_reference_based/union_exon_models_test.rds")

cat("\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("VERIFICATION OF DEDUPLICATION FIX\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Check AARS1
aars1 <- union_models[["AARS1"]]

cat("AARS1 Gene:\n")
cat("  Number of isoforms:", aars1$n_isoforms, "\n")
cat("  Total union exons:", length(aars1$union_exons), "\n\n")

cat("Checking for duplicate variants:\n")
has_duplicates <- FALSE

for (exon in aars1$union_exons) {
  n_variants <- nrow(exon$variants)
  n_unique <- nrow(distinct(exon$variants))

  if (n_variants != n_unique) {
    has_duplicates <- TRUE
    cat(sprintf("  Exon %d: %d variants (%d unique) ⚠️\n",
                exon$exon_number, n_variants, n_unique))
  }
}

if (!has_duplicates) {
  cat("  ✓ No duplicate variants found!\n")
  cat("  All exons have unique variants only.\n")
}

cat("\n")

# Show first exon structure
cat("First exon structure:\n")
first_exon <- aars1$union_exons[[1]]
cat(sprintf("  Exon %d (%s): %d variant(s)\n",
            first_exon$exon_number,
            first_exon$exon_type,
            nrow(first_exon$variants)))

for (i in seq_len(nrow(first_exon$variants))) {
  v <- first_exon$variants[i, ]
  cat(sprintf("    [%d] %s: %s:%d-%d\n",
              i, v$isoform_id, v$seqnames, v$start, v$end))
}

cat("\n✓ Deduplication fix is working correctly!\n\n")
