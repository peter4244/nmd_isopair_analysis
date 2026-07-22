#!/usr/bin/env Rscript

# Quick check of filter results

library(tidyverse)

# Load filter statistics
filter_summary <- read_tsv("../tmp/isoform_filter_summary_2026.1.3.tsv", show_col_types = FALSE)
sig_coverage <- read_tsv("../tmp/isoform_filter_sig_coverage_2026.1.3.tsv", show_col_types = FALSE)

cat("=== ISOFORM FILTER SUMMARY ===\n\n")
cat("Filter: >= 0.25 CPM in EITHER DMSO or Smg1i per cell type\n\n")

print(filter_summary)
cat("\n")

cat("=== SIGNIFICANT ISOFORM COVERAGE ===\n\n")
print(sig_coverage)
cat("\n")

cat("Total significant isoforms retained:",
    sum(sig_coverage$n_in_dictionary), "/", sum(sig_coverage$n_significant),
    paste0("(", round(100 * sum(sig_coverage$n_in_dictionary) / sum(sig_coverage$n_significant), 2), "%)"), "\n")
