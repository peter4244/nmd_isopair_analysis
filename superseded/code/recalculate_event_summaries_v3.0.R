#!/usr/bin/env Rscript
# Recalculate Event Summaries from Fixed Event Vectors v3.0
# Generates summary statistics after removing invalid events

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   RECALCULATE EVENT SUMMARIES v3.0 - From Fixed Vectors      ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

output_dir <- "results/isoform_transitions/v3.0_reference_based"

# ============================================================================
# Process FILTERED version
# ============================================================================

cat("═══ Processing FILTERED Version ═══\n\n")

cat("Loading filtered event vectors...\n")
filtered <- readRDS(file.path(output_dir, "reference_event_vectors_v3.0_filtered.rds"))
cat("  Loaded:", nrow(filtered), "transitions\n\n")

# Calculate summary statistics
cat("Calculating summary statistics...\n")
filtered_summary <- filtered %>%
  summarize(
    n_transitions = n(),
    total_alt_tss = sum(n_alt_tss),
    total_alt_tes = sum(n_alt_tes),
    total_se = sum(n_se),
    total_a5ss = sum(n_a5ss),
    total_a3ss = sum(n_a3ss),
    total_ir = sum(n_ir),
    total_constitutive = sum(n_constitutive),
    pct_with_alt_tss = 100 * mean(n_alt_tss > 0),
    pct_with_alt_tes = 100 * mean(n_alt_tes > 0),
    pct_with_se = 100 * mean(n_se > 0),
    pct_with_a5ss = 100 * mean(n_a5ss > 0),
    pct_with_a3ss = 100 * mean(n_a3ss > 0),
    pct_with_ir = 100 * mean(n_ir > 0)
  )

print(filtered_summary)
cat("\n")

# Save filtered summary
write_tsv(filtered_summary, file.path(output_dir, "reference_event_summary_filtered.tsv"))
cat("Saved: reference_event_summary_filtered.tsv\n\n")

# ============================================================================
# Process ALL version
# ============================================================================

cat("═══ Processing ALL Version ═══\n\n")

cat("Loading ALL event vectors...\n")
all_data <- readRDS(file.path(output_dir, "reference_event_vectors_v3.0_all.rds"))
cat("  Loaded:", nrow(all_data), "transitions\n\n")

# Calculate summary statistics
cat("Calculating summary statistics...\n")
all_summary <- all_data %>%
  summarize(
    n_transitions = n(),
    total_alt_tss = sum(n_alt_tss),
    total_alt_tes = sum(n_alt_tes),
    total_se = sum(n_se),
    total_a5ss = sum(n_a5ss),
    total_a3ss = sum(n_a3ss),
    total_ir = sum(n_ir),
    total_constitutive = sum(n_constitutive),
    pct_with_alt_tss = 100 * mean(n_alt_tss > 0),
    pct_with_alt_tes = 100 * mean(n_alt_tes > 0),
    pct_with_se = 100 * mean(n_se > 0),
    pct_with_a5ss = 100 * mean(n_a5ss > 0),
    pct_with_a3ss = 100 * mean(n_a3ss > 0),
    pct_with_ir = 100 * mean(n_ir > 0)
  )

print(all_summary)
cat("\n")

# Save ALL summary
write_tsv(all_summary, file.path(output_dir, "reference_event_summary_all.tsv"))
cat("Saved: reference_event_summary_all.tsv\n\n")

# ============================================================================
# Comparison Summary
# ============================================================================

cat("═══ Comparison: ALL vs FILTERED ═══\n\n")

comparison <- tibble(
  version = c("ALL", "FILTERED"),
  n_transitions = c(all_summary$n_transitions, filtered_summary$n_transitions),
  total_se = c(all_summary$total_se, filtered_summary$total_se),
  total_a5ss = c(all_summary$total_a5ss, filtered_summary$total_a5ss),
  total_a3ss = c(all_summary$total_a3ss, filtered_summary$total_a3ss),
  total_ir = c(all_summary$total_ir, filtered_summary$total_ir)
)

print(comparison)
cat("\n")

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EVENT SUMMARY RECALCULATION COMPLETE                       ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("All summaries calculated from FIXED event vectors (0 invalid events)\n\n")
