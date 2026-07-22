#!/usr/bin/env Rscript
# Analyze A5SS/A3SS Size Categories and Partial IR
# Analysis of v2.1 enhanced event detection

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   SPLICE SITE SIZE ANALYSIS                                    ║\n")
cat("║   Categorizing A5SS/A3SS by shift distance                     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

output_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions"

# Load v2.1 data
cat("Loading enhanced event data (v2.1)...\n")
transitions <- readRDS(file.path(output_dir, "detailed_event_vectors_v2.1.rds"))
cat("Loaded:", nrow(transitions), "transitions\n\n")

# ============================================================================
# 1. A5SS Size Distribution
# ============================================================================

cat("═══ A5SS Size Distribution ═══\n\n")

a5ss_summary <- transitions %>%
  filter(n_a5ss > 0) %>%
  group_by(a5ss_category) %>%
  summarize(
    n_transitions = n(),
    pct_of_a5ss = 100 * n() / sum(transitions$n_a5ss > 0),
    mean_shift = mean(max_a5ss_shift),
    median_shift = median(max_a5ss_shift),
    min_shift = min(max_a5ss_shift),
    max_shift = max(max_a5ss_shift),
    .groups = "drop"
  ) %>%
  arrange(desc(n_transitions))

cat("A5SS events by size category:\n")
print(a5ss_summary, n = 20)
cat("\n")

# ============================================================================
# 2. A3SS Size Distribution
# ============================================================================

cat("═══ A3SS Size Distribution ═══\n\n")

a3ss_summary <- transitions %>%
  filter(n_a3ss > 0) %>%
  group_by(a3ss_category) %>%
  summarize(
    n_transitions = n(),
    pct_of_a3ss = 100 * n() / sum(transitions$n_a3ss > 0),
    mean_shift = mean(max_a3ss_shift),
    median_shift = median(max_a3ss_shift),
    min_shift = min(max_a3ss_shift),
    max_shift = max(max_a3ss_shift),
    .groups = "drop"
  ) %>%
  arrange(desc(n_transitions))

cat("A3SS events by size category:\n")
print(a3ss_summary, n = 20)
cat("\n")

# ============================================================================
# 3. Partial IR vs Full IR
# ============================================================================

cat("═══ Intron Retention Analysis ═══\n\n")

ir_summary <- transitions %>%
  summarize(
    total_transitions = n(),
    n_partial_ir = sum(has_partial_ir),
    pct_partial_ir = 100 * mean(has_partial_ir),
    n_full_ir = sum(has_full_ir),
    pct_full_ir = 100 * mean(has_full_ir),
    n_either_ir = sum(has_partial_ir | has_full_ir),
    pct_either_ir = 100 * mean(has_partial_ir | has_full_ir),
    n_both_ir = sum(has_partial_ir & has_full_ir),
    pct_both_ir = 100 * mean(has_partial_ir & has_full_ir)
  )

cat("Intron retention summary:\n")
print(ir_summary)
cat("\n")

# By cell type
ir_by_celltype <- transitions %>%
  group_by(cell_type) %>%
  summarize(
    n_transitions = n(),
    n_partial_ir = sum(has_partial_ir),
    pct_partial_ir = 100 * mean(has_partial_ir),
    n_full_ir = sum(has_full_ir),
    pct_full_ir = 100 * mean(has_full_ir),
    .groups = "drop"
  )

cat("IR by cell type:\n")
print(ir_by_celltype)
cat("\n")

# ============================================================================
# 4. Combined Size Categories
# ============================================================================

cat("═══ Combined A5SS/A3SS Size Analysis ═══\n\n")

size_combo <- transitions %>%
  filter(n_a5ss > 0 | n_a3ss > 0) %>%
  count(a5ss_category, a3ss_category) %>%
  arrange(desc(n))

cat("A5SS × A3SS category combinations:\n")
print(size_combo, n = 30)
cat("\n")

# ============================================================================
# 5. Shift Distance Distributions
# ============================================================================

cat("═══ Shift Distance Distributions ═══\n\n")

# A5SS distances
a5ss_distances <- transitions %>%
  filter(n_a5ss > 0) %>%
  pull(max_a5ss_shift)

cat("A5SS shift distances (n =", length(a5ss_distances), "):\n")
cat("  Min:", min(a5ss_distances), "bp\n")
cat("  Q1:", quantile(a5ss_distances, 0.25), "bp\n")
cat("  Median:", median(a5ss_distances), "bp\n")
cat("  Q3:", quantile(a5ss_distances, 0.75), "bp\n")
cat("  Max:", max(a5ss_distances), "bp\n")
cat("  Mean:", round(mean(a5ss_distances), 1), "bp\n\n")

# A3SS distances
a3ss_distances <- transitions %>%
  filter(n_a3ss > 0) %>%
  pull(max_a3ss_shift)

cat("A3SS shift distances (n =", length(a3ss_distances), "):\n")
cat("  Min:", min(a3ss_distances), "bp\n")
cat("  Q1:", quantile(a3ss_distances, 0.25), "bp\n")
cat("  Median:", median(a3ss_distances), "bp\n")
cat("  Q3:", quantile(a3ss_distances, 0.75), "bp\n")
cat("  Max:", max(a3ss_distances), "bp\n")
cat("  Mean:", round(mean(a3ss_distances), 1), "bp\n\n")

# ============================================================================
# 6. Save Results
# ============================================================================

write_tsv(a5ss_summary, file.path(output_dir, "a5ss_size_categories.tsv"))
write_tsv(a3ss_summary, file.path(output_dir, "a3ss_size_categories.tsv"))
write_tsv(ir_summary, file.path(output_dir, "ir_summary.tsv"))
write_tsv(ir_by_celltype, file.path(output_dir, "ir_by_celltype.tsv"))
write_tsv(size_combo, file.path(output_dir, "a5ss_a3ss_size_combinations.tsv"))

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║           SIZE ANALYSIS COMPLETE                               ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("Saved results:\n")
cat("  - a5ss_size_categories.tsv\n")
cat("  - a3ss_size_categories.tsv\n")
cat("  - ir_summary.tsv\n")
cat("  - ir_by_celltype.tsv\n")
cat("  - a5ss_a3ss_size_combinations.tsv\n")
cat("\n")
