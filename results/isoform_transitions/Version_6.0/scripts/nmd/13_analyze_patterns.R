#!/usr/bin/env Rscript
#
# Script: 13_analyze_patterns.R
# Version: 6.0
# Purpose: Classify and analyze splicing choice pattern frequencies (complexity-controlled)
#
# Input:
#   - splicing_choice_profiles_with_bins.rds (from Script 09)
#
# Output:
#   - {output_dir}/pattern_frequencies_by_complexity.tsv
#   - {output_dir}/top_patterns.tsv
#   - {output_dir}/figures/pattern_classification.pdf
#
# Flags:
#   --input path.rds       Input profiles with bins (default: data/splicing_choice_profiles_with_bins.rds)
#   --output-dir dir/      Output directory (default: results)
#

library(tidyverse)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 13: Analyze Pattern Classification                     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Paths
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"
setwd(base_dir)

# ═══════════════════════════════════════════════════════════════════
# Parse CLI arguments
# ═══════════════════════════════════════════════════════════════════

args <- commandArgs(trailingOnly = TRUE)

input_path <- "data/splicing_choice_profiles_with_bins.rds"
if ("--input" %in% args) {
  idx <- which(args == "--input")
  input_path <- args[idx + 1]
}

output_dir <- "results"
if ("--output-dir" %in% args) {
  idx <- which(args == "--output-dir")
  output_dir <- args[idx + 1]
}

cat(sprintf("  Input:      %s\n", input_path))
cat(sprintf("  Output dir: %s\n\n", output_dir))

# ═══════════════════════════════════════════════════════════════════
# SECTION 1: Load Data
# ═══════════════════════════════════════════════════════════════════

cat("Loading data...\n")

if (!file.exists(input_path)) {
  stop(sprintf("ERROR: Input file '%s' not found. Run Script 09 first.", input_path))
}

profiles <- readRDS(input_path)
cat(sprintf("  Loaded %d profiles\n", nrow(profiles)))

# Verify complexity_bin column exists
if (!"complexity_bin" %in% names(profiles)) {
  stop("ERROR: Input profiles do not contain 'complexity_bin' column. Run Script 09 first.")
}

cat(sprintf("  Complexity bins: %s\n", paste(sort(unique(profiles$complexity_bin)), collapse = ", ")))

# ═══════════════════════════════════════════════════════════════════
# SECTION 2: Classify Profiles
# ═══════════════════════════════════════════════════════════════════

cat("\nClassifying profiles into pattern types...\n")

# Define event categories
# Terminal: Alt_TSS, Alt_TES
# Boundary: A5SS, A3SS
# Inclusion: SE, Missing_Internal
# Partial Retention: Partial_IR
# Full Retention: IR, IR_diff

profiles <- profiles %>%
  mutate(
    has_terminal = (n_alt_tss > 0) | (n_alt_tes > 0),
    has_boundary = (n_a5ss > 0) | (n_a3ss > 0),
    has_inclusion = (n_se > 0) | (n_missing_internal > 0),
    has_partial_retention = (n_partial_ir > 0),
    has_full_retention = (n_ir > 0) | (n_ir_diff > 0),
    n_categories = has_terminal + has_boundary + has_inclusion +
                   has_partial_retention + has_full_retention,
    profile_type = case_when(
      n_events == 0 ~ "Identical",
      n_categories == 1 & has_terminal ~ "Terminal-only",
      n_categories == 1 & has_boundary ~ "Boundary-only",
      n_categories == 1 & has_inclusion ~ "Inclusion-only",
      n_categories == 1 & has_partial_retention ~ "Partial_Retention-only",
      n_categories == 1 & has_full_retention ~ "Full_Retention-only",
      n_categories >= 2 ~ "Combined",
      TRUE ~ "Unclassified"
    )
  )

# Order profile types for display
# Include all observed types to avoid NA from factor conversion
type_order <- c("Identical", "Terminal-only", "Boundary-only", "Inclusion-only",
                "Partial_Retention-only", "Full_Retention-only", "Combined")
observed_types <- unique(profiles$profile_type)
extra_types <- setdiff(observed_types, type_order)
if (length(extra_types) > 0) type_order <- c(type_order, extra_types)
profiles$profile_type <- factor(profiles$profile_type, levels = type_order)

type_counts <- profiles %>%
  count(profile_type, .drop = FALSE) %>%
  mutate(pct = 100 * n / sum(n))

cat("\nProfile type distribution:\n")
print(type_counts)

# ═══════════════════════════════════════════════════════════════════
# SECTION 3: Pattern Frequencies by Complexity
# ═══════════════════════════════════════════════════════════════════

cat("\nAnalyzing patterns by complexity bin...\n")

pattern_by_complexity <- profiles %>%
  count(complexity_bin, profile_type, .drop = FALSE) %>%
  group_by(complexity_bin) %>%
  mutate(
    bin_total = sum(n),
    proportion = ifelse(bin_total > 0, n / bin_total, 0)
  ) %>%
  ungroup()

cat("\nPattern frequencies by complexity:\n")
pattern_wide <- pattern_by_complexity %>%
  select(complexity_bin, profile_type, proportion) %>%
  pivot_wider(names_from = profile_type, values_from = proportion, values_fill = 0)
print(pattern_wide)

# ═══════════════════════════════════════════════════════════════════
# SECTION 4: Pattern × Complexity Association
# ═══════════════════════════════════════════════════════════════════

cat("\nTesting pattern-complexity association...\n")

# Build contingency table
contingency <- profiles %>%
  count(complexity_bin, profile_type, .drop = FALSE) %>%
  pivot_wider(names_from = profile_type, values_from = n, values_fill = 0) %>%
  column_to_rownames("complexity_bin") %>%
  as.matrix()

# Check expected cell counts
expected <- chisq.test(contingency)$expected
min_expected <- min(expected)
cat(sprintf("  Minimum expected cell count: %.1f\n", min_expected))

# Use Monte Carlo simulation if any expected count < 5
if (min_expected < 5) {
  cat("  Using Monte Carlo chi-square (some expected counts < 5)\n")
  chisq_result <- chisq.test(contingency, simulate.p.value = TRUE, B = 10000)
} else {
  chisq_result <- chisq.test(contingency)
}

# Cramer's V (effect size)
n_total <- sum(contingency)
k <- min(nrow(contingency), ncol(contingency))
cramers_v <- sqrt(chisq_result$statistic / (n_total * (k - 1)))

cat(sprintf("  Chi-square: %.2f\n", chisq_result$statistic))
cat(sprintf("  p-value: %.2e\n", chisq_result$p.value))
cat(sprintf("  Cramer's V: %.3f\n", cramers_v))

# Post-hoc: standardized residuals
if (min_expected >= 5) {
  std_residuals <- chisq_result$stdres
  cat("\nStandardized residuals (|z| > 2 highlighted):\n")

  residuals_df <- as.data.frame(as.table(std_residuals)) %>%
    rename(complexity_bin = Var1, profile_type = Var2, std_residual = Freq) %>%
    mutate(significant = abs(std_residual) > 2)

  sig_residuals <- residuals_df %>% filter(significant) %>% arrange(desc(abs(std_residual)))
  if (nrow(sig_residuals) > 0) {
    print(sig_residuals)
  } else {
    cat("  No standardized residuals with |z| > 2\n")
  }
} else {
  residuals_df <- NULL
}

association_results <- tibble(
  chi_square = as.numeric(chisq_result$statistic),
  df = ifelse(!is.null(chisq_result$parameter), chisq_result$parameter, NA_real_),
  p_value = chisq_result$p.value,
  cramers_v = as.numeric(cramers_v),
  min_expected_count = min_expected,
  method = ifelse(min_expected < 5, "Monte Carlo (B=10000)", "Pearson chi-square")
)

# ═══════════════════════════════════════════════════════════════════
# SECTION 5: Detailed Combination Analysis
# ═══════════════════════════════════════════════════════════════════

cat("\nAnalyzing 'Combined' profile compositions...\n")

combined_profiles <- profiles %>%
  filter(profile_type == "Combined")

cat(sprintf("  Combined profiles: %d (%.1f%%)\n",
            nrow(combined_profiles),
            100 * nrow(combined_profiles) / nrow(profiles)))

if (nrow(combined_profiles) > 0) {
  # Create event signature string
  combined_profiles <- combined_profiles %>%
    mutate(
      signature = paste0(
        ifelse(has_terminal, "Terminal", ""),
        ifelse(has_terminal & (has_boundary | has_inclusion | has_partial_retention | has_full_retention), "+", ""),
        ifelse(has_boundary, "Boundary", ""),
        ifelse(has_boundary & (has_inclusion | has_partial_retention | has_full_retention), "+", ""),
        ifelse(has_inclusion, "Inclusion", ""),
        ifelse(has_inclusion & (has_partial_retention | has_full_retention), "+", ""),
        ifelse(has_partial_retention, "PartialRetention", ""),
        ifelse(has_partial_retention & has_full_retention, "+", ""),
        ifelse(has_full_retention, "FullRetention", "")
      )
    )

  # Top 10 most frequent combinations
  top_patterns <- combined_profiles %>%
    count(signature, sort = TRUE) %>%
    mutate(
      pct = 100 * n / sum(n),
      cumulative_pct = cumsum(pct)
    ) %>%
    head(10)

  cat("\nTop 10 Combined profile signatures:\n")
  print(top_patterns)

  # Break down by complexity bin
  top_by_complexity <- combined_profiles %>%
    count(complexity_bin, signature) %>%
    group_by(complexity_bin) %>%
    mutate(
      bin_total = sum(n),
      proportion = n / bin_total
    ) %>%
    ungroup() %>%
    arrange(complexity_bin, desc(n))

} else {
  top_patterns <- tibble(signature = character(0), n = integer(0),
                          pct = numeric(0), cumulative_pct = numeric(0))
  top_by_complexity <- tibble(complexity_bin = character(0), signature = character(0),
                               n = integer(0), bin_total = integer(0), proportion = numeric(0))
}

# ═══════════════════════════════════════════════════════════════════
# SECTION 6: Visualization
# ═══════════════════════════════════════════════════════════════════

cat("\nCreating pattern plots...\n")

fig_dir <- file.path(output_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(fig_dir, "pattern_classification.pdf"), width = 14, height = 10)

# Define color palette for profile types
type_colors <- c(
  "Identical" = "#999999",
  "Terminal-only" = "#E69F00",
  "Boundary-only" = "#56B4E9",
  "Inclusion-only" = "#009E73",
  "Partial_Retention-only" = "#F0E442",
  "Full_Retention-only" = "#D55E00",
  "Combined" = "#CC79A7"
)

# Plot 1: Overall profile type distribution
p1 <- ggplot(type_counts, aes(x = reorder(profile_type, -n), y = n, fill = profile_type)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = type_colors) +
  labs(
    title = "Profile Type Distribution",
    subtitle = sprintf("n = %d profiles", nrow(profiles)),
    x = "Profile type",
    y = "Count"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(face = "bold"))

print(p1)

# Plot 2: Stacked barplot — profile type proportions by complexity bin
p2 <- ggplot(pattern_by_complexity, aes(x = complexity_bin, y = proportion, fill = profile_type)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = type_colors) +
  labs(
    title = "Profile Type Proportions by Complexity Bin",
    subtitle = sprintf("Chi-sq = %.1f, p = %.2e, Cramer's V = %.3f",
                       chisq_result$statistic, chisq_result$p.value, cramers_v),
    x = "Complexity bin",
    y = "Proportion",
    fill = "Profile type"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

print(p2)

# Plot 3: Mosaic-style — standardized residuals heatmap
if (!is.null(residuals_df)) {
  p3 <- ggplot(residuals_df, aes(x = complexity_bin, y = profile_type, fill = std_residual)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", std_residual)), size = 3) +
    scale_fill_gradient2(
      low = "steelblue", mid = "white", high = "firebrick",
      midpoint = 0,
      name = "Std.\nResidual"
    ) +
    labs(
      title = "Pattern-Complexity Association: Standardized Residuals",
      subtitle = "Blue = depleted, Red = enriched (|z| > 2 = significant)",
      x = "Complexity bin",
      y = "Profile type"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank())

  print(p3)
}

# Plot 4: Top-10 combination frequency barplot
if (nrow(top_patterns) > 0) {
  p4 <- ggplot(top_patterns, aes(x = reorder(signature, n), y = n)) +
    geom_col(fill = "#CC79A7") +
    geom_text(aes(label = sprintf("%.1f%%", pct)), hjust = -0.1, size = 3.5) +
    coord_flip() +
    labs(
      title = "Top 10 Combined Profile Signatures",
      subtitle = sprintf("Among %d combined profiles", nrow(combined_profiles)),
      x = "Event combination signature",
      y = "Count"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"))

  print(p4)
}

dev.off()

cat(sprintf("  ✓ %s\n", file.path(fig_dir, "pattern_classification.pdf")))

# ═══════════════════════════════════════════════════════════════════
# SECTION 7: Save Outputs
# ═══════════════════════════════════════════════════════════════════

cat("\nSaving results...\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save pattern frequencies by complexity
out_path <- file.path(output_dir, "pattern_frequencies_by_complexity.tsv")
write_tsv(pattern_by_complexity, out_path)
cat(sprintf("  ✓ %s\n", out_path))

# Save association test separately
out_path <- file.path(output_dir, "pattern_complexity_association.tsv")
write_tsv(association_results, out_path)
cat(sprintf("  ✓ %s\n", out_path))

# Save top patterns
out_path <- file.path(output_dir, "top_patterns.tsv")
write_tsv(top_patterns, out_path)
cat(sprintf("  ✓ %s\n", out_path))

# Save top patterns by complexity
out_path <- file.path(output_dir, "top_patterns_by_complexity.tsv")
write_tsv(top_by_complexity, out_path)
cat(sprintf("  ✓ %s\n", out_path))

cat("\n✓ Step 13 complete\n\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Profiles classified: %d\n", nrow(profiles)))
cat(sprintf("Pattern types: %d\n", n_distinct(profiles$profile_type)))
cat(sprintf("Combined profiles: %d (%.1f%%)\n",
            nrow(combined_profiles),
            100 * nrow(combined_profiles) / nrow(profiles)))
cat(sprintf("Association: Chi-sq = %.1f, p = %.2e, V = %.3f\n",
            chisq_result$statistic, chisq_result$p.value, cramers_v))
cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")
