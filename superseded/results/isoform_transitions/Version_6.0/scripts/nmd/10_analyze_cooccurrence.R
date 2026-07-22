#!/usr/bin/env Rscript
#
# Script: 10_analyze_cooccurrence.R
# Version: 6.0
# Purpose: Test co-occurrence of splicing events
#
# Input:
#   - splicing_choice_profiles_with_bins.rds (from Script 09)
#
# Output:
#   - {output_dir}/cooccurrence_crude_results.tsv
#   - {output_dir}/cooccurrence_stratified_results.tsv
#   - {output_dir}/figures/cooccurrence_heatmap.pdf
#
# Flags:
#   --input path.rds       Input profiles with bins (default: data/splicing_choice_profiles_with_bins.rds)
#   --output-dir dir/      Output directory (default: results)
#

library(tidyverse)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 10: Analyze Event Co-occurrence                        ║\n")
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

cat("Loading splicing profiles...\n")

if (!file.exists(input_path)) {
  stop(sprintf("ERROR: Input file '%s' not found. Run Script 09 first.", input_path))
}

profiles <- tryCatch({
  readRDS(input_path)
}, error = function(e) {
  stop(sprintf("ERROR: Failed to load profiles: %s", e$message))
})

cat(sprintf("  Loaded %d profiles\n", nrow(profiles)))
cat(sprintf("  Covering %d genes\n", n_distinct(profiles$gene_id)))

# ═══════════════════════════════════════════════════════════════════
# SECTION 2: Prepare Event Binary Matrix
# ═══════════════════════════════════════════════════════════════════

cat("\nPreparing binary event matrix...\n")

# Convert event counts to binary (present/absent)
profiles <- profiles %>%
  mutate(
    has_alt_tss = as.integer(tss_changed),
    has_alt_tes = as.integer(tes_changed),
    has_a5ss = as.integer(n_a5ss > 0),
    has_a3ss = as.integer(n_a3ss > 0),
    has_partial_ir = as.integer(n_partial_ir > 0),
    has_ir = as.integer((n_ir + n_ir_diff) > 0),
    has_se = as.integer(n_se > 0),
    has_missing_internal = as.integer(n_missing_internal > 0)
  )

# Event types for analysis
event_types <- c("has_alt_tss", "has_alt_tes", "has_a5ss", "has_a3ss",
                 "has_partial_ir", "has_ir", "has_se", "has_missing_internal")

# Pretty names for output
event_names <- c(
  has_alt_tss = "Alt_TSS",
  has_alt_tes = "Alt_TES",
  has_a5ss = "A5SS",
  has_a3ss = "A3SS",
  has_partial_ir = "Partial_IR",
  has_ir = "IR",
  has_se = "SE",
  has_missing_internal = "Missing_Internal"
)

# Calculate event frequencies
event_freq <- profiles %>%
  summarise(across(all_of(event_types), sum)) %>%
  pivot_longer(everything(), names_to = "event", values_to = "count") %>%
  mutate(
    freq = count / nrow(profiles),
    event_name = event_names[event]
  )

cat("\nEvent frequencies:\n")
print(event_freq %>% select(event_name, count, freq), n = Inf)

# ═══════════════════════════════════════════════════════════════════
# SECTION 3: Crude Co-occurrence Analysis (All Profiles)
# ═══════════════════════════════════════════════════════════════════

cat("\nRunning crude co-occurrence tests...\n")

# Function to test co-occurrence between two events
test_cooccurrence <- function(data, event1, event2) {

  # Build 2x2 contingency table with explicit levels to ensure 2x2
  tbl <- table(
    factor(data[[event1]], levels = c(0, 1)),
    factor(data[[event2]], levels = c(0, 1))
  )

  # Fisher's exact test (works for all sample sizes)
  test <- fisher.test(tbl)

  # Calculate observed vs expected
  n <- nrow(data)
  p1 <- sum(data[[event1]]) / n
  p2 <- sum(data[[event2]]) / n

  observed <- sum(data[[event1]] & data[[event2]]) / n
  expected_indep <- p1 * p2

  # Return results
  tibble(
    event1_name = event_names[event1],
    event2_name = event_names[event2],
    n_both = sum(data[[event1]] & data[[event2]]),
    n_event1_only = sum(data[[event1]] & !data[[event2]]),
    n_event2_only = sum(!data[[event1]] & data[[event2]]),
    n_neither = sum(!data[[event1]] & !data[[event2]]),
    freq_both_observed = observed,
    freq_both_expected = expected_indep,
    enrichment = if (expected_indep > 0) observed / expected_indep else NA_real_,
    odds_ratio = test$estimate,
    or_ci_lower = test$conf.int[1],
    or_ci_upper = test$conf.int[2],
    p_value = test$p.value
  )
}

# Test all pairs
crude_results <- list()
counter <- 1

for (i in 1:(length(event_types) - 1)) {
  for (j in (i + 1):length(event_types)) {
    event1 <- event_types[i]
    event2 <- event_types[j]

    crude_results[[counter]] <- test_cooccurrence(profiles, event1, event2)
    counter <- counter + 1
  }
}

crude_results <- bind_rows(crude_results) %>%
  mutate(
    adj_p_value = p.adjust(p_value, method = "fdr"),
    significant = adj_p_value < 0.05,
    interpretation = case_when(
      !significant ~ "No association",
      enrichment > 1.2 ~ "Co-occur",
      enrichment < 0.8 ~ "Mutually exclusive",
      TRUE ~ "Weak association"
    )
  )

cat(sprintf("  Tested %d event pairs\n", nrow(crude_results)))
cat(sprintf("  Significant associations: %d (%.1f%%)\n",
            sum(crude_results$significant),
            100 * sum(crude_results$significant) / nrow(crude_results)))

# ═══════════════════════════════════════════════════════════════════
# SECTION 4: Stratified Co-occurrence (By Complexity Bin)
# ═══════════════════════════════════════════════════════════════════

cat("\nRunning stratified co-occurrence tests (by complexity)...\n")

# Test within each complexity bin
stratified_results <- list()
counter <- 1

for (bin in unique(profiles$complexity_bin)) {

  bin_data <- profiles %>% filter(complexity_bin == bin)
  cat(sprintf("  Testing bin: %s (n = %d)\n", bin, nrow(bin_data)))

  for (i in 1:(length(event_types) - 1)) {
    for (j in (i + 1):length(event_types)) {
      event1 <- event_types[i]
      event2 <- event_types[j]

      result <- test_cooccurrence(bin_data, event1, event2)
      result$complexity_bin <- bin

      stratified_results[[counter]] <- result
      counter <- counter + 1
    }
  }
}

stratified_results <- bind_rows(stratified_results) %>%
  group_by(complexity_bin) %>%
  mutate(
    adj_p_value = p.adjust(p_value, method = "fdr"),
    significant = adj_p_value < 0.05,
    interpretation = case_when(
      !significant ~ "No association",
      enrichment > 1.2 ~ "Co-occur",
      enrichment < 0.8 ~ "Mutually exclusive",
      TRUE ~ "Weak association"
    )
  ) %>%
  ungroup()

cat(sprintf("  Tested %d event pairs across %d bins\n",
            nrow(stratified_results),
            n_distinct(stratified_results$complexity_bin)))

# ═══════════════════════════════════════════════════════════════════
# SECTION 5: Summary Statistics
# ═══════════════════════════════════════════════════════════════════

cat("\nSummary of co-occurrence patterns:\n")

# Overall patterns
cat("\nCrude analysis (all profiles):\n")
crude_summary <- crude_results %>%
  group_by(interpretation) %>%
  summarise(n_pairs = n(), .groups = "drop") %>%
  arrange(desc(n_pairs))
print(crude_summary)

# Stratified patterns
cat("\nStratified analysis (by complexity):\n")
strat_summary <- stratified_results %>%
  group_by(complexity_bin, interpretation) %>%
  summarise(n_pairs = n(), .groups = "drop") %>%
  arrange(complexity_bin, desc(n_pairs))
print(strat_summary, n = Inf)

# ═══════════════════════════════════════════════════════════════════
# SECTION 6: Visualization
# ═══════════════════════════════════════════════════════════════════

cat("\nCreating heatmaps...\n")

fig_dir <- file.path(output_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(fig_dir, "cooccurrence_heatmap.pdf"), width = 14, height = 10)

# Helper function to create matrix for heatmap
create_or_matrix <- function(results_df) {
  n_events <- length(event_names)
  or_matrix <- matrix(1, nrow = n_events, ncol = n_events)
  rownames(or_matrix) <- colnames(or_matrix) <- event_names

  for (i in 1:nrow(results_df)) {
    row <- results_df[i, ]
    or_matrix[row$event1_name, row$event2_name] <- row$odds_ratio
    or_matrix[row$event2_name, row$event1_name] <- row$odds_ratio
  }

  return(or_matrix)
}

# Plot 1: Crude co-occurrence heatmap
or_matrix <- create_or_matrix(crude_results)

# Convert to long format for ggplot
plot_data <- as.data.frame(or_matrix) %>%
  rownames_to_column("event1") %>%
  pivot_longer(-event1, names_to = "event2", values_to = "odds_ratio") %>%
  left_join(
    crude_results %>%
      select(event1_name, event2_name, significant) %>%
      rename(event1 = event1_name, event2 = event2_name),
    by = c("event1", "event2")
  ) %>%
  mutate(
    log_or = log2(odds_ratio),
    log_or_capped = pmax(pmin(log_or, 3), -3),  # Cap at ±3 for visualization
    sig_marker = ifelse(significant %in% TRUE, "*", "")
  )

p1 <- ggplot(plot_data, aes(x = event1, y = event2, fill = log_or_capped)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_marker), size = 6, vjust = 0.75) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    limits = c(-3, 3),
    name = "log2(OR)",
    breaks = c(-3, -1.5, 0, 1.5, 3),
    labels = c("≤-3", "-1.5", "0", "1.5", "≥3")
  ) +
  labs(
    title = "Event Co-occurrence: Crude Analysis",
    subtitle = "* indicates FDR < 0.05",
    x = NULL, y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    panel.grid = element_blank()
  )

print(p1)

# Plot 2: Stratified co-occurrence (faceted by complexity bin)
for (bin in unique(stratified_results$complexity_bin)) {

  bin_results <- stratified_results %>% filter(complexity_bin == bin)
  or_matrix_bin <- create_or_matrix(bin_results)

  plot_data_bin <- as.data.frame(or_matrix_bin) %>%
    rownames_to_column("event1") %>%
    pivot_longer(-event1, names_to = "event2", values_to = "odds_ratio") %>%
    left_join(
      bin_results %>%
        select(event1_name, event2_name, significant) %>%
        rename(event1 = event1_name, event2 = event2_name),
      by = c("event1", "event2")
    ) %>%
    mutate(
      log_or = log2(odds_ratio),
      log_or_capped = pmax(pmin(log_or, 3), -3),
      sig_marker = ifelse(significant %in% TRUE, "*", "")
    )

  p <- ggplot(plot_data_bin, aes(x = event1, y = event2, fill = log_or_capped)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sig_marker), size = 6, vjust = 0.75) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = c(-3, 3),
      name = "log2(OR)"
    ) +
    labs(
      title = sprintf("Event Co-occurrence: %s", bin),
      subtitle = "* indicates FDR < 0.05",
      x = NULL, y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 14),
      panel.grid = element_blank()
    )

  print(p)
}

dev.off()

cat(sprintf("  ✓ %s\n", file.path(fig_dir, "cooccurrence_heatmap.pdf")))

# ═══════════════════════════════════════════════════════════════════
# SECTION 7: Save Outputs
# ═══════════════════════════════════════════════════════════════════

cat("\nSaving results...\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save crude results
out_path <- file.path(output_dir, "cooccurrence_crude_results.tsv")
write_tsv(crude_results, out_path)
cat(sprintf("  ✓ %s\n", out_path))

# Save stratified results
out_path <- file.path(output_dir, "cooccurrence_stratified_results.tsv")
write_tsv(stratified_results, out_path)
cat(sprintf("  ✓ %s\n", out_path))

cat("\n✓ Step 10 complete\n\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Event pairs tested: %d\n", nrow(crude_results)))
cat(sprintf("Significant crude associations: %d (%.1f%%)\n",
            sum(crude_results$significant),
            100 * sum(crude_results$significant) / nrow(crude_results)))
cat(sprintf("Stratified tests: %d pairs × %d bins = %d tests\n",
            nrow(crude_results),
            n_distinct(stratified_results$complexity_bin),
            nrow(stratified_results)))
cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")
