#!/usr/bin/env Rscript
#
# Script: 09_analyze_complexity_relationship.R
# Version: 6.0
# Purpose: Test relationship between isoform complexity and event count (foundational for stratification)
#
# Input:
#   - data/splicing_choice_profiles.rds (from Script 07)
#
# Output:
#   - results/complexity_relationship_results.rds
#   - results/complexity_bins_definition.rds
#   - figures/complexity_vs_events.pdf
#

library(tidyverse)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 9: Analyze Complexity vs Event Count                   ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Paths
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"
setwd(base_dir)

# ═══════════════════════════════════════════════════════════════════
# SECTION 1: Load Data
# ═══════════════════════════════════════════════════════════════════

cat("Loading splicing profiles...\n")

# Check if input file exists
if (!file.exists("data/splicing_choice_profiles.rds")) {
  stop("ERROR: Input file 'data/splicing_choice_profiles.rds' not found. Run Script 07 first.")
}

profiles <- tryCatch({
  readRDS("data/splicing_choice_profiles.rds")
}, error = function(e) {
  stop(sprintf("ERROR: Failed to load profiles: %s", e$message))
})

cat(sprintf("  Loaded %d profiles\n", nrow(profiles)))
cat(sprintf("  Covering %d genes\n", n_distinct(profiles$gene_id)))

# ═══════════════════════════════════════════════════════════════════
# SECTION 2: Test Each Complexity Metric
# ═══════════════════════════════════════════════════════════════════

cat("\nTesting complexity metrics...\n")

# Define complexity metrics to test
complexity_metrics <- c(
  "n_exons_non_dom",           # Number of exons in non-dominant isoform
  "n_junctions_non_dom",       # Number of junctions
  "length_non_dom",            # Transcript length
  "n_union_exons_total"         # Gene-level complexity (union exons in gene)
)

# Event count: total splicing events
# Calculate from actual event counts in the data
profiles <- profiles %>%
  mutate(
    total_events = n_a5ss + n_a3ss + n_ir + n_ir_diff + n_partial_ir + n_se +
                   n_missing_internal + n_alt_tss + n_alt_tes
  )

# Store results
model_results <- list()

for (metric in complexity_metrics) {

  if (!metric %in% names(profiles)) {
    cat(sprintf("  ⚠ Skipping %s (not found in profiles)\n", metric))
    next
  }

  cat(sprintf("  Testing %s...\n", metric))

  # Prepare data
  data_for_model <- profiles %>%
    filter(!is.na(!!sym(metric)), !is.na(total_events))

  # Fit multiple models
  # 1. Linear
  fit_linear <- lm(total_events ~ get(metric), data = data_for_model)
  r2_linear <- summary(fit_linear)$r.squared
  aic_linear <- AIC(fit_linear)

  # 2. Polynomial (degree 2)
  fit_poly <- lm(total_events ~ poly(get(metric), 2), data = data_for_model)
  r2_poly <- summary(fit_poly)$r.squared
  aic_poly <- AIC(fit_poly)

  # 3. Log-linear (if metric > 0)
  if (all(data_for_model[[metric]] > 0)) {
    fit_log <- lm(total_events ~ log(get(metric)), data = data_for_model)
    r2_log <- summary(fit_log)$r.squared
    aic_log <- AIC(fit_log)
  } else {
    r2_log <- NA
    aic_log <- NA
  }

  # Select best model (highest R²)
  r2_values <- c(linear = r2_linear, poly = r2_poly, log = r2_log)
  aic_values <- c(linear = aic_linear, poly = aic_poly, log = aic_log)

  best_model <- names(which.max(r2_values))

  # Store results
  model_results[[metric]] <- list(
    metric = metric,
    n_obs = nrow(data_for_model),
    r2_linear = r2_linear,
    r2_poly = r2_poly,
    r2_log = r2_log,
    aic_linear = aic_linear,
    aic_poly = aic_poly,
    aic_log = aic_log,
    best_model = best_model,
    best_r2 = r2_values[best_model],
    best_aic = aic_values[best_model]
  )

  cat(sprintf("    Best model: %s (R² = %.3f)\n", best_model, r2_values[best_model]))
}

# Convert to data frame
results_df <- bind_rows(model_results)

# ═══════════════════════════════════════════════════════════════════
# SECTION 3: Select Primary Complexity Metric
# ═══════════════════════════════════════════════════════════════════

cat("\nSelecting primary complexity metric...\n")

# Select metric with highest R²
primary_metric <- results_df %>%
  filter(!is.na(best_r2)) %>%
  arrange(desc(best_r2)) %>%
  slice(1)

cat(sprintf("  Primary metric: %s\n", primary_metric$metric))
cat(sprintf("  Best model: %s\n", primary_metric$best_model))
cat(sprintf("  R²: %.3f\n", primary_metric$best_r2))
cat(sprintf("  AIC: %.1f\n", primary_metric$best_aic))

# ═══════════════════════════════════════════════════════════════════
# SECTION 4: Define Complexity Bins
# ═══════════════════════════════════════════════════════════════════

cat("\nDefining complexity bins...\n")

# Use quartiles of primary metric
metric_name <- primary_metric$metric
metric_values <- profiles[[metric_name]]
metric_values <- metric_values[!is.na(metric_values)]

# Calculate quartiles
quartiles <- quantile(metric_values, probs = c(0, 0.25, 0.5, 0.75, 1))

cat("  Quartiles:\n")
cat(sprintf("    Q1 (0-25%%): %.1f - %.1f\n", quartiles[1], quartiles[2]))
cat(sprintf("    Q2 (25-50%%): %.1f - %.1f\n", quartiles[2], quartiles[3]))
cat(sprintf("    Q3 (50-75%%): %.1f - %.1f\n", quartiles[3], quartiles[4]))
cat(sprintf("    Q4 (75-100%%): %.1f - %.1f\n", quartiles[4], quartiles[5]))

# Add complexity bin to profiles
profiles <- profiles %>%
  mutate(
    complexity_bin = case_when(
      !!sym(metric_name) <= quartiles[2] ~ "Q1_low",
      !!sym(metric_name) <= quartiles[3] ~ "Q2_medium_low",
      !!sym(metric_name) <= quartiles[4] ~ "Q3_medium_high",
      TRUE ~ "Q4_high"
    )
  )

# Summary by bin
bin_summary <- profiles %>%
  group_by(complexity_bin) %>%
  summarise(
    n_profiles = n(),
    mean_metric = mean(!!sym(metric_name), na.rm = TRUE),
    mean_events = mean(total_events, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n  Profiles per bin:\n")
print(bin_summary)

# ═══════════════════════════════════════════════════════════════════
# SECTION 5: Visualization
# ═══════════════════════════════════════════════════════════════════

cat("\nCreating plots...\n")

dir.create("figures", showWarnings = FALSE)

pdf("figures/complexity_vs_events.pdf", width = 12, height = 10)

# Plot for each metric
for (metric in names(model_results)) {

  if (!metric %in% names(profiles)) next

  result <- model_results[[metric]]

  p <- ggplot(profiles, aes(x = !!sym(metric), y = total_events)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_smooth(method = "lm", color = "blue", se = TRUE, linetype = "dashed", alpha = 0.2) +
    labs(
      title = sprintf("Complexity vs Event Count: %s", metric),
      subtitle = sprintf("Best model: %s | R² = %.3f | AIC = %.1f",
                        result$best_model, result$best_r2, result$best_aic),
      x = metric,
      y = "Total splicing changes"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"))

  print(p)
}

# Primary metric with bins
p_primary <- ggplot(profiles, aes(x = !!sym(metric_name), y = total_events, color = complexity_bin)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = sprintf("Primary Complexity Metric: %s", metric_name),
    subtitle = "Colored by complexity quartile bins",
    x = metric_name,
    y = "Total splicing changes",
    color = "Complexity Bin"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

print(p_primary)

dev.off()

cat("  ✓ figures/complexity_vs_events.pdf\n")

# ═══════════════════════════════════════════════════════════════════
# SECTION 6: Save Outputs
# ═══════════════════════════════════════════════════════════════════

cat("\nSaving results...\n")

dir.create("results", showWarnings = FALSE)

# Save model results
saveRDS(results_df, "results/complexity_relationship_results.rds")
cat("  ✓ results/complexity_relationship_results.rds\n")

# Save complexity bins definition
complexity_bins <- list(
  primary_metric = metric_name,
  best_model = primary_metric$best_model,
  r_squared = primary_metric$best_r2,
  quartiles = quartiles,
  bin_summary = bin_summary
)

saveRDS(complexity_bins, "results/complexity_bins_definition.rds")
cat("  ✓ results/complexity_bins_definition.rds\n")

# Update profiles with complexity bins
saveRDS(profiles, "data/splicing_choice_profiles_with_bins.rds")
cat("  ✓ data/splicing_choice_profiles_with_bins.rds\n")

cat("\n✓ Step 9 complete\n\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Primary complexity metric: %s\n", metric_name))
cat(sprintf("Best model: %s (R² = %.3f)\n", primary_metric$best_model, primary_metric$best_r2))
cat(sprintf("Complexity bins defined: %d quartiles\n", nrow(bin_summary)))
cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")
