#!/usr/bin/env Rscript
# Quick analysis: GAIN/LOSS direction stratified by event type, NMD vs baseline
# Focused on Alt_TES (C1 signal) and SE (C2 signal)

suppressPackageStartupMessages(library(tidyverse))

base_dir <- "comparisons/nonNMD_0.50"

# Discover per-comparison profile files
comparisons <- c("C1", "C2", "C4")
runs <- c("all_samples", "at", "dd", "dd_ali", "fb", "mv")

# Load and unnest events from each comparison x run
all_events <- list()

for (comp in comparisons) {
  for (run in runs) {
    profile_path <- file.path(base_dir, comp, run, "splicing_choice_profiles.rds")
    if (!file.exists(profile_path)) next

    profiles <- readRDS(profile_path)

    events <- profiles %>%
      select(gene_id, dominant_isoform_id, non_dominant_isoform_id, detailed_events) %>%
      unnest(detailed_events, names_sep = "_") %>%
      select(
        gene_id,
        event_type = detailed_events_event_type,
        direction = detailed_events_direction
      ) %>%
      mutate(comparison = comp, run = run)

    all_events[[paste(comp, run, sep = "_")]] <- events
    cat(sprintf("  Loaded %s / %s: %d events\n", comp, run, nrow(events)))
  }
}

all_events <- bind_rows(all_events)
cat(sprintf("\nTotal events: %d\n\n", nrow(all_events)))

# Compute per event-type direction proportions
direction_by_type <- all_events %>%
  filter(direction %in% c("GAIN", "LOSS")) %>%
  group_by(comparison, run, event_type) %>%
  summarise(
    n_events = n(),
    n_loss = sum(direction == "LOSS"),
    n_gain = sum(direction == "GAIN"),
    loss_prop = n_loss / n_events,
    .groups = "drop"
  )

# Focus on Alt_TES and SE: compare NMD (C1, C2) vs baseline (C4)
cat("═══════════════════════════════════════════════════════════════\n")
cat("GAIN/LOSS DIRECTION BY EVENT TYPE: NMD vs BASELINE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

focal_events <- c("Alt_TES", "SE", "Alt_TSS", "A5SS", "A3SS",
                   "Missing_Internal", "IR", "Partial_IR_5", "Partial_IR_3")

results <- list()

for (nmd_comp in c("C1", "C2")) {
  cat(sprintf("\n─── %s vs C4 ───\n\n", nmd_comp))

  for (r in runs) {
    nmd_data <- direction_by_type %>%
      filter(comparison == nmd_comp, run == r)
    baseline_data <- direction_by_type %>%
      filter(comparison == "C4", run == r)

    if (nrow(nmd_data) == 0 || nrow(baseline_data) == 0) next

    cat(sprintf("  Run: %s\n", r))

    for (et in focal_events) {
      nmd_row <- nmd_data %>% filter(event_type == et)
      base_row <- baseline_data %>% filter(event_type == et)

      if (nrow(nmd_row) == 0 || nrow(base_row) == 0) next
      if (nmd_row$n_events < 10 || base_row$n_events < 10) next

      pt <- prop.test(
        c(nmd_row$n_loss, base_row$n_loss),
        c(nmd_row$n_events, base_row$n_events)
      )

      results[[length(results) + 1]] <- tibble(
        nmd_comparison = nmd_comp,
        run = r,
        event_type = et,
        nmd_n = nmd_row$n_events,
        nmd_loss_prop = nmd_row$loss_prop,
        baseline_n = base_row$n_events,
        baseline_loss_prop = base_row$loss_prop,
        diff = nmd_row$loss_prop - base_row$loss_prop,
        p_value = pt$p.value
      )
    }
  }
}

results_df <- bind_rows(results)

# Add FDR correction within each comparison
results_df <- results_df %>%
  group_by(nmd_comparison) %>%
  mutate(fdr_p = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Print focused results for Alt_TES and SE
cat("\n\n═══ Alt_TES Direction (C1 signal) ═══\n\n")
results_df %>%
  filter(event_type == "Alt_TES") %>%
  arrange(nmd_comparison, run) %>%
  mutate(across(c(nmd_loss_prop, baseline_loss_prop, diff), ~round(., 3)),
         across(c(p_value, fdr_p), ~signif(., 3))) %>%
  print(n = 20, width = 120)

cat("\n\n═══ SE Direction (C2 signal) ═══\n\n")
results_df %>%
  filter(event_type == "SE") %>%
  arrange(nmd_comparison, run) %>%
  mutate(across(c(nmd_loss_prop, baseline_loss_prop, diff), ~round(., 3)),
         across(c(p_value, fdr_p), ~signif(., 3))) %>%
  print(n = 20, width = 120)

# Full table of all event types
cat("\n\n═══ All Event Types — Summary Across Runs ═══\n\n")
summary_by_event <- results_df %>%
  group_by(nmd_comparison, event_type) %>%
  summarise(
    n_runs = n(),
    mean_nmd_loss = mean(nmd_loss_prop),
    mean_baseline_loss = mean(baseline_loss_prop),
    mean_diff = mean(diff),
    n_significant = sum(fdr_p < 0.05),
    min_fdr = min(fdr_p),
    .groups = "drop"
  ) %>%
  arrange(nmd_comparison, event_type)

print(summary_by_event, n = 40, width = 120)

# Save
output_path <- file.path(base_dir, "cross_comparison", "direction_by_event_type.tsv")
write_tsv(results_df, output_path)
cat(sprintf("\n\nSaved: %s\n", output_path))

# Also save the summary
summary_path <- file.path(base_dir, "cross_comparison", "direction_by_event_type_summary.tsv")
write_tsv(summary_by_event, summary_path)
cat(sprintf("Saved: %s\n", summary_path))
