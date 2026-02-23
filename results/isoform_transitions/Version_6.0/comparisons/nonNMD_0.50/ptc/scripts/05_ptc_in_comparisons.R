#!/usr/bin/env Rscript
# 05_ptc_in_comparisons.R
#
# PTC prevalence and frame disruption analysis in C1/C2/C4 comparison pairs.
#
# Part A: PTC prevalence in comparator isoforms — NMD vs baseline
# Part B: Frame disruption — events that shift the dominant's reading frame
#
# Requires: ptc_status.rds from 01_compute_ptc_status.R
#
# Outputs:
#   ptc_prevalence_nmd_vs_baseline.tsv
#   ptc_rates_by_comparison_run.tsv
#   ptc_within_pair_asymmetry.tsv
#   frame_disruption_nmd_vs_baseline.tsv
#   frame_disruption_loss_nmd_vs_baseline.tsv
#
# Run from: results/isoform_transitions/Version_6.0/

suppressPackageStartupMessages(library(tidyverse))

cat("\n")
cat("==================================================================\n")
cat("   05: PTC in Comparison Pairs + Frame Disruption\n")
cat("==================================================================\n\n")

# ─── 1. Load data ────────────────────────────────────────────────────────────

cat("Loading data...\n")
ptc_status <- readRDS("comparisons/nonNMD_0.50/ptc/results/ptc_status.rds")
cds_meta   <- readRDS("data/isoform_cds_metadata.rds") %>%
  filter(coding_status == "coding")
profiles   <- readRDS("comparisons/nonNMD_0.50/deduplicated/splicing_choice_profiles.rds")

ptc_lookup <- ptc_status %>%
  filter(!is.na(has_ptc)) %>%
  select(isoform_id, n_exons, ptc_distance, has_ptc, orf_length)

# Load comparison pairs
base_dir <- "comparisons/nonNMD_0.50"
comparisons <- c("C1", "C2", "C4")
runs <- c("all_samples", "at", "dd", "dd_ali", "fb", "mv")

all_pairs <- list()
for (comp in comparisons) {
  for (run in runs) {
    path <- file.path(base_dir, comp, run, "pairs.tsv")
    if (!file.exists(path)) next
    all_pairs[[paste(comp, run)]] <- read_tsv(path, show_col_types = FALSE) %>%
      mutate(comparison = comp, run = run) %>%
      rename(non_dominant_isoform_id = comparator_isoform_id)
  }
}
pair_map <- bind_rows(all_pairs)
cat(sprintf("  Total pair-run combinations: %d\n", nrow(pair_map)))

# ═══════════════════════════════════════════════════════════════════════════════
# PART A: PTC Prevalence in Comparison Pairs
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("==================================================================\n")
cat("   PART A: PTC Prevalence in Comparator Isoforms\n")
cat("==================================================================\n\n")

# Annotate both dominant and comparator with PTC status
annotated <- pair_map %>%
  left_join(ptc_lookup %>% rename_with(~paste0("dom_", .), -isoform_id),
            by = c("dominant_isoform_id" = "isoform_id")) %>%
  left_join(ptc_lookup %>% rename_with(~paste0("comp_", .), -isoform_id),
            by = c("non_dominant_isoform_id" = "isoform_id"))

annotated_clean <- annotated %>%
  filter(!is.na(dom_has_ptc), !is.na(comp_has_ptc))

cat(sprintf("  Pairs with PTC status for both: %d / %d\n",
            nrow(annotated_clean), nrow(annotated)))

# ─── Comparator PTC rates by comparison x run ────────────────────────────────

cat("\n--- Comparator PTC rates by comparison x run ---\n\n")

comp_ptc_by_run <- annotated_clean %>%
  group_by(comparison, run) %>%
  summarise(
    n_pairs = n(),
    n_comp_ptc = sum(comp_has_ptc),
    comp_ptc_rate = mean(comp_has_ptc),
    n_dom_ptc = sum(dom_has_ptc),
    dom_ptc_rate = mean(dom_has_ptc),
    .groups = "drop"
  )

print(comp_ptc_by_run, n = 30)

# ─── Statistical tests: NMD vs C4 ────────────────────────────────────────────

cat("\n--- Statistical tests: NMD comparators vs C4 baseline ---\n\n")

results <- list()
for (nmd_comp in c("C1", "C2")) {
  for (r in runs) {
    nmd_data <- annotated_clean %>% filter(comparison == nmd_comp, run == r)
    base_data <- annotated_clean %>% filter(comparison == "C4", run == r)
    if (nrow(nmd_data) < 10 || nrow(base_data) < 10) next

    nmd_ptc <- sum(nmd_data$comp_has_ptc)
    nmd_n <- nrow(nmd_data)
    base_ptc <- sum(base_data$comp_has_ptc)
    base_n <- nrow(base_data)

    ft <- fisher.test(matrix(c(nmd_ptc, nmd_n - nmd_ptc,
                                base_ptc, base_n - base_ptc), nrow = 2))
    results[[length(results) + 1]] <- tibble(
      nmd_comparison = nmd_comp, run = r,
      nmd_n = nmd_n, nmd_ptc_n = nmd_ptc,
      nmd_ptc_rate = nmd_ptc / nmd_n,
      baseline_n = base_n, baseline_ptc_n = base_ptc,
      baseline_ptc_rate = base_ptc / base_n,
      diff = nmd_ptc / nmd_n - base_ptc / base_n,
      odds_ratio = ft$estimate,
      or_ci_lower = ft$conf.int[1], or_ci_upper = ft$conf.int[2],
      p_value = ft$p.value
    )
  }
}

results_df <- bind_rows(results) %>%
  mutate(fdr_p = p.adjust(p_value, method = "BH"))

results_df %>%
  mutate(across(c(nmd_ptc_rate, baseline_ptc_rate, diff), ~round(., 3)),
         across(c(odds_ratio, or_ci_lower, or_ci_upper), ~round(., 2)),
         across(c(p_value, fdr_p), ~signif(., 3))) %>%
  print(n = 20, width = 140)

# ─── Within-pair PTC asymmetry ───────────────────────────────────────────────

cat("\n--- Within-pair PTC asymmetry ---\n\n")

dedup_pairs <- annotated_clean %>%
  distinct(comparison, dominant_isoform_id, non_dominant_isoform_id,
           .keep_all = TRUE)

pair_asym <- dedup_pairs %>%
  mutate(
    comp_ptc_not_dom = comp_has_ptc & !dom_has_ptc,
    dom_ptc_not_comp = dom_has_ptc & !comp_has_ptc,
    both_ptc = comp_has_ptc & dom_has_ptc,
    neither_ptc = !comp_has_ptc & !dom_has_ptc
  ) %>%
  group_by(comparison) %>%
  summarise(
    n = n(),
    comp_only_ptc = sum(comp_ptc_not_dom),
    dom_only_ptc = sum(dom_ptc_not_comp),
    both_ptc = sum(both_ptc),
    neither_ptc = sum(neither_ptc),
    pct_comp_only = 100 * sum(comp_ptc_not_dom) / n(),
    pct_dom_only = 100 * sum(dom_ptc_not_comp) / n(),
    .groups = "drop"
  )

print(pair_asym)

cat("\nMcNemar's test for PTC asymmetry within pairs:\n")
for (comp in c("C1", "C2", "C4")) {
  d <- dedup_pairs %>% filter(comparison == comp)
  b <- sum(d$comp_has_ptc & !d$dom_has_ptc)
  cc <- sum(!d$comp_has_ptc & d$dom_has_ptc)
  if (b + cc > 0) {
    mt <- mcnemar.test(matrix(c(
      sum(d$comp_has_ptc & d$dom_has_ptc), b, cc,
      sum(!d$comp_has_ptc & !d$dom_has_ptc)
    ), nrow = 2))
    cat(sprintf("  %s: comp_only=%d, dom_only=%d, chi2=%.2f, p=%g\n",
                comp, b, cc, mt$statistic, mt$p.value))
  }
}

# ─── PTC distance distribution comparison ────────────────────────────────────

cat("\n--- PTC distance distributions (comparator isoforms) ---\n\n")

for (nmd_comp in c("C1", "C2")) {
  nmd_dists <- dedup_pairs %>%
    filter(comparison == nmd_comp, comp_has_ptc == TRUE) %>%
    pull(comp_ptc_distance)
  base_dists <- dedup_pairs %>%
    filter(comparison == "C4", comp_has_ptc == TRUE) %>%
    pull(comp_ptc_distance)

  if (length(nmd_dists) > 5 && length(base_dists) > 5) {
    cat(sprintf("%s PTC+: median=%.0f (n=%d)\n", nmd_comp, median(nmd_dists), length(nmd_dists)))
    cat(sprintf("C4 PTC+: median=%.0f (n=%d)\n", median(base_dists), length(base_dists)))
    wt <- wilcox.test(nmd_dists, base_dists)
    cat(sprintf("  Wilcoxon p = %g\n\n", wt$p.value))
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# PART B: Frame Disruption Analysis
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("==================================================================\n")
cat("   PART B: Frame Disruption Analysis\n")
cat("==================================================================\n\n")

# Unnest events and annotate with dominant CDS
cat("Unnesting events and computing CDS overlap...\n")

events <- profiles %>%
  select(gene_id, dominant_isoform_id, non_dominant_isoform_id,
         detailed_events) %>%
  unnest(detailed_events, names_sep = "_") %>%
  select(
    gene_id, dominant_isoform_id, non_dominant_isoform_id,
    event_type = detailed_events_event_type,
    direction = detailed_events_direction,
    five_prime = detailed_events_five_prime,
    three_prime = detailed_events_three_prime,
    bp_diff = detailed_events_bp_diff
  )

cat(sprintf("  %d total events\n", nrow(events)))

# Join dominant CDS coordinates
events <- events %>%
  left_join(
    cds_meta %>% select(isoform_id, cds_start, cds_stop, orf_length),
    by = c("dominant_isoform_id" = "isoform_id")
  )

events_cds <- events %>%
  filter(!is.na(cds_start)) %>%
  mutate(
    event_start = pmin(five_prime, three_prime),
    event_end = pmax(five_prime, three_prime),
    overlap_start = pmax(event_start, cds_start),
    overlap_end = pmin(event_end, cds_stop),
    cds_overlap_bp = pmax(0, overlap_end - overlap_start),
    affects_cds = cds_overlap_bp > 0,
    frame_disrupting = affects_cds & (cds_overlap_bp %% 3 != 0)
  )

# Event types that can cause frameshifts
frameshift_events <- c("SE", "Missing_Internal", "A5SS", "A3SS",
                        "Partial_IR_5", "Partial_IR_3")

# ─── Event-level summary ─────────────────────────────────────────────────────

cat("\n--- CDS overlap and frame disruption by event type ---\n\n")

event_summary <- events_cds %>%
  group_by(event_type) %>%
  summarise(
    n_events = n(),
    n_affects_cds = sum(affects_cds),
    pct_affects_cds = round(100 * mean(affects_cds), 1),
    n_frame_disrupting = sum(frame_disrupting),
    pct_frame_disrupting_of_cds = round(
      100 * sum(frame_disrupting) / max(sum(affects_cds), 1), 1),
    median_cds_overlap = median(cds_overlap_bp[affects_cds]),
    .groups = "drop"
  ) %>%
  arrange(desc(n_events))

print(event_summary, n = 15, width = 120)

# ─── Profile-level frame disruption ─────────────────────────────────────────

cat("\n--- Profile-level frame disruption ---\n\n")

profile_frame <- events_cds %>%
  filter(event_type %in% frameshift_events) %>%
  group_by(gene_id, dominant_isoform_id, non_dominant_isoform_id) %>%
  summarise(
    n_cds_events = sum(affects_cds),
    n_frame_disrupting = sum(frame_disrupting),
    has_frame_disruption = any(frame_disrupting),
    .groups = "drop"
  )

profiles_annotated <- profiles %>%
  select(gene_id, dominant_isoform_id, non_dominant_isoform_id) %>%
  left_join(
    cds_meta %>% select(isoform_id, cds_start),
    by = c("dominant_isoform_id" = "isoform_id")
  ) %>%
  mutate(dom_has_cds = !is.na(cds_start)) %>%
  select(-cds_start) %>%
  left_join(profile_frame, by = c("gene_id", "dominant_isoform_id",
                                   "non_dominant_isoform_id")) %>%
  mutate(
    has_frame_disruption = replace_na(has_frame_disruption, FALSE),
    n_frame_disrupting = replace_na(n_frame_disrupting, 0)
  )

cds_profiles <- profiles_annotated %>% filter(dom_has_cds)
cat(sprintf("Profiles with dominant CDS: %d / %d (%.1f%%)\n",
            sum(profiles_annotated$dom_has_cds), nrow(profiles_annotated),
            100 * mean(profiles_annotated$dom_has_cds)))
cat(sprintf("Profiles with frame disruption: %d (%.1f%%)\n",
            sum(cds_profiles$has_frame_disruption),
            100 * mean(cds_profiles$has_frame_disruption)))

# ─── NMD vs baseline frame disruption rates ──────────────────────────────────

cat("\n--- NMD vs baseline: frame disruption rates ---\n\n")

# Map profiles to comparisons (use distinct pair_map for cleaner join)
pair_map_distinct <- pair_map %>%
  distinct(comparison, run, gene_id, dominant_isoform_id, non_dominant_isoform_id)

profile_comparison <- profiles_annotated %>%
  inner_join(pair_map_distinct, by = c("gene_id", "dominant_isoform_id",
                                        "non_dominant_isoform_id"))

frame_results <- list()
for (nmd_comp in c("C1", "C2")) {
  for (r in runs) {
    nmd_data <- profile_comparison %>%
      filter(comparison == nmd_comp, run == r, dom_has_cds)
    base_data <- profile_comparison %>%
      filter(comparison == "C4", run == r, dom_has_cds)
    if (nrow(nmd_data) < 10 || nrow(base_data) < 10) next

    nmd_fd <- sum(nmd_data$has_frame_disruption)
    base_fd <- sum(base_data$has_frame_disruption)

    ft <- fisher.test(matrix(c(
      nmd_fd, nrow(nmd_data) - nmd_fd,
      base_fd, nrow(base_data) - base_fd
    ), nrow = 2))

    frame_results[[length(frame_results) + 1]] <- tibble(
      nmd_comparison = nmd_comp, run = r,
      nmd_n = nrow(nmd_data), nmd_fd_n = nmd_fd,
      nmd_fd_rate = nmd_fd / nrow(nmd_data),
      baseline_n = nrow(base_data), baseline_fd_n = base_fd,
      baseline_fd_rate = base_fd / nrow(base_data),
      diff = nmd_fd / nrow(nmd_data) - base_fd / nrow(base_data),
      odds_ratio = ft$estimate,
      or_ci_lower = ft$conf.int[1], or_ci_upper = ft$conf.int[2],
      p_value = ft$p.value
    )
  }
}

frame_results_df <- bind_rows(frame_results) %>%
  mutate(fdr_p = p.adjust(p_value, method = "BH"))

frame_results_df %>%
  mutate(across(c(nmd_fd_rate, baseline_fd_rate, diff), ~round(., 3)),
         across(c(odds_ratio, or_ci_lower, or_ci_upper), ~round(., 2)),
         across(c(p_value, fdr_p), ~signif(., 3))) %>%
  print(n = 20, width = 140)

# ─── LOSS + frame-disrupting events (most likely PTC) ─────────────────────────

cat("\n--- LOSS + frame-disrupting events (comparator has frameshift -> PTC) ---\n\n")

loss_fd_profiles <- events_cds %>%
  filter(event_type %in% frameshift_events,
         direction == "LOSS",
         frame_disrupting) %>%
  distinct(gene_id, dominant_isoform_id, non_dominant_isoform_id) %>%
  mutate(has_loss_fd = TRUE)

profiles_loss_fd <- profile_comparison %>%
  filter(dom_has_cds) %>%
  left_join(loss_fd_profiles,
            by = c("gene_id", "dominant_isoform_id", "non_dominant_isoform_id")) %>%
  mutate(has_loss_fd = replace_na(has_loss_fd, FALSE))

# Aggregate
loss_fd_agg <- profiles_loss_fd %>%
  distinct(comparison, gene_id, dominant_isoform_id, non_dominant_isoform_id,
           .keep_all = TRUE) %>%
  group_by(comparison) %>%
  summarise(
    n_pairs = n(), n_loss_fd = sum(has_loss_fd),
    loss_fd_rate = round(mean(has_loss_fd), 3),
    .groups = "drop"
  )
print(loss_fd_agg)

# Test LOSS frame-disruption rates
cat("\n--- Tests for LOSS frame-disruption: NMD vs C4 ---\n\n")

loss_results <- list()
for (nmd_comp in c("C1", "C2")) {
  for (r in runs) {
    nmd_data <- profiles_loss_fd %>%
      filter(comparison == nmd_comp, run == r)
    base_data <- profiles_loss_fd %>%
      filter(comparison == "C4", run == r)
    if (nrow(nmd_data) < 10 || nrow(base_data) < 10) next

    ft <- fisher.test(matrix(c(
      sum(nmd_data$has_loss_fd), sum(!nmd_data$has_loss_fd),
      sum(base_data$has_loss_fd), sum(!base_data$has_loss_fd)
    ), nrow = 2))

    loss_results[[length(loss_results) + 1]] <- tibble(
      nmd_comparison = nmd_comp, run = r,
      nmd_n = nrow(nmd_data), nmd_rate = mean(nmd_data$has_loss_fd),
      baseline_n = nrow(base_data), baseline_rate = mean(base_data$has_loss_fd),
      diff = mean(nmd_data$has_loss_fd) - mean(base_data$has_loss_fd),
      odds_ratio = ft$estimate, p_value = ft$p.value
    )
  }
}

loss_results_df <- bind_rows(loss_results) %>%
  mutate(fdr_p = p.adjust(p_value, method = "BH"))

loss_results_df %>%
  mutate(across(c(nmd_rate, baseline_rate, diff), ~round(., 3)),
         odds_ratio = round(odds_ratio, 2),
         across(c(p_value, fdr_p), ~signif(., 3))) %>%
  print(n = 20, width = 120)

# ─── Save all outputs ────────────────────────────────────────────────────────

output_dir <- "comparisons/nonNMD_0.50/ptc/results"

write_tsv(results_df, file.path(output_dir, "ptc_prevalence_nmd_vs_baseline.tsv"))
write_tsv(comp_ptc_by_run, file.path(output_dir, "ptc_rates_by_comparison_run.tsv"))
write_tsv(pair_asym, file.path(output_dir, "ptc_within_pair_asymmetry.tsv"))
write_tsv(frame_results_df, file.path(output_dir, "frame_disruption_nmd_vs_baseline.tsv"))
write_tsv(loss_results_df, file.path(output_dir, "frame_disruption_loss_nmd_vs_baseline.tsv"))

cat(sprintf("\nSaved 5 output files to: %s\n", output_dir))
cat("\nDone.\n")
