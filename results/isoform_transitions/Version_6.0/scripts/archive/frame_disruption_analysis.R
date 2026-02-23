#!/usr/bin/env Rscript
# Frame Disruption Analysis
#
# Question: Do NMD-triggering transitions more often disrupt the dominant
# isoform's reading frame than baseline transitions?
#
# Approach: For each splicing event between a pair, compute the overlap
# with the dominant isoform's CDS. If the CDS-affecting portion of an
# SE, Missing_Internal, A5SS, or A3SS event is not divisible by 3, the
# event introduces a frameshift in the comparator → likely PTC downstream.
#
# This avoids the self-referential problem of using each isoform's own
# CDS annotation — we always use the DOMINANT's reading frame as reference.

suppressPackageStartupMessages(library(tidyverse))

cat("Loading data...\n")

cds_meta <- readRDS("data/isoform_cds_metadata.rds") %>%
  filter(coding_status == "coding")
profiles <- readRDS("comparisons/nonNMD_0.50/deduplicated/splicing_choice_profiles.rds")

cat(sprintf("  %d coding isoforms with CDS annotations\n", nrow(cds_meta)))
cat(sprintf("  %d splicing profiles\n", nrow(profiles)))

# Load comparison pairs to map profiles to comparisons
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
pair_map <- bind_rows(all_pairs) %>%
  distinct(comparison, run, gene_id, dominant_isoform_id, non_dominant_isoform_id)

# ─── Unnest events and annotate with dominant CDS ───────────────────────────

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

cat(sprintf("  %d total events across all profiles\n", nrow(events)))

# Join dominant CDS coordinates
events <- events %>%
  left_join(
    cds_meta %>% select(isoform_id, cds_start, cds_stop, orf_length),
    by = c("dominant_isoform_id" = "isoform_id")
  )

cat(sprintf("  %d events with dominant CDS annotation (%.1f%%)\n",
            sum(!is.na(events$cds_start)),
            100 * mean(!is.na(events$cds_start))))

# Filter to events with CDS annotation on dominant
events_cds <- events %>% filter(!is.na(cds_start))

# ─── Compute CDS overlap for each event ────────────────────────────────────

# Convert five_prime/three_prime to genomic coordinates (they are strand-aware)
events_cds <- events_cds %>%
  mutate(
    event_start = pmin(five_prime, three_prime),
    event_end = pmax(five_prime, three_prime)
  )

# Compute overlap with CDS region
# CDS overlap = max(0, min(event_end, cds_stop) - max(event_start, cds_start))
events_cds <- events_cds %>%
  mutate(
    overlap_start = pmax(event_start, cds_start),
    overlap_end = pmin(event_end, cds_stop),
    cds_overlap_bp = pmax(0, overlap_end - overlap_start),
    affects_cds = cds_overlap_bp > 0,
    # Frame disruption: CDS overlap not divisible by 3
    frame_disrupting = affects_cds & (cds_overlap_bp %% 3 != 0)
  )

# ─── Classify events by frame impact ───────────────────────────────────────

# Focus on events that can cause frameshifts:
# - SE, Missing_Internal: exon inclusion/exclusion
# - A5SS, A3SS: splice site shifts
# - IR events typically introduce stop codons directly (intronic sequence)

frameshift_events <- c("SE", "Missing_Internal", "A5SS", "A3SS",
                        "Partial_IR_5", "Partial_IR_3")

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("FRAME DISRUPTION ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Summary by event type
cat("\n─── CDS overlap and frame disruption by event type ───\n\n")
event_summary <- events_cds %>%
  group_by(event_type) %>%
  summarise(
    n_events = n(),
    n_affects_cds = sum(affects_cds),
    pct_affects_cds = round(100 * mean(affects_cds), 1),
    n_frame_disrupting = sum(frame_disrupting),
    pct_frame_disrupting_of_cds = round(
      100 * sum(frame_disrupting) / max(sum(affects_cds), 1), 1
    ),
    median_cds_overlap = median(cds_overlap_bp[affects_cds]),
    .groups = "drop"
  ) %>%
  arrange(desc(n_events))

print(event_summary, n = 15, width = 120)

# ─── Profile-level frame disruption ────────────────────────────────────────

cat("\n─── Profile-level: does any event disrupt the reading frame? ───\n\n")

# For each profile, determine if ANY event disrupts the CDS reading frame
profile_frame <- events_cds %>%
  filter(event_type %in% frameshift_events) %>%
  group_by(gene_id, dominant_isoform_id, non_dominant_isoform_id) %>%
  summarise(
    n_cds_events = sum(affects_cds),
    n_frame_disrupting = sum(frame_disrupting),
    has_frame_disruption = any(frame_disrupting),
    # Also track if any event is frame-preserving but CDS-affecting
    has_cds_preserving = any(affects_cds & !frame_disrupting),
    .groups = "drop"
  )

# Join back to all profiles (including those without CDS annotation)
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

cat(sprintf("Profiles with dominant CDS: %d / %d (%.1f%%)\n",
            sum(profiles_annotated$dom_has_cds),
            nrow(profiles_annotated),
            100 * mean(profiles_annotated$dom_has_cds)))

# Among profiles with CDS-annotated dominants
cds_profiles <- profiles_annotated %>% filter(dom_has_cds)

cat(sprintf("Of those, profiles with frame-disrupting events: %d (%.1f%%)\n",
            sum(cds_profiles$has_frame_disruption),
            100 * mean(cds_profiles$has_frame_disruption)))

# ─── Compare NMD vs baseline ───────────────────────────────────────────────

cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("NMD vs BASELINE: FRAME DISRUPTION RATES\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Map profiles to comparisons
profile_comparison <- profiles_annotated %>%
  inner_join(pair_map, by = c("gene_id", "dominant_isoform_id",
                               "non_dominant_isoform_id"))

cat(sprintf("\nProfile-comparison mappings: %d\n", nrow(profile_comparison)))

# Per comparison x run
cat("\n─── Frame disruption rate by comparison × run ───\n\n")

frame_by_run <- profile_comparison %>%
  filter(dom_has_cds) %>%
  group_by(comparison, run) %>%
  summarise(
    n_pairs = n(),
    n_frame_disrupted = sum(has_frame_disruption),
    frame_disruption_rate = mean(has_frame_disruption),
    .groups = "drop"
  )

frame_by_run %>%
  mutate(frame_disruption_rate = round(frame_disruption_rate, 3)) %>%
  pivot_wider(names_from = comparison, values_from = c(n_pairs, frame_disruption_rate)) %>%
  print(width = 120)

# Aggregate by comparison (deduplicated)
cat("\n─── Aggregate by comparison (deduplicated) ───\n\n")

dedup_comp <- profile_comparison %>%
  filter(dom_has_cds) %>%
  distinct(comparison, gene_id, dominant_isoform_id, non_dominant_isoform_id,
           .keep_all = TRUE)

agg <- dedup_comp %>%
  group_by(comparison) %>%
  summarise(
    n_pairs = n(),
    n_frame_disrupted = sum(has_frame_disruption),
    frame_disruption_rate = round(mean(has_frame_disruption), 3),
    .groups = "drop"
  )
print(agg)

# ─── Statistical tests ─────────────────────────────────────────────────────

cat("\n─── Statistical tests: NMD vs C4 baseline ───\n\n")

results <- list()

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

    results[[length(results) + 1]] <- tibble(
      nmd_comparison = nmd_comp,
      run = r,
      nmd_n = nrow(nmd_data),
      nmd_fd_n = nmd_fd,
      nmd_fd_rate = nmd_fd / nrow(nmd_data),
      baseline_n = nrow(base_data),
      baseline_fd_n = base_fd,
      baseline_fd_rate = base_fd / nrow(base_data),
      diff = nmd_fd / nrow(nmd_data) - base_fd / nrow(base_data),
      odds_ratio = ft$estimate,
      or_ci_lower = ft$conf.int[1],
      or_ci_upper = ft$conf.int[2],
      p_value = ft$p.value
    )
  }
}

results_df <- bind_rows(results) %>%
  mutate(fdr_p = p.adjust(p_value, method = "BH"))

results_df %>%
  mutate(across(c(nmd_fd_rate, baseline_fd_rate, diff), ~round(., 3)),
         across(c(odds_ratio, or_ci_lower, or_ci_upper), ~round(., 2)),
         across(c(p_value, fdr_p), ~signif(., 3))) %>%
  print(n = 20, width = 140)

# ─── Breakdown by event type: which events cause frame disruption? ──────────

cat("\n\n─── Frame-disrupting events by type: NMD vs baseline ───\n\n")

# Per event type, compare frame disruption rate
event_level <- events_cds %>%
  filter(event_type %in% frameshift_events, affects_cds) %>%
  inner_join(
    pair_map,
    by = c("gene_id", "dominant_isoform_id", "non_dominant_isoform_id")
  ) %>%
  distinct(comparison, gene_id, dominant_isoform_id, non_dominant_isoform_id,
           event_type, direction, frame_disrupting, cds_overlap_bp)

event_fd_rates <- event_level %>%
  group_by(comparison, event_type) %>%
  summarise(
    n_cds_events = n(),
    n_frame_disrupting = sum(frame_disrupting),
    fd_rate = round(mean(frame_disrupting), 3),
    median_cds_bp = median(cds_overlap_bp),
    .groups = "drop"
  )

event_fd_rates %>%
  pivot_wider(names_from = comparison, values_from = c(n_cds_events, fd_rate)) %>%
  print(width = 120)

# ─── Direction analysis: LOSS events that disrupt frame ─────────────────────

cat("\n─── Frame disruption by direction (LOSS = comparator missing sequence) ───\n\n")

# LOSS events are where the comparator is MISSING sequence relative to dominant
# These are the most likely to cause PTCs in the comparator
direction_fd <- event_level %>%
  group_by(comparison, direction) %>%
  summarise(
    n_events = n(),
    n_fd = sum(frame_disrupting),
    fd_rate = round(mean(frame_disrupting), 3),
    .groups = "drop"
  )

print(direction_fd, n = 10)

# Specifically LOSS + frame-disrupting = comparator missing a non-3n CDS chunk
cat("\n─── LOSS + frame-disrupting events (comparator has frameshift → PTC) ───\n\n")

loss_fd <- event_level %>%
  filter(direction == "LOSS", frame_disrupting) %>%
  group_by(comparison, event_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = comparison, values_from = n, values_fill = 0)

print(loss_fd, width = 100)

# Profile-level: any LOSS frame-disrupting event?
cat("\n─── Profiles with LOSS frame-disrupting events (most likely PTC) ───\n\n")

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

loss_fd_agg <- profiles_loss_fd %>%
  distinct(comparison, gene_id, dominant_isoform_id, non_dominant_isoform_id,
           .keep_all = TRUE) %>%
  group_by(comparison) %>%
  summarise(
    n_pairs = n(),
    n_loss_fd = sum(has_loss_fd),
    loss_fd_rate = round(mean(has_loss_fd), 3),
    .groups = "drop"
  )

print(loss_fd_agg)

# Test LOSS frame-disruption rates
cat("\n─── Tests for LOSS frame-disruption: NMD vs C4 ───\n\n")

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
      nmd_n = nrow(nmd_data),
      nmd_rate = mean(nmd_data$has_loss_fd),
      baseline_n = nrow(base_data),
      baseline_rate = mean(base_data$has_loss_fd),
      diff = mean(nmd_data$has_loss_fd) - mean(base_data$has_loss_fd),
      odds_ratio = ft$estimate,
      p_value = ft$p.value
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

# ─── Save results ──────────────────────────────────────────────────────────

output_dir <- file.path(base_dir, "cross_comparison")

write_tsv(results_df, file.path(output_dir, "frame_disruption_nmd_vs_baseline.tsv"))
write_tsv(loss_results_df, file.path(output_dir, "frame_disruption_loss_nmd_vs_baseline.tsv"))
cat(sprintf("\nSaved to: %s\n", output_dir))

cat("\nDone.\n")
