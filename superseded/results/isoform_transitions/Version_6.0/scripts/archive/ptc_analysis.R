#!/usr/bin/env Rscript
# PTC (Premature Termination Codon) Analysis
# Question: Do NMD-triggering transitions involve isoforms with PTCs
# more often than baseline transitions?
#
# PTC definition: stop codon >50nt upstream of the last exon-exon junction
# (the canonical NMD "50-nucleotide rule")

suppressPackageStartupMessages(library(tidyverse))

cat("Loading data...\n")
cds_meta <- readRDS("data/isoform_cds_metadata.rds")
structs <- readRDS("data/isoform_structures.rds")

# ─── Compute PTC status for all coding isoforms ─────────────────────────────

cat("Computing PTC status for all coding isoforms...\n")

coding <- cds_meta %>%
  filter(coding_status == "coding") %>%
  inner_join(structs, by = "isoform_id")

cat(sprintf("  %d coding isoforms with exon structures\n", nrow(coding)))

# For each isoform, compute:
# 1. Stop codon position in mRNA coordinates
# 2. Last EJC position in mRNA coordinates
# 3. Distance = last_EJC - stop_codon_mRNA
#
# Coordinate conventions:
#   cds_start/cds_stop: always min/max genomic coordinate
#   + strand: stop codon biological position = cds_stop
#   - strand: stop codon biological position = cds_start

compute_ptc_distance <- function(cds_start, cds_stop, strand,
                                  exon_starts, exon_ends, n_exons) {
  if (n_exons == 1) return(NA_real_)  # No EJCs in single-exon transcripts

  # Determine biological stop codon position (genomic)
  stop_pos <- if (strand == "+") cds_stop else cds_start

  # Order exons biologically (5' to 3')
  if (strand == "+") {
    ord <- order(exon_starts)
  } else {
    ord <- order(exon_starts, decreasing = TRUE)
  }
  e_starts <- exon_starts[ord]
  e_ends <- exon_ends[ord]

  # Compute exon lengths and cumulative mRNA positions
  exon_lengths <- e_ends - e_starts  # 0-based lengths (end is exclusive in BED, but GFF/GRanges uses inclusive)
  # Actually for GRanges-style coordinates (1-based inclusive): length = end - start + 1
  # But these are from rtracklayer which uses 1-based start, so length = end - start + 1
  # Let's check: the first exon for our + strand example was 2794970-2795244
  # That's 2795244 - 2794970 + 1 = 275 bp. But GFF/BED conventions vary.
  # Using the convention that start is 0-based and end is 1-based (BED-like):
  # length = end - start. This seems to match the data (e.g. exon 2794970-2795244 = 274bp).
  # Let's use end - start (works correctly for distance calculations regardless of off-by-one)
  exon_lengths <- e_ends - e_starts

  # Cumulative start positions in mRNA space (0-based)
  cum_starts <- c(0, cumsum(exon_lengths[-length(exon_lengths)]))

  # Last EJC position in mRNA = cumulative length through penultimate exon
  last_ejc_mRNA <- sum(exon_lengths[-length(exon_lengths)])

  # Find which exon contains the stop codon
  if (strand == "+") {
    # Stop codon in exon where start <= stop_pos <= end
    exon_idx <- which(e_starts <= stop_pos & stop_pos <= e_ends)
  } else {
    exon_idx <- which(e_starts <= stop_pos & stop_pos <= e_ends)
  }

  if (length(exon_idx) == 0) return(NA_real_)  # Stop codon not in any exon
  exon_idx <- exon_idx[1]  # Take first match if multiple

  # Stop codon position in mRNA coordinates
  if (strand == "+") {
    offset_in_exon <- stop_pos - e_starts[exon_idx]
  } else {
    offset_in_exon <- e_ends[exon_idx] - stop_pos
  }
  stop_mRNA <- cum_starts[exon_idx] + offset_in_exon

  # PTC distance = last_EJC - stop_codon in mRNA space
  ptc_distance <- last_ejc_mRNA - stop_mRNA

  return(ptc_distance)
}

# Apply to all coding isoforms
ptc_results <- coding %>%
  rowwise() %>%
  mutate(
    ptc_distance = compute_ptc_distance(
      cds_start, cds_stop, as.character(strand),
      exon_starts, exon_ends, n_exons
    )
  ) %>%
  ungroup() %>%
  mutate(
    has_ptc = case_when(
      n_exons == 1 ~ FALSE,     # Single-exon: no EJCs possible
      is.na(ptc_distance) ~ NA, # Stop codon not found in exons
      ptc_distance > 50 ~ TRUE, # PTC: >50nt from last EJC
      TRUE ~ FALSE              # Not PTC: stop within 50nt of last EJC
    )
  )

cat(sprintf("  Single-exon (no EJCs): %d\n",
            sum(ptc_results$n_exons == 1, na.rm = TRUE)))
cat(sprintf("  Stop not in exon (NA): %d\n",
            sum(is.na(ptc_results$ptc_distance) & ptc_results$n_exons > 1)))
cat(sprintf("  PTC (>50nt): %d\n", sum(ptc_results$has_ptc == TRUE, na.rm = TRUE)))
cat(sprintf("  Not PTC (<=50nt): %d\n",
            sum(ptc_results$has_ptc == FALSE & ptc_results$n_exons > 1, na.rm = TRUE)))

# Quick sanity: distribution of PTC distances
cat("\nPTC distance distribution (multi-exon coding isoforms):\n")
multi_exon <- ptc_results %>% filter(n_exons > 1, !is.na(ptc_distance))
cat(sprintf("  n = %d\n", nrow(multi_exon)))
cat(sprintf("  Median: %.0f\n", median(multi_exon$ptc_distance)))
cat(sprintf("  Mean: %.0f\n", mean(multi_exon$ptc_distance)))
q <- quantile(multi_exon$ptc_distance, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1))
cat("  Quantiles:\n")
print(q)

# Trim to just what we need for the comparison
ptc_lookup <- ptc_results %>%
  select(isoform_id, n_exons, ptc_distance, has_ptc, orf_length)

# ─── Load comparison pairs and annotate with PTC status ─────────────────────

cat("\n\nLoading comparison pairs...\n")

base_dir <- "comparisons/nonNMD_0.50"
comparisons <- c("C1", "C2", "C4")
runs <- c("all_samples", "at", "dd", "dd_ali", "fb", "mv")

all_pairs <- list()
for (comp in comparisons) {
  for (run in runs) {
    pairs_path <- file.path(base_dir, comp, run, "pairs.tsv")
    if (!file.exists(pairs_path)) next
    pairs <- read_tsv(pairs_path, show_col_types = FALSE) %>%
      mutate(comparison = comp, run = run) %>%
      rename(non_dominant_isoform_id = comparator_isoform_id)
    all_pairs[[paste(comp, run, sep = "_")]] <- pairs
  }
}
all_pairs <- bind_rows(all_pairs)
cat(sprintf("  Total pair-run combinations: %d\n", nrow(all_pairs)))

# Annotate both dominant and comparator with PTC status
annotated <- all_pairs %>%
  left_join(ptc_lookup %>% rename_with(~paste0("dom_", .), -isoform_id),
            by = c("dominant_isoform_id" = "isoform_id")) %>%
  left_join(ptc_lookup %>% rename_with(~paste0("comp_", .), -isoform_id),
            by = c("non_dominant_isoform_id" = "isoform_id"))

# Focus on pairs where both isoforms are coding with known PTC status
annotated_clean <- annotated %>%
  filter(!is.na(dom_has_ptc), !is.na(comp_has_ptc))

cat(sprintf("  Pairs with PTC status for both isoforms: %d / %d\n",
            nrow(annotated_clean), nrow(annotated)))

# ─── Analysis: PTC prevalence by comparison ─────────────────────────────────

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("PTC PREVALENCE IN COMPARATOR ISOFORMS: NMD vs BASELINE\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Deduplicate: each unique pair counted once per comparison x run
# (pairs may appear across multiple runs)

cat("\n─── Comparator PTC rates by comparison × run ───\n\n")

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

# Aggregate by comparison (pool across runs, deduplicate pairs first)
cat("\n─── Aggregate by comparison (deduplicated pairs) ───\n\n")

dedup_pairs <- annotated_clean %>%
  distinct(comparison, dominant_isoform_id, non_dominant_isoform_id,
           .keep_all = TRUE)

comp_ptc_agg <- dedup_pairs %>%
  group_by(comparison) %>%
  summarise(
    n_pairs = n(),
    n_comp_ptc = sum(comp_has_ptc),
    comp_ptc_rate = mean(comp_has_ptc),
    n_dom_ptc = sum(dom_has_ptc),
    dom_ptc_rate = mean(dom_has_ptc),
    .groups = "drop"
  )

print(comp_ptc_agg)

# ─── Statistical tests: NMD comparators vs baseline comparators ─────────────

cat("\n─── Statistical tests ───\n\n")

# Per-run tests: C1 vs C4, C2 vs C4
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

    # Two-proportion test
    ft <- fisher.test(matrix(c(nmd_ptc, nmd_n - nmd_ptc,
                                base_ptc, base_n - base_ptc), nrow = 2))

    results[[length(results) + 1]] <- tibble(
      nmd_comparison = nmd_comp,
      run = r,
      nmd_n = nmd_n,
      nmd_ptc_n = nmd_ptc,
      nmd_ptc_rate = nmd_ptc / nmd_n,
      baseline_n = base_n,
      baseline_ptc_n = base_ptc,
      baseline_ptc_rate = base_ptc / base_n,
      diff = nmd_ptc / nmd_n - base_ptc / base_n,
      odds_ratio = ft$estimate,
      or_ci_lower = ft$conf.int[1],
      or_ci_upper = ft$conf.int[2],
      p_value = ft$p.value
    )
  }
}

results_df <- bind_rows(results) %>%
  mutate(fdr_p = p.adjust(p_value, method = "BH"))

cat("Comparator PTC rates: NMD vs Baseline (C4)\n\n")
results_df %>%
  mutate(across(c(nmd_ptc_rate, baseline_ptc_rate, diff), ~round(., 3)),
         across(c(odds_ratio, or_ci_lower, or_ci_upper), ~round(., 2)),
         across(c(p_value, fdr_p), ~signif(., 3))) %>%
  print(n = 20, width = 140)

# ─── Also check dominant isoform PTC rates ──────────────────────────────────

cat("\n\n─── Dominant isoform PTC rates (sanity check) ───\n\n")

dom_results <- list()
for (nmd_comp in c("C1", "C2")) {
  for (r in runs) {
    nmd_data <- annotated_clean %>% filter(comparison == nmd_comp, run == r)
    base_data <- annotated_clean %>% filter(comparison == "C4", run == r)
    if (nrow(nmd_data) < 10 || nrow(base_data) < 10) next

    dom_results[[length(dom_results) + 1]] <- tibble(
      nmd_comparison = nmd_comp,
      run = r,
      nmd_dom_ptc_rate = mean(nmd_data$dom_has_ptc),
      baseline_dom_ptc_rate = mean(base_data$dom_has_ptc)
    )
  }
}
dom_df <- bind_rows(dom_results)
dom_df %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print(n = 20)

# ─── PTC distance distribution comparison ───────────────────────────────────

cat("\n\n─── PTC distance distributions (comparator isoforms) ───\n\n")

# Compare PTC distance distributions for PTC+ isoforms
for (nmd_comp in c("C1", "C2")) {
  nmd_ptc_dists <- dedup_pairs %>%
    filter(comparison == nmd_comp, comp_has_ptc == TRUE) %>%
    pull(comp_ptc_distance)

  base_ptc_dists <- dedup_pairs %>%
    filter(comparison == "C4", comp_has_ptc == TRUE) %>%
    pull(comp_ptc_distance)

  if (length(nmd_ptc_dists) > 5 && length(base_ptc_dists) > 5) {
    cat(sprintf("%s PTC+ comparators: median distance = %.0f (n = %d)\n",
                nmd_comp, median(nmd_ptc_dists), length(nmd_ptc_dists)))
    cat(sprintf("C4 PTC+ comparators: median distance = %.0f (n = %d)\n",
                median(base_ptc_dists), length(base_ptc_dists)))
    wt <- wilcox.test(nmd_ptc_dists, base_ptc_dists)
    cat(sprintf("  Wilcoxon p = %g\n\n", wt$p.value))
  }
}

# ─── Pair-level analysis: is the comparator more likely to have PTC than dom? ─

cat("\n─── Within-pair PTC asymmetry ───\n\n")

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

# McNemar's test: is there asymmetry in PTC status within pairs?
cat("\nMcNemar's test for PTC asymmetry within pairs:\n")
for (comp in c("C1", "C2", "C4")) {
  d <- dedup_pairs %>% filter(comparison == comp)
  b <- sum(d$comp_has_ptc & !d$dom_has_ptc)  # comp PTC, dom not
  cc <- sum(!d$comp_has_ptc & d$dom_has_ptc)   # dom PTC, comp not
  if (b + cc > 0) {
    mt <- mcnemar.test(matrix(c(
      sum(d$comp_has_ptc & d$dom_has_ptc),
      b, cc,
      sum(!d$comp_has_ptc & !d$dom_has_ptc)
    ), nrow = 2))
    cat(sprintf("  %s: comp_only=%d, dom_only=%d, chi2=%.2f, p=%g\n",
                comp, b, cc, mt$statistic, mt$p.value))
  }
}

# ─── Save results ───────────────────────────────────────────────────────────

output_dir <- file.path(base_dir, "cross_comparison")

write_tsv(results_df, file.path(output_dir, "ptc_prevalence_nmd_vs_baseline.tsv"))
cat(sprintf("\n\nSaved: %s\n", file.path(output_dir, "ptc_prevalence_nmd_vs_baseline.tsv")))

write_tsv(comp_ptc_by_run, file.path(output_dir, "ptc_rates_by_comparison_run.tsv"))
cat(sprintf("Saved: %s\n", file.path(output_dir, "ptc_rates_by_comparison_run.tsv")))

write_tsv(pair_asym, file.path(output_dir, "ptc_within_pair_asymmetry.tsv"))
cat(sprintf("Saved: %s\n", file.path(output_dir, "ptc_within_pair_asymmetry.tsv")))

cat("\nDone.\n")
