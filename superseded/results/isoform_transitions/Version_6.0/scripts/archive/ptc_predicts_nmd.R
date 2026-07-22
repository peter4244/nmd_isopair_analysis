#!/usr/bin/env Rscript
# PTC Predicts NMD Responsiveness
#
# Question: Among expressed isoforms passing the 5% filter, does PTC presence
# predict NMD responsiveness (adj.P.Val < 0.05, logFC > 0)?
#
# Sensitivity: Stratify by distance from stop codon to last EJC.
# Biology predicts stronger NMD effect for larger distances.

suppressPackageStartupMessages(library(tidyverse))

cat("\n")
cat("==================================================================\n")
cat("   PTC → NMD Responsiveness Analysis\n")
cat("==================================================================\n\n")

# ─── 1. Load data ────────────────────────────────────────────────────────────

cat("Loading data...\n")
cds_meta <- readRDS("data/isoform_cds_metadata.rds")
structs <- readRDS("data/isoform_structures.rds")
expression_data <- readRDS("data/expression_data_filtered.rds")
sample_metadata <- readRDS("data/sample_metadata.rds")

# Join treatment info
expression_data <- expression_data %>%
  left_join(sample_metadata %>% select(sample_id, treatment, ct), by = "sample_id")

# DE files
de_dir <- "/Users/petecastaldi/claude_projects/nmd/longread_dge"
de_date <- "2026.1.18"
ct_to_de_suffix <- c(
  "AT"     = "at2",
  "DD"     = "dd",
  "DD_ALI" = "ddali",
  "DO"     = "doali",
  "FB"     = "fb",
  "MV"     = "mv"
)

cell_types <- sort(unique(as.character(sample_metadata$ct)))

# ─── 2. Compute PTC status for all coding isoforms ──────────────────────────

cat("Computing PTC status...\n")

coding <- cds_meta %>%
  filter(coding_status == "coding") %>%
  inner_join(structs, by = "isoform_id")

cat(sprintf("  %d coding isoforms with exon structures\n", nrow(coding)))

compute_ptc_distance <- function(cds_start, cds_stop, strand,
                                  exon_starts, exon_ends, n_exons) {
  if (n_exons == 1) return(NA_real_)

  stop_pos <- if (strand == "+") cds_stop else cds_start

  if (strand == "+") {
    ord <- order(exon_starts)
  } else {
    ord <- order(exon_starts, decreasing = TRUE)
  }
  e_starts <- exon_starts[ord]
  e_ends <- exon_ends[ord]

  exon_lengths <- e_ends - e_starts
  cum_starts <- c(0, cumsum(exon_lengths[-length(exon_lengths)]))
  last_ejc_mRNA <- sum(exon_lengths[-length(exon_lengths)])

  if (strand == "+") {
    exon_idx <- which(e_starts <= stop_pos & stop_pos <= e_ends)
  } else {
    exon_idx <- which(e_starts <= stop_pos & stop_pos <= e_ends)
  }

  if (length(exon_idx) == 0) return(NA_real_)
  exon_idx <- exon_idx[1]

  if (strand == "+") {
    offset_in_exon <- stop_pos - e_starts[exon_idx]
  } else {
    offset_in_exon <- e_ends[exon_idx] - stop_pos
  }
  stop_mRNA <- cum_starts[exon_idx] + offset_in_exon

  last_ejc_mRNA - stop_mRNA
}

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
      n_exons == 1 ~ FALSE,
      is.na(ptc_distance) ~ NA,
      ptc_distance > 50 ~ TRUE,
      TRUE ~ FALSE
    ),
    ptc_distance_bin = case_when(
      n_exons == 1 ~ "single_exon",
      is.na(ptc_distance) ~ "unknown",
      ptc_distance <= 50 ~ "0-50 (no PTC)",
      ptc_distance <= 100 ~ "51-100",
      ptc_distance <= 200 ~ "101-200",
      ptc_distance <= 500 ~ "201-500",
      TRUE ~ ">500"
    )
  )

ptc_lookup <- ptc_results %>%
  select(isoform_id, gene_id, n_exons, ptc_distance, has_ptc, ptc_distance_bin)

cat(sprintf("  PTC (>50nt): %d\n", sum(ptc_lookup$has_ptc == TRUE, na.rm = TRUE)))
cat(sprintf("  Not PTC (<=50nt): %d\n",
            sum(ptc_lookup$has_ptc == FALSE & ptc_lookup$n_exons > 1, na.rm = TRUE)))
cat(sprintf("  Single-exon: %d\n", sum(ptc_lookup$n_exons == 1)))

# ─── 3. Apply 5% expression filter per cell type ────────────────────────────

cat("\nApplying 5% expression filter and classifying per cell type...\n")

apply_5pct_filter <- function(expr_df, sample_ids) {
  expr_subset <- expr_df %>% filter(sample_id %in% sample_ids)

  gene_totals <- expr_subset %>%
    group_by(gene_id, sample_id) %>%
    summarize(gene_cpm = sum(cpm), .groups = "drop")

  proportions <- expr_subset %>%
    left_join(gene_totals, by = c("gene_id", "sample_id")) %>%
    mutate(proportion = if_else(gene_cpm == 0, 0, cpm / gene_cpm))

  proportions %>%
    group_by(isoform_id) %>%
    summarize(max_prop = max(proportion), .groups = "drop") %>%
    filter(max_prop >= 0.05) %>%
    pull(isoform_id)
}

# ─── 4. Per-cell-type analysis ───────────────────────────────────────────────

all_results <- list()
all_stratified <- list()

for (ct in cell_types) {
  cat(sprintf("\n  --- %s ---\n", ct))

  # Load DE results
  de_file <- file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix[ct], "_", de_date, ".csv"))
  de <- read_csv(de_file, show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC) %>%
    rename(isoform_id = txid)

  # Get samples for this cell type
  ct_samples <- sample_metadata %>%
    filter(ct == !!ct) %>%
    pull(sample_id)

  # Apply 5% filter using all samples for this cell type
  passing_5pct <- apply_5pct_filter(expression_data, ct_samples)
  cat(sprintf("    Isoforms passing 5%% filter: %d\n", length(passing_5pct)))

  # Join PTC status with DE results, restricted to 5%-passing coding isoforms
  analysis_df <- ptc_lookup %>%
    filter(isoform_id %in% passing_5pct) %>%
    inner_join(de, by = "isoform_id") %>%
    filter(!is.na(has_ptc))  # exclude unknown PTC status

  cat(sprintf("    Coding isoforms with PTC status + DE results: %d\n", nrow(analysis_df)))

  # Classify NMD responsiveness
  analysis_df <- analysis_df %>%
    mutate(
      nmd_responsive = adj.P.Val < 0.05 & logFC > 0
    )

  # ── Cross-tabulation: PTC status × NMD responsiveness ──
  tab <- table(PTC = analysis_df$has_ptc, NMD_responsive = analysis_df$nmd_responsive)
  cat("    PTC × NMD responsiveness:\n")
  print(tab)

  if (all(dim(tab) == 2)) {
    ft <- fisher.test(tab)
    cat(sprintf("    Fisher OR = %.2f (%.2f–%.2f), p = %g\n",
                ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))

    nmd_rate_ptc <- tab["TRUE", "TRUE"] / sum(tab["TRUE", ])
    nmd_rate_no_ptc <- tab["FALSE", "TRUE"] / sum(tab["FALSE", ])
    cat(sprintf("    NMD rate with PTC: %.1f%% (%d/%d)\n",
                100 * nmd_rate_ptc, tab["TRUE", "TRUE"], sum(tab["TRUE", ])))
    cat(sprintf("    NMD rate without PTC: %.1f%% (%d/%d)\n",
                100 * nmd_rate_no_ptc, tab["FALSE", "TRUE"], sum(tab["FALSE", ])))

    all_results[[ct]] <- tibble(
      cell_type = ct,
      n_total = nrow(analysis_df),
      n_ptc = sum(analysis_df$has_ptc),
      n_no_ptc = sum(!analysis_df$has_ptc),
      n_nmd_responsive = sum(analysis_df$nmd_responsive),
      nmd_rate_ptc = nmd_rate_ptc,
      nmd_rate_no_ptc = nmd_rate_no_ptc,
      odds_ratio = ft$estimate,
      or_ci_lower = ft$conf.int[1],
      or_ci_upper = ft$conf.int[2],
      p_value = ft$p.value
    )
  }

  # ── Stratified by PTC distance bin ──
  cat("\n    NMD responsiveness by PTC distance bin:\n")

  strat <- analysis_df %>%
    filter(n_exons > 1, !is.na(ptc_distance)) %>%
    group_by(ptc_distance_bin) %>%
    summarise(
      n = n(),
      n_nmd = sum(nmd_responsive),
      nmd_rate = mean(nmd_responsive),
      median_distance = median(ptc_distance),
      .groups = "drop"
    ) %>%
    arrange(median_distance)

  print(strat, n = 10)

  strat$cell_type <- ct
  all_stratified[[ct]] <- strat
}

# ─── 5. Summary across cell types ────────────────────────────────────────────

cat("\n\n")
cat("==================================================================\n")
cat("   SUMMARY: PTC → NMD Responsiveness\n")
cat("==================================================================\n\n")

summary_df <- bind_rows(all_results) %>%
  mutate(fdr = p.adjust(p_value, method = "BH"))

cat("Per-cell-type results:\n\n")
summary_df %>%
  mutate(
    across(c(nmd_rate_ptc, nmd_rate_no_ptc), ~ sprintf("%.1f%%", 100 * .)),
    across(c(odds_ratio, or_ci_lower, or_ci_upper), ~ round(., 2)),
    p_value = signif(p_value, 3),
    fdr = signif(fdr, 3)
  ) %>%
  select(cell_type, n_total, n_ptc, nmd_rate_ptc, nmd_rate_no_ptc,
         odds_ratio, or_ci_lower, or_ci_upper, p_value, fdr) %>%
  print(n = 10, width = 140)

# ── Meta-analysis (fixed-effects, inverse-variance) ──
cat("\n\nFixed-effects meta-analysis (log-OR):\n")

meta_df <- bind_rows(all_results) %>%
  mutate(
    log_or = log(odds_ratio),
    se_log_or = (log(or_ci_upper) - log(or_ci_lower)) / (2 * 1.96)
  )

weights <- 1 / meta_df$se_log_or^2
pooled_log_or <- sum(weights * meta_df$log_or) / sum(weights)
pooled_se <- sqrt(1 / sum(weights))
pooled_or <- exp(pooled_log_or)
pooled_ci_lower <- exp(pooled_log_or - 1.96 * pooled_se)
pooled_ci_upper <- exp(pooled_log_or + 1.96 * pooled_se)
pooled_z <- pooled_log_or / pooled_se
pooled_p <- 2 * pnorm(-abs(pooled_z))

# Heterogeneity
Q <- sum(weights * (meta_df$log_or - pooled_log_or)^2)
k <- nrow(meta_df)
I2 <- max(0, (Q - (k - 1)) / Q) * 100
Q_p <- pchisq(Q, df = k - 1, lower.tail = FALSE)

cat(sprintf("  Pooled OR = %.2f (%.2f–%.2f), p = %g\n",
            pooled_or, pooled_ci_lower, pooled_ci_upper, pooled_p))
cat(sprintf("  Heterogeneity: Q = %.1f (p = %g), I² = %.0f%%\n", Q, Q_p, I2))

# ─── 6. Distance-stratified summary ─────────────────────────────────────────

cat("\n\n")
cat("==================================================================\n")
cat("   NMD RESPONSIVENESS BY PTC DISTANCE (all cell types pooled)\n")
cat("==================================================================\n\n")

strat_all <- bind_rows(all_stratified)

# Pool across cell types per distance bin
strat_pooled <- strat_all %>%
  group_by(ptc_distance_bin) %>%
  summarise(
    n_cell_types = n(),
    total_n = sum(n),
    total_nmd = sum(n_nmd),
    pooled_nmd_rate = sum(n_nmd) / sum(n),
    mean_nmd_rate = mean(nmd_rate),
    sd_nmd_rate = sd(nmd_rate),
    median_distance = median(median_distance),
    .groups = "drop"
  ) %>%
  arrange(median_distance)

cat("Pooled across cell types:\n")
strat_pooled %>%
  mutate(
    pooled_nmd_rate = sprintf("%.1f%%", 100 * pooled_nmd_rate),
    mean_nmd_rate = sprintf("%.1f%%", 100 * mean_nmd_rate)
  ) %>%
  print(n = 10, width = 120)

# ── Trend test: does NMD rate increase with PTC distance? ──
cat("\n\nTrend test: NMD rate vs PTC distance (Spearman correlation on isoform-level data)\n")

for (ct in cell_types) {
  de_file <- file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix[ct], "_", de_date, ".csv"))
  de <- read_csv(de_file, show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC) %>%
    rename(isoform_id = txid)

  ct_samples <- sample_metadata %>%
    filter(ct == !!ct) %>%
    pull(sample_id)
  passing_5pct <- apply_5pct_filter(expression_data, ct_samples)

  trend_df <- ptc_lookup %>%
    filter(isoform_id %in% passing_5pct, n_exons > 1, !is.na(ptc_distance)) %>%
    inner_join(de, by = "isoform_id") %>%
    mutate(nmd_responsive = as.integer(adj.P.Val < 0.05 & logFC > 0))

  if (nrow(trend_df) > 50) {
    st <- cor.test(trend_df$ptc_distance, trend_df$nmd_responsive, method = "spearman")
    cat(sprintf("  %s: rho = %.4f, p = %g (n = %d)\n",
                ct, st$estimate, st$p.value, nrow(trend_df)))
  }
}

# ── Logistic regression: PTC distance predicting NMD responsiveness ──
cat("\n\nLogistic regression: log(ptc_distance) → NMD responsiveness\n")

for (ct in cell_types) {
  de_file <- file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix[ct], "_", de_date, ".csv"))
  de <- read_csv(de_file, show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC) %>%
    rename(isoform_id = txid)

  ct_samples <- sample_metadata %>%
    filter(ct == !!ct) %>%
    pull(sample_id)
  passing_5pct <- apply_5pct_filter(expression_data, ct_samples)

  glm_df <- ptc_lookup %>%
    filter(isoform_id %in% passing_5pct, n_exons > 1, !is.na(ptc_distance),
           ptc_distance > 0) %>%
    inner_join(de, by = "isoform_id") %>%
    mutate(
      nmd_responsive = as.integer(adj.P.Val < 0.05 & logFC > 0),
      log_ptc_distance = log(ptc_distance)
    )

  if (nrow(glm_df) > 50 && sum(glm_df$nmd_responsive) > 10) {
    fit <- glm(nmd_responsive ~ log_ptc_distance, data = glm_df, family = binomial)
    s <- summary(fit)
    coef_row <- s$coefficients["log_ptc_distance", ]
    or <- exp(coef_row["Estimate"])
    or_ci <- exp(coef_row["Estimate"] + c(-1.96, 1.96) * coef_row["Std. Error"])
    cat(sprintf("  %s: OR per log-unit = %.3f (%.3f–%.3f), p = %g (n = %d, n_nmd = %d)\n",
                ct, or, or_ci[1], or_ci[2], coef_row["Pr(>|z|)"],
                nrow(glm_df), sum(glm_df$nmd_responsive)))
  }
}

# ─── 7. Save outputs ────────────────────────────────────────────────────────

output_dir <- "comparisons/nonNMD_0.50/cross_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv(summary_df, file.path(output_dir, "ptc_predicts_nmd_by_celltype.tsv"))
cat(sprintf("\n\nSaved: %s\n", file.path(output_dir, "ptc_predicts_nmd_by_celltype.tsv")))

write_tsv(strat_pooled, file.path(output_dir, "ptc_distance_nmd_rate.tsv"))
cat(sprintf("Saved: %s\n", file.path(output_dir, "ptc_distance_nmd_rate.tsv")))

# Per-cell-type stratified results
strat_all_out <- strat_all %>%
  select(cell_type, ptc_distance_bin, n, n_nmd, nmd_rate, median_distance)
write_tsv(strat_all_out, file.path(output_dir, "ptc_distance_nmd_rate_by_celltype.tsv"))
cat(sprintf("Saved: %s\n", file.path(output_dir, "ptc_distance_nmd_rate_by_celltype.tsv")))

cat("\nDone.\n")
