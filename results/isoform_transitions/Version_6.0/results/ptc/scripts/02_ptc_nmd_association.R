#!/usr/bin/env Rscript
# 02_ptc_nmd_association.R
#
# Tests whether PTC presence and/or GENCODE NMD biotype predicts
# NMD responsiveness (adj.P.Val < 0.05, logFC > 0).
#
# Analyses:
#   A) PTC → NMD responsiveness (5%-filtered + unfiltered)
#   B) GENCODE NMD biotype → NMD responsiveness (5%-filtered + ENST-only)
#   C) GENCODE vs PacBio isoform stratification
#   D) PTC distance stratification + trend tests
#   E) Logistic regression + meta-analysis
#
# Requires: ptc_status.rds from 01_compute_ptc_status.R
#
# Outputs:
#   ptc_predicts_nmd_by_celltype.tsv
#   ptc_predicts_nmd_unfiltered_by_celltype.tsv
#   ptc_predicts_nmd_by_source.tsv         (GENCODE vs PacBio)
#   ptc_distance_nmd_rate.tsv
#   ptc_distance_nmd_rate_by_celltype.tsv
#   nmd_biotype_predicts_nmd_by_celltype.tsv
#   nmd_biotype_predicts_nmd_enst_only.tsv
#   nmd_rate_by_biotype.tsv
#
# Run from: results/isoform_transitions/Version_6.0/

suppressPackageStartupMessages(library(tidyverse))

cat("\n")
cat("==================================================================\n")
cat("   02: PTC / NMD Biotype → NMD Responsiveness\n")
cat("==================================================================\n\n")

# ─── 1. Load data ────────────────────────────────────────────────────────────

cat("Loading data...\n")
ptc_status      <- readRDS("results/ptc/results/ptc_status.rds")
expression_data <- readRDS("data/expression_data_filtered.rds")
sample_metadata <- readRDS("data/sample_metadata.rds")

expression_data <- expression_data %>%
  left_join(sample_metadata %>% select(sample_id, treatment, ct), by = "sample_id")

de_dir <- "/Users/petecastaldi/claude_projects/nmd/longread_dge"
de_date <- "2026.1.18"
ct_to_de_suffix <- c(
  "AT" = "at2", "DD" = "dd", "DD_ALI" = "ddali",
  "DO" = "doali", "FB" = "fb", "MV" = "mv"
)
cell_types <- sort(unique(as.character(sample_metadata$ct)))

# ─── 2. 5% expression filter ────────────────────────────────────────────────

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

# ─── Helper: fixed-effects meta-analysis ─────────────────────────────────────

run_meta_analysis <- function(results_list, label = "") {
  meta_df <- bind_rows(results_list) %>%
    filter(!is.na(odds_ratio), odds_ratio > 0,
           !is.na(or_ci_lower), or_ci_lower > 0,
           !is.na(or_ci_upper), or_ci_upper > 0) %>%
    mutate(
      log_or = log(odds_ratio),
      se_log_or = (log(or_ci_upper) - log(or_ci_lower)) / (2 * 1.96)
    ) %>%
    filter(is.finite(se_log_or), se_log_or > 0)

  if (nrow(meta_df) < 2) {
    cat(sprintf("  %s: too few cell types for meta-analysis\n", label))
    return(invisible(NULL))
  }

  weights <- 1 / meta_df$se_log_or^2
  pooled_log_or <- sum(weights * meta_df$log_or) / sum(weights)
  pooled_se <- sqrt(1 / sum(weights))
  pooled_or <- exp(pooled_log_or)
  pooled_ci <- exp(pooled_log_or + c(-1.96, 1.96) * pooled_se)
  pooled_z <- pooled_log_or / pooled_se
  pooled_p <- 2 * pnorm(-abs(pooled_z))

  Q <- sum(weights * (meta_df$log_or - pooled_log_or)^2)
  k <- nrow(meta_df)
  I2 <- max(0, (Q - (k - 1)) / Q) * 100
  Q_p <- pchisq(Q, df = k - 1, lower.tail = FALSE)

  cat(sprintf("  %s: Pooled OR = %.2f (%.2f-%.2f), p = %g\n",
              label, pooled_or, pooled_ci[1], pooled_ci[2], pooled_p))
  cat(sprintf("    Heterogeneity: Q = %.1f (p = %g), I^2 = %.0f%%\n", Q, Q_p, I2))
}

# ═══════════════════════════════════════════════════════════════════════════════
# PART A: PTC → NMD Responsiveness
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("==================================================================\n")
cat("   PART A: PTC → NMD Responsiveness (5%-filtered)\n")
cat("==================================================================\n\n")

ptc_lookup <- ptc_status %>% filter(!is.na(has_ptc))

all_results_filt <- list()
all_results_unfilt <- list()
all_stratified <- list()
all_results_gencode <- list()
all_results_pacbio <- list()

for (ct in cell_types) {
  cat(sprintf("  --- %s ---\n", ct))

  # Load DE results
  de_file <- file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix[ct], "_", de_date, ".csv"))
  de <- read_csv(de_file, show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC) %>%
    rename(isoform_id = txid)

  ct_samples <- sample_metadata %>% filter(ct == !!ct) %>% pull(sample_id)
  passing_5pct <- apply_5pct_filter(expression_data, ct_samples)

  # ── Filtered analysis ──
  analysis_filt <- ptc_lookup %>%
    filter(isoform_id %in% passing_5pct) %>%
    inner_join(de, by = "isoform_id") %>%
    mutate(nmd_responsive = adj.P.Val < 0.05 & logFC > 0)

  tab <- table(PTC = analysis_filt$has_ptc, NMD = analysis_filt$nmd_responsive)
  if (all(dim(tab) == 2)) {
    ft <- fisher.test(tab)
    nmd_rate_ptc <- tab["TRUE", "TRUE"] / sum(tab["TRUE", ])
    nmd_rate_no_ptc <- tab["FALSE", "TRUE"] / sum(tab["FALSE", ])
    cat(sprintf("    [5%% filt] OR = %.2f (%.2f-%.2f), p = %g\n",
                ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
    all_results_filt[[ct]] <- tibble(
      cell_type = ct, filter = "5pct",
      n_total = nrow(analysis_filt), n_ptc = sum(analysis_filt$has_ptc),
      n_nmd_responsive = sum(analysis_filt$nmd_responsive),
      nmd_rate_ptc = nmd_rate_ptc, nmd_rate_no_ptc = nmd_rate_no_ptc,
      odds_ratio = ft$estimate,
      or_ci_lower = ft$conf.int[1], or_ci_upper = ft$conf.int[2],
      p_value = ft$p.value
    )
  }

  # ── Unfiltered analysis ──
  analysis_unfilt <- ptc_lookup %>%
    inner_join(de, by = "isoform_id") %>%
    mutate(nmd_responsive = adj.P.Val < 0.05 & logFC > 0)

  tab_u <- table(PTC = analysis_unfilt$has_ptc, NMD = analysis_unfilt$nmd_responsive)
  if (all(dim(tab_u) == 2)) {
    ft_u <- fisher.test(tab_u)
    all_results_unfilt[[ct]] <- tibble(
      cell_type = ct, filter = "none",
      n_total = nrow(analysis_unfilt), n_ptc = sum(analysis_unfilt$has_ptc),
      n_nmd_responsive = sum(analysis_unfilt$nmd_responsive),
      nmd_rate_ptc = tab_u["TRUE", "TRUE"] / sum(tab_u["TRUE", ]),
      nmd_rate_no_ptc = tab_u["FALSE", "TRUE"] / sum(tab_u["FALSE", ]),
      odds_ratio = ft_u$estimate,
      or_ci_lower = ft_u$conf.int[1], or_ci_upper = ft_u$conf.int[2],
      p_value = ft_u$p.value
    )
  }

  # ── GENCODE-only analysis (5% filtered) ──
  analysis_gc <- analysis_filt %>% filter(is_gencode)
  tab_gc <- table(PTC = analysis_gc$has_ptc, NMD = analysis_gc$nmd_responsive)
  if (all(dim(tab_gc) == 2)) {
    ft_gc <- fisher.test(tab_gc)
    all_results_gencode[[ct]] <- tibble(
      cell_type = ct, source = "GENCODE", filter = "5pct",
      n_total = nrow(analysis_gc), n_ptc = sum(analysis_gc$has_ptc),
      n_nmd_responsive = sum(analysis_gc$nmd_responsive),
      nmd_rate_ptc = tab_gc["TRUE", "TRUE"] / sum(tab_gc["TRUE", ]),
      nmd_rate_no_ptc = tab_gc["FALSE", "TRUE"] / sum(tab_gc["FALSE", ]),
      odds_ratio = ft_gc$estimate,
      or_ci_lower = ft_gc$conf.int[1], or_ci_upper = ft_gc$conf.int[2],
      p_value = ft_gc$p.value
    )
    cat(sprintf("    [GENCODE] OR = %.2f (%.2f-%.2f), p = %g (n=%d)\n",
                ft_gc$estimate, ft_gc$conf.int[1], ft_gc$conf.int[2],
                ft_gc$p.value, nrow(analysis_gc)))
  }

  # ── PacBio-only analysis (5% filtered) ──
  analysis_pb <- analysis_filt %>% filter(is_novel)
  tab_pb <- table(PTC = analysis_pb$has_ptc, NMD = analysis_pb$nmd_responsive)
  if (all(dim(tab_pb) == 2)) {
    ft_pb <- fisher.test(tab_pb)
    all_results_pacbio[[ct]] <- tibble(
      cell_type = ct, source = "PacBio", filter = "5pct",
      n_total = nrow(analysis_pb), n_ptc = sum(analysis_pb$has_ptc),
      n_nmd_responsive = sum(analysis_pb$nmd_responsive),
      nmd_rate_ptc = tab_pb["TRUE", "TRUE"] / sum(tab_pb["TRUE", ]),
      nmd_rate_no_ptc = tab_pb["FALSE", "TRUE"] / sum(tab_pb["FALSE", ]),
      odds_ratio = ft_pb$estimate,
      or_ci_lower = ft_pb$conf.int[1], or_ci_upper = ft_pb$conf.int[2],
      p_value = ft_pb$p.value
    )
    cat(sprintf("    [PacBio]  OR = %.2f (%.2f-%.2f), p = %g (n=%d)\n",
                ft_pb$estimate, ft_pb$conf.int[1], ft_pb$conf.int[2],
                ft_pb$p.value, nrow(analysis_pb)))
  }

  # ── Distance-stratified analysis ──
  strat <- analysis_filt %>%
    filter(n_exons > 1, !is.na(ptc_distance)) %>%
    group_by(ptc_distance_bin) %>%
    summarise(
      n = n(), n_nmd = sum(nmd_responsive),
      nmd_rate = mean(nmd_responsive),
      median_distance = median(ptc_distance),
      .groups = "drop"
    ) %>%
    arrange(median_distance) %>%
    mutate(cell_type = ct)

  all_stratified[[ct]] <- strat
}

# ─── Summary: 5%-filtered ────────────────────────────────────────────────────

cat("\n\n")
cat("==================================================================\n")
cat("   SUMMARY: PTC → NMD Responsiveness\n")
cat("==================================================================\n\n")

summary_filt <- bind_rows(all_results_filt) %>%
  mutate(fdr = p.adjust(p_value, method = "BH"))

cat("5%-filtered results:\n\n")
summary_filt %>%
  mutate(
    across(c(nmd_rate_ptc, nmd_rate_no_ptc), ~ sprintf("%.1f%%", 100 * .)),
    across(c(odds_ratio, or_ci_lower, or_ci_upper), ~ round(., 2)),
    p_value = signif(p_value, 3), fdr = signif(fdr, 3)
  ) %>%
  select(cell_type, n_total, n_ptc, nmd_rate_ptc, nmd_rate_no_ptc,
         odds_ratio, or_ci_lower, or_ci_upper, p_value, fdr) %>%
  print(n = 10, width = 140)

cat("\nMeta-analysis (5%-filtered):\n")
run_meta_analysis(all_results_filt, "PTC 5%-filtered")

# ─── Summary: unfiltered ─────────────────────────────────────────────────────

summary_unfilt <- bind_rows(all_results_unfilt) %>%
  mutate(fdr = p.adjust(p_value, method = "BH"))

cat("\n\nUnfiltered results:\n\n")
summary_unfilt %>%
  mutate(
    across(c(nmd_rate_ptc, nmd_rate_no_ptc), ~ sprintf("%.1f%%", 100 * .)),
    across(c(odds_ratio, or_ci_lower, or_ci_upper), ~ round(., 2)),
    p_value = signif(p_value, 3), fdr = signif(fdr, 3)
  ) %>%
  select(cell_type, n_total, n_ptc, nmd_rate_ptc, nmd_rate_no_ptc,
         odds_ratio, or_ci_lower, or_ci_upper, p_value, fdr) %>%
  print(n = 10, width = 140)

cat("\nMeta-analysis (unfiltered):\n")
run_meta_analysis(all_results_unfilt, "PTC unfiltered")

# ─── Summary: GENCODE vs PacBio ──────────────────────────────────────────────

cat("\n\n")
cat("==================================================================\n")
cat("   GENCODE vs PacBio: PTC → NMD Responsiveness\n")
cat("==================================================================\n\n")

source_results <- bind_rows(
  bind_rows(all_results_gencode),
  bind_rows(all_results_pacbio)
) %>%
  group_by(source) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup()

source_results %>%
  mutate(
    across(c(nmd_rate_ptc, nmd_rate_no_ptc), ~ sprintf("%.1f%%", 100 * .)),
    across(c(odds_ratio, or_ci_lower, or_ci_upper), ~ round(., 2)),
    p_value = signif(p_value, 3), fdr = signif(fdr, 3)
  ) %>%
  select(source, cell_type, n_total, n_ptc, nmd_rate_ptc, nmd_rate_no_ptc,
         odds_ratio, or_ci_lower, or_ci_upper, p_value, fdr) %>%
  print(n = 20, width = 150)

cat("\nMeta-analysis by source:\n")
run_meta_analysis(all_results_gencode, "GENCODE only")
run_meta_analysis(all_results_pacbio, "PacBio only")

# ─── Distance-stratified summary ─────────────────────────────────────────────

cat("\n\n")
cat("==================================================================\n")
cat("   NMD Rate by PTC Distance (pooled)\n")
cat("==================================================================\n\n")

strat_all <- bind_rows(all_stratified)

strat_pooled <- strat_all %>%
  group_by(ptc_distance_bin) %>%
  summarise(
    n_cell_types = n(),
    total_n = sum(n), total_nmd = sum(n_nmd),
    pooled_nmd_rate = sum(n_nmd) / sum(n),
    mean_nmd_rate = mean(nmd_rate),
    sd_nmd_rate = sd(nmd_rate),
    median_distance = median(median_distance),
    .groups = "drop"
  ) %>%
  arrange(median_distance)

strat_pooled %>%
  mutate(pooled_nmd_rate = sprintf("%.1f%%", 100 * pooled_nmd_rate)) %>%
  print(n = 10, width = 120)

# ─── Trend tests ─────────────────────────────────────────────────────────────

cat("\nSpearman correlation: PTC distance vs NMD responsiveness\n")
for (ct in cell_types) {
  de_file <- file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix[ct], "_", de_date, ".csv"))
  de <- read_csv(de_file, show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC) %>%
    rename(isoform_id = txid)
  ct_samples <- sample_metadata %>% filter(ct == !!ct) %>% pull(sample_id)
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

# ─── Logistic regression ─────────────────────────────────────────────────────

cat("\nLogistic regression: log(ptc_distance) -> NMD responsiveness\n")
for (ct in cell_types) {
  de_file <- file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix[ct], "_", de_date, ".csv"))
  de <- read_csv(de_file, show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC) %>%
    rename(isoform_id = txid)
  ct_samples <- sample_metadata %>% filter(ct == !!ct) %>% pull(sample_id)
  passing_5pct <- apply_5pct_filter(expression_data, ct_samples)

  glm_df <- ptc_lookup %>%
    filter(isoform_id %in% passing_5pct, n_exons > 1, !is.na(ptc_distance),
           ptc_distance > 0) %>%
    inner_join(de, by = "isoform_id") %>%
    mutate(nmd_responsive = as.integer(adj.P.Val < 0.05 & logFC > 0),
           log_ptc_distance = log(ptc_distance))

  if (nrow(glm_df) > 50 && sum(glm_df$nmd_responsive) > 10) {
    fit <- glm(nmd_responsive ~ log_ptc_distance, data = glm_df, family = binomial)
    s <- summary(fit)
    coef_row <- s$coefficients["log_ptc_distance", ]
    or <- exp(coef_row["Estimate"])
    or_ci <- exp(coef_row["Estimate"] + c(-1.96, 1.96) * coef_row["Std. Error"])
    cat(sprintf("  %s: OR per log-unit = %.3f (%.3f-%.3f), p = %g (n = %d)\n",
                ct, or, or_ci[1], or_ci[2], coef_row["Pr(>|z|)"], nrow(glm_df)))
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# PART B: GENCODE NMD Biotype → NMD Responsiveness
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\n")
cat("==================================================================\n")
cat("   PART B: GENCODE NMD Biotype → NMD Responsiveness\n")
cat("==================================================================\n\n")

all_biotype_results <- list()
all_biotype_tables <- list()
enst_results <- list()

for (ct in cell_types) {
  cat(sprintf("  --- %s ---\n", ct))

  de_file <- file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix[ct], "_", de_date, ".csv"))
  de <- read_csv(de_file, show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC, biotype) %>%
    rename(isoform_id = txid)

  ct_samples <- sample_metadata %>% filter(ct == !!ct) %>% pull(sample_id)
  passing_5pct <- apply_5pct_filter(expression_data, ct_samples)

  # ── 5%-filtered, all isoforms ──
  analysis_df <- de %>%
    filter(isoform_id %in% passing_5pct, !is.na(biotype)) %>%
    mutate(nmd_biotype = biotype == "nonsense_mediated_decay",
           nmd_responsive = adj.P.Val < 0.05 & logFC > 0)

  tab <- table(NMD_biotype = analysis_df$nmd_biotype,
               NMD_responsive = analysis_df$nmd_responsive)

  if (all(dim(tab) == 2)) {
    ft <- fisher.test(tab)
    nmd_rate_biotype <- tab["TRUE", "TRUE"] / sum(tab["TRUE", ])
    nmd_rate_other <- tab["FALSE", "TRUE"] / sum(tab["FALSE", ])

    cat(sprintf("    OR = %.2f (%.2f-%.2f), p = %g\n",
                ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))

    all_biotype_results[[ct]] <- tibble(
      cell_type = ct,
      n_total = nrow(analysis_df), n_nmd_biotype = sum(analysis_df$nmd_biotype),
      n_nmd_responsive = sum(analysis_df$nmd_responsive),
      nmd_rate_biotype = nmd_rate_biotype, nmd_rate_other = nmd_rate_other,
      rate_ratio = nmd_rate_biotype / nmd_rate_other,
      odds_ratio = ft$estimate,
      or_ci_lower = ft$conf.int[1], or_ci_upper = ft$conf.int[2],
      p_value = ft$p.value
    )
  }

  # ── Biotype breakdown ──
  bt_table <- analysis_df %>%
    group_by(biotype) %>%
    summarise(n = n(), n_nmd = sum(nmd_responsive),
              nmd_rate = mean(nmd_responsive), .groups = "drop") %>%
    arrange(desc(n_nmd)) %>%
    mutate(cell_type = ct)
  all_biotype_tables[[ct]] <- bt_table

  # ── ENST-only sensitivity analysis ──
  analysis_enst <- de %>%
    filter(isoform_id %in% passing_5pct,
           startsWith(isoform_id, "ENST"),
           !is.na(biotype)) %>%
    mutate(nmd_biotype = biotype == "nonsense_mediated_decay",
           nmd_responsive = adj.P.Val < 0.05 & logFC > 0)

  tab_e <- table(NMD_biotype = analysis_enst$nmd_biotype,
                  NMD_responsive = analysis_enst$nmd_responsive)
  if (all(dim(tab_e) == 2)) {
    ft_e <- fisher.test(tab_e)
    enst_results[[ct]] <- tibble(
      cell_type = ct,
      n_total = nrow(analysis_enst), n_nmd_biotype = sum(analysis_enst$nmd_biotype),
      nmd_rate_biotype = tab_e["TRUE", "TRUE"] / sum(tab_e["TRUE", ]),
      nmd_rate_other = tab_e["FALSE", "TRUE"] / sum(tab_e["FALSE", ]),
      odds_ratio = ft_e$estimate,
      or_ci_lower = ft_e$conf.int[1], or_ci_upper = ft_e$conf.int[2],
      p_value = ft_e$p.value
    )
  }
}

# ─── Biotype summary ─────────────────────────────────────────────────────────

biotype_summary <- bind_rows(all_biotype_results) %>%
  mutate(fdr = p.adjust(p_value, method = "BH"))

cat("\n\nNMD biotype results (5%-filtered):\n\n")
biotype_summary %>%
  mutate(
    across(c(nmd_rate_biotype, nmd_rate_other), ~ sprintf("%.1f%%", 100 * .)),
    rate_ratio = round(rate_ratio, 2),
    across(c(odds_ratio, or_ci_lower, or_ci_upper), ~ round(., 2)),
    p_value = signif(p_value, 3), fdr = signif(fdr, 3)
  ) %>%
  select(cell_type, n_total, n_nmd_biotype, nmd_rate_biotype, nmd_rate_other,
         rate_ratio, odds_ratio, or_ci_lower, or_ci_upper, p_value, fdr) %>%
  print(n = 10, width = 140)

cat("\nMeta-analysis (NMD biotype):\n")
run_meta_analysis(all_biotype_results, "NMD biotype")

# Biotype pooled table
bt_all <- bind_rows(all_biotype_tables)
bt_pooled <- bt_all %>%
  group_by(biotype) %>%
  summarise(
    n_cell_types = n(), total_n = sum(n), total_nmd = sum(n_nmd),
    pooled_nmd_rate = sum(n_nmd) / sum(n),
    mean_nmd_rate = mean(nmd_rate), .groups = "drop"
  ) %>%
  filter(total_n >= 100) %>%
  arrange(desc(pooled_nmd_rate))

cat("\nNMD rate by biotype (pooled, n >= 100):\n")
bt_pooled %>%
  mutate(pooled_nmd_rate = sprintf("%.2f%%", 100 * pooled_nmd_rate)) %>%
  print(n = 30, width = 120)

# ENST-only meta
cat("\nENST-only meta-analysis:\n")
run_meta_analysis(enst_results, "ENST-only NMD biotype")

# ─── Save all outputs ────────────────────────────────────────────────────────

output_dir <- "results/ptc/results"

write_tsv(summary_filt, file.path(output_dir, "ptc_predicts_nmd_by_celltype.tsv"))
write_tsv(summary_unfilt, file.path(output_dir, "ptc_predicts_nmd_unfiltered_by_celltype.tsv"))
write_tsv(source_results, file.path(output_dir, "ptc_predicts_nmd_by_source.tsv"))
write_tsv(strat_pooled, file.path(output_dir, "ptc_distance_nmd_rate.tsv"))

strat_all_out <- strat_all %>%
  select(cell_type, ptc_distance_bin, n, n_nmd, nmd_rate, median_distance)
write_tsv(strat_all_out, file.path(output_dir, "ptc_distance_nmd_rate_by_celltype.tsv"))

write_tsv(biotype_summary, file.path(output_dir, "nmd_biotype_predicts_nmd_by_celltype.tsv"))
write_tsv(bt_pooled, file.path(output_dir, "nmd_rate_by_biotype.tsv"))

if (length(enst_results) > 0) {
  write_tsv(bind_rows(enst_results),
            file.path(output_dir, "nmd_biotype_predicts_nmd_enst_only.tsv"))
}

cat(sprintf("\nSaved %d output files to: %s\n", 8, output_dir))
cat("\nDone.\n")
