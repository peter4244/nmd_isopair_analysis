#!/usr/bin/env Rscript
# GENCODE NMD Biotype Predicts NMD Responsiveness
#
# Question: Among expressed isoforms passing the 5% filter, does GENCODE
# "nonsense_mediated_decay" biotype predict NMD responsiveness
# (adj.P.Val < 0.05, logFC > 0)?
#
# Mirrors ptc_predicts_nmd.R but uses annotation-based NMD status instead
# of PTC distance.

suppressPackageStartupMessages(library(tidyverse))

cat("\n")
cat("==================================================================\n")
cat("   GENCODE NMD Biotype → NMD Responsiveness Analysis\n")
cat("==================================================================\n\n")

# ─── 1. Load data ────────────────────────────────────────────────────────────

cat("Loading data...\n")
expression_data <- readRDS("data/expression_data_filtered.rds")
sample_metadata <- readRDS("data/sample_metadata.rds")

expression_data <- expression_data %>%
  left_join(sample_metadata %>% select(sample_id, treatment, ct), by = "sample_id")

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

# ─── 3. Per-cell-type analysis ───────────────────────────────────────────────

all_results <- list()
all_biotype_tables <- list()

for (ct in cell_types) {
  cat(sprintf("\n  --- %s ---\n", ct))

  # Load DE results with biotype
  de_file <- file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix[ct], "_", de_date, ".csv"))
  de <- read_csv(de_file, show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC, biotype) %>%
    rename(isoform_id = txid)

  # Get samples and apply 5% filter
  ct_samples <- sample_metadata %>%
    filter(ct == !!ct) %>%
    pull(sample_id)

  passing_5pct <- apply_5pct_filter(expression_data, ct_samples)
  cat(sprintf("    Isoforms passing 5%% filter: %d\n", length(passing_5pct)))

  # Filter to 5%-passing isoforms with DE results
  analysis_df <- de %>%
    filter(isoform_id %in% passing_5pct, !is.na(biotype))

  cat(sprintf("    Isoforms with biotype annotation: %d\n", nrow(analysis_df)))

  # Classify
  analysis_df <- analysis_df %>%
    mutate(
      nmd_biotype = biotype == "nonsense_mediated_decay",
      nmd_responsive = adj.P.Val < 0.05 & logFC > 0
    )

  cat(sprintf("    GENCODE NMD biotype: %d\n", sum(analysis_df$nmd_biotype)))
  cat(sprintf("    NMD-responsive (FDR<0.05, logFC>0): %d\n", sum(analysis_df$nmd_responsive)))

  # ── Cross-tabulation ──
  tab <- table(NMD_biotype = analysis_df$nmd_biotype,
               NMD_responsive = analysis_df$nmd_responsive)
  cat("    NMD biotype × NMD responsiveness:\n")
  print(tab)

  if (all(dim(tab) == 2)) {
    ft <- fisher.test(tab)
    nmd_rate_biotype <- tab["TRUE", "TRUE"] / sum(tab["TRUE", ])
    nmd_rate_other <- tab["FALSE", "TRUE"] / sum(tab["FALSE", ])

    cat(sprintf("    Fisher OR = %.2f (%.2f–%.2f), p = %g\n",
                ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
    cat(sprintf("    NMD rate with NMD biotype: %.1f%% (%d/%d)\n",
                100 * nmd_rate_biotype, tab["TRUE", "TRUE"], sum(tab["TRUE", ])))
    cat(sprintf("    NMD rate without NMD biotype: %.1f%% (%d/%d)\n",
                100 * nmd_rate_other, tab["FALSE", "TRUE"], sum(tab["FALSE", ])))

    all_results[[ct]] <- tibble(
      cell_type = ct,
      n_total = nrow(analysis_df),
      n_nmd_biotype = sum(analysis_df$nmd_biotype),
      n_other_biotype = sum(!analysis_df$nmd_biotype),
      n_nmd_responsive = sum(analysis_df$nmd_responsive),
      nmd_rate_biotype = nmd_rate_biotype,
      nmd_rate_other = nmd_rate_other,
      rate_ratio = nmd_rate_biotype / nmd_rate_other,
      odds_ratio = ft$estimate,
      or_ci_lower = ft$conf.int[1],
      or_ci_upper = ft$conf.int[2],
      p_value = ft$p.value
    )
  }

  # ── Biotype breakdown for NMD-responsive isoforms ──
  bt_table <- analysis_df %>%
    group_by(biotype) %>%
    summarise(
      n = n(),
      n_nmd = sum(nmd_responsive),
      nmd_rate = mean(nmd_responsive),
      .groups = "drop"
    ) %>%
    arrange(desc(n_nmd)) %>%
    mutate(cell_type = ct)

  all_biotype_tables[[ct]] <- bt_table
}

# ─── 4. Summary across cell types ────────────────────────────────────────────

cat("\n\n")
cat("==================================================================\n")
cat("   SUMMARY: GENCODE NMD Biotype → NMD Responsiveness\n")
cat("==================================================================\n\n")

summary_df <- bind_rows(all_results) %>%
  mutate(fdr = p.adjust(p_value, method = "BH"))

cat("Per-cell-type results:\n\n")
summary_df %>%
  mutate(
    across(c(nmd_rate_biotype, nmd_rate_other), ~ sprintf("%.1f%%", 100 * .)),
    rate_ratio = round(rate_ratio, 2),
    across(c(odds_ratio, or_ci_lower, or_ci_upper), ~ round(., 2)),
    p_value = signif(p_value, 3),
    fdr = signif(fdr, 3)
  ) %>%
  select(cell_type, n_total, n_nmd_biotype, nmd_rate_biotype, nmd_rate_other,
         rate_ratio, odds_ratio, or_ci_lower, or_ci_upper, p_value, fdr) %>%
  print(n = 10, width = 140)

# ── Meta-analysis ──
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

Q <- sum(weights * (meta_df$log_or - pooled_log_or)^2)
k <- nrow(meta_df)
I2 <- max(0, (Q - (k - 1)) / Q) * 100
Q_p <- pchisq(Q, df = k - 1, lower.tail = FALSE)

cat(sprintf("  Pooled OR = %.2f (%.2f–%.2f), p = %g\n",
            pooled_or, pooled_ci_lower, pooled_ci_upper, pooled_p))
cat(sprintf("  Heterogeneity: Q = %.1f (p = %g), I² = %.0f%%\n", Q, Q_p, I2))

# ─── 5. Biotype breakdown (top biotypes across all cell types) ───────────────

cat("\n\n")
cat("==================================================================\n")
cat("   NMD RESPONSIVENESS BY BIOTYPE (pooled across cell types)\n")
cat("==================================================================\n\n")

bt_all <- bind_rows(all_biotype_tables)

bt_pooled <- bt_all %>%
  group_by(biotype) %>%
  summarise(
    n_cell_types = n(),
    total_n = sum(n),
    total_nmd = sum(n_nmd),
    pooled_nmd_rate = sum(n_nmd) / sum(n),
    mean_nmd_rate = mean(nmd_rate),
    .groups = "drop"
  ) %>%
  filter(total_n >= 100) %>%
  arrange(desc(pooled_nmd_rate))

bt_pooled %>%
  mutate(pooled_nmd_rate = sprintf("%.2f%%", 100 * pooled_nmd_rate),
         mean_nmd_rate = sprintf("%.2f%%", 100 * mean_nmd_rate)) %>%
  print(n = 30, width = 120)

# ─── 6. Restrict to GENCODE-annotated (ENST) isoforms only ──────────────────

cat("\n\n")
cat("==================================================================\n")
cat("   SENSITIVITY: GENCODE (ENST) isoforms only\n")
cat("==================================================================\n\n")

cat("PacBio novel isoforms lack GENCODE biotype. Restricting to ENST isoforms:\n\n")

enst_results <- list()

for (ct in cell_types) {
  de_file <- file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix[ct], "_", de_date, ".csv"))
  de <- read_csv(de_file, show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC, biotype) %>%
    rename(isoform_id = txid)

  ct_samples <- sample_metadata %>%
    filter(ct == !!ct) %>%
    pull(sample_id)
  passing_5pct <- apply_5pct_filter(expression_data, ct_samples)

  # Restrict to ENST isoforms only
  analysis_df <- de %>%
    filter(isoform_id %in% passing_5pct,
           startsWith(isoform_id, "ENST"),
           !is.na(biotype)) %>%
    mutate(
      nmd_biotype = biotype == "nonsense_mediated_decay",
      nmd_responsive = adj.P.Val < 0.05 & logFC > 0
    )

  tab <- table(NMD_biotype = analysis_df$nmd_biotype,
               NMD_responsive = analysis_df$nmd_responsive)

  if (all(dim(tab) == 2)) {
    ft <- fisher.test(tab)
    nmd_rate_biotype <- tab["TRUE", "TRUE"] / sum(tab["TRUE", ])
    nmd_rate_other <- tab["FALSE", "TRUE"] / sum(tab["FALSE", ])

    cat(sprintf("  %s: NMD biotype %.1f%% (%d/%d) vs other %.1f%% (%d/%d), OR=%.2f (%.2f–%.2f), p=%g\n",
                ct,
                100 * nmd_rate_biotype, tab["TRUE", "TRUE"], sum(tab["TRUE", ]),
                100 * nmd_rate_other, tab["FALSE", "TRUE"], sum(tab["FALSE", ]),
                ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))

    enst_results[[ct]] <- tibble(
      cell_type = ct,
      n_total = nrow(analysis_df),
      n_nmd_biotype = sum(analysis_df$nmd_biotype),
      nmd_rate_biotype = nmd_rate_biotype,
      nmd_rate_other = nmd_rate_other,
      odds_ratio = ft$estimate,
      or_ci_lower = ft$conf.int[1],
      or_ci_upper = ft$conf.int[2],
      p_value = ft$p.value
    )
  }
}

# Meta-analysis for ENST-only
cat("\n  ENST-only meta-analysis:\n")
enst_df <- bind_rows(enst_results) %>%
  mutate(
    log_or = log(odds_ratio),
    se_log_or = (log(or_ci_upper) - log(or_ci_lower)) / (2 * 1.96)
  )

w <- 1 / enst_df$se_log_or^2
p_log_or <- sum(w * enst_df$log_or) / sum(w)
p_se <- sqrt(1 / sum(w))
cat(sprintf("  Pooled OR = %.2f (%.2f–%.2f), p = %g\n",
            exp(p_log_or), exp(p_log_or - 1.96 * p_se), exp(p_log_or + 1.96 * p_se),
            2 * pnorm(-abs(p_log_or / p_se))))

Q_e <- sum(w * (enst_df$log_or - p_log_or)^2)
I2_e <- max(0, (Q_e - (nrow(enst_df) - 1)) / Q_e) * 100
cat(sprintf("  Heterogeneity: Q = %.1f, I² = %.0f%%\n", Q_e, I2_e))

# ─── 7. Save outputs ────────────────────────────────────────────────────────

output_dir <- "comparisons/nonNMD_0.50/cross_comparison"

write_tsv(summary_df, file.path(output_dir, "nmd_biotype_predicts_nmd_by_celltype.tsv"))
cat(sprintf("\n\nSaved: %s\n", file.path(output_dir, "nmd_biotype_predicts_nmd_by_celltype.tsv")))

write_tsv(bt_pooled, file.path(output_dir, "nmd_rate_by_biotype.tsv"))
cat(sprintf("Saved: %s\n", file.path(output_dir, "nmd_rate_by_biotype.tsv")))

if (length(enst_results) > 0) {
  write_tsv(bind_rows(enst_results), file.path(output_dir, "nmd_biotype_predicts_nmd_enst_only.tsv"))
  cat(sprintf("Saved: %s\n", file.path(output_dir, "nmd_biotype_predicts_nmd_enst_only.tsv")))
}

cat("\nDone.\n")
