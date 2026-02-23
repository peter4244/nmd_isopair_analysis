#!/usr/bin/env Rscript
# 04_ptc_signal_exploration.R
#
# Systematic search for subsets where PTC predicts NMD responsiveness.
# Explores 10 angles: extreme tails, GENCODE NMD biotype, large PTC distances,
# zero-DMSO isoforms, within-gene paired analysis, multiple downstream EJCs,
# novel vs annotated isoforms, combined strong PTC, quantile enrichment,
# and very low DMSO expression.
#
# Requires: ptc_status.rds from 01_compute_ptc_status.R
#
# Outputs:
#   ptc_signal_exploration_all_tests.tsv
#
# Run from: results/isoform_transitions/Version_6.0/

suppressPackageStartupMessages(library(tidyverse))

cat("\n")
cat("======================================================================\n")
cat("   04: PTC Signal Exploration\n")
cat("======================================================================\n\n")

# ─── 1. Load data ────────────────────────────────────────────────────────────

cat("Loading data...\n")
ptc_status      <- readRDS("comparisons/nonNMD_0.50/ptc/results/ptc_status.rds")
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

ptc_lookup <- ptc_status %>% filter(!is.na(has_ptc), n_exons > 1)

# ─── 2. DMSO CPM + 5% filter ────────────────────────────────────────────────

dmso_cpm <- expression_data %>%
  filter(treatment == "DMSO") %>%
  group_by(isoform_id, ct) %>%
  summarise(mean_dmso_cpm = mean(cpm), .groups = "drop")

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

# ─── 3. Load DE + build merged datasets ─────────────────────────────────────

cat("Loading DE results...\n")

all_de <- list()
for (ct in cell_types) {
  de_file <- file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix[ct], "_", de_date, ".csv"))
  de <- read_csv(de_file, show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC, biotype) %>%
    rename(isoform_id = txid) %>%
    mutate(cell_type = ct)
  all_de[[ct]] <- de
}
all_de <- bind_rows(all_de)

# Unfiltered
merged_unf <- all_de %>%
  inner_join(ptc_lookup, by = "isoform_id") %>%
  left_join(dmso_cpm, by = c("isoform_id", "cell_type" = "ct")) %>%
  mutate(
    mean_dmso_cpm = replace_na(mean_dmso_cpm, 0),
    nmd_responsive = adj.P.Val < 0.05 & logFC > 0,
    nmd_biotype = !is.na(biotype) & biotype == "nonsense_mediated_decay"
  )

# 5%-filtered
passing_5pct <- list()
for (ct in cell_types) {
  ct_samples <- sample_metadata %>% filter(ct == !!ct) %>% pull(sample_id)
  passing_5pct[[ct]] <- tibble(isoform_id = apply_5pct_filter(expression_data, ct_samples),
                                cell_type = ct)
}
passing_5pct <- bind_rows(passing_5pct)

merged_filt <- merged_unf %>%
  semi_join(passing_5pct, by = c("isoform_id", "cell_type"))

cat(sprintf("  Unfiltered: %d observations\n", nrow(merged_unf)))
cat(sprintf("  5%%-filtered: %d observations\n", nrow(merged_filt)))

# ─── Helpers ─────────────────────────────────────────────────────────────────

run_logfc_test <- function(df, label) {
  if (sum(df$has_ptc) < 20 || sum(!df$has_ptc) < 20) return(NULL)
  wt <- wilcox.test(logFC ~ has_ptc, data = df)
  tibble(
    subset = label, n = nrow(df), n_ptc = sum(df$has_ptc),
    median_logFC_ptc = median(df$logFC[df$has_ptc]),
    median_logFC_no = median(df$logFC[!df$has_ptc]),
    median_shift = median(df$logFC[df$has_ptc]) - median(df$logFC[!df$has_ptc]),
    mean_shift = mean(df$logFC[df$has_ptc]) - mean(df$logFC[!df$has_ptc]),
    p = wt$p.value
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat("======================================================================\n")
cat("  ANGLE 1: Extreme Tail Enrichment\n")
cat("======================================================================\n\n")

for (dataset_label in c("unfiltered", "5pct_filtered")) {
  df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
  cat(sprintf("--- %s ---\n", dataset_label))
  for (ct in cell_types) {
    ct_df <- df %>% filter(cell_type == ct)
    n_total <- nrow(ct_df)
    for (top_n in c(100, 500, 1000)) {
      if (n_total < top_n * 2) next
      top <- ct_df %>% slice_max(logFC, n = top_n)
      rest <- ct_df %>% filter(!isoform_id %in% top$isoform_id)
      ptc_top <- mean(top$has_ptc)
      ptc_rest <- mean(rest$has_ptc)
      tab <- matrix(c(sum(top$has_ptc), sum(!top$has_ptc),
                       sum(rest$has_ptc), sum(!rest$has_ptc)), nrow = 2)
      ft <- fisher.test(tab)
      if (ft$p.value < 0.05) {
        cat(sprintf("  %s top%d: PTC=%.1f%% vs rest=%.1f%%, OR=%.2f, p=%g\n",
                    ct, top_n, 100*ptc_top, 100*ptc_rest, ft$estimate, ft$p.value))
      }
    }
  }
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 2: GENCODE NMD Biotype Isoforms Only\n")
cat("======================================================================\n\n")

results_a2 <- list()
for (ct in cell_types) {
  for (dataset_label in c("unfiltered", "5pct_filtered")) {
    df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
    ct_df <- df %>% filter(cell_type == ct, nmd_biotype == TRUE)
    r <- run_logfc_test(ct_df, sprintf("NMD_biotype_%s_%s", ct, dataset_label))
    if (!is.null(r)) results_a2[[length(results_a2) + 1]] <- r
  }
}
if (length(results_a2) > 0) {
  bind_rows(results_a2) %>%
    mutate(across(where(is.numeric), ~round(., 4))) %>%
    print(n = 20, width = 140)
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat("======================================================================\n")
cat("  ANGLE 3: Large PTC Distances (>200nt, >500nt, >1000nt)\n")
cat("======================================================================\n\n")

for (dataset_label in c("unfiltered", "5pct_filtered")) {
  df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
  cat(sprintf("--- %s ---\n", dataset_label))
  results_a3 <- list()
  for (ct in cell_types) {
    for (dist_thresh in c(200, 500, 1000)) {
      ct_df <- df %>% filter(cell_type == ct) %>%
        mutate(has_ptc = ptc_distance > dist_thresh)
      r <- run_logfc_test(ct_df, sprintf("%s_PTC>%d", ct, dist_thresh))
      if (!is.null(r)) results_a3[[length(results_a3) + 1]] <- r
    }
  }
  if (length(results_a3) > 0) {
    bind_rows(results_a3) %>%
      mutate(across(where(is.numeric), ~round(., 4))) %>%
      print(n = 30, width = 140)
  }
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 4: Multiple Downstream EJCs (>= 2, >= 3, >= 5)\n")
cat("======================================================================\n\n")

for (dataset_label in c("unfiltered", "5pct_filtered")) {
  df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
  cat(sprintf("--- %s ---\n", dataset_label))
  results_a4 <- list()
  for (ct in cell_types) {
    for (min_ejcs in c(2, 3, 5)) {
      ct_df <- df %>% filter(cell_type == ct) %>%
        mutate(has_ptc = n_downstream_ejcs >= min_ejcs)
      r <- run_logfc_test(ct_df, sprintf("%s_EJC>=%d", ct, min_ejcs))
      if (!is.null(r)) results_a4[[length(results_a4) + 1]] <- r
    }
  }
  if (length(results_a4) > 0) {
    bind_rows(results_a4) %>%
      mutate(across(where(is.numeric), ~round(., 4))) %>%
      print(n = 30, width = 140)
  }
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 5: Zero DMSO Expression\n")
cat("======================================================================\n\n")

results_a5 <- list()
for (ct in cell_types) {
  ct_df <- merged_unf %>% filter(cell_type == ct, mean_dmso_cpm == 0)
  r <- run_logfc_test(ct_df, sprintf("zero_dmso_%s", ct))
  if (!is.null(r)) results_a5[[length(results_a5) + 1]] <- r
}
if (length(results_a5) > 0) {
  bind_rows(results_a5) %>%
    mutate(across(where(is.numeric), ~round(., 4))) %>%
    print(n = 10, width = 140)
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat("======================================================================\n")
cat("  ANGLE 6: Novel (PacBio) vs GENCODE Isoforms\n")
cat("======================================================================\n\n")

for (dataset_label in c("unfiltered", "5pct_filtered")) {
  df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
  cat(sprintf("--- %s ---\n", dataset_label))
  results_a6 <- list()
  for (ct in cell_types) {
    for (iso_type in c("GENCODE", "PacBio")) {
      ct_df <- df %>%
        filter(cell_type == ct,
               if (iso_type == "GENCODE") is_gencode else is_novel)
      r <- run_logfc_test(ct_df, sprintf("%s_%s", ct, iso_type))
      if (!is.null(r)) results_a6[[length(results_a6) + 1]] <- r
    }
  }
  if (length(results_a6) > 0) {
    bind_rows(results_a6) %>%
      mutate(across(where(is.numeric), ~round(., 4))) %>%
      print(n = 20, width = 140)
  }
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 7: Within-Gene Paired Analysis\n")
cat("======================================================================\n\n")

for (dataset_label in c("unfiltered", "5pct_filtered")) {
  df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
  cat(sprintf("--- %s ---\n", dataset_label))
  for (ct in cell_types) {
    ct_df <- df %>% filter(cell_type == ct)
    gene_ptc_status <- ct_df %>%
      group_by(gene_id) %>%
      summarise(has_both = any(has_ptc) & any(!has_ptc), .groups = "drop") %>%
      filter(has_both)
    paired_df <- ct_df %>% filter(gene_id %in% gene_ptc_status$gene_id)
    gene_means <- paired_df %>%
      group_by(gene_id, has_ptc) %>%
      summarise(mean_logFC = mean(logFC), .groups = "drop") %>%
      pivot_wider(names_from = has_ptc, values_from = mean_logFC,
                  names_prefix = "ptc_") %>%
      mutate(diff = ptc_TRUE - ptc_FALSE)
    if (nrow(gene_means) >= 20) {
      wt <- wilcox.test(gene_means$diff, mu = 0)
      cat(sprintf("  %s: %d genes, mean diff=%.4f, median=%.4f, p=%g, %%PTC+ higher=%.1f%%\n",
                  ct, nrow(gene_means), mean(gene_means$diff, na.rm = TRUE),
                  median(gene_means$diff, na.rm = TRUE), wt$p.value,
                  100 * mean(gene_means$diff > 0, na.rm = TRUE)))
    }
  }
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 8: Combined Strong PTC Features\n")
cat("  PTC distance > 500 AND >= 2 downstream EJCs AND GENCODE\n")
cat("======================================================================\n\n")

for (dataset_label in c("unfiltered", "5pct_filtered")) {
  df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
  cat(sprintf("--- %s ---\n", dataset_label))
  results_a8 <- list()
  for (ct in cell_types) {
    ct_df <- df %>% filter(cell_type == ct) %>%
      mutate(strong_ptc = ptc_distance > 500 & n_downstream_ejcs >= 2 & is_gencode)
    n_strong <- sum(ct_df$strong_ptc)
    if (n_strong >= 20) {
      r <- run_logfc_test(ct_df %>% mutate(has_ptc = strong_ptc),
                          sprintf("%s_strong_PTC", ct))
      if (!is.null(r)) {
        r$n_strong_ptc <- n_strong
        results_a8[[length(results_a8) + 1]] <- r
      }
    }
  }
  if (length(results_a8) > 0) {
    bind_rows(results_a8) %>%
      mutate(across(where(is.numeric), ~round(., 4))) %>%
      print(n = 10, width = 140)
  }
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 9: Quantile Enrichment\n")
cat("======================================================================\n\n")

for (dataset_label in c("unfiltered", "5pct_filtered")) {
  df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
  cat(sprintf("--- %s ---\n", dataset_label))
  for (ct in cell_types) {
    ct_df <- df %>% filter(cell_type == ct) %>%
      mutate(logFC_decile = ntile(logFC, 10))
    decile_ptc <- ct_df %>%
      group_by(logFC_decile) %>%
      summarise(ptc_frac = mean(has_ptc), .groups = "drop")
    st <- cor.test(decile_ptc$logFC_decile, decile_ptc$ptc_frac, method = "spearman")
    if (st$p.value < 0.1) {
      cat(sprintf("  %s: rho = %.3f, p = %g\n", ct, st$estimate, st$p.value))
    }
  }
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 10: Very Low DMSO Expression (bottom 10%)\n")
cat("======================================================================\n\n")

for (ct in cell_types) {
  ct_df <- merged_unf %>% filter(cell_type == ct, mean_dmso_cpm > 0)
  q10 <- quantile(ct_df$mean_dmso_cpm, 0.10)
  low_df <- ct_df %>% filter(mean_dmso_cpm <= q10)
  r <- run_logfc_test(low_df, sprintf("bottom10pct_dmso_%s", ct))
  if (!is.null(r)) {
    cat(sprintf("  %s (DMSO CPM <= %.2f): n=%d, median shift=%.4f, p=%g\n",
                ct, q10, r$n, r$median_shift, r$p))
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat("======================================================================\n")
cat("  SUMMARY: All Tests Collected\n")
cat("======================================================================\n\n")

all_tests <- list()

for (ct in cell_types) {
  for (dataset_label in c("unfiltered", "5pct_filtered")) {
    df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
    ct_df <- df %>% filter(cell_type == ct)

    # Basic PTC
    r <- run_logfc_test(ct_df, sprintf("%s_%s_basic_PTC", ct, dataset_label))
    if (!is.null(r)) all_tests[[length(all_tests) + 1]] <- r

    # Large PTC distance
    for (d in c(500, 1000)) {
      r <- run_logfc_test(ct_df %>% mutate(has_ptc = ptc_distance > d),
                          sprintf("%s_%s_PTC>%d", ct, dataset_label, d))
      if (!is.null(r)) all_tests[[length(all_tests) + 1]] <- r
    }

    # Multiple EJCs
    r <- run_logfc_test(ct_df %>% mutate(has_ptc = n_downstream_ejcs >= 3),
                        sprintf("%s_%s_EJC>=3", ct, dataset_label))
    if (!is.null(r)) all_tests[[length(all_tests) + 1]] <- r

    # GENCODE only
    r <- run_logfc_test(ct_df %>% filter(is_gencode),
                        sprintf("%s_%s_GENCODE_only", ct, dataset_label))
    if (!is.null(r)) all_tests[[length(all_tests) + 1]] <- r

    # Novel only
    r <- run_logfc_test(ct_df %>% filter(is_novel),
                        sprintf("%s_%s_PacBio_only", ct, dataset_label))
    if (!is.null(r)) all_tests[[length(all_tests) + 1]] <- r

    # Strong PTC
    r <- run_logfc_test(
      ct_df %>% mutate(has_ptc = ptc_distance > 500 & n_downstream_ejcs >= 2 & is_gencode),
      sprintf("%s_%s_strong_PTC", ct, dataset_label))
    if (!is.null(r)) all_tests[[length(all_tests) + 1]] <- r
  }

  # Zero DMSO (unfiltered only)
  r <- run_logfc_test(merged_unf %>% filter(cell_type == ct, mean_dmso_cpm == 0),
                      sprintf("%s_zero_dmso", ct))
  if (!is.null(r)) all_tests[[length(all_tests) + 1]] <- r

  # Bottom 10% DMSO
  ct_expr <- merged_unf %>% filter(cell_type == ct, mean_dmso_cpm > 0)
  q10 <- quantile(ct_expr$mean_dmso_cpm, 0.10)
  r <- run_logfc_test(ct_expr %>% filter(mean_dmso_cpm <= q10),
                      sprintf("%s_bottom10pct_dmso", ct))
  if (!is.null(r)) all_tests[[length(all_tests) + 1]] <- r
}

summary_all <- bind_rows(all_tests) %>%
  mutate(
    direction = case_when(median_shift > 0 ~ "PTC higher (expected)",
                          median_shift < 0 ~ "PTC lower (reversed)",
                          TRUE ~ "equal"),
    fdr = p.adjust(p, method = "BH")
  ) %>%
  arrange(p)

cat("Top 30 most significant PTC logFC shifts:\n\n")
summary_all %>%
  head(30) %>%
  mutate(across(where(is.numeric), ~signif(., 3))) %>%
  print(n = 30, width = 150)

cat("\n\nTop 20 with EXPECTED direction (PTC higher logFC):\n\n")
summary_all %>%
  filter(direction == "PTC higher (expected)") %>%
  head(20) %>%
  mutate(across(where(is.numeric), ~signif(., 3))) %>%
  print(n = 20, width = 150)

# Save
output_dir <- "comparisons/nonNMD_0.50/ptc/results"
write_tsv(summary_all, file.path(output_dir, "ptc_signal_exploration_all_tests.tsv"))
cat(sprintf("\nSaved: %s\n", file.path(output_dir, "ptc_signal_exploration_all_tests.tsv")))

cat("\nDone.\n")
