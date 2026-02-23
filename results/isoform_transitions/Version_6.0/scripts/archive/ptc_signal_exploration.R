#!/usr/bin/env Rscript
# PTC Signal Exploration
#
# Systematic search for subsets where PTC predicts NMD responsiveness.
# Explores multiple angles: extreme tails, GENCODE NMD biotype, large PTC
# distances, zero-DMSO isoforms, within-gene paired analysis, number of
# downstream EJCs, and novel vs annotated isoforms.

suppressPackageStartupMessages(library(tidyverse))

cat("\n")
cat("======================================================================\n")
cat("   PTC Signal Exploration: Searching for a Detectable PTC Effect\n")
cat("======================================================================\n\n")

# ─── 1. Load data ────────────────────────────────────────────────────────────

cat("Loading data...\n")
cds_meta <- readRDS("data/isoform_cds_metadata.rds")
structs <- readRDS("data/isoform_structures.rds")
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

# ─── 2. Compute PTC status with additional features ─────────────────────────

cat("Computing PTC status...\n")

coding <- cds_meta %>%
  filter(coding_status == "coding") %>%
  inner_join(structs, by = "isoform_id")

compute_ptc_features <- function(cds_start, cds_stop, strand,
                                  exon_starts, exon_ends, n_exons) {
  result <- list(ptc_distance = NA_real_, n_downstream_ejcs = NA_integer_,
                 stop_in_last_exon = NA)
  if (n_exons == 1) return(result)

  stop_pos <- if (strand == "+") cds_stop else cds_start
  if (strand == "+") { ord <- order(exon_starts) } else { ord <- order(exon_starts, decreasing = TRUE) }
  e_starts <- exon_starts[ord]; e_ends <- exon_ends[ord]
  exon_lengths <- e_ends - e_starts
  cum_starts <- c(0, cumsum(exon_lengths[-length(exon_lengths)]))
  last_ejc_mRNA <- sum(exon_lengths[-length(exon_lengths)])

  exon_idx <- which(e_starts <= stop_pos & stop_pos <= e_ends)
  if (length(exon_idx) == 0) return(result)
  exon_idx <- exon_idx[1]

  if (strand == "+") { offset_in_exon <- stop_pos - e_starts[exon_idx] }
  else { offset_in_exon <- e_ends[exon_idx] - stop_pos }
  stop_mRNA <- cum_starts[exon_idx] + offset_in_exon

  result$ptc_distance <- last_ejc_mRNA - stop_mRNA
  result$n_downstream_ejcs <- as.integer(n_exons - exon_idx)
  result$stop_in_last_exon <- (exon_idx == n_exons)
  result
}

ptc_results <- coding %>%
  filter(n_exons > 1) %>%
  rowwise() %>%
  mutate(
    features = list(compute_ptc_features(
      cds_start, cds_stop, as.character(strand),
      exon_starts, exon_ends, n_exons
    )),
    ptc_distance = features$ptc_distance,
    n_downstream_ejcs = features$n_downstream_ejcs,
    stop_in_last_exon = features$stop_in_last_exon
  ) %>%
  ungroup() %>%
  select(-features) %>%
  filter(!is.na(ptc_distance)) %>%
  mutate(
    has_ptc = ptc_distance > 50,
    is_novel = startsWith(isoform_id, "PB."),
    is_gencode = startsWith(isoform_id, "ENST")
  )

ptc_lookup <- ptc_results %>%
  select(isoform_id, gene_id, n_exons, ptc_distance, has_ptc,
         n_downstream_ejcs, stop_in_last_exon, is_novel, is_gencode, orf_length)

cat(sprintf("  %d multi-exon coding isoforms\n", nrow(ptc_lookup)))
cat(sprintf("  PTC+: %d, PTC-: %d\n", sum(ptc_lookup$has_ptc), sum(!ptc_lookup$has_ptc)))
cat(sprintf("  GENCODE: %d, PacBio novel: %d\n",
            sum(ptc_lookup$is_gencode), sum(ptc_lookup$is_novel)))

# ─── 3. Compute DMSO CPM ────────────────────────────────────────────────────

dmso_cpm <- expression_data %>%
  filter(treatment == "DMSO") %>%
  group_by(isoform_id, ct) %>%
  summarise(mean_dmso_cpm = mean(cpm), .groups = "drop")

# 5% filter function
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

# ─── 4. Load DE + biotype for all cell types ────────────────────────────────

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

# ─── 5. Build merged dataset ────────────────────────────────────────────────

cat("Building merged datasets...\n\n")

# Unfiltered: all DE results × PTC status
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

cat(sprintf("  Unfiltered: %d isoform-celltype observations\n", nrow(merged_unf)))
cat(sprintf("  5%%-filtered: %d isoform-celltype observations\n", nrow(merged_filt)))

# ═══════════════════════════════════════════════════════════════════════════
# EXPLORATION ANALYSES
# ═══════════════════════════════════════════════════════════════════════════

run_fisher <- function(df, label) {
  tab <- table(PTC = df$has_ptc, NMD = df$nmd_responsive)
  if (!all(dim(tab) == 2) || any(tab == 0)) return(NULL)
  ft <- fisher.test(tab)
  tibble(
    subset = label,
    n = nrow(df),
    n_ptc = sum(df$has_ptc),
    n_nmd = sum(df$nmd_responsive),
    nmd_rate_ptc = tab["TRUE","TRUE"] / sum(tab["TRUE",]),
    nmd_rate_no_ptc = tab["FALSE","TRUE"] / sum(tab["FALSE",]),
    OR = ft$estimate,
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    p = ft$p.value
  )
}

run_logfc_test <- function(df, label) {
  if (sum(df$has_ptc) < 20 || sum(!df$has_ptc) < 20) return(NULL)
  wt <- wilcox.test(logFC ~ has_ptc, data = df)
  tibble(
    subset = label,
    n = nrow(df),
    n_ptc = sum(df$has_ptc),
    median_logFC_ptc = median(df$logFC[df$has_ptc]),
    median_logFC_no = median(df$logFC[!df$has_ptc]),
    median_shift = median(df$logFC[df$has_ptc]) - median(df$logFC[!df$has_ptc]),
    mean_shift = mean(df$logFC[df$has_ptc]) - mean(df$logFC[!df$has_ptc]),
    p = wt$p.value
  )
}

# ═══════════════════════════════════════════════════════════════════════════
cat("\n")
cat("======================================================================\n")
cat("  ANGLE 1: Extreme Tail Enrichment\n")
cat("  Are PTC isoforms overrepresented among the most upregulated?\n")
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
      bottom <- ct_df %>% slice_min(logFC, n = top_n)
      rest <- ct_df %>% filter(!isoform_id %in% c(top$isoform_id, bottom$isoform_id))

      ptc_top <- mean(top$has_ptc)
      ptc_rest <- mean(rest$has_ptc)
      ptc_bottom <- mean(bottom$has_ptc)

      # Fisher: PTC enrichment in top vs rest
      tab <- matrix(c(sum(top$has_ptc), sum(!top$has_ptc),
                       sum(rest$has_ptc), sum(!rest$has_ptc)), nrow = 2)
      ft <- fisher.test(tab)

      if (ft$p.value < 0.05) {
        cat(sprintf("  %s top%d: PTC=%.1f%% vs rest=%.1f%% vs bottom=%.1f%%, OR=%.2f, p=%g\n",
                    ct, top_n, 100*ptc_top, 100*ptc_rest, 100*ptc_bottom,
                    ft$estimate, ft$p.value))
      }
    }
  }
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 2: GENCODE NMD Biotype Isoforms Only\n")
cat("  Restrict to gold-standard annotated NMD transcripts\n")
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

# ═══════════════════════════════════════════════════════════════════════════
cat("\n")
cat("======================================================================\n")
cat("  ANGLE 3: Large PTC Distances (>200nt, >500nt, >1000nt)\n")
cat("  Stronger NMD trigger → bigger effect?\n")
cat("======================================================================\n\n")

for (dataset_label in c("unfiltered", "5pct_filtered")) {
  df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
  cat(sprintf("--- %s ---\n", dataset_label))

  results_a3 <- list()
  for (ct in cell_types) {
    for (dist_thresh in c(200, 500, 1000)) {
      ct_df <- df %>%
        filter(cell_type == ct) %>%
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

# ═══════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 4: Multiple Downstream EJCs (>= 2, >= 3)\n")
cat("  More EJCs downstream of stop = stronger NMD\n")
cat("======================================================================\n\n")

for (dataset_label in c("unfiltered", "5pct_filtered")) {
  df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
  cat(sprintf("--- %s ---\n", dataset_label))

  results_a4 <- list()
  for (ct in cell_types) {
    for (min_ejcs in c(2, 3, 5)) {
      ct_df <- df %>%
        filter(cell_type == ct) %>%
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

# ═══════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 5: Zero DMSO Expression (Completely Silenced)\n")
cat("  True NMD targets might have zero baseline expression\n")
cat("======================================================================\n\n")

results_a5 <- list()
for (ct in cell_types) {
  ct_df <- merged_unf %>%
    filter(cell_type == ct, mean_dmso_cpm == 0)
  r <- run_logfc_test(ct_df, sprintf("zero_dmso_%s", ct))
  if (!is.null(r)) results_a5[[length(results_a5) + 1]] <- r
}
if (length(results_a5) > 0) {
  cat("Isoforms with zero DMSO expression:\n")
  bind_rows(results_a5) %>%
    mutate(across(where(is.numeric), ~round(., 4))) %>%
    print(n = 10, width = 140)
}

# ═══════════════════════════════════════════════════════════════════════════
cat("\n")
cat("======================================================================\n")
cat("  ANGLE 6: Novel (PacBio) vs GENCODE Isoforms\n")
cat("  PTC predictions may be more reliable for GENCODE\n")
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

# ═══════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 7: Within-Gene Paired Analysis\n")
cat("  For genes with both PTC+ and PTC- isoforms, compare logFC\n")
cat("======================================================================\n\n")

for (dataset_label in c("unfiltered", "5pct_filtered")) {
  df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
  cat(sprintf("--- %s ---\n", dataset_label))

  for (ct in cell_types) {
    ct_df <- df %>% filter(cell_type == ct)

    # Find genes with at least one PTC+ and one PTC- isoform
    gene_ptc_status <- ct_df %>%
      group_by(gene_id) %>%
      summarise(has_both = any(has_ptc) & any(!has_ptc), .groups = "drop") %>%
      filter(has_both)

    paired_df <- ct_df %>% filter(gene_id %in% gene_ptc_status$gene_id)

    # For each gene, compute mean logFC for PTC+ vs PTC-
    gene_means <- paired_df %>%
      group_by(gene_id, has_ptc) %>%
      summarise(mean_logFC = mean(logFC), .groups = "drop") %>%
      pivot_wider(names_from = has_ptc, values_from = mean_logFC,
                  names_prefix = "ptc_") %>%
      mutate(diff = ptc_TRUE - ptc_FALSE)

    if (nrow(gene_means) >= 20) {
      wt <- wilcox.test(gene_means$diff, mu = 0)
      cat(sprintf("  %s: %d genes with both PTC+/PTC- isoforms\n", ct, nrow(gene_means)))
      cat(sprintf("    Mean within-gene logFC diff (PTC+ minus PTC-): %.4f\n", mean(gene_means$diff, na.rm = TRUE)))
      cat(sprintf("    Median: %.4f, Wilcoxon p = %g\n", median(gene_means$diff, na.rm = TRUE), wt$p.value))
      cat(sprintf("    %% genes where PTC+ has higher logFC: %.1f%%\n",
                  100 * mean(gene_means$diff > 0, na.rm = TRUE)))
    }
  }
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 8: Combined Strong PTC Features\n")
cat("  PTC distance > 500 AND >= 2 downstream EJCs AND GENCODE isoform\n")
cat("======================================================================\n\n")

for (dataset_label in c("unfiltered", "5pct_filtered")) {
  df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
  cat(sprintf("--- %s ---\n", dataset_label))

  results_a8 <- list()
  for (ct in cell_types) {
    ct_df <- df %>%
      filter(cell_type == ct) %>%
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

# ═══════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 9: Quantile Enrichment (QQ-style)\n")
cat("  PTC fraction in sliding quantiles of logFC\n")
cat("======================================================================\n\n")

cat("PTC fraction by logFC decile (pooled across cell types, unfiltered):\n\n")

# Pool all cell types
for (dataset_label in c("unfiltered", "5pct_filtered")) {
  df <- if (dataset_label == "unfiltered") merged_unf else merged_filt
  cat(sprintf("--- %s ---\n", dataset_label))

  for (ct in cell_types) {
    ct_df <- df %>% filter(cell_type == ct)
    ct_df <- ct_df %>%
      mutate(logFC_decile = ntile(logFC, 10))

    decile_ptc <- ct_df %>%
      group_by(logFC_decile) %>%
      summarise(
        n = n(),
        n_ptc = sum(has_ptc),
        ptc_frac = mean(has_ptc),
        mean_logFC = mean(logFC),
        .groups = "drop"
      )

    # Is there a trend? Spearman correlation between decile and PTC fraction
    st <- cor.test(decile_ptc$logFC_decile, decile_ptc$ptc_frac, method = "spearman")

    if (st$p.value < 0.1) {
      cat(sprintf("  %s: Spearman rho = %.3f, p = %g\n", ct, st$estimate, st$p.value))
      cat("    Decile PTC fractions: ")
      cat(paste(sprintf("%.1f%%", 100 * decile_ptc$ptc_frac), collapse = " → "))
      cat("\n")
    }
  }
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════
cat("======================================================================\n")
cat("  ANGLE 10: Very Low DMSO Expression (bottom 10% of expressed)\n")
cat("  NMD-degraded isoforms should be low but detectable in DMSO\n")
cat("======================================================================\n\n")

for (ct in cell_types) {
  ct_df <- merged_unf %>%
    filter(cell_type == ct, mean_dmso_cpm > 0)

  q10 <- quantile(ct_df$mean_dmso_cpm, 0.10)
  low_df <- ct_df %>% filter(mean_dmso_cpm <= q10)
  r <- run_logfc_test(low_df, sprintf("bottom10pct_dmso_%s", ct))
  if (!is.null(r)) {
    cat(sprintf("  %s (DMSO CPM <= %.2f): n=%d, median shift=%.4f, p=%g\n",
                ct, q10, r$n, r$median_shift, r$p))
  }
}

# ═══════════════════════════════════════════════════════════════════════════
cat("\n")
cat("======================================================================\n")
cat("  SUMMARY: Best Signals Found\n")
cat("======================================================================\n\n")

cat("Collecting all logFC shift tests across all angles...\n\n")

# Re-run everything and collect into one summary
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

cat("Top 30 most significant PTC logFC shifts (any direction):\n\n")
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
output_dir <- "comparisons/nonNMD_0.50/cross_comparison"
write_tsv(summary_all, file.path(output_dir, "ptc_signal_exploration_all_tests.tsv"))
cat(sprintf("\nSaved: %s\n", file.path(output_dir, "ptc_signal_exploration_all_tests.tsv")))

cat("\nDone.\n")
