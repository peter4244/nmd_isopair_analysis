#!/usr/bin/env Rscript
# PTC logFC Distributions
#
# Show distribution of beta-coefficients (logFC) for PTC+ vs PTC- isoforms,
# stratified by cell type. Then repeat stratified by DMSO expression level.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

cat("\n")
cat("==================================================================\n")
cat("   PTC logFC Distributions\n")
cat("==================================================================\n\n")

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
  "AT"     = "at2",
  "DD"     = "dd",
  "DD_ALI" = "ddali",
  "DO"     = "doali",
  "FB"     = "fb",
  "MV"     = "mv"
)

cell_types <- sort(unique(as.character(sample_metadata$ct)))

# ─── 2. Compute PTC status ──────────────────────────────────────────────────

cat("Computing PTC status...\n")

coding <- cds_meta %>%
  filter(coding_status == "coding") %>%
  inner_join(structs, by = "isoform_id")

compute_ptc_distance <- function(cds_start, cds_stop, strand,
                                  exon_starts, exon_ends, n_exons) {
  if (n_exons == 1) return(NA_real_)
  stop_pos <- if (strand == "+") cds_stop else cds_start
  if (strand == "+") { ord <- order(exon_starts) } else { ord <- order(exon_starts, decreasing = TRUE) }
  e_starts <- exon_starts[ord]; e_ends <- exon_ends[ord]
  exon_lengths <- e_ends - e_starts
  cum_starts <- c(0, cumsum(exon_lengths[-length(exon_lengths)]))
  last_ejc_mRNA <- sum(exon_lengths[-length(exon_lengths)])
  exon_idx <- which(e_starts <= stop_pos & stop_pos <= e_ends)
  if (length(exon_idx) == 0) return(NA_real_)
  exon_idx <- exon_idx[1]
  if (strand == "+") { offset_in_exon <- stop_pos - e_starts[exon_idx] }
  else { offset_in_exon <- e_ends[exon_idx] - stop_pos }
  stop_mRNA <- cum_starts[exon_idx] + offset_in_exon
  last_ejc_mRNA - stop_mRNA
}

ptc_results <- coding %>%
  filter(n_exons > 1) %>%
  rowwise() %>%
  mutate(
    ptc_distance = compute_ptc_distance(
      cds_start, cds_stop, as.character(strand),
      exon_starts, exon_ends, n_exons
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(ptc_distance)) %>%
  mutate(
    has_ptc = ptc_distance > 50,
    ptc_label = if_else(has_ptc, "PTC (>50nt)", "No PTC (≤50nt)")
  )

ptc_lookup <- ptc_results %>%
  select(isoform_id, ptc_distance, has_ptc, ptc_label)

cat(sprintf("  %d multi-exon coding isoforms with PTC status\n", nrow(ptc_lookup)))

# ─── 3. Compute mean DMSO CPM per isoform per cell type ─────────────────────

cat("Computing DMSO expression levels...\n")

dmso_cpm <- expression_data %>%
  filter(treatment == "DMSO") %>%
  group_by(isoform_id, ct) %>%
  summarise(mean_dmso_cpm = mean(cpm), .groups = "drop")

# ─── 4. Apply 5% filter per cell type ───────────────────────────────────────

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

# ─── 5. Build analysis data per cell type ────────────────────────────────────

cat("Building analysis datasets...\n\n")

all_data <- list()

for (ct in cell_types) {
  de_file <- file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix[ct], "_", de_date, ".csv"))
  de <- read_csv(de_file, show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC) %>%
    rename(isoform_id = txid)

  ct_samples <- sample_metadata %>%
    filter(ct == !!ct) %>%
    pull(sample_id)
  passing_5pct <- apply_5pct_filter(expression_data, ct_samples)

  ct_dmso <- dmso_cpm %>% filter(ct == !!ct) %>% select(isoform_id, mean_dmso_cpm)

  df <- ptc_lookup %>%
    filter(isoform_id %in% passing_5pct) %>%
    inner_join(de, by = "isoform_id") %>%
    left_join(ct_dmso, by = "isoform_id") %>%
    mutate(cell_type = ct)

  cat(sprintf("  %s: %d PTC+, %d PTC- (passing 5%% filter)\n",
              ct, sum(df$has_ptc), sum(!df$has_ptc)))

  all_data[[ct]] <- df
}

combined <- bind_rows(all_data)

# Create expression tertiles within each cell type
combined <- combined %>%
  group_by(cell_type) %>%
  mutate(
    dmso_cpm_tertile = case_when(
      is.na(mean_dmso_cpm) | mean_dmso_cpm == 0 ~ "Zero/NA",
      mean_dmso_cpm <= quantile(mean_dmso_cpm[mean_dmso_cpm > 0], 1/3, na.rm = TRUE) ~ "Low",
      mean_dmso_cpm <= quantile(mean_dmso_cpm[mean_dmso_cpm > 0], 2/3, na.rm = TRUE) ~ "Mid",
      TRUE ~ "High"
    )
  ) %>%
  ungroup()

# ─── 6. Plot 1: logFC distributions by PTC status, per cell type ────────────

cat("\n\nGenerating plots...\n")

output_dir <- "comparisons/nonNMD_0.50/cross_comparison/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Summary stats for annotation
summary_stats <- combined %>%
  group_by(cell_type, ptc_label) %>%
  summarise(
    n = n(),
    median_logFC = median(logFC),
    mean_logFC = mean(logFC),
    .groups = "drop"
  )

cat("\nlogFC summary by PTC status and cell type:\n")
summary_stats %>%
  mutate(across(c(median_logFC, mean_logFC), ~round(., 4))) %>%
  print(n = 20)

# Wilcoxon tests
cat("\nWilcoxon rank-sum tests (PTC+ vs PTC-):\n")
for (ct in cell_types) {
  d <- combined %>% filter(cell_type == ct)
  wt <- wilcox.test(logFC ~ has_ptc, data = d)
  cat(sprintf("  %s: p = %g\n", ct, wt$p.value))
}

# Plot
p1 <- ggplot(combined, aes(x = logFC, fill = ptc_label)) +
  geom_density(alpha = 0.4, linewidth = 0.3) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("PTC (>50nt)" = "#d62728", "No PTC (≤50nt)" = "#1f77b4")) +
  labs(
    title = "logFC Distribution: PTC vs No PTC (5% filter, multi-exon coding)",
    subtitle = "logFC = Smg1i vs DMSO; positive = upregulated under NMD inhibition",
    x = "logFC (Smg1i vs DMSO)",
    y = "Density",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))

ggsave(file.path(output_dir, "ptc_logfc_by_celltype.pdf"), p1, width = 10, height = 7)
cat(sprintf("Saved: %s\n", file.path(output_dir, "ptc_logfc_by_celltype.pdf")))

# ─── 7. Plot 2: Stratified by expression level ──────────────────────────────

# Exclude Zero/NA for cleaner plots
combined_expr <- combined %>%
  filter(dmso_cpm_tertile != "Zero/NA") %>%
  mutate(dmso_cpm_tertile = factor(dmso_cpm_tertile, levels = c("Low", "Mid", "High")))

# Summary stats stratified
summary_strat <- combined_expr %>%
  group_by(cell_type, dmso_cpm_tertile, ptc_label) %>%
  summarise(
    n = n(),
    median_logFC = median(logFC),
    mean_logFC = mean(logFC),
    .groups = "drop"
  )

cat("\nlogFC summary stratified by expression tertile:\n")
summary_strat %>%
  mutate(across(c(median_logFC, mean_logFC), ~round(., 4))) %>%
  print(n = 40)

# Wilcoxon tests stratified
cat("\nWilcoxon tests stratified by expression tertile:\n")
for (ct in cell_types) {
  for (tert in c("Low", "Mid", "High")) {
    d <- combined_expr %>% filter(cell_type == ct, dmso_cpm_tertile == tert)
    if (sum(d$has_ptc) >= 10 && sum(!d$has_ptc) >= 10) {
      wt <- wilcox.test(logFC ~ has_ptc, data = d)
      median_ptc <- median(d$logFC[d$has_ptc])
      median_no <- median(d$logFC[!d$has_ptc])
      cat(sprintf("  %s / %s: median PTC=%.3f, no PTC=%.3f, diff=%.3f, p=%g (n_ptc=%d, n_no=%d)\n",
                  ct, tert, median_ptc, median_no, median_ptc - median_no,
                  wt$p.value, sum(d$has_ptc), sum(!d$has_ptc)))
    }
  }
}

p2 <- ggplot(combined_expr, aes(x = logFC, fill = ptc_label)) +
  geom_density(alpha = 0.4, linewidth = 0.3) +
  facet_grid(dmso_cpm_tertile ~ cell_type, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("PTC (>50nt)" = "#d62728", "No PTC (≤50nt)" = "#1f77b4")) +
  labs(
    title = "logFC Distribution by PTC Status and DMSO Expression Level",
    subtitle = "Rows = DMSO expression tertile (Low/Mid/High); Columns = Cell type",
    x = "logFC (Smg1i vs DMSO)",
    y = "Density",
    fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))

ggsave(file.path(output_dir, "ptc_logfc_by_celltype_expression.pdf"), p2, width = 12, height = 8)
cat(sprintf("Saved: %s\n", file.path(output_dir, "ptc_logfc_by_celltype_expression.pdf")))

# ─── 8. Also show quantile-quantile style comparison ────────────────────────

# Mean logFC shift by PTC status × expression tertile × cell type
shift_summary <- combined_expr %>%
  group_by(cell_type, dmso_cpm_tertile, has_ptc) %>%
  summarise(
    n = n(),
    mean_logFC = mean(logFC),
    median_logFC = median(logFC),
    pct_positive = 100 * mean(logFC > 0),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = has_ptc,
    values_from = c(n, mean_logFC, median_logFC, pct_positive),
    names_sep = "_"
  ) %>%
  mutate(
    mean_shift = mean_logFC_TRUE - mean_logFC_FALSE,
    median_shift = median_logFC_TRUE - median_logFC_FALSE,
    pct_positive_shift = pct_positive_TRUE - pct_positive_FALSE
  )

cat("\n\nPTC vs non-PTC shift summary:\n")
shift_summary %>%
  select(cell_type, dmso_cpm_tertile,
         n_TRUE, n_FALSE,
         mean_shift, median_shift, pct_positive_shift) %>%
  mutate(across(c(mean_shift, median_shift), ~round(., 4)),
         pct_positive_shift = round(pct_positive_shift, 1)) %>%
  rename(n_ptc = n_TRUE, n_no_ptc = n_FALSE,
         mean_logFC_shift = mean_shift,
         median_logFC_shift = median_shift,
         pct_pos_shift = pct_positive_shift) %>%
  print(n = 30, width = 120)

write_tsv(shift_summary,
          file.path(output_dir, "ptc_logfc_shift_summary.tsv"))
cat(sprintf("\nSaved: %s\n", file.path(output_dir, "ptc_logfc_shift_summary.tsv")))

cat("\nDone.\n")
