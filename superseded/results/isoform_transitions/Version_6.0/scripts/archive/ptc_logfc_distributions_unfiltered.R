#!/usr/bin/env Rscript
# PTC logFC Distributions — No 5% filter
#
# Same as ptc_logfc_distributions.R but using the complete DE results
# without applying the 5% expression proportion filter.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

cat("\n")
cat("==================================================================\n")
cat("   PTC logFC Distributions (No 5% Filter)\n")
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
    ptc_label = if_else(has_ptc, "PTC (>50nt)", "No PTC (<=50nt)")
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

# ─── 4. Build analysis data per cell type (NO 5% filter) ────────────────────

cat("Building analysis datasets (no 5%% filter)...\n\n")

all_data <- list()

for (ct in cell_types) {
  de_file <- file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix[ct], "_", de_date, ".csv"))
  de <- read_csv(de_file, show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC) %>%
    rename(isoform_id = txid)

  ct_dmso <- dmso_cpm %>% filter(ct == !!ct) %>% select(isoform_id, mean_dmso_cpm)

  # Join PTC status with DE — no expression filter
  df <- ptc_lookup %>%
    inner_join(de, by = "isoform_id") %>%
    left_join(ct_dmso, by = "isoform_id") %>%
    mutate(
      cell_type = ct,
      mean_dmso_cpm = replace_na(mean_dmso_cpm, 0)
    )

  cat(sprintf("  %s: %d PTC+, %d PTC- (all DE results)\n",
              ct, sum(df$has_ptc), sum(!df$has_ptc)))

  all_data[[ct]] <- df
}

combined <- bind_rows(all_data)

# Create expression tertiles within each cell type (among non-zero DMSO)
combined <- combined %>%
  group_by(cell_type) %>%
  mutate(
    dmso_cpm_tertile = case_when(
      mean_dmso_cpm == 0 ~ "Zero",
      mean_dmso_cpm <= quantile(mean_dmso_cpm[mean_dmso_cpm > 0], 1/3) ~ "Low",
      mean_dmso_cpm <= quantile(mean_dmso_cpm[mean_dmso_cpm > 0], 2/3) ~ "Mid",
      TRUE ~ "High"
    )
  ) %>%
  ungroup()

# ─── 5. Plot 1: logFC distributions by PTC status, per cell type ────────────

cat("\n\nGenerating plots...\n")

output_dir <- "comparisons/nonNMD_0.50/cross_comparison/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summary_stats <- combined %>%
  group_by(cell_type, ptc_label) %>%
  summarise(
    n = n(),
    median_logFC = median(logFC),
    mean_logFC = mean(logFC),
    .groups = "drop"
  )

cat("\nlogFC summary by PTC status and cell type (unfiltered):\n")
summary_stats %>%
  mutate(across(c(median_logFC, mean_logFC), ~round(., 4))) %>%
  print(n = 20)

cat("\nWilcoxon rank-sum tests (PTC+ vs PTC-):\n")
for (ct in cell_types) {
  d <- combined %>% filter(cell_type == ct)
  wt <- wilcox.test(logFC ~ has_ptc, data = d)
  cat(sprintf("  %s: p = %g (n_ptc=%d, n_no=%d)\n", ct, wt$p.value,
              sum(d$has_ptc), sum(!d$has_ptc)))
}

# Compute Wilcoxon p-values for facet labels
wilcox_labels <- combined %>%
  group_by(cell_type) %>%
  summarise(
    p = wilcox.test(logFC ~ has_ptc)$p.value,
    n_ptc = sum(has_ptc),
    n_no = sum(!has_ptc),
    .groups = "drop"
  ) %>%
  mutate(
    sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"),
    label = sprintf("%s (p=%s)", cell_type, ifelse(p < 0.001, format(p, digits = 2, scientific = TRUE), format(round(p, 3), nsmall = 3))),
    is_sig = p < 0.05
  )

# Create named vector for facet labeling
ct_labels <- setNames(wilcox_labels$label, wilcox_labels$cell_type)
combined$cell_type_label <- ct_labels[combined$cell_type]
# Preserve ordering
combined$cell_type_label <- factor(combined$cell_type_label, levels = ct_labels[cell_types])

# Significance indicator for shading
sig_cts <- wilcox_labels$cell_type[wilcox_labels$is_sig]
sig_bg <- combined %>%
  distinct(cell_type_label, cell_type) %>%
  mutate(is_sig = cell_type %in% sig_cts)

p1 <- ggplot(combined, aes(x = logFC, color = ptc_label)) +
  geom_rect(data = sig_bg %>% filter(is_sig),
            aes(fill = is_sig), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            alpha = 0.07, inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "gold"), guide = "none") +
  geom_density(linewidth = 0.7, fill = NA) +
  facet_wrap(~cell_type_label, scales = "free_y", ncol = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("PTC (>50nt)" = "#d62728", "No PTC (<=50nt)" = "#1f77b4")) +
  labs(
    title = "logFC Distribution: PTC vs No PTC (NO 5% filter, multi-exon coding)",
    subtitle = "Density curves only; yellow background = significant (Wilcoxon p < 0.05)",
    x = "logFC (Smg1i vs DMSO)",
    y = "Density",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))

ggsave(file.path(output_dir, "ptc_logfc_by_celltype_unfiltered.pdf"), p1, width = 10, height = 7)
cat(sprintf("\nSaved: %s\n", file.path(output_dir, "ptc_logfc_by_celltype_unfiltered.pdf")))

# ─── 6. Plot 2: Stratified by expression level ──────────────────────────────

combined_expr <- combined %>%
  mutate(dmso_cpm_tertile = factor(dmso_cpm_tertile,
                                    levels = c("Zero", "Low", "Mid", "High")))

# Summary stats stratified
cat("\nlogFC summary stratified by expression tertile (unfiltered):\n")
summary_strat <- combined_expr %>%
  group_by(cell_type, dmso_cpm_tertile, ptc_label) %>%
  summarise(
    n = n(),
    median_logFC = median(logFC),
    mean_logFC = mean(logFC),
    .groups = "drop"
  )

summary_strat %>%
  mutate(across(c(median_logFC, mean_logFC), ~round(., 4))) %>%
  print(n = 50)

# Wilcoxon tests stratified
cat("\nWilcoxon tests stratified by expression tertile:\n")
for (ct in cell_types) {
  for (tert in c("Zero", "Low", "Mid", "High")) {
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

# Compute Wilcoxon p-values per cell_type × expression tertile
wilcox_strat <- combined_expr %>%
  group_by(cell_type, dmso_cpm_tertile) %>%
  summarise(
    p = tryCatch(wilcox.test(logFC ~ has_ptc)$p.value, error = function(e) NA_real_),
    n_ptc = sum(has_ptc),
    n_no = sum(!has_ptc),
    .groups = "drop"
  ) %>%
  mutate(is_sig = !is.na(p) & p < 0.05)

# Create background rectangles for significant panels
sig_panels <- wilcox_strat %>%
  filter(is_sig) %>%
  mutate(
    p_label = ifelse(p < 0.001, format(p, digits = 2, scientific = TRUE), format(round(p, 3), nsmall = 3))
  )

# Add p-value annotations
p_annotations <- wilcox_strat %>%
  filter(!is.na(p)) %>%
  mutate(
    p_label = ifelse(p < 0.001, format(p, digits = 2, scientific = TRUE), format(round(p, 3), nsmall = 3)),
    sig_star = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
  )

p2 <- ggplot(combined_expr, aes(x = logFC, color = ptc_label)) +
  geom_rect(data = sig_panels,
            aes(fill = is_sig), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            alpha = 0.07, inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "gold"), guide = "none") +
  geom_density(linewidth = 0.5, fill = NA) +
  facet_grid(dmso_cpm_tertile ~ cell_type, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_text(data = p_annotations,
            aes(label = sig_star),
            x = Inf, y = Inf, hjust = 1.2, vjust = 1.5,
            size = 5, color = "black", fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_manual(values = c("PTC (>50nt)" = "#d62728", "No PTC (<=50nt)" = "#1f77b4")) +
  labs(
    title = "logFC Distribution by PTC Status and DMSO Expression Level (NO 5% filter)",
    subtitle = "Density curves; yellow background / stars = significant (Wilcoxon p < 0.05); * p<.05, ** p<.01, *** p<.001",
    x = "logFC (Smg1i vs DMSO)",
    y = "Density",
    color = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))

ggsave(file.path(output_dir, "ptc_logfc_by_celltype_expression_unfiltered.pdf"), p2, width = 12, height = 9)
cat(sprintf("Saved: %s\n", file.path(output_dir, "ptc_logfc_by_celltype_expression_unfiltered.pdf")))

cat("\nDone.\n")
