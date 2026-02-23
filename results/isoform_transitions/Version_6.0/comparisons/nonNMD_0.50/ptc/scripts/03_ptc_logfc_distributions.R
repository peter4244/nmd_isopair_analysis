#!/usr/bin/env Rscript
# 03_ptc_logfc_distributions.R
#
# Show distribution of logFC (beta-coefficients) for PTC+ vs PTC- isoforms,
# stratified by cell type and DMSO expression level. Produces both 5%-filtered
# and unfiltered versions, plus GENCODE vs PacBio stratification.
#
# Requires: ptc_status.rds from 01_compute_ptc_status.R
#
# Outputs (figures/):
#   ptc_logfc_by_celltype.pdf                          — 5%-filtered, by cell type
#   ptc_logfc_by_celltype_expression.pdf               — 5%-filtered, by cell type x expression
#   ptc_logfc_by_celltype_unfiltered.pdf               — unfiltered, by cell type
#   ptc_logfc_by_celltype_expression_unfiltered.pdf    — unfiltered, by cell type x expression
#   ptc_logfc_by_source.pdf                            — GENCODE vs PacBio, 5%-filtered
#   ptc_logfc_by_source_unfiltered.pdf                 — GENCODE vs PacBio, unfiltered
# Outputs (results/):
#   ptc_logfc_shift_summary.tsv
#
# Run from: results/isoform_transitions/Version_6.0/

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

cat("\n")
cat("==================================================================\n")
cat("   03: PTC logFC Distributions\n")
cat("==================================================================\n\n")

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

ptc_lookup <- ptc_status %>%
  filter(!is.na(has_ptc), n_exons > 1) %>%
  mutate(ptc_label = if_else(has_ptc, "PTC (>50nt)", "No PTC (<=50nt)"))

fig_dir <- "comparisons/nonNMD_0.50/ptc/figures"
res_dir <- "comparisons/nonNMD_0.50/ptc/results"

# ─── 2. 5% filter ────────────────────────────────────────────────────────────

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

# ─── 3. Compute DMSO CPM ────────────────────────────────────────────────────

cat("Computing DMSO expression levels...\n")
dmso_cpm <- expression_data %>%
  filter(treatment == "DMSO") %>%
  group_by(isoform_id, ct) %>%
  summarise(mean_dmso_cpm = mean(cpm), .groups = "drop")

# ─── 4. Build datasets per cell type (filtered + unfiltered) ─────────────────

cat("Building analysis datasets...\n\n")

data_filt <- list()
data_unfilt <- list()

for (ct in cell_types) {
  de_file <- file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix[ct], "_", de_date, ".csv"))
  de <- read_csv(de_file, show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC) %>%
    rename(isoform_id = txid)

  ct_samples <- sample_metadata %>% filter(ct == !!ct) %>% pull(sample_id)
  passing_5pct <- apply_5pct_filter(expression_data, ct_samples)
  ct_dmso <- dmso_cpm %>% filter(ct == !!ct) %>% select(isoform_id, mean_dmso_cpm)

  # Unfiltered
  df_u <- ptc_lookup %>%
    inner_join(de, by = "isoform_id") %>%
    left_join(ct_dmso, by = "isoform_id") %>%
    mutate(cell_type = ct, mean_dmso_cpm = replace_na(mean_dmso_cpm, 0))
  data_unfilt[[ct]] <- df_u

  # 5%-filtered
  df_f <- df_u %>% filter(isoform_id %in% passing_5pct)
  data_filt[[ct]] <- df_f

  cat(sprintf("  %s: %d filtered / %d unfiltered\n", ct, nrow(df_f), nrow(df_u)))
}

combined_filt <- bind_rows(data_filt)
combined_unfilt <- bind_rows(data_unfilt)

# Add expression tertiles
add_tertiles <- function(df, include_zero = FALSE) {
  df %>%
    group_by(cell_type) %>%
    mutate(
      dmso_cpm_tertile = if (include_zero) {
        case_when(
          mean_dmso_cpm == 0 ~ "Zero",
          mean_dmso_cpm <= quantile(mean_dmso_cpm[mean_dmso_cpm > 0], 1/3) ~ "Low",
          mean_dmso_cpm <= quantile(mean_dmso_cpm[mean_dmso_cpm > 0], 2/3) ~ "Mid",
          TRUE ~ "High"
        )
      } else {
        case_when(
          is.na(mean_dmso_cpm) | mean_dmso_cpm == 0 ~ "Zero/NA",
          mean_dmso_cpm <= quantile(mean_dmso_cpm[mean_dmso_cpm > 0], 1/3, na.rm = TRUE) ~ "Low",
          mean_dmso_cpm <= quantile(mean_dmso_cpm[mean_dmso_cpm > 0], 2/3, na.rm = TRUE) ~ "Mid",
          TRUE ~ "High"
        )
      }
    ) %>%
    ungroup()
}

combined_filt <- add_tertiles(combined_filt, include_zero = FALSE)
combined_unfilt <- add_tertiles(combined_unfilt, include_zero = TRUE)

# ─── Helper: compute Wilcoxon p-values per cell type ─────────────────────────

get_wilcox_labels <- function(df) {
  df %>%
    group_by(cell_type) %>%
    summarise(
      p = wilcox.test(logFC ~ has_ptc)$p.value,
      n_ptc = sum(has_ptc), n_no = sum(!has_ptc),
      .groups = "drop"
    ) %>%
    mutate(
      sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"),
      label = sprintf("%s (p=%s)", cell_type,
                       ifelse(p < 0.001, format(p, digits = 2, scientific = TRUE),
                              format(round(p, 3), nsmall = 3))),
      is_sig = p < 0.05
    )
}

# ═══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\nGenerating plots...\n")

# ─── Plot 1: 5%-filtered, by cell type ───────────────────────────────────────

p1 <- ggplot(combined_filt, aes(x = logFC, fill = ptc_label)) +
  geom_density(alpha = 0.4, linewidth = 0.3) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("PTC (>50nt)" = "#d62728", "No PTC (<=50nt)" = "#1f77b4")) +
  labs(
    title = "logFC Distribution: PTC vs No PTC (5% filter, multi-exon coding)",
    subtitle = "logFC = Smg1i vs DMSO; positive = upregulated under NMD inhibition",
    x = "logFC (Smg1i vs DMSO)", y = "Density", fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "ptc_logfc_by_celltype.pdf"), p1, width = 10, height = 7)
cat(sprintf("  Saved: ptc_logfc_by_celltype.pdf\n"))

# ─── Plot 2: 5%-filtered, by cell type x expression ──────────────────────────

combined_filt_expr <- combined_filt %>%
  filter(dmso_cpm_tertile != "Zero/NA") %>%
  mutate(dmso_cpm_tertile = factor(dmso_cpm_tertile, levels = c("Low", "Mid", "High")))

p2 <- ggplot(combined_filt_expr, aes(x = logFC, fill = ptc_label)) +
  geom_density(alpha = 0.4, linewidth = 0.3) +
  facet_grid(dmso_cpm_tertile ~ cell_type, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("PTC (>50nt)" = "#d62728", "No PTC (<=50nt)" = "#1f77b4")) +
  labs(
    title = "logFC Distribution by PTC Status and DMSO Expression Level",
    subtitle = "Rows = DMSO expression tertile (Low/Mid/High); Columns = Cell type",
    x = "logFC (Smg1i vs DMSO)", y = "Density", fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "ptc_logfc_by_celltype_expression.pdf"), p2, width = 12, height = 8)
cat(sprintf("  Saved: ptc_logfc_by_celltype_expression.pdf\n"))

# ─── Plot 3: Unfiltered, by cell type ────────────────────────────────────────

wilcox_labels <- get_wilcox_labels(combined_unfilt)
ct_labels <- setNames(wilcox_labels$label, wilcox_labels$cell_type)
combined_unfilt$cell_type_label <- factor(
  ct_labels[combined_unfilt$cell_type], levels = ct_labels[cell_types]
)

sig_cts <- wilcox_labels$cell_type[wilcox_labels$is_sig]
sig_bg <- combined_unfilt %>%
  distinct(cell_type_label, cell_type) %>%
  mutate(is_sig = cell_type %in% sig_cts)

p3 <- ggplot(combined_unfilt, aes(x = logFC, color = ptc_label)) +
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
    subtitle = "Yellow background = significant (Wilcoxon p < 0.05)",
    x = "logFC (Smg1i vs DMSO)", y = "Density", color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "ptc_logfc_by_celltype_unfiltered.pdf"), p3, width = 10, height = 7)
cat(sprintf("  Saved: ptc_logfc_by_celltype_unfiltered.pdf\n"))

# ─── Plot 4: Unfiltered, by cell type x expression ──────────────────────────

combined_unfilt_expr <- combined_unfilt %>%
  mutate(dmso_cpm_tertile = factor(dmso_cpm_tertile,
                                    levels = c("Zero", "Low", "Mid", "High")))

wilcox_strat <- combined_unfilt_expr %>%
  group_by(cell_type, dmso_cpm_tertile) %>%
  summarise(
    p = tryCatch(wilcox.test(logFC ~ has_ptc)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(is_sig = !is.na(p) & p < 0.05)

sig_panels <- wilcox_strat %>% filter(is_sig)
p_annotations <- wilcox_strat %>%
  filter(!is.na(p)) %>%
  mutate(sig_star = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ ""))

p4 <- ggplot(combined_unfilt_expr, aes(x = logFC, color = ptc_label)) +
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
    title = "logFC Distribution by PTC Status and DMSO Expression (NO 5% filter)",
    subtitle = "Yellow / stars = significant (Wilcoxon p < 0.05)",
    x = "logFC (Smg1i vs DMSO)", y = "Density", color = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "ptc_logfc_by_celltype_expression_unfiltered.pdf"), p4,
       width = 12, height = 9)
cat(sprintf("  Saved: ptc_logfc_by_celltype_expression_unfiltered.pdf\n"))

# ─── Plot 5: GENCODE vs PacBio, 5%-filtered ──────────────────────────────────

combined_filt_source <- combined_filt %>%
  mutate(source = case_when(is_gencode ~ "GENCODE", is_novel ~ "PacBio", TRUE ~ NA_character_)) %>%
  filter(!is.na(source))

p5 <- ggplot(combined_filt_source, aes(x = logFC, fill = ptc_label)) +
  geom_density(alpha = 0.4, linewidth = 0.3) +
  facet_grid(source ~ cell_type, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("PTC (>50nt)" = "#d62728", "No PTC (<=50nt)" = "#1f77b4")) +
  labs(
    title = "logFC Distribution by PTC Status: GENCODE vs PacBio (5% filter)",
    subtitle = "Rows = isoform source; Columns = Cell type",
    x = "logFC (Smg1i vs DMSO)", y = "Density", fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "ptc_logfc_by_source.pdf"), p5, width = 12, height = 6)
cat(sprintf("  Saved: ptc_logfc_by_source.pdf\n"))

# ─── Plot 6: GENCODE vs PacBio, unfiltered ──────────────────────────────────

combined_unfilt_source <- combined_unfilt %>%
  mutate(source = case_when(is_gencode ~ "GENCODE", is_novel ~ "PacBio", TRUE ~ NA_character_)) %>%
  filter(!is.na(source))

# Wilcoxon tests per source x cell_type
wilcox_source <- combined_unfilt_source %>%
  group_by(cell_type, source) %>%
  summarise(
    p = tryCatch(wilcox.test(logFC ~ has_ptc)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    is_sig = !is.na(p) & p < 0.05,
    sig_star = case_when(
      is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ ""
    )
  )

sig_source <- wilcox_source %>% filter(is_sig)

p6 <- ggplot(combined_unfilt_source, aes(x = logFC, color = ptc_label)) +
  geom_rect(data = sig_source,
            aes(fill = is_sig), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            alpha = 0.07, inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "gold"), guide = "none") +
  geom_density(linewidth = 0.5, fill = NA) +
  facet_grid(source ~ cell_type, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_text(data = wilcox_source,
            aes(label = sig_star),
            x = Inf, y = Inf, hjust = 1.2, vjust = 1.5,
            size = 5, color = "black", fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_manual(values = c("PTC (>50nt)" = "#d62728", "No PTC (<=50nt)" = "#1f77b4")) +
  labs(
    title = "logFC Distribution by PTC Status: GENCODE vs PacBio (NO 5% filter)",
    subtitle = "Yellow / stars = significant (Wilcoxon p < 0.05)",
    x = "logFC (Smg1i vs DMSO)", y = "Density", color = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "ptc_logfc_by_source_unfiltered.pdf"), p6, width = 12, height = 6)
cat(sprintf("  Saved: ptc_logfc_by_source_unfiltered.pdf\n"))

# ─── Shift summary statistics ────────────────────────────────────────────────

cat("\nComputing shift summary...\n")

# 5%-filtered, expression-stratified
shift_summary <- combined_filt %>%
  filter(dmso_cpm_tertile != "Zero/NA") %>%
  mutate(dmso_cpm_tertile = factor(dmso_cpm_tertile, levels = c("Low", "Mid", "High"))) %>%
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

shift_summary %>%
  select(cell_type, dmso_cpm_tertile, n_TRUE, n_FALSE,
         mean_shift, median_shift, pct_positive_shift) %>%
  mutate(across(c(mean_shift, median_shift), ~round(., 4)),
         pct_positive_shift = round(pct_positive_shift, 1)) %>%
  rename(n_ptc = n_TRUE, n_no_ptc = n_FALSE) %>%
  print(n = 30, width = 120)

write_tsv(shift_summary, file.path(res_dir, "ptc_logfc_shift_summary.tsv"))
cat(sprintf("  Saved: ptc_logfc_shift_summary.tsv\n"))

# ─── Print Wilcoxon summary ──────────────────────────────────────────────────

cat("\nWilcoxon rank-sum tests (PTC+ vs PTC-):\n")
cat("\n  5%-filtered:\n")
for (ct in cell_types) {
  d <- combined_filt %>% filter(cell_type == ct)
  wt <- wilcox.test(logFC ~ has_ptc, data = d)
  cat(sprintf("    %s: p = %g\n", ct, wt$p.value))
}

cat("\n  Unfiltered:\n")
for (ct in cell_types) {
  d <- combined_unfilt %>% filter(cell_type == ct)
  wt <- wilcox.test(logFC ~ has_ptc, data = d)
  cat(sprintf("    %s: p = %g\n", ct, wt$p.value))
}

cat("\n  By source (unfiltered):\n")
for (src in c("GENCODE", "PacBio")) {
  for (ct in cell_types) {
    d <- combined_unfilt_source %>% filter(cell_type == ct, source == src)
    if (sum(d$has_ptc) >= 10 && sum(!d$has_ptc) >= 10) {
      wt <- wilcox.test(logFC ~ has_ptc, data = d)
      if (wt$p.value < 0.05) {
        shift <- median(d$logFC[d$has_ptc]) - median(d$logFC[!d$has_ptc])
        cat(sprintf("    %s / %s: median shift = %.3f, p = %g *\n", src, ct, shift, wt$p.value))
      }
    }
  }
}

cat("\nDone.\n")
