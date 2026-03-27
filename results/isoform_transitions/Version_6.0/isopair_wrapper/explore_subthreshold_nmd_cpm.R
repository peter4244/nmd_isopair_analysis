#!/usr/bin/env Rscript
# ==============================================================================
# Explore: Sub-Threshold NMD via CPM ratio (SMG1i / DMSO)
#
# Same isoform sets as report Section 3c "Sub-Threshold NMD" analysis:
#   - Non-NMD isoforms from NMD-harboring genes vs non-NMD genes
#   - Dominant non-NMD isoform per gene (same as ORF landscape analysis)
#
# Instead of mashr logFC, report the ratio of:
#   mean CPM (SMG1i) / mean CPM (DMSO)
# where CPM is counts / library_size * 1e6 (NOT TMM-normalized).
#
# A ratio > 1 means the isoform is more abundant under SMG1i (rescued from NMD).
# A ratio = 1 means no change.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
fmt <- function(x) format(x, big.mark = ",")

cat("=== Sub-Threshold NMD: CPM Ratio Analysis ===\n\n")

# === Load data ===
count_file <- "/Users/petecastaldi/claude_projects/nmd/isocall/nmd_lungcells/results/call/nmd_isocall.count_matrix.txt"
cat("Loading count matrix...\n")
count_df <- read.csv(count_file, check.names = FALSE)
count_mat <- as.matrix(count_df[, -1])
rownames(count_mat) <- count_df$id

sample_meta <- readRDS("data_mashr/sample_metadata.rds")
nmd_class <- readRDS("data_mashr/nmd_classification.rds")
structures <- readRDS("data_mashr/structures.rds")
cds <- readRDS("data_mashr/cds.rds")
gene_map <- readRDS("data_mashr/gene_map.rds")

# Exclude DO 029T (same as 01_prepare_data_mashr.R)
exclude_samples <- c("Sample25_DO_029T_DMSO", "Sample26_DO_029T_Smg1i")
count_mat <- count_mat[, !colnames(count_mat) %in% exclude_samples]
sample_meta <- sample_meta[!sample_meta$sample_id %in% exclude_samples, ]

# === Compute raw CPM (NOT TMM) ===
cat("Computing raw CPM (no TMM normalization)...\n")
lib_sizes <- colSums(count_mat)
cpm_mat <- sweep(count_mat, 2, lib_sizes, "/") * 1e6

dmso_samples <- sample_meta$sample_id[sample_meta$treatment == "DMSO"]
smg1i_samples <- sample_meta$sample_id[sample_meta$treatment == "Smg1i"]

mean_cpm_dmso <- rowMeans(cpm_mat[, dmso_samples])
mean_cpm_smg1i <- rowMeans(cpm_mat[, smg1i_samples])

# === Identify the same isoform/gene sets as the report ===
nmd_ids <- nmd_class$all_samples$nmd
non_nmd_ids <- nmd_class$all_samples$non_nmd
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]

# NMD-harboring genes
nmd_genes <- unique(gene_map$gene_id[gene_map$isoform_id %in% nmd_ids])

# Dominant non-NMD coding isoform per gene (highest DMSO CPM)
non_nmd_coding <- intersect(non_nmd_ids, coding_ids)
non_nmd_in_cpm <- intersect(non_nmd_coding, rownames(cpm_mat))

dom_iso <- data.frame(
  isoform_id = non_nmd_in_cpm,
  gene_id = gene_map$gene_id[match(non_nmd_in_cpm, gene_map$isoform_id)],
  mean_dmso = mean_cpm_dmso[non_nmd_in_cpm],
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>%
  slice_max(mean_dmso, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Dominant non-NMD isoforms:", nrow(dom_iso), "\n")

# ORFik-based gene set (same filter as report: genes with ORFik data)
orfik_data <- readRDS("data_mashr/analysis_cache/orfik_scan.rds")
orfik_genes <- unique(gene_map$gene_id[gene_map$isoform_id %in% orfik_data$tx_summary$isoform_id])
dom_iso_orfik <- dom_iso %>% filter(gene_id %in% orfik_genes)
cat("Dominant isoforms in ORFik gene set:", nrow(dom_iso_orfik), "\n")

# === Compute CPM ratio for all non-NMD isoforms ===
cat("\nComputing CPM ratios...\n")

# All non-NMD isoforms (not just dominant)
non_nmd_expr <- data.frame(
  isoform_id = non_nmd_in_cpm,
  gene_id = gene_map$gene_id[match(non_nmd_in_cpm, gene_map$isoform_id)],
  cpm_dmso = mean_cpm_dmso[non_nmd_in_cpm],
  cpm_smg1i = mean_cpm_smg1i[non_nmd_in_cpm],
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(gene_id)) %>%
  # Restrict to the same gene set as report (ORFik genes)
  filter(gene_id %in% orfik_genes) %>%
  mutate(
    cpm_ratio = cpm_smg1i / pmax(cpm_dmso, 0.01),  # avoid division by zero
    log2_ratio = log2(cpm_ratio),
    gene_type = ifelse(gene_id %in% nmd_genes, "NMD-harboring gene", "Non-NMD gene")
  )

cat("Non-NMD isoforms with CPM ratio:", nrow(non_nmd_expr), "\n")
cat("  NMD-harboring genes:", sum(non_nmd_expr$gene_type == "NMD-harboring gene"), "\n")
cat("  Non-NMD genes:", sum(non_nmd_expr$gene_type == "Non-NMD gene"), "\n")

# === Per-isoform summary ===
cat("\n=== Per-Isoform CPM Ratio Summary ===\n")
for (gt in c("NMD-harboring gene", "Non-NMD gene")) {
  sub <- non_nmd_expr[non_nmd_expr$gene_type == gt, ]
  cat(sprintf("\n%s (n = %s isoforms):\n", gt, fmt(nrow(sub))))
  cat(sprintf("  Median CPM ratio (SMG1i/DMSO): %.3f\n", median(sub$cpm_ratio)))
  cat(sprintf("  Mean CPM ratio: %.3f\n", mean(sub$cpm_ratio)))
  cat(sprintf("  Median log2 ratio: %.3f\n", median(sub$log2_ratio)))
  cat(sprintf("  IQR log2 ratio: [%.3f, %.3f]\n",
              quantile(sub$log2_ratio, 0.25), quantile(sub$log2_ratio, 0.75)))
  cat(sprintf("  %% with ratio > 1: %.1f%%\n", 100 * mean(sub$cpm_ratio > 1)))
}

wt_iso <- wilcox.test(
  non_nmd_expr$log2_ratio[non_nmd_expr$gene_type == "NMD-harboring gene"],
  non_nmd_expr$log2_ratio[non_nmd_expr$gene_type == "Non-NMD gene"])
cat(sprintf("\nWilcoxon p (isoform-level): %s\n", signif(wt_iso$p.value, 3)))

# === Per-gene summary (average across non-NMD isoforms per gene) ===
cat("\n=== Per-Gene CPM Ratio Summary ===\n")
gene_ratio <- non_nmd_expr %>%
  group_by(gene_id, gene_type) %>%
  summarise(
    n_isoforms = n(),
    mean_cpm_ratio = mean(cpm_ratio),
    mean_log2_ratio = mean(log2_ratio),
    .groups = "drop"
  )

for (gt in c("NMD-harboring gene", "Non-NMD gene")) {
  sub <- gene_ratio[gene_ratio$gene_type == gt, ]
  cat(sprintf("\n%s (n = %s genes):\n", gt, fmt(nrow(sub))))
  cat(sprintf("  Median mean CPM ratio: %.3f\n", median(sub$mean_cpm_ratio)))
  cat(sprintf("  Median mean log2 ratio: %.3f\n", median(sub$mean_log2_ratio)))
  cat(sprintf("  IQR mean log2 ratio: [%.3f, %.3f]\n",
              quantile(sub$mean_log2_ratio, 0.25), quantile(sub$mean_log2_ratio, 0.75)))
  cat(sprintf("  %% genes with mean ratio > 1: %.1f%%\n", 100 * mean(sub$mean_cpm_ratio > 1)))
}

wt_gene <- wilcox.test(
  gene_ratio$mean_log2_ratio[gene_ratio$gene_type == "NMD-harboring gene"],
  gene_ratio$mean_log2_ratio[gene_ratio$gene_type == "Non-NMD gene"])
cat(sprintf("\nWilcoxon p (gene-level): %s\n", signif(wt_gene$p.value, 3)))

# === Figures ===
theme_set(theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank()))

# Isoform-level boxplot
fig1 <- ggplot(non_nmd_expr, aes(x = gene_type, y = log2_ratio, fill = gene_type)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("NMD-harboring gene" = "#ef8a62", "Non-NMD gene" = "#67a9cf")) +
  annotate("text", x = 1.5, y = max(non_nmd_expr$log2_ratio, na.rm = TRUE) * 0.8,
           label = sprintf("p = %s", signif(wt_iso$p.value, 3)), size = 4.5) +
  labs(title = "Non-NMD Isoform CPM Ratio by Gene NMD Status",
       subtitle = sprintf("log2(mean SMG1i CPM / mean DMSO CPM), raw CPM (no TMM) | %s isoforms",
                          fmt(nrow(non_nmd_expr))),
       x = NULL, y = "log2(SMG1i CPM / DMSO CPM)") +
  theme(legend.position = "none")

ggsave("figures/explore_subthreshold_isoform_cpm_ratio.pdf", fig1, width = 10, height = 6)
ggsave("figures/explore_subthreshold_isoform_cpm_ratio.png", fig1, width = 10, height = 6, dpi = 300)

# Gene-level boxplot
fig2 <- ggplot(gene_ratio, aes(x = gene_type, y = mean_log2_ratio, fill = gene_type)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("NMD-harboring gene" = "#ef8a62", "Non-NMD gene" = "#67a9cf")) +
  annotate("text", x = 1.5, y = max(gene_ratio$mean_log2_ratio, na.rm = TRUE) * 0.8,
           label = sprintf("p = %s", signif(wt_gene$p.value, 3)), size = 4.5) +
  labs(title = "Non-NMD Isoform CPM Ratio by Gene NMD Status (Gene-Level)",
       subtitle = sprintf("Mean log2(SMG1i/DMSO) per gene, raw CPM (no TMM) | %s genes",
                          fmt(nrow(gene_ratio))),
       x = NULL, y = "Mean log2(SMG1i CPM / DMSO CPM)") +
  theme(legend.position = "none")

ggsave("figures/explore_subthreshold_gene_cpm_ratio.pdf", fig2, width = 10, height = 6)
ggsave("figures/explore_subthreshold_gene_cpm_ratio.png", fig2, width = 10, height = 6, dpi = 300)

cat("\nFigures saved to figures/explore_subthreshold_*.{pdf,png}\n")
cat("Done.\n")
