#!/usr/bin/env Rscript
# Paper Figure: Short-Read vs Long-Read Concordance (Manuscript Section 1)
# Panels: A) Per-cell-type SR-vs-LR expression correlation (gene-level log2 CPM+1)
#         B) Per-cell-type SR-vs-LR logFC concordance (mashr posterior, gene-level)
#
# Notes on data construction:
#  - SR is gene-level natively (DGEList counts -> CPM -> mean across DMSO per CT).
#  - LR has only isoform-level mashr; gene aggregation is performed here:
#      Expression: sum isoform CPMs per gene_id (gene CPM = sum-of-isoform-CPM).
#      logFC:      median of isoform mashr posterior logFC per gene_id.
#  - Gene IDs are matched on the unversioned ENSG (SR DGEList vs LR gene_map use
#    different version suffixes).
#  - SR DGEList CT label "AT2" corresponds to mashr column "AT" (paper scope).

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(edgeR)
})

paper_theme <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        plot.margin = margin(5, 8, 5, 5),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold", size = 11))

# Cell types in paper scope; SR DGEList uses "AT2" but both mashr files use "AT"
ct_paper <- c("AT", "DD", "FB", "MV")
sr_ct_map <- c(AT = "AT2", DD = "DD", FB = "FB", MV = "MV")

# ============================================================================
# Load data
# ============================================================================

# SR gene-level DGEList
sr_dge <- readRDS("/Users/petecastaldi/claude_projects/nmd/rds/dge_gene_2026.1.3.rds")
sr_cpm <- edgeR::cpm(sr_dge, normalized.lib.sizes = TRUE)
sr_meta <- sr_dge$samples
sr_meta$ct <- as.character(sr_meta$ct)
sr_meta$treatment <- as.character(sr_meta$treatment)

# LR transcript-level CPM matrix and metadata
lr_cpm_iso <- readRDS("data_mashr/expression_data.rds")
lr_meta <- readRDS("data_mashr/sample_metadata.rds")
lr_meta$ct <- as.character(lr_meta$ct)
lr_meta$treatment <- as.character(lr_meta$treatment)
gene_map <- readRDS("data_mashr/gene_map.rds")  # isoform_id -> gene_id

# SR mashr posterior means (gene-level)
sr_pm <- read.csv("/Users/petecastaldi/claude_projects/nmd/shortread_dge/mashr/mashr_posterior_means_2026.3.10.csv",
                  stringsAsFactors = FALSE)
# LR mashr posterior means (isoform-level)
lr_pm <- read.csv("/Users/petecastaldi/claude_projects/nmd/isocall_dge/mashr/mashr_isoform_posterior_means_2026.3.10.csv",
                  stringsAsFactors = FALSE)

# ============================================================================
# Aggregate LR isoform -> gene
# ============================================================================
# Helper to strip ENSG version suffix (handles novel suffix too)
strip_ver <- function(x) sub("[.][0-9]+$", "", x)

# Build LR gene-level CPM matrix by summing isoform CPMs within each gene
iso_to_gene <- setNames(gene_map$gene_id, gene_map$isoform_id)
lr_iso_ids <- rownames(lr_cpm_iso)
lr_gene_for_iso <- iso_to_gene[lr_iso_ids]
keep <- !is.na(lr_gene_for_iso)
if (sum(!keep) > 0) {
  message(sprintf("Dropping %d isoforms not in gene_map", sum(!keep)))
}
lr_cpm_iso_kept <- lr_cpm_iso[keep, , drop = FALSE]
lr_gene_for_iso_kept <- lr_gene_for_iso[keep]

# Aggregate (sum) isoform CPMs to gene level
lr_cpm_gene <- rowsum(lr_cpm_iso_kept, group = lr_gene_for_iso_kept, reorder = TRUE)
rownames_unv_lr <- strip_ver(rownames(lr_cpm_gene))

# ============================================================================
# Panel A: Per-CT expression correlation (DMSO mean log2(CPM+1))
# ============================================================================

build_expr_pair <- function(ct) {
  sr_ct <- sr_ct_map[[ct]]
  # SR samples: this CT, DMSO
  sr_keep <- sr_meta$ct == sr_ct & sr_meta$treatment == "DMSO"
  sr_mean <- rowMeans(sr_cpm[, sr_keep, drop = FALSE])
  sr_df <- data.frame(gene_unv = strip_ver(names(sr_mean)),
                      sr_log2cpm = log2(sr_mean + 1),
                      stringsAsFactors = FALSE)
  # If duplicates after stripping versions, keep highest-expressed
  sr_df <- sr_df[order(-sr_df$sr_log2cpm), ]
  sr_df <- sr_df[!duplicated(sr_df$gene_unv), ]

  # LR samples: this CT, DMSO
  lr_keep <- lr_meta$ct == ct & lr_meta$treatment == "DMSO"
  lr_sample_ids <- lr_meta$sample_id[lr_keep]
  lr_sample_ids <- intersect(lr_sample_ids, colnames(lr_cpm_gene))
  lr_mean <- rowMeans(lr_cpm_gene[, lr_sample_ids, drop = FALSE])
  lr_df <- data.frame(gene_unv = strip_ver(names(lr_mean)),
                      lr_log2cpm = log2(lr_mean + 1),
                      stringsAsFactors = FALSE)
  lr_df <- lr_df[order(-lr_df$lr_log2cpm), ]
  lr_df <- lr_df[!duplicated(lr_df$gene_unv), ]

  m <- merge(sr_df, lr_df, by = "gene_unv")
  # Drop genes with zero on both platforms (uninformative)
  m <- m[m$sr_log2cpm > 0 | m$lr_log2cpm > 0, ]
  m$ct <- ct
  m
}

expr_list <- lapply(ct_paper, build_expr_pair)
expr_df <- do.call(rbind, expr_list)

# Compute per-CT Pearson r
expr_stats <- expr_df %>%
  group_by(ct) %>%
  summarise(
    n = n(),
    r = cor(sr_log2cpm, lr_log2cpm, use = "complete.obs"),
    .groups = "drop"
  )
expr_stats$label <- sprintf("r = %.2f, n = %s",
                            expr_stats$r,
                            format(expr_stats$n, big.mark = ","))
# Match facet ordering
expr_df$ct <- factor(expr_df$ct, levels = ct_paper)
expr_stats$ct <- factor(expr_stats$ct, levels = ct_paper)

# Common axis range so facets are comparable
xrange <- range(c(expr_df$sr_log2cpm, expr_df$lr_log2cpm), finite = TRUE)
xpad <- diff(xrange) * 0.02
xlim_expr <- c(xrange[1] - xpad, xrange[2] + xpad)
label_x <- xlim_expr[1] + 0.05 * diff(xlim_expr)
label_y <- xlim_expr[2] - 0.05 * diff(xlim_expr)

pA <- ggplot(expr_df, aes(x = sr_log2cpm, y = lr_log2cpm)) +
  geom_point(alpha = 0.08, size = 0.4, color = "grey20") +
  geom_abline(slope = 1, intercept = 0, color = "red",
              linetype = "dashed", linewidth = 0.4) +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue",
              linewidth = 0.6) +
  geom_text(data = expr_stats, aes(x = label_x, y = label_y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.4,
            fontface = "bold", color = "black") +
  facet_wrap(~ ct, nrow = 1) +
  coord_cartesian(xlim = xlim_expr, ylim = xlim_expr) +
  labs(title = "Per-cell-type expression correlation: SR vs LR (gene-level)",
       subtitle = "DMSO mean log2(CPM+1); LR gene CPM = sum of isoform CPMs; red dashed = y=x",
       x = "SR gene log2(CPM+1)",
       y = "LR gene log2(CPM+1)") +
  paper_theme

# ============================================================================
# Panel B: Per-CT logFC concordance (gene-level mashr posterior)
# ============================================================================

# LR isoform -> gene logFC: median across isoforms per gene per CT
lr_pm$gene_unv <- strip_ver(lr_pm$gene_id)
sr_pm$gene_unv <- strip_ver(sr_pm$gene_id)

logfc_long_list <- list()
for (ct in ct_paper) {
  sr_col <- paste0("Smg1i_in_", ct)
  lr_col <- paste0("Smg1i_in_", ct)
  # SR: deduplicate to one row per unversioned gene (take row with max abs logFC)
  sr_sub <- sr_pm[, c("gene_unv", sr_col)]
  names(sr_sub)[2] <- "sr_logFC"
  sr_sub <- sr_sub[order(-abs(sr_sub$sr_logFC)), ]
  sr_sub <- sr_sub[!duplicated(sr_sub$gene_unv), ]

  # LR: median across isoforms per gene
  lr_sub <- lr_pm %>%
    group_by(gene_unv) %>%
    summarise(lr_logFC = median(.data[[lr_col]], na.rm = TRUE),
              n_iso = dplyr::n(),
              .groups = "drop")

  m <- merge(sr_sub, lr_sub, by = "gene_unv")
  m$ct <- ct
  logfc_long_list[[ct]] <- m
}
logfc_df <- do.call(rbind, logfc_long_list)

logfc_stats <- logfc_df %>%
  group_by(ct) %>%
  summarise(
    n = n(),
    r = cor(sr_logFC, lr_logFC, use = "complete.obs"),
    .groups = "drop"
  )
logfc_stats$label <- sprintf("r = %.2f, n = %s",
                              logfc_stats$r,
                              format(logfc_stats$n, big.mark = ","))
logfc_df$ct <- factor(logfc_df$ct, levels = ct_paper)
logfc_stats$ct <- factor(logfc_stats$ct, levels = ct_paper)

# Clip extreme outliers for plotting (1st-99th percentile within each axis), but
# DISCLOSE clipped count and leave actual values used for r computation untouched.
logfc_q <- quantile(c(logfc_df$sr_logFC, logfc_df$lr_logFC),
                    c(0.005, 0.995), na.rm = TRUE)
xlim_lfc <- c(logfc_q[1], logfc_q[2])
ylim_lfc <- xlim_lfc
n_clipped <- sum(logfc_df$sr_logFC < xlim_lfc[1] | logfc_df$sr_logFC > xlim_lfc[2] |
                 logfc_df$lr_logFC < ylim_lfc[1] | logfc_df$lr_logFC > ylim_lfc[2],
                 na.rm = TRUE)
clip_pct <- 100 * n_clipped / nrow(logfc_df)

label_x_b <- xlim_lfc[1] + 0.05 * diff(xlim_lfc)
label_y_b <- ylim_lfc[2] - 0.05 * diff(ylim_lfc)

pB <- ggplot(logfc_df, aes(x = sr_logFC, y = lr_logFC)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(alpha = 0.10, size = 0.4, color = "grey20") +
  geom_abline(slope = 1, intercept = 0, color = "red",
              linetype = "dashed", linewidth = 0.4) +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue",
              linewidth = 0.6) +
  geom_text(data = logfc_stats, aes(x = label_x_b, y = label_y_b, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.4,
            fontface = "bold", color = "black") +
  facet_wrap(~ ct, nrow = 1) +
  coord_cartesian(xlim = xlim_lfc, ylim = ylim_lfc) +
  labs(title = "Per-cell-type logFC concordance: SR (gene mashr) vs LR (gene-aggregated isoform mashr)",
       subtitle = sprintf("LR gene logFC = median of isoform mashr posterior logFC; axes clipped to 0.5-99.5%% (%d points / %.2f%% off-panel); red dashed = y=x",
                          n_clipped, clip_pct),
       x = "SR gene mashr logFC (Smg1i)",
       y = "LR gene mashr logFC (median across isoforms)") +
  paper_theme

# ============================================================================
# Compose & save
# ============================================================================

combined <- (pA / pB) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.margin = margin(10, 10, 10, 10))) &
  theme(plot.tag = element_text(size = 18, face = "bold"))

# Validate composite layout (skip-on-error so build doesn't fail silently)
val_path <- "~/.claude/utils/validate_composite_layout.R"
if (file.exists(path.expand(val_path))) {
  source(path.expand(val_path))
  tryCatch(
    validate_composite_layout(combined, fig_width = 14, fig_height = 9),
    error = function(e) message("Validator note: ", conditionMessage(e))
  )
}

out_pdf <- "paper_figures/FigX-SRLRConcordance_generated.pdf"
out_png <- "paper_figures/FigX-SRLRConcordance_generated.png"
ggsave(out_png, combined, width = 14, height = 9, dpi = 300, bg = "white")
ggsave(out_pdf, combined, width = 14, height = 9, bg = "white")

cat(sprintf("Saved: %s and %s\n", out_png, out_pdf))

# Print summary stats for the report
cat("\n=== Summary statistics ===\n")
cat("\nPanel A (expression):\n"); print(as.data.frame(expr_stats))
cat("\nPanel B (logFC):\n"); print(as.data.frame(logfc_stats))
