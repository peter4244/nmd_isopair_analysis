#!/usr/bin/env Rscript
# Paper Figure: Residual NMD Investigation
#
# Focus: ~15% of NMD coding pairs that are not explained by PTCs after
# reference-ATG reclassification. After ref-ATG analysis the residual splits
# into three subgroups:
#   - no_downstream_ejc : ref ATG retained, ORF terminates without downstream EJCs
#   - truncated_no_ejc  : ref ATG retained but ORF is truncated, no downstream EJCs
#   - ref_atg_lost      : ref ATG missing from comparator (different category)
#
# Panels:
#   A) Subgroup composition (3 NMD subgroups + Control)
#   B) 5'UTR length per subgroup (vs Control, Wilcoxon)
#   C) Longest 5'UTR ORF length per subgroup (vs Control, Wilcoxon)
#   D) NMD response magnitude (mashr posterior logFC averaged across AT/DD/FB/MV)
#
# Output: paper_figures/FigX-ResidualNMD_generated.{pdf,png}

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
source("visualization_functions.R")
source("analysis_functions.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(data.table)
})

paper_theme <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        plot.margin = margin(5, 8, 5, 5))

fmt <- function(x) format(x, big.mark = ",")
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 2.2e-16) return("< 2.2e-16")
  sprintf("%.1e", p)
}
sig_label <- function(p) {
  ifelse(is.na(p), "NA",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", "ns"))))
}

cat("=== Building Residual NMD Investigation figure ===\n")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
profiles_c2 <- readRDS("data_mashr/profiles_c2_allsamples.rds")
profiles_c4 <- readRDS("data_mashr/profiles_c4_allsamples.rds")
cds <- readRDS("data_mashr/cds.rds")

cache_dir <- "data_mashr/analysis_cache"
ref_atg <- readRDS(file.path(cache_dir, "ref_atg_analysis.rds"))
utr5_feat_all <- readRDS(file.path(cache_dir, "utr5_features_all.rds"))
div_c2 <- readRDS(file.path(cache_dir, "div_c2_allsamples.rds"))
div_c4 <- readRDS(file.path(cache_dir, "div_c4_allsamples.rds"))

mashr_csv <- "/Users/petecastaldi/claude_projects/nmd/isocall_dge/mashr/mashr_isoform_posterior_means_2026.3.10.csv"
mashr_pm <- fread(mashr_csv)

# ---------------------------------------------------------------------------
# Gene-match C2 and C4 (same convention as PTC accounting figure)
# ---------------------------------------------------------------------------
shared <- merge(
  unique(profiles_c2[, c("gene_id", "reference_isoform_id")]),
  unique(profiles_c4[, c("gene_id", "reference_isoform_id")]),
  by = c("gene_id", "reference_isoform_id"))
gm_c2 <- merge(profiles_c2, shared, by = c("gene_id", "reference_isoform_id"))
gm_c4 <- merge(profiles_c4, shared, by = c("gene_id", "reference_isoform_id"))

coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
coding_c2 <- gm_c2[gm_c2$reference_isoform_id %in% coding_ids &
                     gm_c2$comparator_isoform_id %in% coding_ids, ]
coding_c4 <- gm_c4[gm_c4$reference_isoform_id %in% coding_ids &
                     gm_c4$comparator_isoform_id %in% coding_ids, ]

# ---------------------------------------------------------------------------
# Subgroup definitions
# ---------------------------------------------------------------------------
c2_ref <- ref_atg$c2

# PTC-negative residual subgroups (the focus of this figure)
nde_rows   <- c2_ref[c2_ref$category == "no_downstream_ejc", ]
trunc_rows <- c2_ref[c2_ref$category == "truncated_no_ejc", ]
lost_rows  <- c2_ref[c2_ref$category == "ref_atg_lost", ]

no_ejc_ids   <- nde_rows$comparator_isoform_id
trunc_ids    <- trunc_rows$comparator_isoform_id
ref_lost_ids <- lost_rows$comparator_isoform_id
ctrl_ids     <- coding_c4$comparator_isoform_id

cat(sprintf("Subgroup counts (NMD): no_ejc=%d, truncated=%d, ref_atg_lost=%d | Control=%d\n",
            length(no_ejc_ids), length(trunc_ids), length(ref_lost_ids), length(ctrl_ids)))

# ---------------------------------------------------------------------------
# Shared aesthetics
# ---------------------------------------------------------------------------
group_levels <- c("No downstream\nEJC",
                  "Truncated ORF,\nno EJC",
                  "Ref ATG\nmissing",
                  "Control")
group_pal <- c("No downstream\nEJC"   = "#d95f02",
               "Truncated ORF,\nno EJC" = "#984ea3",
               "Ref ATG\nmissing"     = "#e6ab02",
               "Control"              = "#67a9cf")

# ---------------------------------------------------------------------------
# Panel A: Subgroup composition
# ---------------------------------------------------------------------------
total_nmd <- nrow(c2_ref)
bar_a <- data.frame(
  group = factor(group_levels,
                 levels = group_levels),
  n = c(length(no_ejc_ids), length(trunc_ids), length(ref_lost_ids), length(ctrl_ids)),
  is_nmd = c(TRUE, TRUE, TRUE, FALSE))

# Percent: NMD subgroups expressed as % of all NMD pairs; Control as % of C4
bar_a$pct <- ifelse(bar_a$is_nmd,
                    round(100 * bar_a$n / total_nmd, 1),
                    round(100 * bar_a$n / nrow(coding_c4), 1))
bar_a$denom_label <- ifelse(bar_a$is_nmd, "% of NMD pairs", "% of Control pairs")

pA <- ggplot(bar_a, aes(x = group, y = n, fill = group)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("%s\n(%s%%)", fmt(n), pct)),
            vjust = -0.3, size = 3.4, lineheight = 0.85) +
  scale_fill_manual(values = group_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "Residual PTC-negative NMD subgroups",
       subtitle = sprintf("NMD pairs: %s coding (C2) | Control: %s coding (C4)",
                          fmt(total_nmd), fmt(nrow(coding_c4))),
       x = NULL, y = "Number of pairs") +
  paper_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 9))

# ---------------------------------------------------------------------------
# Panel B: 5'UTR length comparison
# ---------------------------------------------------------------------------
# For no_downstream_ejc and truncated_no_ejc: use ref-ATG-derived 5'UTR length
#   5utr_len = (comp_stop_tx_pos - comp_orf_length_from_ref_atg) - 1
# For ref_atg_lost and Control: fall back to SQANTI 5'UTR length from div_*
utr5_lk_c2 <- setNames(div_c2$utr5_bp_comp, div_c2$comparator_isoform_id)
utr5_lk_c4 <- setNames(div_c4$utr5_bp_comp, div_c4$comparator_isoform_id)

nde_utr5 <- (nde_rows$comp_stop_tx_pos - nde_rows$comp_orf_length_from_ref_atg) - 1L
names(nde_utr5) <- nde_rows$comparator_isoform_id

trunc_utr5 <- (trunc_rows$comp_stop_tx_pos - trunc_rows$comp_orf_length_from_ref_atg) - 1L
names(trunc_utr5) <- trunc_rows$comparator_isoform_id

utr5_df <- rbind(
  data.frame(group = "No downstream\nEJC",      utr5 = unname(nde_utr5),
             source = "ref-ATG"),
  data.frame(group = "Truncated ORF,\nno EJC",   utr5 = unname(trunc_utr5),
             source = "ref-ATG"),
  data.frame(group = "Ref ATG\nmissing",        utr5 = unname(utr5_lk_c2[ref_lost_ids]),
             source = "SQANTI"),
  data.frame(group = "Control",                 utr5 = unname(utr5_lk_c4[ctrl_ids]),
             source = "SQANTI"))
utr5_df <- utr5_df[!is.na(utr5_df$utr5) & utr5_df$utr5 > 0, ]
utr5_df$group <- factor(utr5_df$group, levels = group_levels)

# Wilcoxon vs Control
utr5_split <- split(utr5_df$utr5, utr5_df$group)
utr5_pvals <- sapply(group_levels[1:3], function(g) {
  wilcox.test(utr5_split[[g]], utr5_split[["Control"]])$p.value
})
utr5_sig <- sig_label(utr5_pvals)

utr5_y_max <- tapply(log10(utr5_df$utr5), utr5_df$group, function(x)
  quantile(x, 0.75, na.rm = TRUE) + 1.5 * IQR(x, na.rm = TRUE))
sig_utr5 <- data.frame(
  group = factor(group_levels[1:3], levels = group_levels),
  y = 10^(pmax(utr5_y_max[group_levels[1:3]], utr5_y_max["Control"]) + 0.18),
  label = utr5_sig)

med_utr5 <- tapply(utr5_df$utr5, utr5_df$group, median, na.rm = TRUE)

pB <- ggplot(utr5_df, aes(x = group, y = utr5, fill = group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7, width = 0.6) +
  geom_text(data = sig_utr5, aes(x = group, y = y, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  scale_fill_manual(values = group_pal) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "5'UTR length",
       subtitle = sprintf("Median (bp): %s=%s, %s=%s, %s=%s, Ctrl=%s",
                          "no-EJC", fmt(round(med_utr5[group_levels[1]])),
                          "trunc",  fmt(round(med_utr5[group_levels[2]])),
                          "lost",   fmt(round(med_utr5[group_levels[3]])),
                                    fmt(round(med_utr5["Control"]))),
       x = NULL, y = "5'UTR length (bp, log scale)") +
  paper_theme +
  theme(legend.position = "none")

# ---------------------------------------------------------------------------
# Panel C: Longest 5'UTR ORF (uORF) per subgroup
# ---------------------------------------------------------------------------
# utr5_feat_all$isoform_features has per-isoform 5'UTR ORF features.
# longest_orf_nt is the longest 5'UTR ORF length in nt (0 if none).
utr5_feat <- utr5_feat_all$isoform_features %>% filter(!excluded)
uorf_lk <- setNames(utr5_feat$longest_orf_nt, utr5_feat$isoform_id)

uorf_df <- rbind(
  data.frame(group = "No downstream\nEJC",     longest_uorf = unname(uorf_lk[no_ejc_ids])),
  data.frame(group = "Truncated ORF,\nno EJC",  longest_uorf = unname(uorf_lk[trunc_ids])),
  data.frame(group = "Ref ATG\nmissing",       longest_uorf = unname(uorf_lk[ref_lost_ids])),
  data.frame(group = "Control",                longest_uorf = unname(uorf_lk[ctrl_ids])))
uorf_df <- uorf_df[!is.na(uorf_df$longest_uorf) & uorf_df$longest_uorf > 0, ]
uorf_df$group <- factor(uorf_df$group, levels = group_levels)

uorf_split <- split(uorf_df$longest_uorf, uorf_df$group)
uorf_pvals <- sapply(group_levels[1:3], function(g) {
  wilcox.test(uorf_split[[g]], uorf_split[["Control"]])$p.value
})
uorf_sig <- sig_label(uorf_pvals)

uorf_y_max <- tapply(log10(uorf_df$longest_uorf), uorf_df$group, function(x)
  quantile(x, 0.75, na.rm = TRUE) + 1.5 * IQR(x, na.rm = TRUE))
sig_uorf <- data.frame(
  group = factor(group_levels[1:3], levels = group_levels),
  y = 10^(pmax(uorf_y_max[group_levels[1:3]], uorf_y_max["Control"]) + 0.18),
  label = uorf_sig)

med_uorf <- tapply(uorf_df$longest_uorf, uorf_df$group, median, na.rm = TRUE)

pC <- ggplot(uorf_df, aes(x = group, y = longest_uorf, fill = group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7, width = 0.6) +
  geom_text(data = sig_uorf, aes(x = group, y = y, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  scale_fill_manual(values = group_pal) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "Longest 5'UTR ORF (uORF)",
       subtitle = sprintf("Median (nt): %s=%s, %s=%s, %s=%s, Ctrl=%s",
                          "no-EJC", fmt(round(med_uorf[group_levels[1]])),
                          "trunc",  fmt(round(med_uorf[group_levels[2]])),
                          "lost",   fmt(round(med_uorf[group_levels[3]])),
                                    fmt(round(med_uorf["Control"]))),
       x = NULL, y = "Longest uORF (nt, log scale)") +
  paper_theme +
  theme(legend.position = "none")

# ---------------------------------------------------------------------------
# Panel D: NMD response magnitude (mashr posterior logFC, averaged 4 CTs)
# ---------------------------------------------------------------------------
# Average the four cell-type Smg1i posterior means per isoform.
ct_cols <- c("Smg1i_in_AT", "Smg1i_in_DD", "Smg1i_in_FB", "Smg1i_in_MV")
stopifnot(all(ct_cols %in% colnames(mashr_pm)))
mashr_pm[, mean_lfc := rowMeans(.SD, na.rm = TRUE), .SDcols = ct_cols]
lfc_lk <- setNames(mashr_pm$mean_lfc, mashr_pm$txid)

lfc_df <- rbind(
  data.frame(group = "No downstream\nEJC",     lfc = unname(lfc_lk[no_ejc_ids])),
  data.frame(group = "Truncated ORF,\nno EJC",  lfc = unname(lfc_lk[trunc_ids])),
  data.frame(group = "Ref ATG\nmissing",       lfc = unname(lfc_lk[ref_lost_ids])),
  data.frame(group = "Control",                lfc = unname(lfc_lk[ctrl_ids])))
lfc_df <- lfc_df[!is.na(lfc_df$lfc), ]
lfc_df$group <- factor(lfc_df$group, levels = group_levels)

lfc_split <- split(lfc_df$lfc, lfc_df$group)
lfc_pvals <- sapply(group_levels[1:3], function(g) {
  wilcox.test(lfc_split[[g]], lfc_split[["Control"]])$p.value
})
lfc_sig <- sig_label(lfc_pvals)

# Position significance labels above whisker tops
lfc_y_max <- tapply(lfc_df$lfc, lfc_df$group, function(x)
  quantile(x, 0.75, na.rm = TRUE) + 1.5 * IQR(x, na.rm = TRUE))
sig_lfc <- data.frame(
  group = factor(group_levels[1:3], levels = group_levels),
  y = pmax(lfc_y_max[group_levels[1:3]], lfc_y_max["Control"]) + 0.15,
  label = lfc_sig)

med_lfc <- tapply(lfc_df$lfc, lfc_df$group, median, na.rm = TRUE)

pD <- ggplot(lfc_df, aes(x = group, y = lfc, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7, width = 0.6) +
  geom_text(data = sig_lfc, aes(x = group, y = y, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  scale_fill_manual(values = group_pal) +
  labs(title = "NMD response magnitude (Smg1i)",
       subtitle = sprintf("Median logFC: no-EJC=%.2f, trunc=%.2f, lost=%.2f, Ctrl=%.2f (mean across AT/DD/FB/MV)",
                          med_lfc[group_levels[1]],
                          med_lfc[group_levels[2]],
                          med_lfc[group_levels[3]],
                          med_lfc["Control"]),
       x = NULL, y = expression("Mean posterior log"[2]*"FC")) +
  paper_theme +
  theme(legend.position = "none")

# ---------------------------------------------------------------------------
# Compose
# ---------------------------------------------------------------------------
combined <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A",
    theme = theme(plot.margin = margin(10, 10, 10, 10))) &
  theme(plot.tag = element_text(size = 18, face = "bold"))

# ---------------------------------------------------------------------------
# Validate before saving
# ---------------------------------------------------------------------------
source("~/.claude/utils/validate_composite_layout.R")
val_result <- validate_composite_layout(combined, fig_width = 13, fig_height = 9)
cat("\nValidator result:\n")
print(val_result)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
ggsave("paper_figures/FigX-ResidualNMD_generated.png", combined,
       width = 13, height = 9, dpi = 300, bg = "white")
ggsave("paper_figures/FigX-ResidualNMD_generated.pdf", combined,
       width = 13, height = 9, device = cairo_pdf, bg = "white")

cat("\nSaved: paper_figures/FigX-ResidualNMD_generated.png and .pdf\n")

# ---------------------------------------------------------------------------
# Summary table for sanity
# ---------------------------------------------------------------------------
summary_tbl <- data.frame(
  group = group_levels,
  n_pairs = c(length(no_ejc_ids), length(trunc_ids), length(ref_lost_ids), length(ctrl_ids)),
  median_utr5_bp = round(med_utr5[group_levels]),
  median_uorf_nt = round(med_uorf[group_levels]),
  median_lfc = round(med_lfc[group_levels], 3),
  utr5_p_vs_ctrl = c(sapply(utr5_pvals, fmt_p), "—"),
  uorf_p_vs_ctrl = c(sapply(uorf_pvals, fmt_p), "—"),
  lfc_p_vs_ctrl  = c(sapply(lfc_pvals, fmt_p),  "—"))
cat("\nPer-subgroup summary:\n")
print(summary_tbl, row.names = FALSE)
