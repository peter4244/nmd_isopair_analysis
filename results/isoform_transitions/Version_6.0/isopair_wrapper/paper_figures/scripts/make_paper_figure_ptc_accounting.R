#!/usr/bin/env Rscript
# Paper Figure: PTC Accounting and Residual PTC-Negative Characterization
# Panels: A) NMD classification group counts
#         B) 3'UTR length (3 residual groups + Control)
#         C) 5'UTR length density (same 4 groups)
#         D) Longest uORF (same 4 groups)

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
library(ggplot2)
library(patchwork)
library(dplyr)

paper_theme <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        plot.margin = margin(5, 8, 5, 5))

fmt <- function(x) format(x, big.mark = ",")

# ---- Load data ----
profiles_c2 <- readRDS("data_mashr/profiles_c2_allsamples.rds")
profiles_c4 <- readRDS("data_mashr/profiles_c4_allsamples.rds")
cds <- readRDS("data_mashr/cds.rds")
ptc_data <- readRDS("data_mashr/ptc.rds")
structures <- readRDS("data_mashr/structures.rds")

cache_dir <- "data_mashr/analysis_cache"
div_c2 <- readRDS(file.path(cache_dir, "div_c2_allsamples.rds"))
div_c4 <- readRDS(file.path(cache_dir, "div_c4_allsamples.rds"))
ref_atg <- readRDS(file.path(cache_dir, "ref_atg_analysis.rds"))
utr5_feat_all <- readRDS(file.path(cache_dir, "utr5_features_all.rds"))

# ---- Gene-match C2 and C4 ----
shared <- merge(
  unique(profiles_c2[, c("gene_id", "reference_isoform_id")]),
  unique(profiles_c4[, c("gene_id", "reference_isoform_id")]),
  by = c("gene_id", "reference_isoform_id"))
gm_c2 <- merge(profiles_c2, shared, by = c("gene_id", "reference_isoform_id"))
gm_c4 <- merge(profiles_c4, shared, by = c("gene_id", "reference_isoform_id"))

# Restrict to coding pairs
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
coding_c2 <- gm_c2[gm_c2$reference_isoform_id %in% coding_ids &
                      gm_c2$comparator_isoform_id %in% coding_ids, ]
coding_c4 <- gm_c4[gm_c4$reference_isoform_id %in% coding_ids &
                      gm_c4$comparator_isoform_id %in% coding_ids, ]

# ---- Classification groups ----
c2_ref <- ref_atg$c2
ptc_lookup <- setNames(ptc_data$has_ptc, ptc_data$isoform_id)

# Original PTC+ (from coding pairs)
coding_c2$comp_has_ptc <- ptc_lookup[coding_c2$comparator_isoform_id]
n_ptc_pos <- sum(coding_c2$comp_has_ptc == TRUE, na.rm = TRUE)

# PTC-negative categories from ref-ATG analysis
ptc_neg_ref <- c2_ref[!c2_ref$original_ptc, ]
n_eff_ptc     <- sum(ptc_neg_ref$category == "effectively_ptc")
n_no_ejc      <- sum(ptc_neg_ref$category == "no_downstream_ejc")
n_trunc_noejc <- sum(ptc_neg_ref$category == "truncated_no_ejc")
n_ref_lost    <- sum(ptc_neg_ref$category == "ref_atg_lost")

# Residual group isoform IDs
no_ejc_ids     <- ptc_neg_ref$comparator_isoform_id[ptc_neg_ref$category == "no_downstream_ejc"]
trunc_ids      <- ptc_neg_ref$comparator_isoform_id[ptc_neg_ref$category == "truncated_no_ejc"]
ref_lost_ids   <- ptc_neg_ref$comparator_isoform_id[ptc_neg_ref$category == "ref_atg_lost"]
ctrl_comp_ids  <- coding_c4$comparator_isoform_id

# ---- Panel A: Classification barplot ----
panel_a_pal <- c("NMD+/PTC+" = "#d6604d",
                  "Reclassified\nNMD+/PTC+" = "#ef8a62",
                  "Ref ATG intact,\nno main ORF PTC" = "#d95f02",
                  "Ref ATG intact,\ntruncated ORF" = "#984ea3",
                  "Ref ATG\nmissing" = "#e6ab02")

bar_a <- data.frame(
  group = factor(c("NMD+/PTC+", "Reclassified\nNMD+/PTC+",
                    "Ref ATG intact,\nno main ORF PTC",
                    "Ref ATG intact,\ntruncated ORF",
                    "Ref ATG\nmissing"),
                 levels = c("NMD+/PTC+", "Reclassified\nNMD+/PTC+",
                            "Ref ATG intact,\nno main ORF PTC",
                            "Ref ATG intact,\ntruncated ORF",
                            "Ref ATG\nmissing")),
  n = c(n_ptc_pos, n_eff_ptc, n_no_ejc, n_trunc_noejc, n_ref_lost))

bar_a$pct <- round(100 * bar_a$n / sum(bar_a$n), 1)

pA <- ggplot(bar_a, aes(x = group, y = n, fill = group)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("%s\n(%s%%)", fmt(n), pct)),
            vjust = -0.3, size = 3.2, lineheight = 0.85) +
  scale_fill_manual(values = panel_a_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "NMD Pair Classification After Reference-ATG Reclassification",
       subtitle = sprintf("%s total coding NMD pairs", fmt(sum(bar_a$n))),
       x = NULL, y = "Number of pairs") +
  paper_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 9))

# ---- Panels B-D: 4-group comparisons ----
# Consistent group labels and colors for residual + Control
group_levels <- c("No main ORF\nPTC", "Truncated\nORF", "Ref ATG\nmissing", "Control")
group_pal <- c("No main ORF\nPTC" = "#d95f02",
               "Truncated\nORF" = "#984ea3",
               "Ref ATG\nmissing" = "#e6ab02",
               "Control" = "#67a9cf")

# ---- Panel B: 3'UTR length boxplot ----
# Use SQANTI-based 3'UTR from divergence (appropriate for PTC-negative groups)
utr3_lk_c2 <- setNames(div_c2$utr3_bp_comp, div_c2$comparator_isoform_id)
utr3_lk_c4 <- setNames(div_c4$utr3_bp_comp, div_c4$comparator_isoform_id)

utr3_df <- rbind(
  data.frame(group = "No main ORF\nPTC",  utr3 = utr3_lk_c2[no_ejc_ids]),
  data.frame(group = "Truncated\nORF",     utr3 = utr3_lk_c2[trunc_ids]),
  data.frame(group = "Ref ATG\nmissing",   utr3 = utr3_lk_c2[ref_lost_ids]),
  data.frame(group = "Control",            utr3 = utr3_lk_c4[ctrl_comp_ids]))
utr3_df <- utr3_df[!is.na(utr3_df$utr3) & utr3_df$utr3 > 0, ]
utr3_df$group <- factor(utr3_df$group, levels = group_levels)

# Wilcoxon tests vs Control
utr3_split <- split(utr3_df$utr3, utr3_df$group)
utr3_tests <- sapply(group_levels[1:3], function(g) {
  wt <- wilcox.test(utr3_split[[g]], utr3_split[["Control"]])
  wt$p.value
})
utr3_sig <- ifelse(utr3_tests < 0.001, "***",
              ifelse(utr3_tests < 0.01, "**",
                ifelse(utr3_tests < 0.05, "*", "ns")))

# Annotation positions
utr3_y_max <- tapply(log10(utr3_df$utr3), utr3_df$group, function(x)
  quantile(x, 0.75, na.rm = TRUE) + 1.5 * IQR(x, na.rm = TRUE))
sig_utr3 <- data.frame(
  group = factor(group_levels[1:3], levels = group_levels),
  y = 10^(pmax(utr3_y_max[group_levels[1:3]], utr3_y_max["Control"]) + 0.15),
  label = utr3_sig)

pB <- ggplot(utr3_df, aes(x = group, y = utr3, fill = group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7, width = 0.6) +
  geom_text(data = sig_utr3, aes(x = group, y = y, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  scale_fill_manual(values = group_pal) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "3'UTR Length",
       subtitle = "Significance vs Control (Wilcoxon)",
       x = NULL, y = "3'UTR length (bp, log scale)") +
  paper_theme +
  theme(legend.position = "none")

# ---- Panel C: 5'UTR length density ----
utr5_lk_c2 <- setNames(div_c2$utr5_bp_comp, div_c2$comparator_isoform_id)
utr5_lk_c4 <- setNames(div_c4$utr5_bp_comp, div_c4$comparator_isoform_id)

utr5_df <- rbind(
  data.frame(group = "No main ORF\nPTC",  utr5 = utr5_lk_c2[no_ejc_ids]),
  data.frame(group = "Truncated\nORF",     utr5 = utr5_lk_c2[trunc_ids]),
  data.frame(group = "Ref ATG\nmissing",   utr5 = utr5_lk_c2[ref_lost_ids]),
  data.frame(group = "Control",            utr5 = utr5_lk_c4[ctrl_comp_ids]))
utr5_df <- utr5_df[!is.na(utr5_df$utr5) & utr5_df$utr5 > 0, ]
utr5_df$group <- factor(utr5_df$group, levels = group_levels)

# Medians and Wilcoxon tests vs Control for subtitle
utr5_med <- tapply(utr5_df$utr5, utr5_df$group, median, na.rm = TRUE)
utr5_split <- split(utr5_df$utr5, utr5_df$group)
utr5_pvals <- sapply(group_levels[1:3], function(g) {
  wilcox.test(utr5_split[[g]], utr5_split[["Control"]])$p.value
})
fmt_p <- function(p) ifelse(p < 2.2e-16, "< 2.2e-16", sprintf("%.1e", p))

# Single-line legend labels for density plot
legend_labels <- c("No main ORF\nPTC" = "No main ORF PTC",
                   "Truncated\nORF" = "Truncated ORF",
                   "Ref ATG\nmissing" = "Ref ATG missing",
                   "Control" = "Control")

pC <- ggplot(utr5_df, aes(x = utr5, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values = group_pal, labels = legend_labels) +
  labs(title = "5'UTR Length",
       subtitle = sprintf("vs Control (Wilcoxon): No PTC p %s | Trunc p = %s | Lost p %s",
                          fmt_p(utr5_pvals[1]), fmt_p(utr5_pvals[2]),
                          fmt_p(utr5_pvals[3])),
       x = "5'UTR length (bp, log scale)", y = "Density", fill = NULL) +
  paper_theme +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.4, "cm"))

# ---- Panel D: Longest uORF boxplot ----
utr5_feat <- utr5_feat_all$isoform_features %>% filter(!excluded)

uorf_df <- rbind(
  data.frame(group = "No main ORF\nPTC",
             longest_uorf = utr5_feat$longest_orf_nt[utr5_feat$isoform_id %in% no_ejc_ids]),
  data.frame(group = "Truncated\nORF",
             longest_uorf = utr5_feat$longest_orf_nt[utr5_feat$isoform_id %in% trunc_ids]),
  data.frame(group = "Ref ATG\nmissing",
             longest_uorf = utr5_feat$longest_orf_nt[utr5_feat$isoform_id %in% ref_lost_ids]),
  data.frame(group = "Control",
             longest_uorf = utr5_feat$longest_orf_nt[utr5_feat$isoform_id %in% ctrl_comp_ids]))
uorf_df <- uorf_df[!is.na(uorf_df$longest_uorf) & uorf_df$longest_uorf > 0, ]
uorf_df$group <- factor(uorf_df$group, levels = group_levels)

# Wilcoxon tests vs Control
uorf_split <- split(uorf_df$longest_uorf, uorf_df$group)
uorf_tests <- sapply(group_levels[1:3], function(g) {
  wt <- wilcox.test(uorf_split[[g]], uorf_split[["Control"]])
  wt$p.value
})
uorf_sig <- ifelse(uorf_tests < 0.001, "***",
              ifelse(uorf_tests < 0.01, "**",
                ifelse(uorf_tests < 0.05, "*", "ns")))

uorf_y_max <- tapply(log10(uorf_df$longest_uorf), uorf_df$group, function(x)
  quantile(x, 0.75, na.rm = TRUE) + 1.5 * IQR(x, na.rm = TRUE))
sig_uorf <- data.frame(
  group = factor(group_levels[1:3], levels = group_levels),
  y = 10^(pmax(uorf_y_max[group_levels[1:3]], uorf_y_max["Control"]) + 0.15),
  label = uorf_sig)

pD <- ggplot(uorf_df, aes(x = group, y = longest_uorf, fill = group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7, width = 0.6) +
  geom_text(data = sig_uorf, aes(x = group, y = y, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  scale_fill_manual(values = group_pal) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "Longest 5'UTR ORF",
       subtitle = "Significance vs Control (Wilcoxon)",
       x = NULL, y = "Longest uORF (nt, log scale)") +
  paper_theme +
  theme(legend.position = "none")

# ---- Compose ----
combined <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A",
    theme = theme(plot.margin = margin(10, 10, 10, 10))) &
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave("paper_figures/FigX-PTCAccounting_generated.png", combined,
       width = 16, height = 10, dpi = 300, bg = "white")
ggsave("paper_figures/FigX-PTCAccounting_generated.pdf", combined,
       width = 16, height = 10, bg = "white")

cat("Saved: paper_figures/FigX-PTCAccounting_generated.png and .pdf\n")
