#!/usr/bin/env Rscript
# Paper Figure: Attribution Dissociation in Occult-PTC Isoforms
#
# Panels:
#   A) TD2 model SHAP — top features for occult-PTC vs non-occult NMD
#   B) Ref-CDS model SHAP — top features for occult-PTC vs non-occult NMD
#   C) Three-way ORF caller agreement (TD2 vs ref-ATG vs ORFik) by NMD subset
#   D) PTC-bearing rate among ORFik-ref-ATG agreements
#
# Population: holdout test isoforms (NMD only) split into:
#   - Occult-PTC: NMD pairs where ref-ATG analysis reclassified the
#     comparator from TD2-PTC-negative to effectively_ptc
#   - Non-occult NMD: all other NMD test isoforms
#
# Output: paper_figures/FigX-AttributionDissociation_generated.{pdf,png}

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
source("visualization_functions.R")
source("analysis_functions.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(Isopair)
})

paper_theme <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        plot.margin = margin(5, 8, 5, 5))

pal_subset <- c("Occult-PTC" = "#d73027", "Non-occult NMD" = "#4575b4")
pal_orfik <- c(
  "ORFik = ref-ATG only"   = "#1f78b4",
  "All three agree"        = "#666666",
  "ORFik picks third ATG"  = "#cccccc",
  "ORFik = TD2 only"       = "#e08214"
)

cat("=== Building Attribution Dissociation figure ===\n")

# ---- Load data ----
cache_dir <- "data_mashr/analysis_cache"
mc <- readRDS(file.path(cache_dir, "model_comparison.rds"))
ra <- readRDS(file.path(cache_dir, "ref_atg_analysis.rds"))
orfik <- readRDS(file.path(cache_dir, "orfik_scan.rds"))
structures <- readRDS("data_mashr/structures.rds")
cds <- readRDS("data_mashr/cds.rds")

# Identify occult-PTC isoforms (NMD pairs reclassified TD2-PTC- → effectively PTC+)
occult_ids <- ra$c2$comparator_isoform_id[
  ra$c2$category == "effectively_ptc" & !ra$c2$original_ptc]
cat(sprintf("Occult-PTC NMD pairs: %d\n", length(occult_ids)))

# ---- Panel A: TD2 model SHAP top features ----
td2_shap <- mc$td2_models$step4$shap_test
td2_test_df <- mc$test
td2_occult_idx <- which(td2_test_df$is_nmd == 1 &
                          td2_test_df$isoform_id %in% occult_ids)
td2_nonocc_idx <- setdiff(which(td2_test_df$is_nmd == 1), td2_occult_idx)

mean_abs <- function(mat, idx) colMeans(abs(mat[idx, , drop = FALSE]))
td2_occult_mean <- mean_abs(td2_shap, td2_occult_idx)
td2_nonocc_mean <- mean_abs(td2_shap, td2_nonocc_idx)

td2_shap_df <- data.frame(
  feature  = names(td2_occult_mean),
  occult   = as.numeric(td2_occult_mean),
  nonocc   = as.numeric(td2_nonocc_mean),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(pmax(occult, nonocc))) %>%
  head(8)

td2_shap_long <- td2_shap_df %>%
  pivot_longer(cols = c(occult, nonocc),
               names_to = "subset", values_to = "mean_abs_shap") %>%
  mutate(subset = factor(subset,
                         levels = c("occult", "nonocc"),
                         labels = c("Occult-PTC", "Non-occult NMD")),
         feature = factor(feature, levels = rev(td2_shap_df$feature)))

pA <- ggplot(td2_shap_long,
             aes(y = feature, x = mean_abs_shap, fill = subset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = pal_subset, name = NULL) +
  labs(title = sprintf("TD2 model SHAP (n=%s occult, n=%s non-occult)",
                       format(length(td2_occult_idx), big.mark = ","),
                       format(length(td2_nonocc_idx), big.mark = ",")),
       subtitle = "Length-priority CDS — attribution shifts to spurious 5'UTR/uORF features for occult-PTC",
       x = "Mean |SHAP|", y = NULL) +
  paper_theme +
  theme(legend.position = "top")

# ---- Panel B: Ref-CDS model SHAP top features ----
ref_shap <- mc$ref_models$step4$shap_test
ref_filter_idx <- which(!is.na(mc$test$ref_utr5_excluded) &
                          mc$test$ref_utr5_excluded == FALSE)
ref_test_df <- mc$test[ref_filter_idx, ]
ref_occult_idx <- which(ref_test_df$is_nmd == 1 &
                          ref_test_df$isoform_id %in% occult_ids)
ref_nonocc_idx <- setdiff(which(ref_test_df$is_nmd == 1), ref_occult_idx)

ref_occult_mean <- mean_abs(ref_shap, ref_occult_idx)
ref_nonocc_mean <- mean_abs(ref_shap, ref_nonocc_idx)

ref_shap_df <- data.frame(
  feature  = names(ref_occult_mean),
  occult   = as.numeric(ref_occult_mean),
  nonocc   = as.numeric(ref_nonocc_mean),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(pmax(occult, nonocc))) %>%
  head(8)

ref_shap_long <- ref_shap_df %>%
  pivot_longer(cols = c(occult, nonocc),
               names_to = "subset", values_to = "mean_abs_shap") %>%
  mutate(subset = factor(subset,
                         levels = c("occult", "nonocc"),
                         labels = c("Occult-PTC", "Non-occult NMD")),
         feature = factor(feature, levels = rev(ref_shap_df$feature)))

pB <- ggplot(ref_shap_long,
             aes(y = feature, x = mean_abs_shap, fill = subset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = pal_subset, name = NULL) +
  labs(title = sprintf("Ref-CDS model SHAP (n=%s occult, n=%s non-occult)",
                       format(length(ref_occult_idx), big.mark = ","),
                       format(length(ref_nonocc_idx), big.mark = ",")),
       subtitle = "Splicing-aware ATG — attribution concentrated on actual PTC for occult-PTC",
       x = "Mean |SHAP|", y = NULL) +
  paper_theme +
  theme(legend.position = "top")

# ---- Panel C: ORFik 3-way agreement stacked bar ----
struct_lookup <- setNames(
  lapply(seq_len(nrow(structures)), function(i) {
    list(starts = structures$exon_starts[[i]],
         ends   = structures$exon_ends[[i]],
         strand = structures$strand[i])
  }), structures$isoform_id)

orfik_features <- orfik$orf_features
orfik_top <- orfik_features %>%
  group_by(isoform_id) %>%
  arrange(desc(kozak_score), orf_start) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(isoform_id,
         orfik_atg_tx     = orf_start,
         orfik_kozak      = kozak_score,
         orfik_has_dejc   = has_downstream_ejc)

td2_atg_genomic_for <- function(iso_id) {
  cd <- cds[match(iso_id, cds$isoform_id), ]
  if (length(cd$strand) == 0 || is.na(cd$strand)) return(NA_real_)
  ifelse(cd$strand == "+", cd$cds_start, cd$cds_stop)
}
g2t_safe <- function(genomic, iso_id) {
  s <- struct_lookup[[iso_id]]
  if (is.null(s) || is.na(genomic)) return(NA_real_)
  Isopair::genomicToTranscript(genomic, s$starts, s$ends, s$strand)
}

ra_c2 <- ra$c2
ra_c2$td2_atg_tx <- vapply(seq_len(nrow(ra_c2)), function(i) {
  cid <- ra_c2$comparator_isoform_id[i]
  g2t_safe(td2_atg_genomic_for(cid), cid)
}, numeric(1))
ra_c2$ref_atg_tx <- vapply(seq_len(nrow(ra_c2)), function(i) {
  cid <- ra_c2$comparator_isoform_id[i]
  g2t_safe(ra_c2$ref_atg_genomic[i], cid)
}, numeric(1))
ra_c2 <- merge(ra_c2, orfik_top,
               by.x = "comparator_isoform_id", by.y = "isoform_id",
               all.x = TRUE)

ra_c2$subset <- dplyr::case_when(
  ra_c2$category == "effectively_ptc" & !ra_c2$original_ptc ~ "Occult-PTC",
  ra_c2$original_ptc ~ "Genuinely PTC+",
  ra_c2$category == "no_downstream_ejc" ~ "PTC- no_dEJC",
  ra_c2$category == "truncated_no_ejc"  ~ "PTC- truncated",
  TRUE ~ "Other"
)

resolvable <- !is.na(ra_c2$td2_atg_tx) & !is.na(ra_c2$ref_atg_tx) &
               !is.na(ra_c2$orfik_atg_tx)
TOL <- 2L
ra_c2$agree_orfik_td2 <- resolvable &
  abs(ra_c2$orfik_atg_tx - ra_c2$td2_atg_tx) <= TOL
ra_c2$agree_orfik_ref <- resolvable &
  abs(ra_c2$orfik_atg_tx - ra_c2$ref_atg_tx) <= TOL

ra_c2$bucket <- dplyr::case_when(
  !resolvable                                   ~ "Unresolvable",
  ra_c2$agree_orfik_td2 & ra_c2$agree_orfik_ref ~ "All three agree",
  ra_c2$agree_orfik_td2                         ~ "ORFik = TD2 only",
  ra_c2$agree_orfik_ref                         ~ "ORFik = ref-ATG only",
  TRUE                                          ~ "ORFik picks third ATG"
)

agreement_long <- ra_c2 %>%
  filter(subset != "Other", bucket != "Unresolvable") %>%
  count(subset, bucket) %>%
  group_by(subset) %>%
  mutate(pct = 100 * n / sum(n),
         total_n = sum(n)) %>%
  ungroup()

agreement_long$bucket <- factor(agreement_long$bucket,
  levels = c("ORFik = TD2 only", "ORFik picks third ATG",
             "All three agree", "ORFik = ref-ATG only"))
agreement_long$subset <- factor(agreement_long$subset,
  levels = c("Occult-PTC", "Genuinely PTC+", "PTC- no_dEJC", "PTC- truncated"))

subset_labels <- agreement_long %>%
  distinct(subset, total_n) %>%
  mutate(label = sprintf("%s\n(n=%s)", subset, format(total_n, big.mark = ","))) %>%
  arrange(subset)

agreement_long$subset_label <- factor(
  agreement_long$subset,
  levels = subset_labels$subset,
  labels = subset_labels$label
)

pC <- ggplot(agreement_long,
             aes(x = subset_label, y = pct, fill = bucket)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
  scale_fill_manual(values = pal_orfik, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Three-way ORF caller agreement",
       subtitle = "ORFik (Kozak-aware, no splicing context) corroborates ref-ATG in occult-PTC",
       x = NULL, y = "Percent of pairs") +
  paper_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0))

# ---- Panel D: PTC-bearing rate among ORFik-ref-ATG agreements ----
ptc_bearing_summary <- ra_c2 %>%
  filter(bucket == "ORFik = ref-ATG only", subset != "Other") %>%
  group_by(subset) %>%
  summarise(
    n_total = n(),
    n_dejc  = sum(orfik_has_dejc, na.rm = TRUE),
    pct_dejc = 100 * n_dejc / pmax(n_total, 1),
    .groups = "drop"
  )

ptc_bearing_summary$subset <- factor(ptc_bearing_summary$subset,
  levels = c("Occult-PTC", "Genuinely PTC+", "PTC- no_dEJC", "PTC- truncated"))

pD <- ggplot(ptc_bearing_summary,
             aes(x = subset, y = pct_dejc, fill = subset)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%\n(n=%d)", pct_dejc, n_total)),
            vjust = -0.3, size = 3) +
  scale_y_continuous(limits = c(0, 115),
                     breaks = seq(0, 100, 25),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c(
    "Occult-PTC"     = "#d73027",
    "Genuinely PTC+" = "#fc8d59",
    "PTC- no_dEJC"   = "#91bfdb",
    "PTC- truncated" = "#4575b4"
  ), guide = "none") +
  labs(title = "ORFik-ref-ATG agreements bear PTCs",
       subtitle = "Of pairs where ORFik agrees with ref-ATG, % of those ORFs with downstream EJC",
       x = NULL, y = "% with downstream EJC") +
  paper_theme

# ---- Compose multipanel ----
combined <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# ---- Validate before saving ----
cat("\nRunning composite layout validator...\n")
source("~/.claude/utils/validate_composite_layout.R")
fig_w <- 14; fig_h <- 10
val_result <- tryCatch(
  validate_composite_layout(combined, fig_width = fig_w, fig_height = fig_h),
  error = function(e) { cat("Validator error:", e$message, "\n"); NULL }
)

# ---- Save ----
out_dir <- "paper_figures"
ggsave(file.path(out_dir, "FigX-AttributionDissociation_generated.pdf"),
       combined, width = fig_w, height = fig_h, device = cairo_pdf)
ggsave(file.path(out_dir, "FigX-AttributionDissociation_generated.png"),
       combined, width = fig_w, height = fig_h, dpi = 300, bg = "white")

cat("\nSaved to:", out_dir, "\n")
cat(sprintf("  FigX-AttributionDissociation_generated.{pdf,png} (%g x %g in)\n",
            fig_w, fig_h))
