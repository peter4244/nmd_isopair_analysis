#!/usr/bin/env Rscript
# ==============================================================================
# Phase 0: Disambiguate the 5'UTR signal in residual PTC-negative NMD pairs
#
# Before investigating uORF-mediated NMD, we must determine whether the long
# 5'UTRs in ref_atg_lost (n=240) and no_downstream_ejc (n=187) groups are:
#   (a) Genuine 5'UTR extensions (biological), or
#   (b) Artifacts of TD2 CDS misprediction (same artifact as reclassified PTCs)
#
# Approach:
#   0a. Compare SQANTI 5'UTR vs ref-ATG 5'UTR for no_downstream_ejc
#   0b. Characterize TSS position for ref_atg_lost (Alt_TSS vs exon loss)
#   0c. Compute ref-ATG 3'UTR for no_downstream_ejc
#   0d. Include PTC+ as artifact control
#
# Inputs:
#   - data_mashr/analysis_cache/ref_atg_analysis.rds
#   - data_mashr/cds.rds, structures.rds
#   - data_mashr/profiles_c2_allsamples.rds
#   - mashr DE files (for lfsr)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
fmt <- function(x) format(x, big.mark = ",")

cat("=== Phase 0: Disambiguate 5'UTR Signal ===\n\n")

# Load data
ref_atg <- readRDS("data_mashr/analysis_cache/ref_atg_analysis.rds")
cds <- readRDS("data_mashr/cds.rds")
structures <- readRDS("data_mashr/structures.rds")
profiles_c2 <- readRDS("data_mashr/profiles_c2_allsamples.rds")
profiles_c4 <- readRDS("data_mashr/profiles_c4_allsamples.rds")

# Transcript lengths
tx_lengths <- setNames(
  sapply(seq_len(nrow(structures)), function(i)
    sum(structures$exon_ends[[i]] - structures$exon_starts[[i]] + 1L)),
  structures$isoform_id)

# ==============================================================================
# Helper: compute SQANTI 5'UTR from TD2 CDS
# ==============================================================================
compute_sqanti_5utr <- function(iso_ids) {
  sapply(iso_ids, function(iso_id) {
    cds_row <- cds[cds$isoform_id == iso_id, ]
    if (nrow(cds_row) == 0 || cds_row$coding_status != "coding") return(NA_real_)
    str_row <- structures[structures$isoform_id == iso_id, ]
    if (nrow(str_row) == 0) return(NA_real_)
    strand <- cds_row$strand
    cds5 <- if (strand == "+") cds_row$cds_start else cds_row$cds_stop
    exon_s <- str_row$exon_starts[[1]]; exon_e <- str_row$exon_ends[[1]]
    utr5 <- 0
    for (j in seq_along(exon_s)) {
      if (strand == "+") {
        if (exon_s[j] < cds5) utr5 <- utr5 + min(exon_e[j], cds5 - 1) - exon_s[j] + 1
      } else {
        if (exon_e[j] > cds5) utr5 <- utr5 + exon_e[j] - max(exon_s[j], cds5 + 1) + 1
      }
    }
    utr5
  })
}

# ==============================================================================
# Define populations
# ==============================================================================
c2_res <- ref_atg$c2
ptc_neg <- c2_res[!c2_res$original_ptc, ]

# Residual groups
ref_atg_lost <- ptc_neg[ptc_neg$category == "ref_atg_lost", ]
no_dejc <- ptc_neg[ptc_neg$category == "no_downstream_ejc", ]
trunc_no_ejc <- ptc_neg[ptc_neg$category == "truncated_no_ejc", ]
eff_ptc <- ptc_neg[ptc_neg$category == "effectively_ptc", ]

# PTC+ pairs (artifact control)
ptc_pos <- c2_res[c2_res$original_ptc, ]

# Control
c4_res <- ref_atg$c4
c4_avail <- c4_res[c4_res$ref_atg_exonic_in_comp == TRUE &
                     !is.na(c4_res$comp_stop_tx_pos), ]

cat("Populations:\n")
cat("  ref_atg_lost:", nrow(ref_atg_lost), "\n")
cat("  no_downstream_ejc:", nrow(no_dejc), "\n")
cat("  truncated_no_ejc:", nrow(trunc_no_ejc), "\n")
cat("  effectively_ptc (reclassified):", nrow(eff_ptc), "\n")
cat("  original PTC+:", nrow(ptc_pos), "\n")
cat("  Control (C4, ref ATG avail):", nrow(c4_avail), "\n")

# ==============================================================================
# 0a. SQANTI 5'UTR vs ref-ATG 5'UTR for no_downstream_ejc
# ==============================================================================
cat("\n=== 0a. SQANTI vs Ref-ATG 5'UTR (no_downstream_ejc) ===\n")

# SQANTI 5'UTR (from TD2 CDS)
no_dejc$sqanti_5utr <- compute_sqanti_5utr(no_dejc$comparator_isoform_id)

# Ref-ATG 5'UTR: ATG tx position = comp_stop_tx_pos - comp_orf_length_from_ref_atg
# 5'UTR = ATG_tx_pos - 1
no_dejc$ref_atg_tx_pos <- no_dejc$comp_stop_tx_pos - no_dejc$comp_orf_length_from_ref_atg
no_dejc$ref_5utr <- no_dejc$ref_atg_tx_pos - 1L

cat("  SQANTI 5'UTR: median =", round(median(no_dejc$sqanti_5utr, na.rm = TRUE)),
    "bp, IQR = [", round(quantile(no_dejc$sqanti_5utr, 0.25, na.rm = TRUE)), ",",
    round(quantile(no_dejc$sqanti_5utr, 0.75, na.rm = TRUE)), "]\n")
cat("  Ref-ATG 5'UTR: median =", round(median(no_dejc$ref_5utr, na.rm = TRUE)),
    "bp, IQR = [", round(quantile(no_dejc$ref_5utr, 0.25, na.rm = TRUE)), ",",
    round(quantile(no_dejc$ref_5utr, 0.75, na.rm = TRUE)), "]\n")

# If ref-ATG 5'UTR is much shorter than SQANTI, it's a TD2 artifact
cat("  Ratio SQANTI/Ref-ATG: median =",
    round(median(no_dejc$sqanti_5utr / pmax(no_dejc$ref_5utr, 1), na.rm = TRUE), 2), "\n")
cat("  Cases where SQANTI > 2x Ref-ATG:",
    sum(no_dejc$sqanti_5utr > 2 * no_dejc$ref_5utr, na.rm = TRUE),
    "/", sum(!is.na(no_dejc$ref_5utr)), "\n")

# Control comparison
c4_avail$sqanti_5utr <- compute_sqanti_5utr(c4_avail$comparator_isoform_id)
c4_avail$ref_atg_tx_pos <- c4_avail$comp_stop_tx_pos - c4_avail$comp_orf_length_from_ref_atg
c4_avail$ref_5utr <- c4_avail$ref_atg_tx_pos - 1L

cat("\n  Control SQANTI 5'UTR: median =",
    round(median(c4_avail$sqanti_5utr, na.rm = TRUE)), "bp\n")
cat("  Control Ref-ATG 5'UTR: median =",
    round(median(c4_avail$ref_5utr, na.rm = TRUE)), "bp\n")

# ==============================================================================
# 0b. Characterize ref_atg_lost: what removed the ATG?
# ==============================================================================
cat("\n=== 0b. What removed the reference ATG? (ref_atg_lost) ===\n")

# Match to profiles to get splicing events
ref_lost_profiles <- profiles_c2 %>%
  filter(comparator_isoform_id %in% ref_atg_lost$comparator_isoform_id)

cat("  Pairs with profiles:", nrow(ref_lost_profiles), "\n")
cat("  TSS changed:", sum(ref_lost_profiles$tss_changed),
    sprintf("(%.1f%%)\n", 100 * mean(ref_lost_profiles$tss_changed)))
cat("  Has SE:", sum(ref_lost_profiles$n_se > 0),
    sprintf("(%.1f%%)\n", 100 * mean(ref_lost_profiles$n_se > 0)))
cat("  Has Missing_Internal:", sum(ref_lost_profiles$n_missing_internal > 0),
    sprintf("(%.1f%%)\n", 100 * mean(ref_lost_profiles$n_missing_internal > 0)))

# SQANTI 5'UTR for ref_atg_lost
ref_atg_lost$sqanti_5utr <- compute_sqanti_5utr(ref_atg_lost$comparator_isoform_id)
cat("\n  SQANTI 5'UTR: median =", round(median(ref_atg_lost$sqanti_5utr, na.rm = TRUE)), "bp\n")

# ==============================================================================
# 0c. Ref-ATG 3'UTR for no_downstream_ejc
# ==============================================================================
cat("\n=== 0c. Ref-ATG 3'UTR (no_downstream_ejc) ===\n")

# 3'UTR = transcript_length - (stop_tx_pos + 2)
no_dejc$tx_len <- tx_lengths[no_dejc$comparator_isoform_id]
no_dejc$ref_3utr <- no_dejc$tx_len - (no_dejc$comp_stop_tx_pos + 2L)

cat("  Ref-ATG 3'UTR: median =", round(median(no_dejc$ref_3utr, na.rm = TRUE)),
    "bp, IQR = [", round(quantile(no_dejc$ref_3utr, 0.25, na.rm = TRUE)), ",",
    round(quantile(no_dejc$ref_3utr, 0.75, na.rm = TRUE)), "]\n")

c4_avail$tx_len <- tx_lengths[c4_avail$comparator_isoform_id]
c4_avail$ref_3utr <- c4_avail$tx_len - (c4_avail$comp_stop_tx_pos + 2L)

cat("  Control ref-ATG 3'UTR: median =",
    round(median(c4_avail$ref_3utr, na.rm = TRUE)), "bp\n")

# Is ref-ATG 3'UTR significantly longer than control?
wt_3utr <- wilcox.test(no_dejc$ref_3utr, c4_avail$ref_3utr)
cat("  Wilcoxon p (no_dejc vs Control):", signif(wt_3utr$p.value, 3), "\n")

# ==============================================================================
# 0d. PTC+ artifact control: SQANTI 5'UTR
# ==============================================================================
cat("\n=== 0d. PTC+ artifact control ===\n")

ptc_pos$sqanti_5utr <- compute_sqanti_5utr(ptc_pos$comparator_isoform_id)
cat("  PTC+ SQANTI 5'UTR: median =",
    round(median(ptc_pos$sqanti_5utr, na.rm = TRUE)), "bp\n")

eff_ptc$sqanti_5utr <- compute_sqanti_5utr(eff_ptc$comparator_isoform_id)
cat("  Reclassified PTC SQANTI 5'UTR: median =",
    round(median(eff_ptc$sqanti_5utr, na.rm = TRUE)), "bp\n")

cat("  Control SQANTI 5'UTR: median =",
    round(median(c4_avail$sqanti_5utr, na.rm = TRUE)), "bp\n")

# ==============================================================================
# 0e. Load lfsr for residual classification confidence
# ==============================================================================
cat("\n=== 0e. Classification confidence (lfsr) ===\n")

mashr_model <- readRDS("/Users/petecastaldi/claude_projects/nmd/isocall_dge/mashr/mashr_isoform_model_2026.3.10.rds")
lfsr_mat <- ashr::get_lfsr(mashr_model)
# Use minimum lfsr across conditions (most confident condition)
main_cols <- c("Smg1i_in_DD", "Smg1i_in_AT", "Smg1i_in_DO_ALI",
               "Smg1i_in_FB", "Smg1i_in_MV")
min_lfsr <- apply(lfsr_mat[, main_cols], 1, min)

get_lfsr <- function(iso_ids) {
  idx <- match(iso_ids, rownames(lfsr_mat))
  min_lfsr[idx]
}

groups <- list(
  "ref_atg_lost" = ref_atg_lost$comparator_isoform_id,
  "no_downstream_ejc" = no_dejc$comparator_isoform_id,
  "truncated_no_ejc" = trunc_no_ejc$comparator_isoform_id,
  "reclassified_PTC" = eff_ptc$comparator_isoform_id,
  "original_PTC+" = ptc_pos$comparator_isoform_id
)

for (grp in names(groups)) {
  lf <- get_lfsr(groups[[grp]])
  lf <- lf[!is.na(lf)]
  cat(sprintf("  %s (n=%d): median lfsr = %.4f, %% near boundary (>0.01) = %.1f%%\n",
              grp, length(lf), median(lf),
              100 * mean(lf > 0.01)))
}

# ==============================================================================
# 0f. NMD response magnitude (logFC) by group
# ==============================================================================
cat("\n=== 0f. NMD response magnitude ===\n")

meta_logfc <- {
  post_mean <- ashr::get_pm(mashr_model)
  data.frame(
    txid = rownames(post_mean),
    mean_logFC = rowMeans(post_mean[, main_cols]),
    stringsAsFactors = FALSE)
}

for (grp in names(groups)) {
  lfc <- meta_logfc$mean_logFC[match(groups[[grp]], meta_logfc$txid)]
  lfc <- lfc[!is.na(lfc)]
  cat(sprintf("  %s: median logFC = %.3f, IQR = [%.3f, %.3f]\n",
              grp, median(lfc), quantile(lfc, 0.25), quantile(lfc, 0.75)))
}

# ==============================================================================
# Summary figure
# ==============================================================================
cat("\n=== Building summary figure ===\n")

theme_set(theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank()))

# Panel 1: SQANTI vs Ref-ATG 5'UTR for no_downstream_ejc
utr_compare <- rbind(
  data.frame(group = "no_dejc", source = "SQANTI (TD2)", utr5 = no_dejc$sqanti_5utr),
  data.frame(group = "no_dejc", source = "Ref-ATG", utr5 = no_dejc$ref_5utr),
  data.frame(group = "Control", source = "SQANTI (TD2)", utr5 = c4_avail$sqanti_5utr),
  data.frame(group = "Control", source = "Ref-ATG", utr5 = c4_avail$ref_5utr)
)

p1 <- ggplot(utr_compare, aes(x = interaction(group, source, sep = "\n"),
                                y = utr5, fill = source)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = c("SQANTI (TD2)" = "#ef8a62", "Ref-ATG" = "#67a9cf")) +
  labs(title = "5'UTR: SQANTI (TD2) vs Reference ATG",
       subtitle = "no_downstream_ejc vs Control",
       x = NULL, y = "5'UTR length (bp, log scale)") +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 10))

ggsave("figures/explore_phase0_5utr_comparison.png", p1, width = 10, height = 6, dpi = 300)
ggsave("figures/explore_phase0_5utr_comparison.pdf", p1, width = 10, height = 6)

# Panel 2: logFC by group
logfc_data <- do.call(rbind, lapply(names(groups), function(grp) {
  lfc <- meta_logfc$mean_logFC[match(groups[[grp]], meta_logfc$txid)]
  data.frame(group = grp, logFC = lfc[!is.na(lfc)])
}))
logfc_data$group <- factor(logfc_data$group,
  levels = c("original_PTC+", "reclassified_PTC", "no_downstream_ejc",
             "ref_atg_lost", "truncated_no_ejc"))

p2 <- ggplot(logfc_data, aes(x = group, y = logFC, fill = group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "NMD Response by Residual Category",
       subtitle = "mashr posterior mean logFC (5 conditions)",
       x = NULL, y = "Mean logFC") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figures/explore_phase0_logfc_by_group.png", p2, width = 10, height = 6, dpi = 300)
ggsave("figures/explore_phase0_logfc_by_group.pdf", p2, width = 10, height = 6)

cat("\nFigures saved to figures/explore_phase0_*.{pdf,png}\n")
cat("Done.\n")
