#!/usr/bin/env Rscript
# ==============================================================================
# Phase 1: Specific ORF vs 5'UTR landscape investigation
#
# Q1: Can any specific ORF explain NMD susceptibility?
#     → For each no_downstream_ejc isoform, find ORFs with downstream EJCs
#       and strong Kozak context. If one dominant ORF explains most cases,
#       it's a specific-ORF mechanism.
#
# Q2: What's different about these 5'UTRs compared to Control and PTC+?
#     → Compare 5'UTR features (ATG density, uORF counts, Kozak, coverage)
#       across groups WITHOUT length-matching. The long 5'UTR is the phenotype.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
fmt <- function(x) format(x, big.mark = ",")

cat("=== Phase 1: Specific ORF vs 5'UTR Landscape ===\n\n")

# Load data
ref_atg <- readRDS("data_mashr/analysis_cache/ref_atg_analysis.rds")
orfik_data <- readRDS("data_mashr/analysis_cache/orfik_scan.rds")
utr5_all <- readRDS("data_mashr/analysis_cache/utr5_features_all.rds")
cds <- readRDS("data_mashr/cds.rds")

orf_feat <- orfik_data$orf_features
orf_tx <- orfik_data$tx_summary

# Define populations
c2_res <- ref_atg$c2
ptc_neg <- c2_res[!c2_res$original_ptc, ]
no_dejc <- ptc_neg[ptc_neg$category == "no_downstream_ejc", ]
eff_ptc <- ptc_neg[ptc_neg$category == "effectively_ptc", ]
ptc_pos <- c2_res[c2_res$original_ptc, ]

c4_res <- ref_atg$c4
c4_avail <- c4_res[c4_res$ref_atg_exonic_in_comp == TRUE &
                     !is.na(c4_res$comp_stop_tx_pos), ]

cat("Populations:\n")
cat("  no_downstream_ejc:", nrow(no_dejc), "\n")
cat("  original PTC+:", nrow(ptc_pos), "\n")
cat("  reclassified PTC:", nrow(eff_ptc), "\n")
cat("  Control (C4):", nrow(c4_avail), "\n\n")

# ==============================================================================
# Q1: Per-ORF analysis for no_downstream_ejc
# ==============================================================================
cat("=== Q1: Can specific ORFs explain NMD? ===\n\n")

# Get all ORFs for no_downstream_ejc isoforms
no_dejc_ids <- no_dejc$comparator_isoform_id
no_dejc_orfs <- orf_feat[orf_feat$isoform_id %in% no_dejc_ids, ]

cat("ORFs in no_downstream_ejc isoforms:", nrow(no_dejc_orfs), "\n")
cat("  with downstream EJC:", sum(no_dejc_orfs$has_downstream_ejc), "\n")
cat("  with strong Kozak (score=2) + dEJC:",
    sum(no_dejc_orfs$has_downstream_ejc & no_dejc_orfs$kozak_score == 2), "\n")
cat("  with any Kozak (score>=1) + dEJC:",
    sum(no_dejc_orfs$has_downstream_ejc & no_dejc_orfs$kozak_score >= 1), "\n\n")

# Per-isoform: how many NMD-susceptible ORFs?
no_dejc_nmd_orfs <- no_dejc_orfs %>%
  filter(has_downstream_ejc) %>%
  group_by(isoform_id) %>%
  summarise(
    n_nmd_orfs = n(),
    n_strong_kozak = sum(kozak_score == 2),
    n_any_kozak = sum(kozak_score >= 1),
    max_kozak = max(kozak_score),
    .groups = "drop"
  )

cat("Per-isoform NMD-susceptible ORF counts:\n")
cat("  Isoforms with >= 1 NMD-susceptible ORF:",
    nrow(no_dejc_nmd_orfs), "/", length(no_dejc_ids),
    sprintf("(%.1f%%)\n", 100 * nrow(no_dejc_nmd_orfs) / length(no_dejc_ids)))
cat("  Median NMD-susceptible ORFs per isoform:",
    median(no_dejc_nmd_orfs$n_nmd_orfs), "\n")
cat("  Isoforms with strong Kozak NMD ORF:",
    sum(no_dejc_nmd_orfs$n_strong_kozak > 0),
    sprintf("(%.1f%%)\n", 100 * mean(no_dejc_nmd_orfs$n_strong_kozak > 0)), "\n")

# Is there typically ONE dominant NMD ORF or many?
cat("Distribution of NMD-susceptible ORF count per isoform:\n")
print(table(pmin(no_dejc_nmd_orfs$n_nmd_orfs, 10)))

# Where are these ORFs relative to the main CDS?
# Check: are they upstream of the CDS start (true uORFs) or alternative CDS?
no_dejc_orfs_with_ejc <- no_dejc_orfs[no_dejc_orfs$has_downstream_ejc, ]

# Get CDS start position for each isoform
cds_starts <- cds %>%
  filter(isoform_id %in% no_dejc_ids, coding_status == "coding") %>%
  transmute(isoform_id,
            cds_5prime_tx = NA_integer_)  # we need transcript-space position

# Actually, orf_feat has orf_start (transcript position) — compare to CDS start
# The ORFik orf_start is in transcript coordinates. The CDS start from TD2 is
# also expressible in transcript coordinates via the 5'UTR length.
# If orf_start < 5'UTR length, it's a uORF; if >= 5'UTR length, it's in/past CDS.

utr5_lookup <- utr5_all$isoform_features %>%
  filter(isoform_id %in% no_dejc_ids) %>%
  select(isoform_id, utr5_length)

no_dejc_orfs_with_ejc <- no_dejc_orfs_with_ejc %>%
  left_join(utr5_lookup, by = "isoform_id") %>%
  mutate(is_upstream = orf_start < utr5_length)

cat("\nNMD-susceptible ORFs: upstream vs downstream of CDS start:\n")
cat("  Upstream (true uORFs):", sum(no_dejc_orfs_with_ejc$is_upstream, na.rm = TRUE),
    sprintf("(%.1f%%)\n", 100 * mean(no_dejc_orfs_with_ejc$is_upstream, na.rm = TRUE)))
cat("  At/past CDS start:", sum(!no_dejc_orfs_with_ejc$is_upstream, na.rm = TRUE), "\n")

# Among upstream NMD-susceptible ORFs, Kozak distribution
upstream_orfs <- no_dejc_orfs_with_ejc[no_dejc_orfs_with_ejc$is_upstream == TRUE, ]
cat("\nUpstream NMD-susceptible ORF Kozak distribution:\n")
print(table(upstream_orfs$kozak_score))

# ==============================================================================
# Q2: 5'UTR features across groups
# ==============================================================================
cat("\n=== Q2: 5'UTR Feature Comparison ===\n\n")

# Get 5'UTR features for each group
utr5_feat <- utr5_all$isoform_features %>%
  filter(!excluded)

groups <- list(
  "no_downstream_ejc" = no_dejc$comparator_isoform_id,
  "original_PTC+" = ptc_pos$comparator_isoform_id,
  "reclassified_PTC" = eff_ptc$comparator_isoform_id,
  "Control" = c4_avail$comparator_isoform_id
)

# Features to compare
features <- c("utr5_length", "n_atg", "atg_density", "n_strong_kozak_atg",
              "n_atg_without_orf", "n_orfs_overlapping", "n_orfs_inframe",
              "n_orfs_outframe", "longest_orf_nt", "pct_utr5_in_orfs",
              "stop_density")

feature_labels <- c(
  utr5_length = "5'UTR length (bp)",
  n_atg = "ATG count",
  atg_density = "ATG density (per 100 bp)",
  n_strong_kozak_atg = "Strong Kozak ATG count",
  n_atg_without_orf = "Orphan ATG count",
  n_orfs_overlapping = "Overlapping uORF count",
  n_orfs_inframe = "In-frame uORF count",
  n_orfs_outframe = "Out-of-frame uORF count",
  longest_orf_nt = "Longest uORF (nt)",
  pct_utr5_in_orfs = "5'UTR ORF coverage (%)",
  stop_density = "Stop density (per 100 bp)"
)

# Build comparison table
comparison_rows <- list()
for (grp_name in names(groups)) {
  sub <- utr5_feat[utr5_feat$isoform_id %in% groups[[grp_name]], ]
  row <- data.frame(Group = grp_name, n = nrow(sub), stringsAsFactors = FALSE)
  for (feat in features) {
    vals <- sub[[feat]]
    row[[feature_labels[feat]]] <- sprintf("%.1f [%.1f, %.1f]",
      median(vals, na.rm = TRUE),
      quantile(vals, 0.25, na.rm = TRUE),
      quantile(vals, 0.75, na.rm = TRUE))
  }
  comparison_rows[[grp_name]] <- row
}
comparison_df <- do.call(rbind, comparison_rows)

cat("5'UTR Feature Comparison (median [IQR]):\n\n")
for (i in seq_len(nrow(comparison_df))) {
  cat(sprintf("--- %s (n=%s) ---\n", comparison_df$Group[i], fmt(comparison_df$n[i])))
  for (feat in names(feature_labels)) {
    cat(sprintf("  %s: %s\n", feature_labels[feat],
                comparison_df[[feature_labels[feat]]][i]))
  }
  cat("\n")
}

# Wilcoxon tests: no_downstream_ejc vs Control
cat("=== Wilcoxon tests: no_downstream_ejc vs Control ===\n")
no_dejc_utr5 <- utr5_feat[utr5_feat$isoform_id %in% groups[["no_downstream_ejc"]], ]
ctrl_utr5 <- utr5_feat[utr5_feat$isoform_id %in% groups[["Control"]], ]

for (feat in features) {
  wt <- wilcox.test(no_dejc_utr5[[feat]], ctrl_utr5[[feat]])
  effect <- median(no_dejc_utr5[[feat]], na.rm = TRUE) - median(ctrl_utr5[[feat]], na.rm = TRUE)
  cat(sprintf("  %s: p = %s, median diff = %.1f\n",
              feature_labels[feat], signif(wt$p.value, 3), effect))
}

# Wilcoxon tests: no_downstream_ejc vs PTC+
cat("\n=== Wilcoxon tests: no_downstream_ejc vs PTC+ ===\n")
ptcpos_utr5 <- utr5_feat[utr5_feat$isoform_id %in% groups[["original_PTC+"]], ]

for (feat in features) {
  wt <- wilcox.test(no_dejc_utr5[[feat]], ptcpos_utr5[[feat]])
  effect <- median(no_dejc_utr5[[feat]], na.rm = TRUE) - median(ptcpos_utr5[[feat]], na.rm = TRUE)
  cat(sprintf("  %s: p = %s, median diff = %.1f\n",
              feature_labels[feat], signif(wt$p.value, 3), effect))
}

# ==============================================================================
# Summary figure: key distinguishing features
# ==============================================================================
cat("\n=== Building figures ===\n")

theme_set(theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank()))

# Build long-format data for selected features
key_feats <- c("utr5_length", "n_atg", "atg_density", "n_strong_kozak_atg",
               "n_orfs_overlapping", "longest_orf_nt", "pct_utr5_in_orfs")

plot_data <- do.call(rbind, lapply(names(groups), function(grp_name) {
  sub <- utr5_feat[utr5_feat$isoform_id %in% groups[[grp_name]], ]
  sub$group <- grp_name
  sub[, c("group", key_feats)]
}))
plot_data$group <- factor(plot_data$group,
  levels = c("Control", "original_PTC+", "reclassified_PTC", "no_downstream_ejc"))

plot_long <- plot_data %>%
  pivot_longer(cols = all_of(key_feats), names_to = "feature", values_to = "value") %>%
  mutate(feature = feature_labels[feature],
         feature = factor(feature, levels = feature_labels[key_feats]))

fig <- ggplot(plot_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
  facet_wrap(~feature, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Control" = "#67a9cf", "original_PTC+" = "#fddbc7",
                                "reclassified_PTC" = "#ef8a62",
                                "no_downstream_ejc" = "#b2182b")) +
  labs(title = "5'UTR Features: no_downstream_ejc vs Other NMD Groups",
       subtitle = "TD2 CDS-defined 5'UTR. no_downstream_ejc has confirmed genuine long 5'UTR.",
       x = NULL) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        strip.text = element_text(size = 10))

ggsave("figures/explore_phase1_utr5_features.png", fig, width = 14, height = 8, dpi = 300)
ggsave("figures/explore_phase1_utr5_features.pdf", fig, width = 14, height = 8)

cat("Saved figures/explore_phase1_utr5_features.{pdf,png}\n")
cat("Done.\n")
