#!/usr/bin/env Rscript
# Paper Figure: Isoform Expression Landscape (Manuscript Section 2 main figure)
#
# Panels:
#   A) Isoform length distribution by SQANTI structural category
#   B) Number of expressed isoforms per cell type (5% gene-proportion filter)
#   C) Cell-type sharing of NMD-responsive isoforms
#   D) % transcriptional output lost to NMD per cell type
#
# Output: paper_figures/FigX-IsoformLandscape_generated.{pdf,png}

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
source("visualization_functions.R")
source("analysis_functions.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

paper_theme <- theme_bw(base_size = 12) +
  theme(plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.title    = element_text(size = 11),
        axis.text     = element_text(size = 10),
        legend.text   = element_text(size = 10),
        legend.title  = element_text(size = 11),
        plot.margin   = margin(5, 8, 5, 5))

# Cell types (paper scope, 2026-04-29)
CTS <- c("AT", "DD", "FB", "MV")

ct_pal <- c("AT" = "#1b9e77",
            "DD" = "#d95f02",
            "FB" = "#7570b3",
            "MV" = "#e7298a")

# ---- Load shared inputs ------------------------------------------------------
expr_mat <- readRDS("data_mashr/expression_data.rds")            # CPM matrix
sample_metadata <- readRDS("data_mashr/sample_metadata.rds")
gene_map <- readRDS("data_mashr/gene_map.rds")
nmd_classification <- readRDS("data_mashr/nmd_classification.rds")

# SQANTI classification table (one row per isoform, providing length + category)
sqanti_path <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_classification.txt"
sq <- read.table(sqanti_path, sep = "\t", header = TRUE,
                 stringsAsFactors = FALSE, quote = "", comment.char = "")

# ============================================================================
# Panel A: Isoform length distribution by SQANTI structural category
# ============================================================================
# Use SQANTI `length` column (transcript length in nt) for isoforms in the
# 5%-filtered expression matrix. Map ugly SQANTI category names to short labels.

cat_map <- c("full-splice_match"       = "FSM",
             "incomplete-splice_match" = "ISM",
             "novel_in_catalog"        = "NIC",
             "novel_not_in_catalog"    = "NNC",
             "fusion"                  = "Fusion",
             "genic"                   = "Genic")

panelA_df <- sq %>%
  filter(isoform %in% rownames(expr_mat)) %>%
  transmute(isoform_id  = isoform,
            length      = length,
            category    = unname(cat_map[structural_category])) %>%
  filter(!is.na(category), length > 0) %>%
  mutate(log10_length = log10(length))

# Order categories by descending count for legend readability
cat_levels <- panelA_df %>%
  count(category, sort = TRUE) %>%
  pull(category)
panelA_df$category <- factor(panelA_df$category, levels = cat_levels)

# Color palette (qualitative, distinct for 6 SQANTI categories)
sq_pal <- c("FSM"    = "#1f78b4",
            "ISM"    = "#a6cee3",
            "NIC"    = "#33a02c",
            "NNC"    = "#b2df8a",
            "Fusion" = "#e31a1c",
            "Genic"  = "#fdbf6f")

panelA_counts <- panelA_df %>%
  count(category) %>%
  mutate(label = sprintf("%s (n=%s)", category, format(n, big.mark = ",")))

panelA_legend_labels <- setNames(panelA_counts$label, panelA_counts$category)

pA <- ggplot(panelA_df, aes(x = log10_length, color = category)) +
  geom_density(linewidth = 0.8, alpha = 0.85) +
  scale_color_manual(values = sq_pal, labels = panelA_legend_labels) +
  labs(title = "Isoform Length Distribution by SQANTI Category",
       subtitle = sprintf("n = %s expressed isoforms",
                          format(nrow(panelA_df), big.mark = ",")),
       x = expression(log[10]~"(transcript length, nt)"),
       y = "Density",
       color = "SQANTI category") +
  paper_theme +
  theme(legend.position  = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.key.height = unit(0.6, "lines"))

# ============================================================================
# Panel B: Expressed isoforms per cell type (5% within-CT gene proportion filter)
# ============================================================================
# Per CT, isoform passes if its max gene-proportion in either CT-DMSO or
# CT-Smg1i samples is >= 5%. This mirrors the global filter in
# 01_prepare_data_mashr.R but scoped to within-CT samples.

gene_for_iso <- setNames(gene_map$gene_id, gene_map$isoform_id)[rownames(expr_mat)]

compute_max_prop_ct <- function(samples) {
  cpm_sub <- expr_mat[, samples, drop = FALSE]
  gene_totals <- rowsum(cpm_sub, gene_for_iso)
  props <- cpm_sub / gene_totals[gene_for_iso, , drop = FALSE]
  props[is.nan(props)] <- 0
  apply(props, 1, max, na.rm = TRUE)
}

count_expressed_per_ct <- sapply(CTS, function(ct) {
  d_s <- sample_metadata$sample_id[sample_metadata$ct == ct &
                                     sample_metadata$treatment == "DMSO"]
  s_s <- sample_metadata$sample_id[sample_metadata$ct == ct &
                                     sample_metadata$treatment == "Smg1i"]
  max_d <- compute_max_prop_ct(d_s)
  max_s <- compute_max_prop_ct(s_s)
  sum(max_d >= 0.05 | max_s >= 0.05)
})

panelB_df <- data.frame(
  ct   = factor(CTS, levels = CTS),
  n    = unname(count_expressed_per_ct)
)

pB <- ggplot(panelB_df, aes(x = ct, y = n, fill = ct)) +
  geom_col(width = 0.7, alpha = 0.9, color = "grey20") +
  geom_text(aes(label = format(n, big.mark = ",")),
            vjust = -0.4, size = 3.6, fontface = "bold") +
  scale_fill_manual(values = ct_pal, guide = "none") +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Expressed Isoforms per Cell Type",
       subtitle = "5% gene-proportion filter (DMSO or Smg1i, within CT)",
       x = "Cell type",
       y = "Number of expressed isoforms") +
  paper_theme

# ============================================================================
# Panel C: Cell-type sharing of NMD-responsive isoforms
# ============================================================================
# For each isoform, count how many of the 4 CTs flag it as NMD-responsive
# (mashr nmd_responsive == TRUE). Plot a stacked bar showing the per-CT
# composition by sharing class (CT-specific, 2 CTs, 3 CTs, all 4 CTs).

de_dir <- "/Users/petecastaldi/claude_projects/nmd/isocall_dge/mashr"
ct_to_de_suffix <- c(AT = "at", DD = "dd", FB = "fb", MV = "mv")

per_ct_nmd <- lapply(CTS, function(ct) {
  fp <- file.path(de_dir,
                  sprintf("nmd_mashr_die_%s_2026.3.10.csv", ct_to_de_suffix[ct]))
  d <- read.csv(fp, stringsAsFactors = FALSE)
  d$txid[d$nmd_responsive == TRUE]
})
names(per_ct_nmd) <- CTS

# Matrix of presence/absence in each CT for the union of NMD isoforms
all_nmd <- sort(unique(unlist(per_ct_nmd)))
share_mat <- sapply(per_ct_nmd, function(ids) all_nmd %in% ids)
rownames(share_mat) <- all_nmd
n_cts_per_iso <- rowSums(share_mat)

share_label <- function(k) {
  if (k == 1) "CT-specific (1 CT)"
  else if (k == length(CTS)) sprintf("All %d CTs", length(CTS))
  else sprintf("Shared in %d CTs", k)
}
share_levels <- c(sapply(seq_len(length(CTS)), share_label))

# For each CT, count its NMD isoforms by sharing class
panelC_df <- do.call(rbind, lapply(CTS, function(ct) {
  ids <- per_ct_nmd[[ct]]
  k   <- n_cts_per_iso[ids]
  data.frame(ct       = ct,
             share    = factor(sapply(k, share_label), levels = share_levels),
             stringsAsFactors = FALSE)
})) %>%
  count(ct, share) %>%
  group_by(ct) %>%
  mutate(total = sum(n),
         pct   = 100 * n / total) %>%
  ungroup() %>%
  mutate(ct = factor(ct, levels = CTS))

share_pal <- c("CT-specific (1 CT)" = "#fdae61",
               "Shared in 2 CTs"    = "#abdda4",
               "Shared in 3 CTs"    = "#66c2a5",
               "All 4 CTs"          = "#3288bd")

# Total NMD count for each CT (annotation)
ct_totals <- panelC_df %>%
  group_by(ct) %>%
  summarise(total = first(total), .groups = "drop")

pC <- ggplot(panelC_df, aes(x = ct, y = n, fill = share)) +
  geom_col(width = 0.7, color = "grey20", alpha = 0.9) +
  geom_text(data = ct_totals,
            aes(x = ct, y = total, label = format(total, big.mark = ",")),
            inherit.aes = FALSE,
            vjust = -0.4, size = 3.6, fontface = "bold") +
  scale_fill_manual(values = share_pal) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Cell-Type Sharing of NMD-Responsive Isoforms",
       subtitle = sprintf("Union across CTs: %s isoforms",
                          format(length(all_nmd), big.mark = ",")),
       x = "Cell type",
       y = "NMD-responsive isoforms",
       fill = "Sharing") +
  paper_theme +
  theme(legend.position  = "bottom",
        legend.title     = element_text(size = 10),
        legend.key.size  = unit(0.7, "lines"))

# ============================================================================
# Panel D: % transcriptional output lost to NMD per cell type
# ============================================================================
# For each CT:  pct = sum_{iso in NMD_ct} (mean Smg1i CPM - mean DMSO CPM)
#               -------------------------------------------------------------
#                              sum_{all iso} mean Smg1i CPM
#
# Numerator = transcriptional rescue when NMD is blocked = the output that
# NMD removes in DMSO. Denominator = total Smg1i transcriptional output.

panelD_rows <- lapply(CTS, function(ct) {
  d_s <- sample_metadata$sample_id[sample_metadata$ct == ct &
                                     sample_metadata$treatment == "DMSO"]
  s_s <- sample_metadata$sample_id[sample_metadata$ct == ct &
                                     sample_metadata$treatment == "Smg1i"]
  m_d <- rowMeans(expr_mat[, d_s, drop = FALSE])
  m_s <- rowMeans(expr_mat[, s_s, drop = FALSE])
  nmd_ids <- intersect(nmd_classification[[ct]]$nmd, rownames(expr_mat))
  num <- sum(m_s[nmd_ids] - m_d[nmd_ids])
  denom <- sum(m_s)
  data.frame(ct = ct,
             pct_lost = 100 * num / denom,
             n_nmd    = length(nmd_ids),
             stringsAsFactors = FALSE)
})
panelD_df <- do.call(rbind, panelD_rows)
panelD_df$ct <- factor(panelD_df$ct, levels = CTS)

pD <- ggplot(panelD_df, aes(x = ct, y = pct_lost, fill = ct)) +
  geom_col(width = 0.7, alpha = 0.9, color = "grey20") +
  geom_text(aes(label = sprintf("%.1f%%", pct_lost)),
            vjust = -0.4, size = 3.6, fontface = "bold") +
  scale_fill_manual(values = ct_pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "Transcriptional Output Lost to NMD",
       subtitle = expression(("Smg1i CPM" - "DMSO CPM")[NMD] / "Total Smg1i CPM"),
       x = "Cell type",
       y = "% of total Smg1i output") +
  paper_theme

# ============================================================================
# Compose
# ============================================================================
combined <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.margin = margin(10, 10, 10, 10))) &
  theme(plot.tag = element_text(size = 18, face = "bold"))

# ---- Validate ---------------------------------------------------------------
source("~/.claude/utils/validate_composite_layout.R")
val <- validate_composite_layout(combined, fig_width = 14, fig_height = 10)
if (val$summary$n_errors > 0) {
  stop("Composite layout validation failed; see messages above.")
}

# ---- Save -------------------------------------------------------------------
ggsave("paper_figures/FigX-IsoformLandscape_generated.png", combined,
       width = 14, height = 10, dpi = 300, bg = "white")
ggsave("paper_figures/FigX-IsoformLandscape_generated.pdf", combined,
       width = 14, height = 10, bg = "white")

cat("\nSaved: paper_figures/FigX-IsoformLandscape_generated.{png,pdf}\n")

# ---- Summary print-out ------------------------------------------------------
cat("\n=== Panel data summary ===\n")
cat("Panel A (length distribution):\n")
print(panelA_counts)
cat("\nPanel B (expressed isoforms per CT):\n"); print(panelB_df)
cat("\nPanel C (CT sharing — totals per CT and shared union):\n")
print(ct_totals); cat(sprintf("Union across CTs: %d\n", length(all_nmd)))
cat("\nPanel D (% transcriptional output lost):\n"); print(panelD_df)
