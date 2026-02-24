#!/usr/bin/env Rscript
# 08_poison_exon_structure.R
#
# Gene structure diagrams for SRSF1 and SRSF7 showing:
# - All annotated isoforms with exon structures
# - rMATS poison exon event region highlighted
# - Which isoforms contain/overlap/skip the poison exon
# - Expression status and logFC from PacBio DGE
#
# Run from: results/isoform_transitions/Version_6.0/

suppressPackageStartupMessages({
  library(tidyverse)
  library(maser)
  library(patchwork)
})

cat("\n")
cat("==================================================================\n")
cat("   08: Poison Exon Structure Diagrams\n")
cat("==================================================================\n\n")

fig_dir <- "results/ptc/figures"

# ── Load data ──
iso_str  <- readRDS("data/isoform_structures.rds")
rmats_jc <- readRDS("../../../results/rds/allrmats_jc_2025.8.20.rds")
ptc      <- readRDS("results/ptc/results/ptc_status.rds")

de_files <- list(
  DD = "../../../longread_dge/nmd_dge_dd_2026.1.18.csv",
  MV = "../../../longread_dge/nmd_dge_mv_2026.1.18.csv",
  DD_ALI = "../../../longread_dge/nmd_dge_ddali_2026.1.18.csv"
)

# ── Helper: build exon tibble for all isoforms of a gene ──
get_exon_df <- function(gene_id_val, de_gene_symbol) {
  isos <- iso_str %>% filter(gene_id == gene_id_val)

  # Get expression/DE status from DD
  dd_de <- read.csv(de_files[["DD"]])
  gene_de <- dd_de %>%
    filter(hgnc_id == de_gene_symbol | hgnc_id.gc == de_gene_symbol |
           hgnc_id.sq == de_gene_symbol)

  exon_rows <- list()
  for (j in seq_len(nrow(isos))) {
    iso <- isos[j, ]
    starts <- iso$exon_starts[[1]]
    ends   <- iso$exon_ends[[1]]
    # Check if expressed
    de_row <- gene_de %>% filter(txid == iso$isoform_id)
    expressed <- nrow(de_row) > 0
    logfc <- if (expressed) de_row$logFC[1] else NA_real_
    adjp  <- if (expressed) de_row$adj.P.Val[1] else NA_real_
    biotype <- if (expressed && "biotype" %in% colnames(de_row)) de_row$biotype[1] else ""
    # PTC status
    ptc_row <- ptc %>% filter(isoform_id == iso$isoform_id)
    has_ptc_val <- if (nrow(ptc_row) > 0 && !is.na(ptc_row$has_ptc[1])) ptc_row$has_ptc[1] else NA

    source <- if_else(grepl("^ENST", iso$isoform_id), "GENCODE", "PacBio")

    for (k in seq_along(starts)) {
      exon_rows[[length(exon_rows) + 1]] <- tibble(
        isoform_id = iso$isoform_id,
        source = source,
        strand = iso$strand,
        exon_num = k,
        n_exons = iso$n_exons,
        start = starts[k],
        end = ends[k],
        width = ends[k] - starts[k] + 1,
        expressed = expressed,
        logFC = logfc,
        adj.P.Val = adjp,
        biotype = biotype,
        has_ptc = has_ptc_val,
        tx_start = iso$tx_start,
        tx_end = iso$tx_end
      )
    }
  }
  bind_rows(exon_rows)
}

# ── Helper: get rMATS SE event coordinates for a gene ──
get_rmats_se <- function(gene_symbol) {
  results <- list()
  for (ct_rmats in names(rmats_jc)) {
    obj <- rmats_jc[[ct_rmats]]
    se_ev <- slot(obj, "SE_events") %>% as_tibble()
    se_st <- slot(obj, "SE_stats") %>% as_tibble()
    se_gr <- slot(obj, "SE_gr")

    idx <- which(se_ev$geneSymbol == gene_symbol)
    if (length(idx) == 0) next

    for (i in idx) {
      target <- se_gr[["exon_target"]][i]
      up     <- se_gr[["exon_upstream"]][i]
      down   <- se_gr[["exon_downstream"]][i]
      stats  <- se_st[i, ]

      results[[length(results) + 1]] <- tibble(
        cell_type = ct_rmats,
        event_id = se_ev$ID[i],
        chr = as.character(seqnames(target)),
        strand = as.character(strand(target)),
        target_start = start(target),
        target_end = end(target),
        target_width = width(target),
        up_start = start(up), up_end = end(up),
        down_start = start(down), down_end = end(down),
        dPSI = stats$IncLevelDifference,
        FDR = stats$FDR,
        frame_disrupting = (width(target) %% 3) != 0
      )
    }
  }
  bind_rows(results)
}

# ═══════════════════════════════════════════════════════════════════════════════
# SRSF1
# ═══════════════════════════════════════════════════════════════════════════════
cat("Building SRSF1 structure diagram...\n")

srsf1_exons <- get_exon_df("ENSG00000136450", "SRSF1")
srsf1_rmats <- get_rmats_se("SRSF1")

# The significant poison exon event in DD
srsf1_poison <- srsf1_rmats %>%
  filter(cell_type == "DD", FDR < 0.05, frame_disrupting) %>%
  arrange(FDR) %>%
  dplyr::slice(1)

cat(sprintf("  Poison exon: %d-%d (width %d)\n",
            srsf1_poison$target_start, srsf1_poison$target_end,
            srsf1_poison$target_width))

# Create isoform labels with expression info
iso_labels <- srsf1_exons %>%
  distinct(isoform_id, expressed, logFC, adj.P.Val, biotype, has_ptc, n_exons) %>%
  mutate(
    expr_label = case_when(
      !expressed ~ "not detected",
      adj.P.Val < 0.05 & logFC > 0 ~ sprintf("logFC=%+.2f *", logFC),
      adj.P.Val < 0.05 ~ sprintf("logFC=%+.2f *", logFC),
      TRUE ~ sprintf("logFC=%+.2f (ns)", logFC)
    ),
    ptc_label = case_when(
      is.na(has_ptc) ~ "",
      has_ptc ~ " [PTC+]",
      TRUE ~ ""
    ),
    label = paste0(isoform_id, "\n", biotype, ptc_label, "\n", expr_label)
  )

srsf1_exons <- srsf1_exons %>%
  left_join(iso_labels %>% select(isoform_id, label), by = "isoform_id") %>%
  mutate(
    # Check if this exon overlaps the poison exon
    overlaps_poison = start <= srsf1_poison$target_end &
                      end >= srsf1_poison$target_start,
    contains_poison = start <= srsf1_poison$target_start &
                      end >= srsf1_poison$target_end,
    is_poison_exon = start == srsf1_poison$target_start &
                     end == srsf1_poison$target_end
  )

# Order isoforms: expressed first, then by n_exons
iso_order <- srsf1_exons %>%
  distinct(isoform_id, label, expressed, n_exons) %>%
  arrange(desc(expressed), desc(n_exons)) %>%
  pull(label)

srsf1_exons$label <- factor(srsf1_exons$label, levels = rev(iso_order))

# Gene region for x-axis
gene_min <- min(srsf1_exons$start) - 200
gene_max <- max(srsf1_exons$end) + 200

p_srsf1 <- ggplot(srsf1_exons) +
  # Poison exon highlight region
  annotate("rect",
           xmin = srsf1_poison$target_start, xmax = srsf1_poison$target_end,
           ymin = -Inf, ymax = Inf,
           fill = "#E74C3C", alpha = 0.12) +
  # Intron lines (tx_start to tx_end)
  geom_segment(
    data = srsf1_exons %>% distinct(label, tx_start, tx_end, expressed),
    aes(x = tx_start, xend = tx_end, y = label, yend = label,
        linetype = expressed),
    linewidth = 0.3, color = "grey40"
  ) +
  # Exon rectangles
  geom_rect(
    aes(xmin = start, xmax = end,
        ymin = as.numeric(label) - 0.3,
        ymax = as.numeric(label) + 0.3,
        fill = case_when(
          is_poison_exon ~ "Exact poison exon",
          overlaps_poison & !is_poison_exon ~ "Overlaps poison region",
          TRUE ~ "Other exon"
        ),
        alpha = expressed)
  ) +
  # rMATS flanking exons (upstream and downstream)
  annotate("point",
           x = (srsf1_poison$up_start + srsf1_poison$up_end) / 2,
           y = 0.3, shape = 25, size = 2, fill = "#2ECC71") +
  annotate("point",
           x = (srsf1_poison$down_start + srsf1_poison$down_end) / 2,
           y = 0.3, shape = 25, size = 2, fill = "#2ECC71") +
  annotate("text",
           x = (srsf1_poison$target_start + srsf1_poison$target_end) / 2,
           y = 0.3,
           label = sprintf("rMATS SE: dPSI=%+.2f\nFDR=%.3f (DD)",
                           srsf1_poison$dPSI, srsf1_poison$FDR),
           size = 3, color = "#C0392B", fontface = "bold") +
  scale_fill_manual(
    values = c("Exact poison exon" = "#E74C3C",
               "Overlaps poison region" = "#F39C12",
               "Other exon" = "#3498DB"),
    name = "Exon type"
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4),
                     name = "Detected in PacBio",
                     labels = c("TRUE" = "Yes", "FALSE" = "No")) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
                        guide = "none") +
  scale_x_continuous(labels = scales::comma) +
  coord_cartesian(xlim = c(gene_min, gene_max)) +
  labs(title = "SRSF1 (chr17, minus strand) - poison exon at 58,005,398-58,005,600",
       subtitle = "Red shading = rMATS poison exon region (203 bp, not divisible by 3). Only 1 of 6 isoforms detected by PacBio.",
       x = "Genomic position (chr17)", y = "") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 7, family = "mono"),
    plot.title = element_text(face = "bold")
  )

# ═══════════════════════════════════════════════════════════════════════════════
# SRSF7
# ═══════════════════════════════════════════════════════════════════════════════
cat("Building SRSF7 structure diagram...\n")

srsf7_exons <- get_exon_df("ENSG00000115875", "SRSF7")
srsf7_rmats <- get_rmats_se("SRSF7")

# The strongest significant poison exon event in DD
srsf7_poison <- srsf7_rmats %>%
  filter(cell_type == "DD", FDR < 0.05, frame_disrupting) %>%
  arrange(FDR) %>%
  dplyr::slice(1)

cat(sprintf("  Poison exon: %d-%d (width %d)\n",
            srsf7_poison$target_start, srsf7_poison$target_end,
            srsf7_poison$target_width))

iso_labels7 <- srsf7_exons %>%
  distinct(isoform_id, expressed, logFC, adj.P.Val, biotype, has_ptc, n_exons) %>%
  mutate(
    expr_label = case_when(
      !expressed ~ "not detected",
      adj.P.Val < 0.05 & logFC > 0 ~ sprintf("logFC=%+.2f *", logFC),
      adj.P.Val < 0.05 ~ sprintf("logFC=%+.2f *", logFC),
      TRUE ~ sprintf("logFC=%+.2f (ns)", logFC)
    ),
    ptc_label = case_when(
      is.na(has_ptc) ~ "",
      has_ptc ~ " [PTC+]",
      TRUE ~ ""
    ),
    label = paste0(isoform_id, "\n", biotype, ptc_label, "\n", expr_label)
  )

srsf7_exons <- srsf7_exons %>%
  left_join(iso_labels7 %>% select(isoform_id, label), by = "isoform_id") %>%
  mutate(
    overlaps_poison = start <= srsf7_poison$target_end &
                      end >= srsf7_poison$target_start,
    contains_poison = start <= srsf7_poison$target_start &
                      end >= srsf7_poison$target_end,
    is_poison_exon = start == srsf7_poison$target_start &
                     end == srsf7_poison$target_end
  )

iso_order7 <- srsf7_exons %>%
  distinct(isoform_id, label, expressed, n_exons) %>%
  arrange(desc(expressed), desc(n_exons)) %>%
  pull(label)

srsf7_exons$label <- factor(srsf7_exons$label, levels = rev(iso_order7))

gene_min7 <- min(srsf7_exons$start) - 200
gene_max7 <- max(srsf7_exons$end) + 200

p_srsf7 <- ggplot(srsf7_exons) +
  annotate("rect",
           xmin = srsf7_poison$target_start, xmax = srsf7_poison$target_end,
           ymin = -Inf, ymax = Inf,
           fill = "#E74C3C", alpha = 0.12) +
  geom_segment(
    data = srsf7_exons %>% distinct(label, tx_start, tx_end, expressed),
    aes(x = tx_start, xend = tx_end, y = label, yend = label,
        linetype = expressed),
    linewidth = 0.3, color = "grey40"
  ) +
  geom_rect(
    aes(xmin = start, xmax = end,
        ymin = as.numeric(label) - 0.3,
        ymax = as.numeric(label) + 0.3,
        fill = case_when(
          is_poison_exon ~ "Exact poison exon",
          overlaps_poison & !is_poison_exon ~ "Overlaps poison region",
          TRUE ~ "Other exon"
        ),
        alpha = expressed)
  ) +
  annotate("text",
           x = (srsf7_poison$target_start + srsf7_poison$target_end) / 2,
           y = 0.3,
           label = sprintf("rMATS SE: dPSI=%+.2f\nFDR=%.3f (DD)",
                           srsf7_poison$dPSI, srsf7_poison$FDR),
           size = 3, color = "#C0392B", fontface = "bold") +
  scale_fill_manual(
    values = c("Exact poison exon" = "#E74C3C",
               "Overlaps poison region" = "#F39C12",
               "Other exon" = "#3498DB"),
    name = "Exon type"
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4),
                     name = "Detected in PacBio",
                     labels = c("TRUE" = "Yes", "FALSE" = "No")) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
                        guide = "none") +
  scale_x_continuous(labels = scales::comma) +
  coord_cartesian(xlim = c(gene_min7, gene_max7)) +
  labs(title = "SRSF7 (chr2, minus strand) - poison exon at 38,748,898-38,749,346",
       subtitle = "Red shading = rMATS poison exon region (449 bp, not divisible by 3). No isoform contains the full poison exon.",
       x = "Genomic position (chr2)", y = "") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 7, family = "mono"),
    plot.title = element_text(face = "bold")
  )

# ═══════════════════════════════════════════════════════════════════════════════
# Combine
# ═══════════════════════════════════════════════════════════════════════════════
cat("Assembling figure...\n")

combined <- p_srsf1 / p_srsf7 +
  plot_annotation(
    title = "Poison exon events in isoform context",
    subtitle = paste0(
      "rMATS detects poison exon inclusion changes not represented in PacBio isoform models.\n",
      "Blue = normal exons, Orange = overlaps poison region, Red = exact poison exon match.\n",
      "Faded = not detected in PacBio expression data. * = significant (adj.P < 0.05, DD)."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
  )

ggsave(file.path(fig_dir, "poison_exon_structures.pdf"),
       combined, width = 14, height = 12)
cat("  Wrote poison_exon_structures.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Console summary: compatibility assessment
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("==================================================================\n")
cat("   Compatibility Assessment\n")
cat("==================================================================\n\n")

cat("SRSF1 (ENSG00000136450):\n")
cat("  Total isoforms: 6 GENCODE, 0 PacBio novel\n")
cat("  rMATS poison exon: 58,005,398-58,005,600 (203 bp)\n")
cat("  Detected isoform: ENST00000583741.1 (GENCODE, NMD biotype)\n")
cat("    -> Contains poison region within larger exon (58,004,849-58,005,600)\n")
cat("    -> Very low expression (AveExpr ~ -3.7), not significant in any cell type\n")
cat("    -> The canonical poison exon isoform (ENST00000584668.5) is NOT detected\n")
cat("  No PacBio novel isoforms exist in the isoform catalog.\n")
cat("  ENST00000583741.1 (NMD biotype): AveExpr=-3.65 (~1 read/sample), not significant.\n")

cat("\nSRSF7 (ENSG00000115875):\n")
cat("  Total isoforms: 7 GENCODE, 0 PacBio novel\n")
cat("  rMATS poison exon: 38,748,898-38,749,346 (449 bp)\n")
cat("  No isoform (GENCODE or PacBio) contains the full poison exon.\n")
cat("    ENST00000415527.1 has a partial overlap (276 bp) but is NOT expressed.\n")
cat("  Detected isoforms: 6 (all GENCODE), all either skip the region or end before it.\n")
cat("  Conclusion: The poison exon-containing transcript is absent from\n")
cat("    both GENCODE annotations and PacBio novel isoform discovery.\n")

cat("\nSequencing depth context:\n")
cat("  Long-read: ~14M reads/sample (DD median), ~2000bp each = ~28 Gbp total bases\n")
cat("  Short-read: ~49M reads/sample (DD median), ~150bp each = ~7.4 Gbp total bases\n")
cat("  PacBio generates MORE total bases, but distributes across ~300K isoforms\n")
cat("  Short-read quantifies ~21K genes → ~58x more reads per feature at median\n")
cat("  SRSF1 poison exon event: ~1,400-2,100 junction reads/sample (rMATS)\n")
cat("    vs ~1 read/sample for NMD isoform (PacBio)\n")

cat("\n  Done.\n")
