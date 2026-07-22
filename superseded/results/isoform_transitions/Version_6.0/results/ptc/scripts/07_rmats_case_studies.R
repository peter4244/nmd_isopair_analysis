#!/usr/bin/env Rscript
# 07_rmats_case_studies.R
#
# Case study figures for canonical NMD-coupled poison exon genes,
# showing the disconnect between rMATS event-level signal and
# isoform-level PTC analysis.
#
# Run from: results/isoform_transitions/Version_6.0/

suppressPackageStartupMessages({
  library(tidyverse)
  library(maser)
  library(patchwork)
})

cat("\n")
cat("==================================================================\n")
cat("   07: Case Study Figures - rMATS vs Isoform-Level\n")
cat("==================================================================\n\n")

fig_dir <- "results/ptc/figures"

# ═══════════════════════════════════════════════════════════════════════════════
# Load data
# ═══════════════════════════════════════════════════════════════════════════════

cat("Loading data...\n")

rmats_jc <- readRDS("../../../results/rds/allrmats_jc_2025.8.20.rds")
ptc <- readRDS("results/ptc/results/ptc_status.rds")

de_files <- list(
  AT = "../../../longread_dge/nmd_dge_at2_2026.1.18.csv",
  DD = "../../../longread_dge/nmd_dge_dd_2026.1.18.csv",
  DD_ALI = "../../../longread_dge/nmd_dge_ddali_2026.1.18.csv",
  DO = "../../../longread_dge/nmd_dge_doali_2026.1.18.csv",
  FB = "../../../longread_dge/nmd_dge_fb_2026.1.18.csv",
  MV = "../../../longread_dge/nmd_dge_mv_2026.1.18.csv"
)

CT_MAP <- c(AT2 = "AT", DD = "DD", DD_ALI = "DD_ALI",
            DO_ALI = "DO", FB = "FB", MV = "MV")

# ═══════════════════════════════════════════════════════════════════════════════
# Helper: extract rMATS PSI data for a gene
# ═══════════════════════════════════════════════════════════════════════════════

get_rmats_psi <- function(gene_symbol, event_type = "SE") {
  results <- list()
  for (ct_rmats in names(rmats_jc)) {
    obj <- rmats_jc[[ct_rmats]]
    ct <- CT_MAP[ct_rmats]

    events <- slot(obj, paste0(event_type, "_events")) %>% as_tibble()
    stats  <- slot(obj, paste0(event_type, "_stats"))  %>% as_tibble()
    psi    <- slot(obj, paste0(event_type, "_PSI"))

    gene_idx <- which(events$geneSymbol == gene_symbol)
    if (length(gene_idx) == 0) next

    merged <- inner_join(events[gene_idx, ], stats[gene_idx, ], by = "ID")

    # Get PSI values per sample
    n1 <- slot(obj, "n_cond1")
    n2 <- slot(obj, "n_cond2")
    psi_mat <- psi[gene_idx, , drop = FALSE]

    for (i in seq_len(nrow(merged))) {
      row <- merged[i, ]
      psi_vals <- as.numeric(psi_mat[i, ])
      # cond1 = Smg1i, cond2 = DMSO
      smg1i_psi <- psi_vals[1:n1]
      dmso_psi  <- psi_vals[(n1 + 1):(n1 + n2)]

      smg1i_df <- tibble(
        cell_type = ct, event_id = row$ID, treatment = "Smg1i",
        PSI = smg1i_psi[!is.na(smg1i_psi)],
        dPSI = row$IncLevelDifference, FDR = row$FDR
      )
      dmso_df <- tibble(
        cell_type = ct, event_id = row$ID, treatment = "DMSO",
        PSI = dmso_psi[!is.na(dmso_psi)],
        dPSI = row$IncLevelDifference, FDR = row$FDR
      )
      results[[paste(ct, row$ID)]] <- bind_rows(smg1i_df, dmso_df)
    }
  }
  bind_rows(results)
}

# ═══════════════════════════════════════════════════════════════════════════════
# Helper: get isoform-level data for a gene
# ═══════════════════════════════════════════════════════════════════════════════

get_isoform_data <- function(gene_symbol) {
  results <- list()
  for (ct_name in names(de_files)) {
    de <- read.csv(de_files[[ct_name]])
    gene_de <- de %>%
      filter(hgnc_id == gene_symbol | hgnc_id.gc == gene_symbol |
             hgnc_id.sq == gene_symbol)
    if (nrow(gene_de) == 0) next

    gene_de <- gene_de %>%
      left_join(ptc %>% select(isoform_id, has_ptc, ptc_distance),
                by = c("txid" = "isoform_id")) %>%
      mutate(
        cell_type = ct_name,
        source = if_else(grepl("^ENST", txid), "GENCODE", "PacBio"),
        significant = adj.P.Val < 0.05,
        nmd_responsive = adj.P.Val < 0.05 & logFC > 0,
        ptc_label = case_when(
          has_ptc == TRUE ~ "PTC+",
          has_ptc == FALSE ~ "PTC-",
          TRUE ~ "unknown"
        )
      )
    results[[ct_name]] <- gene_de
  }
  bind_rows(results)
}

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE: SRSF1 case study
# ═══════════════════════════════════════════════════════════════════════════════

cat("Generating SRSF1 case study...\n")

# Get SRSF1 SE events - filter to the poison exon region
srsf1_psi <- get_rmats_psi("SRSF1", "SE")

# Find the significant events in DD (the poison exon)
srsf1_sig <- srsf1_psi %>%
  filter(FDR < 0.05) %>%
  distinct(cell_type, event_id, dPSI, FDR)
cat("  SRSF1 significant SE events:\n")
print(as.data.frame(srsf1_sig))

# Pick the strongest event per cell type for the plot
srsf1_top <- srsf1_psi %>%
  filter(event_id %in% srsf1_sig$event_id) %>%
  # Use only the first significant event per cell type to avoid overplotting
  group_by(cell_type, treatment) %>%
  filter(event_id == min(event_id)) %>%
  ungroup()

# Panel A: rMATS PSI by treatment for SRSF1 poison exon
p_srsf1_rmats <- ggplot(srsf1_top,
                         aes(x = treatment, y = PSI, color = treatment)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, linewidth = 0.5) +
  facet_wrap(~cell_type, nrow = 1) +
  scale_color_manual(values = c("DMSO" = "#3498DB", "Smg1i" = "#E74C3C")) +
  labs(title = "SRSF1 poison exon - rMATS SE event",
       subtitle = "Significant SE events (FDR < 0.05). dPSI < 0 = less inclusion with Smg1i.",
       y = "PSI (inclusion level)", x = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# Panel B: Isoform-level logFC for SRSF1
srsf1_iso <- get_isoform_data("SRSF1")

p_srsf1_iso <- ggplot(srsf1_iso,
                        aes(x = cell_type, y = logFC, color = ptc_label)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, lty = 2, color = "grey50") +
  scale_color_manual(values = c("PTC+" = "#E74C3C", "PTC-" = "#3498DB",
                                 "unknown" = "grey60")) +
  labs(title = "SRSF1 - PacBio isoform-level DGE",
       subtitle = sprintf("Only %d isoform detected. None significant (all adj.P.Val > 0.05).",
                           n_distinct(srsf1_iso$txid)),
       y = "logFC (Smg1i vs DMSO)", x = "", color = "PTC status") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE: SRSF7 case study
# ═══════════════════════════════════════════════════════════════════════════════

cat("Generating SRSF7 case study...\n")

srsf7_psi <- get_rmats_psi("SRSF7", "SE")
srsf7_sig <- srsf7_psi %>%
  filter(FDR < 0.05) %>%
  distinct(cell_type, event_id, dPSI, FDR)
cat("  SRSF7 significant SE events:\n")
print(as.data.frame(srsf7_sig))

srsf7_top <- srsf7_psi %>%
  filter(event_id %in% srsf7_sig$event_id) %>%
  group_by(cell_type, treatment) %>%
  filter(event_id == min(event_id)) %>%
  ungroup()

p_srsf7_rmats <- ggplot(srsf7_top,
                          aes(x = treatment, y = PSI, color = treatment)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, linewidth = 0.5) +
  facet_wrap(~cell_type, nrow = 1) +
  scale_color_manual(values = c("DMSO" = "#3498DB", "Smg1i" = "#E74C3C")) +
  labs(title = "SRSF7 poison exon - rMATS SE event",
       subtitle = "Significant SE events (FDR < 0.05). dPSI > 0 = more inclusion with Smg1i.",
       y = "PSI (inclusion level)", x = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

srsf7_iso <- get_isoform_data("SRSF7")

# Show all isoforms, highlight significance
p_srsf7_iso <- ggplot(srsf7_iso,
                        aes(x = cell_type, y = logFC, color = ptc_label,
                            shape = significant)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, lty = 2, color = "grey50") +
  scale_color_manual(values = c("PTC+" = "#E74C3C", "PTC-" = "#3498DB",
                                 "unknown" = "grey60")) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16)) +
  labs(title = "SRSF7 - PacBio isoform-level DGE",
       subtitle = sprintf("%d isoforms per cell type. Triangles = significant (adj.P.Val < 0.05).",
                           n_distinct(srsf7_iso$txid)),
       y = "logFC (Smg1i vs DMSO)", x = "",
       color = "PTC status", shape = "Significant") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ═══════════════════════════════════════════════════════════════════════════════
# Combine and save
# ═══════════════════════════════════════════════════════════════════════════════

cat("Assembling figure...\n")

combined <- (p_srsf1_rmats + p_srsf1_iso) /
            (p_srsf7_rmats + p_srsf7_iso) +
  plot_annotation(
    title = "NMD-coupled poison exons: rMATS event-level vs PacBio isoform-level",
    subtitle = "rMATS detects clear splicing changes at known poison exons; PacBio isoform-level analysis misses them",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
  )

ggsave(file.path(fig_dir, "rmats_case_studies.pdf"),
       combined, width = 14, height = 10)
cat("  Wrote rmats_case_studies.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Summary table for multiple genes
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Summary: rMATS vs isoform-level for NMD autoregulatory genes ──\n\n")

target_genes <- c("SRSF1", "SRSF3", "SRSF5", "SRSF6", "SRSF7",
                   "HNRNPL", "PTBP1", "TRA2B")

# Load pre-computed rMATS classification
rmats_class <- read_tsv("results/ptc/results/rmats_event_classification.tsv",
                         show_col_types = FALSE)

summary_rows <- list()
for (gene in target_genes) {
  # rMATS
  gene_rmats <- rmats_class %>% filter(geneSymbol == gene)
  n_sig_fd <- gene_rmats %>%
    filter(frame_class == "cds_frame_disrupting", FDR < 0.05) %>%
    nrow()
  max_dpsi <- gene_rmats %>%
    filter(frame_class == "cds_frame_disrupting", FDR < 0.05) %>%
    pull(IncLevelDifference) %>%
    {if (length(.) > 0) .[which.max(abs(.))] else NA_real_}

  # Isoform-level
  gene_iso <- tryCatch(get_isoform_data(gene), error = function(e) tibble())
  n_isoforms <- if (nrow(gene_iso) > 0 && "txid" %in% colnames(gene_iso))
    n_distinct(gene_iso$txid) else 0L
  n_sig_iso <- if (nrow(gene_iso) > 0 && "significant" %in% colnames(gene_iso))
    sum(gene_iso$significant, na.rm = TRUE) else 0L
  n_nmd_iso <- if (nrow(gene_iso) > 0 && "nmd_responsive" %in% colnames(gene_iso))
    sum(gene_iso$nmd_responsive, na.rm = TRUE) else 0L
  n_ptc_pos <- if (nrow(gene_iso) > 0 && "has_ptc" %in% colnames(gene_iso))
    sum(gene_iso$has_ptc == TRUE, na.rm = TRUE) else 0L

  summary_rows[[gene]] <- tibble(
    gene = gene,
    rmats_sig_fd = n_sig_fd,
    rmats_max_dpsi = round(max_dpsi, 3),
    iso_total = n_isoforms,
    iso_ptc_pos = n_ptc_pos,
    iso_sig = n_sig_iso,
    iso_nmd = n_nmd_iso
  )
}

summary_tbl <- bind_rows(summary_rows)
cat("Gene | rMATS sig FD | Max |dPSI| | Isoforms | PTC+ | Sig iso | NMD-resp\n")
cat("-----|-------------|-----------|----------|------|---------|--------\n")
for (i in seq_len(nrow(summary_tbl))) {
  r <- summary_tbl[i, ]
  cat(sprintf("%-7s| %11d | %9s | %8d | %4d | %7d | %7d\n",
              r$gene, r$rmats_sig_fd,
              if_else(is.na(r$rmats_max_dpsi), "   --", sprintf("%+.3f", r$rmats_max_dpsi)),
              r$iso_total, r$iso_ptc_pos, r$iso_sig, r$iso_nmd))
}

cat("\n  Done.\n")
