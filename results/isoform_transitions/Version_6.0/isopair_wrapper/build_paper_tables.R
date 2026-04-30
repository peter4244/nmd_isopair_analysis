#!/usr/bin/env Rscript
# ==============================================================================
# build_paper_tables.R
#
# Produces paper-ready tables for the manuscript at non-NMD threshold = 0.30,
# 4-CT scope (AT, DD, FB, MV). Output: paper_outputs/manuscript_tables.md
# (pandoc-convertible to .docx).
#
# Tables produced:
#   1. SR gene-level mashr DGE — per-CT counts, % NMD-responsive, mean logFC
#   2. LR isoform-level mashr DIE — per-CT counts, % NMD-responsive, mean logFC
#   3. Pair set composition — per-CT C2/C4 counts (total, gene-matched, coding)
#   4. Top mashr meta-analysis NMD isoforms (by posterior logFC) — for Top 10
#      table reference. Not the supplement's full STable.
#
# Run from isopair_wrapper directory:
#   Rscript build_paper_tables.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(knitr)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
out_dir <- "paper_outputs"
dir.create(out_dir, showWarnings = FALSE)
out_md  <- file.path(out_dir, "manuscript_tables.md")

cell_types <- c("AT", "DD", "FB", "MV")
ct_to_label <- function(ct) tolower(gsub("_", "", ct))

cat("=== Building paper tables ===\n\n")

# ── Table 1: SR gene-level mashr DGE ───────────────────────────────────────────
cat("Building Table 1 (SR gene mashr DGE)...\n")
sr_dir <- "/Users/petecastaldi/claude_projects/nmd/shortread_dge/mashr"
sr_rows <- lapply(cell_types, function(ct) {
  f <- file.path(sr_dir, sprintf("nmd_mashr_dge_%s_2026.3.10.csv", tolower(ct)))
  if (!file.exists(f)) return(NULL)
  de <- read.csv(f, stringsAsFactors = FALSE)
  n_tested <- nrow(de)
  n_sig <- sum(de$adj.P.Val < 0.05, na.rm = TRUE)
  n_nmd <- sum(de$nmd_responsive, na.rm = TRUE)
  pct_nmd <- 100 * n_nmd / pmax(n_sig, 1)
  mean_logfc <- mean(de$logFC[de$nmd_responsive], na.rm = TRUE)
  data.frame(
    `Cell Type` = ct,
    `N Genes Tested` = format(n_tested, big.mark = ","),
    `N Significant` = format(n_sig, big.mark = ","),
    `N NMD-Responsive` = format(n_nmd, big.mark = ","),
    `% NMD-Responsive` = sprintf("%.1f%%", pct_nmd),
    `Mean logFC (NMD)` = sprintf("%.2f", mean_logfc),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})
sr_table <- do.call(rbind, sr_rows)
cat("  ", nrow(sr_table), "cell type rows\n")

# ── Table 2: LR isoform-level mashr DIE ────────────────────────────────────────
cat("Building Table 2 (LR isoform mashr DIE)...\n")
lr_dir <- "/Users/petecastaldi/claude_projects/nmd/isocall_dge/mashr"
lr_rows <- lapply(cell_types, function(ct) {
  f <- file.path(lr_dir, sprintf("nmd_mashr_die_%s_2026.3.10.csv", tolower(ct)))
  if (!file.exists(f)) return(NULL)
  de <- read.csv(f, stringsAsFactors = FALSE)
  n_tested <- nrow(de)
  n_sig <- sum(de$adj.P.Val < 0.05, na.rm = TRUE)
  n_nmd <- sum(de$nmd_responsive, na.rm = TRUE)
  pct_nmd <- 100 * n_nmd / pmax(n_sig, 1)
  mean_logfc <- mean(de$logFC[de$nmd_responsive], na.rm = TRUE)
  data.frame(
    `Cell Type` = ct,
    `N Isoforms Tested` = format(n_tested, big.mark = ","),
    `N Significant` = format(n_sig, big.mark = ","),
    `N NMD-Responsive` = format(n_nmd, big.mark = ","),
    `% NMD-Responsive` = sprintf("%.1f%%", pct_nmd),
    `Mean logFC (NMD)` = sprintf("%.2f", mean_logfc),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})
lr_table <- do.call(rbind, lr_rows)
cat("  ", nrow(lr_table), "cell type rows\n")

# ── Table 3: Pair set composition ──────────────────────────────────────────────
cat("Building Table 3 (Isopair pair set composition)...\n")
data_dir <- "data_mashr"
cds <- readRDS(file.path(data_dir, "cds.rds"))
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]

ps_rows <- lapply(c(cell_types, "all_samples"), function(ct) {
  ct_lab <- ct_to_label(ct)
  c2_file <- file.path(data_dir, sprintf("profiles_c2_%s.rds", ct_lab))
  c4_file <- file.path(data_dir, sprintf("profiles_c4_%s.rds", ct_lab))
  if (!file.exists(c2_file) || !file.exists(c4_file)) return(NULL)
  p_c2 <- readRDS(c2_file)
  p_c4 <- readRDS(c4_file)

  shared <- inner_join(
    unique(p_c2[, c("gene_id", "reference_isoform_id")]),
    unique(p_c4[, c("gene_id", "reference_isoform_id")]),
    by = c("gene_id", "reference_isoform_id"))
  gm_c2 <- semi_join(p_c2, shared, by = c("gene_id", "reference_isoform_id"))
  gm_c4 <- semi_join(p_c4, shared, by = c("gene_id", "reference_isoform_id"))

  cod_c2 <- gm_c2[gm_c2$reference_isoform_id %in% coding_ids &
                    gm_c2$comparator_isoform_id %in% coding_ids, ]
  cod_c4 <- gm_c4[gm_c4$reference_isoform_id %in% coding_ids &
                    gm_c4$comparator_isoform_id %in% coding_ids, ]

  data.frame(
    `Cell Type` = if (ct == "all_samples") "all_samples" else ct,
    `C2 (NMD) Pairs` = format(nrow(p_c2), big.mark = ","),
    `C4 (Control) Pairs` = format(nrow(p_c4), big.mark = ","),
    `C2 Gene-Matched` = format(nrow(gm_c2), big.mark = ","),
    `C4 Gene-Matched` = format(nrow(gm_c4), big.mark = ","),
    `C2 Coding (gene-matched)` = format(nrow(cod_c2), big.mark = ","),
    `C4 Coding (gene-matched)` = format(nrow(cod_c4), big.mark = ","),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})
ps_table <- do.call(rbind, ps_rows)
cat("  ", nrow(ps_table), "cell type rows\n")

# ── Table 4: meta-mashr summary ────────────────────────────────────────────────
cat("Building Table 4 (mashr meta-analysis summary)...\n")
mashr_lr_pm <- read.csv(file.path(lr_dir, "mashr_isoform_posterior_means_2026.3.10.csv"),
                         stringsAsFactors = FALSE)
mashr_lr_lfsr <- read.csv(file.path(lr_dir, "mashr_isoform_lfsr_2026.3.10.csv"),
                           stringsAsFactors = FALSE)
ct_cols <- c("Smg1i_in_AT", "Smg1i_in_DD", "Smg1i_in_FB", "Smg1i_in_MV")
mashr_lr_pm$mean_logFC <- rowMeans(mashr_lr_pm[, ct_cols, drop = FALSE], na.rm = TRUE)
mashr_lr_pm$min_lfsr <- apply(mashr_lr_lfsr[, ct_cols, drop = FALSE], 1, min, na.rm = TRUE)

mashr_summary <- data.frame(
  Metric = c(
    "Isoforms tested",
    "NMD-responsive in any CT (lfsr < 0.05, logFC > 0)",
    "NMD-responsive in all 4 CTs",
    "Median posterior logFC (NMD-responsive in any)",
    "Mean posterior logFC (NMD-responsive in any)"
  ),
  Value = c(
    format(nrow(mashr_lr_pm), big.mark = ","),
    format(sum(mashr_lr_pm$min_lfsr < 0.05 & mashr_lr_pm$mean_logFC > 0,
               na.rm = TRUE), big.mark = ","),
    format(sum(rowSums(mashr_lr_lfsr[, ct_cols] < 0.05) == 4 &
                 mashr_lr_pm$mean_logFC > 0, na.rm = TRUE), big.mark = ","),
    sprintf("%.3f", median(mashr_lr_pm$mean_logFC[mashr_lr_pm$min_lfsr < 0.05 &
                                                    mashr_lr_pm$mean_logFC > 0],
                            na.rm = TRUE)),
    sprintf("%.3f", mean(mashr_lr_pm$mean_logFC[mashr_lr_pm$min_lfsr < 0.05 &
                                                  mashr_lr_pm$mean_logFC > 0],
                          na.rm = TRUE))
  ),
  stringsAsFactors = FALSE
)

# ── Tables 5 & 6: Top 10 NMD-responsive genes and isoforms ────────────────────
cat("Building Tables 5 & 6 (Top 10 NMD-responsive genes / isoforms)...\n")

# Top 10 NMD-responsive genes (SR mashr): consistent hits across 4 CTs,
# ranked by min posterior logFC across CTs (most consistently strong)
sr_pm <- read.csv(file.path(sr_dir, "mashr_posterior_means_2026.3.10.csv"),
                   stringsAsFactors = FALSE)
sr_lfsr <- read.csv(file.path(sr_dir, "mashr_lfsr_2026.3.10.csv"),
                     stringsAsFactors = FALSE)
sr_ct_cols <- c("Smg1i_in_AT", "Smg1i_in_DD", "Smg1i_in_FB", "Smg1i_in_MV")
# Significant in all 4 CTs (lfsr<0.05) and positive logFC in all 4
sr_sig_all <- rowSums(sr_lfsr[, sr_ct_cols] < 0.05, na.rm = TRUE) == 4 &
              rowSums(sr_pm[, sr_ct_cols] > 0, na.rm = TRUE) == 4
sr_pm$min_logFC <- apply(sr_pm[, sr_ct_cols], 1, min, na.rm = TRUE)
sr_pm$mean_logFC <- rowMeans(sr_pm[, sr_ct_cols], na.rm = TRUE)
sr_top10_genes <- sr_pm[sr_sig_all, ] %>%
  arrange(desc(min_logFC)) %>%
  head(10) %>%
  transmute(
    `Gene ID` = gene_id,
    Symbol = ifelse(external_gene_name == "" | is.na(external_gene_name),
                    "—", external_gene_name),
    `logFC AT` = sprintf("%.2f", Smg1i_in_AT),
    `logFC DD` = sprintf("%.2f", Smg1i_in_DD),
    `logFC FB` = sprintf("%.2f", Smg1i_in_FB),
    `logFC MV` = sprintf("%.2f", Smg1i_in_MV),
    `Min logFC` = sprintf("%.2f", min_logFC),
    `Mean logFC` = sprintf("%.2f", mean_logFC)
  )
sr_top10_genes <- as.data.frame(sr_top10_genes, check.names = FALSE)

# Top 10 NMD-responsive isoforms (LR isoform mashr): consistent across 4 CTs
lr_pm <- mashr_lr_pm  # already loaded above
lr_lfsr <- mashr_lr_lfsr
lr_sig_all <- rowSums(lr_lfsr[, ct_cols] < 0.05, na.rm = TRUE) == 4 &
              rowSums(lr_pm[, ct_cols] > 0, na.rm = TRUE) == 4
lr_pm$min_logFC <- apply(lr_pm[, ct_cols], 1, min, na.rm = TRUE)
lr_top10_iso <- lr_pm[lr_sig_all, ] %>%
  arrange(desc(min_logFC)) %>%
  head(10) %>%
  transmute(
    `Transcript ID` = txid,
    Symbol = ifelse(hgnc_symbol == "" | is.na(hgnc_symbol), "—", hgnc_symbol),
    `Gene ID` = gene_id,
    `logFC AT` = sprintf("%.2f", Smg1i_in_AT),
    `logFC DD` = sprintf("%.2f", Smg1i_in_DD),
    `logFC FB` = sprintf("%.2f", Smg1i_in_FB),
    `logFC MV` = sprintf("%.2f", Smg1i_in_MV),
    `Min logFC` = sprintf("%.2f", min_logFC),
    `Mean logFC` = sprintf("%.2f", mean_logFC)
  )
lr_top10_iso <- as.data.frame(lr_top10_iso, check.names = FALSE)

# ── Write markdown output ──────────────────────────────────────────────────────
cat("Writing output to:", out_md, "\n")

con <- file(out_md, open = "w")
writeLines(c(
  "# NMD Manuscript Tables — 4-CT scope, non-NMD threshold = 0.30",
  "",
  sprintf("Generated: %s", format(Sys.time())),
  "",
  "Cell types: AT, DD, FB, MV (paper scope, non-ALI).",
  "Non-NMD threshold: adj.P.Val > 0.30.",
  "NMD-responsive: nmd_responsive == TRUE (mashr lfsr < 0.05 AND posterior logFC > 0).",
  "",
  "---",
  "",
  "## Table 1. Short-Read Gene-Level mashr DGE",
  "",
  "Differential gene expression at the gene level, mashr-shrunken across 4 cell types.",
  "",
  kable(sr_table, format = "markdown"),
  "",
  "**Column definitions:** Cell Type = primary lung cell type; N Genes Tested = genes passing expression filter; ",
  "N Significant = genes with adj.P.Val < 0.05; N NMD-Responsive = significant genes with mashr posterior ",
  "logFC > 0 (rescued by Smg1i); % NMD-Responsive = N NMD-Responsive / N Significant; Mean logFC (NMD) = ",
  "mean mashr posterior logFC across NMD-responsive genes.",
  "",
  "---",
  "",
  "## Table 2. Long-Read Isoform-Level mashr DIE (Main Text)",
  "",
  "Differential isoform expression at the isoform level, mashr-shrunken across 4 cell types.",
  "",
  kable(lr_table, format = "markdown"),
  "",
  "**Column definitions:** as above, applied at isoform level. N Isoforms Tested = isoforms passing the ",
  "5% expression filter and filterByExpr; mean logFC computed across mashr-NMD-responsive isoforms in each cell type.",
  "",
  "---",
  "",
  "## Table 3. Isopair Pair Set Composition",
  "",
  "Per-cell-type pair counts at each construction stage. The all_samples row aggregates across the 4 ",
  "non-ALI cell types using union NMD / intersection non-NMD classifications.",
  "",
  kable(ps_table, format = "markdown"),
  "",
  "**Column definitions:** C2 (NMD) Pairs = top non-NMD vs top NMD-responsive isoform per gene by CPM; ",
  "C4 (Control) Pairs = top two non-NMD isoforms per gene by CPM; Gene-Matched = pairs restricted to ",
  "(gene_id, reference_isoform_id) pairs shared between C2 and C4; Coding (gene-matched) = both ",
  "reference and comparator isoforms classified as coding by SQANTI3.",
  "",
  "---",
  "",
  "## Table 4. mashr Meta-Analysis Summary",
  "",
  "Mashr posterior summaries integrating across the 4 cell types.",
  "",
  kable(mashr_summary, format = "markdown"),
  "",
  "**Column definitions:** Metric = summary statistic; Value = computed value. NMD-responsive in any CT ",
  "uses min lfsr across CTs < 0.05 AND mean posterior logFC > 0. NMD-responsive in all 4 CTs uses lfsr < ",
  "0.05 in every CT plus mean posterior logFC > 0.",
  "",
  "---",
  "",
  "## Table 5. Top 10 Most Consistently NMD-Responsive Genes (SR mashr)",
  "",
  "Top genes ranked by minimum posterior logFC across the 4 cell types, restricted to genes that are ",
  "NMD-responsive (lfsr < 0.05 and posterior logFC > 0) in **all four** cell types.",
  "",
  kable(sr_top10_genes, format = "markdown"),
  "",
  "**Column definitions:** Gene ID = Ensembl gene ID with version; Symbol = HGNC gene symbol; ",
  "logFC {AT,DD,FB,MV} = mashr posterior logFC in each cell type (Smg1i vs DMSO); ",
  "Min logFC = minimum across the 4 cell types (the ranking metric — captures consistently strong response); ",
  "Mean logFC = mean across the 4 cell types.",
  "",
  "---",
  "",
  "## Table 6. Top 10 Most Consistently NMD-Responsive Isoforms (LR mashr)",
  "",
  "Top isoforms ranked by minimum posterior logFC across the 4 cell types, restricted to isoforms that ",
  "are NMD-responsive in all four cell types. Includes both reference (Ensembl) and novel (PacBio) isoforms.",
  "",
  kable(lr_top10_iso, format = "markdown"),
  "",
  "**Column definitions:** Transcript ID = Ensembl transcript ID (with .version) or PacBio novel ID ",
  "(`ENSG...novelN`); Symbol = parent gene's HGNC symbol; Gene ID = parent Ensembl gene ID; ",
  "logFC and Min/Mean as in Table 5.",
  "",
  "---",
  "",
  "## Notes on Drop-from-Paper Tables",
  "",
  "These tables previously appeared in supplement and are dropped from the paper:",
  "",
  "- **SR gene limma** (`shortread_dge/nmd_sr_dge_*_2026.1.2.csv`) — superseded by SR gene mashr (Table 1).",
  "- **LR isoform limma** (`isocall_dge/limma/yul/nmd_dge_*_2026.3.3.csv`) — superseded by LR isoform mashr (Table 2).",
  "- **DD_ALI / DO_ALI / DO results** — paper restricted to 4 non-ALI cell types (AT, DD, FB, MV); ALI conditions and DO are not used.",
  ""
), con)
close(con)

cat("\n=== Done ===\n")
cat("Tables written to:", out_md, "\n")
cat("To convert to .docx: pandoc", out_md, "-o", file.path(out_dir, "manuscript_tables.docx"), "\n")
