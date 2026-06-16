#!/usr/bin/env Rscript
# GSEA on March 2026 mashr gene-level results (all 6 cell types incl. DO_ALI)
# Ranking: signed mashr posterior-mean logFC
# Pathways: MSigDB Hallmark, KEGG, Reactome, GO-BP
# Match manuscript Methods: minSize 15, maxSize 500, FDR < 0.05 = significant
# Output: long-format TSV across cell types

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(data.table)
})

cat("fgsea:", as.character(packageVersion("fgsea")),
    "  msigdbr:", as.character(packageVersion("msigdbr")), "\n")

mashr_dir <- "/Users/petecastaldi/claude_projects/nmd/shortread_dge/mashr"
out_dir   <- "/Users/petecastaldi/claude_projects/nmd/tmp"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cts <- list(
  AT     = "nmd_mashr_dge_at_2026.3.10.csv",
  DD     = "nmd_mashr_dge_dd_2026.3.10.csv",
  DD_ALI = "nmd_mashr_dge_ddali_2026.3.10.csv",
  DO_ALI = "nmd_mashr_dge_doali_2026.3.10.csv",
  FB     = "nmd_mashr_dge_fb_2026.3.10.csv",
  MV     = "nmd_mashr_dge_mv_2026.3.10.csv"
)

# --- Pathway database (handle msigdbr API drift) ---
fetch <- function(coll, subcat = NULL) {
  args_new <- list(species = "Homo sapiens", collection = coll)
  args_old <- list(species = "Homo sapiens", category   = coll)
  if (!is.null(subcat)) {
    args_new$subcollection <- subcat
    args_old$subcategory   <- subcat
  }
  out <- tryCatch(do.call(msigdbr, args_new), error = function(e) NULL)
  if (is.null(out) || nrow(out) == 0)
    out <- tryCatch(do.call(msigdbr, args_old), error = function(e) NULL)
  out
}

cat("Loading MSigDB collections...\n")
parts <- list(
  H        = fetch("H"),
  KEGG     = fetch("C2", "CP:KEGG_LEGACY"),
  REACTOME = fetch("C2", "CP:REACTOME"),
  GOBP     = fetch("C5", "GO:BP")
)
# Fallback for KEGG naming variants across msigdbr versions
if (is.null(parts$KEGG) || nrow(parts$KEGG) == 0) parts$KEGG <- fetch("C2", "CP:KEGG")
if (is.null(parts$KEGG) || nrow(parts$KEGG) == 0) parts$KEGG <- fetch("C2", "CP:KEGG_MEDICUS")
parts <- parts[!vapply(parts, function(x) is.null(x) || nrow(x) == 0, logical(1))]

for (nm in names(parts))
  cat("  ", nm, ":", length(unique(parts[[nm]]$gs_name)), "pathways\n")

m_df <- rbindlist(parts, fill = TRUE)
sym_col <- if ("gene_symbol" %in% colnames(m_df)) "gene_symbol" else "human_gene_symbol"
pathways <- split(m_df[[sym_col]], m_df$gs_name)
pathways <- lapply(pathways, function(g) unique(g[!is.na(g) & g != ""]))
cat("Total pathways:", length(pathways), "\n\n")

# --- Per-cell-type GSEA ---
do_gsea <- function(ct, fname) {
  cat("=== ", ct, " ===\n", sep = "")
  d <- fread(file.path(mashr_dir, fname))
  d <- d[!is.na(hgnc_symbol) & hgnc_symbol != ""]
  d <- d[is.finite(logFC)]
  # Dedupe per symbol: keep strongest |logFC|
  d <- d[order(-abs(logFC))]
  d <- d[!duplicated(hgnc_symbol)]
  ranks <- setNames(d$logFC, d$hgnc_symbol)
  cat("  genes with rank:", length(ranks),
      "  range:", round(min(ranks), 3), "to", round(max(ranks), 3), "\n")

  set.seed(42)
  res <- fgsea(pathways, ranks, minSize = 15, maxSize = 500, eps = 0)
  res[, cell_type := ct]
  res
}

all_res <- rbindlist(lapply(names(cts), function(ct) do_gsea(ct, cts[[ct]])))

# Flatten leadingEdge list-column for TSV writing
all_res[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]

out <- file.path(out_dir, "gsea_mashr_gene_2026.3.10_run2026-05-18.tsv")
fwrite(all_res, out, sep = "\t")
cat("\nSaved:", out, "\n  rows:", nrow(all_res), "\n")

# --- Per-cell-type summary ---
cat("\n=== Significant (padj < 0.05) per cell type ===\n")
print(all_res[, .(
  tested = .N,
  sig    = sum(padj < 0.05, na.rm = TRUE),
  up     = sum(padj < 0.05 & NES > 0, na.rm = TRUE),
  down   = sum(padj < 0.05 & NES < 0, na.rm = TRUE)
), by = cell_type])

# --- Top 10 by padj per cell type (incl. DO_ALI focus) ---
cat("\n=== Top 10 pathways by padj per cell type ===\n")
for (ct in names(cts)) {
  cat("\n-- ", ct, " --\n", sep = "")
  top <- all_res[cell_type == ct][order(padj, -abs(NES))][1:10,
         .(pathway, padj = signif(padj, 3), NES = signif(NES, 3), size)]
  print(top, row.names = FALSE)
}

cat("\nDone.\n")
