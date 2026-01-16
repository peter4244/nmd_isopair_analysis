#!/usr/bin/env Rscript
# Test script: GSEA Approach 2 for DD cell type with 100,000 permutations
# Purpose: Check if higher permutation count resolves p-value = 0.0001 issue
#
# IMPORTANT NOTE: Initial version used P.Value instead of adj.P.Val for selecting
# most significant isoform, which did not match main report. This has been corrected.

library(tidyverse)
library(fgsea)
library(msigdbr)

cat("==================================================\n")
cat("GSEA Test: DD Cell Type with 100,000 Permutations\n")
cat("==================================================\n\n")

# Configuration
LONGREAD_DIR <- "../longread_dge"
SHORTREAD_DIR <- "../shortread_dge"
DATE_STAMP <- "2026.1.15"
SIGNIFICANCE_THRESHOLD <- 0.05
GSEA_MIN_SIZE <- 15
GSEA_MAX_SIZE <- 500

# Test with 10k vs 100k permutations for comparison
NPERM_TEST <- c(10000, 100000)

# Load long-read DIE results for DD
cat("Loading DD isoform-level DIE results...\n")
dd_die <- read_csv(file.path(LONGREAD_DIR, "nmd_dge_dd_2026.1.2.csv"), show_col_types = FALSE)
cat("Loaded", nrow(dd_die), "isoforms\n\n")

# Load short-read DGE for gene symbols
cat("Loading DD gene-level DGE for symbols...\n")
dd_gene <- read_csv(file.path(SHORTREAD_DIR, "nmd_sr_dge_dd_2026.1.2.csv"), show_col_types = FALSE)

# Create gene symbol mapping
gene_symbol_map <- dd_gene %>%
  mutate(gene_id_no_version = str_replace(ensembl_gene_id_version, "\\.\\d+$", "")) %>%
  select(gene_id_no_version, hgnc_symbol) %>%
  distinct() %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "")

cat("Created gene symbol map with", nrow(gene_symbol_map), "genes\n\n")

# Add gene symbols to DIE results
dd_die <- dd_die %>%
  mutate(gene_id_no_version = str_replace(hgnc_id, "\\.\\d+$", "")) %>%
  left_join(gene_symbol_map, by = "gene_id_no_version") %>%
  # Match main report logic: use hgnc_symbol if available, otherwise fall back to hgnc_id
  mutate(gene_symbol = coalesce(hgnc_symbol, hgnc_id))

# Approach 2: Use most significant isoform per gene
cat("Preparing Approach 2 ranking (most significant isoform per gene)...\n")
cat("IMPORTANT: Using adj.P.Val to match main report implementation\n")
cat("IMPORTANT: Using coalesce(hgnc_symbol, hgnc_id) for gene_symbol like main report\n")
ranked_data <- dd_die %>%
  filter(!is.na(gene_symbol), !is.na(logFC), gene_symbol != "", !str_detect(gene_symbol, "^ENSG")) %>%
  group_by(gene_symbol) %>%
  slice_min(order_by = adj.P.Val, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(logFC))

ranked_vector <- setNames(ranked_data$logFC, ranked_data$gene_symbol)
cat("Ranked genes:", length(ranked_vector), "\n\n")

# Load MSigDB gene sets
cat("Loading MSigDB gene sets...\n")
all_genesets <- bind_rows(
  msigdbr(species = "Homo sapiens", collection = "H"),
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY"),
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
  msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
)

pathways <- split(all_genesets$gene_symbol, all_genesets$gs_name)

# Filter pathways by size
pathways_filtered <- pathways[sapply(pathways, length) >= GSEA_MIN_SIZE &
                               sapply(pathways, length) <= GSEA_MAX_SIZE]
cat("Using", length(pathways_filtered), "pathways (size", GSEA_MIN_SIZE, "-", GSEA_MAX_SIZE, ")\n\n")

# Run GSEA with different permutation counts
results_list <- list()

for (nperm in NPERM_TEST) {
  cat("=========================================\n")
  cat("Running GSEA with", format(nperm, big.mark = ","), "permutations...\n")
  cat("=========================================\n")

  start_time <- Sys.time()

  gsea_res <- fgsea(
    pathways = pathways_filtered,
    stats = ranked_vector,
    minSize = GSEA_MIN_SIZE,
    maxSize = GSEA_MAX_SIZE,
    nperm = nperm
  )

  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  cat("Completed in", round(elapsed, 1), "seconds\n\n")

  # Show top results
  top_10 <- gsea_res %>%
    arrange(pval) %>%
    head(10) %>%
    select(pathway, pval, padj, ES, NES, size) %>%
    as.data.frame()

  cat("Top 10 pathways by p-value:\n")
  print(top_10)
  cat("\n")

  # Count how many hit the p-value floor
  min_detectable_pval <- 1 / nperm
  n_at_floor <- sum(gsea_res$pval == min_detectable_pval)

  cat("P-values at floor (p =", min_detectable_pval, "):", n_at_floor, "\n")
  cat("P-values < 0.05:", sum(gsea_res$pval < 0.05), "\n")
  cat("Adjusted p-values < 0.05:", sum(gsea_res$padj < 0.05), "\n\n")

  results_list[[as.character(nperm)]] <- gsea_res
}

# Compare the top pathways between runs
cat("=========================================\n")
cat("Comparison of P-values (Top 20 pathways)\n")
cat("=========================================\n\n")

top_pathways <- results_list[["10000"]] %>%
  arrange(pval) %>%
  head(20) %>%
  pull(pathway)

comparison <- map_dfr(NPERM_TEST, function(nperm) {
  results_list[[as.character(nperm)]] %>%
    filter(pathway %in% top_pathways) %>%
    select(pathway, pval, padj) %>%
    mutate(nperm = nperm)
}) %>%
  pivot_wider(
    id_cols = pathway,
    names_from = nperm,
    values_from = c(pval, padj),
    names_glue = "{.value}_{nperm}"
  ) %>%
  arrange(pval_10000)

cat("Showing p-values for top 20 pathways from 10K run:\n\n")
print(comparison, n = 20)

cat("\n=========================================\n")
cat("Summary\n")
cat("=========================================\n\n")

cat("With 10,000 permutations:\n")
cat("- Minimum detectable p-value: 0.0001\n")
cat("- Pathways at p-value floor:", sum(results_list[["10000"]]$pval == 0.0001), "\n\n")

cat("With 100,000 permutations:\n")
cat("- Minimum detectable p-value: 0.00001\n")
cat("- Pathways at p-value floor:", sum(results_list[["100000"]]$pval == 0.00001), "\n\n")

cat("Recommendation: ")
if (sum(results_list[["100000"]]$pval == 0.00001) > 0) {
  cat("Even with 100K permutations, some pathways hit the floor.\n")
  cat("Consider 1M permutations for highest-significance pathways.\n")
} else {
  cat("100K permutations sufficient - no pathways at p-value floor.\n")
}
cat("\n")

cat("Test completed successfully!\n")
