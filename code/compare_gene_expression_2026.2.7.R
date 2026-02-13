#!/usr/bin/env Rscript
# Compare gene-level expression quantifications between isoform-aggregated and gene-level data
# Created: 2026-02-07
# Purpose: Modular comparison supporting DGEList or matrix inputs for gene-level data

library(edgeR)
library(tidyverse)

# Main comparison function -------------------------------------------------------
compare_gene_expression <- function(dge_isoform,
                                     gene_data,
                                     gene_id_col = "ensembl_gene_id",
                                     output_prefix = "gene_expression_comparison") {

  cat("\n=== Gene-Level Expression Comparison ===\n\n")

  # Check gene-level data type and prepare DGEList
  cat("Preparing gene-level data...\n")
  if (is(gene_data, "DGEList")) {
    dge_gene <- gene_data
    cat("  Input: DGEList with", nrow(dge_gene), "genes\n")
  } else if (is.matrix(gene_data) || is.data.frame(gene_data)) {
    # Create DGEList from counts matrix using isoform sample info
    cat("  Input: Counts matrix with", nrow(gene_data), "genes\n")
    cat("  Creating DGEList using sample info from isoform data...\n")
    dge_gene <- DGEList(counts = as.matrix(gene_data),
                        samples = dge_isoform$samples)
  } else {
    stop("gene_data must be a DGEList or counts matrix")
  }

  # Inspect dimensions
  cat("\nInput dimensions:\n")
  cat("  Isoform-level:", nrow(dge_isoform$counts), "transcripts x",
      ncol(dge_isoform$counts), "samples\n")
  cat("  Gene-level:", nrow(dge_gene$counts), "genes x",
      ncol(dge_gene$counts), "samples\n")

  # Filter and aggregate isoforms to gene level ---------------------------------
  cat("\nFiltering and aggregating isoforms...\n")

  # Get Ensembl gene IDs from isoform data
  if (!gene_id_col %in% colnames(dge_isoform$genes)) {
    stop(paste("Column", gene_id_col, "not found in dge_isoform$genes"))
  }

  gene_ids <- dge_isoform$genes[[gene_id_col]]

  # Filter to transcripts with valid Ensembl gene IDs
  has_gene_id <- !is.na(gene_ids) & gene_ids != "" & gene_ids != "NA"
  cat("  Transcripts with Ensembl gene IDs:", sum(has_gene_id), "/",
      length(gene_ids), "\n")

  # Subset to transcripts with gene IDs
  isoform_counts_filtered <- dge_isoform$counts[has_gene_id, ]
  gene_ids_filtered <- gene_ids[has_gene_id]

  # Aggregate by summing transcript counts per gene
  gene_counts_from_isoforms <- rowsum(isoform_counts_filtered,
                                      group = gene_ids_filtered)
  cat("  Aggregated to", nrow(gene_counts_from_isoforms), "genes\n")

  # Match datasets by Ensembl gene ID --------------------------------------------
  cat("\nMatching datasets by Ensembl gene ID...\n")

  # Get gene IDs from gene-level data
  gene_ids_gene <- rownames(dge_gene$counts)

  # Find common genes
  common_genes <- intersect(rownames(gene_counts_from_isoforms), gene_ids_gene)
  cat("  Common genes:", length(common_genes), "\n")

  if (length(common_genes) == 0) {
    stop("No common genes found between datasets!")
  }

  # Find common samples
  common_samples <- intersect(colnames(gene_counts_from_isoforms),
                              colnames(dge_gene$counts))
  cat("  Common samples:", length(common_samples), "\n")

  if (length(common_samples) == 0) {
    stop("No common samples found between datasets!")
  }

  # Subset both datasets to common genes and samples
  counts_iso_matched <- gene_counts_from_isoforms[common_genes, common_samples]
  counts_gene_matched <- dge_gene$counts[common_genes, common_samples]

  # Create matched DGEList objects -----------------------------------------------
  cat("\nCreating matched DGEList objects...\n")

  # Get sample metadata (from isoform DGEList)
  samples_matched <- dge_isoform$samples[common_samples, ]

  dge_iso_matched <- DGEList(counts = counts_iso_matched,
                             samples = samples_matched)
  dge_gene_matched <- DGEList(counts = counts_gene_matched,
                              samples = samples_matched)

  # TMM normalize separately -----------------------------------------------------
  cat("\nApplying TMM normalization...\n")

  dge_iso_matched <- calcNormFactors(dge_iso_matched, method = "TMM")
  dge_gene_matched <- calcNormFactors(dge_gene_matched, method = "TMM")

  cat("  Isoform norm factors range: [",
      round(min(dge_iso_matched$samples$norm.factors), 3), ", ",
      round(max(dge_iso_matched$samples$norm.factors), 3), "]\n", sep = "")
  cat("  Gene norm factors range: [",
      round(min(dge_gene_matched$samples$norm.factors), 3), ", ",
      round(max(dge_gene_matched$samples$norm.factors), 3), "]\n", sep = "")

  # Calculate CPM
  cpm_iso <- cpm(dge_iso_matched, log = FALSE)
  cpm_gene <- cpm(dge_gene_matched, log = FALSE)

  # Overall correlations ---------------------------------------------------------
  cat("\n=== Overall Correlations ===\n")

  overall_cor_spearman <- cor(as.vector(cpm_iso), as.vector(cpm_gene),
                              method = "spearman")
  overall_cor_pearson <- cor(as.vector(cpm_iso), as.vector(cpm_gene),
                             method = "pearson")

  cat("  Spearman:", round(overall_cor_spearman, 4), "\n")
  cat("  Pearson:", round(overall_cor_pearson, 4), "\n")

  # Per-sample correlations ------------------------------------------------------
  cat("\n=== Per-Sample Correlations ===\n")

  sample_cors <- data.frame(
    sample = common_samples,
    spearman = sapply(common_samples, function(s) {
      cor(cpm_iso[, s], cpm_gene[, s], method = "spearman")
    }),
    pearson = sapply(common_samples, function(s) {
      cor(cpm_iso[, s], cpm_gene[, s], method = "pearson")
    })
  )

  # Add sample metadata
  sample_cors <- cbind(sample_cors, samples_matched[common_samples, ])

  cat("  Mean Spearman:", round(mean(sample_cors$spearman), 4), "\n")
  cat("  Range: [", round(min(sample_cors$spearman), 4), ", ",
      round(max(sample_cors$spearman), 4), "]\n", sep = "")

  # Per-cell-type correlations ---------------------------------------------------
  if ("ct" %in% colnames(samples_matched)) {
    cat("\n=== Per-Cell-Type Correlations ===\n")

    cell_types <- unique(samples_matched$ct)
    ct_cors <- lapply(cell_types, function(ct) {
      ct_samples <- common_samples[samples_matched$ct == ct]
      ct_cpm_iso <- as.vector(cpm_iso[, ct_samples])
      ct_cpm_gene <- as.vector(cpm_gene[, ct_samples])

      data.frame(
        cell_type = ct,
        n_samples = length(ct_samples),
        spearman = cor(ct_cpm_iso, ct_cpm_gene, method = "spearman"),
        pearson = cor(ct_cpm_iso, ct_cpm_gene, method = "pearson")
      )
    })
    ct_cors <- do.call(rbind, ct_cors)

    print(ct_cors)
  }

  # Gene-level statistics --------------------------------------------------------
  cat("\n=== Gene-Level Analysis ===\n")

  # Mean expression across all samples
  mean_cpm_iso <- rowMeans(cpm_iso)
  mean_cpm_gene <- rowMeans(cpm_gene)

  # Log2 fold difference
  fold_diff <- log2((mean_cpm_iso + 1) / (mean_cpm_gene + 1))

  # Expressed genes (mean CPM > 1)
  expressed_genes <- (mean_cpm_iso > 1) | (mean_cpm_gene > 1)
  cat("  Expressed genes (mean CPM > 1):", sum(expressed_genes), "\n")

  # Discrepant genes (>2-fold difference)
  discrepant_genes <- abs(fold_diff) > 1 & expressed_genes
  cat("  Discrepant genes (>2-fold diff):", sum(discrepant_genes), "\n")

  # Summary statistics
  cat("\n  Summary statistics (expressed genes):\n")
  cat("    Median absolute CPM difference:",
      round(median(abs(mean_cpm_iso - mean_cpm_gene)[expressed_genes]), 2), "\n")
  cat("    Mean absolute CPM difference:",
      round(mean(abs(mean_cpm_iso - mean_cpm_gene)[expressed_genes]), 2), "\n")
  cat("    Median |log2 fold difference|:",
      round(median(abs(fold_diff[expressed_genes])), 3), "\n")

  # Create comparison data frame
  comparison_df <- data.frame(
    gene = common_genes,
    mean_cpm_isoform = mean_cpm_iso,
    mean_cpm_gene = mean_cpm_gene,
    log2_fold_diff = fold_diff,
    abs_log2_fold_diff = abs(fold_diff),
    expressed = expressed_genes,
    discrepant = discrepant_genes,
    stringsAsFactors = FALSE
  )

  # Top discrepancies
  top_disc <- comparison_df %>%
    filter(expressed) %>%
    arrange(desc(abs_log2_fold_diff)) %>%
    head(20)

  cat("\n  Top 20 discrepant genes:\n")
  print(top_disc[, c("gene", "log2_fold_diff", "mean_cpm_isoform", "mean_cpm_gene")])

  # Return results ---------------------------------------------------------------
  results <- list(
    comparison_df = comparison_df,
    sample_correlations = sample_cors,
    cell_type_correlations = if(exists("ct_cors")) ct_cors else NULL,
    overall_correlations = c(spearman = overall_cor_spearman,
                            pearson = overall_cor_pearson),
    cpm_isoform = cpm_iso,
    cpm_gene = cpm_gene,
    dge_isoform_matched = dge_iso_matched,
    dge_gene_matched = dge_gene_matched
  )

  # Save results
  cat("\nSaving results...\n")
  saveRDS(results, file = paste0("results/rds/", output_prefix, "_2026.2.7.rds"))
  write.csv(comparison_df,
            file = paste0("results/", output_prefix, "_2026.2.7.csv"),
            row.names = FALSE)

  cat("  - results/rds/", output_prefix, "_2026.2.7.rds\n", sep = "")
  cat("  - results/", output_prefix, "_2026.2.7.csv\n", sep = "")

  return(results)
}

# Run Analysis 1: DGEList vs DGEList --------------------------------------------
cat("\n########################################\n")
cat("# Analysis 1: DGEList vs DGEList       #\n")
cat("########################################\n")

# Load data
dge_isoform <- readRDS("rds/dge_isoform_nofilter_2026.2.7.rds")
dge_gene <- readRDS("rds/dge_gene_nofilter_2026.1.3.rds")

# Check for gene ID column name
cat("\nAvailable columns in dge_isoform$genes:\n")
print(colnames(dge_isoform$genes))

# Determine gene ID column
gene_id_col <- if ("ensembl_gene_id" %in% colnames(dge_isoform$genes)) {
  "ensembl_gene_id"
} else if ("gene_id" %in% colnames(dge_isoform$genes)) {
  "gene_id"
} else if ("hgnc_id" %in% colnames(dge_isoform$genes)) {
  "hgnc_id"
} else {
  stop("Could not find Ensembl gene ID column")
}

cat("Using gene ID column:", gene_id_col, "\n")

# Run comparison
results1 <- compare_gene_expression(
  dge_isoform = dge_isoform,
  gene_data = dge_gene,
  gene_id_col = gene_id_col,
  output_prefix = "comparison_nofilter_dgelist"
)

cat("\n=== Analysis 1 Complete ===\n\n")

# Placeholder for Analysis 2 (counts matrix) -----------------------------------
cat("# To run Analysis 2 with counts matrix:\n")
cat("# counts_matrix <- readRDS('path/to/counts_matrix.rds')\n")
cat("# results2 <- compare_gene_expression(\n")
cat("#   dge_isoform = dge_isoform,\n")
cat("#   gene_data = counts_matrix,\n")
cat("#   gene_id_col = gene_id_col,\n")
cat("#   output_prefix = 'comparison_nofilter_lengthcorrected'\n")
cat("# )\n")
