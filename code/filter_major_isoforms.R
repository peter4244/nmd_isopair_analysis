#!/usr/bin/env Rscript
# Filter to Major Isoforms and Identify Dominant Isoform
# Part of v4.0 Expression-Based Event Detection Pipeline
#
# Purpose:
#   1. Subset to DMSO samples only (19 samples)
#   2. Calculate TMM normalization factors
#   3. Filter to major isoforms (≥5% of gene expression in ≥1 DMSO sample)
#      - Gene expression = sum of TMM-normalized CPM for all isoforms in gene
#      - Isoform proportion = isoform CPM / gene total CPM
#   4. Filter to genes with median total gene expression > 1 CPM across DMSO samples
#   5. Filter to genes with 2-10 major isoforms per gene
#   6. Identify dominant isoform per gene (highest median expression across DMSO samples)
#
# Inputs:
#   - rds/dge_isoform_nofilter_2026.2.7.rds (DGEList, 1.7M isoforms, 100% ENSG ID coverage)
#   - pheno/nmd_pheno_longreadbamids_2026.1.18.csv (sample metadata)
#
# Outputs:
#   - major_isoforms_dmso.rds: Major isoforms with proportions
#   - dominant_isoforms_dmso.rds: Dominant isoform per gene
#   - major_isoform_filtering_summary.txt: Filtering statistics

library(tidyverse)
library(edgeR)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   FILTER TO MAJOR ISOFORMS - v4.0                              ║\n")
cat("║   Expression-Based Filtering (≥5% threshold)                   ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Output directory
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v4.0_reference_based"
output_dir <- file.path(base_dir, "major_isoforms")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Load Data
# ============================================================================

cat("═══ Loading Data ═══\n\n")

# Load DGEList (new version with better ENSG ID mapping)
cat("Loading DGEList...\n")
dge <- readRDS("/Users/petecastaldi/claude_projects/nmd/rds/dge_isoform_nofilter_2026.2.7.rds")
cat("  Total isoforms:", nrow(dge), "\n")
cat("  Total samples:", ncol(dge), "\n\n")

# Load phenotype data
cat("Loading phenotype data...\n")
pheno <- read_csv("/Users/petecastaldi/claude_projects/nmd/pheno/nmd_pheno_longreadbamids_2026.1.18.csv",
                  show_col_types = FALSE)

# Add BAM basename as sample identifier
pheno <- pheno %>%
  mutate(
    bam_id = str_remove(basename(bam), "\\.bam$"),
    subject = id,
    cell_type = ct
  )

# Filter to DMSO samples
pheno_dmso <- pheno %>%
  filter(treatment == "DMSO")

cat("  DMSO samples:", nrow(pheno_dmso), "\n")
cat("  Cell types:", paste(sort(unique(pheno_dmso$cell_type)), collapse = ", "), "\n")
cat("\nSamples per cell type:\n")
print(table(pheno_dmso$cell_type))
cat("\n")

# Subset DGEList to DMSO samples
dge_dmso <- dge[, pheno_dmso$bam_id]
cat("DGEList subsetted to DMSO samples:", ncol(dge_dmso), "samples\n\n")

# Calculate TMM normalization factors
cat("Calculating TMM normalization factors...\n")
dge_dmso <- calcNormFactors(dge_dmso, method = "TMM")
cat("  Norm factors range:", round(min(dge_dmso$samples$norm.factors), 3), "to",
    round(max(dge_dmso$samples$norm.factors), 3), "\n")
cat("  Mean norm factor:", round(mean(dge_dmso$samples$norm.factors), 3), "\n")
cat("  ✓ TMM normalization applied\n\n")

# ============================================================================
# Extract Gene-Isoform Mapping
# ============================================================================

cat("═══ Building Gene-Isoform Mapping ═══\n\n")

# Extract isoform-gene mapping from DGEList
# The genes data frame contains annotation information
# Use gene_id_ens115_sqanti for ENSG IDs (100% coverage)
isoform_info <- dge_dmso$genes %>%
  as_tibble() %>%
  mutate(
    # Use ENSG ID from Ensembl 115 / SQANTI annotation
    gene_id = gene_id_ens115_sqanti
  ) %>%
  select(
    isoform_id = txid,
    gene_id,
    biotype = any_of("biotype")
  ) %>%
  filter(!is.na(gene_id), gene_id != "")

cat("Isoform-gene mapping:\n")
cat("  Total isoforms:", nrow(isoform_info), "\n")
cat("  Unique genes:", length(unique(isoform_info$gene_id)), "\n\n")

# ============================================================================
# Calculate Isoform Proportions
# ============================================================================

cat("═══ Calculating Isoform Proportions ═══\n\n")

# Calculate CPM (counts per million)
cat("Calculating CPM...\n")
cpm_dmso <- cpm(dge_dmso, log = FALSE)

# Convert to data frame with isoform IDs
cpm_df <- as.data.frame(cpm_dmso) %>%
  rownames_to_column("row_index") %>%
  mutate(row_index = as.integer(row_index)) %>%
  left_join(
    isoform_info %>% mutate(row_index = row_number()),
    by = "row_index"
  ) %>%
  select(-row_index) %>%
  filter(!is.na(isoform_id))

cat("  CPM calculated for", nrow(cpm_df), "isoforms\n\n")

# Calculate isoform proportions per sample
cat("Calculating isoform proportions per sample...\n")

# Get sample columns (exclude metadata columns)
sample_cols <- pheno_dmso$bam_id

# Calculate gene totals per sample
gene_totals <- cpm_df %>%
  select(gene_id, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "cpm"
  ) %>%
  group_by(gene_id, sample) %>%
  summarise(gene_total_cpm = sum(cpm, na.rm = TRUE), .groups = "drop")

# Calculate proportions
isoform_proportions <- cpm_df %>%
  select(isoform_id, gene_id, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "cpm"
  ) %>%
  left_join(gene_totals, by = c("gene_id", "sample")) %>%
  mutate(
    proportion = if_else(gene_total_cpm > 0, cpm / gene_total_cpm, 0)
  ) %>%
  select(isoform_id, gene_id, sample, cpm, gene_total_cpm, proportion)

cat("  Proportions calculated for", length(unique(isoform_proportions$isoform_id)), "isoforms\n")
cat("  Across", length(unique(isoform_proportions$sample)), "samples\n\n")

# ============================================================================
# Filter to Major Isoforms (≥5% in at least one sample)
# ============================================================================

cat("═══ Filtering to Major Isoforms (≥5% threshold) ═══\n\n")

# Identify major isoforms
major_isoforms_list <- isoform_proportions %>%
  group_by(isoform_id, gene_id) %>%
  summarise(
    max_proportion = max(proportion, na.rm = TRUE),
    mean_proportion = mean(proportion, na.rm = TRUE),
    median_proportion = median(proportion, na.rm = TRUE),
    n_samples_above_5pct = sum(proportion >= 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(max_proportion >= 0.05)  # Major isoform criterion

cat("Major isoforms identified:\n")
cat("  Starting isoforms:", length(unique(isoform_proportions$isoform_id)), "\n")
cat("  Major isoforms (≥5%):", nrow(major_isoforms_list), "\n")
cat("  Genes with major isoforms:", length(unique(major_isoforms_list$gene_id)), "\n\n")

# Add back all proportion data for major isoforms
major_isoforms <- isoform_proportions %>%
  semi_join(major_isoforms_list, by = c("isoform_id", "gene_id"))

cat("Major isoform proportions dataset:\n")
cat("  Rows (isoform × sample):", nrow(major_isoforms), "\n")
cat("  Unique isoforms:", length(unique(major_isoforms$isoform_id)), "\n\n")

# ============================================================================
# Filter by Median Gene Expression
# ============================================================================

cat("═══ Filtering by Median Gene Expression ═══\n\n")

# Calculate median total gene CPM across DMSO samples
gene_expression <- major_isoforms %>%
  group_by(gene_id, sample) %>%
  summarise(gene_total_cpm = first(gene_total_cpm), .groups = "drop") %>%
  group_by(gene_id) %>%
  summarise(
    median_gene_cpm = median(gene_total_cpm),
    mean_gene_cpm = mean(gene_total_cpm),
    .groups = "drop"
  )

# Filter to genes with median > 1 CPM
genes_expressed <- gene_expression %>%
  filter(median_gene_cpm > 1)

cat("Gene expression filtering (median > 1 CPM):\n")
cat("  Genes before filter:", nrow(gene_expression), "\n")
cat("  Genes with median > 1 CPM:", nrow(genes_expressed), "\n")
cat("  Genes removed:", nrow(gene_expression) - nrow(genes_expressed), "\n")
cat("  Retention rate:", sprintf("%.1f%%", 100 * nrow(genes_expressed) / nrow(gene_expression)), "\n\n")

# Filter major isoforms to only expressed genes
major_isoforms_list <- major_isoforms_list %>%
  semi_join(genes_expressed, by = "gene_id")

major_isoforms <- major_isoforms %>%
  semi_join(genes_expressed, by = "gene_id")

cat("After expression filtering:\n")
cat("  Major isoforms:", nrow(major_isoforms_list), "\n")
cat("  Genes:", length(unique(major_isoforms_list$gene_id)), "\n\n")

# ============================================================================
# Filter to Genes with 2-10 Major Isoforms
# ============================================================================

cat("═══ Filtering to Genes with 2-10 Major Isoforms ═══\n\n")

# Count major isoforms per gene
isoforms_per_gene <- major_isoforms_list %>%
  group_by(gene_id) %>%
  summarise(n_major_isoforms = n(), .groups = "drop")

cat("Major isoforms per gene distribution:\n")
print(table(isoforms_per_gene$n_major_isoforms))
cat("\n")

# Filter to genes with 2-10 major isoforms
genes_2_10_isoforms <- isoforms_per_gene %>%
  filter(n_major_isoforms >= 2 & n_major_isoforms <= 10)

cat("Genes with 2-10 major isoforms:\n")
for (i in 2:10) {
  count <- sum(genes_2_10_isoforms$n_major_isoforms == i)
  if (count > 0) {
    cat(sprintf("  %d isoforms: %d genes\n", i, count))
  }
}
cat("  Total:", nrow(genes_2_10_isoforms), "genes\n\n")

# Filter major isoforms to only these genes
major_isoforms_filtered <- major_isoforms %>%
  semi_join(genes_2_10_isoforms, by = "gene_id")

major_isoforms_list_filtered <- major_isoforms_list %>%
  semi_join(genes_2_10_isoforms, by = "gene_id")

cat("Filtered major isoforms:\n")
cat("  Isoforms:", nrow(major_isoforms_list_filtered), "\n")
cat("  Genes:", length(unique(major_isoforms_filtered$gene_id)), "\n\n")

# ============================================================================
# Identify Dominant Isoform per Gene
# ============================================================================

cat("═══ Identifying Dominant Isoform per Gene ═══\n\n")
cat("Definition: Isoform with highest MEDIAN proportion across all DMSO samples\n\n")

# Calculate median proportion per isoform across all samples
dominant_isoforms <- major_isoforms_filtered %>%
  group_by(gene_id, isoform_id) %>%
  summarise(
    median_proportion = median(proportion, na.rm = TRUE),
    mean_proportion = mean(proportion, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  group_by(gene_id) %>%
  arrange(desc(median_proportion), .by_group = TRUE) %>%
  slice(1) %>%  # Take isoform with highest median
  ungroup() %>%
  select(
    gene_id,
    dominant_isoform_id = isoform_id,
    dominant_median_proportion = median_proportion,
    dominant_mean_proportion = mean_proportion
  )

cat("Dominant isoforms identified:\n")
cat("  Genes with dominant isoform:", nrow(dominant_isoforms), "\n")
cat("  Coverage: 100% (one dominant per gene)\n\n")

# Statistics on dominance
cat("Dominant isoform characteristics:\n")
cat("  Mean median proportion:", round(mean(dominant_isoforms$dominant_median_proportion), 3), "\n")
cat("  Median median proportion:", round(median(dominant_isoforms$dominant_median_proportion), 3), "\n")
cat("  Genes with dominant >50%:", sum(dominant_isoforms$dominant_median_proportion > 0.5),
    sprintf("(%.1f%%)", 100 * mean(dominant_isoforms$dominant_median_proportion > 0.5)), "\n")
cat("  Genes with dominant >75%:", sum(dominant_isoforms$dominant_median_proportion > 0.75),
    sprintf("(%.1f%%)", 100 * mean(dominant_isoforms$dominant_median_proportion > 0.75)), "\n\n")

# ============================================================================
# Save Results
# ============================================================================

cat("═══ Saving Results ═══\n\n")

# Save major isoforms (long format with all samples)
saveRDS(major_isoforms_filtered,
        file.path(output_dir, "major_isoforms_dmso.rds"))
cat("Saved: major_isoforms_dmso.rds\n")
cat("  Contains: Isoform proportions per sample for all major isoforms\n")
cat("  Dimensions:", nrow(major_isoforms_filtered), "rows (isoform × sample)\n\n")

# Save dominant isoforms
saveRDS(dominant_isoforms,
        file.path(output_dir, "dominant_isoforms_dmso.rds"))
cat("Saved: dominant_isoforms_dmso.rds\n")
cat("  Contains: Dominant isoform per gene (highest median proportion)\n")
cat("  Dimensions:", nrow(dominant_isoforms), "genes\n\n")

# Save summary statistics
summary_file <- file.path(output_dir, "major_isoform_filtering_summary.txt")
sink(summary_file)
cat("MAJOR ISOFORM FILTERING SUMMARY - v4.0\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

cat("INPUT DATA:\n")
cat("  DGEList: results/rds/dge_isoform_2026.1.20.rds\n")
cat("  Total isoforms:", nrow(dge), "\n")
cat("  Total samples:", ncol(dge), "\n")
cat("  DMSO samples:", nrow(pheno_dmso), "\n")
cat("  Cell types:", paste(sort(unique(pheno_dmso$cell_type)), collapse = ", "), "\n\n")

cat("FILTERING CASCADE:\n\n")

cat("Step 1: Major Isoform Filtering (≥5% threshold)\n")
cat("  Starting isoforms:", length(unique(isoform_proportions$isoform_id)), "\n")
cat("  Major isoforms (≥5% in ≥1 sample):", nrow(major_isoforms_list), "\n")
cat("  Genes with major isoforms:", length(unique(major_isoforms_list$gene_id)), "\n\n")

cat("Step 2: Gene Filtering (2-10 major isoforms)\n")
cat("  Starting genes:", length(unique(major_isoforms_list$gene_id)), "\n")
cat("  Genes with 2-10 major isoforms:", nrow(genes_2_10_isoforms), "\n")
cat("    - 2 isoforms:", sum(genes_2_10_isoforms$n_major_isoforms == 2), "\n")
cat("    - 3 isoforms:", sum(genes_2_10_isoforms$n_major_isoforms == 3), "\n")
cat("    - 4 isoforms:", sum(genes_2_10_isoforms$n_major_isoforms == 4), "\n")
cat("  Final major isoforms:", nrow(major_isoforms_list_filtered), "\n\n")

cat("DOMINANT ISOFORM IDENTIFICATION:\n")
cat("  Method: Highest median proportion across all DMSO samples\n")
cat("  Dominant isoforms identified:", nrow(dominant_isoforms), "\n")
cat("  Mean dominance (median proportion):", round(mean(dominant_isoforms$dominant_median_proportion), 3), "\n")
cat("  Median dominance:", round(median(dominant_isoforms$dominant_median_proportion), 3), "\n")
cat("  Genes with clear dominant (>50%):", sum(dominant_isoforms$dominant_median_proportion > 0.5),
    sprintf("(%.1f%%)", 100 * mean(dominant_isoforms$dominant_median_proportion > 0.5)), "\n\n")

cat("OUTPUT FILES:\n")
cat("  major_isoforms_dmso.rds - Proportions per sample (long format)\n")
cat("  dominant_isoforms_dmso.rds - Dominant isoform per gene\n")
cat("  major_isoform_filtering_summary.txt - This file\n\n")

cat("NEXT STEP:\n")
cat("  Run: code/extract_exon_structures.R\n")
cat("  Purpose: Extract exon coordinates for major isoforms from GFF files\n")

sink()

cat("Saved: major_isoform_filtering_summary.txt\n\n")

# ============================================================================
# Summary
# ============================================================================

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   MAJOR ISOFORM FILTERING COMPLETE                             ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("RESULTS SUMMARY:\n")
cat("  Major isoforms identified:", nrow(major_isoforms_list_filtered), "\n")
cat("  Genes with 2-10 major isoforms:", nrow(genes_2_10_isoforms), "\n")
cat("  Dominant isoforms defined:", nrow(dominant_isoforms), "\n\n")

cat("OUTPUT LOCATION:\n")
cat(" ", output_dir, "\n\n")

cat("NEXT STEP:\n")
cat("  Rscript code/extract_exon_structures.R\n\n")
