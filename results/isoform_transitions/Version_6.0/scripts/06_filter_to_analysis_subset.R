#!/usr/bin/env Rscript
#
# Script: 06_filter_to_analysis_subset.R
# Version: 6.0
# Purpose: Filter to analysis-ready dataset matching DGE analysis criteria
#
# Filtering Strategy:
#   1. Apply filterByExpr (matches long-read DGE analysis)
#   2. Exclude fusion/chimeric genes (keep GENCODE + PacBio novel genes)
#
# Input:
#   - /Users/petecastaldi/claude_projects/nmd/rds/dge_isoform_nofilter_2026.2.7.rds
#   - data/*.rds (all unfiltered data from scripts 01-05)
#
# Output:
#   - data/*_filtered.rds (7 filtered datasets, no profiles yet)
#   - results/filtering_report_script06.txt
#

library(tidyverse)
library(edgeR)

# ==============================================================================
# Input Validation Helpers
# ==============================================================================

validate_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required input file not found: %s\nHave you run the prerequisite scripts?", path))
  }
}

validate_columns <- function(df, required_cols, name) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing required columns: %s", name, paste(missing, collapse = ", ")))
  }
}

validate_overlap <- function(ids_a, ids_b, name_a, name_b, min_pct = 10) {
  overlap <- length(intersect(ids_a, ids_b))
  pct <- 100 * overlap / length(ids_a)
  cat(sprintf("  Cross-check: %d/%d %s found in %s (%.1f%%)\n",
              overlap, length(ids_a), name_a, name_b, pct))
  if (pct < min_pct) {
    stop(sprintf("FATAL: Only %.1f%% of %s found in %s. Files may be from different pipeline runs.",
                 pct, name_a, name_b))
  }
}

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 6: Filter to Analysis-Ready Subset                    ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Paths
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"
dge_file <- "/Users/petecastaldi/claude_projects/nmd/rds/dge_isoform_nofilter_2026.2.7.rds"

# ==============================================================================
# 0. Validate Inputs
# ==============================================================================

cat("Validating inputs...\n")
validate_file(dge_file)
validate_file(file.path(base_dir, "data/expression_data.rds"))
validate_file(file.path(base_dir, "data/dominant_isoforms.rds"))
validate_file(file.path(base_dir, "data/isoform_structures.rds"))
validate_file(file.path(base_dir, "data/union_exons.rds"))
validate_file(file.path(base_dir, "data/isoform_union_mapping.rds"))
validate_file(file.path(base_dir, "data/isoform_cds_metadata.rds"))
validate_file(file.path(base_dir, "data/isoform_union_exons_annotated.rds"))
cat("  All required input files found.\n\n")

# ═══════════════════════════════════════════════════════════════════
# SECTION 1: Load Original DGEList and Apply filterByExpr
# ═══════════════════════════════════════════════════════════════════

cat("Loading original DGEList...\n")
dge <- readRDS(dge_file)

cat(sprintf("  Original: %s isoforms, %s genes\n",
            format(nrow(dge), big.mark=","),
            format(length(unique(dge$genes$gene_id_ens115_sqanti)), big.mark=",")))

cat("\nBuilding design matrix...\n")
design_all <- model.matrix(formula("~treatment+ct+ct*treatment"), data = dge$samples)
cat(sprintf("  Design: %d samples × %d coefficients\n", nrow(design_all), ncol(design_all)))

cat("\nApplying filterByExpr (matching DGE analysis)...\n")
cat("  Parameters: min.count=5, min.total.count=10, min.prop=0\n")
keep_expr <- filterByExpr(dge, design=design_all, min.count=5, min.total.count=10, min.prop=0)

cat(sprintf("  Isoforms passing filter: %s (%.1f%%)\n",
            format(sum(keep_expr), big.mark=","),
            100 * mean(keep_expr)))

# ═══════════════════════════════════════════════════════════════════
# SECTION 2: Exclude Fusion/Chimeric Genes
# ═══════════════════════════════════════════════════════════════════

cat("\nExcluding fusion/chimeric genes...\n")
cat("  KEEP: GENCODE (ENSG[0-9]+) + PacBio Novel (novelGene_*)\n")
cat("  EXCLUDE: Fusion/Chimeric (e.g., ENSG_ENSG)\n")

# Subset to isoforms passing expression filter
dge_expr <- dge[keep_expr, , keep.lib.sizes=FALSE]
gene_ids_expr <- dge_expr$genes$gene_id_ens115_sqanti

# Categorize genes
is_gencode <- grepl("^ENSG[0-9]+$", gene_ids_expr)
is_novel <- grepl("^novelGene_", gene_ids_expr)
is_fusion <- !is_gencode & !is_novel

cat(sprintf("\n  Gene categories after filterByExpr:\n"))
cat(sprintf("    GENCODE: %s genes\n", format(length(unique(gene_ids_expr[is_gencode])), big.mark=",")))
cat(sprintf("    PacBio Novel: %s genes\n", format(length(unique(gene_ids_expr[is_novel])), big.mark=",")))
cat(sprintf("    Fusion/Chimeric: %s genes (will be excluded)\n",
            format(length(unique(gene_ids_expr[is_fusion])), big.mark=",")))

# Keep only GENCODE + Novel
keep_genes <- is_gencode | is_novel
dge_filtered <- dge_expr[keep_genes, , keep.lib.sizes=FALSE]

cat(sprintf("\n  Final filtered dataset:\n"))
cat(sprintf("    Isoforms: %s\n", format(nrow(dge_filtered), big.mark=",")))
cat(sprintf("    Genes: %s\n", format(length(unique(dge_filtered$genes$gene_id_ens115_sqanti)), big.mark=",")))

# Create isoform and gene ID lists for filtering downstream data
isoform_ids_to_keep <- dge_filtered$genes$txid
gene_ids_to_keep <- unique(dge_filtered$genes$gene_id_ens115_sqanti)

cat(sprintf("\n  Created filter lists:\n"))
cat(sprintf("    %s isoform IDs\n", format(length(isoform_ids_to_keep), big.mark=",")))
cat(sprintf("    %s gene IDs\n", format(length(gene_ids_to_keep), big.mark=",")))

# ═══════════════════════════════════════════════════════════════════
# SECTION 3: Filter Downstream Data Structures
# ═══════════════════════════════════════════════════════════════════

cat("\n\nFiltering downstream data structures...\n")

# Cross-check: verify DGE isoform IDs overlap with pipeline expression data
cat("\n  Cross-file consistency checks:\n")
expr_check <- readRDS(file.path(base_dir, "data/expression_data.rds"))
validate_overlap(isoform_ids_to_keep, unique(expr_check$isoform_id),
                 "DGE filter isoform_ids", "expression_data")
rm(expr_check)

# ---- 3.1: expression_data ----
cat("\n  [1/7] Filtering expression_data.rds...\n")
expression_data <- readRDS(file.path(base_dir, "data/expression_data.rds"))
cat(sprintf("    Before: %s rows\n", format(nrow(expression_data), big.mark=",")))

expression_data_filtered <- expression_data %>%
  filter(isoform_id %in% isoform_ids_to_keep)

cat(sprintf("    After: %s rows\n", format(nrow(expression_data_filtered), big.mark=",")))
saveRDS(expression_data_filtered, file.path(base_dir, "data/expression_data_filtered.rds"))
cat("    ✓ Saved: data/expression_data_filtered.rds\n")

# ---- 3.2: dominant_isoforms ----
cat("\n  [2/7] Filtering dominant_isoforms.rds...\n")
dominant_isoforms <- readRDS(file.path(base_dir, "data/dominant_isoforms.rds"))
cat(sprintf("    Before: %s genes\n", format(nrow(dominant_isoforms), big.mark=",")))

dominant_isoforms_filtered <- dominant_isoforms %>%
  filter(gene_id %in% gene_ids_to_keep)

cat(sprintf("    After: %s genes\n", format(nrow(dominant_isoforms_filtered), big.mark=",")))
saveRDS(dominant_isoforms_filtered, file.path(base_dir, "data/dominant_isoforms_filtered.rds"))
cat("    ✓ Saved: data/dominant_isoforms_filtered.rds\n")

# ---- 3.3: isoform_structures ----
cat("\n  [3/7] Filtering isoform_structures.rds...\n")
isoform_structures <- readRDS(file.path(base_dir, "data/isoform_structures.rds"))
cat(sprintf("    Before: %s rows\n", format(nrow(isoform_structures), big.mark=",")))

isoform_structures_filtered <- isoform_structures %>%
  filter(isoform_id %in% isoform_ids_to_keep)

cat(sprintf("    After: %s rows\n", format(nrow(isoform_structures_filtered), big.mark=",")))
saveRDS(isoform_structures_filtered, file.path(base_dir, "data/isoform_structures_filtered.rds"))
cat("    ✓ Saved: data/isoform_structures_filtered.rds\n")

# ---- 3.4: union_exons ----
cat("\n  [4/7] Filtering union_exons.rds...\n")
union_exons <- readRDS(file.path(base_dir, "data/union_exons.rds"))
cat(sprintf("    Before: %s union exons from %s genes\n",
            format(nrow(union_exons), big.mark=","),
            format(length(unique(union_exons$gene_id)), big.mark=",")))

union_exons_filtered <- union_exons %>%
  filter(gene_id %in% gene_ids_to_keep)

cat(sprintf("    After: %s union exons from %s genes\n",
            format(nrow(union_exons_filtered), big.mark=","),
            format(length(unique(union_exons_filtered$gene_id)), big.mark=",")))
saveRDS(union_exons_filtered, file.path(base_dir, "data/union_exons_filtered.rds"))
cat("    ✓ Saved: data/union_exons_filtered.rds\n")

# ---- 3.5: isoform_union_mapping ----
cat("\n  [5/7] Filtering isoform_union_mapping.rds...\n")
isoform_union_mapping <- readRDS(file.path(base_dir, "data/isoform_union_mapping.rds"))
cat(sprintf("    Before: %s rows\n", format(nrow(isoform_union_mapping), big.mark=",")))

# Filter by both isoform_id and gene_id
isoform_union_mapping_filtered <- isoform_union_mapping %>%
  filter(isoform_id %in% isoform_ids_to_keep,
         gene_id %in% gene_ids_to_keep)

cat(sprintf("    After: %s rows\n", format(nrow(isoform_union_mapping_filtered), big.mark=",")))
saveRDS(isoform_union_mapping_filtered, file.path(base_dir, "data/isoform_union_mapping_filtered.rds"))
cat("    ✓ Saved: data/isoform_union_mapping_filtered.rds\n")

# ---- 3.6: isoform_cds_metadata ----
cat("\n  [6/7] Filtering isoform_cds_metadata.rds...\n")
isoform_cds_metadata <- readRDS(file.path(base_dir, "data/isoform_cds_metadata.rds"))
cat(sprintf("    Before: %s isoforms\n", format(nrow(isoform_cds_metadata), big.mark=",")))

isoform_cds_metadata_filtered <- isoform_cds_metadata %>%
  filter(isoform_id %in% isoform_ids_to_keep)

cat(sprintf("    After: %s isoforms\n", format(nrow(isoform_cds_metadata_filtered), big.mark=",")))
saveRDS(isoform_cds_metadata_filtered, file.path(base_dir, "data/isoform_cds_metadata_filtered.rds"))
cat("    ✓ Saved: data/isoform_cds_metadata_filtered.rds\n")

# ---- 3.7: isoform_union_exons_annotated ----
cat("\n  [7/7] Filtering isoform_union_exons_annotated.rds...\n")
isoform_union_exons_annotated <- readRDS(file.path(base_dir, "data/isoform_union_exons_annotated.rds"))
cat(sprintf("    Before: %s rows\n", format(nrow(isoform_union_exons_annotated), big.mark=",")))

isoform_union_exons_annotated_filtered <- isoform_union_exons_annotated %>%
  filter(isoform_id %in% isoform_ids_to_keep,
         gene_id %in% gene_ids_to_keep)

cat(sprintf("    After: %s rows\n", format(nrow(isoform_union_exons_annotated_filtered), big.mark=",")))
saveRDS(isoform_union_exons_annotated_filtered, file.path(base_dir, "data/isoform_union_exons_annotated_filtered.rds"))
cat("    ✓ Saved: data/isoform_union_exons_annotated_filtered.rds\n")

# ═══════════════════════════════════════════════════════════════════
# SECTION 4: Generate Filtering Report
# ═══════════════════════════════════════════════════════════════════

cat("\n\nGenerating filtering report...\n")

# Calculate statistics
original_stats <- list(
  isoforms = nrow(dge),
  genes = length(unique(dge$genes$gene_id_ens115_sqanti))
)

filtered_stats <- list(
  isoforms = nrow(dge_filtered),
  genes = length(gene_ids_to_keep)
)

# Gene category breakdown
gene_categories_filtered <- tibble(gene_id = gene_ids_to_keep) %>%
  mutate(
    category = case_when(
      grepl("^ENSG[0-9]+$", gene_id) ~ "GENCODE",
      grepl("^novelGene_", gene_id) ~ "PacBio Novel",
      TRUE ~ "Other"
    )
  ) %>%
  count(category)

# Write report
report_file <- file.path(base_dir, "results/filtering_report_script06.txt")
sink(report_file)

cat("═══════════════════════════════════════════════════════════════════\n")
cat("  FILTERING REPORT - Script 06\n")
cat("  Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("FILTERING STRATEGY:\n")
cat("  1. Apply filterByExpr (min.count=5, min.total.count=10, min.prop=0)\n")
cat("  2. Exclude fusion/chimeric genes (keep GENCODE + PacBio Novel)\n\n")

cat("SUMMARY STATISTICS:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-30s %15s %15s %10s\n", "Metric", "Original", "Filtered", "% Kept"))
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-30s %15s %15s %9.1f%%\n",
            "Isoforms",
            format(original_stats$isoforms, big.mark=","),
            format(filtered_stats$isoforms, big.mark=","),
            100 * filtered_stats$isoforms / original_stats$isoforms))
cat(sprintf("%-30s %15s %15s %9.1f%%\n",
            "Genes",
            format(original_stats$genes, big.mark=","),
            format(filtered_stats$genes, big.mark=","),
            100 * filtered_stats$genes / original_stats$genes))
cat("─────────────────────────────────────────────────────────────────\n\n")

cat("GENE CATEGORY BREAKDOWN (FILTERED):\n")
cat("─────────────────────────────────────────────────────────────────\n")
for (i in 1:nrow(gene_categories_filtered)) {
  cat(sprintf("  %-20s %15s genes\n",
              gene_categories_filtered$category[i],
              format(gene_categories_filtered$n[i], big.mark=",")))
}
cat("─────────────────────────────────────────────────────────────────\n\n")

cat("OUTPUT FILES CREATED:\n")
cat("  ✓ data/expression_data_filtered.rds\n")
cat("  ✓ data/dominant_isoforms_filtered.rds\n")
cat("  ✓ data/isoform_structures_filtered.rds\n")
cat("  ✓ data/union_exons_filtered.rds\n")
cat("  ✓ data/isoform_union_mapping_filtered.rds\n")
cat("  ✓ data/isoform_cds_metadata_filtered.rds\n")
cat("  ✓ data/isoform_union_exons_annotated_filtered.rds\n\n")

cat("NOTE: Splicing profiles will be created in Script 07\n\n")

cat("═══════════════════════════════════════════════════════════════════\n")

sink()

cat(sprintf("  ✓ Report saved: results/filtering_report_script06.txt\n"))

# ═══════════════════════════════════════════════════════════════════
# SECTION 5: Summary
# ═══════════════════════════════════════════════════════════════════

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Filtered dataset: %s isoforms from %s genes\n",
            format(filtered_stats$isoforms, big.mark=","),
            format(filtered_stats$genes, big.mark=",")))
cat(sprintf("Excluded: %s fusion/chimeric genes\n",
            format(sum(is_fusion), big.mark=",")))
cat("\nAll filtered data files saved with '_filtered.rds' suffix\n")
cat("Next: Script 07 will create splicing choice profiles from filtered data\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("✓ Step 6 complete\n\n")
