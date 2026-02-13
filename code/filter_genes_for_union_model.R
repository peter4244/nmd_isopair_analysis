#!/usr/bin/env Rscript
# Filter Genes for Union Exon Model Analysis
# Criteria: Gene expression ≥1 CPM AND >1 isoform

library(tidyverse)
library(edgeR)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   GENE FILTERING FOR UNION EXON MODEL ANALYSIS               ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Configuration
output_dir <- "results/isoform_transitions/v3.0_reference_based"

# ============================================================================
# Step 1: Load Data
# ============================================================================

cat("Step 1: Loading data...\n")

# Load DGE data
dge <- readRDS("results/rds/dge_isoform_2026.1.20.rds")
cat("  Loaded DGEList:\n")
cat("    Isoforms:", nrow(dge$counts), "\n")
cat("    Samples:", ncol(dge$counts), "\n")

# Filter to DMSO samples only
dmso_samples <- dge$samples %>% filter(treatment == "DMSO")
cat("    DMSO samples:", nrow(dmso_samples), "\n\n")

# Subset DGE to DMSO
dge_dmso <- dge[, rownames(dmso_samples)]
cat("  Subsetted to DMSO samples\n\n")

# Calculate CPM
cpm_data <- cpm(dge_dmso)
cat("  Calculated CPM values\n\n")

# ============================================================================
# Step 2: Map Isoforms to Genes
# ============================================================================

cat("Step 2: Mapping isoforms to genes...\n")

# Get gene information from DGE object
gene_isoform_map <- dge_dmso$genes %>%
  as_tibble() %>%
  mutate(
    # Use GENCODE gene ID if available, otherwise SQANTI
    gene_id = if_else(!is.na(hgnc_id.gc), hgnc_id.gc, hgnc_id.sq),
    isoform_id = txid
  ) %>%
  filter(!is.na(gene_id)) %>%  # Remove isoforms with no gene mapping
  select(gene_id, isoform_id) %>%
  distinct()

cat("  Mapped isoforms:", nrow(gene_isoform_map), "\n")
cat("  Unique genes:", n_distinct(gene_isoform_map$gene_id), "\n\n")

# ============================================================================
# Step 3: Calculate Gene-Level Expression
# ============================================================================

cat("Step 3: Calculating gene-level expression...\n")

# For each isoform, get its mean CPM across DMSO samples
# Use txid from genes annotation, not rownames (which are numeric indices)
isoform_expression <- tibble(
  isoform_id = dge_dmso$genes$txid,
  mean_cpm = rowMeans(cpm_data)
)

# Join with gene mapping
isoform_gene_expr <- isoform_expression %>%
  inner_join(gene_isoform_map, by = "isoform_id")

# Calculate gene-level expression (sum of isoform CPMs)
gene_expression <- isoform_gene_expr %>%
  group_by(gene_id) %>%
  summarise(
    gene_cpm = sum(mean_cpm),
    n_isoforms = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(gene_cpm))

cat("  Total genes with expression data:", nrow(gene_expression), "\n")
cat("  Gene CPM range:", sprintf("%.2f - %.2f",
                                  min(gene_expression$gene_cpm),
                                  max(gene_expression$gene_cpm)), "\n\n")

# ============================================================================
# Step 4: Filter by Gene Expression ≥1 CPM
# ============================================================================

cat("Step 4: Filtering by gene expression ≥1 CPM...\n")

genes_expr_filtered <- gene_expression %>%
  filter(gene_cpm >= 1)

cat("  Genes before filter:", nrow(gene_expression), "\n")
cat("  Genes after filter (≥1 CPM):", nrow(genes_expr_filtered), "\n")
cat("  Genes removed:", nrow(gene_expression) - nrow(genes_expr_filtered), "\n\n")

# ============================================================================
# Step 5: Filter by >1 Isoform and ≤10 Isoforms
# ============================================================================

cat("Step 5: Filtering by >1 isoform and ≤10 isoforms...\n")

genes_filtered_iso <- genes_expr_filtered %>%
  filter(n_isoforms > 1)

cat("  Genes before filter:", nrow(genes_expr_filtered), "\n")
cat("  Genes with >1 isoform:", nrow(genes_filtered_iso), "\n")

# Apply cap at 10 isoforms
genes_final <- genes_filtered_iso %>%
  filter(n_isoforms <= 10)

cat("  Genes after cap (≤10 isoforms):", nrow(genes_final), "\n")
cat("  Genes removed (>10 isoforms):", nrow(genes_filtered_iso) - nrow(genes_final), "\n")

cat("  Genes before filter:", nrow(genes_expr_filtered), "\n")
cat("  Genes after filter (>1 isoform):", nrow(genes_final), "\n")
cat("  Genes removed (single isoform):", nrow(genes_expr_filtered) - nrow(genes_final), "\n\n")

# Distribution of isoform counts
isoform_dist <- genes_final %>%
  count(n_isoforms, name = "n_genes") %>%
  arrange(n_isoforms)

cat("  Distribution of isoform counts:\n")
print(isoform_dist)
cat("\n")

# Summary statistics
cat("  Summary statistics:\n")
cat("    Mean isoforms per gene:", sprintf("%.1f", mean(genes_final$n_isoforms)), "\n")
cat("    Median isoforms per gene:", median(genes_final$n_isoforms), "\n")
cat("    Max isoforms per gene:", max(genes_final$n_isoforms), "\n")
cat("    Genes with 2 isoforms:", sum(genes_final$n_isoforms == 2), "\n")
cat("    Genes with 3-4 isoforms:", sum(genes_final$n_isoforms %in% 3:4), "\n")
cat("    Genes with 5+ isoforms:", sum(genes_final$n_isoforms >= 5), "\n\n")

# ============================================================================
# Step 6: Create Isoform List for These Genes
# ============================================================================

cat("Step 6: Creating isoform list for union model construction...\n")

# Get all isoforms for the filtered genes
isoforms_for_union <- isoform_gene_expr %>%
  filter(gene_id %in% genes_final$gene_id) %>%
  arrange(gene_id, desc(mean_cpm))

cat("  Total isoforms:", nrow(isoforms_for_union), "\n")
cat("  Unique genes:", n_distinct(isoforms_for_union$gene_id), "\n\n")

# ============================================================================
# Step 7: Generate All Pairwise Transitions
# ============================================================================

cat("Step 7: Generating pairwise transitions...\n")

# For each gene, create all pairwise combinations of isoforms
all_transitions <- isoforms_for_union %>%
  group_by(gene_id) %>%
  summarise(
    isoforms = list(isoform_id),
    n_isoforms = n(),
    .groups = "drop"
  ) %>%
  mutate(
    transitions = map(isoforms, function(isos) {
      if (length(isos) < 2) return(NULL)

      # Generate all pairs (A, B) where A != B
      pairs <- expand.grid(
        isoform_ref = isos,
        isoform_alt = isos,
        stringsAsFactors = FALSE
      ) %>%
        filter(isoform_ref != isoform_alt) %>%
        as_tibble()

      pairs
    })
  ) %>%
  select(gene_id, n_isoforms, transitions) %>%
  unnest(transitions)

cat("  Total transitions:", nrow(all_transitions), "\n")
cat("  Genes:", n_distinct(all_transitions$gene_id), "\n\n")

# Transitions per gene
trans_per_gene <- all_transitions %>%
  count(gene_id, name = "n_transitions") %>%
  left_join(genes_final %>% select(gene_id, n_isoforms), by = "gene_id")

cat("  Transitions by isoform count:\n")
summary_trans <- trans_per_gene %>%
  group_by(n_isoforms) %>%
  summarise(
    n_genes = n(),
    mean_transitions = mean(n_transitions),
    total_transitions = sum(n_transitions),
    .groups = "drop"
  ) %>%
  arrange(n_isoforms)

print(summary_trans)
cat("\n")

# ============================================================================
# Step 8: Save Results
# ============================================================================

cat("Step 8: Saving filtered gene and isoform lists...\n")

# Save gene list with expression info
saveRDS(genes_final, file.path(output_dir, "genes_for_union_model.rds"))
write_tsv(genes_final, file.path(output_dir, "genes_for_union_model.tsv"))
cat("  Saved: genes_for_union_model.rds\n")
cat("  Saved: genes_for_union_model.tsv\n\n")

# Save isoform list
saveRDS(isoforms_for_union, file.path(output_dir, "isoforms_for_union_model.rds"))
write_tsv(isoforms_for_union, file.path(output_dir, "isoforms_for_union_model.tsv"))
cat("  Saved: isoforms_for_union_model.rds\n")
cat("  Saved: isoforms_for_union_model.tsv\n\n")

# Save transition list (will be filtered later by uncategorizable genes)
saveRDS(all_transitions, file.path(output_dir, "transitions_for_union_model_unfiltered.rds"))
cat("  Saved: transitions_for_union_model_unfiltered.rds\n\n")

# ============================================================================
# Summary Report
# ============================================================================

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   GENE FILTERING COMPLETE                                     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("═══ FILTERING SUMMARY ═══\n\n")
cat(sprintf("  Total genes in dataset:           %6d\n", nrow(gene_expression)))
cat(sprintf("  Genes with ≥1 CPM:                %6d\n", nrow(genes_expr_filtered)))
cat(sprintf("  Genes with >1 isoform:            %6d ← FINAL\n", nrow(genes_final)))
cat("\n")
cat(sprintf("  Total isoforms for these genes:   %6d\n", nrow(isoforms_for_union)))
cat(sprintf("  Total pairwise transitions:       %6d\n", nrow(all_transitions)))
cat("\n")

cat("═══ ISOFORM DISTRIBUTION ═══\n\n")
cat("  Genes by isoform count:\n")
for (i in seq_len(min(10, nrow(isoform_dist)))) {
  cat(sprintf("    %2d isoforms: %5d genes\n",
              isoform_dist$n_isoforms[i],
              isoform_dist$n_genes[i]))
}
if (nrow(isoform_dist) > 10) {
  cat(sprintf("    ... and %d more categories\n", nrow(isoform_dist) - 10))
}
cat("\n")

cat("═══ NEXT STEP ═══\n\n")
cat("  Run union exon model construction on these", nrow(genes_final), "genes\n")
cat("  Expected output: ~", nrow(all_transitions), "transitions with event vectors\n\n")

cat("  Command: Rscript code/build_union_exon_model_full.R\n\n")
