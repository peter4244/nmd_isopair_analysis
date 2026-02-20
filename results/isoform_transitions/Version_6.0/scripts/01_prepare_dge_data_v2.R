#!/usr/bin/env Rscript
#
# Script: 01_prepare_dge_data_v2.R
# Version: 6.0
# Purpose: Memory-optimized version - Filter to major isoforms FIRST, then identify dominant
#
# Strategy: Calculate proportions per sample, filter to major isoforms (≥5% in any sample),
#           then identify dominant. Avoids massive join on full 1.7M isoform dataset.
#
# Input:
#   - /Users/petecastaldi/claude_projects/nmd/rds/dge_isoform_nofilter_2026.2.7.rds
#
# Output:
#   - data/expression_data.rds (FILTERED to major isoforms only)
#   - data/dominant_isoforms.rds
#   - data/sample_metadata.rds
#   - data/filtering_stats.rds (statistics on filtering process)
#

library(tidyverse)
library(edgeR)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 1.1: Prepare DGEList Data (Memory-Optimized v2)       ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Paths
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"
dge_file <- "/Users/petecastaldi/claude_projects/nmd/rds/dge_isoform_nofilter_2026.2.7.rds"

# ═══════════════════════════════════════════════════════════════════
# SECTION 1: Load DGEList
# ═══════════════════════════════════════════════════════════════════

cat("Loading DGEList...\n")
dge <- readRDS(dge_file)

cat(sprintf("  Loaded DGEList with %d isoforms and %d samples\n",
            nrow(dge), ncol(dge)))

# ═══════════════════════════════════════════════════════════════════
# SECTION 2: Extract Sample Metadata
# ═══════════════════════════════════════════════════════════════════

cat("Extracting sample metadata...\n")
sample_metadata <- dge$samples %>%
  as_tibble(rownames = "sample_id") %>%
  mutate(
    sample = if (!"sample" %in% names(.)) sample_id else sample,
    treatment = factor(treatment, levels = c("DMSO", "Smg1i")),
    ct = factor(ct)
  )

cat(sprintf("  Extracted metadata for %d samples\n", nrow(sample_metadata)))
cat(sprintf("  Cell types: %s\n", paste(unique(sample_metadata$ct), collapse = ", ")))
cat(sprintf("  Treatments: %s\n", paste(unique(sample_metadata$treatment), collapse = ", ")))

# ═══════════════════════════════════════════════════════════════════
# SECTION 3: Calculate CPM (All Isoforms)
# ═══════════════════════════════════════════════════════════════════

cat("\nCalculating CPM for all isoforms...\n")

# Calculate CPM using edgeR
cpm_matrix <- cpm(dge, log = FALSE)

# Create isoform-to-gene mapping from genes slot
isoform_gene_map <- dge$genes %>%
  as_tibble() %>%
  select(isoform_id = txid, gene_id = gene_id_ens115_sqanti)

cat(sprintf("  Calculated CPM for %d isoforms\n", nrow(cpm_matrix)))

# ═══════════════════════════════════════════════════════════════════
# SECTION 4: Filter to Major Isoforms (MEMORY-EFFICIENT APPROACH)
# ═══════════════════════════════════════════════════════════════════

cat("\nFiltering to major isoforms (≥5% in at least one sample)...\n")
cat("  This avoids memory-intensive operations on full dataset\n")

# Process per gene to calculate proportions and filter
# This is memory-efficient: process one gene at a time

major_isoform_indices <- integer(0)
filtering_stats <- list()

n_genes <- length(unique(isoform_gene_map$gene_id))
genes_processed <- 0
genes_with_major <- 0

pb <- txtProgressBar(min = 0, max = n_genes, style = 3)

for (gene in unique(isoform_gene_map$gene_id)) {
  genes_processed <- genes_processed + 1

  # Get isoform indices for this gene
  gene_isoform_idx <- which(isoform_gene_map$gene_id == gene)

  if (length(gene_isoform_idx) == 0) next

  # Get CPM for this gene's isoforms
  gene_cpm <- cpm_matrix[gene_isoform_idx, , drop = FALSE]

  # Calculate proportions per sample
  gene_totals <- colSums(gene_cpm)

  # Avoid division by zero
  gene_totals[gene_totals == 0] <- 1

  gene_proportions <- sweep(gene_cpm, 2, gene_totals, "/")

  # Check which isoforms are ≥5% in at least one sample
  max_prop_per_isoform <- apply(gene_proportions, 1, max)
  is_major <- max_prop_per_isoform >= 0.05

  # Keep major isoforms for this gene
  if (any(is_major)) {
    major_isoform_indices <- c(major_isoform_indices,
                                gene_isoform_idx[is_major])
    genes_with_major <- genes_with_major + 1
  }

  # Update progress every 100 genes
  if (genes_processed %% 100 == 0) {
    setTxtProgressBar(pb, genes_processed)
  }
}

close(pb)

cat(sprintf("\n  Filtering complete:\n"))
cat(sprintf("    Genes processed: %d\n", genes_processed))
cat(sprintf("    Genes with major isoforms (≥5%%): %d (%.1f%%)\n",
            genes_with_major, 100 * genes_with_major / genes_processed))
cat(sprintf("    Major isoforms kept: %d / %d (%.1f%%)\n",
            length(major_isoform_indices), nrow(cpm_matrix),
            100 * length(major_isoform_indices) / nrow(cpm_matrix)))

# ═══════════════════════════════════════════════════════════════════
# SECTION 5: Create Filtered Expression Data
# ═══════════════════════════════════════════════════════════════════

cat("\nCreating filtered expression dataset...\n")

# Subset to major isoforms only
cpm_filtered <- cpm_matrix[major_isoform_indices, ]
isoform_gene_map_filtered <- isoform_gene_map[major_isoform_indices, ]

# Add proper row names
rownames(cpm_filtered) <- isoform_gene_map_filtered$isoform_id

# Convert to tidy format
expression_data <- cpm_filtered %>%
  as_tibble(rownames = "isoform_id") %>%
  pivot_longer(
    cols = -isoform_id,
    names_to = "sample_id",
    values_to = "cpm"
  ) %>%
  left_join(isoform_gene_map_filtered, by = "isoform_id")

cat(sprintf("  Filtered expression data:\n"))
cat(sprintf("    Isoforms: %d\n", length(unique(expression_data$isoform_id))))
cat(sprintf("    Genes: %d\n", length(unique(expression_data$gene_id))))
cat(sprintf("    Samples: %d\n", length(unique(expression_data$sample_id))))

# ═══════════════════════════════════════════════════════════════════
# SECTION 6: Identify Dominant Isoforms
# ═══════════════════════════════════════════════════════════════════

cat("\nIdentifying dominant isoforms (from major isoforms)...\n")

# Now this is fast because we're only working with major isoforms
dominant_isoforms <- expression_data %>%
  group_by(gene_id, isoform_id) %>%
  summarise(mean_cpm = mean(cpm), .groups = "drop") %>%
  group_by(gene_id) %>%
  arrange(desc(mean_cpm)) %>%
  mutate(
    rank = row_number(),
    is_dominant = (rank == 1)
  ) %>%
  ungroup() %>%
  filter(is_dominant) %>%
  select(gene_id, dominant_isoform_id = isoform_id, mean_expression = mean_cpm)

cat(sprintf("  Identified dominant isoforms for %d genes\n", nrow(dominant_isoforms)))

# Calculate dominant proportions (now fast with filtered data)
dominant_proportions <- expression_data %>%
  group_by(gene_id, sample_id) %>%
  mutate(gene_total_cpm = sum(cpm)) %>%
  ungroup() %>%
  inner_join(
    dominant_isoforms %>% select(gene_id, dominant_isoform_id),
    by = "gene_id"
  ) %>%
  filter(isoform_id == dominant_isoform_id) %>%
  mutate(dominant_proportion = cpm / gene_total_cpm) %>%
  group_by(gene_id) %>%
  summarise(
    mean_dominant_proportion = mean(dominant_proportion, na.rm = TRUE),
    .groups = "drop"
  )

dominant_isoforms <- dominant_isoforms %>%
  left_join(dominant_proportions, by = "gene_id")

cat(sprintf("  Mean proportion of gene expression from dominant isoform: %.1f%%\n",
            mean(dominant_isoforms$mean_dominant_proportion, na.rm = TRUE) * 100))

# ═══════════════════════════════════════════════════════════════════
# SECTION 7: Filtering Statistics
# ═══════════════════════════════════════════════════════════════════

cat("\nCalculating filtering statistics...\n")

# Isoforms per gene (before and after filtering)
isoforms_per_gene_before <- isoform_gene_map %>%
  group_by(gene_id) %>%
  summarise(n_isoforms_total = n(), .groups = "drop")

isoforms_per_gene_after <- expression_data %>%
  group_by(gene_id) %>%
  summarise(n_isoforms_major = n_distinct(isoform_id), .groups = "drop")

filtering_stats_summary <- isoforms_per_gene_before %>%
  left_join(isoforms_per_gene_after, by = "gene_id") %>%
  mutate(
    n_isoforms_major = replace_na(n_isoforms_major, 0),
    pct_kept = 100 * n_isoforms_major / n_isoforms_total
  )

cat(sprintf("  Statistics saved for %d genes\n", nrow(filtering_stats_summary)))

# ═══════════════════════════════════════════════════════════════════
# SECTION 8: Save Outputs
# ═══════════════════════════════════════════════════════════════════

cat("\nSaving outputs...\n")

# Create data directory
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# Save expression data (FILTERED)
saveRDS(expression_data, file.path(data_dir, "expression_data.rds"))
cat(sprintf("  ✓ %s\n", "expression_data.rds"))

# Save dominant isoforms
saveRDS(dominant_isoforms, file.path(data_dir, "dominant_isoforms.rds"))
cat(sprintf("  ✓ %s\n", "dominant_isoforms.rds"))

# Save sample metadata
saveRDS(sample_metadata, file.path(data_dir, "sample_metadata.rds"))
cat(sprintf("  ✓ %s\n", "sample_metadata.rds"))

# Save filtering statistics
saveRDS(filtering_stats_summary, file.path(data_dir, "filtering_stats.rds"))
cat(sprintf("  ✓ %s\n", "filtering_stats.rds"))

cat("\n✓ Step 1.1 complete (v2 - memory-optimized)\n\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Total isoforms in DGEList: %d\n", nrow(cpm_matrix)))
cat(sprintf("Major isoforms (≥5%% in ≥1 sample): %d (%.1f%%)\n",
            length(unique(expression_data$isoform_id)),
            100 * length(unique(expression_data$isoform_id)) / nrow(cpm_matrix)))
cat(sprintf("Genes with major isoforms: %d\n", length(unique(expression_data$gene_id))))
cat(sprintf("Samples: %d (%d cell types, %d treatments)\n",
            nrow(sample_metadata),
            length(unique(sample_metadata$ct)),
            length(unique(sample_metadata$treatment))))
cat(sprintf("Dominant isoforms identified: %d\n", nrow(dominant_isoforms)))
cat(sprintf("Mean dominant proportion: %.1f%%\n",
            mean(dominant_isoforms$mean_dominant_proportion, na.rm = TRUE) * 100))
cat("═══════════════════════════════════════════════════════════════════\n\n")
