#!/usr/bin/env Rscript
#
# Script: 01_prepare_dge_data.R
# Version: 6.0
# Purpose: Load DGEList, extract sample metadata, calculate expression, identify dominant isoforms
#
# Input:
#   - /Users/petecastaldi/claude_projects/nmd/rds/dge_isoform_nofilter_2026.2.7.rds
#
# Output:
#   - data/expression_data.rds (isoform × sample CPM/TPM matrix)
#   - data/dominant_isoforms.rds (gene_id, dominant_isoform_id, mean_expression)
#   - data/sample_metadata.rds (sample info from DGEList $samples slot)
#

library(tidyverse)
library(edgeR)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 1.1: Prepare DGEList Data                              ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Paths
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"
dge_file <- "/Users/petecastaldi/claude_projects/nmd/rds/dge_isoform_nofilter_2026.2.7.rds"

# ═══════════════════════════════════════════════════════════════════
# SECTION 1: Load DGEList
# ═══════════════════════════════════════════════════════════════════

cat("Loading DGEList...\n")
dge <- readRDS(dge_file)

# Verify structure
cat(sprintf("  Loaded DGEList with %d isoforms and %d samples\n",
            nrow(dge), ncol(dge)))

# ═══════════════════════════════════════════════════════════════════
# SECTION 2: Extract Sample Metadata
# ═══════════════════════════════════════════════════════════════════

cat("Extracting sample metadata...\n")
sample_metadata <- dge$samples %>%
  as_tibble(rownames = "sample_id") %>%
  mutate(
    # Ensure key columns are present
    sample = if (!"sample" %in% names(.)) sample_id else sample,
    treatment = factor(treatment, levels = c("DMSO", "Smg1i")),
    ct = factor(ct)
  )

cat(sprintf("  Extracted metadata for %d samples\n", nrow(sample_metadata)))
cat(sprintf("  Cell types: %s\n", paste(unique(sample_metadata$ct), collapse = ", ")))
cat(sprintf("  Treatments: %s\n", paste(unique(sample_metadata$treatment), collapse = ", ")))

# ═══════════════════════════════════════════════════════════════════
# SECTION 3: Calculate Expression Levels
# ═══════════════════════════════════════════════════════════════════

cat("Calculating CPM...\n")

# Calculate CPM using edgeR (accounts for normalization factors)
cpm_matrix <- cpm(dge, log = FALSE)

# Create isoform-to-gene mapping from genes slot
isoform_gene_map <- dge$genes %>%
  as_tibble() %>%
  select(isoform_id = txid, gene_id = gene_id_ens115_sqanti)

# Add proper row names to CPM matrix
rownames(cpm_matrix) <- isoform_gene_map$isoform_id

# Convert to tidy format
expression_data <- cpm_matrix %>%
  as_tibble(rownames = "isoform_id") %>%
  pivot_longer(
    cols = -isoform_id,
    names_to = "sample_id",
    values_to = "cpm"
  ) %>%
  # Join with gene_id mapping
  left_join(isoform_gene_map, by = "isoform_id")

cat(sprintf("  Calculated CPM for %d isoforms\n", length(unique(expression_data$isoform_id))))
cat(sprintf("  Across %d samples\n", length(unique(expression_data$sample_id))))
cat(sprintf("  Covering %d genes\n", length(unique(expression_data$gene_id))))

# ═══════════════════════════════════════════════════════════════════
# SECTION 4: Identify Dominant Isoforms per Gene
# ═══════════════════════════════════════════════════════════════════

cat("Identifying dominant isoforms...\n")

# Calculate mean expression per isoform across all samples
dominant_isoforms <- expression_data %>%
  group_by(gene_id, isoform_id) %>%
  summarise(mean_cpm = mean(cpm), .groups = "drop") %>%
  # Rank isoforms within each gene
  group_by(gene_id) %>%
  arrange(desc(mean_cpm)) %>%
  mutate(
    rank = row_number(),
    is_dominant = (rank == 1)
  ) %>%
  ungroup() %>%
  # Keep only dominant isoforms
  filter(is_dominant) %>%
  select(gene_id, dominant_isoform_id = isoform_id, mean_expression = mean_cpm)

cat(sprintf("  Identified dominant isoforms for %d genes\n", nrow(dominant_isoforms)))

# Calculate proportion of expression captured by dominant isoform
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
# SECTION 5: Save Outputs
# ═══════════════════════════════════════════════════════════════════

cat("Saving outputs...\n")

# Create data directory if needed
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# Save expression data
saveRDS(expression_data, file.path(data_dir, "expression_data.rds"))
cat(sprintf("  Saved: %s\n", file.path(data_dir, "expression_data.rds")))

# Save dominant isoforms
saveRDS(dominant_isoforms, file.path(data_dir, "dominant_isoforms.rds"))
cat(sprintf("  Saved: %s\n", file.path(data_dir, "dominant_isoforms.rds")))

# Save sample metadata
saveRDS(sample_metadata, file.path(data_dir, "sample_metadata.rds"))
cat(sprintf("  Saved: %s\n", file.path(data_dir, "sample_metadata.rds")))

cat("\n✓ Step 1.1 complete\n\n")
cat(sprintf("Summary:\n"))
cat(sprintf("  - %d isoforms across %d genes\n",
            length(unique(expression_data$isoform_id)),
            length(unique(expression_data$gene_id))))
cat(sprintf("  - %d samples (%d cell types, %d treatments)\n",
            nrow(sample_metadata),
            length(unique(sample_metadata$ct)),
            length(unique(sample_metadata$treatment))))
cat(sprintf("  - %d dominant isoforms identified\n", nrow(dominant_isoforms)))
cat("\n")
