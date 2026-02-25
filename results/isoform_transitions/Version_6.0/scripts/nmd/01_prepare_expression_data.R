#!/usr/bin/env Rscript
#
# Script: 01_prepare_dge_data_v2.R
# Version: 6.0
# Purpose: Memory-optimized version - Filter to major isoforms FIRST, then identify dominant
#
# Strategy: Calculate proportions per sample, filter to major isoforms (≥5% in any sample
#           within either DMSO or Smg1i condition), then identify dominant isoform per gene
#           based on DMSO mean CPM. Stores condition-specific proportions and expression.
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

# Pre-compute condition-specific sample column indices (for CPM matrix)
dmso_cols  <- which(sample_metadata$treatment == "DMSO")
smg1i_cols <- which(sample_metadata$treatment == "Smg1i")
cat(sprintf("  DMSO samples: %d, Smg1i samples: %d\n",
            length(dmso_cols), length(smg1i_cols)))

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

# Free DGEList — we have cpm_matrix and isoform_gene_map now
rm(dge)
gc()
cat("  (Freed DGEList from memory)\n")

# ═══════════════════════════════════════════════════════════════════
# SECTION 4: Filter to Major Isoforms (MEMORY-EFFICIENT APPROACH)
# ═══════════════════════════════════════════════════════════════════

cat("\nFiltering to major isoforms (≥5% in at least one sample)...\n")
cat("  This avoids memory-intensive operations on full dataset\n")

# Process per gene to calculate proportions and filter
# Condition-stratified: pass if ≥5% in any sample within EITHER condition

major_isoform_indices <- integer(0)

n_genes <- length(unique(isoform_gene_map$gene_id))
genes_processed <- 0
genes_with_major <- 0
n_pass_dmso_only  <- 0L
n_pass_smg1i_only <- 0L
n_pass_both       <- 0L

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

  # Condition-stratified 5% filter
  max_prop_dmso  <- apply(gene_proportions[, dmso_cols, drop = FALSE], 1, max)
  max_prop_smg1i <- apply(gene_proportions[, smg1i_cols, drop = FALSE], 1, max)
  pass_dmso  <- max_prop_dmso >= 0.05
  pass_smg1i <- max_prop_smg1i >= 0.05
  is_major   <- pass_dmso | pass_smg1i

  # Track which condition caused each isoform to pass
  n_pass_dmso_only  <- n_pass_dmso_only  + sum(pass_dmso & !pass_smg1i)
  n_pass_smg1i_only <- n_pass_smg1i_only + sum(!pass_dmso & pass_smg1i)
  n_pass_both       <- n_pass_both       + sum(pass_dmso & pass_smg1i)

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
cat(sprintf("    Passed in DMSO only:  %d\n", n_pass_dmso_only))
cat(sprintf("    Passed in Smg1i only: %d\n", n_pass_smg1i_only))
cat(sprintf("    Passed in both:       %d\n", n_pass_both))

# ═══════════════════════════════════════════════════════════════════
# SECTION 5: Create Filtered Expression Data
# ═══════════════════════════════════════════════════════════════════

cat("\nCreating filtered expression dataset...\n")

# Subset to major isoforms only
cpm_filtered <- cpm_matrix[major_isoform_indices, ]
isoform_gene_map_filtered <- isoform_gene_map[major_isoform_indices, ]

# Capture counts before freeing (needed for stats/summary later)
n_total_isoforms <- nrow(cpm_matrix)
isoforms_per_gene_before <- isoform_gene_map %>%
  group_by(gene_id) %>%
  summarise(n_isoforms_total = n(), .groups = "drop")

# Free full CPM matrix — we only need the filtered subset now
rm(cpm_matrix, isoform_gene_map)
gc()
cat("  (Freed full CPM matrix from memory)\n")

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

# Free filtered CPM matrix — tidy expression_data has everything now
rm(cpm_filtered, isoform_gene_map_filtered)
gc()
cat("  (Freed filtered CPM matrix from memory)\n")

# ═══════════════════════════════════════════════════════════════════
# SECTION 6: Identify Dominant Isoforms
# ═══════════════════════════════════════════════════════════════════

cat("\nIdentifying dominant isoforms (from major isoforms, DMSO-based)...\n")

# Join treatment info for condition-specific calculations
expr_with_trt <- expression_data %>%
  left_join(sample_metadata %>% select(sample_id, treatment),
            by = "sample_id")

# Dominant = highest mean CPM in DMSO samples
dominant_isoforms <- expr_with_trt %>%
  filter(treatment == "DMSO") %>%
  group_by(gene_id, isoform_id) %>%
  summarise(mean_expression_dmso = mean(cpm), .groups = "drop") %>%
  group_by(gene_id) %>%
  arrange(desc(mean_expression_dmso)) %>%
  mutate(rank = row_number(), is_dominant = (rank == 1)) %>%
  ungroup() %>%
  filter(is_dominant) %>%
  select(gene_id, dominant_isoform_id = isoform_id, mean_expression_dmso)

cat(sprintf("  Identified dominant isoforms for %d genes\n", nrow(dominant_isoforms)))

# Add overall and Smg1i mean expression
dom_expr_overall <- expression_data %>%
  semi_join(dominant_isoforms, by = c("gene_id", "isoform_id" = "dominant_isoform_id")) %>%
  group_by(gene_id, isoform_id) %>%
  summarise(mean_expression = mean(cpm), .groups = "drop") %>%
  select(gene_id, mean_expression)

dom_expr_smg1i <- expr_with_trt %>%
  filter(treatment == "Smg1i") %>%
  semi_join(dominant_isoforms, by = c("gene_id", "isoform_id" = "dominant_isoform_id")) %>%
  group_by(gene_id, isoform_id) %>%
  summarise(mean_expression_smg1i = mean(cpm), .groups = "drop") %>%
  select(gene_id, mean_expression_smg1i)

dominant_isoforms <- dominant_isoforms %>%
  left_join(dom_expr_overall, by = "gene_id") %>%
  left_join(dom_expr_smg1i, by = "gene_id")

# Calculate condition-specific dominant proportions
dominant_proportions <- expr_with_trt %>%
  group_by(gene_id, sample_id, treatment) %>%
  mutate(gene_total_cpm = sum(cpm)) %>%
  ungroup() %>%
  inner_join(
    dominant_isoforms %>% select(gene_id, dominant_isoform_id),
    by = "gene_id"
  ) %>%
  filter(isoform_id == dominant_isoform_id) %>%
  mutate(dominant_proportion = cpm / gene_total_cpm)

# Overall mean proportion
dom_prop_overall <- dominant_proportions %>%
  group_by(gene_id) %>%
  summarise(mean_dominant_proportion = mean(dominant_proportion, na.rm = TRUE),
            .groups = "drop")

# DMSO mean proportion
dom_prop_dmso <- dominant_proportions %>%
  filter(treatment == "DMSO") %>%
  group_by(gene_id) %>%
  summarise(mean_dominant_proportion_dmso = mean(dominant_proportion, na.rm = TRUE),
            .groups = "drop")

# Smg1i mean proportion
dom_prop_smg1i <- dominant_proportions %>%
  filter(treatment == "Smg1i") %>%
  group_by(gene_id) %>%
  summarise(mean_dominant_proportion_smg1i = mean(dominant_proportion, na.rm = TRUE),
            .groups = "drop")

dominant_isoforms <- dominant_isoforms %>%
  left_join(dom_prop_overall, by = "gene_id") %>%
  left_join(dom_prop_dmso, by = "gene_id") %>%
  left_join(dom_prop_smg1i, by = "gene_id")

rm(expr_with_trt, dominant_proportions, dom_expr_overall, dom_expr_smg1i,
   dom_prop_overall, dom_prop_dmso, dom_prop_smg1i)

cat(sprintf("  Mean dominant proportion (overall): %.1f%%\n",
            mean(dominant_isoforms$mean_dominant_proportion, na.rm = TRUE) * 100))
cat(sprintf("  Mean dominant proportion (DMSO):    %.1f%%\n",
            mean(dominant_isoforms$mean_dominant_proportion_dmso, na.rm = TRUE) * 100))
cat(sprintf("  Mean dominant proportion (Smg1i):   %.1f%%\n",
            mean(dominant_isoforms$mean_dominant_proportion_smg1i, na.rm = TRUE) * 100))

# ═══════════════════════════════════════════════════════════════════
# SECTION 7: Filtering Statistics
# ═══════════════════════════════════════════════════════════════════

cat("\nCalculating filtering statistics...\n")

# isoforms_per_gene_before was computed in Section 5 before freeing isoform_gene_map
isoforms_per_gene_after <- expression_data %>%
  group_by(gene_id) %>%
  summarise(n_isoforms_major = n_distinct(isoform_id), .groups = "drop")

filtering_stats_summary <- isoforms_per_gene_before %>%
  left_join(isoforms_per_gene_after, by = "gene_id") %>%
  mutate(
    n_isoforms_major = replace_na(n_isoforms_major, 0),
    pct_kept = 100 * n_isoforms_major / n_isoforms_total
  )

# Per-condition pass counts (gene-level aggregation not needed; stored as scalars)
condition_pass_counts <- list(
  dmso_only  = n_pass_dmso_only,
  smg1i_only = n_pass_smg1i_only,
  both       = n_pass_both
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

# Save filtering statistics (gene-level summary + condition pass counts)
saveRDS(list(per_gene = filtering_stats_summary,
             condition_pass_counts = condition_pass_counts),
        file.path(data_dir, "filtering_stats.rds"))
cat(sprintf("  ✓ %s\n", "filtering_stats.rds"))

cat("\n✓ Step 1.1 complete (v2 - memory-optimized)\n\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Total isoforms in DGEList: %d\n", n_total_isoforms))
cat(sprintf("Major isoforms (≥5%% in ≥1 sample): %d (%.1f%%)\n",
            length(unique(expression_data$isoform_id)),
            100 * length(unique(expression_data$isoform_id)) / n_total_isoforms))
cat(sprintf("Genes with major isoforms: %d\n", length(unique(expression_data$gene_id))))
cat(sprintf("Samples: %d (%d cell types, %d treatments)\n",
            nrow(sample_metadata),
            length(unique(sample_metadata$ct)),
            length(unique(sample_metadata$treatment))))
cat(sprintf("Dominant isoforms identified: %d (based on DMSO mean CPM)\n", nrow(dominant_isoforms)))
cat(sprintf("Mean dominant proportion (overall): %.1f%%\n",
            mean(dominant_isoforms$mean_dominant_proportion, na.rm = TRUE) * 100))
cat(sprintf("Mean dominant proportion (DMSO):    %.1f%%\n",
            mean(dominant_isoforms$mean_dominant_proportion_dmso, na.rm = TRUE) * 100))
cat(sprintf("Mean dominant proportion (Smg1i):   %.1f%%\n",
            mean(dominant_isoforms$mean_dominant_proportion_smg1i, na.rm = TRUE) * 100))
cat(sprintf("5%% filter: %d DMSO-only, %d Smg1i-only, %d both\n",
            n_pass_dmso_only, n_pass_smg1i_only, n_pass_both))
cat("═══════════════════════════════════════════════════════════════════\n\n")
