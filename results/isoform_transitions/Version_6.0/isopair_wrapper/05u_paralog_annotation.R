#!/usr/bin/env Rscript
# ==============================================================================
# 05u_paralog_annotation.R
#
# Query Ensembl for human gene paralogs via biomaRt.
#
# For predicting NMD from isoform structure, genes with high-similarity
# expressed paralogs pose a data leakage risk when one gene is in the training
# set (non-holdout chromosomes) and its paralog is in the holdout test set.
# Structurally similar isoforms could inflate test performance.
#
# We flag genes that meet ALL three criteria:
#   1. Paralog protein sequence identity >= 80%
#   2. Both the gene and its paralog are expressed in our dataset
#   3. The gene and its paralog are on DIFFERENT sides of the train/test split
#      (one on a holdout chromosome, one on a training chromosome)
#
# Only these cross-split high-similarity pairs represent actual leakage risk.
#
# Inputs:
#   - data_mashr/structures.rds (gene IDs + chromosomes)
#   - data_mashr/nmd_classification.rds (to identify expressed genes)
#
# Output:
#   - data_mashr/analysis_cache/paralog_genes.rds
#     $leakage_genes: versioned gene_ids to remove from test set
#     $leakage_pairs: the cross-split high-similarity pairs (for reporting)
#     $all_expressed_pairs: all expressed paralog pairs (before filtering)
#     $paralog_pairs: full biomaRt result with percent identity
#     $metadata: timestamp, Ensembl version, thresholds, counts
#
# Usage:
#   Rscript 05u_paralog_annotation.R
# ==============================================================================

suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

OUTPUT_PATH <- "data_mashr/analysis_cache/paralog_genes.rds"

# Thresholds
MIN_IDENTITY <- 80     # minimum % protein sequence identity
HOLDOUT_CHRS <- c("chr1", "chr3", "chr5", "chr7")

cat("=== Paralog Annotation ===\n")
cat("Started:", format(Sys.time()), "\n\n")
cat("Thresholds:\n")
cat("  Min protein identity:", MIN_IDENTITY, "%\n")
cat("  Holdout chromosomes:", paste(HOLDOUT_CHRS, collapse = ", "), "\n\n")

# ==============================================================================
# 1. Extract gene IDs and chromosome assignments
# ==============================================================================
cat("Loading gene IDs...\n")

structures <- readRDS("data_mashr/structures.rds")
nmd_class <- readRDS("data_mashr/nmd_classification.rds")

# All classified genes (NMD + non-NMD from all_samples)
classified_isoforms <- union(nmd_class$all_samples$nmd, nmd_class$all_samples$non_nmd)
classified_structures <- structures[structures$isoform_id %in% classified_isoforms, ]

# Gene → chromosome mapping (one gene can only be on one chromosome)
gene_chr <- classified_structures %>%
  distinct(gene_id, chr) %>%
  mutate(in_holdout = chr %in% HOLDOUT_CHRS)

classified_genes_versioned <- unique(gene_chr$gene_id)

# Strip version for biomaRt query
gene_id_map <- data.frame(
  gene_id_versioned = classified_genes_versioned,
  gene_id_unversioned = sub("\\.[0-9]+$", "", classified_genes_versioned),
  stringsAsFactors = FALSE
)

cat("  Classified genes:", length(classified_genes_versioned), "\n")
cat("  On holdout chromosomes:", sum(gene_chr$in_holdout), "\n")
cat("  On training chromosomes:", sum(!gene_chr$in_holdout), "\n")

# ==============================================================================
# 2. Query Ensembl for paralogs WITH percent identity
# ==============================================================================
cat("\nConnecting to Ensembl...\n")

ensembl <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
ens_version <- ensembl@host

cat("  Ensembl host:", ens_version, "\n")
cat("  Querying paralogs for", length(unique(gene_id_map$gene_id_unversioned)), "genes...\n")

# Query in batches to avoid timeout
batch_size <- 2000L
uniq_ids <- unique(gene_id_map$gene_id_unversioned)
n_batches <- ceiling(length(uniq_ids) / batch_size)

paralog_list <- vector("list", n_batches)
t0 <- proc.time()

for (b in seq_len(n_batches)) {
  start_idx <- (b - 1L) * batch_size + 1L
  end_idx <- min(b * batch_size, length(uniq_ids))
  batch_ids <- uniq_ids[start_idx:end_idx]

  # Retry up to 3 times on transient biomaRt failures
  for (attempt in 1:3) {
    paralog_list[[b]] <- tryCatch(
      getBM(
        attributes = c("ensembl_gene_id",
                        "hsapiens_paralog_ensembl_gene",
                        "hsapiens_paralog_perc_id"),
        filters = "ensembl_gene_id",
        values = batch_ids,
        mart = ensembl
      ),
      error = function(e) {
        if (attempt < 3) {
          cat(sprintf("  Batch %d attempt %d failed: %s. Retrying in 10s...\n",
                      b, attempt, conditionMessage(e)))
          Sys.sleep(10)
          NULL
        } else {
          stop(sprintf("biomaRt query failed after 3 attempts (batch %d): %s",
                       b, conditionMessage(e)))
        }
      }
    )
    if (!is.null(paralog_list[[b]])) break
  }

  cat(sprintf("  Batch %d/%d: %d genes, %d paralog rows (%.0f sec)\n",
              b, n_batches, length(batch_ids), nrow(paralog_list[[b]]),
              (proc.time() - t0)[3]))
}

paralog_pairs <- do.call(rbind, paralog_list)
stopifnot("biomaRt returned no paralog pairs" = nrow(paralog_pairs) > 0)

# Remove empty paralog entries
paralog_pairs <- paralog_pairs[paralog_pairs$hsapiens_paralog_ensembl_gene != "", ]

cat("\n  Total paralog pairs:", nrow(paralog_pairs), "\n")
cat("  Genes with >= 1 paralog:", length(unique(paralog_pairs$ensembl_gene_id)), "\n")

# Percent identity distribution
cat("\n  Percent identity distribution:\n")
print(summary(paralog_pairs$hsapiens_paralog_perc_id))

# ==============================================================================
# 3. Filter to high-similarity expressed cross-split pairs
# ==============================================================================
cat("\n=== Filtering paralog pairs ===\n")

expressed_unversioned <- unique(gene_id_map$gene_id_unversioned)

# Step 1: Both genes expressed in our dataset
expressed_pairs <- paralog_pairs %>%
  filter(
    ensembl_gene_id %in% expressed_unversioned,
    hsapiens_paralog_ensembl_gene %in% expressed_unversioned
  )
cat("  Step 1 — Both expressed:", nrow(expressed_pairs), "pairs,",
    length(unique(expressed_pairs$ensembl_gene_id)), "genes\n")

# Step 2: High protein sequence identity
high_sim_pairs <- expressed_pairs %>%
  filter(hsapiens_paralog_perc_id >= MIN_IDENTITY)
cat("  Step 2 — Identity >=", MIN_IDENTITY, "%:", nrow(high_sim_pairs), "pairs,",
    length(unique(high_sim_pairs$ensembl_gene_id)), "genes\n")

# Step 3: Cross-split (one in holdout, one in training)
# Map unversioned IDs to chromosomes via versioned → gene_chr
gene_chr_unver <- gene_chr %>%
  mutate(gene_id_unversioned = sub("\\.[0-9]+$", "", gene_id)) %>%
  distinct(gene_id_unversioned, in_holdout)

cross_split_pairs <- high_sim_pairs %>%
  inner_join(gene_chr_unver, by = c("ensembl_gene_id" = "gene_id_unversioned")) %>%
  rename(gene_in_holdout = in_holdout) %>%
  inner_join(gene_chr_unver, by = c("hsapiens_paralog_ensembl_gene" = "gene_id_unversioned")) %>%
  rename(paralog_in_holdout = in_holdout) %>%
  filter(gene_in_holdout != paralog_in_holdout)

cat("  Step 3 — Cross-split:", nrow(cross_split_pairs), "pairs,",
    length(unique(cross_split_pairs$ensembl_gene_id)), "genes\n")

# Genes to remove from test set: holdout-side genes in cross-split pairs
# (Remove the holdout gene, not the training gene)
leakage_genes_unversioned <- unique(c(
  cross_split_pairs$ensembl_gene_id[cross_split_pairs$gene_in_holdout],
  cross_split_pairs$hsapiens_paralog_ensembl_gene[cross_split_pairs$paralog_in_holdout]
))

# Map back to versioned IDs
leakage_genes_versioned <- gene_id_map$gene_id_versioned[
  gene_id_map$gene_id_unversioned %in% leakage_genes_unversioned
]
# Keep only those actually on holdout chromosomes
leakage_genes_versioned <- leakage_genes_versioned[
  leakage_genes_versioned %in% gene_chr$gene_id[gene_chr$in_holdout]
]

cat("\n  Leakage genes to remove from test set:", length(leakage_genes_versioned), "\n")

# ==============================================================================
# 4. Summary
# ==============================================================================
cat("\n=== SUMMARY ===\n")
cat("  Total classified genes:", length(classified_genes_versioned), "\n")
cat("  Genes with any paralog:", length(unique(paralog_pairs$ensembl_gene_id)), "\n")
cat("  High-similarity (>=", MIN_IDENTITY, "%) expressed cross-split pairs:",
    nrow(cross_split_pairs), "\n")
cat("  Holdout genes to remove (leakage risk):", length(leakage_genes_versioned),
    sprintf("(%.1f%% of holdout genes)\n",
            100 * length(leakage_genes_versioned) / sum(gene_chr$in_holdout)))

# Identity distribution of the leakage pairs
if (nrow(cross_split_pairs) > 0) {
  cat("\n  Identity distribution of leakage pairs:\n")
  print(summary(cross_split_pairs$hsapiens_paralog_perc_id))
}

# ==============================================================================
# 5. Save
# ==============================================================================
results <- list(
  # Primary output: genes to remove from test set
  leakage_genes = leakage_genes_versioned,

  # For reporting
  leakage_pairs = cross_split_pairs,
  all_expressed_pairs = expressed_pairs,
  paralog_pairs = paralog_pairs,

  metadata = list(
    run_timestamp = Sys.time(),
    ensembl_host = ens_version,
    min_identity = MIN_IDENTITY,
    holdout_chrs = HOLDOUT_CHRS,
    n_classified_genes = length(classified_genes_versioned),
    n_holdout_genes = sum(gene_chr$in_holdout),
    n_leakage_genes = length(leakage_genes_versioned),
    n_all_paralog_pairs = nrow(paralog_pairs),
    n_expressed_pairs = nrow(expressed_pairs),
    n_high_sim_pairs = nrow(high_sim_pairs),
    n_cross_split_pairs = nrow(cross_split_pairs)
  )
)

saveRDS(results, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
