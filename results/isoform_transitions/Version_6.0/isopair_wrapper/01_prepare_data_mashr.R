#!/usr/bin/env Rscript
#
# NMD Isoform Transitions — Data Preparation (mashr Classification)
# Version: 7.0-mashr
#
# Purpose: Load isocall count matrix, normalize, filter, classify isoforms as
#          NMD-sensitive or non-NMD using mashr differential isoform expression
#          results. All study-specific logic lives here.
#
# Input:
#   - Isocall count matrix + GTF
#   - 6 mashr DIE CSVs (one per cell type)
#
# Output (to output_dir/):
#   - expression_data.rds    — filtered CPM matrix (isoform x sample)
#   - sample_metadata.rds    — data.frame with sample_id, treatment, ct, donor
#   - gene_map.rds           — tibble with isoform_id, gene_id
#   - nmd_classification.rds — named list: per cell type + all_samples
#   - dmso_samples.rds       — named list of DMSO sample IDs per cell type
#   - smg1i_samples.rds      — named list of Smg1i sample IDs per cell type
#
# Usage:
#   Rscript 01_prepare_data_mashr.R [--output-dir DIR]

library(edgeR)

# ==============================================================================
# 0. Configuration
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
output_dir <- "data_mashr"
if ("--output-dir" %in% args) {
  idx <- which(args == "--output-dir")
  if (idx < length(args)) output_dir <- args[idx + 1]
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Input paths ---
count_matrix_file <- "/Users/petecastaldi/claude_projects/nmd/isocall/nmd_lungcells/results/call/nmd_isocall.count_matrix.txt"
gtf_file <- "/Users/petecastaldi/claude_projects/nmd/isocall/nmd_lungcells/results/call/nmd_isocall.isoforms.gtf.gz"

# mashr DE results
de_dir <- "/Users/petecastaldi/claude_projects/nmd/isocall_dge/mashr"
de_date <- "2026.3.10"
# mashr file naming: nmd_mashr_die_{suffix}_{date}.csv
# Paper scope (2026-04-29): 4 non-ALI cell types only — AT, DD, FB, MV.
# DD_ALI / DO_ALI / DO dropped from the paper.
ct_to_de_suffix <- c(AT = "at", DD = "dd", FB = "fb", MV = "mv")

# --- Thresholds ---
# mashr has pre-computed nmd_responsive column; non-NMD uses adj.P.Val.
# Threshold raised from 0.50 → 0.30 (2026-04-29) because the 4-CT mashr rerun
# tightened adj.P.Val distributions (the per-CT limma fits feeding mashr were
# rerun on the 4-CT subset). At 0.50 only ~7-12% of isoforms qualified as
# non-NMD vs ~71% historically; 0.30 yields ~52-60% non-NMD per CT, restoring
# a usable pair set for downstream gene-matching.
NON_NMD_ADJ_P <- 0.30
MIN_PROP      <- 0.05

cat("=== NMD Data Preparation — mashr Classification ===\n\n")

# ==============================================================================
# 1. Load count matrix and build DGEList
# ==============================================================================

cat("Loading count matrix...\n")
count_df <- read.csv(count_matrix_file, check.names = FALSE)
isoform_ids <- count_df$id
count_mat <- as.matrix(count_df[, -1])
rownames(count_mat) <- isoform_ids
rm(count_df); gc(verbose = FALSE)
cat(sprintf("  %d isoforms x %d samples\n", nrow(count_mat), ncol(count_mat)))

# Parse sample metadata from column names
# Format: Sample{N}_{CellType}_{DonorID}_{Treatment}
parse_sample <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  treatment <- parts[length(parts)]
  donor     <- parts[length(parts) - 1]
  ct_parts  <- parts[2:(length(parts) - 2)]
  ct        <- paste(ct_parts, collapse = "_")
  data.frame(sample_id = x, ct = ct, donor = donor, treatment = treatment,
             stringsAsFactors = FALSE)
}

sample_metadata <- do.call(rbind, lapply(colnames(count_mat), parse_sample))
cat(sprintf("  %d samples: %s\n", nrow(sample_metadata),
            paste(names(table(sample_metadata$ct)), collapse = ", ")))

# --- Sample exclusion ---
# DO donor 029T is a PCA outlier (DMSO sample dist=252 from DO centroid, >3x next).
# Remove both treatment partners to maintain balanced design.
exclude_samples <- c("Sample25_DO_029T_DMSO", "Sample26_DO_029T_Smg1i")
excluded <- colnames(count_mat) %in% exclude_samples
if (sum(excluded) != length(exclude_samples)) {
  stop("Expected to find ", length(exclude_samples), " samples to exclude but found ", sum(excluded))
}
count_mat <- count_mat[, !excluded]
sample_metadata <- sample_metadata[!sample_metadata$sample_id %in% exclude_samples, ]
cat(sprintf("  Excluded %d samples (DO 029T PCA outlier): %s\n",
            length(exclude_samples), paste(exclude_samples, collapse = ", ")))
cat(sprintf("  Remaining: %d samples\n", ncol(count_mat)))

# Build gene map from GTF
cat("Parsing gene map from GTF...\n")
gtf_gr <- rtracklayer::import(gtf_file)
tx_rows <- as.data.frame(gtf_gr[gtf_gr$type == "transcript", ])
gene_map <- data.frame(
  isoform_id = tx_rows$transcript_id,
  gene_id    = tx_rows$gene_id,
  stringsAsFactors = FALSE
)
gene_map <- gene_map[!duplicated(gene_map$isoform_id), ]
cat(sprintf("  %d isoforms mapped to %d genes\n",
            nrow(gene_map), length(unique(gene_map$gene_id))))

# DGEList + TMM normalization
dge <- DGEList(counts = count_mat,
               genes = data.frame(gene_id = gene_map$gene_id[
                 match(rownames(count_mat), gene_map$isoform_id)]))
dge <- calcNormFactors(dge, method = "TMM")
cpm_mat <- cpm(dge, log = FALSE)

# ==============================================================================
# 2. Condition-stratified 5% isoform filter
# ==============================================================================

cat("\nApplying 5% condition-stratified filter...\n")
dmso_cols  <- sample_metadata$sample_id[sample_metadata$treatment == "DMSO"]
smg1i_cols <- sample_metadata$sample_id[sample_metadata$treatment == "Smg1i"]

# Per-gene proportions in each condition
compute_max_prop <- function(cpm, samples, gene_ids) {
  cpm_sub <- cpm[, samples, drop = FALSE]
  # Gene totals per sample
  gene_totals <- rowsum(cpm_sub, gene_ids)
  # Isoform proportions: cpm / gene_total for its gene
  props <- cpm_sub / gene_totals[gene_ids, , drop = FALSE]
  props[is.nan(props)] <- 0
  # Max proportion across samples within this condition
  apply(props, 1, max, na.rm = TRUE)
}

gene_ids_all <- gene_map$gene_id[match(rownames(cpm_mat), gene_map$isoform_id)]
max_prop_dmso  <- compute_max_prop(cpm_mat, dmso_cols, gene_ids_all)
max_prop_smg1i <- compute_max_prop(cpm_mat, smg1i_cols, gene_ids_all)

keep_5pct <- max_prop_dmso >= MIN_PROP | max_prop_smg1i >= MIN_PROP
cat(sprintf("  5%% filter: %d / %d isoforms pass\n", sum(keep_5pct), length(keep_5pct)))

# filterByExpr with design (cell type + treatment)
design <- model.matrix(~ 0 + ct + treatment, data = sample_metadata)
keep_expr <- filterByExpr(dge[keep_5pct, ], design = design,
                          min.count = 2, min.total.count = 4)
keep_final <- names(keep_5pct)[keep_5pct]
keep_final <- keep_final[keep_expr]

# Exclude fusion genes (contain "--")
fusion_genes <- unique(gene_map$gene_id[grepl("--", gene_map$gene_id)])
fusion_isoforms <- gene_map$isoform_id[gene_map$gene_id %in% fusion_genes]
keep_final <- setdiff(keep_final, fusion_isoforms)

cpm_mat <- cpm_mat[keep_final, ]
cat(sprintf("  After filterByExpr + fusion exclusion: %d isoforms\n", nrow(cpm_mat)))

# ==============================================================================
# 3. NMD classification per cell type (mashr)
# ==============================================================================

cat("\nClassifying isoforms using mashr results...\n")
nmd_classification <- list()
dmso_samples_list  <- list()
smg1i_samples_list <- list()

for (ct in names(ct_to_de_suffix)) {
  de_file <- file.path(de_dir, sprintf("nmd_mashr_die_%s_%s.csv",
                                        ct_to_de_suffix[ct], de_date))
  de <- read.csv(de_file)

  # mashr uses 'txid' column (not 'transcript_id')
  # Filter to isoforms in our expression data
  de <- de[de$txid %in% rownames(cpm_mat), ]

  # mashr provides pre-computed nmd_responsive column
  nmd_ids     <- de$txid[de$nmd_responsive == TRUE]
  non_nmd_ids <- de$txid[de$adj.P.Val > NON_NMD_ADJ_P]

  nmd_classification[[ct]] <- list(nmd = nmd_ids, non_nmd = non_nmd_ids)

  ct_samples <- sample_metadata$sample_id[sample_metadata$ct == ct]
  dmso_samples_list[[ct]]  <- ct_samples[ct_samples %in% dmso_cols]
  smg1i_samples_list[[ct]] <- ct_samples[ct_samples %in% smg1i_cols]

  cat(sprintf("  %s: %d NMD, %d non-NMD\n", ct, length(nmd_ids), length(non_nmd_ids)))
}

# all_samples: NMD = union across AT, DD, FB, MV (paper scope: 4 non-ALI cell types)
# non-NMD = intersection across AT, DD, FB, MV
main_cts <- c("AT", "DD", "FB", "MV")
all_nmd <- unique(unlist(lapply(nmd_classification[main_cts], `[[`, "nmd")))
all_non_nmd <- Reduce(intersect, lapply(nmd_classification[main_cts], `[[`, "non_nmd"))
nmd_classification[["all_samples"]] <- list(nmd = all_nmd, non_nmd = all_non_nmd)
dmso_samples_list[["all_samples"]]  <- dmso_cols
smg1i_samples_list[["all_samples"]] <- smg1i_cols
cat(sprintf("  all_samples (AT, DD, FB, MV): %d NMD (union), %d non-NMD (intersection)\n",
            length(all_nmd), length(all_non_nmd)))

# ==============================================================================
# 4. Save outputs
# ==============================================================================

cat(sprintf("\nSaving to %s/...\n", output_dir))
saveRDS(cpm_mat,              file.path(output_dir, "expression_data.rds"))
saveRDS(sample_metadata,      file.path(output_dir, "sample_metadata.rds"))
saveRDS(gene_map,             file.path(output_dir, "gene_map.rds"))
saveRDS(nmd_classification,   file.path(output_dir, "nmd_classification.rds"))
saveRDS(dmso_samples_list,    file.path(output_dir, "dmso_samples.rds"))
saveRDS(smg1i_samples_list,   file.path(output_dir, "smg1i_samples.rds"))

cat("\nDone.\n")
