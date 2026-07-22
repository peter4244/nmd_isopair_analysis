#!/usr/bin/env Rscript
#
# NMD Isoform Transitions — Profile Building (Isopair Wrapper)
# Version: 7.0
#
# Purpose: Parse structures, build union exons, extract CDS, generate C2/C4
#          pairs, build splicing choice profiles using Isopair. Heavy compute;
#          uses checkpointing for robustness.
#
# Input:
#   - RDS files from 01_prepare_data.R
#   - Isocall GTF (isoform structures)
#   - SQANTI corrected CDS GFF3 (CDS annotations for known + novel isoforms)
#
# Output (to output_dir/):
#   - structures.rds, union_exons.rds, isoform_union_mapping.rds
#   - cds.rds, annotated_ue.rds, ptc.rds
#   - Per cell type: pairs_c2_{ct}.rds, pairs_c4_{ct}.rds,
#                     profiles_c2_{ct}.rds, profiles_c4_{ct}.rds
#
# Usage:
#   Rscript 02_build_profiles.R [--data-dir DIR] [--output-dir DIR]

library(Isopair)
library(dplyr)

# ==============================================================================
# 0. Configuration
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
data_dir   <- "data"
output_dir <- "data"
if ("--data-dir" %in% args) {
  idx <- which(args == "--data-dir")
  if (idx < length(args)) data_dir <- args[idx + 1]
}
if ("--output-dir" %in% args) {
  idx <- which(args == "--output-dir")
  if (idx < length(args)) output_dir <- args[idx + 1]
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# GTF paths
isocall_gtf <- "/Users/petecastaldi/claude_projects/nmd/isocall/nmd_lungcells/results/call/nmd_isocall.isoforms.gtf.gz"
sqanti_cds_gff3 <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.cds.gff3"

# Cell types to process (exclude DO — too few NMD pairs)
cell_types <- c("all_samples", "AT", "DD", "DD_ALI", "FB", "MV")

# Min pairs threshold
MIN_PAIRS <- 50

cat("=== NMD Profile Building (Isopair Wrapper) ===\n\n")

# ==============================================================================
# 1. Load prepared data
# ==============================================================================

cat("Loading prepared data...\n")
expr_mat  <- readRDS(file.path(data_dir, "expression_data.rds"))
gene_map  <- readRDS(file.path(data_dir, "gene_map.rds"))
nmd_class <- readRDS(file.path(data_dir, "nmd_classification.rds"))
dmso_samp <- readRDS(file.path(data_dir, "dmso_samples.rds"))
smg1i_samp <- readRDS(file.path(data_dir, "smg1i_samples.rds"))

cat(sprintf("  Expression: %d isoforms x %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))

# ==============================================================================
# 2. Isopair Infrastructure: structures, union exons, CDS, annotations
# ==============================================================================

# Parse isoform structures from isocall GTF
cat("\nParsing isoform structures...\n")
structures <- parseIsoformStructures(isocall_gtf,
  isoform_ids = rownames(expr_mat))
saveRDS(structures, file.path(output_dir, "structures.rds"))

# Build union exons
cat("Building union exons...\n")
ue <- buildUnionExons(structures)
saveRDS(ue$union_exons, file.path(output_dir, "union_exons.rds"))
saveRDS(ue$isoform_union_mapping, file.path(output_dir, "isoform_union_mapping.rds"))

# Extract CDS annotations from SQANTI corrected CDS GFF3
# This file has CDS features for both GENCODE and novel isoforms, using the
# same versioned IDs as isocall — no version stripping needed.
# Note: despite .gff3 extension, attributes are GTF format; specify format.
cat("Extracting CDS annotations...\n")
cat("  Loading SQANTI CDS GFF3...\n")
sqanti_gr <- rtracklayer::import(sqanti_cds_gff3, format = "gtf")
cds <- extractCdsAnnotations(sqanti_gr, isoform_ids = structures$isoform_id)

cat(sprintf("  CDS result: %d coding, %d unknown\n",
            sum(cds$coding_status == "coding"),
            sum(cds$coding_status == "unknown")))
saveRDS(cds, file.path(output_dir, "cds.rds"))

# Annotate region types
cat("Annotating region types...\n")
annotated_ue <- annotateRegionTypes(ue$isoform_union_mapping, cds)
saveRDS(annotated_ue, file.path(output_dir, "annotated_ue.rds"))

# Compute PTC status
cat("Computing PTC status...\n")
ptc <- computePtcStatus(structures, cds)
saveRDS(ptc, file.path(output_dir, "ptc.rds"))

# Gene map as tibble for Isopair
gene_map_tbl <- tibble::tibble(
  isoform_id = gene_map$isoform_id,
  gene_id = gene_map$gene_id
)
# Restrict to isoforms in structures
gene_map_tbl <- gene_map_tbl[gene_map_tbl$isoform_id %in% structures$isoform_id, ]

# ==============================================================================
# 3. Per cell type: Generate pairs and build profiles
# ==============================================================================

cat("\n--- Generating pairs and building profiles ---\n")

for (ct in cell_types) {
  cat(sprintf("\n=== %s ===\n", ct))

  nmd_ids     <- nmd_class[[ct]]$nmd
  non_nmd_ids <- nmd_class[[ct]]$non_nmd

  # Filter to isoforms present in structures
  nmd_ids     <- intersect(nmd_ids, structures$isoform_id)
  non_nmd_ids <- intersect(non_nmd_ids, structures$isoform_id)

  cat(sprintf("  NMD isoforms: %d, non-NMD: %d\n",
              length(nmd_ids), length(non_nmd_ids)))

  if (length(non_nmd_ids) < 2) {
    cat("  Skipping: too few non-NMD isoforms.\n")
    next
  }

  # --- C4: Top two non-NMD by CPM in DMSO ---
  expr_non_nmd <- expr_mat[intersect(non_nmd_ids, rownames(expr_mat)),
                            dmso_samp[[ct]], drop = FALSE]
  gene_map_non_nmd <- gene_map_tbl[gene_map_tbl$isoform_id %in% rownames(expr_non_nmd), ]

  pairs_c4 <- tryCatch(
    generatePairsExpression(expr_non_nmd, gene_map_non_nmd,
                             dmso_samp[[ct]], method = "top_two"),
    error = function(e) {
      cat(sprintf("  C4 pair generation failed: %s\n", e$message))
      NULL
    }
  )

  if (is.null(pairs_c4) || nrow(pairs_c4) < MIN_PAIRS) {
    cat(sprintf("  Skipping C4: %d pairs (< %d minimum)\n",
                if (is.null(pairs_c4)) 0L else nrow(pairs_c4), MIN_PAIRS))
    pairs_c4 <- NULL
  }

  # --- C2: Same reference as C4, paired with top NMD by CPM in Smg1i ---
  pairs_c2 <- NULL
  if (!is.null(pairs_c4) && length(nmd_ids) > 0) {
    expr_nmd <- expr_mat[intersect(nmd_ids, rownames(expr_mat)),
                          smg1i_samp[[ct]], drop = FALSE]
    gene_map_nmd <- gene_map_tbl[gene_map_tbl$isoform_id %in% rownames(expr_nmd), ]

    if (nrow(expr_nmd) > 0) {
      top_nmd <- identifyDominantIsoforms(
        expr_nmd, gene_map_nmd, smg1i_samp[[ct]], threshold = 0)

      pairs_c2 <- tryCatch({
        inner_join(
          pairs_c4[, c("gene_id", "reference_isoform_id")] |> distinct(),
          rename(top_nmd, comparator_isoform_id = dominant_isoform_id),
          by = "gene_id"
        )
      }, error = function(e) {
        cat(sprintf("  C2 pair generation failed: %s\n", e$message))
        NULL
      })
    }

    if (!is.null(pairs_c2) && nrow(pairs_c2) < MIN_PAIRS) {
      cat(sprintf("  Skipping C2: %d pairs (< %d minimum)\n",
                  nrow(pairs_c2), MIN_PAIRS))
      pairs_c2 <- NULL
    }
  }

  # --- Build profiles ---
  ct_label <- tolower(gsub("_", "", ct))
  ckpt_base <- file.path(output_dir, "checkpoints")

  if (!is.null(pairs_c4)) {
    cat(sprintf("  Building C4 profiles (%d pairs)...\n", nrow(pairs_c4)))
    saveRDS(pairs_c4, file.path(output_dir, sprintf("pairs_c4_%s.rds", ct_label)))

    profiles_c4 <- buildProfiles(
      pairs_c4, structures, ue$union_exons, ue$isoform_union_mapping,
      verify = TRUE,
      checkpoint_dir = file.path(ckpt_base, sprintf("c4_%s", ct_label)),
      checkpoint_interval = 500L
    )
    saveRDS(profiles_c4, file.path(output_dir, sprintf("profiles_c4_%s.rds", ct_label)))
    cat(sprintf("  C4: %d profiles, %d genes\n",
                nrow(profiles_c4), n_distinct(profiles_c4$gene_id)))
  }

  if (!is.null(pairs_c2)) {
    cat(sprintf("  Building C2 profiles (%d pairs)...\n", nrow(pairs_c2)))
    saveRDS(pairs_c2, file.path(output_dir, sprintf("pairs_c2_%s.rds", ct_label)))

    profiles_c2 <- buildProfiles(
      pairs_c2, structures, ue$union_exons, ue$isoform_union_mapping,
      verify = TRUE,
      checkpoint_dir = file.path(ckpt_base, sprintf("c2_%s", ct_label)),
      checkpoint_interval = 500L
    )
    saveRDS(profiles_c2, file.path(output_dir, sprintf("profiles_c2_%s.rds", ct_label)))
    cat(sprintf("  C2: %d profiles, %d genes\n",
                nrow(profiles_c2), n_distinct(profiles_c2$gene_id)))
  }
}

cat("\n=== Profile building complete ===\n")
