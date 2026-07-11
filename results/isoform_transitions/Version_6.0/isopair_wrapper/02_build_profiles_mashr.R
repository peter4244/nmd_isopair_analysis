#!/usr/bin/env Rscript
#
# NMD Isoform Transitions — Profile Building (mashr Classification)
# Version: 7.0-mashr
#
# Purpose: Generate C2/C4 pairs and build splicing choice profiles using mashr
#          NMD classifications. Reuses isoform-intrinsic infrastructure
#          (structures, union exons, CDS, PTC) from data/ and caches profiles
#          by pair identity to avoid recomputing already-profiled pairs.
#
# Input:
#   - RDS files from 01_prepare_data_mashr.R (data_mashr/)
#   - Isoform infrastructure (data/) — DE-method-independent
#
# Output (to output_dir/):
#   - Per cell type: pairs_c2_{ct}.rds, pairs_c4_{ct}.rds,
#                     profiles_c2_{ct}.rds, profiles_c4_{ct}.rds
#   - Symlinked/copied infrastructure files for self-contained data_mashr/
#
# Usage:
#   Rscript 02_build_profiles_mashr.R [--data-dir DIR] [--output-dir DIR] [--infra-dir DIR]

library(Isopair)
library(dplyr)

# ==============================================================================
# 0. Configuration
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
data_dir   <- "data_mashr"
output_dir <- "data_mashr"
infra_dir  <- "data_mashr"  # isoform-intrinsic infra was copied here; legacy data/ archived
if ("--data-dir" %in% args) {
  idx <- which(args == "--data-dir")
  if (idx < length(args)) data_dir <- args[idx + 1]
}
if ("--output-dir" %in% args) {
  idx <- which(args == "--output-dir")
  if (idx < length(args)) output_dir <- args[idx + 1]
}
if ("--infra-dir" %in% args) {
  idx <- which(args == "--infra-dir")
  if (idx < length(args)) infra_dir <- args[idx + 1]
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Cell types to process — paper scope (2026-04-29): 4 non-ALI cell types
cell_types <- c("all_samples", "AT", "DD", "FB", "MV")

# Min pairs threshold
MIN_PAIRS <- 50

# Reference-share floor (Family A, 2026-07-10): retain only genes whose selected
# reference isoform accounts for >= this fraction of TOTAL gene DMSO expression.
# Reference *selection* is unchanged (still rank-1 of the strict non-NMD pool);
# this is a post-selection retention filter. See REFERENCE_FLOOR_RATIONALE.md.
REF_SHARE_FLOOR <- 0.25

cat("=== NMD Profile Building — mashr Classification ===\n\n")

# ==============================================================================
# 1. Load prepared data and infrastructure
# ==============================================================================

cat("Loading mashr-classified data...\n")
expr_mat   <- readRDS(file.path(data_dir, "expression_data.rds"))
gene_map   <- readRDS(file.path(data_dir, "gene_map.rds"))
nmd_class  <- readRDS(file.path(data_dir, "nmd_classification.rds"))
dmso_samp  <- readRDS(file.path(data_dir, "dmso_samples.rds"))
smg1i_samp <- readRDS(file.path(data_dir, "smg1i_samples.rds"))

cat(sprintf("  Expression: %d isoforms x %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))

# 4-CT re-scope guard (2026-07-11, decision b — REFERENCE_FLOOR_PLAN.md R2):
# the all_samples sample basis must be the 4 manuscript cell types only
# (13 DMSO + 13 Smg1i = 26; no DD_ALI/DO), so reference selection, the 25%
# floor denominator, and the C2 NMD partner are all 4-CT by construction.
stopifnot(
  ncol(expr_mat) == 26,
  length(dmso_samp[["all_samples"]])  == 13,
  length(smg1i_samp[["all_samples"]]) == 13,
  !any(grepl("DD_ALI|DO_ALI|_DO_",
             c(dmso_samp[["all_samples"]], smg1i_samp[["all_samples"]])))
)

# Load isoform-intrinsic infrastructure (DE-method-independent)
cat("Loading isoform infrastructure...\n")
structures   <- readRDS(file.path(infra_dir, "structures.rds"))
ue_union     <- readRDS(file.path(infra_dir, "union_exons.rds"))
ue_mapping   <- readRDS(file.path(infra_dir, "isoform_union_mapping.rds"))
cds          <- readRDS(file.path(infra_dir, "cds.rds"))
annotated_ue <- readRDS(file.path(infra_dir, "annotated_ue.rds"))
ptc          <- readRDS(file.path(infra_dir, "ptc.rds"))

# Copy infrastructure to output dir for self-contained data_mashr/
infra_files <- c("structures.rds", "union_exons.rds", "isoform_union_mapping.rds",
                 "cds.rds", "cds_exons.rds", "annotated_ue.rds", "ptc.rds")
for (f in infra_files) {
  src <- file.path(infra_dir, f)
  dst <- file.path(output_dir, f)
  if (!file.exists(dst) && file.exists(src)) {
    file.copy(src, dst)
    cat(sprintf("  Copied %s to %s\n", f, output_dir))
  }
}

# Gene map as tibble for Isopair
gene_map_tbl <- tibble::tibble(
  isoform_id = gene_map$isoform_id,
  gene_id = gene_map$gene_id
)
# Restrict to isoforms in structures
gene_map_tbl <- gene_map_tbl[gene_map_tbl$isoform_id %in% structures$isoform_id, ]

# ==============================================================================
# 2. Build profile lookup from existing cached profiles
# ==============================================================================

cat("\nBuilding profile lookup from cached profiles...\n")

# Load all existing profiles into a single lookup keyed by ref|comp
# (cache is by pair identity — cell-type label is irrelevant for reuse)
profile_lookup <- list()
cached_profile_files <- list.files(infra_dir, pattern = "^profiles_c[24]_.*\\.rds$",
                                    full.names = TRUE)
n_cached <- 0L
for (pf in cached_profile_files) {
  profiles <- readRDS(pf)
  for (i in seq_len(nrow(profiles))) {
    key <- paste(profiles$reference_isoform_id[i],
                 profiles$comparator_isoform_id[i], sep = "|")
    profile_lookup[[key]] <- profiles[i, ]
  }
  n_cached <- n_cached + nrow(profiles)
}
cat(sprintf("  Loaded %d cached profiles from %d profile files\n",
            n_cached, length(cached_profile_files)))

# ==============================================================================
# 3. Per cell type: Generate pairs and build profiles (with reuse)
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

  # --- Reference-share floor (Family A): drop genes whose selected reference
  #     is < REF_SHARE_FLOOR of the gene's total DMSO expression across ALL
  #     isoforms (gene_map_tbl = 5%-condition-filtered set), same DMSO basis as
  #     selection. Applied BEFORE the MIN_PAIRS gate so the gate sees the
  #     floored count. C2 inherits the surviving references via its gene_id join.
  if (!is.null(pairs_c4) && nrow(pairs_c4) > 0) {
    dmso_mean_all <- rowMeans(expr_mat[, dmso_samp[[ct]], drop = FALSE])
    gm_all <- gene_map_tbl[gene_map_tbl$isoform_id %in% names(dmso_mean_all), ]
    gm_all$mean_dmso <- dmso_mean_all[gm_all$isoform_id]
    gene_total_dmso <- tapply(gm_all$mean_dmso, gm_all$gene_id, sum)
    ref_dmso  <- dmso_mean_all[pairs_c4$reference_isoform_id]
    gene_tot  <- gene_total_dmso[pairs_c4$gene_id]
    ref_share <- ifelse(!is.na(gene_tot) & gene_tot > 0, ref_dmso / gene_tot, NA_real_)
    keep <- !is.na(ref_share) & ref_share >= REF_SHARE_FLOOR
    n_before <- nrow(pairs_c4)
    pairs_c4 <- pairs_c4[keep, , drop = FALSE]
    cat(sprintf("  Reference-share floor (>=%.0f%% of gene DMSO): %d -> %d C4 pairs (dropped %d)\n",
                100 * REF_SHARE_FLOOR, n_before, nrow(pairs_c4), n_before - nrow(pairs_c4)))
  }

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

  # --- Build profiles with reuse ---
  ct_label <- tolower(gsub("_", "", ct))
  ckpt_base <- file.path(output_dir, "checkpoints")

  build_with_reuse <- function(pairs, comp_label) {
    # Split into cached and novel pairs
    pair_keys <- paste(pairs$reference_isoform_id,
                       pairs$comparator_isoform_id, sep = "|")
    cached_mask <- pair_keys %in% names(profile_lookup)
    n_cached_hit <- sum(cached_mask)
    n_novel <- sum(!cached_mask)

    cat(sprintf("  %s: %d pairs total, %d cached (%.1f%%), %d novel\n",
                comp_label, nrow(pairs), n_cached_hit,
                100 * n_cached_hit / nrow(pairs), n_novel))

    # Retrieve cached profiles
    cached_profiles <- NULL
    if (n_cached_hit > 0) {
      cached_profiles <- dplyr::bind_rows(
        profile_lookup[pair_keys[cached_mask]]
      )
    }

    # Build novel profiles
    novel_profiles <- NULL
    if (n_novel > 0) {
      novel_pairs <- pairs[!cached_mask, ]
      novel_profiles <- buildProfiles(
        novel_pairs, structures, ue_union, ue_mapping,
        verify = TRUE,
        checkpoint_dir = file.path(ckpt_base,
                                    sprintf("%s_%s", comp_label, ct_label)),
        checkpoint_interval = 500L
      )
    }

    # Combine
    result <- dplyr::bind_rows(cached_profiles, novel_profiles)

    # Add any novel profiles to the global lookup for reuse within this run
    if (!is.null(novel_profiles) && nrow(novel_profiles) > 0) {
      for (i in seq_len(nrow(novel_profiles))) {
        key <- paste(novel_profiles$reference_isoform_id[i],
                     novel_profiles$comparator_isoform_id[i], sep = "|")
        profile_lookup[[key]] <<- novel_profiles[i, ]
      }
    }

    result
  }

  if (!is.null(pairs_c4)) {
    cat(sprintf("  Building C4 profiles (%d pairs)...\n", nrow(pairs_c4)))
    saveRDS(pairs_c4, file.path(output_dir, sprintf("pairs_c4_%s.rds", ct_label)))

    profiles_c4 <- build_with_reuse(pairs_c4, "c4")
    saveRDS(profiles_c4, file.path(output_dir, sprintf("profiles_c4_%s.rds", ct_label)))
    cat(sprintf("  C4: %d profiles, %d genes\n",
                nrow(profiles_c4), n_distinct(profiles_c4$gene_id)))
  }

  if (!is.null(pairs_c2)) {
    cat(sprintf("  Building C2 profiles (%d pairs)...\n", nrow(pairs_c2)))
    saveRDS(pairs_c2, file.path(output_dir, sprintf("pairs_c2_%s.rds", ct_label)))

    profiles_c2 <- build_with_reuse(pairs_c2, "c2")
    saveRDS(profiles_c2, file.path(output_dir, sprintf("profiles_c2_%s.rds", ct_label)))
    cat(sprintf("  C2: %d profiles, %d genes\n",
                nrow(profiles_c2), n_distinct(profiles_c2$gene_id)))
  }
}

cat("\n=== Profile building complete ===\n")
