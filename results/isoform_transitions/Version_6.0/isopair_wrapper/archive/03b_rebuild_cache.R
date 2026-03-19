#!/usr/bin/env Rscript
#
# NMD Isoform Transitions — Rebuild Analysis Cache
#
# Purpose: Regenerate all analysis_cache/ objects from gene-matched profiles.
#          Standalone script that replicates the caching logic from
#          03_nmd_analysis_mashr.Rmd without rendering the full exploratory report.
#
# Prerequisite: 01_prepare_data_mashr.R and 02_build_profiles_mashr.R must have
#               been run first to produce data_mashr/*.rds files.
#
# Output: data_mashr/analysis_cache/{div,fc,fw,ptc,ri,er,cooc,cps}_*.rds
#
# Usage:
#   Rscript 03b_rebuild_cache.R [--force]

library(Isopair)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args

data_dir  <- "data_mashr"
cache_dir <- file.path(data_dir, "analysis_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

ct_to_label <- function(ct) tolower(gsub("_", "", ct))

cached_compute <- function(key, expr, force_recompute = force) {
  cache_file <- file.path(cache_dir, paste0(key, ".rds"))
  if (!force_recompute && file.exists(cache_file)) {
    message(sprintf("  [cache hit] %s", key))
    return(readRDS(cache_file))
  }
  message(sprintf("  [computing] %s ...", key))
  result <- eval(expr)
  saveRDS(result, cache_file)
  result
}

cat("=== Rebuilding Analysis Cache ===\n\n")

# ==============================================================================
# 1. Load data
# ==============================================================================

cat("Loading data...\n")
structures   <- readRDS(file.path(data_dir, "structures.rds"))
cds          <- readRDS(file.path(data_dir, "cds.rds"))
ptc          <- readRDS(file.path(data_dir, "ptc.rds"))
union_exons  <- readRDS(file.path(data_dir, "union_exons.rds"))
ue_mapping   <- readRDS(file.path(data_dir, "isoform_union_mapping.rds"))
annotated_ue <- readRDS(file.path(data_dir, "annotated_ue.rds"))
cds_exons    <- readRDS(file.path(data_dir, "cds_exons.rds"))

cell_types <- c("all_samples", "AT", "DD", "DO", "FB", "MV")
comparisons <- c("c2", "c4")

# ==============================================================================
# 2. Load profiles and apply gene-matching
# ==============================================================================

cat("Loading profiles...\n")
profiles_list <- list()
for (ct in cell_types) {
  ct_label <- ct_to_label(ct)
  for (comp in comparisons) {
    key <- paste(comp, ct_label, sep = "_")
    prof_file <- file.path(data_dir, sprintf("profiles_%s.rds", key))
    if (file.exists(prof_file)) profiles_list[[key]] <- readRDS(prof_file)
  }
}

cat("Applying gene-matching...\n")
for (ct in cell_types) {
  ct_label <- ct_to_label(ct)
  c2_key <- paste0("c2_", ct_label)
  c4_key <- paste0("c4_", ct_label)
  if (!c2_key %in% names(profiles_list) || !c4_key %in% names(profiles_list)) next

  n_c2_before <- nrow(profiles_list[[c2_key]])
  n_c4_before <- nrow(profiles_list[[c4_key]])

  c2_gene_refs <- unique(profiles_list[[c2_key]][, c("gene_id", "reference_isoform_id")])
  c4_gene_refs <- unique(profiles_list[[c4_key]][, c("gene_id", "reference_isoform_id")])
  shared_refs  <- inner_join(c2_gene_refs, c4_gene_refs,
                              by = c("gene_id", "reference_isoform_id"))

  profiles_list[[c2_key]] <- semi_join(profiles_list[[c2_key]], shared_refs,
                                        by = c("gene_id", "reference_isoform_id"))
  profiles_list[[c4_key]] <- semi_join(profiles_list[[c4_key]], shared_refs,
                                        by = c("gene_id", "reference_isoform_id"))

  cat(sprintf("  %s: C2 %d→%d, C4 %d→%d (matched genes: %d)\n",
    ct, n_c2_before, nrow(profiles_list[[c2_key]]),
    n_c4_before, nrow(profiles_list[[c4_key]]), nrow(shared_refs)))
}

# ==============================================================================
# 3. Compute and cache all analysis objects
# ==============================================================================

cat("\nComputing analyses...\n")
for (key in names(profiles_list)) {
  profiles <- profiles_list[[key]]
  cat(sprintf("\n--- %s (%d pairs) ---\n", key, nrow(profiles)))

  cached_compute(paste0("div_", key),
    tryCatch(quantifyPairDivergence(profiles, structures, cds), error = function(e) { message(e$message); NULL }))
  cached_compute(paste0("ri_", key),
    tryCatch(quantifyRegionalImpact(profiles, annotated_ue, cds), error = function(e) { message(e$message); NULL }))
  cached_compute(paste0("cooc_", key),
    tryCatch(testCooccurrence(profiles), error = function(e) { message(e$message); NULL }))
  cached_compute(paste0("er_", key),
    tryCatch(mapEventsToRegions(profiles, annotated_ue, cds), error = function(e) { message(e$message); NULL }))
  cached_compute(paste0("ptc_", key),
    tryCatch(testPtcAssociation(ptc, profiles), error = function(e) { message(e$message); NULL }))
  cached_compute(paste0("fw_", key),
    tryCatch(analyzeFrameWalk(profiles, cds, ptc), error = function(e) { message(e$message); NULL }))
  cached_compute(paste0("fc_", key),
    tryCatch(compareIsoformFrames(profiles, cds_exons, cds), error = function(e) { message(e$message); NULL }))
}

# ==============================================================================
# 4. comparePairSets per cell type
# ==============================================================================

cat("\nComputing NMD vs Control cross-comparisons...\n")
for (ct in cell_types) {
  ct_label <- ct_to_label(ct)
  c2_key <- paste0("c2_", ct_label)
  c4_key <- paste0("c4_", ct_label)
  if (!c2_key %in% names(profiles_list) || !c4_key %in% names(profiles_list)) next

  ptc_c2 <- ptc[ptc$isoform_id %in% profiles_list[[c2_key]]$comparator_isoform_id, ]
  ptc_c4 <- ptc[ptc$isoform_id %in% profiles_list[[c4_key]]$comparator_isoform_id, ]

  cached_compute(paste0("cps_", ct_label),
    tryCatch(
      comparePairSets(profiles_list[[c2_key]], profiles_list[[c4_key]],
                       ptc_a = ptc_c2, ptc_b = ptc_c4,
                       paired = TRUE, verbose = FALSE),
      error = function(e) { message(e$message); NULL }))
}

cat("\n=== Cache rebuild complete ===\n")
