#!/usr/bin/env Rscript
# ==============================================================================
# 05m_orf_landscape.R
#
# ORF landscape feature extraction for NMD prediction.
# Backend: Isopair::enumerateOrfs (Kozak-filtered, MANE-calibrated threshold).
#
# For each classified coding isoform, enumerates plausibly-translated ORFs
# (ATGs above the MANE 5th-percentile Kozak threshold), traces each to the
# first in-frame stop, counts downstream exon-exon junctions, and aggregates
# to per-isoform features.
#
# Key change from prior categorical implementation (2026-04-29):
#   - Kozak filtering ON (kozak_filter=TRUE): only plausibly-translated
#     ATGs are enumerated, eliminating low-quality noise from random ATGs.
#   - Continuous PWM Kozak scores replace the categorical 0/1/2 score.
#   - n_atgs (count of all ATGs) is no longer reported — it was a diagnostic
#     feature that lost meaning under Kozak filtering.
#   - n_nmd_competent_strong_kozak (categorical) replaced by best_kozak_nmd
#     and mean_kozak_nmd (both continuous PWM scores).
#
# Inputs:
#   - data_mashr/nmd_classification.rds (mashr NMD/non-NMD labels)
#   - data_mashr/structures.rds (exon coords → junctions)
#   - data_mashr/cds.rds (population filter and is_annotated_cds annotation)
#   - SQANTI corrected FASTA (spliced transcript sequences)
#
# Output:
#   - data_mashr/analysis_cache/orf_landscape.rds (per-isoform aggregates)
#   - data_mashr/analysis_cache/orf_landscape_per_orf.rds (raw enumerateOrfs output)
#
# Usage:
#   Rscript 05m_orf_landscape.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Isopair)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

FASTA_PATH <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta"
OUTPUT_PATH <- "data_mashr/analysis_cache/orf_landscape.rds"
PER_ORF_PATH <- "data_mashr/analysis_cache/orf_landscape_per_orf.rds"

cat("=== ORF Landscape Feature Extraction (Isopair::enumerateOrfs) ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Load data and define population
# ==============================================================================
cat("Loading data...\n")

nmd_class <- readRDS("data_mashr/nmd_classification.rds")
cds <- readRDS("data_mashr/cds.rds")
structures <- readRDS("data_mashr/structures.rds")

nmd_ids <- nmd_class$all_samples$nmd
non_nmd_ids <- nmd_class$all_samples$non_nmd
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]

# Population: all classified coding isoforms (same as unified model)
target_ids <- intersect(c(nmd_ids, non_nmd_ids), coding_ids)
cat("  Target isoforms:", length(target_ids), "\n")

# ==============================================================================
# 2. Extract transcript sequences from FASTA
# ==============================================================================
cat("\nExtracting sequences from FASTA...\n")

id_file <- tempfile(fileext = ".txt")
writeLines(target_ids, id_file)

awk_cmd <- sprintf(
  "awk 'BEGIN{while((getline id < \"%s\") > 0) want[id]=1} /^>/{if(seq && found) print hdr \"\\t\" seq; hdr=substr($1,2); found=(hdr in want); seq=\"\"; next} found{seq=seq$0} END{if(seq && found) print hdr \"\\t\" seq}' '%s'",
  id_file, FASTA_PATH
)

t_awk <- proc.time()
raw_output <- system(awk_cmd, intern = TRUE)
cat("  awk completed in", round((proc.time() - t_awk)[3], 1), "seconds\n")
unlink(id_file)

seq_vec <- character(0)
if (length(raw_output) > 0) {
  split_lines <- strsplit(raw_output, "\t", fixed = TRUE)
  seq_names <- vapply(split_lines, `[`, character(1), 1)
  seq_seqs <- vapply(split_lines, `[`, character(1), 2)
  seq_vec <- setNames(toupper(seq_seqs), seq_names)
}

cat("  Sequences extracted:", length(seq_vec), "\n")
cat("  Missing:", sum(!target_ids %in% names(seq_vec)), "\n")

# Filter to isoforms with sequences
target_ids <- intersect(target_ids, names(seq_vec))
target_structs <- structures[structures$isoform_id %in% target_ids, ]
cat("  Final target:", length(target_ids), "\n")

# ==============================================================================
# 3. Enumerate ORFs (Isopair::enumerateOrfs)
# ==============================================================================
cat("\nEnumerating ORFs...\n")

t_enum <- proc.time()
per_orf <- Isopair::enumerateOrfs(
  structures   = target_structs,
  cds_metadata = cds,
  sequences    = seq_vec[target_ids],
  min_orf_nt   = 30L,         # 10 aa minimum
  ejc_threshold = 50L,        # canonical NMD rule
  include_no_stop = TRUE,     # emit category="no_stop_in_frame" rows
  kozak_filter = TRUE,        # Kozak-gate ATGs (MANE-calibrated default threshold)
  kozak_threshold = NULL      # defaultKozakThreshold(0.05): admits 95% of MANE CDS starts
)
cat(sprintf("  enumerateOrfs completed in %.0f seconds\n", (proc.time() - t_enum)[3]))
cat(sprintf("  Per-ORF rows: %d (%.1f per isoform avg)\n",
            nrow(per_orf), nrow(per_orf) / length(target_ids)))

saveRDS(per_orf, PER_ORF_PATH)

# ==============================================================================
# 4. Roll up to per-isoform features
# ==============================================================================
cat("\nRolling up per-isoform features...\n")

# Per-isoform tx_length and junction count from structures (not in enumerateOrfs output)
tx_meta <- target_structs %>%
  rowwise() %>%
  mutate(tx_length = sum(exon_ends - exon_starts + 1L),
         n_junctions = max(0L, length(exon_starts) - 1L)) %>%
  ungroup() %>%
  select(isoform_id, tx_length, n_junctions)

# Per-isoform ORF features. With include_no_stop=TRUE, no_stop rows have NA
# stop_tx_pos / orf_length. Treat them as a separate count and otherwise
# aggregate over rows that DO have a stop.
orf_df <- per_orf %>%
  group_by(isoform_id) %>%
  summarise(
    # Counts (Kozak-passing ORFs only; n_orfs_nonstop = ATGs with no in-frame stop)
    n_orfs_nonstop  = sum(category == "no_stop_in_frame"),
    n_orfs          = sum(category != "no_stop_in_frame"),
    n_nmd_competent = sum(category == "effectively_ptc"),
    has_any_nmd_competent = any(category == "effectively_ptc"),

    # Magnitudes (NA-safe — empty subsets return -Inf/+Inf without na.rm gating)
    max_downstream_ejc = if (any(category == "effectively_ptc"))
                          max(n_downstream_ejc[category == "effectively_ptc"]) else 0L,
    best_kozak_nmd     = if (any(category == "effectively_ptc"))
                          max(kozak_score[category == "effectively_ptc"]) else NA_real_,
    mean_kozak_nmd     = if (any(category == "effectively_ptc"))
                          mean(kozak_score[category == "effectively_ptc"]) else NA_real_,
    first_nmd_atg_pos  = if (any(category == "effectively_ptc"))
                          min(atg_tx_pos[category == "effectively_ptc"]) else NA_integer_,
    n_annotated_cds_orfs = sum(is_annotated_cds, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  left_join(tx_meta, by = "isoform_id")

# Isoforms with zero Kozak-passing ATGs aren't represented in per_orf — add zero-rows
missing_ids <- setdiff(target_ids, orf_df$isoform_id)
if (length(missing_ids) > 0) {
  zero_rows <- tibble::tibble(
    isoform_id = missing_ids,
    n_orfs_nonstop = 0L, n_orfs = 0L, n_nmd_competent = 0L,
    has_any_nmd_competent = FALSE,
    max_downstream_ejc = 0L,
    best_kozak_nmd = NA_real_, mean_kozak_nmd = NA_real_,
    first_nmd_atg_pos = NA_integer_,
    n_annotated_cds_orfs = 0L
  ) %>% left_join(tx_meta, by = "isoform_id")
  orf_df <- bind_rows(orf_df, zero_rows)
  cat(sprintf("  Added %d zero-ORF rows (isoforms with no Kozak-passing ATGs)\n",
              length(missing_ids)))
}

# Derived features (matches prior schema where applicable)
orf_df$frac_nmd_competent       <- orf_df$n_nmd_competent / pmax(orf_df$n_orfs, 1L)
orf_df$first_nmd_atg_frac       <- orf_df$first_nmd_atg_pos / orf_df$tx_length
orf_df$n_nmd_competent_capped   <- pmin(orf_df$n_nmd_competent, 20L)
orf_df$max_downstream_ejc_capped <- pmin(orf_df$max_downstream_ejc, 5L)

# NMD label
orf_df$is_nmd <- as.integer(orf_df$isoform_id %in% nmd_ids)

# Chromosome for train/test split
chr_map <- structures[, c("isoform_id", "chr")]
orf_df <- merge(orf_df, chr_map, by = "isoform_id", all.x = TRUE)

# ==============================================================================
# 5. Summary
# ==============================================================================
cat("\n=== SUMMARY ===\n\n")
cat("Total isoforms:", nrow(orf_df), "\n")
cat("  NMD:", sum(orf_df$is_nmd), "  Non-NMD:", sum(!orf_df$is_nmd), "\n\n")

cat("ORF landscape (Kozak-filtered):\n")
cat("  Median ORFs per isoform:", median(orf_df$n_orfs), "\n")
cat("  Median NMD-competent ORFs:", median(orf_df$n_nmd_competent), "\n")
cat("  Isoforms with any NMD-competent ORF:",
    sum(orf_df$has_any_nmd_competent),
    sprintf("(%.1f%%)\n", 100 * mean(orf_df$has_any_nmd_competent)))

cat("\nBy NMD status:\n")
for (status in c(1, 0)) {
  sub <- orf_df[orf_df$is_nmd == status, ]
  label <- if (status == 1) "NMD" else "Non-NMD"
  cat(sprintf("  %s (n=%d): median ORFs=%g, median NMD-competent=%g, any NMD-competent=%.1f%%\n",
              label, nrow(sub), median(sub$n_orfs), median(sub$n_nmd_competent),
              100 * mean(sub$has_any_nmd_competent)))
}

# ==============================================================================
# 6. Save
# ==============================================================================
attr(orf_df, "metadata") <- list(
  run_timestamp = Sys.time(),
  fasta_path = FASTA_PATH,
  isopair_version = as.character(packageVersion("Isopair")),
  backend = "Isopair::enumerateOrfs",
  kozak_filter = TRUE,
  kozak_threshold = "defaultKozakThreshold(0.05)",
  schema_change_note = "n_atgs and n_nmd_competent_strong_kozak removed; mean_kozak_nmd added; PWM continuous Kozak scores"
)

saveRDS(orf_df, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Per-ORF table also saved to:", PER_ORF_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
