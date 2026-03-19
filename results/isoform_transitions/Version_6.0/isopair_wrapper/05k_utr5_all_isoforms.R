#!/usr/bin/env Rscript
# ==============================================================================
# 05k_utr5_all_isoforms.R
#
# Compute 5'UTR features for ALL classified coding isoforms (~62K), not just
# the paired subset (~10K). This provides the feature set needed for the
# unified predictive model.
#
# Input files:
#   - data_mashr/structures.rds
#   - data_mashr/cds.rds
#   - data_mashr/nmd_classification.rds
#   - FASTA: sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta
#
# Output:
#   - data_mashr/analysis_cache/utr5_features_all.rds
#
# Usage:
#   Rscript 05k_utr5_all_isoforms.R
# ==============================================================================

suppressPackageStartupMessages({
  library(Isopair)
  library(dplyr)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

FASTA_PATH <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta"
OUTPUT_PATH <- "data_mashr/analysis_cache/utr5_features_all.rds"

cat("=== 5'UTR Feature Scan: All Classified Coding Isoforms ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Load data and identify target isoforms
# ==============================================================================
cat("Loading data files...\n")

structures <- readRDS("data_mashr/structures.rds")
cds <- readRDS("data_mashr/cds.rds")
nmd <- readRDS("data_mashr/nmd_classification.rds")

# All classified isoform IDs from all_samples
nmd_ids <- nmd$all_samples$nmd
non_nmd_ids <- nmd$all_samples$non_nmd
classified_ids <- union(nmd_ids, non_nmd_ids)

cat("  Classified isoforms (all_samples): ", length(classified_ids), "\n")
cat("    NMD: ", length(nmd_ids), "\n")
cat("    non-NMD: ", length(non_nmd_ids), "\n")

# Filter to coding only
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
target_ids <- intersect(classified_ids, coding_ids)

cat("  Coding isoforms in CDS metadata: ", length(coding_ids), "\n")
cat("  Target isoforms (classified & coding): ", length(target_ids), "\n\n")

# Filter structures and cds to target isoforms
structures_filtered <- structures[structures$isoform_id %in% target_ids, ]
cds_filtered <- cds[cds$isoform_id %in% target_ids, ]

n_in_structures <- length(unique(structures_filtered$isoform_id))
n_in_cds <- length(unique(cds_filtered$isoform_id))
cat("  Isoforms in filtered structures: ", n_in_structures, "\n")
cat("  Isoforms in filtered cds: ", n_in_cds, "\n")

# Final target = intersection of all three
target_ids <- intersect(target_ids, unique(structures_filtered$isoform_id))
target_ids <- intersect(target_ids, unique(cds_filtered$isoform_id))
cat("  Final target isoforms (in structures & cds): ", length(target_ids), "\n\n")

# ==============================================================================
# 2. Extract sequences from FASTA
# ==============================================================================
cat("Extracting sequences from FASTA...\n")

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
  seq_vec <- setNames(seq_seqs, seq_names)
}

n_with_seq <- sum(target_ids %in% names(seq_vec))
n_missing_seq <- sum(!target_ids %in% names(seq_vec))
cat("  Found", length(seq_vec), "sequences\n")
cat("  Target IDs with sequence:", n_with_seq, "\n")
cat("  Target IDs missing sequence:", n_missing_seq, "\n\n")

# ==============================================================================
# 3. Scan 5'UTR features
# ==============================================================================
cat("Scanning 5'UTR features for", length(target_ids), "isoforms...\n")

t_scan <- proc.time()
isoform_features <- scan5UtrFeatures(
  structures = structures_filtered,
  cds_metadata = cds_filtered,
  sequences = seq_vec,
  min_orf_nt = 9L,
  verbose = TRUE
)
elapsed <- round((proc.time() - t_scan)[3], 1)
cat("  Feature scan completed in", elapsed, "seconds\n")
cat("  Features computed for", nrow(isoform_features), "isoforms\n")

# ==============================================================================
# 4. ATG validation and exclusion flag
# ==============================================================================
n_validated <- sum(!is.na(isoform_features$atg_validated))
n_pass <- sum(isoform_features$atg_validated, na.rm = TRUE)
n_fail <- sum(!is.na(isoform_features$atg_validated) & !isoform_features$atg_validated)
cat("  ATG validation:", n_pass, "/", n_validated,
    "(", round(100 * n_pass / n_validated, 1), "%)\n")

isoform_features$excluded <- isoform_features$utr5_length == 0L |
  (!is.na(isoform_features$atg_validated) & !isoform_features$atg_validated)

n_zero_utr <- sum(isoform_features$utr5_length == 0L)
n_atg_fail <- sum(!is.na(isoform_features$atg_validated) & !isoform_features$atg_validated)
n_excluded <- sum(isoform_features$excluded)
n_included <- sum(!isoform_features$excluded)

cat("  Excluded:", n_excluded,
    "(", n_zero_utr, "zero-length,",
    n_atg_fail, "ATG validation failed)\n")
cat("  Included (usable):", n_included, "\n\n")

# ==============================================================================
# 5. Add NMD group labels
# ==============================================================================
cat("Adding group labels...\n")

isoform_features$group <- ifelse(
  isoform_features$isoform_id %in% nmd_ids, "NMD",
  ifelse(isoform_features$isoform_id %in% non_nmd_ids, "non_NMD", NA_character_)
)

group_counts <- table(isoform_features$group, useNA = "ifany")
cat("  Group counts:\n")
for (g in names(group_counts)) {
  cat("    ", g, ":", group_counts[g], "\n")
}

# Usable by group
usable <- isoform_features[!isoform_features$excluded, ]
usable_counts <- table(usable$group, useNA = "ifany")
cat("  Usable (non-excluded) by group:\n")
for (g in names(usable_counts)) {
  cat("    ", g, ":", usable_counts[g], "\n")
}
cat("\n")

# ==============================================================================
# 6. Save results
# ==============================================================================
output <- list(
  isoform_features = isoform_features,
  metadata = list(
    run_timestamp = Sys.time(),
    script_path = "05k_utr5_all_isoforms.R",
    fasta_path = FASTA_PATH,
    n_classified = length(classified_ids),
    n_coding = length(coding_ids),
    n_target = length(target_ids),
    n_with_sequence = n_with_seq,
    n_missing_sequence = n_missing_seq,
    n_features_computed = nrow(isoform_features),
    n_excluded = n_excluded,
    n_zero_utr = n_zero_utr,
    n_atg_fail = n_atg_fail,
    n_included = n_included,
    group_counts = group_counts,
    usable_group_counts = usable_counts,
    elapsed_scan_seconds = elapsed
  )
)

saveRDS(output, OUTPUT_PATH)
cat("Saved results to:", OUTPUT_PATH, "\n")
cat("Finished:", format(Sys.time()), "\n")
