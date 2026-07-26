#!/usr/bin/env Rscript
# ==============================================================================
# 05s_orfik_scan.R
#
# Comprehensive ORF detection using ORFik on all expressed isoforms.
# Computes per-ORF features: Kozak context, ORF length, position in transcript,
# upstream ATG count, downstream EJC count.
#
# Inputs:
#   - SQANTI corrected FASTA (spliced transcript sequences)
#   - data_mashr/structures.rds (exon coordinates for EJC positions)
#   - data_mashr/nmd_classification.rds (NMD/non-NMD labels)
#   - data_mashr/cds.rds (for population filtering only)
#
# Output:
#   - data_mashr/analysis_cache/orfik_scan.rds
#     $orfs_per_tx: GRangesList of all ORFs per transcript
#     $orf_features: data.frame with per-ORF features
#     $tx_summary: data.frame with per-transcript summary
#
# Usage:
#   Rscript 05s_orfik_scan.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ORFik)
  library(Biostrings)
  library(GenomicRanges)
  library(dplyr)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

FASTA_PATH <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta"
OUTPUT_PATH <- "data_mashr/analysis_cache/orfik_scan.rds"

cat("=== ORFik Comprehensive ORF Scan ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Define population
# ==============================================================================
cat("Loading metadata...\n")

nmd_class <- readRDS("data_mashr/nmd_classification.rds")
cds <- readRDS("data_mashr/cds.rds")
structures <- readRDS("data_mashr/structures.rds")

nmd_ids <- nmd_class$all_samples$nmd
non_nmd_ids <- nmd_class$all_samples$non_nmd
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
target_ids <- intersect(c(nmd_ids, non_nmd_ids), coding_ids)
cat("  Target isoforms:", length(target_ids), "\n")

# ==============================================================================
# 2. Load transcript sequences
# ==============================================================================
cat("\nLoading FASTA...\n")
t0 <- proc.time()
all_seqs <- readDNAStringSet(FASTA_PATH)
cat("  Loaded", length(all_seqs), "sequences in",
    round((proc.time() - t0)[3], 1), "seconds\n")

# Filter to target isoforms
target_seqs <- all_seqs[names(all_seqs) %in% target_ids]
cat("  Target sequences:", length(target_seqs), "\n")
cat("  Missing:", sum(!target_ids %in% names(target_seqs)), "\n")
rm(all_seqs)

# ==============================================================================
# 3. Find all ORFs with ORFik
# ==============================================================================
cat("\nFinding ORFs with ORFik...\n")
t0 <- proc.time()
orfs <- findORFs(target_seqs, startCodon = "ATG",
                  longestORF = FALSE, minimumLength = 9)
cat("  Completed in", round((proc.time() - t0)[3], 1), "seconds\n")
cat("  Transcripts with ORFs:", length(orfs), "\n")
cat("  Total ORFs:", sum(lengths(orfs)), "\n")
cat("  Median ORFs per transcript:", median(lengths(orfs)), "\n")

# ==============================================================================
# 4. Compute per-ORF features
# ==============================================================================
cat("\nComputing per-ORF features...\n")

# Precompute junction positions in transcript space
cat("  Computing junction positions...\n")
target_structs <- structures[structures$isoform_id %in% names(target_seqs), ]
junction_map <- setNames(vector("list", nrow(target_structs)), target_structs$isoform_id)
for (i in seq_len(nrow(target_structs))) {
  starts <- target_structs$exon_starts[[i]]
  ends <- target_structs$exon_ends[[i]]
  strand <- target_structs$strand[i]
  n_exons <- length(starts)
  if (n_exons <= 1) { junction_map[[i]] <- integer(0); next }
  exon_lengths <- ends - starts + 1L
  if (strand == "-") exon_lengths <- rev(exon_lengths)
  junction_map[[i]] <- cumsum(exon_lengths)[-n_exons]
}

# Build feature table
cat("  Computing features for each ORF...\n")
t0 <- proc.time()

# ORFik findORFs returns a CompressedIRangesList that DROPS transcripts with
# zero ORFs. The names are the original numeric indices (as strings) into the
# input DNAStringSet. Use these indices to map back to transcript names.
# BUG FIX: seq_along(orfs) gives 1:length(orfs), which misaligns when
# transcripts are dropped. Use as.integer(names(orfs)) instead.
names(orfs) <- names(target_seqs)[as.integer(names(orfs))]

feature_rows <- vector("list", length(orfs))
tx_names <- names(orfs)

for (tx_idx in seq_along(orfs)) {
  tx_id <- tx_names[tx_idx]
  tx_orfs <- orfs[[tx_idx]]
  n_orfs <- length(tx_orfs)
  if (n_orfs == 0) next

  if (!tx_id %in% names(target_seqs)) next
  tx_seq <- target_seqs[[tx_id]]
  if (is.null(tx_seq) || length(tx_seq) == 0) next
  tx_len <- length(tx_seq)
  junctions <- junction_map[[tx_id]]
  if (is.null(junctions)) junctions <- integer(0)

  orf_starts <- start(tx_orfs)
  orf_ends <- end(tx_orfs)
  orf_widths <- width(tx_orfs)

  # Per-ORF features
  kozak_m3 <- character(n_orfs)
  kozak_p4 <- character(n_orfs)
  kozak_score <- integer(n_orfs)
  n_upstream_atgs <- integer(n_orfs)
  n_downstream_ejc <- integer(n_orfs)
  start_codon <- character(n_orfs)
  stop_codon <- character(n_orfs)

  # Find all ATG positions in transcript for upstream ATG counting
  atg_matches <- tryCatch(matchPattern("ATG", tx_seq), error = function(e) NULL)
  if (is.null(atg_matches)) next
  all_atg_pos <- start(atg_matches)

  for (j in seq_len(n_orfs)) {
    os <- orf_starts[j]
    oe <- orf_ends[j]

    # Bounds check
    if (os < 1L || oe > tx_len || os > oe) next

    # Start codon (should be ATG)
    if (os + 2L <= tx_len) {
      start_codon[j] <- as.character(subseq(tx_seq, os, os + 2L))
    } else {
      start_codon[j] <- NA_character_
    }

    # Stop codon: last 3 nt of the ORF (ORFik includes the stop codon in the range)
    if (oe >= 3L && oe <= tx_len) {
      stop_codon[j] <- as.character(subseq(tx_seq, oe - 2L, oe))
    } else {
      stop_codon[j] <- NA_character_
    }

    # Kozak context: -3 position (A/G = strong) and +4 position (G = strong)
    if (os >= 4L) {
      kozak_m3[j] <- as.character(subseq(tx_seq, os - 3L, os - 3L))
    } else {
      kozak_m3[j] <- NA_character_
    }
    if (os + 3L >= 1L && os + 3L <= tx_len) {
      kozak_p4[j] <- as.character(subseq(tx_seq, os + 3L, os + 3L))
    } else {
      kozak_p4[j] <- NA_character_
    }

    # Kozak score: 0=weak, 1=moderate, 2=strong
    has_m3 <- !is.na(kozak_m3[j]) && kozak_m3[j] %in% c("A", "G")
    has_p4 <- !is.na(kozak_p4[j]) && kozak_p4[j] == "G"
    kozak_score[j] <- as.integer(has_m3) + as.integer(has_p4)

    # Upstream ATGs (ATG positions before this ORF's start)
    n_upstream_atgs[j] <- sum(all_atg_pos < os)

    # Downstream EJCs: junctions after the stop codon
    # ORFik ORF end = last nt of ORF (including stop codon). Stop ends at oe.
    n_downstream_ejc[j] <- sum(junctions > oe)
  }

  feature_rows[[tx_idx]] <- data.frame(
    isoform_id = tx_id,
    orf_idx = seq_len(n_orfs),
    orf_start = orf_starts,
    orf_end = orf_ends,
    orf_length = orf_widths,
    tx_length = tx_len,
    frac_position = orf_starts / tx_len,
    frac_tx_covered = orf_widths / tx_len,
    start_codon = start_codon,
    stop_codon = stop_codon,
    kozak_m3 = kozak_m3,
    kozak_p4 = kozak_p4,
    kozak_score = kozak_score,
    n_upstream_atgs = n_upstream_atgs,
    n_downstream_ejc = n_downstream_ejc,
    has_downstream_ejc = n_downstream_ejc > 0L,
    n_junctions = length(junctions),
    stringsAsFactors = FALSE
  )

  if (tx_idx %% 10000 == 0) {
    cat(sprintf("    %d / %d transcripts (%.0f sec)\n",
                tx_idx, length(orfs), (proc.time() - t0)[3]))
  }
}

orf_features <- do.call(rbind, feature_rows)
cat("  Completed in", round((proc.time() - t0)[3], 1), "seconds\n")
cat("  Total ORF feature rows:", nrow(orf_features), "\n")

# ==============================================================================
# 5. Per-transcript summary
# ==============================================================================
cat("\nComputing per-transcript summary...\n")

tx_summary <- orf_features %>%
  group_by(isoform_id) %>%
  summarise(
    n_orfs = n(),
    n_nmd_competent = sum(has_downstream_ejc),
    max_orf_length = max(orf_length),
    max_downstream_ejc = max(n_downstream_ejc),
    n_strong_kozak = sum(kozak_score == 2L),
    n_strong_kozak_nmd_comp = sum(kozak_score == 2L & has_downstream_ejc),
    best_kozak_nmd_comp = max(kozak_score[has_downstream_ejc], 0L),
    longest_nmd_comp_orf = ifelse(any(has_downstream_ejc),
                                   max(orf_length[has_downstream_ejc]), 0L),
    tx_length = first(tx_length),
    n_junctions = first(n_junctions),
    .groups = "drop"
  )

# Add NMD status
tx_summary$is_nmd <- as.integer(tx_summary$isoform_id %in% nmd_ids)

# Add chromosome
chr_map <- structures[, c("isoform_id", "chr")]
tx_summary <- merge(tx_summary, chr_map, by = "isoform_id", all.x = TRUE)

# ==============================================================================
# 6. Summary statistics
# ==============================================================================
cat("\n=== SUMMARY ===\n\n")
cat("Transcripts scanned:", nrow(tx_summary), "\n")
cat("  NMD:", sum(tx_summary$is_nmd), "  Non-NMD:", sum(!tx_summary$is_nmd), "\n\n")

cat("ORF features:\n")
cat("  Total ORFs:", nrow(orf_features), "\n")
cat("  Median ORFs per transcript:", median(tx_summary$n_orfs), "\n")
cat("  Median NMD-competent ORFs:", median(tx_summary$n_nmd_competent), "\n")
cat("  Median max ORF length:", median(tx_summary$max_orf_length), "nt\n")
cat("  ORFs with strong Kozak (score=2):",
    sum(orf_features$kozak_score == 2),
    sprintf("(%.1f%%)\n", 100 * mean(orf_features$kozak_score == 2)))
cat("  ORFs with downstream EJC:",
    sum(orf_features$has_downstream_ejc),
    sprintf("(%.1f%%)\n", 100 * mean(orf_features$has_downstream_ejc)))

cat("\nBy NMD status:\n")
for (status in c(1, 0)) {
  sub <- tx_summary[tx_summary$is_nmd == status, ]
  label <- if (status == 1) "NMD" else "Non-NMD"
  cat(sprintf("  %s (n=%d): median ORFs=%d, median NMD-comp=%d, median max_orf=%d nt\n",
              label, nrow(sub), median(sub$n_orfs), median(sub$n_nmd_competent),
              median(sub$max_orf_length)))
}

# ==============================================================================
# 7. Save
# ==============================================================================
results <- list(
  orf_features = orf_features,
  tx_summary = tx_summary,
  metadata = list(
    n_transcripts = nrow(tx_summary),
    n_orfs = nrow(orf_features),
    min_orf_length = 9L,
    start_codon = "ATG",
    run_timestamp = Sys.time()
  )
)

# ==============================================================================
# 8. Verification: confirm all ORFs start with ATG
# ==============================================================================
cat("\n=== Verification ===\n")
n_atg <- sum(orf_features$start_codon == "ATG", na.rm = TRUE)
n_non_atg <- sum(orf_features$start_codon != "ATG" & !is.na(orf_features$start_codon))
n_na <- sum(is.na(orf_features$start_codon))
cat("  ATG starts:", n_atg, sprintf("(%.1f%%)\n", 100 * n_atg / nrow(orf_features)))
cat("  Non-ATG starts:", n_non_atg, "\n")
cat("  NA starts:", n_na, "\n")
if (n_non_atg > 0 || n_na > 0) {
  stop("VERIFICATION FAILED: Found ", n_non_atg, " non-ATG and ", n_na,
       " NA start codons. Name mapping or coordinate extraction is likely wrong.")
}
cat("  PASSED: All ORFs start with ATG\n")

# Also verify transcript IDs are valid
n_valid_tx <- sum(orf_features$isoform_id %in% names(target_seqs))
cat("  Valid transcript IDs:", n_valid_tx, "/", nrow(orf_features), "\n")
if (n_valid_tx != nrow(orf_features)) {
  stop("VERIFICATION FAILED: ", nrow(orf_features) - n_valid_tx,
       " ORFs have isoform_ids not in the input sequences.")
}
cat("  PASSED: All isoform_ids match input sequences\n")

saveRDS(results, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
