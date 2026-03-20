#!/usr/bin/env Rscript
# ==============================================================================
# 05m_orf_landscape.R
#
# ORF landscape feature extraction for NMD prediction.
#
# For each classified coding isoform, scans the full spliced transcript
# sequence for ALL ORFs (ATG → first in-frame stop, all three frames).
# For each ORF stop codon, counts downstream exon-exon junctions (EJCs).
# An ORF whose stop has ≥1 downstream EJC is "NMD-competent."
#
# This analysis is CDS-annotation-free: it reasons from first principles
# about ribosomal encounters with ORFs and EJCs.
#
# Inputs:
#   - nmd_classification.rds (mashr NMD/non-NMD labels)
#   - structures.rds (exon coordinates → EJC positions)
#   - cds.rds (only for population filtering, not for features)
#   - SQANTI corrected FASTA (spliced transcript sequences)
#
# Output:
#   - data_mashr/analysis_cache/orf_landscape.rds
#
# Usage:
#   Rscript 05m_orf_landscape.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

FASTA_PATH <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta"
OUTPUT_PATH <- "data_mashr/analysis_cache/orf_landscape.rds"

cat("=== ORF Landscape Feature Extraction ===\n")
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
# 2. Compute exon junction positions in transcript space
# ==============================================================================
cat("\nComputing exon junction positions in transcript space...\n")

# For each isoform, compute the transcript-space positions of exon-exon junctions.
# A junction occurs between consecutive exons. In transcript space, junctions
# are at cumulative exon length boundaries.
#
# For + strand: exons are in genomic order (5'→3')
# For - strand: exons are in reverse genomic order (5'→3' in transcript)

target_structs <- structures[structures$isoform_id %in% target_ids, ]

junction_list <- setNames(vector("list", nrow(target_structs)), target_structs$isoform_id)

for (i in seq_len(nrow(target_structs))) {
  starts <- target_structs$exon_starts[[i]]
  ends <- target_structs$exon_ends[[i]]
  strand <- target_structs$strand[i]
  n_exons <- length(starts)

  if (n_exons <= 1) {
    junction_list[[i]] <- integer(0)
    next
  }

  # Exon lengths (1-based inclusive coordinates: length = end - start + 1)
  exon_lengths <- ends - starts + 1L

  # For - strand, reverse exon order (transcript is read 3'→5' genomically)
  if (strand == "-") {
    exon_lengths <- rev(exon_lengths)
  }

  # Junction positions in transcript space: cumulative sum of exon lengths
  # Junction j is at position sum(exon_lengths[1:j]) — the last nt of exon j
  cum_lengths <- cumsum(exon_lengths)
  # Junctions are between exon j and exon j+1, at transcript position cum_lengths[j]
  junction_positions <- cum_lengths[-n_exons]

  junction_list[[i]] <- junction_positions
}

cat("  Computed junctions for", length(junction_list), "isoforms\n")

# Quick check
n_junctions <- sapply(junction_list, length)
cat("  Median junctions per isoform:", median(n_junctions), "\n")

# ==============================================================================
# 3. Extract transcript sequences from FASTA
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
cat("  Final target:", length(target_ids), "\n")

# ==============================================================================
# 4. Scan all ORFs and compute features
# ==============================================================================
cat("\nScanning ORFs...\n")

# Kozak scoring function (simplified, matches Isopair logic)
# Strong Kozak: A/G at -3 AND G at +4 relative to ATG
score_kozak <- function(seq, atg_pos) {
  # atg_pos is 1-based position of the A in ATG
  # -3 position = atg_pos - 3, +4 position = atg_pos + 3 (G after ATG)
  pos_m3 <- atg_pos - 3
  pos_p4 <- atg_pos + 3  # 4th position after A (i.e., position after G of ATG)
  n <- nchar(seq)

  nt_m3 <- if (pos_m3 >= 1) substr(seq, pos_m3, pos_m3) else ""
  nt_p4 <- if (pos_p4 <= n) substr(seq, pos_p4, pos_p4) else ""

  has_m3 <- nt_m3 %in% c("A", "G")
  has_p4 <- nt_p4 == "G"

  if (has_m3 && has_p4) return(2L)  # strong
  if (has_m3 || has_p4) return(1L)  # moderate
  return(0L)  # weak
}

# Stop codons
stop_codons <- c("TAA", "TAG", "TGA")

# Process each isoform
results <- vector("list", length(target_ids))
names(results) <- target_ids

t_scan <- proc.time()
for (idx in seq_along(target_ids)) {
  iso_id <- target_ids[idx]
  seq <- seq_vec[iso_id]
  seq_len <- nchar(seq)
  junctions <- junction_list[[iso_id]]

  if (is.null(junctions)) junctions <- integer(0)

  # Find all ATG positions (1-based)
  atg_matches <- gregexpr("ATG", seq, fixed = TRUE)[[1]]
  if (atg_matches[1] == -1L) {
    # No ATGs — rare but possible
    results[[idx]] <- data.frame(
      isoform_id = iso_id, n_orfs = 0L, n_nmd_competent = 0L,
      max_downstream_ejc = 0L, n_nmd_competent_strong_kozak = 0L,
      best_kozak_nmd = 0L, first_nmd_atg_pos = NA_integer_,
      has_any_nmd_competent = FALSE, n_atgs = 0L, tx_length = seq_len,
      n_junctions = length(junctions), n_orfs_nonstop = 0L,
      stringsAsFactors = FALSE
    )
    next
  }

  atg_positions <- as.integer(atg_matches)
  n_atgs <- length(atg_positions)

  # For each ATG, find first in-frame stop codon
  n_orfs <- 0L
  n_nmd_competent <- 0L
  n_nmd_strong <- 0L
  max_dejc <- 0L
  best_kozak_nmd <- 0L
  first_nmd_atg <- NA_integer_
  n_nonstop <- 0L

  for (atg_pos in atg_positions) {
    # Scan downstream in-frame (every 3 nt) for stop codon
    stop_pos <- NA_integer_
    pos <- atg_pos + 3L  # skip the ATG itself
    while (pos + 2L <= seq_len) {
      codon <- substr(seq, pos, pos + 2L)
      if (codon %in% stop_codons) {
        stop_pos <- pos
        break
      }
      pos <- pos + 3L
    }

    if (is.na(stop_pos)) {
      # No in-frame stop found — nonstop ORF, skip
      n_nonstop <- n_nonstop + 1L
      next
    }

    # We have an ORF: ATG at atg_pos, stop at stop_pos
    n_orfs <- n_orfs + 1L
    orf_length <- stop_pos - atg_pos  # in nt (includes ATG, excludes stop)

    # Count downstream EJCs: junctions after the stop codon position
    # The stop codon occupies positions stop_pos to stop_pos+2
    # EJCs downstream = junctions at positions > stop_pos + 2
    stop_end <- stop_pos + 2L
    n_dejc <- sum(junctions > stop_end)

    if (n_dejc > 0L) {
      n_nmd_competent <- n_nmd_competent + 1L
      if (n_dejc > max_dejc) max_dejc <- n_dejc

      # Kozak score for this ATG
      kozak <- score_kozak(seq, atg_pos)
      if (kozak == 2L) n_nmd_strong <- n_nmd_strong + 1L
      if (kozak > best_kozak_nmd) best_kozak_nmd <- kozak

      # Track first NMD-competent ATG position
      if (is.na(first_nmd_atg) || atg_pos < first_nmd_atg) {
        first_nmd_atg <- atg_pos
      }
    }
  }

  results[[idx]] <- data.frame(
    isoform_id = iso_id,
    n_orfs = n_orfs,
    n_nmd_competent = n_nmd_competent,
    max_downstream_ejc = max_dejc,
    n_nmd_competent_strong_kozak = n_nmd_strong,
    best_kozak_nmd = best_kozak_nmd,
    first_nmd_atg_pos = first_nmd_atg,
    has_any_nmd_competent = n_nmd_competent > 0L,
    n_atgs = n_atgs,
    tx_length = seq_len,
    n_junctions = length(junctions),
    n_orfs_nonstop = n_nonstop,
    stringsAsFactors = FALSE
  )

  if (idx %% 10000 == 0) {
    cat(sprintf("  Processed %d / %d isoforms (%.0f sec)\n",
                idx, length(target_ids), (proc.time() - t_scan)[3]))
  }
}

cat(sprintf("  ORF scanning completed in %.0f seconds\n", (proc.time() - t_scan)[3]))

# Combine results
orf_df <- do.call(rbind, results)

# Add derived features
orf_df$frac_nmd_competent <- orf_df$n_nmd_competent / pmax(orf_df$n_orfs, 1L)
orf_df$first_nmd_atg_frac <- orf_df$first_nmd_atg_pos / orf_df$tx_length
orf_df$n_nmd_competent_capped <- pmin(orf_df$n_nmd_competent, 20L)
orf_df$max_downstream_ejc_capped <- pmin(orf_df$max_downstream_ejc, 5L)

# Add NMD status
orf_df$is_nmd <- as.integer(orf_df$isoform_id %in% nmd_ids)

# ==============================================================================
# 5. Summary
# ==============================================================================
cat("\n=== SUMMARY ===\n\n")
cat("Total isoforms:", nrow(orf_df), "\n")
cat("  NMD:", sum(orf_df$is_nmd), "  Non-NMD:", sum(!orf_df$is_nmd), "\n\n")

cat("ORF landscape:\n")
cat("  Median ORFs per isoform:", median(orf_df$n_orfs), "\n")
cat("  Median ATGs:", median(orf_df$n_atgs), "\n")
cat("  Median NMD-competent ORFs:", median(orf_df$n_nmd_competent), "\n")
cat("  Isoforms with any NMD-competent ORF:",
    sum(orf_df$has_any_nmd_competent),
    sprintf("(%.1f%%)\n", 100 * mean(orf_df$has_any_nmd_competent)))

cat("\nBy NMD status:\n")
for (status in c(1, 0)) {
  sub <- orf_df[orf_df$is_nmd == status, ]
  label <- if (status == 1) "NMD" else "Non-NMD"
  cat(sprintf("  %s (n=%d): median ORFs=%d, median NMD-competent=%d, any NMD-competent=%.1f%%\n",
              label, nrow(sub), median(sub$n_orfs), median(sub$n_nmd_competent),
              100 * mean(sub$has_any_nmd_competent)))
}

# ==============================================================================
# 6. Save
# ==============================================================================
# Add chromosome for train/test split
chr_map <- structures[, c("isoform_id", "chr")]
orf_df <- merge(orf_df, chr_map, by = "isoform_id", all.x = TRUE)

saveRDS(orf_df, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
