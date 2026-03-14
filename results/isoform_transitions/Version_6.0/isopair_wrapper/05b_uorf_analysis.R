#!/usr/bin/env Rscript
# ==============================================================================
# 05b_uorf_analysis.R
#
# Scan 5'UTR sequences for upstream Open Reading Frames (uORFs) and compare
# prevalence across NMD PTC-negative, NMD PTC-positive, and Control (C4)
# comparator isoforms.
#
# uORFs in the 5'UTR are a known mechanism that can trigger NMD independently
# of premature termination codons in the main CDS. This script tests whether
# NMD-responsive isoforms that lack PTCs are enriched for uORFs.
#
# Input files (all relative to working directory):
#   - data_mashr/cds.rds              (CDS boundaries, genomic coords)
#   - data_mashr/structures.rds       (exon structures)
#   - data_mashr/ptc.rds              (PTC status per isoform)
#   - data_mashr/profiles_c2_allsamples.rds  (C2 NMD comparisons)
#   - data_mashr/profiles_c4_allsamples.rds  (C4 control comparisons)
#   - data_mashr/analysis_cache/ptc_c2_allsamples.rds  (PTC by comparison)
#   - FASTA: sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta
#
# Output:
#   - data_mashr/analysis_cache/uorf_analysis.rds
#
# Usage:
#   Rscript 05b_uorf_analysis.R
# ==============================================================================

suppressPackageStartupMessages({
  library(Biostrings)
})

# Set working directory
setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

FASTA_PATH <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta"
OUTPUT_PATH <- "data_mashr/analysis_cache/uorf_analysis.rds"

MIN_UORF_CODONS <- 3  # minimum uORF length in codons (including start codon)
MIN_UORF_NT <- MIN_UORF_CODONS * 3  # = 9 nt

cat("=== 5'UTR uORF Analysis ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Load data and define groups
# ==============================================================================
cat("Loading data files...\n")

cds <- readRDS("data_mashr/cds.rds")
structures <- readRDS("data_mashr/structures.rds")
ptc <- readRDS("data_mashr/ptc.rds")
profiles_c2 <- readRDS("data_mashr/profiles_c2_allsamples.rds")
profiles_c4 <- readRDS("data_mashr/profiles_c4_allsamples.rds")
ptc_c2 <- readRDS("data_mashr/analysis_cache/ptc_c2_allsamples.rds")

cat("  CDS:", nrow(cds), "isoforms\n")
cat("  Structures:", nrow(structures), "isoforms\n")
cat("  PTC:", nrow(ptc), "isoforms\n")
cat("  C2 profiles:", nrow(profiles_c2), "pairs\n")
cat("  C4 profiles:", nrow(profiles_c4), "pairs\n")

# --- Gene-matching: restrict C4 to genes/references shared with C2 ---
shared_keys <- paste(profiles_c2$gene_id, profiles_c2$reference_isoform_id)
profiles_c4_matched <- profiles_c4[
  paste(profiles_c4$gene_id, profiles_c4$reference_isoform_id) %in% shared_keys,
]
cat("  C4 after gene-matching:", nrow(profiles_c4_matched), "pairs\n")

# --- Coding filter: both ref and comparator must be coding ---
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]

c2_coding <- profiles_c2[
  profiles_c2$comparator_isoform_id %in% coding_ids &
  profiles_c2$reference_isoform_id %in% coding_ids,
]
c4_coding <- profiles_c4_matched[
  profiles_c4_matched$comparator_isoform_id %in% coding_ids &
  profiles_c4_matched$reference_isoform_id %in% coding_ids,
]

cat("  C2 after coding filter:", nrow(c2_coding), "pairs\n")
cat("  C4 after coding filter:", nrow(c4_coding), "pairs\n")

# --- Split C2 into PTC-positive and PTC-negative comparators ---
ptc_lookup <- ptc_c2$pair_summary[, c("comparator_isoform_id", "comp_has_ptc")]
c2_with_ptc <- merge(c2_coding, ptc_lookup, by = "comparator_isoform_id", all.x = TRUE)

# Isoforms not in ptc_c2 pair_summary shouldn't exist (both datasets are C2), but handle gracefully
if (any(is.na(c2_with_ptc$comp_has_ptc))) {
  cat("  WARNING:", sum(is.na(c2_with_ptc$comp_has_ptc)),
      "C2 comparators missing PTC status (dropped)\n")
  c2_with_ptc <- c2_with_ptc[!is.na(c2_with_ptc$comp_has_ptc), ]
}

ptc_pos_ids <- unique(c2_with_ptc$comparator_isoform_id[c2_with_ptc$comp_has_ptc == TRUE])
ptc_neg_ids <- unique(c2_with_ptc$comparator_isoform_id[c2_with_ptc$comp_has_ptc == FALSE])
control_ids <- unique(c4_coding$comparator_isoform_id)

cat("\n--- Group sizes ---\n")
cat("  NMD PTC-negative comparators:", length(ptc_neg_ids), "\n")
cat("  NMD PTC-positive comparators:", length(ptc_pos_ids), "\n")
cat("  Control (C4) comparators:    ", length(control_ids), "\n")

# Build lookup: isoform_id -> group
# An isoform could theoretically appear in multiple groups; assign priority:
# PTC_pos > PTC_neg > Control (shouldn't happen given C2 vs C4, but be safe)
all_ids <- unique(c(ptc_neg_ids, ptc_pos_ids, control_ids))
group_map <- rep(NA_character_, length(all_ids))
names(group_map) <- all_ids
group_map[all_ids %in% control_ids] <- "Control"
group_map[all_ids %in% ptc_neg_ids] <- "PTC_neg"
group_map[all_ids %in% ptc_pos_ids] <- "PTC_pos"

cat("  Total unique isoforms to analyze:", length(all_ids), "\n\n")

# ==============================================================================
# 2. Compute 5'UTR length (in transcript coordinates) for each isoform
# ==============================================================================
cat("Computing 5'UTR lengths from genomic CDS + exon structure...\n")

# Merge structures and CDS for our isoforms
# Coordinates are 1-based closed (verified: exon_end - exon_start + 1 = exon length,
# and sum of exon lengths matches FASTA sequence length)

struct_sub <- structures[structures$isoform_id %in% all_ids, ]
cds_sub <- cds[cds$isoform_id %in% all_ids & cds$coding_status == "coding", ]

# Merge
merged <- merge(
  struct_sub[, c("isoform_id", "strand", "exon_starts", "exon_ends")],
  cds_sub[, c("isoform_id", "cds_start", "cds_stop")],
  by = "isoform_id"
)

cat("  Merged structure+CDS for", nrow(merged), "isoforms\n")

# Function to compute 5'UTR length in transcript coordinates
# For + strand: translation start = cds_start (lowest genomic coord)
#   5'UTR = exonic bases with genomic position < cds_start
# For - strand: translation start = cds_stop (highest genomic coord = 5' end)
#   5'UTR = exonic bases with genomic position > cds_stop
#
# Coordinates are 1-based closed, so an exon from S to E has (E - S + 1) bases.

compute_utr5_length <- function(strand, exon_starts, exon_ends, cds_start, cds_stop) {
  utr5_len <- 0L

  if (strand == "+") {
    # Translation starts at cds_start
    tstart <- cds_start
    for (i in seq_along(exon_starts)) {
      es <- exon_starts[i]
      ee <- exon_ends[i]
      if (ee < tstart) {
        # Entire exon is in 5'UTR
        utr5_len <- utr5_len + (ee - es + 1L)
      } else if (es < tstart) {
        # Exon spans the translation start
        utr5_len <- utr5_len + (tstart - es)
        break
      } else {
        # Exon starts at or after translation start — no more 5'UTR
        break
      }
    }
  } else {
    # Minus strand: translation starts at cds_stop (highest genomic coord = 5')
    # Exons are stored low-to-high genomically, but the transcript reads high-to-low.
    # So we iterate from the LAST exon (highest genomic = 5' of transcript).
    tstart <- cds_stop
    for (i in rev(seq_along(exon_starts))) {
      es <- exon_starts[i]
      ee <- exon_ends[i]
      if (es > tstart) {
        # Entire exon is in 5'UTR (all positions > tstart)
        utr5_len <- utr5_len + (ee - es + 1L)
      } else if (ee > tstart) {
        # Exon spans the translation start
        utr5_len <- utr5_len + (ee - tstart)
        break
      } else {
        # Exon ends at or before translation start — no more 5'UTR
        break
      }
    }
  }

  return(utr5_len)
}

# Compute for all isoforms
utr5_lengths <- integer(nrow(merged))
for (i in seq_len(nrow(merged))) {
  utr5_lengths[i] <- compute_utr5_length(
    strand = merged$strand[i],
    exon_starts = merged$exon_starts[[i]],
    exon_ends = merged$exon_ends[[i]],
    cds_start = merged$cds_start[i],
    cds_stop = merged$cds_stop[i]
  )
}
merged$utr5_length <- utr5_lengths

cat("  5'UTR length distribution:\n")
cat("    Min:", min(merged$utr5_length),
    " Median:", median(merged$utr5_length),
    " Mean:", round(mean(merged$utr5_length), 1),
    " Max:", max(merged$utr5_length), "\n")
cat("    Zero-length 5'UTRs:", sum(merged$utr5_length == 0), "\n\n")

# ==============================================================================
# 3. Extract needed sequences from FASTA
# ==============================================================================
cat("Extracting sequences from FASTA for", length(all_ids), "isoforms...\n")

# Strategy: write needed IDs to a temp file, use awk to extract only matching
# sequences from the 1.7 GB FASTA. This is much faster than line-by-line R I/O.

id_file <- tempfile(fileext = ".txt")
writeLines(all_ids, id_file)

# awk script: read IDs into an associative array, then extract matching sequences
# Output format: ID<tab>SEQUENCE (one line per matched record)
awk_cmd <- sprintf(
  "awk 'BEGIN{while((getline id < \"%s\") > 0) want[id]=1} /^>/{if(seq && found) print hdr \"\\t\" seq; hdr=substr($1,2); found=(hdr in want); seq=\"\"; next} found{seq=seq$0} END{if(seq && found) print hdr \"\\t\" seq}' '%s'",
  id_file, FASTA_PATH
)

cat("  Running awk extraction...\n")
t_awk <- proc.time()
raw_output <- system(awk_cmd, intern = TRUE)
cat("  awk completed in", round((proc.time() - t_awk)[3], 1), "seconds\n")

unlink(id_file)

# Parse awk output into named vector
seq_vec <- character(0)
if (length(raw_output) > 0) {
  split_lines <- strsplit(raw_output, "\t", fixed = TRUE)
  seq_names <- vapply(split_lines, `[`, character(1), 1)
  seq_seqs <- vapply(split_lines, `[`, character(1), 2)
  seq_vec <- setNames(seq_seqs, seq_names)
}

cat("  Found", length(seq_vec), "of", length(all_ids), "sequences\n")

missing_seqs <- setdiff(all_ids, names(seq_vec))
if (length(missing_seqs) > 0) {
  cat("  WARNING:", length(missing_seqs), "isoforms not found in FASTA\n")
  cat("  First few:", head(missing_seqs, 5), "\n")
}

# ==============================================================================
# 4. Sanity checks on 5'UTR extraction
# ==============================================================================
cat("\nRunning sanity checks...\n")

# For coding transcripts, the first codon of the CDS should be ATG.
# The CDS starts right after the 5'UTR, so position utr5_length+1 in the
# transcript sequence should be 'A', utr5_length+2 = 'T', utr5_length+3 = 'G'.
n_checked <- 0L
n_atg_ok <- 0L
n_too_short <- 0L

for (i in seq_len(min(nrow(merged), nrow(merged)))) {
  iso <- merged$isoform_id[i]
  ulen <- merged$utr5_length[i]

  if (is.null(seq_vec[iso]) || is.na(seq_vec[iso])) next

  seq_i <- toupper(seq_vec[iso])

  # CDS start codon should be at position utr5_length + 1
  cds_codon_start <- ulen + 1L
  cds_codon_end <- ulen + 3L

  if (nchar(seq_i) < cds_codon_end) {
    n_too_short <- n_too_short + 1L
    next
  }

  start_codon <- substr(seq_i, cds_codon_start, cds_codon_end)
  n_checked <- n_checked + 1L
  if (start_codon == "ATG") n_atg_ok <- n_atg_ok + 1L
}

cat("  Checked", n_checked, "isoforms for ATG at CDS start\n")
cat("  ATG found:", n_atg_ok, "/", n_checked,
    "(", round(100 * n_atg_ok / n_checked, 1), "%)\n")
if (n_too_short > 0) cat("  Sequence too short for check:", n_too_short, "\n")

# Show a few examples where ATG is NOT found (for debugging)
if (n_atg_ok < n_checked) {
  cat("  Examples without ATG at CDS start:\n")
  n_shown <- 0L
  for (i in seq_len(nrow(merged))) {
    if (n_shown >= 3) break
    iso <- merged$isoform_id[i]
    ulen <- merged$utr5_length[i]
    if (is.null(seq_vec[iso]) || is.na(seq_vec[iso])) next
    seq_i <- toupper(seq_vec[iso])
    cds_start_pos <- ulen + 1L
    if (nchar(seq_i) < ulen + 3L) next
    codon <- substr(seq_i, cds_start_pos, cds_start_pos + 2L)
    if (codon != "ATG") {
      cat("    ", iso, ": 5'UTR=", ulen, "bp, codon at CDS start=", codon, "\n")
      n_shown <- n_shown + 1L
    }
  }
}

cat("\n")

# ==============================================================================
# 5. Scan 5'UTRs for uORFs
# ==============================================================================
cat("Scanning 5'UTRs for uORFs...\n")

# A uORF starts with AUG in the 5'UTR and ends at the first in-frame stop codon.
# We track:
#   - uORFs entirely within the 5'UTR
#   - Overlapping uORFs: AUG in 5'UTR but stop codon falls at or beyond main CDS start
# Minimum length: >= MIN_UORF_NT (9 nt, = 3 codons including ATG)

find_uorfs <- function(utr5_seq, full_seq, utr5_length) {
  # utr5_seq: the 5'UTR sequence (uppercase)
  # full_seq: the full transcript sequence (uppercase), needed for overlapping uORFs
  # utr5_length: length of 5'UTR in nt
  #
  # Returns a list with:
  #   n_uorfs: total count of qualifying uORFs
  #   n_overlapping: count of uORFs that extend into/past CDS
  #   max_length: max uORF length in nt (0 if none)
  #   total_length: sum of all uORF lengths in nt

  result <- list(
    n_uorfs = 0L,
    n_overlapping = 0L,
    max_length = 0L,
    total_length = 0L
  )

  if (utr5_length < 3L) return(result)

  full_len <- nchar(full_seq)

  # Use gregexpr to find all ATG positions in the 5'UTR (much faster than loop)
  atg_matches <- gregexpr("ATG", utr5_seq, fixed = TRUE)[[1]]
  if (atg_matches[1] == -1L) return(result)
  atg_positions <- as.integer(atg_matches)

  # For each ATG, walk in-frame through the full transcript using substr
  for (atg_pos in atg_positions) {
    # Walk in-frame from the codon after AUG through the full transcript
    found_stop <- FALSE
    stop_pos <- NA_integer_

    codon_start <- atg_pos + 3L  # skip the ATG itself
    while (codon_start + 2L <= full_len) {
      codon <- substr(full_seq, codon_start, codon_start + 2L)
      if (codon == "TAA" || codon == "TAG" || codon == "TGA") {
        found_stop <- TRUE
        stop_pos <- codon_start + 2L  # end of stop codon
        break
      }
      codon_start <- codon_start + 3L
    }

    if (!found_stop) next  # no in-frame stop found in transcript

    uorf_length <- stop_pos - atg_pos + 1L

    if (uorf_length < MIN_UORF_NT) next

    result$n_uorfs <- result$n_uorfs + 1L
    result$total_length <- result$total_length + uorf_length
    if (uorf_length > result$max_length) result$max_length <- uorf_length

    # Is the stop codon at or beyond the main CDS start?
    # CDS starts at position utr5_length + 1 in the transcript
    if (stop_pos > utr5_length) {
      result$n_overlapping <- result$n_overlapping + 1L
    }
  }

  return(result)
}

# Process all isoforms
results_list <- vector("list", nrow(merged))

t0 <- proc.time()
for (i in seq_len(nrow(merged))) {
  iso <- merged$isoform_id[i]
  ulen <- merged$utr5_length[i]

  if (i %% 1000 == 0) {
    elapsed <- (proc.time() - t0)[3]
    cat("  Processed", i, "of", nrow(merged),
        "(", round(elapsed, 1), "s elapsed)\n")
  }

  # Default result for missing sequences or zero-length UTR
  default_res <- list(
    isoform_id = iso,
    group = group_map[iso],
    utr5_length = ulen,
    n_uorfs = 0L,
    n_overlapping_uorfs = 0L,
    max_uorf_length_nt = 0L,
    total_uorf_length_nt = 0L
  )

  if (!(iso %in% names(seq_vec))) {
    results_list[[i]] <- default_res
    next
  }

  full_seq <- toupper(seq_vec[iso])

  if (ulen == 0L || ulen > nchar(full_seq)) {
    results_list[[i]] <- default_res
    next
  }

  utr5_seq <- substr(full_seq, 1L, ulen)

  uorf_res <- find_uorfs(utr5_seq, full_seq, ulen)

  results_list[[i]] <- list(
    isoform_id = iso,
    group = group_map[iso],
    utr5_length = ulen,
    n_uorfs = uorf_res$n_uorfs,
    n_overlapping_uorfs = uorf_res$n_overlapping,
    max_uorf_length_nt = uorf_res$max_length,
    total_uorf_length_nt = uorf_res$total_length
  )
}

elapsed_total <- (proc.time() - t0)[3]
cat("  Scanning complete in", round(elapsed_total, 1), "seconds\n\n")

# Convert to data frame
per_isoform <- do.call(rbind, lapply(results_list, as.data.frame,
                                      stringsAsFactors = FALSE))

# Drop any isoforms that had no group assignment (shouldn't happen)
per_isoform <- per_isoform[!is.na(per_isoform$group), ]

cat("=== Results Summary ===\n\n")

# Per-group summary
for (g in c("PTC_neg", "PTC_pos", "Control")) {
  sub <- per_isoform[per_isoform$group == g, ]
  cat("Group:", g, "(n =", nrow(sub), ")\n")
  cat("  5'UTR length: median =", median(sub$utr5_length),
      ", mean =", round(mean(sub$utr5_length), 1), "\n")
  cat("  Has any uORF:", sum(sub$n_uorfs > 0), "/", nrow(sub),
      "(", round(100 * mean(sub$n_uorfs > 0), 1), "%)\n")
  cat("  Has overlapping uORF:", sum(sub$n_overlapping_uorfs > 0), "/", nrow(sub),
      "(", round(100 * mean(sub$n_overlapping_uorfs > 0), 1), "%)\n")
  cat("  uORFs per isoform: median =", median(sub$n_uorfs),
      ", mean =", round(mean(sub$n_uorfs), 2), "\n")
  cat("\n")
}

# Statistical test: proportion with uORFs across groups
cat("--- Statistical Tests ---\n")
has_uorf <- per_isoform$n_uorfs > 0
group_factor <- factor(per_isoform$group, levels = c("Control", "PTC_neg", "PTC_pos"))

# Fisher's exact test: PTC_neg vs Control
tbl_neg <- table(
  group = per_isoform$group %in% "PTC_neg",
  uorf = has_uorf
)[c("TRUE", "FALSE"), ]
if (nrow(tbl_neg) == 2 && ncol(tbl_neg) == 2) {
  # Compare PTC_neg vs Control directly
  sub_neg <- per_isoform[per_isoform$group %in% c("PTC_neg", "Control"), ]
  tbl1 <- table(sub_neg$group, sub_neg$n_uorfs > 0)
  ft1 <- fisher.test(tbl1)
  cat("PTC_neg vs Control (any uORF): OR =", round(ft1$estimate, 3),
      ", p =", signif(ft1$p.value, 3), "\n")
}

# Fisher's exact test: PTC_pos vs Control
sub_pos <- per_isoform[per_isoform$group %in% c("PTC_pos", "Control"), ]
tbl2 <- table(sub_pos$group, sub_pos$n_uorfs > 0)
if (nrow(tbl2) == 2 && ncol(tbl2) == 2) {
  ft2 <- fisher.test(tbl2)
  cat("PTC_pos vs Control (any uORF): OR =", round(ft2$estimate, 3),
      ", p =", signif(ft2$p.value, 3), "\n")
}

# Fisher's exact test: PTC_neg vs PTC_pos
sub_np <- per_isoform[per_isoform$group %in% c("PTC_neg", "PTC_pos"), ]
tbl3 <- table(sub_np$group, sub_np$n_uorfs > 0)
if (nrow(tbl3) == 2 && ncol(tbl3) == 2) {
  ft3 <- fisher.test(tbl3)
  cat("PTC_neg vs PTC_pos (any uORF): OR =", round(ft3$estimate, 3),
      ", p =", signif(ft3$p.value, 3), "\n")
}

# Overlapping uORFs: PTC_neg vs Control
sub_neg2 <- per_isoform[per_isoform$group %in% c("PTC_neg", "Control"), ]
tbl4 <- table(sub_neg2$group, sub_neg2$n_overlapping_uorfs > 0)
if (nrow(tbl4) == 2 && ncol(tbl4) == 2) {
  ft4 <- fisher.test(tbl4)
  cat("PTC_neg vs Control (overlapping uORF): OR =", round(ft4$estimate, 3),
      ", p =", signif(ft4$p.value, 3), "\n")
}

cat("\n")

# ==============================================================================
# 6. Save results
# ==============================================================================
output <- list(
  per_isoform = per_isoform,
  metadata = list(
    run_timestamp = Sys.time(),
    group_sizes = table(per_isoform$group),
    script_path = "05b_uorf_analysis.R",
    min_uorf_codons = MIN_UORF_CODONS,
    fasta_path = FASTA_PATH,
    n_atg_check_passed = n_atg_ok,
    n_atg_check_total = n_checked
  )
)

saveRDS(output, OUTPUT_PATH)
cat("Saved results to:", OUTPUT_PATH, "\n")
cat("Finished:", format(Sys.time()), "\n")
