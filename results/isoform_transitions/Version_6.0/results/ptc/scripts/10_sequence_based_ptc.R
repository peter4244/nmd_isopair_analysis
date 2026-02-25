#!/usr/bin/env Rscript
# 10_sequence_based_ptc.R
#
# Sequence-based PTC detection in isoform pairs.
#
# For each isoform pair sharing a start codon, scans the actual transcript
# sequence to determine (1) whether it has a PTC via the 50nt EJC rule, and
# (2) which splicing event(s) cause it. Compares PTC rates between NMD
# (C1/C2) and baseline (C4) pairs.
#
# PTC definition: stop codon >50nt upstream of the last exon-exon junction
# (the canonical NMD "50-nucleotide rule").
#
# Uses dual FASTA sources:
#   - GENCODE: merged reference FASTA + genomic-to-mRNA coordinate mapping
#   - PacBio: SQANTI-corrected FASTA + CDS_start from classification file
#   (SQANTI correction trims PacBio sequences, so CDS coordinates from
#    sqanti3_classification.txt only match sqanti3_corrected.fasta)
#
# Requires:
#   - gencode49_merged_collapsed_2025.12.21.fasta + .fai (GENCODE seqs)
#   - sqanti3_corrected.fasta + .fai (PacBio corrected seqs)
#   - sqanti3_classification.txt (PacBio CDS coordinates in mRNA space)
#   - isoform_cds_metadata.rds + isoform_structures.rds
#   - splicing_choice_profiles.rds + pairs.tsv per comparison
#   - ptc_status.rds from 01_compute_ptc_status.R (for validation)
#   - Biostrings + Rsamtools (Bioconductor)
#
# Outputs (to results/ptc/results/nonNMD_{threshold}/):
#   sequence_ptc_pair_results.rds/.tsv
#   sequence_ptc_prevalence.tsv
#   sequence_ptc_causal_events.tsv
#   sequence_ptc_event_type_enrichment.tsv
#   sequence_ptc_validation.tsv
#   sequence_ptc_frame_shift_summary.tsv
#   sequence_ptc_within_pair_asymmetry.tsv
#
# Run from: results/isoform_transitions/Version_6.0/
# Usage: Rscript results/ptc/scripts/10_sequence_based_ptc.R [--threshold 0.50|0.95]

# Load Biostrings/Rsamtools FIRST so tidyverse/dplyr functions take precedence
# (Biostrings/IRanges mask rename, filter, select, etc.)
suppressPackageStartupMessages({
  library(Biostrings)
  library(Rsamtools)
  library(tidyverse)
})

# ─── Parse arguments ─────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
threshold <- "0.50"  # default
if ("--threshold" %in% args) {
  idx <- which(args == "--threshold")
  if (idx < length(args)) threshold <- args[idx + 1]
}
stopifnot(threshold %in% c("0.50", "0.95"))

cat("\n")
cat("==================================================================\n")
cat(sprintf("   10: Sequence-Based PTC Detection (nonNMD_%s)\n", threshold))
cat("==================================================================\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Load data and build FASTA lookups
# ═══════════════════════════════════════════════════════════════════════════════

cat("--- Section 1: Loading data ---\n\n")

# Paths — dual FASTA sources
gencode_fasta_path <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode49_merged_collapsed_2025.12.21.fasta"
sqanti_fasta_path  <- "/Users/petecastaldi/claude_projects/nmd/results/sqanti_runs/isoseq_sqanti3_filtered/sqanti3_corrected.fasta"
sqanti_cls_path    <- "/Users/petecastaldi/claude_projects/nmd/results/sqanti_runs/isoseq_sqanti3_filtered/sqanti3_classification.txt"
base_dir <- sprintf("comparisons/nonNMD_%s", threshold)

# Core data
cds_meta  <- readRDS("data/isoform_cds_metadata.rds")
structs   <- readRDS("data/isoform_structures.rds")
profiles  <- readRDS(file.path(base_dir, "deduplicated/splicing_choice_profiles.rds"))
ptc_status <- readRDS("results/ptc/results/ptc_status.rds")

cat(sprintf("  CDS metadata: %d isoforms\n", nrow(cds_meta)))
cat(sprintf("  Isoform structures: %d isoforms\n", nrow(structs)))
cat(sprintf("  Profiles: %d pairs\n", nrow(profiles)))

# Join CDS metadata with structures (need strand, gene_id, exon coords)
coding <- cds_meta %>%
  filter(coding_status == "coding") %>%
  inner_join(
    structs %>% select(isoform_id, gene_id, strand, seqnames,
                       exon_starts, exon_ends, n_exons),
    by = "isoform_id"
  ) %>%
  mutate(strand = as.character(strand))

cat(sprintf("  Coding isoforms with structures: %d\n", nrow(coding)))

# ─── GENCODE FASTA (merged reference, for ENST isoforms) ────────────────────
cat("  Building GENCODE FASTA index...\n")
gc_fai <- read_tsv(paste0(gencode_fasta_path, ".fai"),
                   col_names = c("name", "length", "offset", "linebases", "linewidth"),
                   show_col_types = FALSE)
gc_id_to_name   <- setNames(gc_fai$name, sub("[|].*", "", gc_fai$name))
gc_id_to_length <- setNames(gc_fai$length, sub("[|].*", "", gc_fai$name))
gc_fa <- FaFile(gencode_fasta_path)
cat(sprintf("  GENCODE FASTA index: %d sequences\n", nrow(gc_fai)))

# ─── SQANTI FASTA (corrected PacBio sequences) ──────────────────────────────
cat("  Building SQANTI FASTA index...\n")
sq_fai <- read_tsv(paste0(sqanti_fasta_path, ".fai"),
                   col_names = c("name", "length", "offset", "linebases", "linewidth"),
                   show_col_types = FALSE)
sq_id_to_name   <- setNames(sq_fai$name, sq_fai$name)  # PB headers are just PB.X.Y
sq_id_to_length <- setNames(sq_fai$length, sq_fai$name)
sq_fa <- FaFile(sqanti_fasta_path)
cat(sprintf("  SQANTI FASTA index: %d sequences\n", nrow(sq_fai)))

# ─── SQANTI classification (PacBio CDS coordinates in mRNA space) ───────────
cat("  Loading SQANTI classification (CDS columns only)...\n")
sqanti_cls <- read_tsv(sqanti_cls_path, col_types = cols_only(
  isoform = col_character(), coding = col_character(),
  CDS_start = col_integer(), CDS_end = col_integer(),
  CDS_length = col_integer(), ORF_length = col_integer(),
  predicted_NMD = col_character()
))
# Filter to our PacBio coding isoforms
pb_coding_ids <- cds_meta$isoform_id[cds_meta$coding_status == "coding" &
                                       startsWith(cds_meta$isoform_id, "PB.")]
sqanti_cds <- sqanti_cls %>%
  filter(isoform %in% pb_coding_ids, coding == "coding") %>%
  rename(isoform_id = isoform)
cat(sprintf("  SQANTI CDS info for PacBio coding isoforms: %d\n", nrow(sqanti_cds)))

# Helper: extract sequences using appropriate FASTA source
extract_sequences <- function(iso_ids) {
  gencode_ids <- iso_ids[startsWith(iso_ids, "ENST")]
  pacbio_ids  <- iso_ids[startsWith(iso_ids, "PB.")]

  result <- setNames(character(0), character(0))

  # GENCODE from merged reference FASTA
  if (length(gencode_ids) > 0) {
    gc_headers <- gc_id_to_name[gencode_ids]
    gc_lengths <- gc_id_to_length[gencode_ids]
    gc_valid <- !is.na(gc_headers)
    if (sum(gc_valid) > 0) {
      gr <- GRanges(seqnames = gc_headers[gc_valid],
                    ranges = IRanges(start = 1, width = gc_lengths[gc_valid]))
      seqs <- scanFa(gc_fa, gr)
      result <- c(result, setNames(as.character(seqs), gencode_ids[gc_valid]))
    }
  }

  # PacBio from SQANTI-corrected FASTA
  if (length(pacbio_ids) > 0) {
    sq_headers <- sq_id_to_name[pacbio_ids]
    sq_lengths <- sq_id_to_length[pacbio_ids]
    sq_valid <- !is.na(sq_headers)
    if (sum(sq_valid) > 0) {
      gr <- GRanges(seqnames = sq_headers[sq_valid],
                    ranges = IRanges(start = 1, width = sq_lengths[sq_valid]))
      seqs <- scanFa(sq_fa, gr)
      result <- c(result, setNames(as.character(seqs), pacbio_ids[sq_valid]))
    }
  }

  result
}

# Load comparison pairs
comparisons <- c("C1", "C2", "C4")
runs <- c("all_samples", "at", "dd", "dd_ali", "fb", "mv")

all_pairs <- list()
for (comp in comparisons) {
  for (run in runs) {
    path <- file.path(base_dir, comp, run, "pairs.tsv")
    if (!file.exists(path)) next
    all_pairs[[paste(comp, run)]] <- read_tsv(path, show_col_types = FALSE) %>%
      mutate(comparison = comp, run = run) %>%
      rename(non_dominant_isoform_id = comparator_isoform_id)
  }
}
pair_map <- bind_rows(all_pairs)
cat(sprintf("  Total pair-run combinations: %d\n", nrow(pair_map)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Identify same-start-codon pairs
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 2: Identifying same-start-codon pairs ---\n\n")

# Compute start codon genomic position (strand-aware)
coding_starts <- coding %>%
  mutate(
    start_codon_genomic = ifelse(strand == "+", cds_start, cds_stop)
  ) %>%
  select(isoform_id, gene_id, strand, start_codon_genomic,
         cds_start, cds_stop, orf_length)

# Filter profiles to both-coding pairs
profile_pairs <- profiles %>%
  select(gene_id, dominant_isoform_id, non_dominant_isoform_id)

# Annotate dominant and comparator with coding status + start codon
profile_annotated <- profile_pairs %>%
  inner_join(
    coding_starts %>% select(isoform_id, start_codon_genomic) %>%
      rename(dom_start_codon = start_codon_genomic),
    by = c("dominant_isoform_id" = "isoform_id")
  ) %>%
  inner_join(
    coding_starts %>% select(isoform_id, start_codon_genomic) %>%
      rename(comp_start_codon = start_codon_genomic),
    by = c("non_dominant_isoform_id" = "isoform_id")
  )

n_both_coding <- nrow(profile_annotated)
cat(sprintf("  Both-coding pairs: %d / %d\n", n_both_coding, nrow(profiles)))

# Same start codon
same_start <- profile_annotated %>%
  filter(dom_start_codon == comp_start_codon)

cat(sprintf("  Same start codon: %d / %d both-coding (%.1f%%)\n",
            nrow(same_start), n_both_coding,
            100 * nrow(same_start) / n_both_coding))

# Check which have FASTA sequences (dual source: GENCODE from merged, PacBio from SQANTI)
fasta_available_ids <- c(
  names(gc_id_to_name)[startsWith(names(gc_id_to_name), "ENST")],
  names(sq_id_to_name)[startsWith(names(sq_id_to_name), "PB.")]
)
same_start_with_seq <- same_start %>%
  filter(dominant_isoform_id %in% fasta_available_ids,
         non_dominant_isoform_id %in% fasta_available_ids)

cat(sprintf("  With FASTA sequences for both: %d / %d\n",
            nrow(same_start_with_seq), nrow(same_start)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Map start codon to transcript coordinates and scan ORFs
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 3: Mapping start codons and scanning ORFs ---\n\n")

# Get unique isoforms to scan
isoforms_to_scan <- unique(c(same_start_with_seq$dominant_isoform_id,
                              same_start_with_seq$non_dominant_isoform_id))
cat(sprintf("  Unique isoforms to scan: %d\n", length(isoforms_to_scan)))

gencode_to_scan <- isoforms_to_scan[startsWith(isoforms_to_scan, "ENST")]
pacbio_to_scan  <- isoforms_to_scan[startsWith(isoforms_to_scan, "PB.")]
cat(sprintf("    GENCODE: %d, PacBio: %d\n",
            length(gencode_to_scan), length(pacbio_to_scan)))

# Prepare isoform data for scanning
scan_data <- coding %>%
  filter(isoform_id %in% isoforms_to_scan) %>%
  mutate(
    source = ifelse(startsWith(isoform_id, "ENST"), "GENCODE", "PacBio"),
    start_codon_genomic = ifelse(strand == "+", cds_start, cds_stop),
    stop_codon_genomic = ifelse(strand == "+", cds_stop, cds_start)
  )

cat(sprintf("  Isoforms with structure data: %d\n", nrow(scan_data)))

# Function: map genomic position to mRNA coordinate (0-based from 5' end)
# Exon coordinates are 1-based inclusive (GFF-style): length = end - start + 1
genomic_to_mrna <- function(genomic_pos, strand, exon_starts, exon_ends) {
  # Order exons 5' → 3'
  if (strand == "+") {
    ord <- order(exon_starts)
  } else {
    ord <- order(exon_starts, decreasing = TRUE)
  }
  e_starts <- exon_starts[ord]
  e_ends   <- exon_ends[ord]
  exon_lengths <- e_ends - e_starts + 1L  # 1-based inclusive
  cum_starts   <- c(0L, cumsum(exon_lengths[-length(exon_lengths)]))

  # Find exon containing the position
  exon_idx <- which(e_starts <= genomic_pos & genomic_pos <= e_ends)
  if (length(exon_idx) == 0) return(NA_integer_)
  exon_idx <- exon_idx[1]

  if (strand == "+") {
    offset <- genomic_pos - e_starts[exon_idx]
  } else {
    offset <- e_ends[exon_idx] - genomic_pos
  }
  as.integer(cum_starts[exon_idx] + offset)  # 0-based
}

# Function: scan ORF from a given start position in transcript sequence
# start_pos_1based: 1-based position of 'A' in ATG
# Returns: list(valid_atg, orf_length_nt, stop_codon_pos_1based, stop_codon_type, reached_end)
scan_orf <- function(seq_str, start_pos_1based) {
  seq_len <- nchar(seq_str)

  # Validate ATG
  if (start_pos_1based + 2 > seq_len) {
    return(list(valid_atg = FALSE, orf_length_nt = NA_integer_,
                stop_codon_pos_1based = NA_integer_, stop_codon_type = NA_character_,
                reached_end = NA))
  }
  first_codon <- substr(seq_str, start_pos_1based, start_pos_1based + 2)
  valid_atg <- first_codon == "ATG"

  # Scan codons from start
  pos <- start_pos_1based + 3  # skip ATG, start at next codon
  stop_codons <- c("TAA", "TAG", "TGA")

  while (pos + 2 <= seq_len) {
    codon <- substr(seq_str, pos, pos + 2)
    if (codon %in% stop_codons) {
      orf_len <- as.integer(pos - start_pos_1based)  # ATG to stop (exclusive of stop)
      return(list(valid_atg = valid_atg, orf_length_nt = orf_len,
                  stop_codon_pos_1based = as.integer(pos), stop_codon_type = codon,
                  reached_end = FALSE))
    }
    pos <- pos + 3
  }

  # No stop codon found — ORF runs to end of transcript
  orf_len <- as.integer(seq_len - start_pos_1based + 1)
  return(list(valid_atg = valid_atg, orf_length_nt = orf_len,
              stop_codon_pos_1based = NA_integer_, stop_codon_type = NA_character_,
              reached_end = TRUE))
}

# Function: compute last EJC position in mRNA space (1-based)
# Last EJC = junction between penultimate and last exon
# For GENCODE: exon structure lengths == FASTA length, so direct computation works.
# For PacBio: SQANTI correction trims 5'/3' ends, so exon total != FASTA length.
#   We use: last_ejc = actual_seq_length - last_exon_length_from_structure
#   This is accurate because SQANTI preserves internal splice sites; only terminal
#   exons are trimmed. The last exon length from genomic coords may overestimate
#   slightly if 3' polyA is trimmed, making this a conservative PTC estimate.
compute_last_ejc_1based <- function(strand, exon_starts, exon_ends, seq_length = NULL) {
  n <- length(exon_starts)
  if (n <= 1) return(NA_integer_)  # single-exon, no EJCs
  if (strand == "+") { ord <- order(exon_starts) }
  else { ord <- order(exon_starts, decreasing = TRUE) }
  e_starts <- exon_starts[ord]
  e_ends <- exon_ends[ord]
  exon_lengths <- e_ends - e_starts + 1L

  if (!is.null(seq_length)) {
    # Use actual FASTA length - handles SQANTI trimming for PacBio
    last_exon_len <- exon_lengths[n]
    return(as.integer(seq_length - last_exon_len))
  }

  # Default: compute from exon lengths (for GENCODE where exon total == FASTA length)
  as.integer(sum(exon_lengths[-n]))
}

# ─── Compute start positions: dual-source strategy ──────────────────────────
cat("  Computing mRNA start positions...\n")

# GENCODE: use genomic-to-mRNA coordinate mapping (100% accurate)
gencode_scan <- scan_data %>%
  filter(source == "GENCODE") %>%
  rowwise() %>%
  mutate(
    start_pos_1based = genomic_to_mrna(
      start_codon_genomic, strand, exon_starts, exon_ends
    ) + 1L  # convert 0-based to 1-based
  ) %>%
  ungroup()

cat(sprintf("    GENCODE: %d isoforms, %d with valid start mapping\n",
            nrow(gencode_scan), sum(!is.na(gencode_scan$start_pos_1based))))

# PacBio: use CDS_start from SQANTI classification (1-based in SQANTI FASTA)
pacbio_scan <- scan_data %>%
  filter(source == "PacBio") %>%
  left_join(sqanti_cds %>% select(isoform_id, CDS_start), by = "isoform_id") %>%
  mutate(start_pos_1based = CDS_start)

cat(sprintf("    PacBio: %d isoforms, %d with SQANTI CDS_start\n",
            nrow(pacbio_scan), sum(!is.na(pacbio_scan$start_pos_1based))))

# Combine and filter
scan_data <- bind_rows(
  gencode_scan %>% select(-any_of("CDS_start")),
  pacbio_scan %>% select(-any_of("CDS_start"))
) %>%
  filter(!is.na(start_pos_1based))

cat(sprintf("  Total isoforms with valid start positions: %d\n", nrow(scan_data)))

# ─── Extract sequences and scan ORFs ────────────────────────────────────────
cat("  Extracting sequences from dual FASTA sources...\n")
seq_map <- extract_sequences(scan_data$isoform_id)
cat(sprintf("  Extracted %d sequences\n", length(seq_map)))

scan_data <- scan_data %>%
  mutate(seq_str = seq_map[isoform_id]) %>%
  filter(!is.na(seq_str))

cat(sprintf("  Isoforms with sequences: %d\n", nrow(scan_data)))

cat("  Scanning ORFs...\n")
orf_results <- scan_data %>%
  rowwise() %>%
  mutate(
    orf_scan = list(scan_orf(seq_str, start_pos_1based)),
    valid_atg = orf_scan$valid_atg,
    scanned_orf_length_nt = orf_scan$orf_length_nt,
    stop_codon_pos_1based = orf_scan$stop_codon_pos_1based,
    stop_codon_type = orf_scan$stop_codon_type,
    reached_end = orf_scan$reached_end,
    # Compute last EJC for 50nt rule (pass seq length for PacBio coordinate adjustment)
    last_ejc_1based = compute_last_ejc_1based(
      strand, exon_starts, exon_ends,
      seq_length = if (source == "PacBio") nchar(seq_str) else NULL
    ),
    # 50nt EJC rule: stop codon >50nt upstream of last EJC
    # ptc_distance = last_ejc - stop_codon_pos (in mRNA 1-based coords)
    ptc_distance = ifelse(!is.na(stop_codon_pos_1based) & !is.na(last_ejc_1based),
                          last_ejc_1based - stop_codon_pos_1based, NA_integer_),
    has_ptc = case_when(
      n_exons == 1 ~ FALSE,                 # single-exon: no EJCs
      is.na(stop_codon_pos_1based) ~ NA,    # no stop found
      is.na(last_ejc_1based) ~ NA,          # no EJC info
      ptc_distance > 50 ~ TRUE,             # PTC: stop >50nt before last EJC
      TRUE ~ FALSE
    )
  ) %>%
  ungroup() %>%
  select(isoform_id, gene_id, strand, source,
         start_codon_genomic, stop_codon_genomic,
         start_pos_1based, n_exons, valid_atg, scanned_orf_length_nt,
         stop_codon_pos_1based, stop_codon_type, reached_end,
         last_ejc_1based, ptc_distance, has_ptc, orf_length)

cat(sprintf("  Scanned %d isoforms\n", nrow(orf_results)))
cat(sprintf("  Valid ATG: %d (%.1f%%)\n",
            sum(orf_results$valid_atg), 100 * mean(orf_results$valid_atg)))
cat(sprintf("    GENCODE: %d / %d (%.1f%%)\n",
            sum(orf_results$valid_atg & orf_results$source == "GENCODE"),
            sum(orf_results$source == "GENCODE"),
            100 * mean(orf_results$valid_atg[orf_results$source == "GENCODE"])))
cat(sprintf("    PacBio:  %d / %d (%.1f%%)\n",
            sum(orf_results$valid_atg & orf_results$source == "PacBio"),
            sum(orf_results$source == "PacBio"),
            100 * mean(orf_results$valid_atg[orf_results$source == "PacBio"])))
cat(sprintf("  Reached end (no stop): %d\n",
            sum(orf_results$reached_end, na.rm = TRUE)))
cat(sprintf("  Has PTC (50nt EJC rule): %d (%.1f%% of valid ATG multi-exon)\n",
            sum(orf_results$has_ptc == TRUE, na.rm = TRUE),
            100 * mean(orf_results$has_ptc[orf_results$valid_atg &
                                            orf_results$n_exons > 1], na.rm = TRUE)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Build pair-level PTC results
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 4: Building pair-level PTC results ---\n\n")

# Create lookup for ORF results — now includes 50nt EJC rule PTC status per isoform
orf_lookup <- orf_results %>%
  select(isoform_id, valid_atg, scanned_orf_length_nt,
         stop_codon_pos_1based, stop_codon_type, reached_end,
         last_ejc_1based, ptc_distance, has_ptc)

# Annotate pairs with ORF scan results for both sides
pair_results <- same_start_with_seq %>%
  inner_join(
    orf_lookup %>% rename_with(~paste0("dom_", .), -isoform_id),
    by = c("dominant_isoform_id" = "isoform_id")
  ) %>%
  inner_join(
    orf_lookup %>% rename_with(~paste0("comp_", .), -isoform_id),
    by = c("non_dominant_isoform_id" = "isoform_id")
  )

# PTC status is already computed per-isoform using the 50nt EJC rule.
# For pair-level: rename for clarity
pair_results <- pair_results %>%
  mutate(
    both_valid_atg = dom_valid_atg & comp_valid_atg,
    comp_has_ptc = ifelse(both_valid_atg & comp_valid_atg, comp_has_ptc, NA),
    dom_has_ptc  = ifelse(both_valid_atg & dom_valid_atg, dom_has_ptc, NA),
    orf_length_diff = comp_scanned_orf_length_nt - dom_scanned_orf_length_nt
  )

n_valid <- sum(pair_results$both_valid_atg, na.rm = TRUE)
n_comp_ptc <- sum(pair_results$comp_has_ptc == TRUE, na.rm = TRUE)
n_dom_ptc  <- sum(pair_results$dom_has_ptc == TRUE, na.rm = TRUE)
cat(sprintf("  Pairs with both valid ATG: %d / %d\n", n_valid, nrow(pair_results)))
cat(sprintf("  Comparator has PTC (50nt EJC rule): %d (%.1f%%)\n",
            n_comp_ptc, 100 * n_comp_ptc / n_valid))
cat(sprintf("  Dominant has PTC (50nt EJC rule): %d (%.1f%%)\n",
            n_dom_ptc, 100 * n_dom_ptc / n_valid))
cat(sprintf("  Both have PTC: %d\n",
            sum(pair_results$comp_has_ptc == TRUE & pair_results$dom_has_ptc == TRUE,
                na.rm = TRUE)))
cat(sprintf("  Neither has PTC: %d\n",
            sum(pair_results$comp_has_ptc == FALSE & pair_results$dom_has_ptc == FALSE,
                na.rm = TRUE)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5: Validation against existing annotations
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 5: Validation against annotations ---\n\n")

# ─── 5a. ATG validation by source ───────────────────────────────────────────
cat("ATG validation by source:\n")
orf_results %>%
  group_by(source) %>%
  summarise(
    n = n(),
    n_valid_atg = sum(valid_atg),
    atg_pct = round(100 * mean(valid_atg), 1),
    .groups = "drop"
  ) %>%
  print()

# ─── 5b. ORF length: scanned vs annotated ──────────────────────────────────
cat("\nORF length validation (scanned vs annotated orf_length from CDS metadata):\n")
validation <- orf_results %>%
  filter(valid_atg, !reached_end) %>%
  mutate(
    orf_diff = scanned_orf_length_nt - orf_length,
    match_exact = orf_diff == 0,
    match_close = abs(orf_diff) <= 3
  )

validation %>%
  group_by(source) %>%
  summarise(
    n = n(),
    exact_match = sum(match_exact),
    close_match = sum(match_close),
    exact_pct = round(100 * mean(match_exact), 1),
    close_pct = round(100 * mean(match_close), 1),
    median_diff = median(orf_diff),
    .groups = "drop"
  ) %>%
  print()

# ─── 5c. PTC concordance: sequence-based (50nt rule) vs Script 01 ──────────
cat("\nPTC concordance: sequence-based (50nt rule) vs annotation-based (Script 01):\n")

ptc_comparison <- orf_results %>%
  filter(valid_atg, !is.na(has_ptc)) %>%
  inner_join(
    ptc_status %>% select(isoform_id, has_ptc, ptc_distance) %>%
      rename(annot_has_ptc = has_ptc, annot_ptc_distance = ptc_distance),
    by = "isoform_id"
  )

if (nrow(ptc_comparison) > 0) {
  cat(sprintf("  Isoforms with both calls: %d\n", nrow(ptc_comparison)))

  # Concordance table
  conc_tab <- ptc_comparison %>%
    filter(!is.na(annot_has_ptc)) %>%
    count(source, sequence_ptc = has_ptc, annotation_ptc = annot_has_ptc)
  print(conc_tab, n = 20)

  # Overall concordance
  agree <- ptc_comparison %>%
    filter(!is.na(annot_has_ptc)) %>%
    summarise(
      n = n(),
      concordant = sum(has_ptc == annot_has_ptc),
      concordance_pct = round(100 * mean(has_ptc == annot_has_ptc), 1)
    )
  cat(sprintf("  Overall concordance: %d / %d (%.1f%%)\n",
              agree$concordant, agree$n, agree$concordance_pct))

  # By source
  ptc_comparison %>%
    filter(!is.na(annot_has_ptc)) %>%
    group_by(source) %>%
    summarise(
      n = n(),
      concordant = sum(has_ptc == annot_has_ptc),
      concordance_pct = round(100 * mean(has_ptc == annot_has_ptc), 1),
      .groups = "drop"
    ) %>%
    print()
}

# Build validation output for saving
validation_out <- orf_results %>%
  filter(valid_atg) %>%
  group_by(source) %>%
  summarise(
    n_isoforms = n(),
    n_valid_atg = sum(valid_atg),
    atg_pct = round(100 * mean(valid_atg), 1),
    n_with_stop = sum(!reached_end, na.rm = TRUE),
    n_has_ptc = sum(has_ptc == TRUE, na.rm = TRUE),
    ptc_rate = round(100 * mean(has_ptc[!is.na(has_ptc)], na.rm = TRUE), 1),
    median_ptc_distance = median(ptc_distance[!is.na(ptc_distance)]),
    .groups = "drop"
  )
print(validation_out)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6: Attribute PTCs to causal splicing events
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 6: Attributing PTCs to causal events ---\n\n")

# Get pairs where comparator has PTC
ptc_pairs <- pair_results %>%
  filter(comp_has_ptc == TRUE)
cat(sprintf("  PTC-positive pairs: %d\n", nrow(ptc_pairs)))

# Unnest events for these pairs
events <- profiles %>%
  semi_join(ptc_pairs, by = c("gene_id", "dominant_isoform_id",
                               "non_dominant_isoform_id")) %>%
  select(gene_id, dominant_isoform_id, non_dominant_isoform_id,
         detailed_events) %>%
  unnest(detailed_events, names_sep = "_") %>%
  select(
    gene_id, dominant_isoform_id, non_dominant_isoform_id,
    event_type = detailed_events_event_type,
    direction = detailed_events_direction,
    five_prime = detailed_events_five_prime,
    three_prime = detailed_events_three_prime,
    bp_diff = detailed_events_bp_diff,
    strand = detailed_events_strand
  ) %>%
  mutate(strand = as.character(strand))

cat(sprintf("  Events in PTC pairs: %d\n", nrow(events)))

# Filter to CDS-region events only.
# Since both isoforms share the same start codon, events upstream (5'UTR)
# cannot affect the reading frame. Only events between the start codon and
# the dominant's stop codon can cause frame disruption or introduce PTCs.
cds_lookup <- coding_starts %>%
  select(isoform_id, cds_start, cds_stop)

events <- events %>%
  left_join(cds_lookup, by = c("dominant_isoform_id" = "isoform_id")) %>%
  mutate(
    event_start_genomic = pmin(five_prime, three_prime),
    event_end_genomic = pmax(five_prime, three_prime),
    # CDS region: between cds_start and cds_stop (genomic envelope)
    overlaps_cds = event_end_genomic > cds_start & event_start_genomic < cds_stop
  )

n_cds <- sum(events$overlaps_cds, na.rm = TRUE)
cat(sprintf("  Events overlapping CDS: %d / %d (%.1f%%)\n",
            n_cds, nrow(events), 100 * n_cds / nrow(events)))

events_cds <- events %>% filter(overlaps_cds)

# Function: determine if event could directly insert a stop codon
# IR/Partial_IR events = gained intron sequence → likely contains stops
direct_stop_event_types <- c("IR", "IR_diff_5", "IR_diff_3", "IR_diff_5_3",
                              "Partial_IR_5", "Partial_IR_3")

# For each PTC pair, walk through events 5'→3' and track cumulative frame shift
cat("  Attributing causal events...\n")

# Pre-compute: genomic midpoint for ordering
events_positioned <- events_cds %>%
  mutate(
    event_mid_genomic = (event_start_genomic + event_end_genomic) / 2
  )

# For each pair, order events 5'→3' and track cumulative shift
attribute_causal_events <- function(pair_events, pair_strand) {
  # Order events 5'→3' based on strand
  if (pair_strand == "+") {
    pair_events <- pair_events %>% arrange(event_mid_genomic)
  } else {
    pair_events <- pair_events %>% arrange(desc(event_mid_genomic))
  }

  n_ev <- nrow(pair_events)
  cumulative_shift <- integer(n_ev)
  is_causal <- logical(n_ev)
  is_contributing <- logical(n_ev)
  cause_type <- character(n_ev)

  cum <- 0L
  causal_found <- FALSE

  for (i in seq_len(n_ev)) {
    ev_type <- pair_events$event_type[i]
    ev_dir  <- pair_events$direction[i]
    ev_bp   <- pair_events$bp_diff[i]

    # Compute bp change to comparator relative to dominant:
    # LOSS = comparator lost sequence (dominant has MORE)
    #   → comparator is SHORTER by bp_diff
    # GAIN = comparator gained sequence (dominant has LESS)
    #   → comparator is LONGER by bp_diff
    if (is.na(ev_bp)) {
      # Skip events with NA bp_diff (e.g., Alt_TSS/Alt_TES with no size)
      cumulative_shift[i] <- cum
      next
    }

    if (ev_dir == "LOSS") {
      delta <- -abs(ev_bp)
    } else {
      delta <- abs(ev_bp)
    }

    cum <- cum + delta
    cumulative_shift[i] <- cum

    if (!causal_found) {
      is_contributing[i] <- TRUE

      # Check for direct stop insertion (IR/retention events that GAIN sequence)
      if (ev_type %in% direct_stop_event_types && ev_dir == "GAIN") {
        is_causal[i] <- TRUE
        cause_type[i] <- "direct_stop_insertion"
        causal_found <- TRUE
      }
      # Check for frame disruption
      else if (cum %% 3 != 0) {
        is_causal[i] <- TRUE
        cause_type[i] <- "frame_disruption"
        causal_found <- TRUE
      }
    }
  }

  pair_events %>%
    mutate(
      cumulative_shift = cumulative_shift,
      is_causal = is_causal,
      is_contributing = is_contributing,
      cause_type = ifelse(is_causal, cause_type, NA_character_),
      net_shift = cum,
      net_frame_status = ifelse(cum %% 3 == 0, "frame_preserved", "frame_disrupted")
    )
}

# Apply to each pair
causal_results <- events_positioned %>%
  group_by(gene_id, dominant_isoform_id, non_dominant_isoform_id) %>%
  group_modify(~{
    pair_strand <- .x$strand[1]
    attribute_causal_events(.x, pair_strand)
  }) %>%
  ungroup()

cat(sprintf("  Events attributed: %d\n", nrow(causal_results)))
cat(sprintf("  Causal events identified: %d\n", sum(causal_results$is_causal)))

# Summarize causal event types
cat("\nCausal event types:\n")
causal_results %>%
  filter(is_causal) %>%
  count(event_type, cause_type) %>%
  arrange(desc(n)) %>%
  print(n = 20)

# Also process non-PTC pairs for the net shift analysis
cat("\n  Computing net frame shift for all same-start-codon pairs...\n")

# Unnest events for ALL same-start-codon pairs (not just PTC)
all_events <- profiles %>%
  semi_join(pair_results %>% filter(both_valid_atg),
            by = c("gene_id", "dominant_isoform_id", "non_dominant_isoform_id")) %>%
  select(gene_id, dominant_isoform_id, non_dominant_isoform_id,
         detailed_events) %>%
  unnest(detailed_events, names_sep = "_") %>%
  select(
    gene_id, dominant_isoform_id, non_dominant_isoform_id,
    event_type = detailed_events_event_type,
    direction = detailed_events_direction,
    five_prime = detailed_events_five_prime,
    three_prime = detailed_events_three_prime,
    bp_diff = detailed_events_bp_diff
  )

# Filter to CDS-region events only (same logic as Section 6)
all_events <- all_events %>%
  left_join(cds_lookup, by = c("dominant_isoform_id" = "isoform_id")) %>%
  mutate(
    event_start_genomic = pmin(five_prime, three_prime),
    event_end_genomic = pmax(five_prime, three_prime),
    overlaps_cds = event_end_genomic > cds_start & event_start_genomic < cds_stop
  ) %>%
  filter(overlaps_cds)

# Compute net bp change per pair (CDS-region events only)
net_shift_per_pair <- all_events %>%
  mutate(
    delta = ifelse(direction == "LOSS", -abs(bp_diff), abs(bp_diff))
  ) %>%
  group_by(gene_id, dominant_isoform_id, non_dominant_isoform_id) %>%
  summarise(
    n_cds_events = n(),
    net_shift = sum(delta, na.rm = TRUE),
    net_frame_status = ifelse(sum(delta, na.rm = TRUE) %% 3 == 0,
                               "frame_preserved", "frame_disrupted"),
    .groups = "drop"
  )

# Attach to pair results
pair_results <- pair_results %>%
  left_join(net_shift_per_pair, by = c("gene_id", "dominant_isoform_id",
                                        "non_dominant_isoform_id"))

cat(sprintf("  Pairs with CDS events and net shift data: %d\n",
            sum(!is.na(pair_results$net_shift))))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7: Statistical comparisons
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 7: Statistical comparisons ---\n\n")

# Map pairs to comparisons
pair_comparison <- pair_results %>%
  filter(both_valid_atg) %>%
  inner_join(
    pair_map %>% distinct(comparison, run, gene_id, dominant_isoform_id,
                           non_dominant_isoform_id),
    by = c("gene_id", "dominant_isoform_id", "non_dominant_isoform_id")
  )

cat(sprintf("  Pairs mapped to comparisons: %d\n", nrow(pair_comparison)))

# ─── 7a. PTC prevalence: NMD vs C4 baseline ─────────────────────────────────

cat("\n--- 7a: PTC prevalence ---\n\n")

ptc_by_run <- pair_comparison %>%
  group_by(comparison, run) %>%
  summarise(
    n_pairs = n(),
    n_comp_ptc = sum(comp_has_ptc, na.rm = TRUE),
    comp_ptc_rate = mean(comp_has_ptc, na.rm = TRUE),
    .groups = "drop"
  )
print(ptc_by_run, n = 30)

prevalence_results <- list()
for (nmd_comp in c("C1", "C2")) {
  for (r in runs) {
    nmd_data <- pair_comparison %>%
      filter(comparison == nmd_comp, run == r, !is.na(comp_has_ptc))
    base_data <- pair_comparison %>%
      filter(comparison == "C4", run == r, !is.na(comp_has_ptc))
    if (nrow(nmd_data) < 5 || nrow(base_data) < 5) next

    nmd_ptc <- sum(nmd_data$comp_has_ptc)
    nmd_n <- nrow(nmd_data)
    base_ptc <- sum(base_data$comp_has_ptc)
    base_n <- nrow(base_data)

    ft <- fisher.test(matrix(c(nmd_ptc, nmd_n - nmd_ptc,
                                base_ptc, base_n - base_ptc), nrow = 2))
    prevalence_results[[length(prevalence_results) + 1]] <- tibble(
      nmd_comparison = nmd_comp, run = r,
      nmd_n = nmd_n, nmd_ptc_n = nmd_ptc,
      nmd_ptc_rate = nmd_ptc / nmd_n,
      baseline_n = base_n, baseline_ptc_n = base_ptc,
      baseline_ptc_rate = base_ptc / base_n,
      diff = nmd_ptc / nmd_n - base_ptc / base_n,
      odds_ratio = ft$estimate,
      or_ci_lower = ft$conf.int[1], or_ci_upper = ft$conf.int[2],
      p_value = ft$p.value
    )
  }
}

prevalence_df <- bind_rows(prevalence_results) %>%
  mutate(fdr_p = p.adjust(p_value, method = "BH"))

cat("\nSequence-based PTC prevalence: NMD vs C4 baseline\n")
prevalence_df %>%
  mutate(across(c(nmd_ptc_rate, baseline_ptc_rate, diff), ~round(., 3)),
         across(c(odds_ratio, or_ci_lower, or_ci_upper), ~round(., 2)),
         across(c(p_value, fdr_p), ~signif(., 3))) %>%
  print(n = 20, width = 140)

# ─── 7b. Event-type attribution enrichment ────────────────────────────────────

cat("\n--- 7b: Event-type enrichment among causal events ---\n\n")

# Map causal events to comparisons
causal_with_comp <- causal_results %>%
  filter(is_causal) %>%
  inner_join(
    pair_map %>% distinct(comparison, run, gene_id, dominant_isoform_id,
                           non_dominant_isoform_id),
    by = c("gene_id", "dominant_isoform_id", "non_dominant_isoform_id")
  )

# Aggregate: for deduplicated pairs, pick one row per pair
causal_dedup <- causal_with_comp %>%
  distinct(comparison, gene_id, dominant_isoform_id, non_dominant_isoform_id,
           .keep_all = TRUE)

cat("Causal event type distribution by comparison:\n")
event_enrichment <- causal_dedup %>%
  count(comparison, event_type, cause_type) %>%
  group_by(comparison) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  arrange(comparison, desc(n))

print(event_enrichment, n = 30)

# Fisher's test per event type: NMD (C1+C2 pooled) vs C4
cat("\nEvent type enrichment: NMD (C1+C2) vs C4\n")
event_type_tests <- list()
all_event_types <- unique(causal_dedup$event_type)

nmd_total <- sum(causal_dedup$comparison %in% c("C1", "C2"))
c4_total  <- sum(causal_dedup$comparison == "C4")

for (et in all_event_types) {
  nmd_et <- sum(causal_dedup$comparison %in% c("C1", "C2") &
                  causal_dedup$event_type == et)
  c4_et  <- sum(causal_dedup$comparison == "C4" &
                  causal_dedup$event_type == et)

  if (nmd_total > 0 && c4_total > 0) {
    ft <- fisher.test(matrix(c(nmd_et, nmd_total - nmd_et,
                                c4_et, c4_total - c4_et), nrow = 2))
    event_type_tests[[length(event_type_tests) + 1]] <- tibble(
      event_type = et,
      nmd_n = nmd_et, nmd_pct = round(100 * nmd_et / nmd_total, 1),
      c4_n = c4_et, c4_pct = round(100 * c4_et / c4_total, 1),
      odds_ratio = ft$estimate, p_value = ft$p.value
    )
  }
}

event_type_df <- bind_rows(event_type_tests) %>%
  mutate(fdr_p = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

print(event_type_df, n = 15, width = 120)

# ─── 7c. Within-pair PTC asymmetry ───────────────────────────────────────────

cat("\n--- 7c: Within-pair PTC asymmetry ---\n\n")

dedup_pairs <- pair_comparison %>%
  filter(!is.na(comp_has_ptc), !is.na(dom_has_ptc)) %>%
  distinct(comparison, dominant_isoform_id, non_dominant_isoform_id,
           .keep_all = TRUE)

pair_asym <- dedup_pairs %>%
  mutate(
    comp_ptc_not_dom = comp_has_ptc & !dom_has_ptc,
    dom_ptc_not_comp = dom_has_ptc & !comp_has_ptc,
    both_ptc = comp_has_ptc & dom_has_ptc,
    neither_ptc = !comp_has_ptc & !dom_has_ptc
  ) %>%
  group_by(comparison) %>%
  summarise(
    n = n(),
    comp_only_ptc = sum(comp_ptc_not_dom),
    dom_only_ptc = sum(dom_ptc_not_comp),
    both_ptc = sum(both_ptc),
    neither_ptc = sum(neither_ptc),
    pct_comp_only = round(100 * sum(comp_ptc_not_dom) / n(), 1),
    pct_dom_only = round(100 * sum(dom_ptc_not_comp) / n(), 1),
    .groups = "drop"
  )

print(pair_asym)

cat("\nMcNemar's test for PTC asymmetry within pairs:\n")
mcnemar_results <- list()
for (comp in c("C1", "C2", "C4")) {
  d <- dedup_pairs %>% filter(comparison == comp)
  b <- sum(d$comp_has_ptc & !d$dom_has_ptc)
  cc <- sum(!d$comp_has_ptc & d$dom_has_ptc)
  if (b + cc > 0) {
    mt <- mcnemar.test(matrix(c(
      sum(d$comp_has_ptc & d$dom_has_ptc), b, cc,
      sum(!d$comp_has_ptc & !d$dom_has_ptc)
    ), nrow = 2))
    cat(sprintf("  %s: comp_only=%d, dom_only=%d, chi2=%.2f, p=%g\n",
                comp, b, cc, mt$statistic, mt$p.value))
    mcnemar_results[[comp]] <- tibble(
      comparison = comp, comp_only_ptc = b, dom_only_ptc = cc,
      chi2 = mt$statistic, p_value = mt$p.value
    )
  }
}

# ─── 7d. Net frame shift analysis ────────────────────────────────────────────

cat("\n--- 7d: Net frame shift analysis ---\n\n")

frame_shift_summary <- pair_comparison %>%
  filter(!is.na(net_frame_status)) %>%
  distinct(comparison, gene_id, dominant_isoform_id, non_dominant_isoform_id,
           .keep_all = TRUE) %>%
  group_by(comparison, net_frame_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(comparison) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

cat("Net frame shift status by comparison:\n")
print(frame_shift_summary)

# Cross-tab: net frame shift × sequence PTC
cat("\nCross-tabulation: net frame shift × sequence-based PTC:\n")
frame_ptc_crosstab <- pair_comparison %>%
  filter(!is.na(net_frame_status), !is.na(comp_has_ptc)) %>%
  distinct(comparison, gene_id, dominant_isoform_id, non_dominant_isoform_id,
           .keep_all = TRUE) %>%
  count(comparison, net_frame_status, comp_has_ptc) %>%
  pivot_wider(names_from = comp_has_ptc, values_from = n,
              names_prefix = "ptc_", values_fill = 0)
print(frame_ptc_crosstab, n = 10)

# Test: NMD vs C4 for frame-disrupted rate
cat("\nFrame disruption rate: NMD vs C4\n")
frame_shift_tests <- list()
for (nmd_comp in c("C1", "C2")) {
  for (r in runs) {
    nmd_data <- pair_comparison %>%
      filter(comparison == nmd_comp, run == r, !is.na(net_frame_status))
    base_data <- pair_comparison %>%
      filter(comparison == "C4", run == r, !is.na(net_frame_status))
    if (nrow(nmd_data) < 5 || nrow(base_data) < 5) next

    nmd_fd <- sum(nmd_data$net_frame_status == "frame_disrupted")
    base_fd <- sum(base_data$net_frame_status == "frame_disrupted")

    ft <- fisher.test(matrix(c(
      nmd_fd, nrow(nmd_data) - nmd_fd,
      base_fd, nrow(base_data) - base_fd
    ), nrow = 2))

    frame_shift_tests[[length(frame_shift_tests) + 1]] <- tibble(
      nmd_comparison = nmd_comp, run = r,
      nmd_n = nrow(nmd_data),
      nmd_fd_n = nmd_fd,
      nmd_fd_rate = nmd_fd / nrow(nmd_data),
      baseline_n = nrow(base_data),
      baseline_fd_n = base_fd,
      baseline_fd_rate = base_fd / nrow(base_data),
      diff = nmd_fd / nrow(nmd_data) - base_fd / nrow(base_data),
      odds_ratio = ft$estimate,
      p_value = ft$p.value
    )
  }
}

frame_shift_test_df <- bind_rows(frame_shift_tests) %>%
  mutate(fdr_p = p.adjust(p_value, method = "BH"))

frame_shift_test_df %>%
  mutate(across(c(nmd_fd_rate, baseline_fd_rate, diff), ~round(., 3)),
         odds_ratio = round(odds_ratio, 2),
         across(c(p_value, fdr_p), ~signif(., 3))) %>%
  print(n = 20, width = 140)

# ─── 7e. Alt_TSS/Alt_TES enrichment (uORF-mediated NMD mechanism) ───────────

cat("\n--- 7e: Alt_TSS/Alt_TES enrichment in NMD pairs ---\n\n")

# Among same-start-codon pairs, test whether Alt_TSS/Alt_TES events are
# enriched in NMD comparisons (C1/C2) vs baseline (C4).
# Rationale: Alt_TSS can introduce upstream ORFs (uORFs) in the 5'UTR,
# which are a known NMD trigger independent of the main ORF.

# Get ALL events for same-start-codon pairs (not just CDS-overlapping)
tss_tes_events <- profiles %>%
  semi_join(pair_results %>% filter(both_valid_atg),
            by = c("gene_id", "dominant_isoform_id", "non_dominant_isoform_id")) %>%
  select(gene_id, dominant_isoform_id, non_dominant_isoform_id,
         detailed_events) %>%
  unnest(detailed_events, names_sep = "_") %>%
  select(
    gene_id, dominant_isoform_id, non_dominant_isoform_id,
    event_type = detailed_events_event_type
  ) %>%
  filter(event_type %in% c("Alt_TSS", "Alt_TES"))

# Flag pairs that have Alt_TSS and/or Alt_TES
tss_tes_per_pair <- tss_tes_events %>%
  group_by(gene_id, dominant_isoform_id, non_dominant_isoform_id) %>%
  summarise(
    has_alt_tss = any(event_type == "Alt_TSS"),
    has_alt_tes = any(event_type == "Alt_TES"),
    .groups = "drop"
  )

# Join with comparison mapping
tss_tes_comparison <- pair_comparison %>%
  left_join(tss_tes_per_pair,
            by = c("gene_id", "dominant_isoform_id", "non_dominant_isoform_id")) %>%
  mutate(
    has_alt_tss = replace_na(has_alt_tss, FALSE),
    has_alt_tes = replace_na(has_alt_tes, FALSE)
  )

# Rates by comparison
cat("Alt_TSS/Alt_TES rates by comparison:\n")
tss_tes_comparison %>%
  distinct(comparison, gene_id, dominant_isoform_id, non_dominant_isoform_id,
           .keep_all = TRUE) %>%
  group_by(comparison) %>%
  summarise(
    n = n(),
    n_alt_tss = sum(has_alt_tss),
    pct_alt_tss = round(100 * mean(has_alt_tss), 1),
    n_alt_tes = sum(has_alt_tes),
    pct_alt_tes = round(100 * mean(has_alt_tes), 1),
    .groups = "drop"
  ) %>%
  print()

# Fisher's exact test: NMD vs C4 for Alt_TSS and Alt_TES presence
cat("\nAlt_TSS/Alt_TES enrichment: NMD vs C4 baseline\n")
tss_tes_tests <- list()
for (event_name in c("Alt_TSS", "Alt_TES")) {
  has_col <- if (event_name == "Alt_TSS") "has_alt_tss" else "has_alt_tes"

  for (nmd_comp in c("C1", "C2")) {
    for (r in runs) {
      nmd_data <- tss_tes_comparison %>%
        filter(comparison == nmd_comp, run == r)
      base_data <- tss_tes_comparison %>%
        filter(comparison == "C4", run == r)
      if (nrow(nmd_data) < 5 || nrow(base_data) < 5) next

      nmd_has <- sum(nmd_data[[has_col]])
      nmd_n <- nrow(nmd_data)
      base_has <- sum(base_data[[has_col]])
      base_n <- nrow(base_data)

      ft <- fisher.test(matrix(c(nmd_has, nmd_n - nmd_has,
                                  base_has, base_n - base_has), nrow = 2))
      tss_tes_tests[[length(tss_tes_tests) + 1]] <- tibble(
        event_type = event_name, nmd_comparison = nmd_comp, run = r,
        nmd_n = nmd_n, nmd_has = nmd_has,
        nmd_rate = nmd_has / nmd_n,
        baseline_n = base_n, baseline_has = base_has,
        baseline_rate = base_has / base_n,
        diff = nmd_has / nmd_n - base_has / base_n,
        odds_ratio = ft$estimate,
        or_ci_lower = ft$conf.int[1], or_ci_upper = ft$conf.int[2],
        p_value = ft$p.value
      )
    }
  }
}

tss_tes_df <- bind_rows(tss_tes_tests) %>%
  mutate(fdr_p = p.adjust(p_value, method = "BH"))

tss_tes_df %>%
  mutate(across(c(nmd_rate, baseline_rate, diff), ~round(., 3)),
         across(c(odds_ratio, or_ci_lower, or_ci_upper), ~round(., 2)),
         across(c(p_value, fdr_p), ~signif(., 3))) %>%
  print(n = 30, width = 140)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 8: Save outputs
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 8: Saving outputs ---\n\n")

output_dir <- sprintf("results/ptc/results/nonNMD_%s", threshold)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Per-pair results (RDS + TSV)
pair_out <- pair_results %>%
  select(gene_id, dominant_isoform_id, non_dominant_isoform_id,
         dom_start_codon, both_valid_atg,
         dom_scanned_orf_length_nt, comp_scanned_orf_length_nt,
         dom_stop_codon_type, comp_stop_codon_type,
         dom_reached_end, comp_reached_end,
         dom_ptc_distance, comp_ptc_distance,
         dom_has_ptc, comp_has_ptc, orf_length_diff,
         n_cds_events, net_shift, net_frame_status)

saveRDS(pair_out, file.path(output_dir, "sequence_ptc_pair_results.rds"))
write_tsv(pair_out, file.path(output_dir, "sequence_ptc_pair_results.tsv"))
cat(sprintf("  Saved pair results: %d pairs\n", nrow(pair_out)))

# Prevalence test results
write_tsv(prevalence_df, file.path(output_dir, "sequence_ptc_prevalence.tsv"))
cat(sprintf("  Saved prevalence tests: %d tests\n", nrow(prevalence_df)))

# Causal event attribution
causal_out <- causal_results %>%
  filter(is_causal | is_contributing) %>%
  select(gene_id, dominant_isoform_id, non_dominant_isoform_id,
         event_type, direction, five_prime, three_prime, bp_diff,
         cumulative_shift, is_causal, is_contributing, cause_type,
         net_shift, net_frame_status)

write_tsv(causal_out, file.path(output_dir, "sequence_ptc_causal_events.tsv"))
cat(sprintf("  Saved causal events: %d events\n", nrow(causal_out)))

# Event type enrichment
write_tsv(event_type_df, file.path(output_dir, "sequence_ptc_event_type_enrichment.tsv"))
cat(sprintf("  Saved event type enrichment: %d event types\n", nrow(event_type_df)))

# Validation
write_tsv(validation_out, file.path(output_dir, "sequence_ptc_validation.tsv"))
cat("  Saved validation results\n")

# Frame shift summary
write_tsv(frame_shift_summary, file.path(output_dir, "sequence_ptc_frame_shift_summary.tsv"))
cat("  Saved frame shift summary\n")

# Within-pair asymmetry
write_tsv(pair_asym, file.path(output_dir, "sequence_ptc_within_pair_asymmetry.tsv"))
cat("  Saved within-pair asymmetry\n")

# Alt_TSS/Alt_TES enrichment
write_tsv(tss_tes_df, file.path(output_dir, "sequence_ptc_alt_tss_tes_enrichment.tsv"))
cat(sprintf("  Saved Alt_TSS/Alt_TES enrichment: %d tests\n", nrow(tss_tes_df)))

cat(sprintf("\nAll outputs saved to: %s/\n", output_dir))
cat("\nDone.\n")
