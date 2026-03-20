#!/usr/bin/env Rscript
# ==============================================================================
# 05p_cnn_data_prep.R
#
# Prepare data for CNN NMD prediction model.
# Extracts full spliced transcript sequence (padded/truncated to 8192 nt)
# with exon-exon junction positions for each classified coding isoform.
#
# Output:
#   - data_mashr/analysis_cache/cnn_data.tsv
#     Columns: isoform_id, is_nmd, chr, tx_length, seq, junctions
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

FASTA_PATH <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta"
OUTPUT_PATH <- "data_mashr/analysis_cache/cnn_data.tsv"
MAX_LEN <- 4096L

cat("=== CNN Data Preparation ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# Load metadata
nmd_class <- readRDS("data_mashr/nmd_classification.rds")
cds <- readRDS("data_mashr/cds.rds")
structures <- readRDS("data_mashr/structures.rds")

nmd_ids <- nmd_class$all_samples$nmd
non_nmd_ids <- nmd_class$all_samples$non_nmd
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
target_ids <- intersect(c(nmd_ids, non_nmd_ids), coding_ids)
cat("Target isoforms:", length(target_ids), "\n")

# Extract sequences
cat("Extracting sequences from FASTA...\n")
id_file <- tempfile(fileext = ".txt")
writeLines(target_ids, id_file)
awk_cmd <- sprintf(
  "awk 'BEGIN{while((getline id < \"%s\") > 0) want[id]=1} /^>/{if(seq && found) print hdr \"\\t\" seq; hdr=substr($1,2); found=(hdr in want); seq=\"\"; next} found{seq=seq$0} END{if(seq && found) print hdr \"\\t\" seq}' '%s'",
  id_file, FASTA_PATH
)
t0 <- proc.time()
raw_output <- system(awk_cmd, intern = TRUE)
cat("  awk completed in", round((proc.time() - t0)[3], 1), "seconds\n")
unlink(id_file)

split_lines <- strsplit(raw_output, "\t", fixed = TRUE)
seq_names <- vapply(split_lines, `[`, character(1), 1)
seq_seqs <- toupper(vapply(split_lines, `[`, character(1), 2))
cat("  Sequences:", length(seq_names), "\n")
cat("  Missing:", sum(!target_ids %in% seq_names), "\n")

# Compute junction positions in transcript space
cat("\nComputing junction positions...\n")
target_structs <- structures[structures$isoform_id %in% seq_names, ]
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

# Build output
cat("\nPreparing output (max", MAX_LEN, "nt)...\n")
chr_lookup <- setNames(structures$chr, structures$isoform_id)

out_lines <- character(length(seq_names))
n_truncated <- 0L
n_padded <- 0L

for (i in seq_along(seq_names)) {
  iso <- seq_names[i]
  seq <- seq_seqs[i]
  tx_len <- nchar(seq)

  # Pad or truncate sequence
  if (tx_len > MAX_LEN) {
    seq_out <- substr(seq, 1, MAX_LEN)
    n_truncated <- n_truncated + 1L
  } else if (tx_len < MAX_LEN) {
    seq_out <- paste0(seq, paste(rep("N", MAX_LEN - tx_len), collapse = ""))
    n_padded <- n_padded + 1L
  } else {
    seq_out <- seq
  }

  is_nmd <- as.integer(iso %in% nmd_ids)
  chr <- chr_lookup[iso]

  # Junction positions within MAX_LEN
  junctions <- junction_map[[iso]]
  if (is.null(junctions)) junctions <- integer(0)
  junctions <- junctions[junctions <= MAX_LEN]
  junc_str <- if (length(junctions) > 0) paste(junctions, collapse = ",") else ""

  out_lines[i] <- paste(iso, is_nmd, chr, tx_len, seq_out, junc_str, sep = "\t")
}

cat("  Full length:", length(seq_names) - n_truncated - n_padded, "\n")
cat("  Truncated (>", MAX_LEN, "):", n_truncated, "\n")
cat("  Padded (<", MAX_LEN, "):", n_padded, "\n")

con <- file(OUTPUT_PATH, "w")
writeLines("isoform_id\tis_nmd\tchr\ttx_length\tseq\tjunctions", con)
writeLines(out_lines, con)
close(con)

cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Rows:", length(out_lines), "\n")
cat("Completed:", format(Sys.time()), "\n")
