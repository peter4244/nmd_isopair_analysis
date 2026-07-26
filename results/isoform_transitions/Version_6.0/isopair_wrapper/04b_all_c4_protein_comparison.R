#!/usr/bin/env Rscript
# 04b_all_c4_protein_comparison.R
#
# Extends protein sequence comparison from the 107 productive frameshifts
# (04_productive_frameshift_precompute.R) to ALL C4 coding comparator isoforms.
# For each C4 ref+comp pair where both are coding, extracts protein sequences
# from the SQANTI FASTA and computes actual_shared_aa via prefix comparison.
#
# Output: data_mashr/analysis_cache/all_c4_protein_comparison_allsamples.rds
#
# Run from: results/isoform_transitions/Version_6.0/isopair_wrapper/
# Usage: Rscript 04b_all_c4_protein_comparison.R [--test]

args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args

cat("\n")
cat("==================================================================\n")
cat("   04b: All C4 Protein Sequence Comparison\n")
if (test_mode) cat("   *** TEST MODE: first 50 pairs only ***\n")
cat("==================================================================\n\n")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 0: Setup & data loading
# ═══════════════════════════════════════════════════════════════════════════════

cat("--- Section 0: Loading data ---\n\n")

data_dir <- "data_mashr"
cache_dir <- file.path(data_dir, "analysis_cache")

fc <- readRDS(file.path(cache_dir, "fc_c4_allsamples.rds"))
cds <- readRDS(file.path(data_dir, "cds.rds"))
ptc <- readRDS(file.path(data_dir, "ptc.rds"))

# All C4 pairs from frame comparison
ps <- fc$pair_summary
cat(sprintf("  Total C4 pairs (frame comparison): %d\n", nrow(ps)))

# Restrict to coding pairs (both ref and comp are coding)
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
coding_pairs <- ps[ps$reference_isoform_id %in% coding_ids &
                     ps$comparator_isoform_id %in% coding_ids, ]
cat(sprintf("  Coding pairs (both ref+comp coding): %d\n", nrow(coding_pairs)))

if (test_mode) {
  coding_pairs <- coding_pairs[1:min(50, nrow(coding_pairs)), ]
  cat(sprintf("  TEST MODE: using %d pairs\n", nrow(coding_pairs)))
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Protein sequence extraction
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 1: Extracting protein sequences ---\n\n")

protein_fasta <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.faa"
stopifnot(file.exists(protein_fasta))

# Collect all isoform IDs to extract
all_ids <- unique(c(coding_pairs$reference_isoform_id,
                    coding_pairs$comparator_isoform_id))
cat(sprintf("  Unique isoforms to extract: %d\n", length(all_ids)))

# Write IDs to temp file, then use awk for single-pass extraction
id_tmpfile <- tempfile(fileext = ".txt")
writeLines(all_ids, id_tmpfile)

awk_cmd <- sprintf(
  "awk 'NR==FNR{ids[$1]; next} /^>/{id=$1; sub(/^>/,\"\",id); sub(/\\t.*/,\"\",id); p=(id in ids)} p' '%s' '%s'",
  id_tmpfile, protein_fasta
)
fasta_lines <- system(awk_cmd, intern = TRUE)
unlink(id_tmpfile)

# Parse FASTA lines into named sequences
seq_list <- list()
current_id <- NULL
current_seq <- character(0)
for (line in fasta_lines) {
  if (startsWith(line, ">")) {
    if (!is.null(current_id)) {
      seq_list[[current_id]] <- paste0(current_seq, collapse = "")
    }
    header <- sub("^>", "", line)
    current_id <- sub("\t.*", "", header)
    current_seq <- character(0)
  } else {
    current_seq <- c(current_seq, line)
  }
}
if (!is.null(current_id)) {
  seq_list[[current_id]] <- paste0(current_seq, collapse = "")
}

cat(sprintf("  Extracted sequences: %d / %d\n", length(seq_list), length(all_ids)))
missing_ids <- setdiff(all_ids, names(seq_list))
if (length(missing_ids) > 0) {
  cat(sprintf("  WARNING: %d IDs not found in protein FASTA\n", length(missing_ids)))
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Protein sequence comparison
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 2: Computing protein comparisons ---\n\n")

# Build orf_length lookup from ptc.rds
orf_lookup <- setNames(ptc$orf_length, ptc$isoform_id)

results <- vector("list", nrow(coding_pairs))

n_both_seq <- 0
n_missing_seq <- 0
n_same_start <- 0
n_diff_start <- 0

for (i in seq_len(nrow(coding_pairs))) {
  if (i %% 1000 == 0) cat(sprintf("  [%d/%d]\n", i, nrow(coding_pairs)))

  row <- coding_pairs[i, ]
  ref_id <- row$reference_isoform_id
  comp_id <- row$comparator_isoform_id

  ref_seq <- seq_list[[ref_id]]
  comp_seq <- seq_list[[comp_id]]

  if (is.null(ref_seq) || is.null(comp_seq)) {
    n_missing_seq <- n_missing_seq + 1
    next
  }
  n_both_seq <- n_both_seq + 1

  # Remove stop codon marker
  ref_clean <- sub("\\*$", "", ref_seq)
  comp_clean <- sub("\\*$", "", comp_seq)

  ref_len <- nchar(ref_clean)
  comp_len <- nchar(comp_clean)

  # Compute actual_shared_aa via prefix comparison
  ref_chars <- strsplit(ref_clean, "")[[1]]
  comp_chars <- strsplit(comp_clean, "")[[1]]
  min_len <- min(length(ref_chars), length(comp_chars))

  shared_aa <- min_len  # default: identical up to shorter
  for (j in seq_len(min_len)) {
    if (ref_chars[j] != comp_chars[j]) {
      shared_aa <- j - 1L
      break
    }
  }

  # Categorize by start codon
  same_start <- shared_aa > 0 || (ref_len == 0 && comp_len == 0)
  # More precise: check if frame_category indicates same start
  frame_cat <- row$frame_category
  is_same_start <- grepl("^same_start", frame_cat)
  if (is_same_start) n_same_start <- n_same_start + 1 else n_diff_start <- n_diff_start + 1

  # Novel protein content
  novel_aa <- comp_len - shared_aa
  pct_novel_protein <- if (comp_len > 0) 100 * novel_aa / comp_len else 0

  # Get orf_length and PTC status
  comp_orf <- orf_lookup[comp_id]
  comp_has_ptc <- comp_id %in% ptc$isoform_id[ptc$has_ptc == TRUE]

  results[[i]] <- data.frame(
    gene_id = row$gene_id,
    reference_isoform_id = ref_id,
    comparator_isoform_id = comp_id,
    frame_category = frame_cat,
    ref_protein_len = ref_len,
    comp_protein_len = comp_len,
    actual_shared_aa = shared_aa,
    novel_aa = novel_aa,
    pct_novel_protein = round(pct_novel_protein, 2),
    comp_has_ptc = comp_has_ptc,
    comp_orf_length = if (!is.na(comp_orf)) as.integer(comp_orf) else NA_integer_,
    stringsAsFactors = FALSE
  )
}

comparison <- do.call(rbind, results)
cat(sprintf("\n  Pairs with both sequences: %d\n", n_both_seq))
cat(sprintf("  Pairs missing sequence: %d\n", n_missing_seq))
cat(sprintf("  Same-start pairs: %d\n", n_same_start))
cat(sprintf("  Diff-start pairs: %d\n", n_diff_start))
cat(sprintf("  Median pct_novel_protein: %.1f%%\n",
            median(comparison$pct_novel_protein, na.rm = TRUE)))
cat(sprintf("  Mean pct_novel_protein: %.1f%%\n",
            mean(comparison$pct_novel_protein, na.rm = TRUE)))

# Summary by frame category
cat("\n  By frame category:\n")
for (fc_cat in unique(comparison$frame_category)) {
  sub <- comparison[comparison$frame_category == fc_cat, ]
  cat(sprintf("    %s: n=%d, median novel=%.1f%%\n",
              fc_cat, nrow(sub), median(sub$pct_novel_protein, na.rm = TRUE)))
}


# ═══════════════════════════════════════════════════════════════════════════════
# SAVE RESULTS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Saving results ---\n\n")

output <- list(
  comparison = comparison,
  metadata = list(
    run_date = Sys.time(),
    n_total_c4_pairs = nrow(ps),
    n_coding_pairs = nrow(coding_pairs),
    n_compared = nrow(comparison),
    n_missing_seq = n_missing_seq,
    fasta_file = protein_fasta,
    test_mode = test_mode
  )
)

out_file <- file.path(cache_dir, "all_c4_protein_comparison_allsamples.rds")
saveRDS(output, out_file)
cat(sprintf("  Saved to: %s\n", out_file))
cat(sprintf("  Comparison rows: %d\n", nrow(comparison)))
cat("\nDone.\n")
