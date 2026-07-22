#!/usr/bin/env Rscript
#
# Extract Failure Cases from Production Data
#
# Reads reconstruction_verification.rds to identify FAIL pairs, then extracts
# their exon structures from isoform_structures_filtered.rds and appends them
# to the test GTF and pairs files.
#
# This is a one-time helper script. Re-run whenever you want to refresh the
# real-data test cases (e.g., after a new full pipeline run).
#
# Usage:
#   Rscript scripts/tests/extract_failure_cases.R
#
# Input:
#   - data/reconstruction_verification.rds
#   - data/isoform_structures_filtered.rds
#   - scripts/tests/test_data/test_cases.gtf   (synthetic base — will be appended to)
#   - scripts/tests/test_data/pairs.tsv         (synthetic base — will be appended to)
#
# Output:
#   - scripts/tests/test_data/test_cases.gtf   (synthetic + real failures)
#   - scripts/tests/test_data/pairs.tsv        (synthetic + real failures)
#

library(tidyverse)

# ==============================================================================
# Paths
# ==============================================================================

base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"

verification_file <- file.path(base_dir, "data/reconstruction_verification.rds")
structures_file <- file.path(base_dir, "data/isoform_structures_filtered.rds")
synthetic_gtf <- file.path(base_dir, "development/reconstruction/synthetic_data/base_events.gtf")
synthetic_pairs <- file.path(base_dir, "development/reconstruction/synthetic_data/pairs.tsv")
output_gtf <- file.path(base_dir, "scripts/tests/test_data/test_cases.gtf")
output_pairs <- file.path(base_dir, "scripts/tests/test_data/pairs.tsv")

for (f in c(verification_file, structures_file, synthetic_gtf, synthetic_pairs)) {
  if (!file.exists(f)) stop(sprintf("Required file not found: %s", f))
}

cat("\n")
cat("======================================================================\n")
cat("  EXTRACT FAILURE CASES FOR TEST DATASET\n")
cat("======================================================================\n\n")

# ==============================================================================
# 1. Load verification results and identify failures
# ==============================================================================

cat("Loading verification results...\n")
verification <- readRDS(verification_file)
cat(sprintf("  Total pairs: %d\n", nrow(verification)))
cat(sprintf("  PASS: %d, FAIL: %d, ERROR: %d\n",
            sum(verification$status == "PASS"),
            sum(verification$status == "FAIL"),
            sum(verification$status == "ERROR")))

failures <- verification %>%
  filter(status %in% c("FAIL", "ERROR"))

cat(sprintf("\n  Extracting %d failure/error cases\n", nrow(failures)))

if (nrow(failures) == 0) {
  cat("\nNo failures to extract. Test data will contain only synthetic cases.\n")
  # Still rebuild from synthetic base
  file.copy(synthetic_gtf, output_gtf, overwrite = TRUE)
  file.copy(synthetic_pairs, output_pairs, overwrite = TRUE)
  cat("Done.\n")
  quit(status = 0)
}

# ==============================================================================
# 2. Load isoform structures (compact format)
# ==============================================================================

cat("\nLoading isoform structures...\n")
structures <- readRDS(structures_file)
cat(sprintf("  %d isoforms loaded\n", nrow(structures)))

# ==============================================================================
# 3. Extract exon coordinates and build GTF + pairs lines
# ==============================================================================

cat("\nExtracting exon coordinates for failure pairs...\n")

# Collect all isoform IDs we need (both dominant and comparator)
needed_isoforms <- unique(c(failures$dominant_isoform_id, failures$comparator_isoform_id))
cat(sprintf("  Need structures for %d unique isoforms\n", length(needed_isoforms)))

# Filter to needed isoforms
needed_structures <- structures %>%
  filter(isoform_id %in% needed_isoforms)

cat(sprintf("  Found %d/%d in isoform_structures\n",
            nrow(needed_structures), length(needed_isoforms)))

missing <- setdiff(needed_isoforms, needed_structures$isoform_id)
if (length(missing) > 0) {
  cat(sprintf("  WARNING: %d isoforms not found: %s\n",
              length(missing), paste(head(missing, 5), collapse = ", ")))
}

# Build GTF lines and pairs entries per failure pair.
# IMPORTANT: Include ALL isoforms for each failing gene so that union exons
# match production (Script 03 builds UEs from all isoforms). Use the pipeline
# gene_id (from verification data), not the native gene_id from isoform_structures,
# since PacBio isoforms have PB-based gene_ids (e.g., "PB.5114") while the
# pipeline uses ENSG gene_ids.

# First, build a mapping from pipeline gene_id -> all isoform_ids for that gene
# We need the gene_id -> isoform_id mapping from the filtered structures
cat("  Building gene -> isoform mapping for all genes with failures...\n")

# Get all pipeline gene_ids that have failures
failure_gene_ids <- unique(failures$gene_id)
cat(sprintf("  %d unique genes with failures\n", length(failure_gene_ids)))

# For each failure gene, find ALL isoforms via the isoform_structures data
# The tricky part: pipeline gene_ids are ENSG-based, but structures may use PB gene_ids
# So we match via isoform_id — find all isoforms that share a gene with any failure isoform
failure_isoform_ids <- unique(c(failures$dominant_isoform_id, failures$comparator_isoform_id))
failure_native_gene_ids <- structures %>%
  filter(isoform_id %in% failure_isoform_ids) %>%
  pull(gene_id) %>%
  unique()

# Get ALL isoforms for these genes (using native gene_ids)
all_gene_isoforms <- structures %>%
  filter(gene_id %in% failure_native_gene_ids)

cat(sprintf("  Total isoforms across %d failure genes: %d\n",
            length(failure_native_gene_ids), nrow(all_gene_isoforms)))

# Build the native-to-pipeline gene_id mapping (via BOTH dom and comp isoforms)
gene_id_map_dom <- failures %>%
  select(gene_id, isoform_id = dominant_isoform_id)
gene_id_map_comp <- failures %>%
  select(gene_id, isoform_id = comparator_isoform_id)
gene_id_map <- bind_rows(gene_id_map_dom, gene_id_map_comp) %>%
  left_join(structures %>% select(isoform_id, native_gene_id = gene_id),
            by = "isoform_id") %>%
  select(pipeline_gene_id = gene_id, native_gene_id) %>%
  distinct()

cat(sprintf("  Gene ID mappings: %d (pipeline -> native)\n", nrow(gene_id_map)))

gtf_lines <- list()
pair_lines <- list()
skipped <- 0
written_isoforms <- character(0)

# Write ALL isoforms for each failure gene
for (j in seq_len(nrow(gene_id_map))) {
  pipeline_gid <- gene_id_map$pipeline_gene_id[j]
  native_gid <- gene_id_map$native_gene_id[j]

  gene_isoforms <- all_gene_isoforms %>% filter(gene_id == native_gid)

  for (k in seq_len(nrow(gene_isoforms))) {
    iso <- gene_isoforms[k, ]
    iso_key <- paste0(pipeline_gid, "::", iso$isoform_id)
    if (iso_key %in% written_isoforms) next
    written_isoforms <- c(written_isoforms, iso_key)

    starts <- iso$exon_starts[[1]]
    ends <- iso$exon_ends[[1]]
    chr <- iso$seqnames
    strand <- iso$strand

    for (m in seq_along(starts)) {
      gtf_lines[[length(gtf_lines) + 1]] <- sprintf(
        '%s\tproduction\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; exon_number "%d";',
        chr, starts[m], ends[m], strand, pipeline_gid, iso$isoform_id, m
      )
    }
  }
}

# Build pairs entries
for (j in seq_len(nrow(failures))) {
  f <- failures[j, ]

  dom_struct <- needed_structures %>% filter(isoform_id == f$dominant_isoform_id)
  comp_struct <- needed_structures %>% filter(isoform_id == f$comparator_isoform_id)

  if (nrow(dom_struct) == 0 || nrow(comp_struct) == 0) {
    skipped <- skipped + 1
    next
  }

  chr <- dom_struct$seqnames[1]
  strand <- dom_struct$strand[1]

  pair_lines[[length(pair_lines) + 1]] <- sprintf(
    "%s\t%s\t%s\t%s\t%s",
    f$gene_id, f$dominant_isoform_id, f$comparator_isoform_id, chr, strand
  )
}

if (skipped > 0) {
  cat(sprintf("  Skipped %d pairs with missing structures\n", skipped))
}
cat(sprintf("  Generated %d GTF exon lines (all isoforms for failure genes)\n", length(gtf_lines)))
cat(sprintf("  Generated %d pairs entries\n", length(pair_lines)))

# ==============================================================================
# 5. Combine with synthetic data and write output
# ==============================================================================

cat("\nCombining with synthetic test data...\n")

# Read synthetic base files
synthetic_gtf_lines <- readLines(synthetic_gtf)
synthetic_pairs_lines <- readLines(synthetic_pairs)

# Add separator comment and real-data lines to GTF
combined_gtf <- c(
  synthetic_gtf_lines,
  "",
  "# =====================================================",
  sprintf("# REAL-DATA FAILURE CASES (extracted %s)", Sys.Date()),
  sprintf("# %d pairs from reconstruction_verification.rds", length(pair_lines)),
  "# =====================================================",
  unlist(gtf_lines)
)

# Add real-data pairs (skip header since synthetic already has it)
combined_pairs <- c(
  synthetic_pairs_lines,
  unlist(pair_lines)
)

# Write output
writeLines(combined_gtf, output_gtf)
writeLines(combined_pairs, output_pairs)

# Count final totals
n_synthetic_pairs <- length(synthetic_pairs_lines) - 1  # minus header
n_real_pairs <- length(pair_lines)

cat(sprintf("\nOutput written:\n"))
cat(sprintf("  GTF:   %s (%d lines)\n", output_gtf, length(combined_gtf)))
cat(sprintf("  Pairs: %s (%d synthetic + %d real = %d total)\n",
            output_pairs, n_synthetic_pairs, n_real_pairs,
            n_synthetic_pairs + n_real_pairs))

cat("\n======================================================================\n")
cat("  EXTRACTION COMPLETE\n")
cat("======================================================================\n\n")
