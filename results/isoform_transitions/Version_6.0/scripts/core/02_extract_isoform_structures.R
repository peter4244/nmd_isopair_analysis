#!/usr/bin/env Rscript

################################################################################
# Script 02: Extract Isoform Structures from GENCODE and SQANTI GFF Files
################################################################################
#
# Purpose:
#   Extract exon coordinates for all major isoforms from genomic annotation files.
#   This creates the structural foundation for downstream union exon construction
#   and splicing profile analysis.
#
# Input:
#   - data/expression_data.rds (major isoforms list)
#   - reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz
#   - isoseq/collapse/merge-collapsed.gff (PacBio isoforms)
#
# Output:
#   - data/isoform_structures.rds
#
# Structure:
#   One row per isoform with nested exon coordinates:
#   - isoform_id, gene_id, seqnames, strand
#   - n_exons, exon_starts (list), exon_ends (list)
#   - tx_start, tx_end, n_junctions
#
################################################################################

library(tidyverse)
library(rtracklayer)

# ==============================================================================
# Input Validation Helpers
# ==============================================================================

validate_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required input file not found: %s\nHave you run the prerequisite scripts?", path))
  }
}

validate_columns <- function(df, required_cols, name) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing required columns: %s", name, paste(missing, collapse = ", ")))
  }
}

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 2: Extract Isoform Structures from GFF Files          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ==============================================================================
# 0. Validate Inputs
# ==============================================================================

cat("Validating inputs...\n")
validate_file("data/expression_data.rds")
validate_file("/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz")
validate_file("/Users/petecastaldi/claude_projects/nmd/isoseq/collapse/merge-collapsed.gff")
cat("  All required input files found.\n\n")

# ==============================================================================
# 1. Load Major Isoforms List
# ==============================================================================

cat("Loading major isoforms list...\n")
expression_data <- readRDS("data/expression_data.rds")
validate_columns(expression_data, c("isoform_id", "gene_id"), "expression_data")
cat("  Required columns verified.\n")

major_isoforms <- expression_data %>%
  distinct(isoform_id, gene_id)

cat(sprintf("  Loaded %d major isoforms from %d genes\n",
            nrow(major_isoforms), n_distinct(major_isoforms$gene_id)))

# Separate GENCODE (ENST) and PacBio (PB) isoforms
# Strip version numbers for matching with GFF
gencode_isoforms <- major_isoforms %>%
  filter(str_detect(isoform_id, "^ENST")) %>%
  mutate(isoform_id_no_version = str_replace(isoform_id, "\\.\\d+$", ""))

pacbio_isoforms <- major_isoforms %>%
  filter(str_detect(isoform_id, "^PB\\."))

cat(sprintf("  GENCODE isoforms: %d\n", nrow(gencode_isoforms)))
cat(sprintf("  PacBio isoforms: %d\n", nrow(pacbio_isoforms)))

# ==============================================================================
# 2. Extract GENCODE Structures
# ==============================================================================

cat("\nExtracting GENCODE isoform structures...\n")
cat("  Loading GENCODE GFF3 (this may take several minutes)...\n")

gencode_gff <- import("/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz")

cat("  Filtering to exon features...\n")
gencode_exons <- gencode_gff %>%
  as_tibble() %>%
  filter(type == "exon") %>%
  select(transcript_id, gene_id, seqnames, start, end, strand) %>%
  # Strip version numbers from IDs for matching
  mutate(
    transcript_id_no_version = str_replace(transcript_id, "\\.\\d+$", ""),
    gene_id = str_replace(gene_id, "\\.\\d+$", "")
  ) %>%
  # Filter to major isoforms only (using version-stripped IDs)
  filter(transcript_id_no_version %in% gencode_isoforms$isoform_id_no_version) %>%
  # Keep original transcript_id with version for final output
  select(-transcript_id_no_version) %>%
  dplyr::rename(isoform_id = transcript_id)

cat(sprintf("  Extracted %d exons for %d GENCODE isoforms\n",
            nrow(gencode_exons), n_distinct(gencode_exons$isoform_id)))

# Build nested structure
cat("  Building nested exon structures...\n")
gencode_structures <- gencode_exons %>%
  group_by(isoform_id) %>%
  arrange(start) %>%
  summarise(
    gene_id = dplyr::first(gene_id),
    seqnames = dplyr::first(seqnames),
    strand = dplyr::first(strand),
    n_exons = n(),
    exon_starts = list(start),
    exon_ends = list(end),
    tx_start = min(start),
    tx_end = max(end),
    .groups = "drop"
  )

# Calculate number of junctions (n_exons - 1)
gencode_structures <- gencode_structures %>%
  mutate(n_junctions = n_exons - 1)

cat(sprintf("  Created structures for %d GENCODE isoforms\n",
            nrow(gencode_structures)))

# ==============================================================================
# 3. Extract PacBio Structures (if any)
# ==============================================================================

if (nrow(pacbio_isoforms) > 0) {
  cat("\nExtracting PacBio isoform structures...\n")
  cat("  Loading PacBio GFF...\n")

  pacbio_gff <- import("/Users/petecastaldi/claude_projects/nmd/isoseq/collapse/merge-collapsed.gff")

  cat("  Filtering to exon features...\n")
  pacbio_exons <- pacbio_gff %>%
    as_tibble() %>%
    filter(type == "exon") %>%
    select(transcript_id, gene_id, seqnames, start, end, strand) %>%
    # Filter to major isoforms only
    filter(transcript_id %in% pacbio_isoforms$isoform_id)

  cat(sprintf("  Extracted %d exons for %d PacBio isoforms\n",
              nrow(pacbio_exons), n_distinct(pacbio_exons$transcript_id)))

  # Build nested structure
  cat("  Building nested exon structures...\n")
  pacbio_structures <- pacbio_exons %>%
    group_by(transcript_id) %>%
    arrange(start) %>%
    summarise(
      gene_id = dplyr::first(gene_id),
      seqnames = dplyr::first(seqnames),
      strand = dplyr::first(strand),
      n_exons = n(),
      exon_starts = list(start),
      exon_ends = list(end),
      tx_start = min(start),
      tx_end = max(end),
      .groups = "drop"
    ) %>%
    dplyr::rename(isoform_id = transcript_id)

  pacbio_structures <- pacbio_structures %>%
    mutate(n_junctions = n_exons - 1)

  cat(sprintf("  Created structures for %d PacBio isoforms\n",
              nrow(pacbio_structures)))
} else {
  cat("\nNo PacBio isoforms to extract (skipping)\n")
  pacbio_structures <- tibble()
}

# ==============================================================================
# 4. Combine and Validate
# ==============================================================================

cat("\nCombining GENCODE and PacBio structures...\n")

if (nrow(pacbio_structures) > 0) {
  all_structures <- bind_rows(gencode_structures, pacbio_structures)
} else {
  all_structures <- gencode_structures
}

cat(sprintf("  Total isoforms: %d\n", nrow(all_structures)))

# Validation: Check coverage
coverage_pct <- (nrow(all_structures) / nrow(major_isoforms)) * 100
cat(sprintf("  Coverage: %.2f%% (%d / %d major isoforms)\n",
            coverage_pct, nrow(all_structures), nrow(major_isoforms)))

if (coverage_pct < 95) {
  cat("  ⚠ WARNING: Coverage is less than 95%\n")

  missing_isoforms <- major_isoforms %>%
    filter(!isoform_id %in% all_structures$isoform_id)

  cat(sprintf("  Missing isoforms: %d\n", nrow(missing_isoforms)))
  cat("  First 10 missing:\n")
  print(head(missing_isoforms, 10))
}

# Validation: Check for duplicates
duplicates <- all_structures %>%
  group_by(isoform_id) %>%
  filter(n() > 1)

if (nrow(duplicates) > 0) {
  cat("  ⚠ WARNING: Found duplicate isoform IDs\n")
  cat(sprintf("  Duplicates: %d\n", n_distinct(duplicates$isoform_id)))
} else {
  cat("  ✓ No duplicate isoform IDs\n")
}

# Summary statistics
cat("\nSummary statistics:\n")
cat(sprintf("  Mean exons per isoform: %.1f\n", mean(all_structures$n_exons)))
cat(sprintf("  Median exons per isoform: %.0f\n", median(all_structures$n_exons)))
cat(sprintf("  Range: %d - %d exons\n", min(all_structures$n_exons), max(all_structures$n_exons)))

# ==============================================================================
# 5. Save Output
# ==============================================================================

cat("\nSaving output...\n")
saveRDS(all_structures, "data/isoform_structures.rds")
cat("  ✓ data/isoform_structures.rds\n")

cat("\n✓ Step 2 complete\n")
cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Major isoforms: %d\n", nrow(major_isoforms)))
cat(sprintf("Structures extracted: %d (%.2f%%)\n", nrow(all_structures), coverage_pct))
cat(sprintf("GENCODE: %d\n", nrow(gencode_structures)))
if (nrow(pacbio_structures) > 0) {
  cat(sprintf("PacBio: %d\n", nrow(pacbio_structures)))
}
cat(sprintf("Mean exons per isoform: %.1f\n", mean(all_structures$n_exons)))
cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")
