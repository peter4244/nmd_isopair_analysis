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
#   --source oarfish (default):
#     - data/expression_data.rds (major isoforms list)
#     - reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz
#     - isoseq/collapse/merge-collapsed.gff (PacBio isoforms)
#   --source isocall:
#     - data/isocall/expression_data.rds
#     - isocall GTF (single source for all isoforms)
#
# Output:
#   - data/isoform_structures.rds (oarfish) or data/isocall/isoform_structures.rds (isocall)
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

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse --source (default: oarfish)
source_type <- "oarfish"
if ("--source" %in% args) {
  src_idx <- which(args == "--source")
  if (src_idx < length(args)) {
    source_type <- args[src_idx + 1]
    if (!source_type %in% c("oarfish", "isocall")) {
      stop("--source must be 'oarfish' or 'isocall'")
    }
  }
}

# Set paths based on source
data_dir <- if (source_type == "isocall") "data/isocall" else "data"
output_file <- file.path(data_dir, "isoform_structures.rds")

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
cat(sprintf("║   STEP 2: Extract Isoform Structures [source: %s]%s║\n",
            source_type, strrep(" ", 14 - nchar(source_type))))
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ==============================================================================
# 0. Validate Inputs
# ==============================================================================

cat("Validating inputs...\n")
validate_file(file.path(data_dir, "expression_data.rds"))

if (source_type == "isocall") {
  isocall_gtf_file <- "/Users/petecastaldi/claude_projects/nmd/isocall/nmd_lungcells/results/call/nmd_isocall.isoforms.gtf.gz"
  validate_file(isocall_gtf_file)
} else {
  validate_file("/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz")
  validate_file("/Users/petecastaldi/claude_projects/nmd/isoseq/collapse/merge-collapsed.gff")
}
cat("  All required input files found.\n\n")

# ==============================================================================
# 1. Load Major Isoforms List
# ==============================================================================

cat("Loading major isoforms list...\n")
expression_data <- readRDS(file.path(data_dir, "expression_data.rds"))
validate_columns(expression_data, c("isoform_id", "gene_id"), "expression_data")
cat("  Required columns verified.\n")

major_isoforms <- expression_data %>%
  distinct(isoform_id, gene_id)

cat(sprintf("  Loaded %d major isoforms from %d genes\n",
            nrow(major_isoforms), n_distinct(major_isoforms$gene_id)))

if (source_type == "isocall") {
  # ==============================================================================
  # 2. Extract Structures from Isocall GTF (single source for all isoforms)
  # CRITICAL: No version-stripping — isocall IDs match count matrix exactly
  # ==============================================================================

  cat("\nExtracting isoform structures from isocall GTF...\n")
  cat("  Loading GTF...\n")

  isocall_gff <- import(isocall_gtf_file)

  cat("  Filtering to exon features...\n")
  isocall_exons <- isocall_gff %>%
    as_tibble() %>%
    filter(type == "exon") %>%
    select(transcript_id, gene_id, seqnames, start, end, strand) %>%
    # Filter to major isoforms only — direct ID matching, no version stripping
    filter(transcript_id %in% major_isoforms$isoform_id)

  rm(isocall_gff); gc()

  cat(sprintf("  Extracted %d exons for %d isoforms\n",
              nrow(isocall_exons), n_distinct(isocall_exons$transcript_id)))

  # Build nested structure
  cat("  Building nested exon structures...\n")
  all_structures <- isocall_exons %>%
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
    dplyr::rename(isoform_id = transcript_id) %>%
    mutate(n_junctions = n_exons - 1)

  cat(sprintf("  Created structures for %d isoforms\n", nrow(all_structures)))

} else {
  # ==============================================================================
  # 2-3. OARFISH: Extract from GENCODE + PacBio GFFs
  # ==============================================================================

  # Separate GENCODE (ENST) and PacBio (PB) isoforms
  # Strip version numbers for matching with GFF
  gencode_isoforms <- major_isoforms %>%
    filter(str_detect(isoform_id, "^ENST")) %>%
    mutate(isoform_id_no_version = str_replace(isoform_id, "\\.\\d+$", ""))

  pacbio_isoforms <- major_isoforms %>%
    filter(str_detect(isoform_id, "^PB\\."))

  cat(sprintf("  GENCODE isoforms: %d\n", nrow(gencode_isoforms)))
  cat(sprintf("  PacBio isoforms: %d\n", nrow(pacbio_isoforms)))

  # --- GENCODE ---
  cat("\nExtracting GENCODE isoform structures...\n")
  cat("  Loading GENCODE GFF3 (this may take several minutes)...\n")

  gencode_gff <- import("/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz")

  cat("  Filtering to exon features...\n")
  gencode_exons <- gencode_gff %>%
    as_tibble() %>%
    filter(type == "exon") %>%
    select(transcript_id, gene_id, seqnames, start, end, strand) %>%
    mutate(
      transcript_id_no_version = str_replace(transcript_id, "\\.\\d+$", ""),
      gene_id = str_replace(gene_id, "\\.\\d+$", "")
    ) %>%
    filter(transcript_id_no_version %in% gencode_isoforms$isoform_id_no_version) %>%
    select(-transcript_id_no_version) %>%
    dplyr::rename(isoform_id = transcript_id)

  cat(sprintf("  Extracted %d exons for %d GENCODE isoforms\n",
              nrow(gencode_exons), n_distinct(gencode_exons$isoform_id)))

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
    ) %>%
    mutate(n_junctions = n_exons - 1)

  cat(sprintf("  Created structures for %d GENCODE isoforms\n",
              nrow(gencode_structures)))

  # --- PacBio ---
  if (nrow(pacbio_isoforms) > 0) {
    cat("\nExtracting PacBio isoform structures...\n")
    cat("  Loading PacBio GFF...\n")

    pacbio_gff <- import("/Users/petecastaldi/claude_projects/nmd/isoseq/collapse/merge-collapsed.gff")

    cat("  Filtering to exon features...\n")
    pacbio_exons <- pacbio_gff %>%
      as_tibble() %>%
      filter(type == "exon") %>%
      select(transcript_id, gene_id, seqnames, start, end, strand) %>%
      filter(transcript_id %in% pacbio_isoforms$isoform_id)

    cat(sprintf("  Extracted %d exons for %d PacBio isoforms\n",
                nrow(pacbio_exons), n_distinct(pacbio_exons$transcript_id)))

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
      dplyr::rename(isoform_id = transcript_id) %>%
      mutate(n_junctions = n_exons - 1)

    cat(sprintf("  Created structures for %d PacBio isoforms\n",
                nrow(pacbio_structures)))
  } else {
    cat("\nNo PacBio isoforms to extract (skipping)\n")
    pacbio_structures <- tibble()
  }

  # Combine
  cat("\nCombining GENCODE and PacBio structures...\n")
  if (nrow(pacbio_structures) > 0) {
    all_structures <- bind_rows(gencode_structures, pacbio_structures)
  } else {
    all_structures <- gencode_structures
  }
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
saveRDS(all_structures, output_file)
cat(sprintf("  ✓ %s\n", output_file))

cat("\n✓ Step 2 complete\n")
cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Major isoforms: %d\n", nrow(major_isoforms)))
cat(sprintf("Structures extracted: %d (%.2f%%)\n", nrow(all_structures), coverage_pct))
if (source_type == "oarfish") {
  cat(sprintf("GENCODE: %d\n", nrow(gencode_structures)))
  if (exists("pacbio_structures") && nrow(pacbio_structures) > 0) {
    cat(sprintf("PacBio: %d\n", nrow(pacbio_structures)))
  }
}
cat(sprintf("Mean exons per isoform: %.1f\n", mean(all_structures$n_exons)))
cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")
