#!/usr/bin/env Rscript
# Extract Exon Structures for Major Isoforms
# Part of v4.0 Expression-Based Event Detection Pipeline
#
# Purpose:
#   Extract exon coordinates for major isoforms from annotation files
#   - GENCODE v49 comprehensive GFF3 (for ENST isoforms)
#   - PacBio GFF (for PB isoforms)
#
# Inputs:
#   - major_isoforms_dmso.rds (list of major isoforms from Step 1)
#   - reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz
#   - isoseq/collapse/merge-collapsed.gff
#
# Outputs:
#   - exon_structures_major_isoforms.rds: Exon coordinates for all major isoforms

library(tidyverse)
library(rtracklayer)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EXTRACT EXON STRUCTURES - v4.0                               ║\n")
cat("║   Major Isoforms from Expression Filtering                     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Directories
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v4.0_reference_based"
results_dir <- file.path(base_dir, "major_isoforms")
ref_dir <- "/Users/petecastaldi/claude_projects/nmd/reference_files"
pacbio_dir <- "/Users/petecastaldi/claude_projects/nmd/isoseq/collapse"

# ============================================================================
# Load Major Isoforms List
# ============================================================================

cat("═══ Loading Major Isoforms ═══\n\n")

major_isoforms <- readRDS(file.path(results_dir, "major_isoforms_dmso.rds"))

# Get unique isoform IDs
major_isoform_ids <- unique(major_isoforms$isoform_id)

cat("Major isoforms to extract:", length(major_isoform_ids), "\n")
cat("Unique genes:", length(unique(major_isoforms$gene_id)), "\n\n")

# Identify ENST vs PB isoforms
enst_ids <- major_isoform_ids[grepl("^ENST", major_isoform_ids)]
pb_ids <- major_isoform_ids[grepl("^PB\\.", major_isoform_ids)]

cat("Isoform sources:\n")
cat("  ENST (GENCODE annotated):", length(enst_ids), "\n")
cat("  PB (PacBio novel):", length(pb_ids), "\n\n")

# ============================================================================
# Load DGEList for Proper Gene ID Mapping
# ============================================================================

cat("═══ Loading Gene ID Mapping from DGEList ═══\n\n")

library(edgeR)
dge <- readRDS("/Users/petecastaldi/claude_projects/nmd/rds/dge_isoform_nofilter_2026.2.7.rds")

# Create isoform_id -> gene_id mapping
# NOTE: IDs already match exactly - no version stripping needed!
isoform_gene_map <- dge$genes %>%
  as_tibble() %>%
  select(isoform_id = txid, gene_id_with_version = gene_id_ens115_sqanti) %>%
  # Strip version from gene_id only (to match exon_structures format from GFF)
  mutate(gene_id = str_remove(gene_id_with_version, "\\.\\d+$")) %>%
  select(isoform_id, gene_id) %>%
  filter(isoform_id %in% major_isoform_ids)

cat("Gene ID mapping created:\n")
cat("  Total mapped:", nrow(isoform_gene_map), "isoforms\n")
cat("  Unique genes:", length(unique(isoform_gene_map$gene_id)), "\n\n")

# ============================================================================
# Extract from GENCODE GFF3 (Comprehensive)
# ============================================================================

cat("═══ Extracting from GENCODE v49 Comprehensive GFF3 ═══\n")
cat("Note: Using comprehensive annotation (includes alt scaffolds)\n\n")

gencode_file <- file.path(ref_dir, "gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz")

if (!file.exists(gencode_file)) {
  stop("GENCODE GFF3 file not found: ", gencode_file)
}

cat("Loading GENCODE GFF3...\n")
gencode_gff <- import(gencode_file)
cat("  Total features:", length(gencode_gff), "\n")

# Extract exons for major isoforms
cat("Extracting exons for ENST isoforms...\n")
gencode_exons <- as.data.frame(gencode_gff) %>%
  filter(type == "exon") %>%
  filter(transcript_id %in% enst_ids) %>%
  select(
    isoform_id = transcript_id,
    gene_id,
    seqnames,
    start,
    end,
    strand
  ) %>%
  as_tibble() %>%
  mutate(
    # Strip version numbers from gene_id to match DGEList format
    gene_id = str_remove(gene_id, "\\.\\d+$")
  )

cat("  GENCODE exons extracted:", nrow(gencode_exons), "\n")
cat("  GENCODE isoforms matched:", length(unique(gencode_exons$isoform_id)), "/", length(enst_ids), "\n")

# Check coverage
enst_missing <- setdiff(enst_ids, unique(gencode_exons$isoform_id))
if (length(enst_missing) > 0) {
  cat("  WARNING:", length(enst_missing), "ENST isoforms not found in GENCODE\n")
  cat("  First 10 missing:", paste(head(enst_missing, 10), collapse = ", "), "\n")
}
cat("\n")

# ============================================================================
# Extract from PacBio GFF
# ============================================================================

cat("═══ Extracting from PacBio GFF ═══\n\n")

pacbio_file <- file.path(pacbio_dir, "merge-collapsed.gff")

if (!file.exists(pacbio_file)) {
  stop("PacBio GFF file not found: ", pacbio_file)
}

cat("Loading PacBio GFF...\n")
pb_gff <- import(pacbio_file)
cat("  Total features:", length(pb_gff), "\n")

# Extract exons
cat("Extracting exons for PB isoforms...\n")
pb_exons <- as.data.frame(pb_gff) %>%
  filter(type == "exon") %>%
  filter(transcript_id %in% pb_ids) %>%
  select(
    isoform_id = transcript_id,
    seqnames,
    start,
    end,
    strand
  ) %>%
  as_tibble() %>%
  # Add gene_id from DGEList mapping (PacBio GFF doesn't have proper gene_id)
  left_join(isoform_gene_map, by = "isoform_id") %>%
  select(isoform_id, gene_id, seqnames, start, end, strand)

cat("  PacBio exons extracted:", nrow(pb_exons), "\n")
cat("  PacBio isoforms matched:", length(unique(pb_exons$isoform_id)), "/", length(pb_ids), "\n")

# Check coverage
pb_missing <- setdiff(pb_ids, unique(pb_exons$isoform_id))
if (length(pb_missing) > 0) {
  cat("  WARNING:", length(pb_missing), "PB isoforms not found in PacBio GFF\n")
  cat("  First 10 missing:", paste(head(pb_missing, 10), collapse = ", "), "\n")
}
cat("\n")

# ============================================================================
# Combine and Process Exon Structures
# ============================================================================

cat("═══ Building Exon Structures ═══\n\n")

# Combine both sources
all_exons <- bind_rows(gencode_exons, pb_exons)

cat("Combined exon data:\n")
cat("  Total exons:", nrow(all_exons), "\n")
cat("  Unique isoforms:", length(unique(all_exons$isoform_id)), "\n")
cat("  Unique genes:", length(unique(all_exons$gene_id)), "\n\n")

# Build exon structures (grouped by isoform)
cat("Processing exon structures...\n")
exon_structures <- all_exons %>%
  group_by(isoform_id, gene_id, seqnames, strand) %>%
  # Sort exons by genomic position (5' to 3' on forward strand)
  arrange(start) %>%
  summarise(
    exon_starts = list(start),
    exon_ends = list(end),
    n_exons = n(),
    .groups = "drop"
  )

cat("  Isoforms with complete exon structures:", nrow(exon_structures), "\n\n")

# Mark first and last exons
cat("Marking first and last exons (transcript order)...\n")

# For each isoform, mark which exons are first/last
# Note: On minus strand, first exon in transcript is last in genomic coordinates
exon_structures_expanded <- all_exons %>%
  group_by(isoform_id) %>%
  arrange(start) %>%
  mutate(
    genomic_exon_number = row_number(),
    # On minus strand, reverse the numbering for transcript order
    transcript_exon_number = if_else(
      strand == "-",
      n() - row_number() + 1,
      row_number()
    ),
    is_first = (transcript_exon_number == 1),
    is_last = (transcript_exon_number == max(transcript_exon_number))
  ) %>%
  ungroup()

cat("  First/last flags assigned\n")
cat("  Verification:\n")
cat("    Isoforms with first exon:", sum(exon_structures_expanded$is_first), "\n")
cat("    Isoforms with last exon:", sum(exon_structures_expanded$is_last), "\n\n")

# ============================================================================
# Final Coverage Check
# ============================================================================

cat("═══ Coverage Verification ═══\n\n")

extracted_isoforms <- unique(exon_structures$isoform_id)
coverage_pct <- 100 * length(extracted_isoforms) / length(major_isoform_ids)

cat("Coverage statistics:\n")
cat("  Major isoforms requested:", length(major_isoform_ids), "\n")
cat("  Isoforms extracted:", length(extracted_isoforms), "\n")
cat("  Coverage:", sprintf("%.2f%%", coverage_pct), "\n\n")

# Breakdown by source
enst_extracted <- sum(grepl("^ENST", extracted_isoforms))
pb_extracted <- sum(grepl("^PB\\.", extracted_isoforms))

cat("By source:\n")
cat("  ENST extracted:", enst_extracted, "/", length(enst_ids),
    sprintf("(%.1f%%)", 100 * enst_extracted / length(enst_ids)), "\n")
cat("  PB extracted:", pb_extracted, "/", length(pb_ids),
    sprintf("(%.1f%%)", 100 * pb_extracted / length(pb_ids)), "\n\n")

# Report missing isoforms
all_missing <- setdiff(major_isoform_ids, extracted_isoforms)
if (length(all_missing) > 0) {
  cat("WARNING: Missing", length(all_missing), "isoforms\n")
  cat("First 20 missing:\n")
  print(head(all_missing, 20))
  cat("\n")
}

# ============================================================================
# Save Results
# ============================================================================

cat("═══ Saving Results ═══\n\n")

# Save exon structures
output_file <- file.path(results_dir, "exon_structures_major_isoforms.rds")
saveRDS(exon_structures, output_file)

cat("Saved: exon_structures_major_isoforms.rds\n")
cat("  Location:", output_file, "\n")
cat("  Dimensions:", nrow(exon_structures), "isoforms\n")
cat("  Columns:\n")
cat("    - isoform_id, gene_id, seqnames, strand\n")
cat("    - exon_starts (list), exon_ends (list)\n")
cat("    - n_exons\n\n")

# Save expanded version with first/last flags (for verification)
expanded_output <- file.path(results_dir, "exon_structures_expanded.rds")
saveRDS(exon_structures_expanded, expanded_output)
cat("Saved: exon_structures_expanded.rds\n")
cat("  Contains: Individual exon rows with is_first/is_last flags\n\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("═══ Summary Statistics ═══\n\n")

exon_count_dist <- exon_structures %>%
  count(n_exons) %>%
  arrange(n_exons)

cat("Exons per isoform distribution:\n")
print(head(exon_count_dist, 20))
cat("\n")

cat("Summary:\n")
cat("  Mean exons per isoform:", round(mean(exon_structures$n_exons), 1), "\n")
cat("  Median exons per isoform:", median(exon_structures$n_exons), "\n")
cat("  Range:", min(exon_structures$n_exons), "-", max(exon_structures$n_exons), "\n\n")

# Strand distribution
strand_dist <- exon_structures %>% count(strand)
cat("Strand distribution:\n")
print(strand_dist)
cat("\n")

# ============================================================================
# Final Summary
# ============================================================================

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EXON STRUCTURE EXTRACTION COMPLETE                           ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("RESULTS:\n")
cat("  Isoforms extracted:", nrow(exon_structures), "\n")
cat("  Coverage:", sprintf("%.2f%%", coverage_pct), "\n")
cat("  Mean exons per isoform:", round(mean(exon_structures$n_exons), 1), "\n\n")

cat("OUTPUT LOCATION:\n")
cat(" ", results_dir, "\n\n")

cat("NEXT STEP:\n")
cat("  Rscript code/build_union_exon_model_expression_based.R\n\n")
