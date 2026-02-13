#!/usr/bin/env Rscript
# Test script to verify GTF coverage fix
# Tests switching from primary assembly GTF to comprehensive GFF3

library(tidyverse)
library(rtracklayer)

output_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions"

cat("=== Testing GTF Coverage Fix ===\n\n")

# Load the filtered major isoforms
cat("Loading major isoforms...\n")
major_isoforms_filtered <- readRDS(file.path(output_dir, "major_isoforms_filtered_dmso.rds"))

major_isoform_ids <- unique(major_isoforms_filtered$isoform_id)
cat("Major isoforms to extract:", length(major_isoform_ids), "\n")

# Identify ENST vs PB IDs
enst_ids <- major_isoform_ids[grepl("^ENST", major_isoform_ids)]
pb_ids <- major_isoform_ids[grepl("^PB\\.", major_isoform_ids)]

cat("  ENST (annotated):", length(enst_ids), "\n")
cat("  PB (novel):", length(pb_ids), "\n\n")

# === Test OLD approach (primary assembly GTF) ===
cat("=== OLD APPROACH: Primary Assembly GTF ===\n")
old_gtf_path <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.chrnamesedited.gtf"

if (file.exists(old_gtf_path)) {
  cat("Loading:", old_gtf_path, "\n")
  old_gencode_gtf <- import(old_gtf_path)

  old_gencode_exons <- as.data.frame(old_gencode_gtf) %>%
    filter(type == "exon") %>%
    filter(transcript_id %in% enst_ids) %>%
    select(
      isoform_id = transcript_id,
      gene_id,
      seqnames,
      start,
      end,
      strand
    )

  old_matched <- length(unique(old_gencode_exons$isoform_id))
  old_coverage <- old_matched / length(enst_ids)

  cat("OLD - GENCODE exons extracted:", nrow(old_gencode_exons), "\n")
  cat("OLD - GENCODE isoforms matched:", old_matched, "/", length(enst_ids), "\n")
  cat("OLD - Coverage:", sprintf("%.2f%%", 100 * old_coverage), "\n\n")

  # Identify missing transcripts
  old_missing <- setdiff(enst_ids, unique(old_gencode_exons$isoform_id))
  cat("OLD - Missing transcripts:", length(old_missing), "\n")
  if (length(old_missing) > 0) {
    cat("First 10 missing:", paste(head(old_missing, 10), collapse = ", "), "\n\n")
  }
} else {
  cat("OLD GTF file not found (expected if already replaced)\n\n")
  old_matched <- NA
  old_coverage <- NA
}

# === Test NEW approach (comprehensive GFF3) ===
cat("=== NEW APPROACH: Comprehensive GFF3 ===\n")
new_gff_path <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz"

cat("Loading:", new_gff_path, "\n")
new_gencode_gff <- import(new_gff_path)

new_gencode_exons <- as.data.frame(new_gencode_gff) %>%
  filter(type == "exon") %>%
  filter(transcript_id %in% enst_ids) %>%
  select(
    isoform_id = transcript_id,
    gene_id,
    seqnames,
    start,
    end,
    strand
  )

new_matched <- length(unique(new_gencode_exons$isoform_id))
new_coverage <- new_matched / length(enst_ids)

cat("NEW - GENCODE exons extracted:", nrow(new_gencode_exons), "\n")
cat("NEW - GENCODE isoforms matched:", new_matched, "/", length(enst_ids), "\n")
cat("NEW - Coverage:", sprintf("%.2f%%", 100 * new_coverage), "\n\n")

# Check if any are still missing
new_missing <- setdiff(enst_ids, unique(new_gencode_exons$isoform_id))
cat("NEW - Missing transcripts:", length(new_missing), "\n")
if (length(new_missing) > 0) {
  cat("First 10 missing:", paste(head(new_missing, 10), collapse = ", "), "\n\n")
}

# === Extract from PacBio GFF (unchanged) ===
cat("=== PacBio GFF (Unchanged) ===\n")
pb_gff <- import("/Users/petecastaldi/claude_projects/nmd/isoseq/collapse/merge-collapsed.gff")

pb_exons <- as.data.frame(pb_gff) %>%
  filter(type == "exon") %>%
  filter(transcript_id %in% pb_ids) %>%
  select(
    isoform_id = transcript_id,
    gene_id = any_of(c("gene_id", "gene")),
    seqnames,
    start,
    end,
    strand
  )

pb_matched <- length(unique(pb_exons$isoform_id))
pb_coverage <- pb_matched / length(pb_ids)

cat("PacBio exons extracted:", nrow(pb_exons), "\n")
cat("PacBio isoforms matched:", pb_matched, "/", length(pb_ids), "\n")
cat("PacBio Coverage:", sprintf("%.2f%%", 100 * pb_coverage), "\n\n")

# === Combined Coverage ===
cat("=== COMBINED COVERAGE ===\n")
total_isoforms <- length(major_isoform_ids)

if (!is.na(old_coverage)) {
  old_total_matched <- old_matched + pb_matched
  old_total_coverage <- old_total_matched / total_isoforms
  cat("OLD TOTAL: ", old_total_matched, "/", total_isoforms,
      sprintf(" (%.2f%%)", 100 * old_total_coverage), "\n")
}

new_total_matched <- new_matched + pb_matched
new_total_coverage <- new_total_matched / total_isoforms
cat("NEW TOTAL: ", new_total_matched, "/", total_isoforms,
    sprintf(" (%.2f%%)", 100 * new_total_coverage), "\n\n")

# === Summary ===
cat("=== SUMMARY ===\n")
if (!is.na(old_coverage)) {
  improvement_enst <- new_matched - old_matched
  improvement_pct <- (new_coverage - old_coverage) * 100
  cat("ENST Improvement: +", improvement_enst, " transcripts (",
      sprintf("+%.2f%%", improvement_pct), ")\n")
}

if (new_total_coverage >= 0.9999) {
  cat("\n✅ SUCCESS: 100% coverage achieved!\n")
} else {
  cat("\n⚠️  Coverage improved but not yet 100%\n")
  cat("   Still missing:", total_isoforms - new_total_matched, "isoforms\n")
}

cat("\n=== Test Complete ===\n")
