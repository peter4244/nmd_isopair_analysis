#!/usr/bin/env Rscript
# Check input data for AARS1 gene

library(tidyverse)

output_dir <- "results/isoform_transitions/v3.0_reference_based"

# Load major isoforms
major_isoforms <- readRDS(file.path(output_dir, "reference_event_vectors_v3.0_filtered.rds"))

# Load exon structures
exon_structures <- readRDS(file.path(output_dir, "../exon_structures_by_isoform.rds"))

cat("\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("INPUT DATA CHECK - AARS1\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Get AARS1 isoforms
aars1_isoforms <- major_isoforms %>% filter(gene_id == "AARS1")

cat("AARS1 isoforms in major_isoforms:\n")
print(aars1_isoforms %>% select(gene_id, isoform_ref))
cat("\nTotal:", nrow(aars1_isoforms), "rows\n\n")

# Show the actual isoform IDs
cat("Isoform IDs:\n")
for (i in seq_len(nrow(aars1_isoforms))) {
  cat(sprintf("  %d: %s\n", i, aars1_isoforms$isoform_ref[i]))
}
cat("\n")

# Get exon structures for these isoforms
cat("Exon structures:\n")
for (i in seq_len(nrow(aars1_isoforms))) {
  iso <- aars1_isoforms$isoform_ref[i]
  exon_data <- exon_structures %>% filter(isoform_id == iso)

  if (nrow(exon_data) > 0) {
    cat(sprintf("\n%s:\n", iso))
    cat(sprintf("  Chromosome: %s, Strand: %s\n", exon_data$seqnames[1], exon_data$strand[1]))
    cat(sprintf("  Number of exons: %d\n", length(exon_data$exon_starts[[1]])))
    cat("  Starts:", paste(head(exon_data$exon_starts[[1]], 5), collapse=", "), "...\n")
    cat("  Ends:", paste(head(exon_data$exon_ends[[1]], 5), collapse=", "), "...\n")
  } else {
    cat(sprintf("\n%s: NO EXON DATA FOUND\n", iso))
  }
}

# Check for duplicates in the input
cat("\n\nChecking for duplicate rows in isoform data:\n")
if (nrow(aars1_isoforms) != length(unique(aars1_isoforms$isoform_ref))) {
  cat("⚠ WARNING: Duplicate isoform IDs found!\n")
  cat("Unique isoforms:", length(unique(aars1_isoforms$isoform_ref)), "\n")
  cat("Total rows:", nrow(aars1_isoforms), "\n\n")

  # Show which are duplicates
  dup_isos <- aars1_isoforms$isoform_ref[duplicated(aars1_isoforms$isoform_ref)]
  cat("Duplicated isoform(s):", paste(unique(dup_isos), collapse=", "), "\n")
} else {
  cat("✓ No duplicate isoform IDs\n")
}

cat("\n")
