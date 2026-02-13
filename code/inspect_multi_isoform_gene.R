#!/usr/bin/env Rscript
# Inspect a gene with multiple isoforms

library(tidyverse)

# Load the union exon models
union_models <- readRDS("results/isoform_transitions/v3.0_reference_based/union_exon_models_test.rds")

cat("\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("MULTI-ISOFORM GENE INSPECTION\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Find genes with multiple isoforms
n_isoforms <- sapply(union_models, function(g) g$n_isoforms)
multi_iso_genes <- names(union_models)[n_isoforms > 1]

cat("Genes with multiple isoforms:", length(multi_iso_genes), "\n")
cat("Gene IDs:", paste(multi_iso_genes, collapse=", "), "\n\n")

# Pick the first multi-isoform gene (or the one with most isoforms)
if (length(multi_iso_genes) > 0) {
  # Get the gene with most isoforms
  max_isoforms <- which.max(n_isoforms)
  target_gene <- names(union_models)[max_isoforms]

  gene_data <- union_models[[target_gene]]

  cat("═══════════════════════════════════════════════════════════\n")
  cat("Inspecting:", target_gene, "\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat("Gene ID:", gene_data$gene_id, "\n")
  cat("Dominant isoform:", gene_data$dominant_isoform, "\n")
  cat("Number of isoforms:", gene_data$n_isoforms, "\n")
  cat("IR exons removed:", gene_data$n_ir_removed, "\n")
  cat("Total union exons:", length(gene_data$union_exons), "\n\n")

  # Get all isoforms in this gene
  all_isoforms <- unique(unlist(lapply(gene_data$union_exons, function(e) e$variants$isoform_id)))
  cat("Isoforms:\n")
  for (iso in all_isoforms) {
    cat("  -", iso, "\n")
  }
  cat("\n")

  cat("Union Exon Structure:\n")
  cat("════════════════════════════════════════════════════════\n\n")

  for (exon in gene_data$union_exons) {
    variants <- exon$variants
    n_variants <- nrow(variants)

    cat(sprintf("Exon %d (%s): %d variant(s)\n",
                exon$exon_number, exon$exon_type, n_variants))

    # Show each variant with details
    for (i in seq_len(nrow(variants))) {
      v <- variants[i, ]
      width <- v$end - v$start + 1
      cat(sprintf("  [%d] %s (original exon %d)\n",
                  i, v$isoform_id, v$exon_index))
      cat(sprintf("      Position: %s:%d-%d (%s) | Width: %d bp\n",
                  v$seqnames, v$start, v$end, v$strand, width))
    }

    # Analysis of boundary sharing
    starts <- unique(variants$start)
    ends <- unique(variants$end)

    cat(sprintf("  Boundary analysis: %d unique start(s), %d unique end(s)\n",
                length(starts), length(ends)))

    if (n_variants > 1) {
      if (length(starts) == 1 && length(ends) == 1) {
        cat("  → All variants IDENTICAL (same start and end)\n")
      } else if (length(starts) == 1) {
        cat("  → Variants share START boundary (alternative 3' splice sites)\n")
        cat(sprintf("    Shared start: %d\n", starts[1]))
        cat("    Alternative ends:", paste(sort(ends), collapse=", "), "\n")
      } else if (length(ends) == 1) {
        cat("  → Variants share END boundary (alternative 5' splice sites)\n")
        cat("    Alternative starts:", paste(sort(starts), collapse=", "), "\n")
        cat(sprintf("    Shared end: %d\n", ends[1]))
      } else {
        cat("  ⚠ WARNING: Variants have DIFFERENT boundaries!\n")
        cat("    Starts:", paste(sort(starts), collapse=", "), "\n")
        cat("    Ends:", paste(sort(ends), collapse=", "), "\n")
      }

      # Check which isoforms have which variants
      cat("  Isoform mapping:\n")
      for (i in seq_len(nrow(variants))) {
        v <- variants[i, ]
        cat(sprintf("    %s uses variant [%d]\n", v$isoform_id, i))
      }
    }

    cat("\n")
  }

  # Additional validation
  cat("Validation:\n")
  cat("════════════════════════════════════════════════════════\n")

  # Check exon numbering
  exon_numbers <- sapply(gene_data$union_exons, function(e) e$exon_number)
  if (all(diff(exon_numbers) == 1) && exon_numbers[1] == 1) {
    cat("✓ Exon numbers are sequential\n")
  } else {
    cat("✗ ERROR: Exon numbers are NOT sequential!\n")
  }

  # Check genomic ordering
  internal_exons <- gene_data$union_exons[sapply(gene_data$union_exons, function(e) e$exon_type == "internal")]
  if (length(internal_exons) > 1) {
    positions <- sapply(internal_exons, function(e) min(e$variants$start))
    if (all(diff(positions) > 0)) {
      cat("✓ Internal exons are in genomic order\n")
    } else {
      cat("✗ ERROR: Internal exons are NOT in genomic order!\n")
      cat("  Positions:", paste(positions, collapse=", "), "\n")
    }
  }

  # Check that all variants in a group share at least one boundary
  all_share_boundary <- TRUE
  for (exon in gene_data$union_exons) {
    if (nrow(exon$variants) > 1) {
      starts <- unique(exon$variants$start)
      ends <- unique(exon$variants$end)
      if (length(starts) > 1 && length(ends) > 1) {
        all_share_boundary <- FALSE
        cat(sprintf("✗ ERROR: Exon %d variants don't share any boundary!\n",
                    exon$exon_number))
      }
    }
  }

  if (all_share_boundary) {
    cat("✓ All variant groups share at least one boundary\n")
  }

} else {
  cat("No multi-isoform genes found in the test set.\n")
}

cat("\n")
