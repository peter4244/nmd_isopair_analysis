#!/usr/bin/env Rscript
# Inspect Union Exon Models

library(tidyverse)

# Load the union exon models
union_models <- readRDS("results/isoform_transitions/v3.0_reference_based/union_exon_models_test.rds")

cat("\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("UNION EXON MODEL INSPECTION\n")
cat("═══════════════════════════════════════════════════════════\n\n")

cat("Total genes processed:", length(union_models), "\n\n")

# Function to display a single gene's union exon model
display_gene <- function(gene_data, gene_id) {
  cat("\n╔══════════════════════════════════════════════════════════╗\n")
  cat(sprintf("║ Gene: %-50s ║\n", gene_id))
  cat("╚══════════════════════════════════════════════════════════╝\n\n")

  cat("Gene ID:", gene_data$gene_id, "\n")
  cat("Dominant isoform:", gene_data$dominant_isoform, "\n")
  cat("Number of isoforms:", gene_data$n_isoforms, "\n")
  cat("IR exons removed:", gene_data$n_ir_removed, "\n")
  cat("Total union exons:", length(gene_data$union_exons), "\n\n")

  cat("Union Exon Structure:\n")
  cat("════════════════════════════════════════════════════════\n\n")

  for (exon in gene_data$union_exons) {
    variants <- exon$variants
    n_variants <- nrow(variants)

    cat(sprintf("Exon %d (%s):\n", exon$exon_number, exon$exon_type))
    cat(sprintf("  Variants: %d\n", n_variants))

    # Show each variant
    for (i in seq_len(nrow(variants))) {
      v <- variants[i, ]
      cat(sprintf("    [%d] %s (exon %d): %s:%d-%d (%s)\n",
                  i,
                  v$isoform_id,
                  v$exon_index,
                  v$seqnames,
                  v$start,
                  v$end,
                  v$strand))
    }

    # Check for shared boundaries
    starts <- unique(variants$start)
    ends <- unique(variants$end)

    cat(sprintf("  Unique starts: %d, Unique ends: %d\n",
                length(starts), length(ends)))

    if (n_variants > 1) {
      if (length(starts) == 1 && length(ends) == 1) {
        cat("  → All variants share BOTH boundaries (identical)\n")
      } else if (length(starts) == 1) {
        cat("  → All variants share START boundary\n")
      } else if (length(ends) == 1) {
        cat("  → All variants share END boundary\n")
      } else {
        cat("  ⚠ WARNING: Variants do NOT share any boundary!\n")
      }
    }

    cat("\n")
  }

  # Check exon number ordering
  cat("Validation:\n")
  cat("════════════════════════════════════════════════════════\n")

  exon_numbers <- sapply(gene_data$union_exons, function(e) e$exon_number)
  if (all(diff(exon_numbers) == 1) && exon_numbers[1] == 1) {
    cat("✓ Exon numbers are sequential (1, 2, 3, ...)\n")
  } else {
    cat("✗ ERROR: Exon numbers are NOT sequential!\n")
  }

  # Check genomic ordering for internal exons
  internal_exons <- gene_data$union_exons[sapply(gene_data$union_exons, function(e) e$exon_type == "internal")]
  if (length(internal_exons) > 1) {
    positions <- sapply(internal_exons, function(e) min(e$variants$start))
    if (all(diff(positions) > 0)) {
      cat("✓ Internal exons are in genomic order\n")
    } else {
      cat("✗ ERROR: Internal exons are NOT in genomic order!\n")
    }
  }

  cat("\n")
}

# Display details for 3 genes
genes_to_inspect <- names(union_models)[1:min(3, length(union_models))]

for (gene_id in genes_to_inspect) {
  display_gene(union_models[[gene_id]], gene_id)
}

# Additional Summary Statistics
cat("\n═══════════════════════════════════════════════════════════\n")
cat("SUMMARY STATISTICS\n")
cat("═══════════════════════════════════════════════════════════\n\n")

n_union_exons <- sapply(union_models, function(g) length(g$union_exons))
n_isoforms <- sapply(union_models, function(g) g$n_isoforms)
n_ir <- sapply(union_models, function(g) g$n_ir_removed)

cat("Union exons per gene:\n")
cat("  Min:", min(n_union_exons), "\n")
cat("  Max:", max(n_union_exons), "\n")
cat("  Mean:", round(mean(n_union_exons), 2), "\n")
cat("  Median:", median(n_union_exons), "\n\n")

cat("Isoforms per gene:\n")
cat("  Min:", min(n_isoforms), "\n")
cat("  Max:", max(n_isoforms), "\n")
cat("  Mean:", round(mean(n_isoforms), 2), "\n")
cat("  Median:", median(n_isoforms), "\n\n")

cat("IR exons removed:\n")
cat("  Total:", sum(n_ir), "\n")
cat("  Genes with IR:", sum(n_ir > 0), "\n\n")
