#!/usr/bin/env Rscript
# Examine Unclassifiable Genes
#
# Purpose: Investigate why certain genes were marked as unclassifiable
#          (exons overlap but share neither start nor end positions)

library(tidyverse)

# Load data
results_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v4.0_reference_based/major_isoforms"

cat("Loading data...\n")
exon_structures <- readRDS(file.path(results_dir, "exon_structures_major_isoforms.rds"))
filtered_genes <- read_tsv(file.path(results_dir, "filtered_genes_major.tsv"), show_col_types = FALSE)

# Get unclassifiable genes
unclassifiable <- filtered_genes %>%
  filter(status == "unclassifiable")

cat("\nUnclassifiable genes:", nrow(unclassifiable), "\n\n")

# Sample some genes with different numbers of isoforms
set.seed(42)
sample_genes <- unclassifiable %>%
  group_by(n_isoforms) %>%
  slice_sample(n = 2) %>%
  ungroup() %>%
  arrange(n_isoforms) %>%
  head(10)

cat("Examining", nrow(sample_genes), "example genes:\n\n")

# Function to visualize exon structure
visualize_gene_structure <- function(gene_id) {
  gene_exons <- exon_structures %>%
    filter(gene_id == !!gene_id)

  if (nrow(gene_exons) == 0) {
    cat("  No exon data found\n\n")
    return()
  }

  cat("Gene:", gene_id, "\n")
  cat("Isoforms:", nrow(gene_exons), "\n")
  cat("Chromosome:", unique(gene_exons$seqnames), "\n")
  cat("Strand:", unique(gene_exons$strand), "\n\n")

  # Expand exon structures for comparison
  expanded <- gene_exons %>%
    select(isoform_id, exon_starts, exon_ends, n_exons) %>%
    mutate(
      exon_data = map2(exon_starts, exon_ends, function(starts, ends) {
        tibble(
          start = starts,
          end = ends,
          exon_number = seq_along(starts)
        )
      })
    ) %>%
    select(isoform_id, exon_data, n_exons) %>%
    unnest(exon_data)

  # Show first few isoforms in detail
  cat("Detailed exon coordinates:\n")
  for (i in 1:min(4, nrow(gene_exons))) {
    iso <- gene_exons$isoform_id[i]
    iso_exons <- expanded %>% filter(isoform_id == iso)

    cat(sprintf("\n  %s (%d exons):\n", iso, nrow(iso_exons)))
    for (j in 1:min(5, nrow(iso_exons))) {
      cat(sprintf("    Exon %d: %d - %d (%d bp)\n",
                  j,
                  iso_exons$start[j],
                  iso_exons$end[j],
                  iso_exons$end[j] - iso_exons$start[j] + 1))
    }
    if (nrow(iso_exons) > 5) {
      cat(sprintf("    ... %d more exons\n", nrow(iso_exons) - 5))
    }
  }

  # Find overlapping but incompatible exons
  cat("\n  Checking for boundary mismatches...\n")

  all_boundaries <- expanded %>%
    select(isoform_id, start, end) %>%
    pivot_longer(cols = c(start, end), names_to = "boundary_type", values_to = "position") %>%
    distinct()

  # Check for exons that overlap but don't share boundaries
  exon_pairs <- expand_grid(
    iso1 = unique(expanded$isoform_id)[1:min(3, length(unique(expanded$isoform_id)))],
    iso2 = unique(expanded$isoform_id)[1:min(3, length(unique(expanded$isoform_id)))]
  ) %>%
    filter(iso1 < iso2)

  for (k in 1:min(3, nrow(exon_pairs))) {
    iso1_exons <- expanded %>% filter(isoform_id == exon_pairs$iso1[k])
    iso2_exons <- expanded %>% filter(isoform_id == exon_pairs$iso2[k])

    # Find overlapping regions
    overlaps <- expand_grid(
      exon1 = 1:nrow(iso1_exons),
      exon2 = 1:nrow(iso2_exons)
    ) %>%
      mutate(
        start1 = iso1_exons$start[exon1],
        end1 = iso1_exons$end[exon1],
        start2 = iso2_exons$start[exon2],
        end2 = iso2_exons$end[exon2],
        overlap = (start1 <= end2) & (end1 >= start2),
        same_start = start1 == start2,
        same_end = end1 == end2,
        incompatible = overlap & !same_start & !same_end
      ) %>%
      filter(incompatible)

    if (nrow(overlaps) > 0) {
      cat(sprintf("\n  ⚠️  Found %d incompatible overlaps between %s and %s:\n",
                  nrow(overlaps),
                  exon_pairs$iso1[k],
                  exon_pairs$iso2[k]))
      for (m in 1:min(3, nrow(overlaps))) {
        cat(sprintf("      Exon %d (%d-%d) overlaps Exon %d (%d-%d) but boundaries differ\n",
                    overlaps$exon1[m], overlaps$start1[m], overlaps$end1[m],
                    overlaps$exon2[m], overlaps$start2[m], overlaps$end2[m]))
      }
    }
  }

  cat("\n" %+% strrep("=", 80) %+% "\n\n")
}

# Examine each sample gene
for (i in 1:nrow(sample_genes)) {
  visualize_gene_structure(sample_genes$gene_id[i])
}

# Summary statistics
cat("\n\n=== SUMMARY STATISTICS ===\n\n")

# Distribution of isoform counts among unclassifiable genes
isoform_dist <- unclassifiable %>%
  count(n_isoforms) %>%
  arrange(n_isoforms)

cat("Distribution of isoform counts (unclassifiable genes):\n")
print(isoform_dist)

cat("\n\nTotal unclassifiable genes by isoform complexity:\n")
cat(sprintf("  2 isoforms: %d genes (%.1f%%)\n",
            sum(isoform_dist$n[isoform_dist$n_isoforms == 2], na.rm = TRUE),
            100 * sum(isoform_dist$n[isoform_dist$n_isoforms == 2], na.rm = TRUE) / nrow(unclassifiable)))
cat(sprintf("  3-5 isoforms: %d genes (%.1f%%)\n",
            sum(isoform_dist$n[isoform_dist$n_isoforms >= 3 & isoform_dist$n_isoforms <= 5], na.rm = TRUE),
            100 * sum(isoform_dist$n[isoform_dist$n_isoforms >= 3 & isoform_dist$n_isoforms <= 5], na.rm = TRUE) / nrow(unclassifiable)))
cat(sprintf("  6-10 isoforms: %d genes (%.1f%%)\n",
            sum(isoform_dist$n[isoform_dist$n_isoforms >= 6 & isoform_dist$n_isoforms <= 10], na.rm = TRUE),
            100 * sum(isoform_dist$n[isoform_dist$n_isoforms >= 6 & isoform_dist$n_isoforms <= 10], na.rm = TRUE) / nrow(unclassifiable)))

cat("\n")
