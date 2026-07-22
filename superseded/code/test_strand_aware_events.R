#!/usr/bin/env Rscript
# Test strand-aware event detection on sample genes

library(tidyverse)
source("code/detect_events_from_union_model_full.R")

output_dir <- "results/isoform_transitions/v3.0_reference_based"

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   TESTING STRAND-AWARE EVENT DETECTION                       ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Load data
union_models <- readRDS(file.path(output_dir, "union_exon_models_full.rds"))
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))

strand_lookup <- exon_structures %>%
  select(gene_id, strand) %>%
  distinct()

cat("Finding test genes (one plus, one minus strand)...\n")

# Find a plus strand gene and a minus strand gene
plus_genes <- strand_lookup %>% filter(strand == "+") %>% pull(gene_id)
minus_genes <- strand_lookup %>% filter(strand == "-") %>% pull(gene_id)

# Get genes that have union models
plus_test <- intersect(plus_genes, names(union_models))[1]
minus_test <- intersect(minus_genes, names(union_models))[1]

cat("  Plus strand test gene:", plus_test, "\n")
cat("  Minus strand test gene:", minus_test, "\n\n")

# Test function on both genes
for (test_gene in c(plus_test, minus_test)) {
  cat("═══ Testing gene:", test_gene, "═══\n")

  gene_strand <- strand_lookup %>%
    filter(gene_id == !!test_gene) %>%
    pull(strand)

  cat("  Strand:", gene_strand, "\n")

  union_model <- union_models[[test_gene]]
  cat("  Union exons:", union_model$n_union_exons, "\n")
  cat("  Isoforms:", union_model$n_isoforms, "\n")

  # Get first two isoforms
  iso_ids <- unique(union_model$union_exons[[1]]$variants$isoform_id)

  if (length(iso_ids) >= 2) {
    iso_A <- iso_ids[1]
    iso_B <- iso_ids[2]

    cat("\n  Comparing:\n")
    cat("    A:", iso_A, "\n")
    cat("    B:", iso_B, "\n\n")

    # Run event detection
    events <- detect_events_from_union(union_model$union_exons, iso_A, iso_B, gene_strand)

    if (nrow(events) > 0) {
      cat("  Events detected:\n")
      print(events %>% select(event_type, exon_number, direction, detail))
    } else {
      cat("  No events (identical isoforms)\n")
    }
  } else {
    cat("  Only one isoform, skipping\n")
  }

  cat("\n")
}

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   TEST COMPLETE                                               ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("If events detected:\n")
cat("  1. Check that Alt_TSS/Alt_TES directions make biological sense\n")
cat("  2. Check that A5SS/A3SS are correctly assigned (strand-dependent)\n")
cat("  3. Verify SE gain/loss is correct\n\n")
