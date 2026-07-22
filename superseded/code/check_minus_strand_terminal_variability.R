#!/usr/bin/env Rscript
# Check if minus strand genes have terminal exon variability >20bp
# This will help determine if low Alt_TSS/Alt_TES on minus strand is genuine or a bug

library(tidyverse)

# Load event detection results (use main file, not checkpoint)
results_file <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v3.0_reference_based/event_vectors_full.rds"
transitions <- readRDS(results_file)

# Load exon structures
exon_file <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v3.0_reference_based/exon_structures_by_isoform_full.rds"
exons <- readRDS(exon_file)

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   MINUS STRAND TERMINAL EXON VARIABILITY CHECK               ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Get strand information from exon structures
gene_strands <- exons %>%
  group_by(gene_id) %>%
  summarise(strand = first(strand), .groups = "drop")

# Add strand to transitions
transitions <- transitions %>%
  left_join(gene_strands, by = "gene_id")

# Get minus strand genes
minus_genes <- transitions %>%
  filter(strand == "-") %>%
  pull(gene_id) %>%
  unique()

cat(sprintf("Minus strand genes in checkpoint: %d\n\n", length(minus_genes)))

# Sample 30 genes for detailed check
set.seed(123)
sample_genes <- sample(minus_genes, min(30, length(minus_genes)))

cat("Checking terminal exon variability for 30 sampled genes...\n\n")

# Function to check terminal variability for one gene
check_gene_terminal_variability <- function(gene_id, exon_data) {
  gene_exons <- filter(exon_data, gene_id == !!gene_id)

  if (nrow(gene_exons) == 0) {
    return(NULL)
  }

  # Get all isoforms
  isoforms <- unique(gene_exons$isoform_id)

  if (length(isoforms) < 2) {
    return(NULL)
  }

  # For minus strand: TSS = end of first exon, TES = start of last exon
  terminal_coords <- gene_exons %>%
    group_by(isoform_id) %>%
    arrange(exon_number) %>%
    summarise(
      first_exon_start = first(start),
      first_exon_end = first(end),     # TSS on minus strand
      last_exon_start = last(start),   # TES on minus strand
      last_exon_end = last(end),
      n_exons = n(),
      .groups = "drop"
    )

  # Calculate TSS variability (first exon END on minus)
  tss_coords <- terminal_coords$first_exon_end
  tss_range <- max(tss_coords) - min(tss_coords)

  # Calculate TES variability (last exon START on minus)
  tes_coords <- terminal_coords$last_exon_start
  tes_range <- max(tes_coords) - min(tes_coords)

  tibble(
    gene_id = gene_id,
    n_isoforms = length(isoforms),
    tss_min = min(tss_coords),
    tss_max = max(tss_coords),
    tss_range = tss_range,
    has_alt_tss = tss_range > 20,
    tes_min = min(tes_coords),
    tes_max = max(tes_coords),
    tes_range = tes_range,
    has_alt_tes = tes_range > 20
  )
}

# Check all sampled genes
variability_results <- map_dfr(sample_genes, function(g) {
  check_gene_terminal_variability(g, exons)
})

# Summary
cat("═══ TERMINAL VARIABILITY SUMMARY ═══\n\n")

n_with_alt_tss <- sum(variability_results$has_alt_tss, na.rm = TRUE)
n_with_alt_tes <- sum(variability_results$has_alt_tes, na.rm = TRUE)

cat(sprintf("Genes with Alt_TSS (>20bp): %d / %d (%.1f%%)\n",
            n_with_alt_tss, nrow(variability_results),
            100 * n_with_alt_tss / nrow(variability_results)))
cat(sprintf("Genes with Alt_TES (>20bp): %d / %d (%.1f%%)\n\n",
            n_with_alt_tes, nrow(variability_results),
            100 * n_with_alt_tes / nrow(variability_results)))

# Show examples of genes that SHOULD have Alt_TSS or Alt_TES
genes_with_variability <- variability_results %>%
  filter(has_alt_tss | has_alt_tes)

if (nrow(genes_with_variability) > 0) {
  cat("═══ GENES WITH TERMINAL VARIABILITY >20bp ═══\n\n")
  print(genes_with_variability, n = 20)

  # Now check if these genes have Alt_TSS/Alt_TES in transitions results
  cat("\n═══ CHECKING IF THESE GENES HAVE Alt_TSS/Alt_TES IN RESULTS ═══\n\n")

  transition_events <- transitions %>%
    filter(gene_id %in% genes_with_variability$gene_id) %>%
    group_by(gene_id) %>%
    summarise(
      n_transitions = n(),
      n_alt_tss = sum(n_alt_tss > 0),
      n_alt_tes = sum(n_alt_tes > 0),
      .groups = "drop"
    )

  comparison <- genes_with_variability %>%
    left_join(transition_events, by = "gene_id") %>%
    mutate(
      tss_detected = coalesce(n_alt_tss, 0L) > 0,
      tes_detected = coalesce(n_alt_tes, 0L) > 0,
      tss_missing = has_alt_tss & !tss_detected,
      tes_missing = has_alt_tes & !tes_detected
    )

  cat("Genes with Alt_TSS variability but NO detection:\n")
  missing_tss <- filter(comparison, tss_missing)
  if (nrow(missing_tss) > 0) {
    print(select(missing_tss, gene_id, tss_range, n_transitions, n_alt_tss), n = 20)
  } else {
    cat("  (none - all detected correctly)\n")
  }

  cat("\nGenes with Alt_TES variability but NO detection:\n")
  missing_tes <- filter(comparison, tes_missing)
  if (nrow(missing_tes) > 0) {
    print(select(missing_tes, gene_id, tes_range, n_transitions, n_alt_tes), n = 20)
  } else {
    cat("  (none - all detected correctly)\n")
  }

  # If we find missing detections, investigate one in detail
  if (nrow(missing_tss) > 0 || nrow(missing_tes) > 0) {
    cat("\n╔════════════════════════════════════════════════════════════════╗\n")
    cat("║   BUG CONFIRMED: Terminal variability exists but not detected ║\n")
    cat("╚════════════════════════════════════════════════════════════════╝\n\n")

    # Pick one gene to investigate
    if (nrow(missing_tss) > 0) {
      test_gene <- missing_tss$gene_id[1]
      cat(sprintf("Example gene with missing Alt_TSS: %s\n", test_gene))
      cat(sprintf("  TSS range: %d bp (should trigger Alt_TSS)\n", missing_tss$tss_range[1]))
      cat(sprintf("  Transitions: %d\n", missing_tss$n_transitions[1]))
      cat(sprintf("  Alt_TSS detected: %d\n\n", coalesce(missing_tss$n_alt_tss[1], 0L)))
    } else if (nrow(missing_tes) > 0) {
      test_gene <- missing_tes$gene_id[1]
      cat(sprintf("Example gene with missing Alt_TES: %s\n", test_gene))
      cat(sprintf("  TES range: %d bp (should trigger Alt_TES)\n", missing_tes$tes_range[1]))
      cat(sprintf("  Transitions: %d\n", missing_tes$n_transitions[1]))
      cat(sprintf("  Alt_TES detected: %d\n\n", coalesce(missing_tes$n_alt_tes[1], 0L)))
    }
  }

} else {
  cat("No genes with terminal variability >20bp found in sample.\n")
  cat("This suggests minus strand genes genuinely lack TSS/TES variability.\n")
}

cat("\n✓ Analysis complete\n")
