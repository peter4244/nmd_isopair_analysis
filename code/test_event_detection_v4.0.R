#!/usr/bin/env Rscript
# Test event detection on small subset of genes (v4.0 bug fixes)

library(tidyverse)

# Configuration
TSS_TES_TOLERANCE <- 20  # Minimum bp difference to call Alt_TSS/Alt_TES
output_dir <- "results/isoform_transitions/v4.0_reference_based"
TEST_GENES <- 10  # Number of genes to test

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EVENT DETECTION TEST - 10 GENES (v4.0 FIXES)              ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Source helper functions
source("code/test_helpers.R", local = TRUE)

cat("Loading data...\n")
union_models <- readRDS(file.path(output_dir, "union_exon_models_full.rds"))
cat(sprintf("  Union models: %d genes\n", length(union_models)))

# Load strand information
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))
strand_lookup <- exon_structures %>%
  select(isoform_id, strand) %>%
  distinct()
cat(sprintf("  Loaded strand information for %d isoforms\n", nrow(strand_lookup)))

# Load filtered isoforms
isoforms_filtered <- readRDS(file.path(output_dir, "isoforms_for_union_model.rds"))
cat(sprintf("  Isoforms in genes with models: %d\n\n", nrow(isoforms_filtered)))

# Filter genes
cat("Filtering genes by union exon count...\n")
genes_before <- length(union_models)
union_models <- union_models[sapply(union_models, function(x) {
  x$n_union_exons > 1 && x$n_union_exons <= 20
})]
cat(sprintf("  After filter (>1 and ≤20 union exons): %d genes\n\n", length(union_models)))

# Filter fusion genes
cat("Filtering fusion genes...\n")
fusion_pattern <- "ENSG[0-9]+_ENSG[0-9]+"
is_fusion <- grepl(fusion_pattern, names(union_models))
fusion_genes <- names(union_models)[is_fusion]
union_models <- union_models[!is_fusion]
cat(sprintf("  After fusion filter: %d genes\n\n", length(union_models)))

# Select test genes - mix of strands
genes_with_models <- names(union_models)

# Get strand for each gene
gene_strands <- sapply(genes_with_models, function(g) {
  gene_isoforms <- isoforms_filtered %>%
    filter(gene_id == g) %>%
    pull(isoform_id)
  if (length(gene_isoforms) == 0) return(NA)

  strand_lookup %>%
    filter(isoform_id == gene_isoforms[1]) %>%
    pull(strand) %>%
    .[1]
})

# Select 5 plus strand and 5 minus strand genes
plus_genes <- names(gene_strands[gene_strands == "+"])[1:5]
minus_genes <- names(gene_strands[gene_strands == "-"])[1:5]
test_genes <- c(plus_genes, minus_genes)

cat(sprintf("Testing %d genes:\n", length(test_genes)))
cat("  Plus strand:", paste(plus_genes, collapse=", "), "\n")
cat("  Minus strand:", paste(minus_genes, collapse=", "), "\n\n")

# Process test genes
all_transitions <- list()

cat("Processing genes...\n")
for (gene_idx in seq_along(test_genes)) {
  gene_id <- test_genes[gene_idx]
  union_model <- union_models[[gene_id]]

  cat(sprintf("[%d/%d] %s (strand: %s)\n", gene_idx, length(test_genes),
              gene_id, gene_strands[gene_id]))

  # Get isoforms for this gene
  gene_isoforms <- isoforms_filtered %>%
    filter(gene_id == !!gene_id) %>%
    pull(isoform_id)

  if (length(gene_isoforms) < 2) next

  # Get strand
  gene_strand <- gene_strands[gene_id]
  if (is.na(gene_strand)) {
    warning(sprintf("No strand info for gene %s, skipping", gene_id))
    next
  }

  # Generate all pairwise transitions
  transitions_for_gene <- expand_grid(
    isoform_A = gene_isoforms,
    isoform_B = gene_isoforms
  ) %>%
    filter(isoform_A != isoform_B)

  # Detect events for each transition
  transition_results <- map_dfr(seq_len(nrow(transitions_for_gene)), function(i) {
    iso_A <- transitions_for_gene$isoform_A[i]
    iso_B <- transitions_for_gene$isoform_B[i]

    events <- tryCatch({
      detect_events_from_union(union_model$union_exons, iso_A, iso_B, gene_strand)
    }, error = function(e) {
      warning(sprintf("Error for %s -> %s: %s", iso_A, iso_B, e$message))
      tibble(event_type = character(),
             exon_number = integer(),
             direction = character(),
             detail = character())
    })

    # Count events by type
    event_counts <- events %>%
      count(event_type, name = "count") %>%
      pivot_wider(names_from = event_type,
                  values_from = count,
                  values_fill = 0,
                  names_prefix = "n_")

    # Ensure all event type columns exist
    all_event_types <- c("n_Alt_TSS", "n_Alt_TES", "n_SE", "n_A5SS", "n_A3SS", "n_CONST")
    for (col in all_event_types) {
      if (!col %in% names(event_counts)) {
        event_counts[[col]] <- 0
      }
    }

    tibble(
      gene_id = gene_id,
      gene_strand = gene_strand,
      isoform_A = iso_A,
      isoform_B = iso_B,
      n_union_exons = union_model$n_union_exons,
      event_vector = list(events),
      !!!event_counts
    )
  })

  all_transitions[[gene_id]] <- transition_results
}

# Combine all results
all_results <- bind_rows(all_transitions)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   TEST RESULTS SUMMARY                                        ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("Total transitions: %d\n", nrow(all_results)))
cat(sprintf("Total genes: %d\n\n", n_distinct(all_results$gene_id)))

# Event distribution
event_summary <- all_results %>%
  summarise(
    transitions = n(),
    with_alt_tss = sum(n_Alt_TSS > 0),
    with_alt_tes = sum(n_Alt_TES > 0),
    with_se = sum(n_SE > 0),
    with_a5ss = sum(n_A5SS > 0),
    with_a3ss = sum(n_A3SS > 0),
    with_const = sum(n_CONST > 0)
  )

cat("Event distribution:\n")
cat(sprintf("  Alt_TSS: %d transitions (%.1f%%)\n",
            event_summary$with_alt_tss,
            100 * event_summary$with_alt_tss / event_summary$transitions))
cat(sprintf("  Alt_TES: %d transitions (%.1f%%)\n",
            event_summary$with_alt_tes,
            100 * event_summary$with_alt_tes / event_summary$transitions))
cat(sprintf("  SE: %d transitions (%.1f%%)\n",
            event_summary$with_se,
            100 * event_summary$with_se / event_summary$transitions))
cat(sprintf("  A5SS: %d transitions (%.1f%%)\n",
            event_summary$with_a5ss,
            100 * event_summary$with_a5ss / event_summary$transitions))
cat(sprintf("  A3SS: %d transitions (%.1f%%)\n",
            event_summary$with_a3ss,
            100 * event_summary$with_a3ss / event_summary$transitions))
cat(sprintf("  CONST: %d transitions (%.1f%%)\n\n",
            event_summary$with_const,
            100 * event_summary$with_const / event_summary$transitions))

# Strand distribution
strand_summary <- all_results %>%
  group_by(gene_strand) %>%
  summarise(
    n_genes = n_distinct(gene_id),
    n_transitions = n(),
    n_alt_tss = sum(n_Alt_TSS > 0),
    n_alt_tes = sum(n_Alt_TES > 0),
    .groups = "drop"
  )

cat("═══ STRAND DISTRIBUTION ═══\n\n")
print(strand_summary)

# Check for minus strand events
minus_events <- filter(strand_summary, gene_strand == "-")
if (nrow(minus_events) > 0 && (minus_events$n_alt_tss > 0 || minus_events$n_alt_tes > 0)) {
  cat("\n✅ SUCCESS: Minus strand Alt_TSS/Alt_TES detected!\n")
} else {
  cat("\n⚠️  WARNING: No minus strand Alt_TSS/Alt_TES detected in test\n")
}

# Save test results
saveRDS(all_results, file.path(output_dir, "test_results_10genes.rds"))
write_tsv(all_results %>% select(-event_vector),
          file.path(output_dir, "test_results_10genes.tsv"))

cat("\n✓ Test complete! Results saved to:\n")
cat(sprintf("  %s/test_results_10genes.rds\n", output_dir))
cat(sprintf("  %s/test_results_10genes.tsv\n", output_dir))
