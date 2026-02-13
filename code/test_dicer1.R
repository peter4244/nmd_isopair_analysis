#!/usr/bin/env Rscript
# Test DICER1 processing time (gene 124)

library(tidyverse, warn.conflicts = FALSE)

output_dir <- "results/isoform_transitions/v3.0_reference_based"

# Load helper functions
source_code <- readLines("code/detect_events_from_union_model_full.R")
helper_start <- which(grepl("# Helper Functions", source_code))[1]
helper_end <- which(grepl("cat.*Helper functions loaded", source_code))[1] - 1
eval(parse(text = source_code[helper_start:helper_end]))

# Load data
cat("Loading data...\n")
union_models <- readRDS(file.path(output_dir, "union_exon_models_full.rds"))
isoforms_to_process <- readRDS(file.path(output_dir, "isoforms_for_union_model.rds"))
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))

strand_lookup <- exon_structures %>%
  select(isoform_id, strand) %>%
  distinct()

genes_with_models <- names(union_models)
isoforms_filtered <- isoforms_to_process %>%
  filter(gene_id %in% genes_with_models)

# Test gene 124 (DICER1)
gene_idx <- 124
gene_id <- genes_with_models[gene_idx]
union_model <- union_models[[gene_id]]

cat(sprintf("\nTesting gene %d: %s\n", gene_idx, gene_id))
cat(sprintf("  Isoforms: %d\n", union_model$n_isoforms))
cat(sprintf("  Union exons: %d\n", union_model$n_union_exons))

gene_isoforms <- isoforms_filtered %>%
  filter(gene_id == !!gene_id) %>%
  pull(isoform_id)

gene_strand <- strand_lookup %>%
  filter(isoform_id == gene_isoforms[1]) %>%
  pull(strand) %>%
  .[1]

transitions_for_gene <- expand_grid(
  isoform_A = gene_isoforms,
  isoform_B = gene_isoforms
) %>%
  filter(isoform_A != isoform_B)

cat(sprintf("  Transitions: %d\n\n", nrow(transitions_for_gene)))

cat("Processing first 10 transitions...\n")
start_time <- Sys.time()

transition_results <- map_dfr(seq_len(min(10, nrow(transitions_for_gene))), function(i) {
  if (i %% 5 == 0) cat("  ", i, "/", min(10, nrow(transitions_for_gene)), "\n")

  iso_A <- transitions_for_gene$isoform_A[i]
  iso_B <- transitions_for_gene$isoform_B[i]

  events <- detect_events_from_union(union_model$union_exons, iso_A, iso_B, gene_strand)

  event_counts <- events %>%
    count(event_type, name = "count") %>%
    pivot_wider(names_from = event_type,
                values_from = count,
                values_fill = 0,
                names_prefix = "n_")

  all_event_types <- c("n_Alt_TSS", "n_Alt_TES", "n_SE", "n_A5SS", "n_A3SS", "n_CONST")
  for (col in all_event_types) {
    if (!col %in% names(event_counts)) {
      event_counts[[col]] <- 0
    }
  }

  tibble(
    gene_id = gene_id,
    isoform_A = iso_A,
    isoform_B = iso_B,
    n_union_exons = union_model$n_union_exons,
    event_vector = list(events),
    !!!event_counts
  )
})

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat(sprintf("\nFirst 10 transitions: %.2f seconds\n", elapsed))
cat(sprintf("Rate: %.2f transitions/sec\n", 10 / elapsed))
cat(sprintf("Est. for all %d transitions: %.1f seconds (%.1f min)\n",
            nrow(transitions_for_gene),
            elapsed / 10 * nrow(transitions_for_gene),
            (elapsed / 10 * nrow(transitions_for_gene)) / 60))
