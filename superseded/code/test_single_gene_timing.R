#!/usr/bin/env Rscript
# Test processing time for a single gene

library(tidyverse, warn.conflicts = FALSE)

output_dir <- "results/isoform_transitions/v3.0_reference_based"

# Load helper functions from main script
source_code <- readLines("code/detect_events_from_union_model_full.R")
# Extract just the helper functions (lines 20-219)
helper_start <- which(grepl("# Helper Functions", source_code))[1]
helper_end <- which(grepl("cat\\(\"Helper functions loaded", source_code))[1] - 1
eval(parse(text = source_code[helper_start:helper_end]))

# Load data
union_models <- readRDS(file.path(output_dir, "union_exon_models_full.rds"))
isoforms_to_process <- readRDS(file.path(output_dir, "isoforms_for_union_model.rds"))
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))

strand_lookup <- exon_structures %>%
  select(isoform_id, strand) %>%
  distinct()

genes_with_models <- names(union_models)
isoforms_filtered <- isoforms_to_process %>%
  filter(gene_id %in% genes_with_models)

# Test SUOX
cat("Testing SUOX processing time...\n\n")

start_time <- Sys.time()
gene_id <- "SUOX"
union_model <- union_models[[gene_id]]

gene_isoforms <- isoforms_filtered %>%
  filter(gene_id == !!gene_id) %>%
  pull(isoform_id)

cat("Isoforms:", length(gene_isoforms), "\n")
cat("Union exons:", union_model$n_union_exons, "\n")

gene_strand <- strand_lookup %>%
  filter(isoform_id == gene_isoforms[1]) %>%
  pull(strand) %>%
  .[1]

cat("Strand:", gene_strand, "\n\n")

transitions_for_gene <- expand_grid(
  isoform_A = gene_isoforms,
  isoform_B = gene_isoforms
) %>%
  filter(isoform_A != isoform_B)

cat("Transitions to process:", nrow(transitions_for_gene), "\n\n")

transition_results <- map_dfr(seq_len(nrow(transitions_for_gene)), function(i) {
  if (i %% 10 == 0) cat("  Transition", i, "/", nrow(transitions_for_gene), "\n")

  iso_A <- transitions_for_gene$isoform_A[i]
  iso_B <- transitions_for_gene$isoform_B[i]

  events <- tryCatch({
    detect_events_from_union(union_model$union_exons, iso_A, iso_B, gene_strand)
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    tibble(event_type = character(),
           exon_number = integer(),
           direction = character(),
           detail = character())
  })

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

cat("\n")
cat("═══════════════════════════════════════\n")
cat("SUOX Processing Complete\n")
cat("═══════════════════════════════════════\n")
cat("Time elapsed:", sprintf("%.2f", elapsed), "seconds\n")
cat("Transitions processed:", nrow(transition_results), "\n")
cat("Rate:", sprintf("%.1f", nrow(transition_results) / elapsed), "transitions/sec\n")
cat("\nExtrapolated full dataset time:\n")

# Estimate based on SUOX
avg_transitions_per_gene <- nrow(transition_results)  # SUOX has 9 isoforms = 72 transitions
time_per_gene <- elapsed
total_genes <- length(union_models)
total_time_hours <- (total_genes * time_per_gene) / 3600

cat("  Genes:", total_genes, "\n")
cat("  Est. time per gene:", sprintf("%.2f", time_per_gene), "sec\n")
cat("  Est. total time:", sprintf("%.1f", total_time_hours), "hours\n")
