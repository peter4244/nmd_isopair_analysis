#!/usr/bin/env Rscript
# Find example isoform transitions for manual verification

library(tidyverse)

output_dir <- "results/isoform_transitions/v3.0_reference_based"

# Load data
union_models <- readRDS(file.path(output_dir, "union_exon_models_full.rds"))
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))

# Load helper functions
source_code <- readLines("code/detect_events_from_union_model_full.R")
helper_start <- which(grepl("# Helper Functions", source_code))[1]
helper_end <- which(grepl("cat.*Helper functions loaded", source_code))[1] - 1
eval(parse(text = source_code[helper_start:helper_end]))

# Configuration
TSS_TES_TOLERANCE <- 20

# Helper function to extract isoform IDs from union exons
get_all_isoforms <- function(model) {
  unique(unlist(lapply(model$union_exons, function(ue) ue$variants$isoform_id)))
}

# Filter genes by number of isoforms and validate data
genes_2_isoforms <- names(union_models)[vapply(union_models, function(x) {
  x$n_isoforms == 2
}, FUN.VALUE = logical(1))]
genes_3_isoforms <- names(union_models)[vapply(union_models, function(x) {
  x$n_isoforms == 3
}, FUN.VALUE = logical(1))]

# Also filter by reasonable union exon count (not too complex)
genes_2_simple <- genes_2_isoforms[vapply(genes_2_isoforms, function(g) {
  n_exons <- union_models[[g]]$n_union_exons
  n_exons >= 3 & n_exons <= 10
}, FUN.VALUE = logical(1))]

genes_3_simple <- genes_3_isoforms[vapply(genes_3_isoforms, function(g) {
  n_exons <- union_models[[g]]$n_union_exons
  n_exons >= 3 & n_exons <= 10
}, FUN.VALUE = logical(1))]

cat("═══════════════════════════════════════════════════════════════\n")
cat("GENES WITH 2 ISOFORMS (3-10 union exons)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat("Total genes:", length(genes_2_simple), "\n")
cat("First 10 genes:\n")
print(head(genes_2_simple, 10))

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("GENES WITH 3 ISOFORMS (3-10 union exons)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat("Total genes:", length(genes_3_simple), "\n")
cat("First 10 genes:\n")
print(head(genes_3_simple, 10))

# Select a few interesting examples - mix of 2 and 3 isoform genes
example_genes <- c(
  genes_2_simple[1:3],  # First 3 genes with 2 isoforms
  genes_3_simple[1:2]   # First 2 genes with 3 isoforms
)

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("SELECTED EXAMPLES FOR MANUAL VERIFICATION\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (gene_id in example_genes) {
  model <- union_models[[gene_id]]

  # Extract isoform IDs from union exons
  isoforms <- get_all_isoforms(model)

  # Skip if can't find enough isoforms
  if (length(isoforms) < 2) {
    cat("Skipping", gene_id, "- insufficient isoforms\n\n")
    next
  }

  cat("───────────────────────────────────────────────────────────────\n")
  cat("GENE:", gene_id, "\n")
  cat("───────────────────────────────────────────────────────────────\n")
  cat("  Isoforms:", model$n_isoforms, "\n")
  cat("  Union exons:", model$n_union_exons, "\n")
  cat("  Isoform IDs:", paste(isoforms, collapse = ", "), "\n\n")

  # Get strand
  strand_lookup <- exon_structures %>%
    filter(isoform_id %in% isoforms) %>%
    select(isoform_id, strand) %>%
    distinct()
  gene_strand <- strand_lookup %>%
    filter(isoform_id == model$dominant_isoform) %>%
    pull(strand) %>%
    .[1]

  cat("  Strand:", ifelse(gene_strand == "+", "Plus (+)", "Minus (-)"), "\n\n")

  # Show union exon structure
  cat("  Union Exon Structure:\n")
  for (i in seq_along(model$union_exons)) {
    ue <- model$union_exons[[i]]
    cat(sprintf("    Exon %d (%s): %d variants\n",
                ue$exon_number, ue$exon_type, nrow(ue$variants)))
  }

  # For first transition only (to keep output manageable)
  isoform_A <- isoforms[1]
  isoform_B <- isoforms[2]

  cat("\n  Example Transition:", isoform_A, "→", isoform_B, "\n\n")

  # Get exon structures
  exons_A <- get_isoform_exons(model$union_exons, isoform_A)
  exons_B <- get_isoform_exons(model$union_exons, isoform_B)

  cat("  Isoform A structure:\n")
  print(exons_A %>% select(exon_number, exon_type, start, end, is_first, is_last))

  cat("\n  Isoform B structure:\n")
  print(exons_B %>% select(exon_number, exon_type, start, end, is_first, is_last))

  # Detect events
  events <- detect_events_from_union(model$union_exons, isoform_A, isoform_B, gene_strand)

  cat("\n  Detected Events:\n")
  if (nrow(events) > 0) {
    print(events)
    cat("\n  Event Summary:\n")
    event_counts <- events %>% count(event_type)
    print(event_counts)
  } else {
    cat("    No events detected\n")
  }

  cat("\n")
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("END OF EXAMPLES\n")
cat("═══════════════════════════════════════════════════════════════\n")
