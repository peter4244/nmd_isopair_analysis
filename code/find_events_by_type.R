#!/usr/bin/env Rscript
# Find example transitions for each event type

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

# Get strand for a gene
get_gene_strand <- function(model, isoforms) {
  strand_lookup <- exon_structures %>%
    filter(isoform_id %in% isoforms) %>%
    select(isoform_id, strand) %>%
    distinct()

  strand_lookup %>%
    filter(isoform_id == model$dominant_isoform) %>%
    pull(strand) %>%
    .[1]
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("SEARCHING FOR EXAMPLES OF EACH EVENT TYPE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Filter to genes with reasonable complexity (2-4 isoforms, 3-15 union exons)
candidate_genes <- names(union_models)[vapply(union_models, function(x) {
  x$n_isoforms >= 2 && x$n_isoforms <= 4 &&
  x$n_union_exons >= 3 && x$n_union_exons <= 15
}, FUN.VALUE = logical(1))]

cat("Candidate genes:", length(candidate_genes), "\n")
cat("Searching through genes to find event examples...\n\n")

# Process more genes to find rare events like Alt_TSS/Alt_TES
sample_genes <- head(candidate_genes, 2000)

# Storage for examples
examples <- list(
  Alt_TSS = list(),
  Alt_TES = list(),
  A5SS = list(),
  A3SS = list(),
  SE = list(),
  CONST = list()
)

# Track how many we've found
found_counts <- c(Alt_TSS = 0, Alt_TES = 0, A5SS = 0, A3SS = 0, SE = 0, CONST = 0)
target_per_type <- 2  # Find 2 examples of each type

for (gene_id in sample_genes) {
  # Stop if we have enough examples of all types
  if (all(found_counts >= target_per_type)) {
    break
  }

  model <- union_models[[gene_id]]
  isoforms <- get_all_isoforms(model)

  if (length(isoforms) < 2) next

  gene_strand <- get_gene_strand(model, isoforms)
  if (is.na(gene_strand)) next

  # Check first transition
  isoform_A <- isoforms[1]
  isoform_B <- isoforms[2]

  # Detect events
  events <- tryCatch({
    detect_events_from_union(model$union_exons, isoform_A, isoform_B, gene_strand)
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(events) || nrow(events) == 0) next

  # Check which event types are present
  event_types_present <- unique(events$event_type)

  # Store examples for each event type we need
  for (event_type in event_types_present) {
    if (found_counts[event_type] < target_per_type) {
      examples[[event_type]][[length(examples[[event_type]]) + 1]] <- list(
        gene_id = gene_id,
        isoform_A = isoform_A,
        isoform_B = isoform_B,
        model = model,
        strand = gene_strand,
        events = events
      )
      found_counts[event_type] <- found_counts[event_type] + 1
    }
  }
}

cat("Found examples:\n")
print(found_counts)
cat("\n")

# Display examples for each event type
event_types_to_show <- c("Alt_TSS", "Alt_TES", "A5SS", "A3SS", "SE")

for (event_type in event_types_to_show) {
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat(event_type, "EXAMPLES\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  if (length(examples[[event_type]]) == 0) {
    cat("  No examples found in sample\n")
    next
  }

  for (i in seq_along(examples[[event_type]])) {
    ex <- examples[[event_type]][[i]]

    cat("───────────────────────────────────────────────────────────────\n")
    cat("Example", i, ":", ex$gene_id, "\n")
    cat("───────────────────────────────────────────────────────────────\n")
    cat("  Isoforms:", ex$model$n_isoforms, "\n")
    cat("  Union exons:", ex$model$n_union_exons, "\n")
    cat("  Strand:", ifelse(ex$strand == "+", "Plus (+)", "Minus (-)"), "\n")
    cat("  Transition:", ex$isoform_A, "→", ex$isoform_B, "\n\n")

    # Get exon structures
    exons_A <- get_isoform_exons(ex$model$union_exons, ex$isoform_A)
    exons_B <- get_isoform_exons(ex$model$union_exons, ex$isoform_B)

    cat("  Isoform A structure:\n")
    print(exons_A %>% select(exon_number, exon_type, start, end, is_first, is_last))

    cat("\n  Isoform B structure:\n")
    print(exons_B %>% select(exon_number, exon_type, start, end, is_first, is_last))

    cat("\n  All Detected Events:\n")
    print(ex$events)

    # Highlight the specific event type
    cat("\n  ", event_type, "event(s):\n")
    specific_events <- ex$events %>% filter(event_type == !!event_type)
    print(specific_events)

    cat("\n")
  }
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("SEARCH COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n")
