#!/usr/bin/env Rscript
# Comprehensive search through ALL genes to find examples of each event type

library(tidyverse)

output_dir <- "results/isoform_transitions/v3.0_reference_based"
union_models <- readRDS(file.path(output_dir, "union_exon_models_full.rds"))
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))

source_code <- readLines("code/detect_events_from_union_model_full.R")
helper_start <- which(grepl("# Helper Functions", source_code))[1]
helper_end <- which(grepl("cat.*Helper functions loaded", source_code))[1] - 1
eval(parse(text = source_code[helper_start:helper_end]))

TSS_TES_TOLERANCE <- 20

get_all_isoforms <- function(model) {
  unique(unlist(lapply(model$union_exons, function(ue) ue$variants$isoform_id)))
}

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

cat("Searching ALL genes for event type examples...\n\n")

examples <- list()
target_per_type <- 2

for (gene_id in names(union_models)) {
  model <- union_models[[gene_id]]

  if (model$n_union_exons > 20 || model$n_isoforms > 4) next

  isoforms <- get_all_isoforms(model)
  if (length(isoforms) < 2) next

  gene_strand <- get_gene_strand(model, isoforms)
  if (is.na(gene_strand)) next

  # Try all pairwise combinations
  for (i in 1:(length(isoforms)-1)) {
    for (j in (i+1):length(isoforms)) {
      isoform_A <- isoforms[i]
      isoform_B <- isoforms[j]

      events <- tryCatch({
        detect_events_from_union(model$union_exons, isoform_A, isoform_B, gene_strand)
      }, error = function(e) NULL)

      if (is.null(events) || nrow(events) == 0) next

      event_types <- unique(events$event_type)

      for (event_type in c("Alt_TSS", "Alt_TES", "A5SS", "A3SS")) {
        if (event_type %in% event_types) {
          key <- paste0(event_type, "_", length(examples[[event_type]]) + 1)
          if (length(examples[[event_type]]) < target_per_type) {
            examples[[event_type]][[key]] <- list(
              gene_id = gene_id,
              isoform_A = isoform_A,
              isoform_B = isoform_B,
              model = model,
              strand = gene_strand,
              events = events
            )
            cat(sprintf("Found %s example in %s (%s strand)\n",
                        event_type, gene_id, ifelse(gene_strand == "+", "plus", "minus")))
          }
        }
      }

      # Stop searching this gene if we have enough examples
      if (all(sapply(c("Alt_TSS", "Alt_TES", "A5SS", "A3SS"),
                     function(t) length(examples[[t]]) >= target_per_type))) {
        break
      }
    }
    if (all(sapply(c("Alt_TSS", "Alt_TES", "A5SS", "A3SS"),
                   function(t) length(examples[[t]]) >= target_per_type))) {
      break
    }
  }

  # Stop if we have enough of all types
  if (all(sapply(c("Alt_TSS", "Alt_TES", "A5SS", "A3SS"),
                 function(t) length(examples[[t]]) >= target_per_type))) {
    cat("\nFound all required examples!\n")
    break
  }
}

cat("\n\nFinal counts:\n")
cat("Alt_TSS:", length(examples$Alt_TSS), "\n")
cat("Alt_TES:", length(examples$Alt_TES), "\n")
cat("A5SS:", length(examples$A5SS), "\n")
cat("A3SS:", length(examples$A3SS), "\n")

# Save examples
saveRDS(examples, "results/isoform_transitions/v3.0_reference_based/event_type_examples.rds")
cat("\nExamples saved to event_type_examples.rds\n")
