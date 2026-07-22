#!/usr/bin/env Rscript
# Display the found event type examples

library(tidyverse)

output_dir <- "results/isoform_transitions/v3.0_reference_based"
examples <- readRDS(file.path(output_dir, "event_type_examples.rds"))

source_code <- readLines("code/detect_events_from_union_model_full.R")
helper_start <- which(grepl("# Helper Functions", source_code))[1]
helper_end <- which(grepl("cat.*Helper functions loaded", source_code))[1] - 1
eval(parse(text = source_code[helper_start:helper_end]))

cat("═══════════════════════════════════════════════════════════════\n")
cat("EVENT TYPE EXAMPLES\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (event_type in c("Alt_TSS", "Alt_TES", "A5SS", "A3SS")) {
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat(event_type, "EXAMPLES\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  if (length(examples[[event_type]]) == 0) {
    cat("  No examples found\n")
    next
  }

  for (i in seq_along(examples[[event_type]])) {
    ex <- examples[[event_type]][[i]]

    cat("───────────────────────────────────────────────────────────────\n")
    cat("Example", i, ":", ex$gene_id, "\n")
    cat("───────────────────────────────────────────────────────────────\n")
    cat("  Strand:", ifelse(ex$strand == "+", "Plus (+)", "Minus (-)"), "\n")
    cat("  Union exons:", ex$model$n_union_exons, "\n")
    cat("  Transition:", ex$isoform_A, "→", ex$isoform_B, "\n\n")

    exons_A <- get_isoform_exons(ex$model$union_exons, ex$isoform_A)
    exons_B <- get_isoform_exons(ex$model$union_exons, ex$isoform_B)

    cat("  Isoform A structure:\n")
    print(exons_A %>% select(exon_number, exon_type, start, end, is_first, is_last) %>% head(8))
    if (nrow(exons_A) > 8) cat("   ", nrow(exons_A) - 8, "more exons...\n")

    cat("\n  Isoform B structure:\n")
    print(exons_B %>% select(exon_number, exon_type, start, end, is_first, is_last) %>% head(8))
    if (nrow(exons_B) > 8) cat("   ", nrow(exons_B) - 8, "more exons...\n")

    cat("\n  All Detected Events:\n")
    print(ex$events)

    cat("\n  ", event_type, "event(s):\n")
    specific <- ex$events %>% filter(event_type == !!event_type)
    print(specific)
    cat("\n")
  }
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("END OF EXAMPLES\n")
cat("═══════════════════════════════════════════════════════════════\n")
