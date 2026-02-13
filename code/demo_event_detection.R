#!/usr/bin/env Rscript
# Demonstrate event detection for SUOX

library(tidyverse)

output_dir <- "results/isoform_transitions/v4.0_reference_based"

# Load helper functions from event detection script
source_code <- readLines("code/detect_events_from_union_model_full.R")
helper_start <- which(grepl("# Helper Functions", source_code))[1]
helper_end <- which(grepl("cat.*Helper functions loaded", source_code))[1] - 1
eval(parse(text = source_code[helper_start:helper_end]))

# Configuration (must match event detection script)
TSS_TES_TOLERANCE <- 20  # Minimum bp difference to call Alt_TSS/Alt_TES

# Load data
union_models <- readRDS(file.path(output_dir, "union_exon_models_full.rds"))
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))

suox_model <- union_models[['SUOX']]

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EVENT DETECTION DEMONSTRATION - SUOX                       ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("SUOX Gene:\n")
cat("  Isoforms:", suox_model$n_isoforms, "\n")
cat("  Union exons:", suox_model$n_union_exons, "\n\n")

# Show union model structure
cat("═══════════════════════════════════════════════════════════════\n")
cat("UNION EXON MODEL\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (i in seq_along(suox_model$union_exons)) {
  ue <- suox_model$union_exons[[i]]
  cat(sprintf("Union Exon %d (%s):\n", ue$exon_number, ue$exon_type))
  print(ue$variants %>% select(isoform_id, start, end, is_first, is_last))
  cat("\n")
}

# Get strand info
strand_lookup <- exon_structures %>% select(isoform_id, strand) %>% distinct()
gene_strand <- strand_lookup %>%
  filter(isoform_id == suox_model$dominant_isoform) %>%
  pull(strand) %>%
  .[1]

cat("Gene strand:", ifelse(gene_strand == "+", "Plus (+)", "Minus (-)"), "\n\n")

# Pick an interesting transition
isoform_A <- "ENST00000886434.1"
isoform_B <- "ENST00000886436.1"

cat("═══════════════════════════════════════════════════════════════\n")
cat("EXAMPLE TRANSITION\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Isoform A:", isoform_A, "\n")
cat("Isoform B:", isoform_B, "\n\n")

# Show which union exons each isoform uses
cat("Isoform A structure:\n")
exons_A <- get_isoform_exons(suox_model$union_exons, isoform_A)
print(exons_A)

cat("\nIsoform B structure:\n")
exons_B <- get_isoform_exons(suox_model$union_exons, isoform_B)
print(exons_B)

# Detect events
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("EVENT DETECTION\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

events <- detect_events_from_union(suox_model$union_exons, isoform_A, isoform_B, gene_strand)

cat("Detected", nrow(events), "events:\n\n")
print(events)

# Explain each event
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("EVENT INTERPRETATION\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (i in seq_len(nrow(events))) {
  event <- events[i, ]
  cat(sprintf("%d. %s at union exon %d (%s):\n",
              i, event$event_type, event$exon_number, event$direction))
  cat("   ", event$detail, "\n")

  # Add biological interpretation
  if (event$event_type == "Alt_TSS") {
    cat("   → Different transcription start sites\n")
  } else if (event$event_type == "Alt_TES") {
    cat("   → Different transcript end sites\n")
  } else if (event$event_type == "SE") {
    if (event$direction == "gain") {
      cat("   → Exon present in A, skipped in B\n")
    } else {
      cat("   → Exon skipped in A, present in B\n")
    }
  } else if (event$event_type == "CONST") {
    cat("   → Constitutive exon (identical in both)\n")
  } else if (event$event_type == "A5SS") {
    cat("   → Alternative 5' splice site (different donor)\n")
  } else if (event$event_type == "A3SS") {
    cat("   → Alternative 3' splice site (different acceptor)\n")
  }
  cat("\n")
}

# Summary
cat("═══════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

event_counts <- events %>% count(event_type)
cat("Event type counts:\n")
print(event_counts)

cat("\nTotal base events (excluding CONST):",
    sum(event_counts$n[event_counts$event_type != "CONST"]), "\n")
