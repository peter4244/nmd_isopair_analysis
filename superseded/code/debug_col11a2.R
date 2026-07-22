#!/usr/bin/env Rscript
# Debug COL11A2 internal exon grouping

library(tidyverse)

output_dir <- "results/isoform_transitions/v3.0_reference_based"

# Load data
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))
isoforms_to_process <- readRDS(file.path(output_dir, "isoforms_for_union_model.rds"))

# Get COL11A2 isoforms
col_isoforms <- isoforms_to_process %>% filter(gene_id == "COL11A2")
cat("COL11A2 isoforms:", nrow(col_isoforms), "\n")
print(col_isoforms$isoform_id)

# Load helper functions with DEBUG
source_code <- readLines("code/build_union_exon_model_full_batched.R")
helper_start <- which(grepl("extract_all_exons <-", source_code))[1]
helper_end <- which(grepl("cat.*Helper functions loaded", source_code))[1] - 1
eval(parse(text = source_code[helper_start:helper_end]))

# Extract exons
col_isoforms_ref <- col_isoforms %>% rename(isoform_ref = isoform_id)
all_exons <- extract_all_exons(col_isoforms_ref, exon_structures)

internal_exons <- all_exons %>% filter(!is_first, !is_last)
cat("\nInternal exons:", nrow(internal_exons), "\n")

# Test grouping with explicit debug
cat("\nTesting group_internal_exons:\n")
cat("First 10 internal exons:\n")
print(internal_exons %>% head(10) %>% select(isoform_id, start, end))

# Manually test the grouping logic for first few iterations
processed <- rep(FALSE, nrow(internal_exons))

for (i in 1:min(5, nrow(internal_exons))) {
  if (processed[i]) {
    cat(sprintf("\nIteration %d: SKIPPED (already processed)\n", i))
    next
  }

  current <- internal_exons[i, ]
  cat(sprintf("\nIteration %d:\n", i))
  cat(sprintf("  Current: %s, start=%d, end=%d\n",
              current$isoform_id, current$start, current$end))

  # Old logic (buggy)
  old_group <- which(internal_exons$start == current$start | internal_exons$end == current$end)

  # New logic (fixed)
  new_group <- which(!processed & (internal_exons$start == current$start | internal_exons$end == current$end))

  cat(sprintf("  Old logic matches: %d rows\n", length(old_group)))
  cat(sprintf("  New logic matches: %d rows\n", length(new_group)))

  if (length(old_group) != length(new_group)) {
    cat("  -> DIFFERENCE DETECTED!\n")
    cat("  Old would include already-processed rows:",
        paste(setdiff(old_group, new_group), collapse=", "), "\n")
  }

  processed[new_group] <- TRUE
}
