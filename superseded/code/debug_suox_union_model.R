#!/usr/bin/env Rscript
# Debug SUOX union model construction

library(tidyverse)

output_dir <- "results/isoform_transitions/v3.0_reference_based"
TSS_TES_TOLERANCE <- 20

# Load helper functions from union model construction script
source_code <- readLines("code/build_union_exon_model_full_batched.R")
helper_start <- which(grepl("extract_all_exons <-", source_code))[1]
helper_end <- which(grepl("cat.*Helper functions loaded", source_code))[1] - 1
eval(parse(text = source_code[helper_start:helper_end]))

# Load data
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))
isoforms_to_process <- readRDS(file.path(output_dir, "isoforms_for_union_model.rds"))
dge <- readRDS("results/rds/dge_isoform_2026.1.20.rds")
cpm_data <- edgeR::cpm(dge)

# Get SUOX isoforms
suox_isoforms <- isoforms_to_process %>% filter(gene_id == "SUOX")
cat("SUOX isoforms:", nrow(suox_isoforms), "\n")
print(suox_isoforms)

# Determine dominant isoform (use mean_cpm from isoforms_to_process)
dominant <- suox_isoforms %>% arrange(desc(mean_cpm)) %>% slice(1) %>% pull(isoform_id)
cat("\nDominant isoform:", dominant, "(mean CPM:", suox_isoforms %>% filter(isoform_id == dominant) %>% pull(mean_cpm), ")\n\n")

# Extract all exons (rename isoform_id to isoform_ref for function compatibility)
suox_isoforms_ref <- suox_isoforms %>% rename(isoform_ref = isoform_id)
all_exons <- extract_all_exons(suox_isoforms_ref, exon_structures)
cat("Total exons extracted:", nrow(all_exons), "\n\n")

# Split by type
first_exons <- all_exons %>% filter(is_first)
last_exons <- all_exons %>% filter(is_last)
internal_exons <- all_exons %>% filter(!is_first, !is_last)

cat("First exons:", nrow(first_exons), "\n")
print(first_exons %>% arrange(isoform_id, start))

cat("\nLast exons:", nrow(last_exons), "\n")
print(last_exons %>% arrange(isoform_id, start))

# Group first exons
cat("\n══════════════════════════════════════════\n")
cat("GROUPING FIRST EXONS\n")
cat("══════════════════════════════════════════\n\n")

first_groups <- group_first_exons(first_exons, dominant)
cat("Number of first exon groups:", length(first_groups), "\n\n")

for (i in seq_along(first_groups)) {
  cat(sprintf("Group %d:\n", i))
  print(first_groups[[i]] %>% select(isoform_id, start, end))
  cat("\n")
}

# Group last exons
cat("══════════════════════════════════════════\n")
cat("GROUPING LAST EXONS\n")
cat("══════════════════════════════════════════\n\n")

last_groups <- group_last_exons(last_exons, dominant)
cat("Number of last exon groups:", length(last_groups), "\n\n")

for (i in seq_along(last_groups)) {
  cat(sprintf("Group %d:\n", i))
  print(last_groups[[i]] %>% select(isoform_id, start, end))
  cat("\n")
}

# Check for duplicates
cat("══════════════════════════════════════════\n")
cat("CHECKING FOR DUPLICATES\n")
cat("══════════════════════════════════════════\n\n")

# Check first exon groups
first_sigs <- sapply(first_groups, function(g) {
  paste(sort(g$isoform_id), collapse="|")
})
cat("First exon group signatures:\n")
print(table(first_sigs))

# Check last exon groups
last_sigs <- sapply(last_groups, function(g) {
  paste(sort(g$isoform_id), collapse="|")
})
cat("\nLast exon group signatures:\n")
print(table(last_sigs))
