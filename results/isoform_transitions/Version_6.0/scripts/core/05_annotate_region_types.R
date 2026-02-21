#!/usr/bin/env Rscript

################################################################################
# Script 05: Annotate Region Types
################################################################################
#
# Purpose:
#   Annotate each (isoform, union_exon) pair with its region type based on
#   CDS boundaries. Region types include: 5'UTR, CDS, 3'UTR, contains_orf_start,
#   contains_orf_stop, non_coding, and unknown.
#
# Input:
#   - data/union_exons.rds
#   - data/isoform_union_mapping.rds
#   - data/isoform_cds_metadata.rds
#
# Output:
#   - data/isoform_union_exons_annotated.rds
#
################################################################################

library(tidyverse)

# ==============================================================================
# Input Validation Helpers
# ==============================================================================

validate_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required input file not found: %s\nHave you run the prerequisite scripts?", path))
  }
}

validate_columns <- function(df, required_cols, name) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing required columns: %s", name, paste(missing, collapse = ", ")))
  }
}

validate_overlap <- function(ids_a, ids_b, name_a, name_b, min_pct = 10) {
  overlap <- length(intersect(ids_a, ids_b))
  pct <- 100 * overlap / length(ids_a)
  cat(sprintf("  Cross-check: %d/%d %s found in %s (%.1f%%)\n",
              overlap, length(ids_a), name_a, name_b, pct))
  if (pct < min_pct) {
    stop(sprintf("FATAL: Only %.1f%% of %s found in %s. Files may be from different pipeline runs.",
                 pct, name_a, name_b))
  }
}

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 5: Annotate Region Types                              ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ==============================================================================
# 0. Validate Inputs
# ==============================================================================

cat("Validating inputs...\n")
validate_file("data/union_exons.rds")
validate_file("data/isoform_union_mapping.rds")
validate_file("data/isoform_cds_metadata.rds")
cat("  All required input files found.\n\n")

# ==============================================================================
# 1. Load Data
# ==============================================================================

cat("Loading data...\n")
union_exons <- readRDS("data/union_exons.rds")
isoform_mapping <- readRDS("data/isoform_union_mapping.rds")
cds_metadata <- readRDS("data/isoform_cds_metadata.rds")

cat(sprintf("  Union exons: %d\n", nrow(union_exons)))
cat(sprintf("  Isoform mappings: %d\n", nrow(isoform_mapping)))
cat(sprintf("  CDS metadata: %d isoforms\n", nrow(cds_metadata)))

# Column validation
validate_columns(union_exons, c("gene_id", "union_exon_id"), "union_exons")
validate_columns(isoform_mapping, c("gene_id", "isoform_id", "union_exon_id"), "isoform_union_mapping")
validate_columns(cds_metadata, c("isoform_id", "coding_status"), "isoform_cds_metadata")

# Cross-file consistency checks
cat("\nCross-file consistency checks:\n")
validate_overlap(unique(union_exons$gene_id), unique(isoform_mapping$gene_id),
                 "union_exons gene_ids", "isoform_union_mapping")
validate_overlap(unique(isoform_mapping$isoform_id), unique(cds_metadata$isoform_id),
                 "mapping isoform_ids", "cds_metadata")
cat("")

# ==============================================================================
# 2. Apply Region Type Logic
# ==============================================================================

cat("\nAnnotating region types...\n")

# Join mapping with CDS metadata
annotated_mapping <- isoform_mapping %>%
  left_join(
    cds_metadata %>% select(isoform_id, coding_status, cds_start, cds_stop),
    by = "isoform_id"
  )

cat("  Applying region type classification...\n")

# Apply region type logic for each (isoform, union_exon) pair
# Based on relationship between isoform exon coordinates and CDS boundaries
#
# IMPORTANT: cds_start and cds_stop are genomic coordinate ranges (min/max)
# - cds_start = minimum genomic coordinate of CDS (works for both strands)
# - cds_stop = maximum genomic coordinate of CDS (works for both strands)
# This allows strand-independent range checking
#
annotated_mapping <- annotated_mapping %>%
  mutate(
    region_type = case_when(
      # Non-coding or unknown isoforms
      coding_status %in% c("non_coding", "unknown") ~ "non_coding",

      # For coding isoforms, check exon position relative to CDS range
      # Use isoform_exon coordinates (actual exon in this isoform)

      # Completely before CDS range (5' UTR on both strands)
      isoform_exon_end < cds_start ~ "5'UTR",

      # Completely after CDS range (3' UTR on both strands)
      isoform_exon_start > cds_stop ~ "3'UTR",

      # Contains BOTH CDS start and stop boundaries (single-CDS-exon genes)
      isoform_exon_start <= cds_start & isoform_exon_end >= cds_stop ~ "contains_orf_start_stop",

      # Contains CDS start boundary (contains first CDS base)
      isoform_exon_start <= cds_start & isoform_exon_end >= cds_start ~ "contains_orf_start",

      # Contains CDS stop boundary (contains last CDS base)
      isoform_exon_start <= cds_stop & isoform_exon_end >= cds_stop ~ "contains_orf_stop",

      # Completely within CDS range
      isoform_exon_start >= cds_start & isoform_exon_end <= cds_stop ~ "CDS",

      # Fallback
      TRUE ~ "unknown"
    )
  )

cat(sprintf("  Annotated %d (isoform, union_exon) pairs\n", nrow(annotated_mapping)))

# ==============================================================================
# 3. Summary Statistics
# ==============================================================================

cat("\nRegion type distribution:\n")

region_counts <- annotated_mapping %>%
  count(region_type, sort = TRUE) %>%
  mutate(percentage = 100 * n / sum(n))

for (i in 1:nrow(region_counts)) {
  cat(sprintf("  %s: %d (%.1f%%)\n",
              region_counts$region_type[i],
              region_counts$n[i],
              region_counts$percentage[i]))
}

# Breakdown by coding status
cat("\nRegion types by coding status:\n")
status_breakdown <- annotated_mapping %>%
  group_by(coding_status, region_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(coding_status, desc(n))

for (status in unique(status_breakdown$coding_status)) {
  cat(sprintf("\n  %s:\n", status))
  subset <- status_breakdown %>% filter(coding_status == status)
  for (i in 1:nrow(subset)) {
    cat(sprintf("    %s: %d\n", subset$region_type[i], subset$n[i]))
  }
}

# ==============================================================================
# 4. Validation
# ==============================================================================

cat("\nValidation checks...\n")

# Check for any unexpected combinations
unexpected <- annotated_mapping %>%
  filter(
    (coding_status == "non_coding" & !region_type %in% c("non_coding", "unknown")) |
    (coding_status == "coding" & region_type == "non_coding")
  )

if (nrow(unexpected) > 0) {
  cat(sprintf("  ⚠ WARNING: %d unexpected region type assignments\n", nrow(unexpected)))
} else {
  cat("  ✓ No unexpected region type assignments\n")
}

# Check coverage
total_mappings <- nrow(isoform_mapping)
annotated_mappings <- sum(!is.na(annotated_mapping$region_type))

cat(sprintf("  ✓ Annotated %d / %d mappings (%.1f%%)\n",
            annotated_mappings,
            total_mappings,
            100 * annotated_mappings / total_mappings))

# ==============================================================================
# 5. Save Output
# ==============================================================================

cat("\nSaving annotated mappings...\n")
saveRDS(annotated_mapping, "data/isoform_union_exons_annotated.rds")
cat("  ✓ data/isoform_union_exons_annotated.rds\n")

cat("\n✓ Step 5 complete\n")
cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Annotated mappings: %d\n", nrow(annotated_mapping)))
cat(sprintf("Region types: %d distinct types\n", n_distinct(annotated_mapping$region_type)))
cat("\nTop region types:\n")
top_regions <- head(region_counts, 5)
for (i in 1:nrow(top_regions)) {
  cat(sprintf("  %s: %.1f%%\n", top_regions$region_type[i], top_regions$percentage[i]))
}
cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")
