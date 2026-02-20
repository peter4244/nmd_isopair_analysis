#!/usr/bin/env Rscript

################################################################################
# Script 03: Build Atomic Union Exon Models
################################################################################
#
# Purpose:
#   Create atomic (non-overlapping) union exon segments by splitting at all
#   unique exon boundaries within each gene. Each atomic segment represents the
#   smallest unit that is either fully included or fully excluded in any given
#   isoform. This provides a common coordinate system for comparing isoform
#   structures and is required by the reconstruction pipeline.
#
# Algorithm:
#   For each gene, collect all unique exon start/end coordinates. Classify each
#   boundary as START, END, or BOTH. Split into atomic segments using
#   boundary-type-aware logic, then filter to segments covered by at least one
#   exon. This produces fine-grained segments at every exon boundary transition.
#
# Input:
#   - data/isoform_structures.rds
#
# Output:
#   - data/union_exons.rds (atomic union exon coordinates per gene)
#   - data/isoform_union_mapping.rds (which union exons each isoform uses)
#
################################################################################

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 3: Build Atomic Union Exon Models                     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ==============================================================================
# 1. Load Isoform Structures
# ==============================================================================

cat("Loading isoform structures...\n")
isoform_structures <- readRDS("data/isoform_structures.rds")

cat(sprintf("  Loaded %d isoforms from %d genes\n",
            nrow(isoform_structures),
            n_distinct(isoform_structures$gene_id)))

# ==============================================================================
# 2. Build Atomic Union Exons Per Gene
# ==============================================================================

cat("\nBuilding atomic union exon models...\n")
cat("  Processing genes with progress tracking...\n")

union_exons_list <- list()
isoform_mapping_list <- list()

genes <- unique(isoform_structures$gene_id)
pb <- txtProgressBar(min = 0, max = length(genes), style = 3)

ue_counter <- 1

for (i in seq_along(genes)) {
  gene <- genes[i]

  # Get all isoforms for this gene
  gene_isoforms <- isoform_structures %>%
    filter(gene_id == gene)

  if (nrow(gene_isoforms) == 0) next

  # Get chromosome and strand
  chr <- gene_isoforms$seqnames[1]
  strand <- gene_isoforms$strand[1]

  # Extract all exon coordinates across all isoforms
  all_exons <- map_dfr(1:nrow(gene_isoforms), function(j) {
    iso <- gene_isoforms[j, ]
    tibble(
      isoform_id = iso$isoform_id,
      exon_start = iso$exon_starts[[1]],
      exon_end = iso$exon_ends[[1]],
      exon_number_in_isoform = seq_along(iso$exon_starts[[1]])
    )
  })

  if (nrow(all_exons) == 0) next

  # --------------------------------------------------------------------------
  # Atomic boundary splitting (from build_atomic_union_exons.R)
  # --------------------------------------------------------------------------

  # Collect all unique boundaries and sort
  all_boundaries <- sort(unique(c(all_exons$exon_start, all_exons$exon_end)))

  if (length(all_boundaries) < 2) next

  # Determine boundary types: START, END, or BOTH
  boundary_types <- sapply(all_boundaries, function(b) {
    is_start <- b %in% all_exons$exon_start
    is_end <- b %in% all_exons$exon_end
    if (is_start && is_end) return("BOTH")
    if (is_start) return("START")
    return("END")
  })

  # Create atomic segments using boundary-type-aware splitting:
  # - START boundary (middle): split before it → [..., b-1] then [b, ...]
  # - END boundary (middle): split after it → [..., b] then [b+1, ...]
  # - BOTH boundary: split at it → [..., b] then [b+1, ...]
  # - First boundary: always segment start
  # - Last boundary: always segment end

  n <- length(all_boundaries)
  seg_starts <- c()
  seg_ends <- c()

  current_start <- all_boundaries[1]

  for (k in seq_len(n - 2) + 1) {
    b <- all_boundaries[k]
    if (boundary_types[k] == "START") {
      # Split before this START boundary
      seg_starts <- c(seg_starts, current_start)
      seg_ends <- c(seg_ends, b - 1)
      current_start <- b
    } else if (boundary_types[k] == "END") {
      # Split after this END boundary
      seg_starts <- c(seg_starts, current_start)
      seg_ends <- c(seg_ends, b)
      current_start <- b + 1
    } else {
      # BOTH: position is both an exon-end and an exon-start (different isoforms).
      # Must produce THREE segments: [..., b-1], [b, b], [b+1, ...]
      # so that both the exon ending at b and the exon starting at b
      # can claim position b through strict containment.
      seg_starts <- c(seg_starts, current_start)
      seg_ends <- c(seg_ends, b - 1)
      seg_starts <- c(seg_starts, b)
      seg_ends <- c(seg_ends, b)
      current_start <- b + 1
    }
  }

  # Add final segment
  seg_starts <- c(seg_starts, current_start)
  seg_ends <- c(seg_ends, all_boundaries[n])

  segments <- tibble(
    seg_start = seg_starts,
    seg_end = seg_ends
  ) %>%
    filter(seg_start <= seg_end)  # Drop phantom segments from adjacent boundaries

  # Filter to segments covered by at least one exon
  covered_segments <- segments %>%
    rowwise() %>%
    filter({
      any(
        all_exons$exon_start <= seg_start &
        all_exons$exon_end >= seg_end
      )
    }) %>%
    ungroup()

  if (nrow(covered_segments) == 0) next

  # --------------------------------------------------------------------------
  # Create union exons from covered segments
  # --------------------------------------------------------------------------

  gene_ue_start <- ue_counter

  atomic_exons <- covered_segments %>%
    mutate(
      gene_id = gene,
      seqnames = chr,
      strand = strand,
      union_exon_id = paste0(gene, "_UE", row_number()),
      union_exon_number = row_number(),
      union_exon_start = seg_start,
      union_exon_end = seg_end,
      union_exon_length = seg_end - seg_start + 1
    ) %>%
    select(union_exon_id, union_exon_number, union_exon_start, union_exon_end,
           union_exon_length, gene_id, seqnames, strand)

  ue_counter <- ue_counter + nrow(atomic_exons)
  union_exons_list[[gene]] <- atomic_exons

  # --------------------------------------------------------------------------
  # Create isoform-to-union-exon mapping
  # Uses containment check: atomic UE is fully within isoform exon
  # --------------------------------------------------------------------------

  for (j in 1:nrow(gene_isoforms)) {
    iso <- gene_isoforms[j, ]
    iso_exons <- tibble(
      exon_start = iso$exon_starts[[1]],
      exon_end = iso$exon_ends[[1]],
      exon_number_in_isoform = seq_along(iso$exon_starts[[1]])
    )

    # Find which atomic UEs are contained within each isoform exon
    mapping <- map_dfr(1:nrow(iso_exons), function(k) {
      ex_start <- iso_exons$exon_start[k]
      ex_end <- iso_exons$exon_end[k]

      # Containment: atomic UE fully within this isoform exon
      contained <- atomic_exons %>%
        filter(
          union_exon_start >= ex_start & union_exon_end <= ex_end
        )

      if (nrow(contained) > 0) {
        contained %>%
          mutate(
            isoform_id = iso$isoform_id,
            exon_number_in_isoform = iso_exons$exon_number_in_isoform[k],
            isoform_exon_start = ex_start,
            isoform_exon_end = ex_end
          )
      } else {
        tibble()
      }
    })

    isoform_mapping_list[[paste0(gene, "_", iso$isoform_id)]] <- mapping
  }

  setTxtProgressBar(pb, i)
}

close(pb)

# Combine results
union_exons <- bind_rows(union_exons_list)
isoform_union_mapping <- bind_rows(isoform_mapping_list)

cat(sprintf("\n  Created %d atomic union exons across %d genes\n",
            nrow(union_exons),
            n_distinct(union_exons$gene_id)))

# ==============================================================================
# 3. Validation
# ==============================================================================

cat("\nValidating atomic union exons...\n")

# Check for overlapping union exons within genes (should be 0)
overlaps_check <- union_exons %>%
  group_by(gene_id) %>%
  arrange(union_exon_start) %>%
  mutate(
    next_start = lead(union_exon_start),
    overlaps_next = !is.na(next_start) & union_exon_end >= next_start
  ) %>%
  ungroup()

n_overlaps <- sum(overlaps_check$overlaps_next, na.rm = TRUE)

if (n_overlaps > 0) {
  cat(sprintf("  WARNING: %d overlapping union exons detected\n", n_overlaps))
} else {
  cat("  No overlapping union exons (correct)\n")
}

# Check mapping completeness
n_isoforms_with_mapping <- n_distinct(isoform_union_mapping$isoform_id)
n_isoforms_total <- nrow(isoform_structures)

cat(sprintf("  Mapping created for %d / %d isoforms (%.1f%%)\n",
            n_isoforms_with_mapping,
            n_isoforms_total,
            100 * n_isoforms_with_mapping / n_isoforms_total))

# Summary statistics
cat("\nAtomic union exon statistics:\n")
cat(sprintf("  Total union exons: %d\n", nrow(union_exons)))
cat(sprintf("  Union exons per gene - mean: %.1f, median: %.0f, range: %d-%d\n",
            mean(table(union_exons$gene_id)),
            median(table(union_exons$gene_id)),
            min(table(union_exons$gene_id)),
            max(table(union_exons$gene_id))))
cat(sprintf("  Union exon length - mean: %.0f bp, median: %.0f bp\n",
            mean(union_exons$union_exon_length),
            median(union_exons$union_exon_length)))

cat("\nMapping statistics:\n")
cat(sprintf("  Total mappings: %d\n", nrow(isoform_union_mapping)))
cat(sprintf("  Union exons per isoform - mean: %.1f\n",
            nrow(isoform_union_mapping) / n_isoforms_with_mapping))

# ==============================================================================
# 4. Save Outputs
# ==============================================================================

cat("\nSaving outputs...\n")
saveRDS(union_exons, "data/union_exons.rds")
cat("  data/union_exons.rds\n")

saveRDS(isoform_union_mapping, "data/isoform_union_mapping.rds")
cat("  data/isoform_union_mapping.rds\n")

cat("\nStep 3 complete\n")
cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Genes processed: %d\n", n_distinct(union_exons$gene_id)))
cat(sprintf("Atomic union exons created: %d\n", nrow(union_exons)))
cat(sprintf("Isoforms mapped: %d\n", n_isoforms_with_mapping))
cat(sprintf("Mapping entries: %d\n", nrow(isoform_union_mapping)))
cat(sprintf("Mean union exons per gene: %.1f\n", mean(table(union_exons$gene_id))))
cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")
