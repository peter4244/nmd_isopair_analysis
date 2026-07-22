#!/usr/bin/env Rscript
# Union Exon Model Construction - Expression-Based (Major Isoforms)
# Part of v4.0 Major Isoforms Pipeline
#
# Purpose:
#   Build union exon models for genes with major isoforms
#   - Uses expression-filtered isoforms only
#   - Includes v4.0 improvements (all first/last exons)
#
# Inputs:
#   - major_isoforms/exon_structures_major_isoforms.rds
#   - major_isoforms/dominant_isoforms_dmso.rds
#
# Outputs:
#   - major_isoforms/union_exon_models_major.rds
#   - major_isoforms/filtered_genes_major.tsv

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   UNION EXON MODEL CONSTRUCTION - MAJOR ISOFORMS               ║\n")
cat("║   Expression-Based Filtering (v4.0)                            ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Configuration
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v4.0_reference_based"
output_dir <- file.path(base_dir, "major_isoforms")
TSS_TES_TOLERANCE <- 20
BATCH_SIZE <- 500  # Save progress every 500 genes

# ============================================================================
# Helper Functions (v4.0 - All First/Last Exons Included)
# ============================================================================

cat("Loading helper functions...\n\n")

# Extract all exons for a gene's isoforms
# Now uses expanded exon structures with strand-aware is_first/is_last flags
extract_all_exons <- function(gene_id, exon_structures_expanded) {
  all_exons <- exon_structures_expanded %>%
    filter(gene_id == !!gene_id) %>%
    select(isoform_id, start, end, is_first, is_last, transcript_exon_number)

  return(all_exons)
}

# Group first exons by shared start coordinate (within tolerance)
# First exons that overlap AND share start (within TSS_TOLERANCE) are grouped together
group_first_exons <- function(first_exons, tolerance = TSS_TES_TOLERANCE) {
  if (nrow(first_exons) == 0) return(list())

  first_groups <- list()
  processed <- rep(FALSE, nrow(first_exons))

  for (i in seq_len(nrow(first_exons))) {
    if (processed[i]) next

    current_start <- first_exons$start[i]
    current_end <- first_exons$end[i]

    # Find exons that share start within tolerance AND haven't been processed yet
    within_tolerance <- abs(first_exons$start - current_start) <= tolerance & !processed

    # Additionally require overlap
    overlaps <- (first_exons$start <= current_end) & (first_exons$end >= current_start)

    # Group = within tolerance AND overlapping
    same_group <- within_tolerance & overlaps
    same_group[i] <- TRUE  # Include current exon

    if (sum(same_group) > 0) {
      first_groups[[length(first_groups) + 1]] <- first_exons[same_group, ]
      processed[same_group] <- TRUE
    }
  }

  return(first_groups)
}

# Group last exons by shared end coordinate (within tolerance)
# Last exons that overlap AND share end (within TES_TOLERANCE) are grouped together
group_last_exons <- function(last_exons, tolerance = TSS_TES_TOLERANCE) {
  if (nrow(last_exons) == 0) return(list())

  last_groups <- list()
  processed <- rep(FALSE, nrow(last_exons))

  for (i in seq_len(nrow(last_exons))) {
    if (processed[i]) next

    current_start <- last_exons$start[i]
    current_end <- last_exons$end[i]

    # Find exons that share end within tolerance AND haven't been processed yet
    within_tolerance <- abs(last_exons$end - current_end) <= tolerance & !processed

    # Additionally require overlap
    overlaps <- (last_exons$start <= current_end) & (last_exons$end >= current_start)

    # Group = within tolerance AND overlapping
    same_group <- within_tolerance & overlaps
    same_group[i] <- TRUE  # Include current exon

    if (sum(same_group) > 0) {
      last_groups[[length(last_groups) + 1]] <- last_exons[same_group, ]
      processed[same_group] <- TRUE
    }
  }

  return(last_groups)
}

# Group monoexonic exons by exact coordinate match
# Monoexonic exons (is_first AND is_last) with identical (start, end) are grouped together
# This ensures each unique monoexonic exon appears in exactly ONE union exon group
group_monoexonic_exons <- function(monoexonic_exons) {
  if (nrow(monoexonic_exons) == 0) return(list())

  monoexonic_groups <- list()
  processed <- rep(FALSE, nrow(monoexonic_exons))

  for (i in seq_len(nrow(monoexonic_exons))) {
    if (processed[i]) next

    current_start <- monoexonic_exons$start[i]
    current_end <- monoexonic_exons$end[i]

    # Group by EXACT coordinate match (both start and end must match)
    # AND hasn't been processed yet
    same_group <- (monoexonic_exons$start == current_start) &
                  (monoexonic_exons$end == current_end) &
                  !processed
    same_group[i] <- TRUE  # Include current exon

    if (sum(same_group) > 0) {
      monoexonic_groups[[length(monoexonic_groups) + 1]] <- monoexonic_exons[same_group, ]
      processed[same_group] <- TRUE
    }
  }

  return(monoexonic_groups)
}

# NOTE: IR detection removed from union exon model building phase (2026-02-12)
#
# RATIONALE:
# - IR is a RELATIONAL property between isoforms, not an intrinsic exon property
# - An exon is only "IR" relative to another specific isoform
# - Detecting IR during model building was incorrectly removing legitimate isoforms
# - Example: Multiple last exons sharing same TES (48026475) but different starts
#   (48026142, 48026167, 48026169, etc.) are A3SS variants, NOT intron retention
# - These were being flagged as "IR" because longer ones "contain" shorter ones
# - But they should be grouped together (shared boundary) and classified during
#   event detection when comparing specific isoform pairs
#
# NEW APPROACH:
# - Keep ALL unique exons (unique start/end coordinates) in union model
# - Group by shared boundaries (position-based)
# - Detect IR during EVENT DETECTION phase when comparing isoform A vs B
# - This preserves all isoform information and correctly classifies events
#
# TODO: Special case for event detection (future implementation)
# - When multiple exons share same TES but have different 5' starts (A3SS pattern)
# - Need to handle these carefully during event detection
# - May require grouping variants that are within small distance of each other
# - Current example: 48026142-475, 48026167-475, 48026169-475 (2-27 bp differences)
# - These should likely be treated as a single A3SS event group, not separate events

# ============================================================================
# UNIFIED GROUPING APPROACH: Position-based with structural metadata
# ============================================================================
# Union exons are defined purely by coordinate matching (overlap + shared boundary)
# Structural role (first/last) is tracked as metadata within each union exon group
# This eliminates cross-type duplicates and simplifies the model

# Group ALL exons by shared boundaries (start OR end)
# Returns list of exon groups, each representing one union exon
group_exons_by_boundaries <- function(all_exons) {
  if (nrow(all_exons) == 0) return(list())

  exon_groups <- list()
  processed <- rep(FALSE, nrow(all_exons))

  for (i in seq_len(nrow(all_exons))) {
    if (processed[i]) next

    current_start <- all_exons$start[i]
    current_end <- all_exons$end[i]

    # Find all exons that share START or END with current exon
    # AND haven't been processed yet (prevents cascading/duplicates)
    same_group <- ((all_exons$start == current_start) | (all_exons$end == current_end)) & !processed
    same_group[i] <- TRUE  # Include current exon

    if (sum(same_group) > 0) {
      exon_groups[[length(exon_groups) + 1]] <- all_exons[same_group, ]
      processed[same_group] <- TRUE
    }
  }

  return(exon_groups)
}

# Add structural metadata to union exon group
# Tracks which isoforms use this exon as first/last
add_structural_metadata <- function(exon_group) {
  # Get isoforms where this union exon is first
  first_isoforms <- exon_group %>%
    filter(is_first) %>%
    pull(isoform_id) %>%
    unique()

  # Get isoforms where this union exon is last
  last_isoforms <- exon_group %>%
    filter(is_last) %>%
    pull(isoform_id) %>%
    unique()

  # Add metadata columns
  exon_group %>%
    mutate(
      first_in_isoforms = if (length(first_isoforms) > 0) list(first_isoforms) else list(NULL),
      last_in_isoforms = if (length(last_isoforms) > 0) list(last_isoforms) else list(NULL)
    )
}

# Build union model for a single gene
build_union_model <- function(gene_id, n_isoforms, exon_structures_expanded, dominant_isoform_id) {
  # Extract all exons for this gene (uses strand-aware is_first/is_last)
  all_exons <- extract_all_exons(gene_id, exon_structures_expanded)

  if (nrow(all_exons) == 0) {
    return(list(
      gene_id = gene_id,
      n_isoforms = n_isoforms,
      dominant_isoform = dominant_isoform_id,
      union_exons = NULL,
      n_union_exons = 0,
      status = "no_exons"
    ))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # NEW UNIFIED APPROACH: Position-based grouping with structural metadata
  # ═══════════════════════════════════════════════════════════════════════════
  # All exons grouped by coordinate matching (overlap + shared boundary)
  # No type-based separation - monoexonic, first, last, internal all treated equally
  # Structural role tracked as metadata (first_in_isoforms, last_in_isoforms)
  # NO IR REMOVAL - all unique exons are kept; IR will be detected during event detection

  # Group ALL exons by shared boundaries (position-based only)
  # Each unique (start, end) coordinate will appear in exactly one union exon group
  exon_groups <- group_exons_by_boundaries(all_exons)

  # NOTE: Unclassifiable check removed - instead of rejecting entire genes,
  # we'll flag individual unclassifiable events during event detection.
  # This maximizes data usage and allows analysis of classifiable parts of genes.

  # Add structural metadata to each union exon group
  union_exons <- lapply(exon_groups, add_structural_metadata)

  # Assign exon numbers and classify union exon groups
  for (i in seq_along(union_exons)) {
    exon_group <- union_exons[[i]]

    # Determine exon type based on what the union exon contains
    # Since groups can now contain mixed types, classify by presence of structural markers
    has_first <- any(exon_group$is_first)
    has_last <- any(exon_group$is_last)
    all_monoexonic <- all(exon_group$is_first & exon_group$is_last)

    if (all_monoexonic) {
      exon_type <- "monoexonic"
    } else if (has_first && has_last) {
      exon_type <- "mixed"  # Contains both first and last exons from different isoforms
    } else if (has_first) {
      exon_type <- "first_containing"
    } else if (has_last) {
      exon_type <- "last_containing"
    } else {
      exon_type <- "internal"
    }

    union_exons[[i]] <- exon_group %>%
      mutate(
        union_exon_number = i,
        exon_type = exon_type
      )
  }

  return(list(
    gene_id = gene_id,
    n_isoforms = n_isoforms,
    dominant_isoform = dominant_isoform_id,
    union_exons = union_exons,
    n_union_exons = length(union_exons),
    status = "success"
  ))
}

# ============================================================================
# Load Data
# ============================================================================

cat("═══ Loading Data ═══\n\n")

# Load exon structures (EXPANDED version with strand-aware is_first/is_last flags)
cat("Loading exon structures (expanded)...\n")
exon_structures_expanded <- readRDS(file.path(output_dir, "exon_structures_expanded.rds"))
cat("  Total exons:", nrow(exon_structures_expanded), "\n")
cat("  Isoforms:", length(unique(exon_structures_expanded$isoform_id)), "\n")
cat("  Genes:", length(unique(exon_structures_expanded$gene_id)), "\n\n")

# Load dominant isoforms
cat("Loading dominant isoforms...\n")
dominant_isoforms <- readRDS(file.path(output_dir, "dominant_isoforms_dmso.rds"))
cat("  Genes with dominant:", nrow(dominant_isoforms), "\n\n")

# ============================================================================
# Filter Monoexonic Genes
# ============================================================================

cat("═══ Filtering Monoexonic Genes ═══\n\n")

# Identify genes where ALL isoforms have only 1 exon
monoexonic_genes <- exon_structures_expanded %>%
  group_by(gene_id, isoform_id) %>%
  summarise(n_exons = n(), .groups = "drop_last") %>%
  summarise(
    max_exons = max(n_exons),
    all_monoexonic = all(n_exons == 1),
    .groups = "drop"
  ) %>%
  filter(all_monoexonic)

cat("Monoexonic genes (all isoforms have 1 exon):", nrow(monoexonic_genes), "\n")

# Filter out monoexonic genes
exon_structures_expanded <- exon_structures_expanded %>%
  filter(!gene_id %in% monoexonic_genes$gene_id)

cat("After filtering monoexonic genes:\n")
cat("  Exons:", nrow(exon_structures_expanded), "\n")
cat("  Isoforms:", length(unique(exon_structures_expanded$isoform_id)), "\n")
cat("  Genes:", length(unique(exon_structures_expanded$gene_id)), "\n\n")

# ============================================================================
# Process Genes
# ============================================================================

cat("═══ Building Union Models ═══\n\n")

# Group isoforms by gene
genes_to_process <- exon_structures_expanded %>%
  group_by(gene_id) %>%
  summarise(
    n_isoforms = n_distinct(isoform_id),
    .groups = "drop"
  ) %>%
  arrange(gene_id)

cat("Genes to process:", nrow(genes_to_process), "\n")
cat("Batch size:", BATCH_SIZE, "genes\n\n")

# Initialize results storage
union_models <- list()
filtered_genes <- tibble(
  gene_id = character(),
  status = character(),
  n_isoforms = integer()
)

# Process genes
start_time <- Sys.time()
n_genes <- nrow(genes_to_process)

for (i in seq_len(n_genes)) {
  gene_id <- genes_to_process$gene_id[i]
  n_isoforms <- genes_to_process$n_isoforms[i]

  # Get dominant isoform
  dominant_info <- dominant_isoforms %>%
    filter(gene_id == !!gene_id)

  dominant_id <- if (nrow(dominant_info) > 0) {
    dominant_info$dominant_isoform_id[1]
  } else {
    NA_character_
  }

  # Build union model (uses strand-aware expanded exon structures)
  union_model <- build_union_model(gene_id, n_isoforms, exon_structures_expanded, dominant_id)

  # Store result
  if (union_model$status == "success") {
    union_models[[gene_id]] <- union_model
  } else {
    filtered_genes <- bind_rows(
      filtered_genes,
      tibble(
        gene_id = gene_id,
        status = union_model$status,
        n_isoforms = union_model$n_isoforms
      )
    )
  }

  # Progress reporting
  if (i %% 100 == 0 || i == n_genes) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    eta <- (elapsed / i) * (n_genes - i)
    cat(sprintf("[%d/%d] %s | Elapsed: %.1f min | ETA: %.1f min\n",
                i, n_genes, gene_id, elapsed, eta))
  }

  # Checkpoint saving
  if (i %% BATCH_SIZE == 0) {
    cat(sprintf("  [CHECKPOINT] Saving progress at gene %d...\n", i))
    saveRDS(union_models, file.path(output_dir, "union_exon_models_major_checkpoint.rds"))
  }
}

cat("\n")

# ============================================================================
# Post-processing: Filter by Union Exon Count
# ============================================================================

cat("═══ Filtering by Union Exon Count ═══\n\n")

# Calculate union exon counts
union_exon_counts <- map_int(union_models, ~.x$n_union_exons)

cat("Before filtering:\n")
cat("  Genes with union models:", length(union_models), "\n")
cat("  Mean union exons:", round(mean(union_exon_counts), 1), "\n")
cat("  Median union exons:", median(union_exon_counts), "\n\n")

# Filter: Keep genes with 2-20 union exons
genes_to_keep <- names(union_models)[union_exon_counts >= 2 & union_exon_counts <= 20]
genes_to_remove <- names(union_models)[union_exon_counts < 2 | union_exon_counts > 20]

cat("Filtering criteria: 2-20 union exons\n")
cat("  Genes to keep:", length(genes_to_keep), "\n")
cat("  Genes to remove:", length(genes_to_remove), "\n")
cat("    - ≤1 exon:", sum(union_exon_counts < 2), "\n")
cat("    - >20 exons:", sum(union_exon_counts > 20), "\n\n")

# Add removed genes to filtered list
if (length(genes_to_remove) > 0) {
  filtered_genes <- bind_rows(
    filtered_genes,
    tibble(
      gene_id = genes_to_remove,
      status = "excessive_complexity",
      n_isoforms = map_int(union_models[genes_to_remove], ~.x$n_isoforms)
    )
  )
}

# Keep only valid genes
union_models_filtered <- union_models[genes_to_keep]

cat("After filtering:\n")
cat("  Final genes:", length(union_models_filtered), "\n")
cat("  Total filtered:", nrow(filtered_genes), "\n\n")

# ============================================================================
# Check for Fusion Genes
# ============================================================================

cat("═══ Checking for Fusion Genes ═══\n\n")

# Identify fusion genes (pattern: ENSG_ENSG)
fusion_pattern <- "_ENSG[0-9]+"
fusion_genes <- names(union_models_filtered)[grepl(fusion_pattern, names(union_models_filtered))]

if (length(fusion_genes) > 0) {
  cat("Fusion genes identified:", length(fusion_genes), "\n")
  cat("First 10:", paste(head(fusion_genes, 10), collapse = ", "), "\n\n")

  # Remove fusion genes
  union_models_final <- union_models_filtered[!names(union_models_filtered) %in% fusion_genes]

  # Add to filtered list
  filtered_genes <- bind_rows(
    filtered_genes,
    tibble(
      gene_id = fusion_genes,
      status = "fusion_gene",
      n_isoforms = map_int(union_models_filtered[fusion_genes], ~.x$n_isoforms)
    )
  )

  cat("After fusion filter:\n")
  cat("  Final genes:", length(union_models_final), "\n\n")
} else {
  cat("No fusion genes found\n\n")
  union_models_final <- union_models_filtered
}

# ============================================================================
# Save Results
# ============================================================================

cat("═══ Saving Results ═══\n\n")

# Save union models
output_file <- file.path(output_dir, "union_exon_models_major.rds")
saveRDS(union_models_final, output_file)
cat("Saved: union_exon_models_major.rds\n")
cat("  Location:", output_file, "\n")
cat("  Genes:", length(union_models_final), "\n\n")

# Save filtered genes list
if (nrow(filtered_genes) > 0) {
  filtered_file <- file.path(output_dir, "filtered_genes_major.tsv")
  write_tsv(filtered_genes, filtered_file)
  cat("Saved: filtered_genes_major.tsv\n")
  cat("  Genes excluded:", nrow(filtered_genes), "\n")
  cat("  Reasons:\n")
  print(count(filtered_genes, status))
  cat("\n")
}

# Remove checkpoint file
checkpoint_file <- file.path(output_dir, "union_exon_models_major_checkpoint.rds")
if (file.exists(checkpoint_file)) {
  file.remove(checkpoint_file)
  cat("Removed checkpoint file\n\n")
}

# ============================================================================
# Summary Statistics
# ============================================================================

cat("═══ Summary Statistics ═══\n\n")

final_exon_counts <- map_int(union_models_final, ~.x$n_union_exons)

cat("Final dataset:\n")
cat("  Genes:", length(union_models_final), "\n")
cat("  Mean union exons:", round(mean(final_exon_counts), 1), "\n")
cat("  Median union exons:", median(final_exon_counts), "\n")
cat("  Range:", min(final_exon_counts), "-", max(final_exon_counts), "\n\n")

cat("Union exon distribution:\n")
exon_dist <- table(final_exon_counts)
print(head(exon_dist, 20))
cat("\n")

# ============================================================================
# Final Summary
# ============================================================================

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   UNION MODEL CONSTRUCTION COMPLETE                            ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("RESULTS:\n")
cat("  Starting genes:", n_genes, "\n")
cat("  Union models built:", length(union_models), "\n")
cat("  After filtering:", length(union_models_final), "\n")
cat("  Genes excluded:", nrow(filtered_genes), "\n\n")

cat("RUNTIME:\n")
cat("  Total time:", round(total_time, 1), "minutes\n\n")

cat("OUTPUT LOCATION:\n")
cat(" ", output_dir, "\n\n")

cat("NEXT STEP:\n")
cat("  Rscript code/detect_events_expression_based.R\n\n")
