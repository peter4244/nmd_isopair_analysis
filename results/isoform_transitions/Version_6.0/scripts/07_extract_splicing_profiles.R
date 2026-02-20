#!/usr/bin/env Rscript

################################################################################
# Script 07: Extract Splicing Choice Profiles
################################################################################
#
# Purpose:
#   Build splicing choice profiles by comparing non-dominant to dominant isoforms.
#   Each profile characterizes the structural differences between isoforms.
#   Performs comprehensive event detection on FILTERED data from Script 06.
#
# Input:
#   - data/isoform_union_exons_annotated_filtered.rds
#   - data/dominant_isoforms_filtered.rds
#   - data/union_exons_filtered.rds
#   - data/isoform_structures_filtered.rds
#
# Output:
#   - data/splicing_choice_profiles.rds
#   - data/splicing_choice_profiles_intermediate.rds (early output for development)
#
################################################################################

library(tidyverse)

# ==============================================================================
# DETECTION THRESHOLDS (adjust these for sensitivity analyses)
# ==============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 7: Extract Splicing Choice Profiles                   ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ==============================================================================
# 0. Load Event Detection Functions (Shared Library)
# ==============================================================================

cat("Loading event detection functions...\n")
source("scripts/event_detection_functions.R")
cat("  ✓ Functions loaded\n\n")

cat("Detection thresholds:\n")
cat(sprintf("  TSS tolerance: %d bp\n", TSS_TOLERANCE))
cat(sprintf("  TES tolerance: %d bp\n", TES_TOLERANCE))
cat(sprintf("  Splice site threshold: %d bp\n\n", SPLICE_SITE_THRESHOLD))

# ==============================================================================
# 0b. Legacy - Function definitions now in event_detection_functions.R
# ==============================================================================
# All detection functions are now sourced from the shared library above

# ==============================================================================
# 1. Load Data
# ==============================================================================

cat("Loading filtered data...\n")
annotated_mapping <- readRDS("data/isoform_union_exons_annotated_filtered.rds")
dominant_isoforms <- readRDS("data/dominant_isoforms_filtered.rds")
union_exons <- readRDS("data/union_exons_filtered.rds")
isoform_structures_compact <- readRDS("data/isoform_structures_filtered.rds")

cat(sprintf("  Annotated mappings: %d\n", nrow(annotated_mapping)))
cat(sprintf("  Dominant isoforms: %d\n", nrow(dominant_isoforms)))
cat(sprintf("  Union exons: %d\n", nrow(union_exons)))
cat(sprintf("  Isoform structures (compact): %d\n", nrow(isoform_structures_compact)))

# Expand isoform_structures to one row per exon (needed for event detection)
cat("\nExpanding isoform structures to one row per exon...\n")
isoform_structures <- isoform_structures_compact %>%
  rowwise() %>%
  mutate(
    exon_data = list(tibble(
      exon_number = seq_along(exon_starts),
      exon_start = exon_starts,
      exon_end = exon_ends
    ))
  ) %>%
  unnest(exon_data) %>%
  select(isoform_id, gene_id, seqnames, strand, exon_number, exon_start, exon_end) %>%
  ungroup()

cat(sprintf("  Expanded to %d exon rows\n", nrow(isoform_structures)))

# ==============================================================================
# 2. Identify Non-Dominant Isoforms
# ==============================================================================

cat("\nIdentifying non-dominant isoforms...\n")

# Get list of all isoforms per gene
all_isoforms <- annotated_mapping %>%
  distinct(gene_id, isoform_id)

# Mark which are dominant
isoforms_classified <- all_isoforms %>%
  left_join(
    dominant_isoforms %>% select(gene_id, dominant_isoform_id),
    by = "gene_id"
  ) %>%
  mutate(
    is_dominant = (isoform_id == dominant_isoform_id),
    role = if_else(is_dominant, "dominant", "non_dominant")
  )

cat(sprintf("  Dominant isoforms: %d\n", sum(isoforms_classified$is_dominant)))
cat(sprintf("  Non-dominant isoforms: %d\n", sum(!isoforms_classified$is_dominant)))

# ==============================================================================
# 3. Build Comparison Profiles with Event Detection
# ==============================================================================

cat("\nBuilding splicing choice profiles with event detection...\n")
cat("  Comparing each non-dominant to dominant isoform...\n")

# For each gene, compare non-dominant to dominant
splicing_profiles <- list()

# Only process genes that have both dominant info AND structure data
genes_with_structures <- unique(isoform_structures$gene_id)
genes_with_dominant <- intersect(
  unique(dominant_isoforms$gene_id),
  genes_with_structures
)

cat(sprintf("  Genes with both dominant and structure data: %d\n", length(genes_with_dominant)))

# Process in batches with progress reporting
batch_size <- 1000
n_batches <- ceiling(length(genes_with_dominant) / batch_size)

for (batch_idx in 1:n_batches) {
  start_idx <- (batch_idx - 1) * batch_size + 1
  end_idx <- min(batch_idx * batch_size, length(genes_with_dominant))
  batch_genes <- genes_with_dominant[start_idx:end_idx]

  for (gene in batch_genes) {
    # Get dominant isoform
    dom_iso <- dominant_isoforms %>%
      filter(gene_id == gene) %>%
      pull(dominant_isoform_id)

    if (length(dom_iso) != 1) next

    # Get gene strand
    gene_strand <- union_exons %>%
      filter(gene_id == gene) %>%
      pull(strand) %>%
      unique()

    if (length(gene_strand) != 1) next

    # Get non-dominant isoforms for this gene
    non_dom_isos <- isoforms_classified %>%
      filter(gene_id == gene, !is_dominant) %>%
      pull(isoform_id)

    if (length(non_dom_isos) == 0) next

    # Get structure for dominant isoform
    dom_structure <- isoform_structures %>%
      filter(isoform_id == dom_iso) %>%
      arrange(exon_number)

    if (nrow(dom_structure) == 0) next

    # Get union exon usage for dominant
    dom_exons <- annotated_mapping %>%
      filter(isoform_id == dom_iso) %>%
      select(union_exon_id, region_type_dom = region_type)

    # Compare each non-dominant to dominant
    for (non_dom_iso in non_dom_isos) {

      # Get structure for non-dominant isoform
      non_dom_structure <- isoform_structures %>%
        filter(isoform_id == non_dom_iso) %>%
        arrange(exon_number)

      if (nrow(non_dom_structure) == 0) next

      # Get union exon usage for non-dominant
      non_dom_exons <- annotated_mapping %>%
        filter(isoform_id == non_dom_iso) %>%
        select(union_exon_id, region_type_non_dom = region_type)

      # Full outer join to find differences
      comparison <- union_exons %>%
        filter(gene_id == gene) %>%
        select(gene_id, union_exon_id, union_exon_number, union_exon_start, union_exon_end) %>%
        left_join(dom_exons, by = "union_exon_id") %>%
        left_join(non_dom_exons, by = "union_exon_id") %>%
        mutate(
          in_dominant = !is.na(region_type_dom),
          in_non_dominant = !is.na(region_type_non_dom),
          exon_status = case_when(
            in_dominant & in_non_dominant ~ "shared",
            in_dominant & !in_non_dominant ~ "dominant_only",
            !in_dominant & in_non_dominant ~ "non_dominant_only",
            TRUE ~ "neither"
          )
        )

      # ===== EVENT DETECTION (validated hierarchical algorithm) =====
      # Uses detect_events_for_pair() which performs hierarchical detection:
      # IR → Boundary (A5SS/A3SS/Partial_IR) → SE → Missing_Internal → Terminal

      # Prepare exon structures for detect_events_for_pair()
      dom_exons_for_detection <- dom_structure %>%
        rename(chr = seqnames, transcript_id = isoform_id)

      comp_exons_for_detection <- non_dom_structure %>%
        rename(chr = seqnames, transcript_id = isoform_id)

      # Run hierarchical event detection
      events <- detect_events_for_pair(
        dom_exons_for_detection, comp_exons_for_detection,
        gene, dom_iso, non_dom_iso, gene_strand
      )

      # Aggregate event counts (backward compatible with Scripts 08-12)
      if (nrow(events) > 0) {
        n_a5ss <- sum(events$event_type == "A5SS")
        n_a3ss <- sum(events$event_type == "A3SS")
        n_partial_ir <- sum(events$event_type %in% c("Partial_IR_5", "Partial_IR_3"))
        n_ir <- sum(events$event_type == "IR")
        n_se <- sum(events$event_type == "SE")
        n_missing_internal <- sum(events$event_type == "Missing_Internal")
        n_ir_diff <- sum(grepl("^IR_diff", events$event_type))
        n_alt_tss <- sum(events$event_type == "Alt_TSS")
        n_alt_tes <- sum(events$event_type == "Alt_TES")
        tss_changed <- n_alt_tss > 0
        tes_changed <- n_alt_tes > 0
      } else {
        n_a5ss <- 0L; n_a3ss <- 0L; n_partial_ir <- 0L
        n_ir <- 0L; n_se <- 0L; n_missing_internal <- 0L
        n_ir_diff <- 0L; n_alt_tss <- 0L; n_alt_tes <- 0L
        tss_changed <- FALSE; tes_changed <- FALSE
      }
      n_dual_boundary <- 0L  # No longer used; decomposed into two events

      # Get complexity metrics from compact isoform_structures
      dom_complexity <- isoform_structures_compact %>%
        filter(isoform_id == dom_iso) %>%
        select(n_exons, n_junctions, tx_start, tx_end) %>%
        mutate(length = tx_end - tx_start + 1)

      non_dom_complexity <- isoform_structures_compact %>%
        filter(isoform_id == non_dom_iso) %>%
        select(n_exons, n_junctions, tx_start, tx_end) %>%
        mutate(length = tx_end - tx_start + 1)

      # Calculate summary statistics for this comparison
      profile <- tibble(
        gene_id = gene,
        dominant_isoform_id = dom_iso,
        non_dominant_isoform_id = non_dom_iso,

        # Exon counts
        n_union_exons_total = nrow(comparison),
        n_exons_shared = sum(comparison$exon_status == "shared"),
        n_exons_dominant_only = sum(comparison$exon_status == "dominant_only"),
        n_exons_non_dominant_only = sum(comparison$exon_status == "non_dominant_only"),

        # Structural info
        n_exons_in_dominant = sum(comparison$in_dominant),
        n_exons_in_non_dominant = sum(comparison$in_non_dominant),

        # Complexity metrics for dominant isoform
        n_exons_dom = dom_complexity$n_exons,
        n_junctions_dom = dom_complexity$n_junctions,
        length_dom = dom_complexity$length,

        # Complexity metrics for non-dominant isoform (used by Script 08)
        n_exons_non_dom = non_dom_complexity$n_exons,
        n_junctions_non_dom = non_dom_complexity$n_junctions,
        length_non_dom = non_dom_complexity$length,

        # TSS/TES changes
        tss_changed = tss_changed,
        tes_changed = tes_changed,

        # Splicing events (event-level counts)
        n_a5ss = n_a5ss,
        n_a3ss = n_a3ss,
        n_partial_ir = n_partial_ir,
        n_ir = n_ir,
        n_se = n_se,
        n_missing_internal = n_missing_internal,
        n_ir_diff = n_ir_diff,
        n_alt_tss = n_alt_tss,
        n_alt_tes = n_alt_tes,

        # Detailed events from hierarchical detection
        n_events = nrow(events),
        detailed_events = list(events),

        # Total differences (needed for Script 08)
        n_differences = n_exons_dominant_only + n_exons_non_dominant_only,

        # Store full comparison as nested tibble
        comparison_detail = list(comparison)
      )

      splicing_profiles[[paste0(gene, "_", non_dom_iso)]] <- profile
    }
  }

  # Progress report
  cat(sprintf("  Processed batch %d/%d (%d genes, %.1f%% complete)\n",
              batch_idx, n_batches, length(batch_genes),
              100 * end_idx / length(genes_with_dominant)))

  # Save intermediate results after batch 1 for downstream development
  if (batch_idx == 1 && length(splicing_profiles) > 0) {
    cat("\n  ═══════════════════════════════════════════════════════════\n")
    cat("  CHECKPOINT: Saving intermediate results for development\n")
    cat("  ═══════════════════════════════════════════════════════════\n")
    intermediate_profiles <- bind_rows(splicing_profiles)
    saveRDS(intermediate_profiles, "data/splicing_choice_profiles_intermediate.rds")
    cat(sprintf("  ✓ Saved %d profiles to:\n", nrow(intermediate_profiles)))
    cat("     data/splicing_choice_profiles_intermediate.rds\n")
    cat(sprintf("  ✓ Covering %d genes\n", n_distinct(intermediate_profiles$gene_id)))
    cat("\n  This file can be used for developing Scripts 08-12\n")
    cat("  while the full processing continues...\n")
    cat("  ═══════════════════════════════════════════════════════════\n\n")
  }
}

# Combine profiles
splicing_profiles_df <- bind_rows(splicing_profiles)

cat(sprintf("\n  Created %d splicing choice profiles\n", nrow(splicing_profiles_df)))
cat(sprintf("  Covering %d genes\n", n_distinct(splicing_profiles_df$gene_id)))

# ==============================================================================
# 4. Summary Statistics
# ==============================================================================

cat("\nProfile statistics:\n")

cat(sprintf("  Mean union exons per gene: %.1f\n",
            mean(splicing_profiles_df$n_union_exons_total)))
cat(sprintf("  Mean shared exons: %.1f (%.1f%%)\n",
            mean(splicing_profiles_df$n_exons_shared),
            100 * mean(splicing_profiles_df$n_exons_shared / splicing_profiles_df$n_union_exons_total)))
cat(sprintf("  Mean dominant-only exons: %.1f (%.1f%%)\n",
            mean(splicing_profiles_df$n_exons_dominant_only),
            100 * mean(splicing_profiles_df$n_exons_dominant_only / splicing_profiles_df$n_union_exons_total)))
cat(sprintf("  Mean non-dominant-only exons: %.1f (%.1f%%)\n",
            mean(splicing_profiles_df$n_exons_non_dominant_only),
            100 * mean(splicing_profiles_df$n_exons_non_dominant_only / splicing_profiles_df$n_union_exons_total)))

cat("\nComplexity metrics:\n")
cat(sprintf("  Dominant isoforms - mean exons: %.1f, junctions: %.1f, length: %.0f bp\n",
            mean(splicing_profiles_df$n_exons_dom, na.rm = TRUE),
            mean(splicing_profiles_df$n_junctions_dom, na.rm = TRUE),
            mean(splicing_profiles_df$length_dom, na.rm = TRUE)))
cat(sprintf("  Non-dominant isoforms - mean exons: %.1f, junctions: %.1f, length: %.0f bp\n",
            mean(splicing_profiles_df$n_exons_non_dom, na.rm = TRUE),
            mean(splicing_profiles_df$n_junctions_non_dom, na.rm = TRUE),
            mean(splicing_profiles_df$length_non_dom, na.rm = TRUE)))

cat("\nTSS/TES changes:\n")
cat(sprintf("  Profiles with TSS change: %d (%.1f%%)\n",
            sum(splicing_profiles_df$tss_changed, na.rm = TRUE),
            100 * mean(splicing_profiles_df$tss_changed, na.rm = TRUE)))
cat(sprintf("  Profiles with TES change: %d (%.1f%%)\n",
            sum(splicing_profiles_df$tes_changed, na.rm = TRUE),
            100 * mean(splicing_profiles_df$tes_changed, na.rm = TRUE)))

cat("\nSplicing events detected:\n")
cat(sprintf("  Total A5SS events: %d (mean %.2f per profile)\n",
            sum(splicing_profiles_df$n_a5ss), mean(splicing_profiles_df$n_a5ss)))
cat(sprintf("  Total A3SS events: %d (mean %.2f per profile)\n",
            sum(splicing_profiles_df$n_a3ss), mean(splicing_profiles_df$n_a3ss)))
cat(sprintf("  Total Partial_IR events: %d (mean %.2f per profile)\n",
            sum(splicing_profiles_df$n_partial_ir), mean(splicing_profiles_df$n_partial_ir)))
cat(sprintf("  Total IR events: %d (mean %.2f per profile)\n",
            sum(splicing_profiles_df$n_ir), mean(splicing_profiles_df$n_ir)))
cat(sprintf("  Total SE events: %d (mean %.2f per profile)\n",
            sum(splicing_profiles_df$n_se), mean(splicing_profiles_df$n_se)))
cat(sprintf("  Total Missing_Internal events: %d (mean %.2f per profile)\n",
            sum(splicing_profiles_df$n_missing_internal), mean(splicing_profiles_df$n_missing_internal)))

# Distribution of profile types
cat("\nProfile complexity:\n")
complexity_bins <- splicing_profiles_df %>%
  mutate(
    complexity = case_when(
      n_differences == 0 ~ "identical",
      n_differences <= 2 ~ "simple (1-2 diffs)",
      n_differences <= 5 ~ "moderate (3-5 diffs)",
      TRUE ~ "complex (6+ diffs)"
    )
  ) %>%
  count(complexity) %>%
  mutate(percentage = 100 * n / sum(n))

for (i in 1:nrow(complexity_bins)) {
  cat(sprintf("  %s: %d (%.1f%%)\n",
              complexity_bins$complexity[i],
              complexity_bins$n[i],
              complexity_bins$percentage[i]))
}

# ==============================================================================
# 5. Save Output
# ==============================================================================

cat("\nSaving splicing choice profiles...\n")
saveRDS(splicing_profiles_df, "data/splicing_choice_profiles.rds")
cat("  ✓ data/splicing_choice_profiles.rds\n")

cat("\n✓ Step 7 complete\n")
cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Splicing profiles created: %d\n", nrow(splicing_profiles_df)))
cat(sprintf("Genes covered: %d\n", n_distinct(splicing_profiles_df$gene_id)))
cat(sprintf("Mean differences per profile: %.1f exons\n",
            mean(splicing_profiles_df$n_differences)))
cat("\n")
cat("Event-level counts (variable, events can co-occur):\n")
cat(sprintf("  A5SS: %d\n", sum(splicing_profiles_df$n_a5ss)))
cat(sprintf("  A3SS: %d\n", sum(splicing_profiles_df$n_a3ss)))
cat(sprintf("  Partial_IR: %d\n", sum(splicing_profiles_df$n_partial_ir)))
cat(sprintf("  IR: %d\n", sum(splicing_profiles_df$n_ir)))
cat(sprintf("  SE: %d\n", sum(splicing_profiles_df$n_se)))
cat(sprintf("  Missing_Internal: %d\n", sum(splicing_profiles_df$n_missing_internal)))
cat(sprintf("  IR_diff: %d\n", sum(splicing_profiles_df$n_ir_diff)))
cat(sprintf("  Alt_TSS: %d\n", sum(splicing_profiles_df$n_alt_tss)))
cat(sprintf("  Alt_TES: %d\n", sum(splicing_profiles_df$n_alt_tes)))
cat("\n")
cat(sprintf("Total events detected: %d\n", sum(splicing_profiles_df$n_events)))
cat(sprintf("  Mean events per profile: %.1f\n",
            mean(splicing_profiles_df$n_events)))
cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")
