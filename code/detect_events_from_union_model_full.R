#!/usr/bin/env Rscript
# Detect Splicing Events from Union Exon Model - Full Dataset
# Processes all 11,379 genes with union models

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   SPLICING EVENT DETECTION - FULL DATASET                    ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Configuration
output_dir <- "results/isoform_transitions/v4.0_reference_based"
BATCH_SIZE <- 1000  # Save every 1000 genes
TSS_TES_TOLERANCE <- 20  # Minimum bp difference to call Alt_TSS/Alt_TES

# ============================================================================
# Helper Functions
# ============================================================================

#' Extract isoform's exons from union model
get_isoform_exons <- function(union_exons, isoform_id) {
  exon_records <- map_dfr(union_exons, function(union_exon) {
    variants_in_iso <- union_exon$variants %>%
      filter(isoform_id == !!isoform_id)

    if (nrow(variants_in_iso) == 0) {
      return(NULL)
    }

    variant <- variants_in_iso[1, ]

    tibble(
      exon_number = union_exon$exon_number,
      exon_type = union_exon$exon_type,
      start = variant$start,
      end = variant$end,
      is_first = variant$is_first,
      is_last = variant$is_last
    )
  })

  exon_records
}

#' Compare two isoforms using union model to detect events
detect_events_from_union <- function(union_exons, isoform_A_id, isoform_B_id, gene_strand) {
  exons_A <- get_isoform_exons(union_exons, isoform_A_id)
  exons_B <- get_isoform_exons(union_exons, isoform_B_id)

  if (nrow(exons_A) == 0 || nrow(exons_B) == 0) {
    return(tibble(event_type = character(),
                  exon_number = integer(),
                  direction = character(),
                  detail = character()))
  }

  events <- list()

  # Get first and last exons from each isoform
  first_A <- exons_A %>% filter(row_number() == 1)
  first_B <- exons_B %>% filter(row_number() == 1)
  last_A <- exons_A %>% filter(row_number() == n())
  last_B <- exons_B %>% filter(row_number() == n())

  # SKIP SINGLE-EXON GENE COMPARISONS
  # If both isoforms are single-exon (first == last), there's no splicing to analyze
  is_single_exon_A <- (first_A$exon_number == last_A$exon_number)
  is_single_exon_B <- (first_B$exon_number == last_B$exon_number)

  if (is_single_exon_A && is_single_exon_B) {
    # Both isoforms are single-exon - no splicing events possible
    return(tibble(event_type = character(),
                  exon_number = integer(),
                  direction = character(),
                  detail = character()))
  }

  # EXPLICIT Alt_TSS/Alt_TES DETECTION (outside union exon iteration)
  # This is necessary because first/last exons that differ by >20bp won't be
  # in the same union exon group, so we need to compare them explicitly

  # Check Alt_TSS (compare first exons)
  if (gene_strand == "+") {
    # Plus strand: TSS = start coordinate
    tss_diff <- abs(first_A$start - first_B$start)
    if (tss_diff > TSS_TES_TOLERANCE) {
      events[[length(events) + 1]] <- tibble(
        event_type = "Alt_TSS",
        exon_number = NA_integer_,  # First exons may be at different union positions
        direction = "none",  # No gain/loss for Alt_TSS - just different positions
        detail = sprintf("TSS diff: %dbp (A: %d, B: %d)",
                        tss_diff, first_A$start, first_B$start)
      )
    }
  } else if (gene_strand == "-") {
    # Minus strand: TSS = end coordinate
    tss_diff <- abs(first_A$end - first_B$end)
    if (tss_diff > TSS_TES_TOLERANCE) {
      events[[length(events) + 1]] <- tibble(
        event_type = "Alt_TSS",
        exon_number = NA_integer_,  # First exons may be at different union positions
        direction = "none",  # No gain/loss for Alt_TSS - just different positions
        detail = sprintf("TSS diff: %dbp (A: %d, B: %d)",
                        tss_diff, first_A$end, first_B$end)
      )
    }
  }

  # Check Alt_TES (compare last exons)
  if (gene_strand == "+") {
    # Plus strand: TES = end coordinate
    tes_diff <- abs(last_A$end - last_B$end)
    if (tes_diff > TSS_TES_TOLERANCE) {
      events[[length(events) + 1]] <- tibble(
        event_type = "Alt_TES",
        exon_number = NA_integer_,  # Last exons may be at different union positions
        direction = "none",  # No gain/loss for Alt_TES - just different positions
        detail = sprintf("TES diff: %dbp (A: %d, B: %d)",
                        tes_diff, last_A$end, last_B$end)
      )
    }
  } else if (gene_strand == "-") {
    # Minus strand: TES = start coordinate
    tes_diff <- abs(last_A$start - last_B$start)
    if (tes_diff > TSS_TES_TOLERANCE) {
      events[[length(events) + 1]] <- tibble(
        event_type = "Alt_TES",
        exon_number = NA_integer_,  # Last exons may be at different union positions
        direction = "none",  # No gain/loss for Alt_TES - just different positions
        detail = sprintf("TES diff: %dbp (A: %d, B: %d)",
                        tes_diff, last_A$start, last_B$start)
      )
    }
  }

  for (union_exon in union_exons) {
    exon_num <- union_exon$exon_number
    exon_type <- union_exon$exon_type

    variants_A <- union_exon$variants %>% filter(isoform_id == isoform_A_id)
    variants_B <- union_exon$variants %>% filter(isoform_id == isoform_B_id)

    has_A <- nrow(variants_A) > 0
    has_B <- nrow(variants_B) > 0

    # Case 1: Both isoforms have this exon
    if (has_A && has_B) {
      var_A <- variants_A[1, ]
      var_B <- variants_B[1, ]

      same_start <- var_A$start == var_B$start
      same_end <- var_A$end == var_B$end

      if (same_start && same_end) {
        # Constitutive exon
        events[[length(events) + 1]] <- tibble(
          event_type = "CONST",
          exon_number = exon_num,
          direction = "none",
          detail = sprintf("Exon %d identical", exon_num)
        )

      } else {
        # NOTE: Alt_TSS and Alt_TES are now handled explicitly above,
        # outside the union exon iteration, to avoid the grouping issue

        # CRITICAL: Skip terminal boundary differences for first/last exons
        # to prevent double-counting with explicit Alt_TSS/Alt_TES detection
        is_first_exon <- var_A$is_first || var_B$is_first
        is_last_exon <- var_A$is_last || var_B$is_last

        # For first exons: only process internal boundary (not TSS)
        # For last exons: only process internal boundary (not TES)
        # For single-exon genes (both first AND last): skip entirely
        if (is_first_exon && is_last_exon) {
          # Single-exon gene - all boundaries handled by Alt_TSS/Alt_TES
          next
        }

        # Adjust which boundaries to process for terminal exons
        # We skip TSS and TES boundaries (handled by explicit Alt_TSS/Alt_TES detection)
        # But we DO process internal splice sites (donors and acceptors)
        if (gene_strand == "+") {
          # Plus strand: TSS at start, TES at end
          process_start <- !is_first_exon  # Skip TSS (first exon start)
          process_end <- !is_last_exon     # Skip TES (last exon end)
        } else {
          # Minus strand: coordinates reversed - TES at start, TSS at end
          process_start <- !is_last_exon   # Skip TES (last exon start on minus)
          process_end <- !is_first_exon    # Skip TSS (first exon end on minus)
        }

        # Apply boundary processing restrictions
        check_start <- !same_start && process_start
        check_end <- !same_end && process_end

        # Internal exon with different coordinates - STRAND-AWARE
        # Key: On minus strand, event types are SWAPPED
        # Plus: different start = A3SS, different end = A5SS
        # Minus: different start = A5SS, different end = A3SS

        if (check_start && check_end) {
          # Both boundaries differ - always A5SS + A3SS
          if (gene_strand == "+") {
            direction_5 <- if_else(var_A$end > var_B$end, "gain", "loss")
            direction_3 <- if_else(var_A$start < var_B$start, "gain", "loss")
          } else {
            direction_5 <- if_else(var_A$start < var_B$start, "gain", "loss")
            direction_3 <- if_else(var_A$end > var_B$end, "gain", "loss")
          }

          events[[length(events) + 1]] <- tibble(
            event_type = "A5SS",
            exon_number = exon_num,
            direction = direction_5,
            detail = sprintf("5'SS diff")
          )

          events[[length(events) + 1]] <- tibble(
            event_type = "A3SS",
            exon_number = exon_num,
            direction = direction_3,
            detail = sprintf("3'SS diff")
          )

        } else if (check_start) {
          # Different start only
          if (gene_strand == "+") {
            # Plus: start = 3'SS (acceptor)
            direction <- if_else(var_A$start < var_B$start, "gain", "loss")
            event_name <- "A3SS"
          } else {
            # Minus: start = 5'SS (donor) - SWAPPED!
            direction <- if_else(var_A$start < var_B$start, "gain", "loss")
            event_name <- "A5SS"
          }

          events[[length(events) + 1]] <- tibble(
            event_type = event_name,
            exon_number = exon_num,
            direction = direction,
            detail = sprintf("%s: %d vs %d", event_name, var_A$start, var_B$start)
          )

        } else if (check_end) {
          # Different end only
          if (gene_strand == "+") {
            # Plus: end = 5'SS (donor)
            direction <- if_else(var_A$end > var_B$end, "gain", "loss")
            event_name <- "A5SS"
          } else {
            # Minus: end = 3'SS (acceptor) - SWAPPED!
            direction <- if_else(var_A$end > var_B$end, "gain", "loss")
            event_name <- "A3SS"
          }

          events[[length(events) + 1]] <- tibble(
            event_type = event_name,
            exon_number = exon_num,
            direction = direction,
            detail = sprintf("%s: %d vs %d", event_name, var_A$end, var_B$end)
          )
        }
      }

    } else if (has_A && !has_B) {
      # Case 2: Only A has this exon
      # Check if it's a first/last exon (already handled as Alt_TSS/TES)
      var_A <- variants_A[1, ]
      if (!var_A$is_first && !var_A$is_last) {
        # Only call SE for internal exons
        events[[length(events) + 1]] <- tibble(
          event_type = "SE",
          exon_number = exon_num,
          direction = "gain",
          detail = sprintf("Exon %d gained", exon_num)
        )
      }

    } else if (!has_A && has_B) {
      # Case 3: Only B has this exon
      # Check if it's a first/last exon (already handled as Alt_TSS/TES)
      var_B <- variants_B[1, ]
      if (!var_B$is_first && !var_B$is_last) {
        # Only call SE for internal exons
        events[[length(events) + 1]] <- tibble(
          event_type = "SE",
          exon_number = exon_num,
          direction = "loss",
          detail = sprintf("Exon %d lost", exon_num)
        )
      }
    }
  }

  if (length(events) == 0) {
    return(tibble(event_type = character(),
                  exon_number = integer(),
                  direction = character(),
                  detail = character()))
  }

  bind_rows(events)
}

cat("Helper functions loaded\n\n")

# ============================================================================
# Load Data
# ============================================================================

cat("Loading data...\n")

# Load union models
union_models <- readRDS(file.path(output_dir, "union_exon_models_full.rds"))
cat("  Union models:", length(union_models), "genes\n")

# Load isoform list
isoforms_to_process <- readRDS(file.path(output_dir, "isoforms_for_union_model.rds"))

# Load exon structures for strand information
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))
strand_lookup <- exon_structures %>%
  select(isoform_id, strand) %>%
  distinct()
cat("  Loaded strand information for", nrow(strand_lookup), "isoforms\n")

# Filter to genes with union models
genes_with_models <- names(union_models)
isoforms_filtered <- isoforms_to_process %>%
  filter(gene_id %in% genes_with_models)

cat("  Isoforms in genes with models:", nrow(isoforms_filtered), "\n\n")

# Filter genes by union exon count (>1 and <=20)
cat("Filtering genes by union exon count...\n")
genes_before_filter <- length(union_models)
genes_to_keep <- names(union_models)[sapply(union_models, function(m) {
  m$n_union_exons > 1 && m$n_union_exons <= 20
})]

union_models <- union_models[genes_to_keep]
genes_with_models <- names(union_models)
isoforms_filtered <- isoforms_filtered %>%
  filter(gene_id %in% genes_with_models)

cat(sprintf("  Before filter: %d genes\n", genes_before_filter))
cat(sprintf("  After filter (>1 and ≤20 union exons): %d genes (%.1f%% retained)\n",
            length(union_models),
            100 * length(union_models) / genes_before_filter))
cat(sprintf("  Genes removed by exon count: %d (%.1f%%)\n\n",
            genes_before_filter - length(union_models),
            100 * (genes_before_filter - length(union_models)) / genes_before_filter))

# Filter out fusion genes (ENSG_ENSG format)
cat("Filtering fusion genes...\n")
genes_after_exon_filter <- length(union_models)
fusion_pattern <- "ENSG[0-9]+_ENSG[0-9]+"
is_fusion <- grepl(fusion_pattern, names(union_models))
fusion_genes <- names(union_models)[is_fusion]

# Save fusion gene info BEFORE removing them
if (length(fusion_genes) > 0) {
  fusion_exclusion <- tibble(
    gene_id = fusion_genes,
    reason = "fusion_gene",
    n_isoforms = sapply(fusion_genes, function(g) union_models[[g]]$n_isoforms),
    n_union_exons = sapply(fusion_genes, function(g) union_models[[g]]$n_union_exons)
  )
  saveRDS(fusion_exclusion, file.path(output_dir, "excluded_fusion_genes.rds"))
  write_tsv(fusion_exclusion, file.path(output_dir, "excluded_fusion_genes.tsv"))
}

# Count isoforms before filtering
fusion_isoforms <- sum(sapply(fusion_genes, function(g) union_models[[g]]$n_isoforms))

# Remove fusion genes
union_models <- union_models[!is_fusion]
genes_with_models <- names(union_models)
isoforms_filtered <- isoforms_filtered %>%
  filter(gene_id %in% genes_with_models)

cat(sprintf("  Fusion genes identified: %d (%.1f%%)\n",
            length(fusion_genes),
            100 * length(fusion_genes) / genes_after_exon_filter))
cat(sprintf("  Fusion isoforms: %d\n", fusion_isoforms))
cat(sprintf("  After fusion filter: %d genes (%.1f%% retained)\n",
            length(union_models),
            100 * length(union_models) / genes_after_exon_filter))
if (length(fusion_genes) > 0) {
  cat("  Saved exclusion list: excluded_fusion_genes.rds/.tsv\n")
}
cat("\n")

# ============================================================================
# Process All Genes
# ============================================================================

cat("Processing genes...\n")
cat("  Progress will be reported every 100 genes\n")
cat("  Saving checkpoints every", BATCH_SIZE, "genes\n\n")

# Initialize results (no resume for this run)
results_file <- file.path(output_dir, "event_vectors_full.rds")

# Remove old results if exists
if (file.exists(results_file)) {
  cat("Removing old results file...\n")
  file.remove(results_file)
}

all_transitions <- list()
start_idx <- 1

start_time <- Sys.time()
genes_in_batch <- 0

for (gene_idx in start_idx:length(genes_with_models)) {
  gene_id <- genes_with_models[gene_idx]
  union_model <- union_models[[gene_id]]

  # Progress reporting
  if (gene_idx %% 100 == 0 || gene_idx == start_idx) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    genes_per_min <- (gene_idx - start_idx + 1) / max(elapsed, 0.01)
    remaining_genes <- length(genes_with_models) - gene_idx
    eta_mins <- remaining_genes / genes_per_min

    cat(sprintf("[%d/%d] %s | Elapsed: %.1f min | ETA: %.1f min\n",
                gene_idx, length(genes_with_models), gene_id,
                elapsed, eta_mins))
  }

  # Get isoforms for this gene
  gene_isoforms <- isoforms_filtered %>%
    filter(gene_id == !!gene_id) %>%
    pull(isoform_id)

  if (length(gene_isoforms) < 2) next

  # Get strand from first isoform (all isoforms in a gene have same strand)
  gene_strand <- strand_lookup %>%
    filter(isoform_id == gene_isoforms[1]) %>%
    pull(strand) %>%
    .[1]

  if (is.na(gene_strand)) {
    warning(sprintf("No strand info for gene %s, skipping", gene_id))
    next
  }

  # Generate all pairwise transitions
  transitions_for_gene <- expand_grid(
    isoform_A = gene_isoforms,
    isoform_B = gene_isoforms
  ) %>%
    filter(isoform_A != isoform_B)

  # Detect events for each transition
  transition_results <- map_dfr(seq_len(nrow(transitions_for_gene)), function(i) {
    iso_A <- transitions_for_gene$isoform_A[i]
    iso_B <- transitions_for_gene$isoform_B[i]

    events <- tryCatch({
      detect_events_from_union(union_model$union_exons, iso_A, iso_B, gene_strand)
    }, error = function(e) {
      tibble(event_type = character(),
             exon_number = integer(),
             direction = character(),
             detail = character())
    })

    # Count events by type
    event_counts <- events %>%
      count(event_type, name = "count") %>%
      pivot_wider(names_from = event_type,
                  values_from = count,
                  values_fill = 0,
                  names_prefix = "n_")

    # Ensure all event type columns exist
    all_event_types <- c("n_Alt_TSS", "n_Alt_TES", "n_SE", "n_A5SS", "n_A3SS", "n_CONST")
    for (col in all_event_types) {
      if (!col %in% names(event_counts)) {
        event_counts[[col]] <- 0
      }
    }

    tibble(
      gene_id = gene_id,
      gene_strand = gene_strand,  # CRITICAL: Save strand info for validation
      isoform_A = iso_A,
      isoform_B = iso_B,
      n_union_exons = union_model$n_union_exons,
      event_vector = list(events),
      !!!event_counts
    )
  })

  all_transitions[[gene_id]] <- transition_results
  genes_in_batch <- genes_in_batch + 1

  # Save checkpoint
  if (genes_in_batch >= BATCH_SIZE) {
    cat(sprintf("  [CHECKPOINT] Saving at gene %d...\n", gene_idx))
    all_results <- bind_rows(all_transitions)
    saveRDS(all_results, results_file)
    genes_in_batch <- 0
  }
}

# Final save
cat("\nSaving final results...\n")
all_results <- bind_rows(all_transitions)
saveRDS(all_results, results_file)
cat("  Saved: event_vectors_full.rds\n")

# Also save as TSV (without event_vector column if it exists)
if ("event_vector" %in% names(all_results)) {
  all_results_tsv <- all_results %>% select(-event_vector)
} else {
  all_results_tsv <- all_results
}
write_tsv(all_results_tsv, file.path(output_dir, "event_vectors_full.tsv"))
cat("  Saved: event_vectors_full.tsv\n")

# ============================================================================
# Summary Statistics
# ============================================================================

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EVENT DETECTION COMPLETE                                    ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("  Total runtime: %.1f minutes\n", total_time))
cat(sprintf("  Genes processed: %d\n", length(unique(all_results$gene_id))))
cat(sprintf("  Transitions analyzed: %d\n", nrow(all_results)))
cat(sprintf("  Mean transitions per gene: %.1f\n\n", nrow(all_results) / length(unique(all_results$gene_id))))

# Event type distribution
cat("Event type distribution:\n")
event_summary <- all_results %>%
  summarise(
    total_transitions = n(),
    transitions_with_Alt_TSS = sum(n_Alt_TSS > 0),
    transitions_with_Alt_TES = sum(n_Alt_TES > 0),
    transitions_with_SE = sum(n_SE > 0),
    transitions_with_A5SS = sum(n_A5SS > 0),
    transitions_with_A3SS = sum(n_A3SS > 0),
    transitions_with_CONST = sum(n_CONST > 0),
    total_Alt_TSS = sum(n_Alt_TSS),
    total_Alt_TES = sum(n_Alt_TES),
    total_SE = sum(n_SE),
    total_A5SS = sum(n_A5SS),
    total_A3SS = sum(n_A3SS),
    total_CONST = sum(n_CONST)
  )

cat(sprintf("  Alt_TSS: %d total, in %.1f%% of transitions\n",
            event_summary$total_Alt_TSS,
            100 * event_summary$transitions_with_Alt_TSS / event_summary$total_transitions))

cat(sprintf("  Alt_TES: %d total, in %.1f%% of transitions\n",
            event_summary$total_Alt_TES,
            100 * event_summary$transitions_with_Alt_TES / event_summary$total_transitions))

cat(sprintf("  SE: %d total, in %.1f%% of transitions\n",
            event_summary$total_SE,
            100 * event_summary$transitions_with_SE / event_summary$total_transitions))

cat(sprintf("  A5SS: %d total, in %.1f%% of transitions\n",
            event_summary$total_A5SS,
            100 * event_summary$transitions_with_A5SS / event_summary$total_transitions))

cat(sprintf("  A3SS: %d total, in %.1f%% of transitions\n",
            event_summary$total_A3SS,
            100 * event_summary$transitions_with_A3SS / event_summary$total_transitions))

cat(sprintf("  CONST: %d total, in %.1f%% of transitions\n\n",
            event_summary$total_CONST,
            100 * event_summary$transitions_with_CONST / event_summary$total_transitions))

# Complexity distribution
cat("Event complexity per transition:\n")
all_results <- all_results %>%
  mutate(n_total_events = n_Alt_TSS + n_Alt_TES + n_SE + n_A5SS + n_A3SS)

complexity_dist <- all_results %>%
  mutate(complexity_bin = case_when(
    n_total_events == 0 ~ "0 (identical)",
    n_total_events == 1 ~ "1",
    n_total_events == 2 ~ "2",
    n_total_events %in% 3:5 ~ "3-5",
    n_total_events %in% 6:10 ~ "6-10",
    n_total_events > 10 ~ "10+"
  )) %>%
  count(complexity_bin)

for (i in seq_len(nrow(complexity_dist))) {
  cat(sprintf("  %s events: %d transitions (%.1f%%)\n",
              complexity_dist$complexity_bin[i],
              complexity_dist$n[i],
              100 * complexity_dist$n[i] / nrow(all_results)))
}

cat("\n═══ NEXT STEPS ═══\n\n")
cat("  1. Run MXE detection: Rscript code/detect_mxe_from_union_model.R\n")
cat("  2. Perform co-occurrence analysis\n")
cat("  3. Analyze topological distances\n\n")
