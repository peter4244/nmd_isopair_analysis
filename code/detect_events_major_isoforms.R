#!/usr/bin/env Rscript
# Event Detection for Major Isoforms - v4.0
#
# Purpose:
#   Detect splicing events by comparing isoform pairs using union exon models
#
# Approach:
#   For each gene, compare all pairs of major isoforms
#   Identify which union exons each isoform contains
#   Classify differences as specific event types
#
# Inputs:
#   - major_isoforms/union_exon_models_major.rds (3,797 genes)
#   - major_isoforms/major_isoforms_dmso.rds (expression data)
#
# Outputs:
#   - major_isoforms/isoform_pairs_events.rds (all pairwise comparisons with events)
#   - major_isoforms/event_summary.tsv (summary statistics)

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EVENT DETECTION - MAJOR ISOFORMS v4.0                       ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Configuration
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v4.0_reference_based"
output_dir <- file.path(base_dir, "major_isoforms")
BATCH_SIZE <- 100  # Process 100 genes per batch
CHECKPOINT_BATCH <- NA  # Set to 1 to stop after first batch, NA to run all genes

# Event detection thresholds
TSS_TES_THRESHOLD <- 20  # Minimum bp difference to report Alt_TSS/Alt_TES (avoids trivial differences)

# ============================================================================
# Helper Functions
# ============================================================================

#' Get union exon numbers present in an isoform
get_isoform_union_exons <- function(union_exons, isoform_id) {
  present_exons <- integer()

  for (i in seq_along(union_exons)) {
    exon_group <- union_exons[[i]]
    if (isoform_id %in% exon_group$isoform_id) {
      present_exons <- c(present_exons, i)
    }
  }

  return(present_exons)
}

#' Get exon details for an isoform from union model
get_isoform_exon_details <- function(union_exons, isoform_id) {
  exon_details <- list()

  for (i in seq_along(union_exons)) {
    exon_group <- union_exons[[i]]
    iso_row <- exon_group %>% filter(isoform_id == !!isoform_id)

    if (nrow(iso_row) > 0) {
      exon_details[[length(exon_details) + 1]] <- list(
        union_exon_number = i,
        exon_type = iso_row$exon_type[1],
        start = iso_row$start[1],
        end = iso_row$end[1],
        is_first = iso_row$is_first[1],
        is_last = iso_row$is_last[1],
        transcript_exon_number = iso_row$transcript_exon_number[1]
      )
    }
  }

  return(bind_rows(exon_details))
}

#' Extract junctions (splice sites) from isoform exons
#' Returns junctions with strand-aware donor/acceptor annotation
extract_junctions <- function(exon_details, gene_strand) {
  if (nrow(exon_details) < 2) {
    return(tibble(
      donor = integer(),
      acceptor = integer(),
      exon_upstream = integer(),
      exon_downstream = integer()
    ))
  }

  # Sort by transcript order
  exons_sorted <- exon_details %>% arrange(transcript_exon_number)

  junctions <- list()

  for (i in 1:(nrow(exons_sorted) - 1)) {
    exon_current <- exons_sorted[i, ]
    exon_next <- exons_sorted[i + 1, ]

    # Junction coordinates depend on strand
    # Plus strand:  donor = end of current exon, acceptor = start of next exon
    # Minus strand: donor = start of next exon, acceptor = end of current exon
    if (gene_strand == "+" || gene_strand == 1) {
      donor <- exon_current$end
      acceptor <- exon_next$start
    } else {
      # Minus strand: donor and acceptor are swapped
      donor <- exon_next$start
      acceptor <- exon_current$end
    }

    junctions[[length(junctions) + 1]] <- tibble(
      donor = donor,
      acceptor = acceptor,
      exon_upstream = exon_current$union_exon_number,
      exon_downstream = exon_next$union_exon_number
    )
  }

  return(bind_rows(junctions))
}

#' Compare junctions between two isoforms to detect A5SS/A3SS
detect_splice_site_events <- function(junctions_A, junctions_B, gene_id, gene_strand,
                                      isoform_A_id, isoform_B_id) {
  events <- list()

  # Find junctions unique to each isoform
  only_A <- anti_join(junctions_A, junctions_B, by = c("donor", "acceptor"))
  only_B <- anti_join(junctions_B, junctions_A, by = c("donor", "acceptor"))

  # For each unique junction in A, check if there's a junction in B with same acceptor (A5SS)
  # or same donor (A3SS)
  for (i in seq_len(nrow(only_A))) {
    junc_A <- only_A[i, ]

    # Check for A5SS: same acceptor, different donor
    a5ss_candidates <- only_B %>%
      filter(acceptor == junc_A$acceptor, donor != junc_A$donor)

    if (nrow(a5ss_candidates) > 0) {
      for (j in seq_len(nrow(a5ss_candidates))) {
        junc_B <- a5ss_candidates[j, ]
        events[[length(events) + 1]] <- tibble(
          gene_id = gene_id,
          strand = gene_strand,
          isoform_A = isoform_A_id,
          isoform_B = isoform_B_id,
          event_type = "A5SS",
          union_exon_A = junc_A$exon_upstream,
          union_exon_B = junc_B$exon_upstream,
          position_A = junc_A$donor,
          position_B = junc_B$donor,
          detail = sprintf("A5SS: donor A=%d, B=%d (same acceptor=%d)",
                          junc_A$donor, junc_B$donor, junc_A$acceptor)
        )
      }
    }

    # Check for A3SS: same donor, different acceptor
    a3ss_candidates <- only_B %>%
      filter(donor == junc_A$donor, acceptor != junc_A$acceptor)

    if (nrow(a3ss_candidates) > 0) {
      for (j in seq_len(nrow(a3ss_candidates))) {
        junc_B <- a3ss_candidates[j, ]
        events[[length(events) + 1]] <- tibble(
          gene_id = gene_id,
          strand = gene_strand,
          isoform_A = isoform_A_id,
          isoform_B = isoform_B_id,
          event_type = "A3SS",
          union_exon_A = junc_A$exon_downstream,
          union_exon_B = junc_B$exon_downstream,
          position_A = junc_A$acceptor,
          position_B = junc_B$acceptor,
          detail = sprintf("A3SS: acceptor A=%d, B=%d (same donor=%d)",
                          junc_A$acceptor, junc_B$acceptor, junc_A$donor)
        )
      }
    }
  }

  # ========================================================================
  # SYMMETRY FIX: Also check from B's perspective
  # ========================================================================
  # The above loop only checks junctions unique to A and looks for matches in B.
  # We also need to check junctions unique to B and look for matches in A.
  # This is critical when isoform A has fewer junctions than B (e.g., monoexonic).

  for (i in seq_len(nrow(only_B))) {
    junc_B <- only_B[i, ]

    # Check for A5SS: same acceptor, different donor
    a5ss_candidates <- only_A %>%
      filter(acceptor == junc_B$acceptor, donor != junc_B$donor)

    if (nrow(a5ss_candidates) > 0) {
      for (j in seq_len(nrow(a5ss_candidates))) {
        junc_A <- a5ss_candidates[j, ]
        # Avoid duplicates: only add if not already detected from A's perspective
        # This happens when both junctions are "unique" (no exact match but share donor/acceptor)
        already_detected <- any(sapply(events, function(e) {
          e$event_type == "A5SS" &&
          e$union_exon_A == junc_A$exon_upstream &&
          e$union_exon_B == junc_B$exon_upstream
        }))

        if (!already_detected) {
          events[[length(events) + 1]] <- tibble(
            gene_id = gene_id,
            strand = gene_strand,
            isoform_A = isoform_A_id,
            isoform_B = isoform_B_id,
            event_type = "A5SS",
            union_exon_A = junc_A$exon_upstream,
            union_exon_B = junc_B$exon_upstream,
            position_A = junc_A$donor,
            position_B = junc_B$donor,
            detail = sprintf("A5SS: donor A=%d, B=%d (same acceptor=%d)",
                            junc_A$donor, junc_B$donor, junc_B$acceptor)
          )
        }
      }
    }

    # Check for A3SS: same donor, different acceptor
    a3ss_candidates <- only_A %>%
      filter(donor == junc_B$donor, acceptor != junc_B$acceptor)

    if (nrow(a3ss_candidates) > 0) {
      for (j in seq_len(nrow(a3ss_candidates))) {
        junc_A <- a3ss_candidates[j, ]
        # Avoid duplicates
        already_detected <- any(sapply(events, function(e) {
          e$event_type == "A3SS" &&
          e$union_exon_A == junc_A$exon_downstream &&
          e$union_exon_B == junc_B$exon_downstream
        }))

        if (!already_detected) {
          events[[length(events) + 1]] <- tibble(
            gene_id = gene_id,
            strand = gene_strand,
            isoform_A = isoform_A_id,
            isoform_B = isoform_B_id,
            event_type = "A3SS",
            union_exon_A = junc_A$exon_downstream,
            union_exon_B = junc_B$exon_downstream,
            position_A = junc_A$acceptor,
            position_B = junc_B$acceptor,
            detail = sprintf("A3SS: acceptor A=%d, B=%d (same donor=%d)",
                            junc_A$acceptor, junc_B$acceptor, junc_B$donor)
          )
        }
      }
    }
  }

  if (length(events) > 0) {
    return(bind_rows(events))
  } else {
    return(tibble(
      gene_id = character(),
      strand = character(),
      isoform_A = character(),
      isoform_B = character(),
      event_type = character(),
      union_exon_A = integer(),
      union_exon_B = integer(),
      position_A = integer(),
      position_B = integer(),
      detail = character()
    ))
  }
}

#' Compare two isoforms and detect events
detect_events <- function(union_exons, isoform_A_id, isoform_B_id, gene_id, gene_strand) {
  # Get union exons present in each isoform
  exons_A_numbers <- get_isoform_union_exons(union_exons, isoform_A_id)
  exons_B_numbers <- get_isoform_union_exons(union_exons, isoform_B_id)

  # Get detailed exon information
  exons_A <- get_isoform_exon_details(union_exons, isoform_A_id)
  exons_B <- get_isoform_exon_details(union_exons, isoform_B_id)

  # Identify exons unique to each isoform
  only_A <- setdiff(exons_A_numbers, exons_B_numbers)
  only_B <- setdiff(exons_B_numbers, exons_A_numbers)
  shared <- intersect(exons_A_numbers, exons_B_numbers)

  events <- list()

  # ========================================================================
  # Terminal events (TSS/TES differences)
  # ========================================================================

  # Check first exons (TSS)
  first_A <- exons_A %>% filter(is_first) %>% slice(1)
  first_B <- exons_B %>% filter(is_first) %>% slice(1)

  if (nrow(first_A) > 0 && nrow(first_B) > 0) {
    # Strand-aware TSS position comparison
    # Plus strand: TSS = start; Minus strand: TSS = end
    if (gene_strand == "+" || gene_strand == 1) {
      pos_A <- first_A$start
      pos_B <- first_B$start
    } else {
      pos_A <- first_A$end
      pos_B <- first_B$end
    }

    # Report Alt_TSS only if positions differ by more than threshold
    if (abs(pos_A - pos_B) > TSS_TES_THRESHOLD) {
      events[[length(events) + 1]] <- tibble(
        gene_id = gene_id,
        strand = gene_strand,
        isoform_A = isoform_A_id,
        isoform_B = isoform_B_id,
        event_type = "Alt_TSS",
        union_exon_A = first_A$union_exon_number,
        union_exon_B = first_B$union_exon_number,
        position_A = pos_A,
        position_B = pos_B,
        detail = sprintf("Different TSS: A=%d, B=%d (diff=%d bp)",
                        pos_A, pos_B, abs(pos_A - pos_B))
      )
    }
  }

  # Check last exons (TES)
  last_A <- exons_A %>% filter(is_last) %>% slice(1)
  last_B <- exons_B %>% filter(is_last) %>% slice(1)

  if (nrow(last_A) > 0 && nrow(last_B) > 0) {
    # Strand-aware TES position comparison
    # Plus strand: TES = end; Minus strand: TES = start
    if (gene_strand == "+" || gene_strand == 1) {
      pos_A <- last_A$end
      pos_B <- last_B$end
    } else {
      pos_A <- last_A$start
      pos_B <- last_B$start
    }

    # Report Alt_TES only if positions differ by more than threshold
    if (abs(pos_A - pos_B) > TSS_TES_THRESHOLD) {
      events[[length(events) + 1]] <- tibble(
        gene_id = gene_id,
        strand = gene_strand,
        isoform_A = isoform_A_id,
        isoform_B = isoform_B_id,
        event_type = "Alt_TES",
        union_exon_A = last_A$union_exon_number,
        union_exon_B = last_B$union_exon_number,
        position_A = pos_A,
        position_B = pos_B,
        detail = sprintf("Different TES: A=%d, B=%d (diff=%d bp)",
                        pos_A, pos_B, abs(pos_A - pos_B))
      )
    }
  }

  # ========================================================================
  # Junction-based events (A5SS/A3SS)
  # ========================================================================

  # Extract junctions (splice sites) for both isoforms
  junctions_A <- extract_junctions(exons_A, gene_strand)
  junctions_B <- extract_junctions(exons_B, gene_strand)

  # Detect A5SS and A3SS by comparing junctions
  splice_site_events <- detect_splice_site_events(
    junctions_A, junctions_B,
    gene_id, gene_strand,
    isoform_A_id, isoform_B_id
  )

  if (nrow(splice_site_events) > 0) {
    for (i in 1:nrow(splice_site_events)) {
      events[[length(events) + 1]] <- splice_site_events[i, ]
    }
  }

  # ========================================================================
  # Internal splicing events (Exon skipping)
  # ========================================================================

  # Exon skipping: union exons unique to one isoform (internal only)
  internal_only_A <- exons_A %>%
    filter(union_exon_number %in% only_A, !is_first, !is_last)

  if (nrow(internal_only_A) > 0) {
    for (i in 1:nrow(internal_only_A)) {
      events[[length(events) + 1]] <- tibble(
        gene_id = gene_id,
        strand = gene_strand,
        isoform_A = isoform_A_id,
        isoform_B = isoform_B_id,
        event_type = "SE",
        union_exon_A = internal_only_A$union_exon_number[i],
        union_exon_B = NA_integer_,
        position_A = internal_only_A$start[i],
        position_B = NA_integer_,
        detail = sprintf("Exon %d skipped in B",
                        internal_only_A$union_exon_number[i])
      )
    }
  }

  internal_only_B <- exons_B %>%
    filter(union_exon_number %in% only_B, !is_first, !is_last)

  if (nrow(internal_only_B) > 0) {
    for (i in 1:nrow(internal_only_B)) {
      events[[length(events) + 1]] <- tibble(
        gene_id = gene_id,
        strand = gene_strand,
        isoform_A = isoform_A_id,
        isoform_B = isoform_B_id,
        event_type = "SE",
        union_exon_A = NA_integer_,
        union_exon_B = internal_only_B$union_exon_number[i],
        position_A = NA_integer_,
        position_B = internal_only_B$start[i],
        detail = sprintf("Exon %d skipped in A",
                        internal_only_B$union_exon_number[i])
      )
    }
  }

  # NOTE: A5SS and A3SS events are now detected by junction-based comparison above.
  # Coordinate-based detection has been removed to prevent double-counting.

  # Return all events
  if (length(events) > 0) {
    return(bind_rows(events))
  } else {
    return(tibble(
      gene_id = gene_id,
      strand = gene_strand,
      isoform_A = isoform_A_id,
      isoform_B = isoform_B_id,
      event_type = character(),
      union_exon_A = integer(),
      union_exon_B = integer(),
      position_A = integer(),
      position_B = integer(),
      detail = character()
    ))
  }
}

# ============================================================================
# Load Data
# ============================================================================

cat("═══ Loading Data ═══\n\n")

# Load union models
cat("Loading union exon models...\n")
union_models <- readRDS(file.path(output_dir, "union_exon_models_major.rds"))
cat("  Genes with union models:", length(union_models), "\n\n")

# Load exon structures to get strand information
cat("Loading exon structures for strand info...\n")
exon_structures <- readRDS(file.path(output_dir, "exon_structures_major_isoforms.rds"))

# Create gene -> strand mapping
gene_strand_map <- exon_structures %>%
  select(gene_id, strand) %>%
  distinct()

cat("  Strand info for", nrow(gene_strand_map), "genes\n\n")

# Load major isoforms expression data
cat("Loading major isoforms expression data...\n")
major_isoforms <- readRDS(file.path(output_dir, "major_isoforms_dmso.rds"))
cat("  Rows:", nrow(major_isoforms), "\n")
cat("  Unique genes:", length(unique(major_isoforms$gene_id)), "\n")
cat("  Unique isoforms:", length(unique(major_isoforms$isoform_id)), "\n\n")

# ============================================================================
# Generate Isoform Pairs
# ============================================================================

cat("═══ Generating Isoform Pairs ═══\n\n")

all_pairs <- list()

for (gene_id in names(union_models)) {
  model <- union_models[[gene_id]]

  # Get all isoforms for this gene from union model
  isoform_ids <- unique(unlist(map(model$union_exons, ~.x$isoform_id)))
  n_isoforms <- length(isoform_ids)

  # Generate all pairwise combinations
  if (n_isoforms >= 2) {
    pairs <- combn(isoform_ids, 2, simplify = FALSE)

    for (pair in pairs) {
      all_pairs[[length(all_pairs) + 1]] <- tibble(
        gene_id = gene_id,
        isoform_A = pair[1],
        isoform_B = pair[2],
        n_isoforms_in_gene = n_isoforms,
        n_union_exons = model$n_union_exons
      )
    }
  }
}

isoform_pairs <- bind_rows(all_pairs)

cat("Total isoform pairs:", nrow(isoform_pairs), "\n")
cat("Genes:", length(unique(isoform_pairs$gene_id)), "\n\n")

# ============================================================================
# Detect Events for All Pairs
# ============================================================================

cat("═══ Detecting Events ═══\n\n")

# Group pairs by gene for efficient processing
pairs_by_gene <- isoform_pairs %>%
  group_by(gene_id) %>%
  nest() %>%
  ungroup()

cat("Processing", nrow(pairs_by_gene), "genes...\n\n")

all_events <- list()
genes_processed <- 0
start_time <- Sys.time()

for (i in 1:nrow(pairs_by_gene)) {
  gene_id <- pairs_by_gene$gene_id[i]
  gene_pairs <- pairs_by_gene$data[[i]]

  model <- union_models[[gene_id]]

  # Get strand for this gene
  gene_strand <- gene_strand_map %>%
    filter(gene_id == !!gene_id) %>%
    pull(strand)

  if (length(gene_strand) == 0) {
    cat(sprintf("WARNING: No strand info for gene %s, skipping\n", gene_id))
    next
  }
  gene_strand <- gene_strand[1]

  # Detect events for all pairs in this gene
  for (j in 1:nrow(gene_pairs)) {
    pair <- gene_pairs[j, ]

    events <- detect_events(
      union_exons = model$union_exons,
      isoform_A_id = pair$isoform_A,
      isoform_B_id = pair$isoform_B,
      gene_id = gene_id,
      gene_strand = gene_strand
    )

    # Add pair metadata
    if (nrow(events) > 0) {
      events <- events %>%
        mutate(
          n_isoforms_in_gene = pair$n_isoforms_in_gene,
          n_union_exons = pair$n_union_exons,
          pair_id = paste(gene_id, pair$isoform_A, pair$isoform_B, sep = "|")
        )

      all_events[[length(all_events) + 1]] <- events
    } else {
      # No events detected - record as no difference
      all_events[[length(all_events) + 1]] <- tibble(
        gene_id = gene_id,
        strand = gene_strand,
        isoform_A = pair$isoform_A,
        isoform_B = pair$isoform_B,
        event_type = "NO_DIFF",
        union_exon_A = NA_integer_,
        union_exon_B = NA_integer_,
        position_A = NA_integer_,
        position_B = NA_integer_,
        detail = "No events detected (differences below thresholds or complex structural variation)",
        n_isoforms_in_gene = pair$n_isoforms_in_gene,
        n_union_exons = pair$n_union_exons,
        pair_id = paste(gene_id, pair$isoform_A, pair$isoform_B, sep = "|")
      )
    }
  }

  genes_processed <- genes_processed + 1

  # Progress reporting
  if (genes_processed %% 100 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    eta <- (elapsed / genes_processed) * (nrow(pairs_by_gene) - genes_processed)
    cat(sprintf("[%d/%d] Processed %s | Elapsed: %.1f min | ETA: %.1f min\n",
                genes_processed, nrow(pairs_by_gene), gene_id, elapsed, eta))
  }

  # CHECKPOINT: Save after first batch
  if (!is.na(CHECKPOINT_BATCH) && genes_processed == CHECKPOINT_BATCH * BATCH_SIZE) {
    cat("\n═══ CHECKPOINT: First Batch Complete ═══\n\n")

    events_checkpoint <- bind_rows(all_events)
    checkpoint_file <- file.path(output_dir, "events_checkpoint_batch1.rds")
    saveRDS(events_checkpoint, checkpoint_file)

    cat("Saved checkpoint:", checkpoint_file, "\n")
    cat("Events detected so far:", nrow(events_checkpoint), "\n")
    cat("Genes processed:", genes_processed, "/", nrow(pairs_by_gene), "\n\n")

    # Event summary for checkpoint
    cat("Event type distribution:\n")
    event_summary <- events_checkpoint %>%
      count(event_type) %>%
      arrange(desc(n))
    print(event_summary)
    cat("\n")

    cat("STOPPING AT CHECKPOINT FOR VERIFICATION\n")
    cat("Resume by commenting out this break or setting CHECKPOINT_BATCH = NA\n\n")
    break
  }
}

cat("\nEvent detection complete\n")
cat("Total events detected:", length(all_events), "\n\n")

# ============================================================================
# Combine and Save Results
# ============================================================================

cat("═══ Saving Results ═══\n\n")

all_events_df <- bind_rows(all_events)

# Save full results
output_file <- file.path(output_dir, "isoform_pairs_events.rds")
saveRDS(all_events_df, output_file)

cat("Saved:", output_file, "\n")
cat("Total rows:", nrow(all_events_df), "\n")
cat("Unique pairs:", length(unique(all_events_df$pair_id)), "\n\n")

# Event summary
cat("═══ Event Summary ═══\n\n")

event_counts <- all_events_df %>%
  count(event_type) %>%
  arrange(desc(n))

print(event_counts)

# Save summary
summary_file <- file.path(output_dir, "event_summary.tsv")
write_tsv(event_counts, summary_file)

cat("\nSaved summary:", summary_file, "\n\n")

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EVENT DETECTION COMPLETE                                     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")
