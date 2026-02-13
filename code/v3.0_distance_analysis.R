#!/usr/bin/env Rscript
# v3.0 Distance Analysis - Structure-Stratified Topological Test
# Tests whether correlated events show enriched/depleted topological configurations
# Uses exact mathematical null distribution with binomial tests

library(tidyverse)

# Load exact topology test function
source("code/exact_topology_test.R")

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   v3.0 DISTANCE ANALYSIS - TOPOLOGICAL PERMUTATION TEST       ║\n")
cat("║   Structure-Stratified Analysis with f2f/b2b Orientation      ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

output_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v3.0_reference_based"

# ============================================================================
# Load Data
# ============================================================================

cat("Loading v3.0 filtered event data...\n")
transitions <- readRDS(file.path(output_dir, "reference_event_vectors_v3.0_filtered.rds"))
cat("Loaded:", nrow(transitions), "transitions\n")

cat("Loading correlation results...\n")
correlations <- read_tsv(file.path(output_dir, "reference_cooccurrence_all.tsv"),
                        show_col_types = FALSE)

# Identify event pairs to test (positive and negative correlations, exclude IR)
event_pairs_to_test <- correlations %>%
  filter(significant == TRUE) %>%  # p < 0.05
  filter(!str_detect(event_A, "IR") & !str_detect(event_B, "IR")) %>%  # No IR
  # Exclude pairs where BOTH events are terminal (at least one must be internal)
  filter(event_A_base %in% c("A5SS", "A3SS", "SE") |
         event_B_base %in% c("A5SS", "A3SS", "SE")) %>%
  select(event_A, event_B, correlation_type, odds_ratio, p_value)

cat("Event pairs to test:", nrow(event_pairs_to_test), "\n")
cat("  Positive correlations:", sum(event_pairs_to_test$correlation_type == "positively correlated"), "\n")
cat("  Negative correlations:", sum(event_pairs_to_test$correlation_type == "negatively correlated"), "\n\n")

# ============================================================================
# Helper Functions
# ============================================================================

# Classify topological relationship for A5SS + A3SS
classify_topology_a5ss_a3ss <- function(a5ss_exon, a3ss_exon) {
  distance <- abs(a5ss_exon - a3ss_exon)

  if (distance == 0) {
    return("same_exon")
  } else if (a3ss_exon < a5ss_exon) {
    # Face to face (A3SS upstream, A5SS downstream)
    return(paste0("dist_", distance, "_f2f"))
  } else {
    # Back to back (A5SS upstream, A3SS downstream)
    return(paste0("dist_", distance, "_b2b"))
  }
}

# Classify topological relationship for events involving SE
# (no f2f/b2b distinction - SE is about exon inclusion, not boundaries)
classify_topology_with_se <- function(exon_a, exon_b) {
  distance <- abs(exon_a - exon_b)

  if (distance == 0) {
    return("same_exon")
  } else {
    return(paste0("dist_", distance))
  }
}

# Extract exon positions for events from event vector
extract_event_positions <- function(event_vector, event_type) {
  if (is.null(event_vector) || length(event_vector) == 0) {
    return(integer(0))
  }

  # Filter to events of this type
  events <- event_vector %>%
    filter(event_type == !!event_type)

  if (nrow(events) == 0) {
    return(integer(0))
  }

  # Return exon numbers
  return(events$exon_number)
}

# Define valid positions for event types
get_valid_positions <- function(event_type, n_exons) {
  base_type <- str_remove(event_type, "_(gain|loss)$")

  if (base_type == "A5SS") {
    # A5SS at downstream boundaries: all exons EXCEPT last (E1, E2, ..., En-1)
    return(1:(n_exons - 1))
  } else if (base_type == "A3SS") {
    # A3SS at upstream boundaries: all exons EXCEPT first (E2, E3, ..., En)
    return(2:n_exons)
  } else if (base_type == "SE") {
    # SE at internal exons only (E2, ..., En-1)
    return(2:(n_exons - 1))
  } else {
    # Terminal events - not permuted
    return(integer(0))
  }
}

# ============================================================================
# Main Analysis Loop
# ============================================================================

cat("═══ Running Structure-Stratified Topological Tests ═══\n\n")

all_results <- list()
all_chisq_results <- list()
MIN_N <- 50  # Minimum observations per structure
# Using exact mathematical null distribution with binomial tests (per topology)
# and chi-square tests (omnibus test per structure)

for (pair_idx in seq_len(nrow(event_pairs_to_test))) {

  event_A <- event_pairs_to_test$event_A[pair_idx]
  event_B <- event_pairs_to_test$event_B[pair_idx]
  corr_type <- event_pairs_to_test$correlation_type[pair_idx]

  cat(sprintf("[%d/%d] Testing: %s + %s (%s)\n",
              pair_idx, nrow(event_pairs_to_test),
              event_A, event_B, corr_type))

  # Filter to transitions with both events
  has_A <- transitions[[paste0("has_", event_A)]]
  has_B <- transitions[[paste0("has_", event_B)]]

  trans_subset <- transitions %>%
    filter(has_A & has_B)

  if (nrow(trans_subset) == 0) {
    cat("  No observations with both events, skipping.\n\n")
    next
  }

  cat(sprintf("  Total observations: %d\n", nrow(trans_subset)))

  # Add n_exons to each transition (from event_vector)
  # IMPORTANT: exon_number means different things for different events:
  #   - Alt_TSS/Alt_TES: actual exon index
  #   - A5SS/A3SS/CONST: junction index (1 to n-1)
  # Correct calculation: use Alt_TES exon_number OR max_junction + 1
  trans_subset <- trans_subset %>%
    mutate(n_exons = map_int(event_vector, function(ev) {
      if (is.null(ev) || nrow(ev) == 0) return(NA_integer_)

      # If Alt_TES exists, use its exon_number (actual last exon)
      tes_events <- ev %>% filter(event_type == "Alt_TES")
      if (nrow(tes_events) > 0) {
        return(as.integer(tes_events$exon_number[1]))
      }

      # Otherwise, use max junction index + 1
      junction_events <- ev %>% filter(event_type %in% c("A5SS", "A3SS", "CONST", "SE"))
      if (nrow(junction_events) > 0) {
        max_junction <- max(junction_events$exon_number, na.rm = TRUE)
        if (is.infinite(max_junction)) return(NA_integer_)
        return(as.integer(max_junction + 1))
      }

      # If only Alt_TSS exists, assume 1 exon (though unusual)
      return(1L)
    })) %>%
    filter(!is.na(n_exons))

  # Stratify by n_exons
  exon_structures <- trans_subset %>%
    count(n_exons) %>%
    filter(n >= MIN_N)

  if (nrow(exon_structures) == 0) {
    cat("  No structures with n >= 50, skipping.\n\n")
    next
  }

  cat(sprintf("  Structures with n>=50: %d\n", nrow(exon_structures)))

  # Process each structure
  for (struct_idx in seq_len(nrow(exon_structures))) {

    n_exons_val <- exon_structures$n_exons[struct_idx]
    n_obs <- exon_structures$n[struct_idx]

    cat(sprintf("    Structure: %d exons (n=%d)\n", n_exons_val, n_obs))

    # Get observations for this structure
    obs_struct <- trans_subset %>%
      filter(n_exons == n_exons_val)

    # Extract base event types once
    event_A_base <- str_remove(event_A, "_(gain|loss)$")
    event_B_base <- str_remove(event_B, "_(gain|loss)$")

    # Extract observed topologies
    # Handle ALL pairwise combinations for multiple events
    observed_topologies <- map(obs_struct$event_vector, function(ev) {
      # Extract positions for event A and B
      pos_A <- extract_event_positions(ev, event_A_base)
      pos_B <- extract_event_positions(ev, event_B_base)

      if (length(pos_A) == 0 || length(pos_B) == 0) {
        return(character(0))
      }

      # ALL PAIRWISE COMBINATIONS (pre-allocate for efficiency)
      n_combinations <- length(pos_A) * length(pos_B)
      topologies <- character(n_combinations)
      idx <- 1

      for (i in seq_along(pos_A)) {
        for (j in seq_along(pos_B)) {
          # Determine which classification function to use
          if (event_A_base == "A5SS" && event_B_base == "A3SS") {
            # A5SS + A3SS: use f2f/b2b classification
            topo <- classify_topology_a5ss_a3ss(pos_A[i], pos_B[j])
          } else if (event_A_base == "A3SS" && event_B_base == "A5SS") {
            # A3SS + A5SS: swap order for correct f2f/b2b
            topo <- classify_topology_a5ss_a3ss(pos_B[j], pos_A[i])
          } else if (event_A_base == "SE" || event_B_base == "SE") {
            # SE + anything: no f2f/b2b, just distance
            topo <- classify_topology_with_se(pos_A[i], pos_B[j])
          } else {
            # Other combinations: default to distance only
            topo <- classify_topology_with_se(pos_A[i], pos_B[j])
          }
          topologies[idx] <- topo
          idx <- idx + 1
        }
      }

      return(topologies)
    })

    # Flatten list of vectors to single vector
    observed_topologies <- unlist(observed_topologies)

    if (length(observed_topologies) == 0) {
      cat("      No valid topologies extracted, skipping.\n")
      next
    }

    obs_counts <- table(observed_topologies)
    n_obs_actual <- length(observed_topologies)  # May be > n_obs due to multiple events

    cat(sprintf("      Extracted %d topology observations from %d transitions\n",
                n_obs_actual, n_obs))

    # Exact topology test using mathematical null distribution
    # Define valid positions for each event type
    valid_A <- get_valid_positions(event_A, n_exons_val)
    valid_B <- get_valid_positions(event_B, n_exons_val)

    # Only skip if BOTH are terminal (no variation possible)
    if (length(valid_A) == 0 && length(valid_B) == 0) {
      cat("      Cannot test (both terminal events, no variation), skipping.\n")
      next
    }

    # Define topology classification function based on event types
    classify_topo <- if (event_A_base == "A5SS" && event_B_base == "A3SS") {
      classify_topology_a5ss_a3ss
    } else if (event_A_base == "A3SS" && event_B_base == "A5SS") {
      function(a, b) classify_topology_a5ss_a3ss(b, a)  # Swap order
    } else {
      classify_topology_with_se
    }

    # Run exact topology test (mathematical null + binomial + chi-square tests)
    cat("      Running exact topology test...\n")
    test_results <- exact_topology_test(
      obs_topologies = observed_topologies,
      obs_struct = obs_struct,
      event_A_base = event_A_base,
      event_B_base = event_B_base,
      valid_A = valid_A,
      valid_B = valid_B,
      extract_event_positions_func = extract_event_positions,
      classify_topology_func = classify_topo
    )

    # Extract topology-specific results (binomial tests)
    topology_results <- test_results$topology_tests %>%
      mutate(
        event_A = event_A,
        event_B = event_B,
        correlation_type = corr_type,
        n_exons = n_exons_val,
        n_transitions = n_obs,  # Original transition count
        n_observations = n_obs_actual,  # Total observations (may be > transitions due to multiple events)
        significance = case_when(
          is.na(p_value) ~ "",
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      ) %>%
      select(event_A, event_B, correlation_type, n_exons, n_transitions, n_observations,
             topology, obs_count, obs_pct, exp_count, exp_pct,
             odds_ratio, p_value, significance)

    # Extract chi-square omnibus test result
    chisq_result <- test_results$chisq_test %>%
      mutate(
        event_A = event_A,
        event_B = event_B,
        correlation_type = corr_type,
        n_exons = n_exons_val,
        n_transitions = n_obs,
        significance_chisq = case_when(
          is.na(p_value_chisq) ~ "",
          p_value_chisq < 0.001 ~ "***",
          p_value_chisq < 0.01 ~ "**",
          p_value_chisq < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      ) %>%
      select(event_A, event_B, correlation_type, n_exons, n_transitions, n_observations,
             n_topologies, chi_squared, df, p_value_chisq, significance_chisq)

    # Store results
    all_results[[length(all_results) + 1]] <- topology_results
    all_chisq_results[[length(all_chisq_results) + 1]] <- chisq_result

    cat(sprintf("      Complete. %d topology categories tested (χ² p=%.4g).\n",
                nrow(topology_results), chisq_result$p_value_chisq))
  }

  cat("\n")
}

# ============================================================================
# Combine and Save Results
# ============================================================================

cat("Combining results...\n")
combined_results <- bind_rows(all_results)
combined_chisq <- bind_rows(all_chisq_results)

cat("Saving results...\n")
write_tsv(combined_results,
          file.path(output_dir, "v3.0_topological_relationships_by_structure.tsv"))
write_tsv(combined_chisq,
          file.path(output_dir, "v3.0_topological_chisquare_tests_by_structure.tsv"))

# Summary statistics
cat("\n═══ Summary ═══\n\n")
cat("Total topology tests:", nrow(combined_results), "\n")
cat("Significant (p<0.05):", sum(combined_results$p_value < 0.05, na.rm = TRUE), "\n")
cat("Enriched (OR>1, p<0.05):", sum(combined_results$odds_ratio > 1 &
                                    combined_results$p_value < 0.05, na.rm = TRUE), "\n")
cat("Depleted (OR<1, p<0.05):", sum(combined_results$odds_ratio < 1 &
                                    combined_results$p_value < 0.05, na.rm = TRUE), "\n\n")

# Top enrichments
cat("Top 10 enriched configurations:\n")
top_enriched <- combined_results %>%
  filter(!is.na(odds_ratio), !is.na(p_value), p_value < 0.05, odds_ratio > 1) %>%
  arrange(desc(odds_ratio)) %>%
  head(10) %>%
  select(event_A, event_B, n_exons, topology, obs_pct, exp_pct, odds_ratio, p_value)

print(top_enriched)

cat("\n")
cat("Chi-square omnibus tests:\n")
cat("  Total structures tested:", nrow(combined_chisq), "\n")
cat("  Significant (p<0.05):", sum(combined_chisq$p_value_chisq < 0.05, na.rm = TRUE), "\n")
cat("  Highly significant (p<0.001):", sum(combined_chisq$p_value_chisq < 0.001, na.rm = TRUE), "\n\n")

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   v3.0 DISTANCE ANALYSIS COMPLETE                             ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("Results saved to:\n")
cat("  - v3.0_topological_relationships_by_structure.tsv (binomial tests)\n")
cat("  - v3.0_topological_chisquare_tests_by_structure.tsv (omnibus tests)\n")
cat("\n")
