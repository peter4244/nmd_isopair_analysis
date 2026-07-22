#!/usr/bin/env Rscript

# ============================================================================
# Event Pattern Analysis - Major Isoforms
# ============================================================================
#
# This script analyzes patterns in splicing events detected between major
# isoform pairs, including:
#   1. Event co-occurrence (correlation analysis)
#   2. Distance relationships (genomic and exon-based)
#   3. Topological relationships (event ordering)
#
# Input:
#   - events_checkpoint_batch1.rds or isoform_pairs_events.rds
#   - union_exon_models_major.rds (for exon position mapping)
#
# Output:
#   - Event co-occurrence matrix
#   - Distance distributions by event pair type
#   - Event ordering patterns
#   - Summary statistics and visualizations

library(tidyverse)

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EVENT PATTERN ANALYSIS - MAJOR ISOFORMS                     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# ============================================================================
# Setup
# ============================================================================

output_dir <- "results/isoform_transitions/v4.0_reference_based/major_isoforms"

# Load event data (use checkpoint for now, will use full data when available)
events_file <- file.path(output_dir, "isoform_pairs_events.rds")
if (!file.exists(events_file)) {
  events_file <- file.path(output_dir, "events_checkpoint_batch1.rds")
  cat("Using checkpoint data:", events_file, "\n")
} else {
  cat("Using full event data:", events_file, "\n")
}

events <- readRDS(events_file)

cat("Loaded", nrow(events), "events from",
    n_distinct(events$gene_id), "genes\n")
cat("Isoform pairs:", n_distinct(events$pair_id), "\n\n")

# Load union exon models for positional context
union_models <- readRDS(file.path(output_dir, "union_exon_models_major.rds"))

# ============================================================================
# Part 1: Event Aggregation by Isoform Pair
# ============================================================================

cat("═══ Part 1: Aggregating Events by Isoform Pair ═══\n\n")

# Create event vectors: one row per pair with all events
pair_events <- events %>%
  group_by(pair_id, gene_id, strand, isoform_A, isoform_B,
           n_isoforms_in_gene, n_union_exons) %>%
  summarise(
    # Count each event type
    n_alt_tss = sum(event_type == "Alt_TSS"),
    n_alt_tes = sum(event_type == "Alt_TES"),
    n_se = sum(event_type == "SE"),
    n_a5ss = sum(event_type == "A5SS"),
    n_a3ss = sum(event_type == "A3SS"),

    # Total events
    n_total_events = n(),

    # Event type combinations (for pattern analysis)
    event_types = paste(sort(unique(event_type)), collapse = ","),

    # Store nested event details for distance calculations
    events = list(tibble(
      event_type = event_type,
      union_exon_A = union_exon_A,
      union_exon_B = union_exon_B,
      position_A = position_A,
      position_B = position_B
    )),

    .groups = "drop"
  )

cat("Created event vectors for", nrow(pair_events), "isoform pairs\n\n")

# Summary statistics
cat("Event complexity distribution:\n")
print(pair_events %>% count(n_total_events) %>% arrange(desc(n)))
cat("\n")

# ============================================================================
# Part 2: Event Co-occurrence Analysis (Q2)
# ============================================================================

cat("═══ Part 2: Event Co-occurrence Analysis ═══\n\n")

# Create binary presence/absence for each event type
pair_events_binary <- pair_events %>%
  mutate(
    has_alt_tss = n_alt_tss > 0,
    has_alt_tes = n_alt_tes > 0,
    has_se = n_se > 0,
    has_a5ss = n_a5ss > 0,
    has_a3ss = n_a3ss > 0
  )

# All pairwise combinations of event types
event_types <- c("alt_tss", "alt_tes", "se", "a5ss", "a3ss")
event_pairs_combinations <- combn(event_types, 2, simplify = FALSE)

# Test co-occurrence for each pair
cooccurrence_results <- map_dfr(event_pairs_combinations, function(pair) {
  event_A <- pair[1]
  event_B <- pair[2]

  has_A <- pair_events_binary[[paste0("has_", event_A)]]
  has_B <- pair_events_binary[[paste0("has_", event_B)]]

  # Build 2x2 contingency table
  #           Event B present  Event B absent
  # Event A present      n_both       n_A_only
  # Event A absent      n_B_only      n_neither
  n_both <- sum(has_A & has_B)
  n_A_only <- sum(has_A & !has_B)
  n_B_only <- sum(!has_A & has_B)
  n_neither <- sum(!has_A & !has_B)

  contingency <- matrix(c(n_both, n_B_only, n_A_only, n_neither),
                        nrow = 2, byrow = TRUE)

  # Fisher's exact test
  test_result <- fisher.test(contingency)

  tibble(
    event_A = event_A,
    event_B = event_B,
    n_both = n_both,
    n_A_only = n_A_only,
    n_B_only = n_B_only,
    n_neither = n_neither,
    n_total = nrow(pair_events_binary),
    prop_both = n_both / nrow(pair_events_binary),
    odds_ratio = as.numeric(test_result$estimate),
    p_value = test_result$p.value
  )
})

# Apply multiple testing correction (Benjamini-Hochberg FDR)
cooccurrence_results <- cooccurrence_results %>%
  mutate(
    adj_p_value = p.adjust(p_value, method = "BH"),
    significant = adj_p_value < 0.05,
    correlation_type = case_when(
      adj_p_value >= 0.05 ~ "non-correlated",
      odds_ratio > 1 ~ "positively correlated",
      odds_ratio < 1 ~ "negatively correlated",
      TRUE ~ "non-correlated"
    )
  )

cat("Co-occurrence test results (with FDR correction):\n")
print(cooccurrence_results %>%
        select(event_A, event_B, n_both, odds_ratio, p_value, adj_p_value, correlation_type))
cat("\n")

# Separate by correlation type
positive_pairs <- cooccurrence_results %>%
  filter(correlation_type == "positively correlated")
negative_pairs <- cooccurrence_results %>%
  filter(correlation_type == "negatively correlated")
noncorr_pairs <- cooccurrence_results %>%
  filter(correlation_type == "non-correlated")

cat("Positively correlated pairs:", nrow(positive_pairs), "\n")
if (nrow(positive_pairs) > 0) {
  print(positive_pairs %>% select(event_A, event_B, odds_ratio, p_value, adj_p_value))
}
cat("\n")

cat("Negatively correlated pairs:", nrow(negative_pairs), "\n")
if (nrow(negative_pairs) > 0) {
  print(negative_pairs %>% select(event_A, event_B, odds_ratio, p_value, adj_p_value))
}
cat("\n")

# Save co-occurrence results
write_tsv(cooccurrence_results,
          file.path(output_dir, "event_cooccurrence_matrix.tsv"))
cat("Saved: event_cooccurrence_matrix.tsv\n\n")

# ============================================================================
# Part 3: Distance Analysis (Q3)
# ============================================================================

cat("═══ Part 3: Distance Analysis ═══\n\n")

# Function to calculate pairwise distances between events
calculate_event_distances <- function(event_A_type, event_B_type) {

  # Filter to pairs with both event types
  col_A <- paste0("n_", event_A_type)
  col_B <- paste0("n_", event_B_type)

  pairs_with_both <- pair_events %>%
    filter(.data[[col_A]] > 0 & .data[[col_B]] > 0)

  if (nrow(pairs_with_both) == 0) {
    return(tibble(
      pair_id = character(),
      gene_id = character(),
      strand = factor(),
      genomic_distance = integer(),
      exon_distance = integer(),
      event_pair = character(),
      event_A = character(),
      event_B = character()
    ))
  }

  # Standardize event type names for matching
  # Convert "alt_tss" -> "Alt_TSS", "se" -> "SE", "a5ss" -> "A5SS", etc.
  standardize_event_name <- function(name) {
    name_upper <- toupper(name)
    # Check if it starts with "ALT"
    if (grepl("^ALT", name_upper)) {
      # Extract the suffix (TSS, TES)
      suffix <- gsub("ALT_?", "", name_upper)
      return(paste0("Alt_", suffix))
    } else {
      # Return as uppercase (SE, A5SS, A3SS, IR)
      return(name_upper)
    }
  }

  event_A_std <- standardize_event_name(event_A_type)
  event_B_std <- standardize_event_name(event_B_type)

  # Calculate distances
  all_dist <- list()

  for (i in 1:nrow(pairs_with_both)) {
    row <- pairs_with_both[i, ]
    ev <- row$events[[1]]

    # Match event types using standardized names
    events_A <- ev %>% filter(event_type == event_A_std)
    events_B <- ev %>% filter(event_type == event_B_std)

    if (nrow(events_A) == 0 || nrow(events_B) == 0) {
      next
    }

    # All pairwise combinations
    for (j in 1:nrow(events_A)) {
      for (k in 1:nrow(events_B)) {
        pos_A <- events_A$position_A[j]
        pos_B <- events_B$position_A[k]

        # Skip if positions are NA
        if (is.na(pos_A) || is.na(pos_B)) next

        # Get exon numbers, using either A or B column (whichever is not NA)
        exon_A_num <- coalesce(events_A$union_exon_A[j], events_A$union_exon_B[j])
        exon_B_num <- coalesce(events_B$union_exon_A[k], events_B$union_exon_B[k])

        # Calculate exon distance only if both are available
        exon_dist <- if (!is.na(exon_A_num) && !is.na(exon_B_num)) {
          abs(exon_A_num - exon_B_num)
        } else {
          NA_integer_
        }

        all_dist[[length(all_dist) + 1]] <- tibble(
          pair_id = row$pair_id,
          gene_id = row$gene_id,
          strand = row$strand,
          genomic_distance = abs(pos_A - pos_B),
          exon_distance = exon_dist,
          event_pair = paste(event_A_type, event_B_type, sep = "-"),
          event_A = event_A_type,
          event_B = event_B_type
        )
      }
    }
  }

  if (length(all_dist) == 0) {
    return(tibble(
      pair_id = character(),
      gene_id = character(),
      strand = factor(),
      genomic_distance = integer(),
      exon_distance = integer(),
      event_pair = character(),
      event_A = character(),
      event_B = character()
    ))
  }

  bind_rows(all_dist)
}

# Calculate distances for all correlated pairs
cat("Calculating distances for correlated event pairs...\n")

all_distances <- bind_rows(
  # Positively correlated pairs
  map_dfr(1:nrow(positive_pairs), function(i) {
    calculate_event_distances(
      positive_pairs$event_A[i],
      positive_pairs$event_B[i]
    ) %>%
      mutate(correlation_type = "positive")
  }),

  # Negatively correlated pairs (only if any exist)
  if (nrow(negative_pairs) > 0) {
    map_dfr(1:nrow(negative_pairs), function(i) {
      calculate_event_distances(
        negative_pairs$event_A[i],
        negative_pairs$event_B[i]
      ) %>%
        mutate(correlation_type = "negative")
    })
  } else {
    tibble()
  },

  # Non-correlated pairs (for comparison)
  map_dfr(1:min(3, nrow(noncorr_pairs)), function(i) {
    calculate_event_distances(
      noncorr_pairs$event_A[i],
      noncorr_pairs$event_B[i]
    ) %>%
      mutate(correlation_type = "non-correlated")
  })
)

cat("Calculated", nrow(all_distances), "pairwise event distances\n\n")

# Summary statistics by correlation type
if (nrow(all_distances) > 0) {
  distance_summary <- all_distances %>%
    group_by(correlation_type, event_pair) %>%
    summarise(
      n_observations = n(),
      mean_genomic_dist = mean(genomic_distance),
      median_genomic_dist = median(genomic_distance),
      mean_exon_dist = mean(exon_distance),
      median_exon_dist = median(exon_distance),
      .groups = "drop"
    )

  cat("Distance summary by correlation type:\n")
  print(distance_summary)
  cat("\n")

  # Save results
  write_tsv(all_distances,
            file.path(output_dir, "event_pair_distances.tsv"))
  write_tsv(distance_summary,
            file.path(output_dir, "distance_summary_by_correlation.tsv"))

  cat("Saved: event_pair_distances.tsv\n")
  cat("Saved: distance_summary_by_correlation.tsv\n\n")
}

# ============================================================================
# Part 4: Topological Relationships (Event Ordering)
# ============================================================================

cat("═══ Part 4: Topological Relationships ═══\n\n")

# For pairs with multiple events, determine ordering
# Use union exon numbers and strand to determine 5' -> 3' order

ordering_analysis <- pair_events %>%
  filter(n_total_events >= 2) %>%
  mutate(
    event_order = map2(events, strand, function(ev, s) {
      # Sort by union exon number (accounting for strand)
      # For minus strand, higher exon number = more 5' (comes first in transcript)
      if (s == "-") {
        ev <- ev %>% arrange(desc(union_exon_A))
      } else {
        ev <- ev %>% arrange(union_exon_A)
      }

      # Create ordered event type sequence
      tibble(
        position = 1:nrow(ev),
        event_type = ev$event_type,
        union_exon = ev$union_exon_A
      )
    })
  ) %>%
  select(pair_id, gene_id, strand, n_total_events, event_order) %>%
  unnest(event_order)

cat("Ordered events for", n_distinct(ordering_analysis$pair_id),
    "pairs with ≥2 events\n\n")

# Most common event orderings
if (nrow(ordering_analysis) > 0) {
  event_patterns <- ordering_analysis %>%
    group_by(pair_id) %>%
    summarise(
      pattern = paste(event_type, collapse = " -> "),
      .groups = "drop"
    ) %>%
    count(pattern, sort = TRUE) %>%
    mutate(proportion = n / sum(n))

  cat("Top 10 event ordering patterns:\n")
  print(event_patterns %>% head(10))
  cat("\n")

  # Save results
  write_tsv(ordering_analysis,
            file.path(output_dir, "event_ordering.tsv"))
  write_tsv(event_patterns,
            file.path(output_dir, "event_patterns_ranked.tsv"))

  cat("Saved: event_ordering.tsv\n")
  cat("Saved: event_patterns_ranked.tsv\n\n")
}

# ============================================================================
# Summary
# ============================================================================

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EVENT PATTERN ANALYSIS COMPLETE                             ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("RESULTS:\n")
cat("  Isoform pairs analyzed:", nrow(pair_events), "\n")
cat("  Co-occurrence tests:", nrow(cooccurrence_results), "\n")
cat("  Positively correlated pairs:", nrow(positive_pairs), "\n")
cat("  Negatively correlated pairs:", nrow(negative_pairs), "\n")
if (nrow(all_distances) > 0) {
  cat("  Distance measurements:", nrow(all_distances), "\n")
}
if (nrow(ordering_analysis) > 0) {
  cat("  Ordered event sequences:", n_distinct(ordering_analysis$pair_id), "\n")
}

cat("\nOUTPUT FILES:\n")
cat("  event_cooccurrence_matrix.tsv\n")
if (nrow(all_distances) > 0) {
  cat("  event_pair_distances.tsv\n")
  cat("  distance_summary_by_correlation.tsv\n")
}
if (nrow(ordering_analysis) > 0) {
  cat("  event_ordering.tsv\n")
  cat("  event_patterns_ranked.tsv\n")
}
cat("\n")
