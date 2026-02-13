#!/usr/bin/env Rscript
# Detailed Event Analysis: Q1-Q3
# Phase 5: Event Dependency Analysis with specific event types

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║        DETAILED EVENT ANALYSIS: Q1-Q3                          ║\n")
cat("║        Event Incidence, Co-occurrence, and Dependencies       ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

output_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions"

# Load detailed event data
cat("Loading detailed event vectors...\n")
transitions <- readRDS(file.path(output_dir, "detailed_event_vectors.rds"))
cat("Loaded:", nrow(transitions), "transitions\n\n")

# ============================================================================
# Q1: Event Incidence by Cell Type
# ============================================================================

cat("═══ Q1: Event Incidence by Cell Type ═══\n\n")

# Overall incidence
overall_incidence <- transitions %>%
  summarize(
    total_transitions = n(),
    total_alt_tss = sum(n_alt_tss, na.rm = TRUE),
    total_alt_tes = sum(n_alt_tes, na.rm = TRUE),
    total_se = sum(n_se, na.rm = TRUE),
    total_a5ss = sum(n_a5ss, na.rm = TRUE),
    total_a3ss = sum(n_a3ss, na.rm = TRUE),
    total_ir = sum(n_ir, na.rm = TRUE),
    pct_with_alt_tss = 100 * mean(n_alt_tss > 0, na.rm = TRUE),
    pct_with_alt_tes = 100 * mean(n_alt_tes > 0, na.rm = TRUE),
    pct_with_se = 100 * mean(n_se > 0, na.rm = TRUE),
    pct_with_a5ss = 100 * mean(n_a5ss > 0, na.rm = TRUE),
    pct_with_a3ss = 100 * mean(n_a3ss > 0, na.rm = TRUE),
    pct_with_ir = 100 * mean(n_ir > 0, na.rm = TRUE)
  )

cat("Overall event incidence:\n")
print(overall_incidence)
cat("\n")

# By cell type
incidence_by_celltype <- transitions %>%
  group_by(cell_type) %>%
  summarize(
    n_transitions = n(),
    # Calculate percentages FIRST (from original columns)
    pct_with_alt_tss = 100 * mean(n_alt_tss > 0),
    pct_with_alt_tes = 100 * mean(n_alt_tes > 0),
    pct_with_se = 100 * mean(n_se > 0),
    pct_with_a5ss = 100 * mean(n_a5ss > 0),
    pct_with_a3ss = 100 * mean(n_a3ss > 0),
    pct_with_ir = 100 * mean(n_ir > 0),
    # Then calculate sums
    n_alt_tss = sum(n_alt_tss),
    n_alt_tes = sum(n_alt_tes),
    n_se = sum(n_se),
    n_a5ss = sum(n_a5ss),
    n_a3ss = sum(n_a3ss),
    n_ir = sum(n_ir),
    .groups = "drop"
  )

cat("Event incidence by cell type:\n")
print(incidence_by_celltype)

write_tsv(incidence_by_celltype,
          file.path(output_dir, "q1_event_incidence_by_celltype.tsv"))
cat("\nSaved: q1_event_incidence_by_celltype.tsv\n\n")

# ============================================================================
# Q2: Event Co-occurrence and Dependencies
# ============================================================================

cat("═══ Q2: Event Co-occurrence and Dependencies ═══\n\n")

# Pairwise co-occurrence testing
event_types <- c("alt_tss", "alt_tes", "se", "a5ss", "a3ss", "ir")
event_pairs <- combn(event_types, 2, simplify = FALSE)

cooccurrence_results <- map_dfr(event_pairs, function(pair) {
  event_A <- pair[1]
  event_B <- pair[2]
  
  has_A <- transitions[[paste0("n_", event_A)]] > 0
  has_B <- transitions[[paste0("n_", event_B)]] > 0
  
  n_both <- sum(has_A & has_B)
  n_A_only <- sum(has_A & !has_B)
  n_B_only <- sum(!has_A & has_B)
  n_neither <- sum(!has_A & !has_B)
  
  contingency <- matrix(c(n_both, n_A_only, n_B_only, n_neither), nrow = 2)
  test_result <- fisher.test(contingency)
  
  tibble(
    event_A = event_A,
    event_B = event_B,
    n_both = n_both,
    n_A_only = n_A_only,
    n_B_only = n_B_only,
    n_neither = n_neither,
    odds_ratio = test_result$estimate,
    p_value = test_result$p.value,
    significant = p_value < 0.05,
    correlation_type = case_when(
      p_value >= 0.05 ~ "non-correlated",
      odds_ratio > 1 ~ "positively correlated",
      odds_ratio < 1 ~ "negatively correlated"
    )
  )
})

cat("Co-occurrence analysis results:\n")
print(cooccurrence_results)

# Separate by correlation type
positive_pairs <- filter(cooccurrence_results, correlation_type == "positively correlated")
negative_pairs <- filter(cooccurrence_results, correlation_type == "negatively correlated")
noncorr_pairs <- filter(cooccurrence_results, correlation_type == "non-correlated")

cat("\nPositively correlated pairs:", nrow(positive_pairs), "\n")
cat("Negatively correlated pairs:", nrow(negative_pairs), "\n")
cat("Non-correlated pairs:", nrow(noncorr_pairs), "\n\n")

# Save results
write_tsv(cooccurrence_results,
          file.path(output_dir, "q2_event_cooccurrence_matrix.tsv"))
write_tsv(positive_pairs,
          file.path(output_dir, "q2_positively_correlated_pairs.tsv"))
write_tsv(negative_pairs,
          file.path(output_dir, "q2_negatively_correlated_pairs.tsv"))

cat("Saved Q2 results\n\n")

# ============================================================================
# Q3: Distance Analysis for Correlated Event Pairs
# ============================================================================

cat("═══ Q3: Distance Analysis for Correlated Pairs ═══\n\n")

# Function to calculate pairwise distances (OPTIMIZED)
calculate_event_distances <- function(event_A_type, event_B_type, transitions) {
  # Convert event type names to match data format
  event_A_match <- case_when(
    event_A_type == "alt_tss" ~ "Alt_TSS",
    event_A_type == "alt_tes" ~ "Alt_TES",
    TRUE ~ toupper(event_A_type)
  )
  event_B_match <- case_when(
    event_B_type == "alt_tss" ~ "Alt_TSS",
    event_B_type == "alt_tes" ~ "Alt_TES",
    TRUE ~ toupper(event_B_type)
  )

  # Filter to transitions with both events
  filtered_trans <- transitions %>%
    filter(.data[[paste0("n_", event_A_type)]] > 0 &
           .data[[paste0("n_", event_B_type)]] > 0)

  cat("    Found", nrow(filtered_trans), "transitions with both events\n")

  # Sample if too many (for performance)
  if (nrow(filtered_trans) > 50000) {
    cat("    Sampling 50000 transitions for performance\n")
    filtered_trans <- sample_n(filtered_trans, 50000)
  }

  # Extract distances
  distances_list <- vector("list", nrow(filtered_trans))

  for (i in seq_len(nrow(filtered_trans))) {
    ev <- filtered_trans$event_vector[[i]]

    if (is.null(ev) || nrow(ev) == 0) next

    events_A <- filter(ev, event_type == event_A_match)
    events_B <- filter(ev, event_type == event_B_match)

    if (nrow(events_A) == 0 || nrow(events_B) == 0) next

    # Calculate distances for all pairs
    n_pairs <- nrow(events_A) * nrow(events_B)
    if (n_pairs > 0) {
      distances_list[[i]] <- expand_grid(
        event_A_idx = seq_len(nrow(events_A)),
        event_B_idx = seq_len(nrow(events_B))
      ) %>%
        mutate(
          gene_id = filtered_trans$gene_id[i],
          cell_type = filtered_trans$cell_type[i],
          genomic_distance = abs(events_A$position[event_A_idx] - events_B$position[event_B_idx]),
          exon_distance = abs(events_A$exon_number[event_A_idx] - events_B$exon_number[event_B_idx]),
          order = if_else(
            events_A$position[event_A_idx] < events_B$position[event_B_idx],
            paste0(event_A_type, "_first"),
            paste0(event_B_type, "_first")
          )
        ) %>%
        select(gene_id, cell_type, genomic_distance, exon_distance, order)
    }
  }

  # Combine all distances
  bind_rows(distances_list) %>%
    filter(!is.na(genomic_distance))
}

# Calculate distances for positive pairs
if (nrow(positive_pairs) > 0) {
  cat("Calculating distances for positively correlated pairs...\n")
  positive_distances <- map_dfr(1:nrow(positive_pairs), function(i) {
    # Extract values from the row
    event_a <- positive_pairs$event_A[i]
    event_b <- positive_pairs$event_B[i]
    or_val <- positive_pairs$odds_ratio[i]
    pval <- positive_pairs$p_value[i]

    cat("  Processing:", event_a, "-", event_b, "\n")

    dists <- calculate_event_distances(event_a, event_b, transitions)
    if (nrow(dists) > 0) {
      dists %>%
        mutate(
          pair = paste(event_a, event_b, sep = "-"),
          correlation_type = "positive",
          odds_ratio = or_val,
          p_value = pval
        )
    } else {
      NULL
    }
  })
  
  if (!is.null(positive_distances) && nrow(positive_distances) > 0) {
    # Summary statistics
    distance_summary <- positive_distances %>%
      group_by(pair, correlation_type) %>%
      summarize(
        n_observations = n(),
        mean_genomic_dist = mean(genomic_distance),
        median_genomic_dist = median(genomic_distance),
        mean_exon_dist = mean(exon_distance),
        median_exon_dist = median(exon_distance),
        .groups = "drop"
      )
    
    cat("\nDistance summary for positively correlated pairs:\n")
    print(distance_summary)
    
    # Save
    write_tsv(positive_distances,
              file.path(output_dir, "q3_positive_pair_distances.tsv"))
    write_tsv(distance_summary,
              file.path(output_dir, "q3_distance_summary.tsv"))
    
    cat("\nSaved Q3 distance results\n")
  } else {
    cat("No distance data available for positively correlated pairs\n")
  }
} else {
  cat("No positively correlated pairs found\n")
}

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║           DETAILED EVENT ANALYSIS COMPLETE (Q1-Q3)            ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")
