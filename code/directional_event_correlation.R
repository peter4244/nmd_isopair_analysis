#!/usr/bin/env Rscript
# Directional Event Correlation Analysis
# Tests co-occurrence of all 12 directional event types

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   DIRECTIONAL EVENT CORRELATION ANALYSIS                       ║\n")
cat("║   Testing all 66 pairwise correlations (12 choose 2)           ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

output_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions"

# Load directional event data
cat("Loading directional event data (v2.2)...\n")
transitions <- readRDS(file.path(output_dir, "detailed_event_vectors_v2.2_directional.rds"))
cat("Loaded:", nrow(transitions), "transitions\n\n")

# ============================================================================
# Directional Event Co-occurrence Testing
# ============================================================================

cat("═══ Testing Directional Event Co-occurrence ═══\n\n")

# All 12 directional event types
event_types <- c(
  "Alt_TSS_gain_A", "Alt_TSS_gain_B",
  "Alt_TES_gain_A", "Alt_TES_gain_B",
  "A5SS_gain_A", "A5SS_gain_B",
  "A3SS_gain_A", "A3SS_gain_B",
  "SE_gain_A", "SE_gain_B",
  "IR_gain_A", "IR_gain_B"
)

# Generate all pairwise combinations (66 total)
event_pairs <- combn(event_types, 2, simplify = FALSE)

cat("Testing", length(event_pairs), "pairwise combinations...\n\n")

# Perform Fisher's exact test for each pair
cooccurrence_results <- map_dfr(event_pairs, function(pair) {
  event_A <- pair[1]
  event_B <- pair[2]

  has_A <- transitions[[paste0("has_", event_A)]]
  has_B <- transitions[[paste0("has_", event_B)]]

  n_both <- sum(has_A & has_B)
  n_A_only <- sum(has_A & !has_B)
  n_B_only <- sum(!has_A & has_B)
  n_neither <- sum(!has_A & !has_B)

  contingency <- matrix(c(n_both, n_A_only, n_B_only, n_neither), nrow = 2)
  test_result <- fisher.test(contingency)

  # Determine relationship type
  correlation_type <- case_when(
    test_result$p.value >= 0.05 ~ "non-correlated",
    test_result$estimate > 1 ~ "positively correlated",
    test_result$estimate < 1 ~ "negatively correlated",
    TRUE ~ "non-correlated"
  )

  # Parse event types to identify patterns
  event_A_parts <- str_split(event_A, "_")[[1]]
  event_B_parts <- str_split(event_B, "_")[[1]]

  event_A_type <- paste(event_A_parts[1:(length(event_A_parts)-2)], collapse="_")
  event_A_iso <- tail(event_A_parts, 1)
  event_B_type <- paste(event_B_parts[1:(length(event_B_parts)-2)], collapse="_")
  event_B_iso <- tail(event_B_parts, 1)

  # Categorize relationship
  same_isoform <- event_A_iso == event_B_iso
  same_event_type <- event_A_type == event_B_type

  relationship_pattern <- case_when(
    same_isoform & same_event_type ~ "same_event_same_iso",
    same_isoform & !same_event_type ~ "coordinated_gain",  # Both gain in same isoform
    !same_isoform & same_event_type ~ "antagonistic",      # Same event, opposite isoforms
    !same_isoform & !same_event_type ~ "cross_isoform"     # Different events, different isoforms
  )

  tibble(
    event_A = event_A,
    event_B = event_B,
    event_A_type = event_A_type,
    event_A_iso = event_A_iso,
    event_B_type = event_B_type,
    event_B_iso = event_B_iso,
    n_both = n_both,
    n_A_only = n_A_only,
    n_B_only = n_B_only,
    n_neither = n_neither,
    odds_ratio = test_result$estimate,
    p_value = test_result$p.value,
    significant = test_result$p.value < 0.05,
    correlation_type = correlation_type,
    relationship_pattern = relationship_pattern
  )
})

# ============================================================================
# Categorize Results
# ============================================================================

cat("Categorizing correlation patterns...\n\n")

# Separate by correlation type
positive_corr <- filter(cooccurrence_results, correlation_type == "positively correlated")
negative_corr <- filter(cooccurrence_results, correlation_type == "negatively correlated")
noncorr <- filter(cooccurrence_results, correlation_type == "non-correlated")

cat("Summary of correlations:\n")
cat("  Positively correlated pairs:", nrow(positive_corr), "\n")
cat("  Negatively correlated pairs:", nrow(negative_corr), "\n")
cat("  Non-correlated pairs:", nrow(noncorr), "\n\n")

# Breakdown by relationship pattern
pattern_summary <- cooccurrence_results %>%
  filter(significant) %>%
  count(relationship_pattern, correlation_type) %>%
  arrange(desc(n))

cat("Correlation patterns:\n")
print(pattern_summary)
cat("\n")

# ============================================================================
# Highlight Key Findings
# ============================================================================

cat("═══ Key Findings ═══\n\n")

# Top coordinated gains (same isoform, positive correlation)
cat("Top COORDINATED GAINS (both events gain in same isoform):\n")
coordinated <- positive_corr %>%
  filter(relationship_pattern == "coordinated_gain") %>%
  arrange(desc(odds_ratio)) %>%
  select(event_A, event_B, n_both, odds_ratio, p_value)
print(head(coordinated, 10))
cat("\n")

# Top antagonistic (same event type, opposite isoforms, negative correlation)
cat("Top ANTAGONISTIC patterns (same event, opposite isoforms):\n")
antagonistic <- negative_corr %>%
  filter(relationship_pattern == "antagonistic") %>%
  arrange(odds_ratio) %>%
  select(event_A, event_B, n_both, odds_ratio, p_value)
print(head(antagonistic, 10))
cat("\n")

# Interesting cross-isoform positive correlations
cat("Interesting CROSS-ISOFORM positive correlations:\n")
cross_iso_pos <- positive_corr %>%
  filter(relationship_pattern == "cross_isoform") %>%
  arrange(desc(odds_ratio)) %>%
  select(event_A, event_B, n_both, odds_ratio, p_value)
print(head(cross_iso_pos, 10))
cat("\n")

# ============================================================================
# Coordinated vs Antagonistic Summary
# ============================================================================

cat("═══ Coordinated vs Antagonistic Events ═══\n\n")

# Summary of coordinated gains by event type combination
coordinated_summary <- positive_corr %>%
  filter(relationship_pattern == "coordinated_gain") %>%
  mutate(event_pair = paste(pmin(event_A_type, event_B_type),
                            pmax(event_A_type, event_B_type), sep = " + ")) %>%
  group_by(event_pair) %>%
  summarize(
    n_pairs = n(),
    mean_OR = mean(odds_ratio),
    max_OR = max(odds_ratio),
    mean_n_both = mean(n_both),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_OR))

cat("Coordinated event pairs (same isoform gains both):\n")
print(coordinated_summary)
cat("\n")

# Summary of antagonistic patterns
antagonistic_summary <- negative_corr %>%
  filter(relationship_pattern == "antagonistic") %>%
  group_by(event_A_type) %>%
  summarize(
    n_antagonistic = n(),
    mean_OR = mean(odds_ratio),
    min_OR = min(odds_ratio),
    .groups = "drop"
  ) %>%
  arrange(mean_OR)

cat("Antagonistic patterns by event type:\n")
print(antagonistic_summary)
cat("\n")

# ============================================================================
# Save Results
# ============================================================================

write_tsv(cooccurrence_results,
          file.path(output_dir, "directional_cooccurrence_all.tsv"))
write_tsv(positive_corr,
          file.path(output_dir, "directional_cooccurrence_positive.tsv"))
write_tsv(negative_corr,
          file.path(output_dir, "directional_cooccurrence_negative.tsv"))
write_tsv(coordinated_summary,
          file.path(output_dir, "coordinated_gains_summary.tsv"))
write_tsv(antagonistic_summary,
          file.path(output_dir, "antagonistic_patterns_summary.tsv"))

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║     DIRECTIONAL CORRELATION ANALYSIS COMPLETE                  ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("Saved results:\n")
cat("  - directional_cooccurrence_all.tsv (all 66 pairs)\n")
cat("  - directional_cooccurrence_positive.tsv\n")
cat("  - directional_cooccurrence_negative.tsv\n")
cat("  - coordinated_gains_summary.tsv\n")
cat("  - antagonistic_patterns_summary.tsv\n")
cat("\n")
