#!/usr/bin/env Rscript
# Reference-Based Event Correlation Analysis v3.0
# Tests co-occurrence of gain/loss events relative to dominant isoform

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   REFERENCE-BASED EVENT CORRELATION ANALYSIS v3.0             ║\n")
cat("║   Testing coordinated gains/losses vs dominant isoform        ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

output_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v3.0_reference_based"

# Load reference-based event data (filtered to SR-DGE genes)
cat("Loading reference-based event data (filtered)...\n")
transitions <- readRDS(file.path(output_dir, "reference_event_vectors_v3.0_filtered.rds"))
cat("Loaded:", nrow(transitions), "transitions\n\n")

# ============================================================================
# Event Co-occurrence Testing
# ============================================================================

cat("═══ Testing Event Co-occurrence ═══\n\n")

# All 12 event types (6 events × 2 directions: gain/loss)
event_types <- c(
  "Alt_TSS_gain", "Alt_TSS_loss",
  "Alt_TES_gain", "Alt_TES_loss",
  "A5SS_gain", "A5SS_loss",
  "A3SS_gain", "A3SS_loss",
  "SE_gain", "SE_loss",
  "IR_gain", "IR_loss"
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

  # Determine correlation type
  correlation_type <- case_when(
    test_result$p.value >= 0.05 ~ "non-correlated",
    test_result$estimate > 1 ~ "positively correlated",
    test_result$estimate < 1 ~ "negatively correlated",
    TRUE ~ "non-correlated"
  )

  # Parse event types to identify patterns
  event_A_parts <- str_split(event_A, "_")[[1]]
  event_B_parts <- str_split(event_B, "_")[[1]]

  # Extract event base and direction
  event_A_base <- paste(event_A_parts[1:(length(event_A_parts)-1)], collapse="_")
  event_A_dir <- tail(event_A_parts, 1)
  event_B_base <- paste(event_B_parts[1:(length(event_B_parts)-1)], collapse="_")
  event_B_dir <- tail(event_B_parts, 1)

  # Categorize relationship
  same_direction <- event_A_dir == event_B_dir
  same_event_type <- event_A_base == event_B_base

  relationship_pattern <- case_when(
    same_event_type & same_direction ~ "same_event_same_direction",
    same_event_type & !same_direction ~ "antagonistic",  # Same event, opposite directions
    !same_event_type & same_direction ~ "coordinated",   # Different events, same direction
    !same_event_type & !same_direction ~ "mixed"         # Different events, different directions
  )

  tibble(
    event_A = event_A,
    event_B = event_B,
    event_A_base = event_A_base,
    event_A_dir = event_A_dir,
    event_B_base = event_B_base,
    event_B_dir = event_B_dir,
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
# Key Findings
# ============================================================================

cat("═══ Key Findings ═══\n\n")

# Top coordinated GAINS (both events gain)
cat("Top COORDINATED GAINS (both events add sequence vs reference):\n")
coordinated_gains <- positive_corr %>%
  filter(relationship_pattern == "coordinated" &
         event_A_dir == "gain" & event_B_dir == "gain") %>%
  arrange(desc(odds_ratio)) %>%
  select(event_A, event_B, n_both, odds_ratio, p_value)
print(head(coordinated_gains, 10))
cat("\n")

# Top coordinated LOSSES (both events lose)
cat("Top COORDINATED LOSSES (both events remove sequence vs reference):\n")
coordinated_losses <- positive_corr %>%
  filter(relationship_pattern == "coordinated" &
         event_A_dir == "loss" & event_B_dir == "loss") %>%
  arrange(desc(odds_ratio)) %>%
  select(event_A, event_B, n_both, odds_ratio, p_value)
print(head(coordinated_losses, 10))
cat("\n")

# Antagonistic patterns (same event, opposite directions)
cat("Top ANTAGONISTIC patterns (same event type, gain vs loss):\n")
antagonistic <- negative_corr %>%
  filter(relationship_pattern == "antagonistic") %>%
  arrange(odds_ratio) %>%
  select(event_A, event_B, n_both, odds_ratio, p_value)
print(head(antagonistic, 10))
cat("\n")

# Mixed patterns (different events, different directions)
cat("Interesting MIXED patterns (gain of one, loss of another):\n")
mixed_pos <- positive_corr %>%
  filter(relationship_pattern == "mixed") %>%
  arrange(desc(odds_ratio)) %>%
  select(event_A, event_B, n_both, odds_ratio, p_value)
print(head(mixed_pos, 10))
cat("\n")

# ============================================================================
# NMD Signature Analysis
# ============================================================================

cat("═══ NMD Signature Candidates ═══\n\n")

# Identify isoforms with coordinated gains (SE + IR = classic NMD signature)
cat("Testing SE_gain + IR_gain correlation (classic NMD signature):\n")
se_ir_gain <- cooccurrence_results %>%
  filter(event_A == "SE_gain" & event_B == "IR_gain")
print(se_ir_gain)
cat("\n")

# Count isoforms with multiple gain events
cat("Isoforms with multiple sequence gains:\n")
multi_gain <- transitions %>%
  mutate(
    n_gains = has_Alt_TSS_gain + has_Alt_TES_gain + has_A5SS_gain +
              has_A3SS_gain + has_SE_gain + has_IR_gain,
    n_losses = has_Alt_TSS_loss + has_Alt_TES_loss + has_A5SS_loss +
               has_A3SS_loss + has_SE_loss + has_IR_loss
  ) %>%
  summarize(
    mean_gains = mean(n_gains),
    mean_losses = mean(n_losses),
    pct_net_gainers = 100 * mean(n_gains > n_losses),
    pct_net_losers = 100 * mean(n_losses > n_gains),
    pct_balanced = 100 * mean(n_gains == n_losses)
  )

print(multi_gain)
cat("\n")

# ============================================================================
# Summary Statistics by Pattern
# ============================================================================

cat("═══ Summary by Relationship Pattern ═══\n\n")

# Coordinated gains summary
if (nrow(coordinated_gains) > 0) {
  coord_gain_summary <- positive_corr %>%
    filter(relationship_pattern == "coordinated" &
           event_A_dir == "gain" & event_B_dir == "gain") %>%
    mutate(event_pair = paste(pmin(event_A_base, event_B_base),
                              pmax(event_A_base, event_B_base), sep = " + ")) %>%
    group_by(event_pair) %>%
    summarize(
      n_pairs = n(),
      mean_OR = mean(odds_ratio),
      max_OR = max(odds_ratio),
      mean_n_both = mean(n_both),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_OR))

  cat("Coordinated GAIN event pairs:\n")
  print(coord_gain_summary)
  cat("\n")
}

# Coordinated losses summary
if (nrow(coordinated_losses) > 0) {
  coord_loss_summary <- positive_corr %>%
    filter(relationship_pattern == "coordinated" &
           event_A_dir == "loss" & event_B_dir == "loss") %>%
    mutate(event_pair = paste(pmin(event_A_base, event_B_base),
                              pmax(event_A_base, event_B_base), sep = " + ")) %>%
    group_by(event_pair) %>%
    summarize(
      n_pairs = n(),
      mean_OR = mean(odds_ratio),
      max_OR = max(odds_ratio),
      mean_n_both = mean(n_both),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_OR))

  cat("Coordinated LOSS event pairs:\n")
  print(coord_loss_summary)
  cat("\n")
}

# ============================================================================
# Save Results
# ============================================================================

cat("Saving results...\n")

write_tsv(cooccurrence_results,
          file.path(output_dir, "reference_cooccurrence_all.tsv"))
write_tsv(positive_corr,
          file.path(output_dir, "reference_cooccurrence_positive.tsv"))
write_tsv(negative_corr,
          file.path(output_dir, "reference_cooccurrence_negative.tsv"))

if (exists("coord_gain_summary")) {
  write_tsv(coord_gain_summary,
            file.path(output_dir, "reference_coordinated_gains_summary.tsv"))
}

if (exists("coord_loss_summary")) {
  write_tsv(coord_loss_summary,
            file.path(output_dir, "reference_coordinated_losses_summary.tsv"))
}

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   REFERENCE-BASED CORRELATION ANALYSIS COMPLETE               ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("Saved results:\n")
cat("  - reference_cooccurrence_all.tsv (all 66 pairs)\n")
cat("  - reference_cooccurrence_positive.tsv\n")
cat("  - reference_cooccurrence_negative.tsv\n")
cat("  - reference_coordinated_gains_summary.tsv\n")
cat("  - reference_coordinated_losses_summary.tsv\n")
cat("\n")
