#!/usr/bin/env Rscript
# Interim Analysis - First 120k Transitions
# QC check to verify directional event classification looks correct

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   INTERIM ANALYSIS - First 200,000 Transitions                ║\n")
cat("║   QC check for v2.2 directional event detection               ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

output_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions"

# Load checkpoint data (first 200k transitions)
cat("Loading checkpoint data (200k transitions)...\n")
checkpoint_file <- file.path(output_dir, "checkpoint_v2.2_200000.rds")

if (!file.exists(checkpoint_file)) {
  cat("ERROR: Checkpoint file not found:", checkpoint_file, "\n")
  cat("Waiting for processing to reach 200,000 transitions...\n")
  quit(status = 1)
}

detailed_events <- readRDS(checkpoint_file)
cat("Loaded:", length(detailed_events), "transitions (raw event data)\n")

# Process into tibble with directional flags (same as main script)
cat("Processing event vectors into directional flags...\n")
event_counts <- map_dfr(seq_along(detailed_events), function(i) {
  ev <- detailed_events[[i]]$event_vector

  # Create binary columns for each directional event type
  tibble(
    row_idx = i,
    # Overall counts
    n_alt_tss = detailed_events[[i]]$n_alt_tss,
    n_alt_tes = detailed_events[[i]]$n_alt_tes,
    n_se = detailed_events[[i]]$n_se,
    n_a5ss = detailed_events[[i]]$n_a5ss,
    n_a3ss = detailed_events[[i]]$n_a3ss,
    n_ir = detailed_events[[i]]$n_ir,
    n_constitutive = detailed_events[[i]]$n_constitutive,

    # Directional event flags (12 types)
    has_Alt_TSS_gain_A = any(ev$event_type == "Alt_TSS" & ev$direction == "gain_A"),
    has_Alt_TSS_gain_B = any(ev$event_type == "Alt_TSS" & ev$direction == "gain_B"),
    has_Alt_TES_gain_A = any(ev$event_type == "Alt_TES" & ev$direction == "gain_A"),
    has_Alt_TES_gain_B = any(ev$event_type == "Alt_TES" & ev$direction == "gain_B"),
    has_A5SS_gain_A = any(ev$event_type == "A5SS" & ev$direction == "gain_A"),
    has_A5SS_gain_B = any(ev$event_type == "A5SS" & ev$direction == "gain_B"),
    has_A3SS_gain_A = any(ev$event_type == "A3SS" & ev$direction == "gain_A"),
    has_A3SS_gain_B = any(ev$event_type == "A3SS" & ev$direction == "gain_B"),
    has_SE_gain_A = any(ev$event_type == "SE" & ev$direction == "gain_A"),
    has_SE_gain_B = any(ev$event_type == "SE" & ev$direction == "gain_B"),
    has_IR_gain_A = any(ev$event_type == "IR" & ev$direction == "gain_A"),
    has_IR_gain_B = any(ev$event_type == "IR" & ev$direction == "gain_B"),

    # Event vector
    event_vector = list(ev)
  )
})

transitions <- event_counts
cat("Processed:", nrow(transitions), "transitions\n\n")

# ============================================================================
# Q1: Event Frequencies
# ============================================================================

cat("═══ Q1: Event Frequencies ═══\n\n")

# Overall counts
event_summary <- transitions %>%
  summarise(
    n_transitions = n(),
    # Calculate percentages FIRST
    pct_Alt_TSS_A = 100 * mean(has_Alt_TSS_gain_A),
    pct_Alt_TSS_B = 100 * mean(has_Alt_TSS_gain_B),
    pct_Alt_TES_A = 100 * mean(has_Alt_TES_gain_A),
    pct_Alt_TES_B = 100 * mean(has_Alt_TES_gain_B),
    pct_A5SS_A = 100 * mean(has_A5SS_gain_A),
    pct_A5SS_B = 100 * mean(has_A5SS_gain_B),
    pct_A3SS_A = 100 * mean(has_A3SS_gain_A),
    pct_A3SS_B = 100 * mean(has_A3SS_gain_B),
    pct_SE_A = 100 * mean(has_SE_gain_A),
    pct_SE_B = 100 * mean(has_SE_gain_B),
    pct_IR_A = 100 * mean(has_IR_gain_A),
    pct_IR_B = 100 * mean(has_IR_gain_B),
    # Then calculate totals
    n_Alt_TSS_A = sum(has_Alt_TSS_gain_A),
    n_Alt_TSS_B = sum(has_Alt_TSS_gain_B),
    n_Alt_TES_A = sum(has_Alt_TES_gain_A),
    n_Alt_TES_B = sum(has_Alt_TES_gain_B),
    n_A5SS_A = sum(has_A5SS_gain_A),
    n_A5SS_B = sum(has_A5SS_gain_B),
    n_A3SS_A = sum(has_A3SS_gain_A),
    n_A3SS_B = sum(has_A3SS_gain_B),
    n_SE_A = sum(has_SE_gain_A),
    n_SE_B = sum(has_SE_gain_B),
    n_IR_A = sum(has_IR_gain_A),
    n_IR_B = sum(has_IR_gain_B)
  )

cat("Directional Event Frequencies:\n")
cat("─────────────────────────────────────────\n")
cat(sprintf("Alt_TSS_gain_A: %6d (%5.1f%%)\n", event_summary$n_Alt_TSS_A, event_summary$pct_Alt_TSS_A))
cat(sprintf("Alt_TSS_gain_B: %6d (%5.1f%%)\n", event_summary$n_Alt_TSS_B, event_summary$pct_Alt_TSS_B))
cat(sprintf("Alt_TES_gain_A: %6d (%5.1f%%)\n", event_summary$n_Alt_TES_A, event_summary$pct_Alt_TES_A))
cat(sprintf("Alt_TES_gain_B: %6d (%5.1f%%)\n", event_summary$n_Alt_TES_B, event_summary$pct_Alt_TES_B))
cat(sprintf("A5SS_gain_A:    %6d (%5.1f%%)\n", event_summary$n_A5SS_A, event_summary$pct_A5SS_A))
cat(sprintf("A5SS_gain_B:    %6d (%5.1f%%)\n", event_summary$n_A5SS_B, event_summary$pct_A5SS_B))
cat(sprintf("A3SS_gain_A:    %6d (%5.1f%%)\n", event_summary$n_A3SS_A, event_summary$pct_A3SS_A))
cat(sprintf("A3SS_gain_B:    %6d (%5.1f%%)\n", event_summary$n_A3SS_B, event_summary$pct_A3SS_B))
cat(sprintf("SE_gain_A:      %6d (%5.1f%%)\n", event_summary$n_SE_A, event_summary$pct_SE_A))
cat(sprintf("SE_gain_B:      %6d (%5.1f%%)\n", event_summary$n_SE_B, event_summary$pct_SE_B))
cat(sprintf("IR_gain_A:      %6d (%5.1f%%)\n", event_summary$n_IR_A, event_summary$pct_IR_A))
cat(sprintf("IR_gain_B:      %6d (%5.1f%%)\n", event_summary$n_IR_B, event_summary$pct_IR_B))
cat("─────────────────────────────────────────\n\n")

# Check balance between A and B
cat("Balance check (should be roughly 50/50):\n")
cat(sprintf("Alt_TSS: A=%5.1f%% vs B=%5.1f%%\n",
    event_summary$pct_Alt_TSS_A, event_summary$pct_Alt_TSS_B))
cat(sprintf("Alt_TES: A=%5.1f%% vs B=%5.1f%%\n",
    event_summary$pct_Alt_TES_A, event_summary$pct_Alt_TES_B))
cat(sprintf("A5SS:    A=%5.1f%% vs B=%5.1f%%\n",
    event_summary$pct_A5SS_A, event_summary$pct_A5SS_B))
cat(sprintf("A3SS:    A=%5.1f%% vs B=%5.1f%%\n",
    event_summary$pct_A3SS_A, event_summary$pct_A3SS_B))
cat(sprintf("SE:      A=%5.1f%% vs B=%5.1f%%\n",
    event_summary$pct_SE_A, event_summary$pct_SE_B))
cat(sprintf("IR:      A=%5.1f%% vs B=%5.1f%%\n\n",
    event_summary$pct_IR_A, event_summary$pct_IR_B))

# ============================================================================
# Q2: Top Event Co-occurrences (Quick Check)
# ============================================================================

cat("═══ Q2: Sample Co-occurrence Tests ═══\n\n")

# Test a few key pairs to verify Fisher's test works
test_pairs <- list(
  c("Alt_TSS_gain_A", "SE_gain_A"),  # Coordinated gain
  c("A5SS_gain_A", "A5SS_gain_B"),   # Antagonistic
  c("SE_gain_A", "SE_gain_B"),       # Antagonistic
  c("Alt_TSS_gain_A", "Alt_TES_gain_A")  # Coordinated gain
)

cat("Testing sample event pairs:\n")
cat("─────────────────────────────────────────────────────────────\n")

for (pair in test_pairs) {
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

  corr_type <- ifelse(test_result$p.value >= 0.05, "non-corr",
                      ifelse(test_result$estimate > 1, "positive", "negative"))

  cat(sprintf("%s + %s:\n", event_A, event_B))
  cat(sprintf("  Both=%d, A_only=%d, B_only=%d, Neither=%d\n",
              n_both, n_A_only, n_B_only, n_neither))
  cat(sprintf("  OR=%.2f, p=%.2e, %s\n\n",
              test_result$estimate, test_result$p.value, corr_type))
}

# ============================================================================
# Quick Sanity Checks
# ============================================================================

cat("═══ Sanity Checks ═══\n\n")

# Check for impossible combinations
impossible <- transitions %>%
  filter(has_Alt_TSS_gain_A & has_Alt_TSS_gain_B)
cat("Transitions with Alt_TSS_gain_A AND Alt_TSS_gain_B (should be 0):", nrow(impossible), "\n")

impossible2 <- transitions %>%
  filter(has_A5SS_gain_A & has_A5SS_gain_B)
cat("Transitions with A5SS_gain_A AND A5SS_gain_B (should be 0):", nrow(impossible2), "\n")

# Check IR detection is working
ir_any <- transitions %>%
  filter(has_IR_gain_A | has_IR_gain_B)
cat("Transitions with any IR event:", nrow(ir_any),
    sprintf("(%.1f%%)\n", 100*nrow(ir_any)/nrow(transitions)))

# Check event vector structure
cat("\nSample event vector (first transition with IR):\n")
if (nrow(ir_any) > 0) {
  sample_events <- ir_any$event_vector[[1]]
  print(head(sample_events, 10))
} else {
  cat("  No IR events detected in first 120k transitions\n")
}

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   INTERIM ANALYSIS COMPLETE                                    ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("Next steps:\n")
cat("  1. Review event frequencies (should be biologically plausible)\n")
cat("  2. Check A/B balance (should be ~50/50)\n")
cat("  3. Verify IR detection is working (>0 events)\n")
cat("  4. Confirm no impossible event combinations\n")
cat("  5. If all looks good, wait for full analysis to complete\n")
cat("\n")
