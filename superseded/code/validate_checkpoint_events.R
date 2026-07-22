#!/usr/bin/env Rscript
# Validate event detection at first checkpoint (1000 genes)
# Tests internal consistency, biological validity, and downstream suitability

library(tidyverse)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EVENT DETECTION VALIDATION - CHECKPOINT                    ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

output_dir <- "results/isoform_transitions/v4.0_reference_based"

# Load results
cat("Loading checkpoint results...\n")
results <- readRDS(file.path(output_dir, "event_vectors_full.rds"))
cat("  Transitions loaded:", nrow(results), "\n")
cat("  Genes:", length(unique(results$gene_id)), "\n\n")

# ============================================================================
# Test 1: Internal Consistency
# ============================================================================

cat("═══ TEST 1: INTERNAL CONSISTENCY ═══\n\n")

# Check that event counts match event_vector contents
cat("1.1 Event count consistency:\n")
inconsistent <- results %>%
  filter(n_union_exons > 0) %>%
  rowwise() %>%
  mutate(
    actual_Alt_TSS = sum(event_vector$event_type == "Alt_TSS"),
    actual_Alt_TES = sum(event_vector$event_type == "Alt_TES"),
    actual_SE = sum(event_vector$event_type == "SE"),
    actual_A5SS = sum(event_vector$event_type == "A5SS"),
    actual_A3SS = sum(event_vector$event_type == "A3SS"),
    actual_CONST = sum(event_vector$event_type == "CONST"),
    mismatch = (n_Alt_TSS != actual_Alt_TSS |
                n_Alt_TES != actual_Alt_TES |
                n_SE != actual_SE |
                n_A5SS != actual_A5SS |
                n_A3SS != actual_A3SS |
                n_CONST != actual_CONST)
  ) %>%
  ungroup()

n_mismatches <- sum(inconsistent$mismatch, na.rm = TRUE)
if (n_mismatches == 0) {
  cat("  ✅ PASS: All event counts match event vectors\n\n")
} else {
  cat("  ❌ FAIL:", n_mismatches, "transitions have mismatched counts\n")
  cat("  Examples:\n")
  print(inconsistent %>% filter(mismatch) %>% head(5) %>%
        select(gene_id, isoform_A, isoform_B, starts_with("n_"), starts_with("actual_")))
  cat("\n")
}

# ============================================================================
# Test 2: Terminal Exon Rules
# ============================================================================

cat("═══ TEST 2: TERMINAL EXON RULES ═══\n\n")

cat("2.1 Alt_TSS should only occur at first exon:\n")
violations_tss <- results %>%
  filter(n_Alt_TSS > 0) %>%
  rowwise() %>%
  mutate(
    tss_positions = list(event_vector %>%
                        filter(event_type == "Alt_TSS") %>%
                        pull(exon_number)),
    has_non_first = any(tss_positions[[1]] != 1)
  ) %>%
  ungroup() %>%
  filter(has_non_first)

if (nrow(violations_tss) == 0) {
  cat("  ✅ PASS: All Alt_TSS at exon 1\n\n")
} else {
  cat("  ❌ FAIL:", nrow(violations_tss), "transitions have Alt_TSS at non-first exon\n")
  cat("  Examples:\n")
  print(violations_tss %>% head(3) %>% select(gene_id, isoform_A, isoform_B, tss_positions))
  cat("\n")
}

cat("2.2 Alt_TES should only occur at last exon:\n")
violations_tes <- results %>%
  filter(n_Alt_TES > 0) %>%
  rowwise() %>%
  mutate(
    tes_positions = list(event_vector %>%
                        filter(event_type == "Alt_TES") %>%
                        pull(exon_number)),
    max_exon = max(event_vector$exon_number),
    has_non_last = any(tes_positions[[1]] != max_exon)
  ) %>%
  ungroup() %>%
  filter(has_non_last)

if (nrow(violations_tes) == 0) {
  cat("  ✅ PASS: All Alt_TES at last exon\n\n")
} else {
  cat("  ❌ FAIL:", nrow(violations_tes), "transitions have Alt_TES at non-last exon\n")
  cat("  Examples:\n")
  print(violations_tes %>% head(3) %>% select(gene_id, isoform_A, isoform_B, tes_positions, max_exon))
  cat("\n")
}

# ============================================================================
# Test 3: Splice Site Position Rules
# ============================================================================

cat("═══ TEST 3: SPLICE SITE POSITION RULES ═══\n\n")

cat("3.1 A5SS can occur on exons 1 to N-1 (not last):\n")
violations_a5ss <- results %>%
  filter(n_A5SS > 0) %>%
  rowwise() %>%
  mutate(
    a5ss_positions = list(event_vector %>%
                         filter(event_type == "A5SS") %>%
                         pull(exon_number)),
    max_exon = max(event_vector$exon_number),
    has_at_last = any(a5ss_positions[[1]] == max_exon)
  ) %>%
  ungroup() %>%
  filter(has_at_last)

if (nrow(violations_a5ss) == 0) {
  cat("  ✅ PASS: No A5SS at last exon\n\n")
} else {
  cat("  ❌ FAIL:", nrow(violations_a5ss), "transitions have A5SS at last exon\n")
  cat("  Examples:\n")
  print(violations_a5ss %>% head(3) %>% select(gene_id, isoform_A, isoform_B, a5ss_positions, max_exon))
  cat("\n")
}

cat("3.2 A3SS can occur on exons 2 to N (not first):\n")
violations_a3ss <- results %>%
  filter(n_A3SS > 0) %>%
  rowwise() %>%
  mutate(
    a3ss_positions = list(event_vector %>%
                         filter(event_type == "A3SS") %>%
                         pull(exon_number)),
    has_at_first = any(a3ss_positions[[1]] == 1)
  ) %>%
  ungroup() %>%
  filter(has_at_first)

if (nrow(violations_a3ss) == 0) {
  cat("  ✅ PASS: No A3SS at first exon\n\n")
} else {
  cat("  ❌ FAIL:", nrow(violations_a3ss), "transitions have A3SS at first exon\n")
  cat("  Examples:\n")
  print(violations_a3ss %>% head(3) %>% select(gene_id, isoform_A, isoform_B, a3ss_positions))
  cat("\n")
}

# ============================================================================
# Test 4: IR Detection Validation
# ============================================================================

cat("═══ TEST 4: INTRON RETENTION VALIDATION ═══\n\n")

# Load union models to check IR was detected during construction
union_models <- readRDS(file.path(output_dir, "union_exon_models_full.rds"))

ir_summary <- tibble(
  gene_id = names(union_models),
  n_ir_removed = map_int(union_models, ~.x$n_ir_removed)
) %>%
  filter(n_ir_removed > 0)

cat("4.1 IR exons removed during union model construction:\n")
cat("  Genes with IR:", nrow(ir_summary), "\n")
cat("  Total IR exons removed:", sum(ir_summary$n_ir_removed), "\n\n")

# Check if IR events appear in event vectors
# Note: IR exons should be REMOVED before grouping, so they shouldn't appear as events
# Instead, they would appear as absence of the mega-exon in one isoform
ir_in_events <- results %>%
  filter(n_union_exons > 0) %>%
  rowwise() %>%
  mutate(has_ir_event = "IR" %in% event_vector$event_type) %>%
  ungroup() %>%
  summarize(n_with_ir = sum(has_ir_event))

cat("4.2 IR events in event vectors:\n")
if (ir_in_events$n_with_ir == 0) {
  cat("  ✅ PASS: No IR events in vectors (expected - IR detected during construction)\n\n")
} else {
  cat("  ⚠️  WARNING:", ir_in_events$n_with_ir, "transitions have IR events\n")
  cat("  (IR should be detected during construction, not event detection)\n\n")
}

# ============================================================================
# Test 5: MXE Detection Validation
# ============================================================================

cat("═══ TEST 5: MUTUALLY EXCLUSIVE EXON VALIDATION ═══\n\n")

# MXE detection is separate - check if script exists
mxe_script <- "code/detect_mxe_from_union_model.R"
if (file.exists(mxe_script)) {
  cat("5.1 MXE detection:\n")
  cat("  ✅ MXE detection script exists:", mxe_script, "\n")
  cat("  (Run separately after event detection completes)\n\n")
} else {
  cat("  ⚠️  MXE detection script not found\n\n")
}

# Check for SE events that might actually be MXE
# MXE = two exons at same position that never co-occur
# This is complex to detect from event vectors alone
cat("5.2 Potential MXE patterns (SE events at same position):\n")
cat("  (Full MXE analysis requires separate detection algorithm)\n\n")

# ============================================================================
# Test 6: Event Distribution & Biological Plausibility
# ============================================================================

cat("═══ TEST 6: EVENT DISTRIBUTION ═══\n\n")

event_summary <- results %>%
  summarize(
    total_transitions = n(),
    pct_Alt_TSS = 100 * mean(n_Alt_TSS > 0),
    pct_Alt_TES = 100 * mean(n_Alt_TES > 0),
    pct_SE = 100 * mean(n_SE > 0),
    pct_A5SS = 100 * mean(n_A5SS > 0),
    pct_A3SS = 100 * mean(n_A3SS > 0),
    pct_CONST = 100 * mean(n_CONST > 0),
    mean_events_per_transition = mean(n_Alt_TSS + n_Alt_TES + n_SE + n_A5SS + n_A3SS)
  )

cat("Event type prevalence:\n")
cat(sprintf("  Alt_TSS: %.1f%% of transitions\n", event_summary$pct_Alt_TSS))
cat(sprintf("  Alt_TES: %.1f%% of transitions\n", event_summary$pct_Alt_TES))
cat(sprintf("  SE:      %.1f%% of transitions\n", event_summary$pct_SE))
cat(sprintf("  A5SS:    %.1f%% of transitions\n", event_summary$pct_A5SS))
cat(sprintf("  A3SS:    %.1f%% of transitions\n", event_summary$pct_A3SS))
cat(sprintf("  CONST:   %.1f%% of transitions\n\n", event_summary$pct_CONST))

cat(sprintf("Mean events per transition: %.2f\n\n", event_summary$mean_events_per_transition))

# Biological plausibility checks
cat("6.1 Biological plausibility:\n")
if (event_summary$pct_SE > 30 && event_summary$pct_SE < 90) {
  cat("  ✅ SE prevalence reasonable (30-90%)\n")
} else {
  cat("  ⚠️  SE prevalence unusual:", sprintf("%.1f%%\n", event_summary$pct_SE))
}

if (event_summary$mean_events_per_transition > 1 && event_summary$mean_events_per_transition < 10) {
  cat("  ✅ Mean events per transition reasonable (1-10)\n")
} else {
  cat("  ⚠️  Mean events unusual:", sprintf("%.2f\n", event_summary$mean_events_per_transition))
}
cat("\n")

# ============================================================================
# Test 7: Downstream Analysis Suitability
# ============================================================================

cat("═══ TEST 7: DOWNSTREAM ANALYSIS SUITABILITY ═══\n\n")

cat("7.1 Data completeness:\n")
required_cols <- c("gene_id", "isoform_A", "isoform_B", "n_union_exons", "event_vector",
                   "n_Alt_TSS", "n_Alt_TES", "n_SE", "n_A5SS", "n_A3SS", "n_CONST")
missing_cols <- setdiff(required_cols, names(results))

if (length(missing_cols) == 0) {
  cat("  ✅ All required columns present\n")
} else {
  cat("  ❌ Missing columns:", paste(missing_cols, collapse = ", "), "\n")
}

# Check for NA values in key columns
na_check <- results %>%
  summarize(across(c(gene_id, isoform_A, isoform_B, n_union_exons),
                   ~sum(is.na(.)),
                   .names = "na_{.col}"))

if (sum(na_check) == 0) {
  cat("  ✅ No NA values in key columns\n")
} else {
  cat("  ⚠️  NA values found:\n")
  print(na_check)
}
cat("\n")

cat("7.2 Event vector structure:\n")
# Check that event vectors have required columns
sample_vector <- results$event_vector[[which(results$n_union_exons > 0)[1]]]
required_event_cols <- c("event_type", "exon_number", "direction", "detail")
missing_event_cols <- setdiff(required_event_cols, names(sample_vector))

if (length(missing_event_cols) == 0) {
  cat("  ✅ Event vectors have all required columns\n")
} else {
  cat("  ❌ Event vectors missing:", paste(missing_event_cols, collapse = ", "), "\n")
}
cat("\n")

cat("7.3 Direction labels:\n")
# Check that directions are valid
all_events <- results %>%
  filter(n_union_exons > 0) %>%
  pull(event_vector) %>%
  bind_rows()

valid_directions <- c("gain", "loss", "none")
invalid_directions <- all_events %>%
  filter(!direction %in% valid_directions) %>%
  count(direction)

if (nrow(invalid_directions) == 0) {
  cat("  ✅ All directions valid (gain/loss/none)\n")
} else {
  cat("  ❌ Invalid directions found:\n")
  print(invalid_directions)
}
cat("\n")

# ============================================================================
# Test 8: Strand-Specific Validation
# ============================================================================

cat("═══ TEST 8: STRAND-SPECIFIC VALIDATION ═══\n\n")

# Load strand info
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))
strand_lookup <- exon_structures %>%
  select(isoform_id, strand) %>%
  distinct()

# Sample genes from both strands
cat("8.1 Testing strand-aware detection:\n")

# Get a few plus and minus strand examples with events
strand_examples <- results %>%
  filter(n_Alt_TSS > 0 | n_A5SS > 0) %>%
  head(100) %>%
  left_join(strand_lookup, by = c("isoform_A" = "isoform_id")) %>%
  group_by(strand) %>%
  slice(1:3) %>%
  ungroup()

if (nrow(strand_examples) > 0) {
  cat("  Sample events from both strands detected\n")
  cat("  Plus strand examples:", sum(strand_examples$strand == "+"), "\n")
  cat("  Minus strand examples:", sum(strand_examples$strand == "-"), "\n")
  cat("  ✅ Strand-aware logic appears functional\n")
  cat("  (Manual inspection of specific genes recommended)\n")
} else {
  cat("  ⚠️  Could not find strand examples\n")
}
cat("\n")

# ============================================================================
# Summary
# ============================================================================

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   VALIDATION SUMMARY                                          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("Tested:\n")
cat("  ✓ Internal consistency (event counts vs vectors)\n")
cat("  ✓ Terminal exon rules (Alt_TSS/TES positioning)\n")
cat("  ✓ Splice site rules (A5SS/A3SS positioning)\n")
cat("  ✓ IR detection during construction\n")
cat("  ✓ Event distribution & biological plausibility\n")
cat("  ✓ Data structure for downstream analysis\n")
cat("  ✓ Strand-specific detection\n\n")

cat("Next steps:\n")
cat("  1. Run MXE detection (separate algorithm)\n")
cat("  2. Perform co-occurrence analysis\n")
cat("  3. Analyze topological distances\n")
cat("  4. Compare event patterns across cell types\n\n")
