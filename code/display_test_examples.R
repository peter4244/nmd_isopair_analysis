#!/usr/bin/env Rscript
# Display examples from test results to verify bug fixes

library(tidyverse)

output_dir <- "results/isoform_transitions/v4.0_reference_based"

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EVENT DETECTION EXAMPLES - v4.0 TEST RESULTS              ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Load test results
results <- readRDS(file.path(output_dir, "test_results_10genes.rds"))

# Load exon structures for coordinate lookup
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))

# Helper function to display event vector
display_event_vector <- function(gene_id, isoform_A, isoform_B, gene_strand, event_vec) {
  cat("─────────────────────────────────────────────────────────────────\n")
  cat(sprintf("Gene: %s (strand: %s)\n", gene_id, gene_strand))
  cat(sprintf("Transition: %s → %s\n", isoform_A, isoform_B))

  # Get exon structures for both isoforms
  exons_A <- exon_structures %>%
    filter(isoform_id == isoform_A) %>%
    select(isoform_id, exon_starts, exon_ends)

  exons_B <- exon_structures %>%
    filter(isoform_id == isoform_B) %>%
    select(isoform_id, exon_starts, exon_ends)

  if (nrow(exons_A) > 0) {
    cat(sprintf("\nIsoform A (%s):\n", isoform_A))
    starts_A <- exons_A$exon_starts[[1]]
    ends_A <- exons_A$exon_ends[[1]]
    for (i in seq_along(starts_A)) {
      cat(sprintf("  Exon %d: %d - %d\n", i, starts_A[i], ends_A[i]))
    }
  }

  if (nrow(exons_B) > 0) {
    cat(sprintf("\nIsoform B (%s):\n", isoform_B))
    starts_B <- exons_B$exon_starts[[1]]
    ends_B <- exons_B$exon_ends[[1]]
    for (i in seq_along(starts_B)) {
      cat(sprintf("  Exon %d: %d - %d\n", i, starts_B[i], ends_B[i]))
    }
  }

  cat("\nDetected Events:\n")
  if (nrow(event_vec) > 0) {
    print(event_vec %>% select(event_type, exon_number, direction, detail), n = 100)
  } else {
    cat("  (no events)\n")
  }
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════
# EXAMPLE 1: Alt_TSS on Plus Strand
# ═══════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════════\n")
cat("EXAMPLE 1: Alt_TSS on Plus Strand\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

alt_tss_plus <- results %>%
  filter(gene_strand == "+", n_Alt_TSS > 0) %>%
  arrange(desc(n_Alt_TSS)) %>%
  head(1)

if (nrow(alt_tss_plus) > 0) {
  display_event_vector(
    alt_tss_plus$gene_id,
    alt_tss_plus$isoform_A,
    alt_tss_plus$isoform_B,
    alt_tss_plus$gene_strand,
    alt_tss_plus$event_vector[[1]]
  )

  cat("✅ VERIFICATION: Alt_TSS detected on plus strand\n")
  cat("   First exons should differ at START coordinate\n\n")
}

# ═══════════════════════════════════════════════════════════════════
# EXAMPLE 2: Alt_TSS on Minus Strand
# ═══════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════════\n")
cat("EXAMPLE 2: Alt_TSS on Minus Strand\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

alt_tss_minus <- results %>%
  filter(gene_strand == "-", n_Alt_TSS > 0) %>%
  arrange(desc(n_Alt_TSS)) %>%
  head(1)

if (nrow(alt_tss_minus) > 0) {
  display_event_vector(
    alt_tss_minus$gene_id,
    alt_tss_minus$isoform_A,
    alt_tss_minus$isoform_B,
    alt_tss_minus$gene_strand,
    alt_tss_minus$event_vector[[1]]
  )

  cat("✅ VERIFICATION: Alt_TSS detected on minus strand\n")
  cat("   First exons should differ at END coordinate (TSS on minus)\n\n")
}

# ═══════════════════════════════════════════════════════════════════
# EXAMPLE 3: Alt_TES on Plus Strand
# ═══════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════════\n")
cat("EXAMPLE 3: Alt_TES on Plus Strand\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

alt_tes_plus <- results %>%
  filter(gene_strand == "+", n_Alt_TES > 0) %>%
  arrange(desc(n_Alt_TES)) %>%
  head(1)

if (nrow(alt_tes_plus) > 0) {
  display_event_vector(
    alt_tes_plus$gene_id,
    alt_tes_plus$isoform_A,
    alt_tes_plus$isoform_B,
    alt_tes_plus$gene_strand,
    alt_tes_plus$event_vector[[1]]
  )

  cat("✅ VERIFICATION: Alt_TES detected on plus strand\n")
  cat("   Last exons should differ at END coordinate\n\n")
}

# ═══════════════════════════════════════════════════════════════════
# EXAMPLE 4: Alt_TES on Minus Strand
# ═══════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════════\n")
cat("EXAMPLE 4: Alt_TES on Minus Strand\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

alt_tes_minus <- results %>%
  filter(gene_strand == "-", n_Alt_TES > 0) %>%
  arrange(desc(n_Alt_TES)) %>%
  head(1)

if (nrow(alt_tes_minus) > 0) {
  display_event_vector(
    alt_tes_minus$gene_id,
    alt_tes_minus$isoform_A,
    alt_tes_minus$isoform_B,
    alt_tes_minus$gene_strand,
    alt_tes_minus$event_vector[[1]]
  )

  cat("✅ VERIFICATION: Alt_TES detected on minus strand\n")
  cat("   Last exons should differ at START coordinate (TES on minus)\n\n")
}

# ═══════════════════════════════════════════════════════════════════
# EXAMPLE 5: SE + A5SS + A3SS
# ═══════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════════\n")
cat("EXAMPLE 5: Complex Transition (SE + A5SS + A3SS)\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

complex <- results %>%
  filter(n_SE > 0, n_A5SS > 0, n_A3SS > 0) %>%
  head(1)

if (nrow(complex) > 0) {
  display_event_vector(
    complex$gene_id,
    complex$isoform_A,
    complex$isoform_B,
    complex$gene_strand,
    complex$event_vector[[1]]
  )

  cat("✅ VERIFICATION: Multiple internal splicing events detected\n\n")
}

# ═══════════════════════════════════════════════════════════════════
# EXAMPLE 6: Verify NO SE at Terminal Exons
# ═══════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════════\n")
cat("EXAMPLE 6: Verify First/Last Exons NOT Called as SE\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

# Check if any SE events have NA exon_number (would indicate terminal exons)
se_events <- results %>%
  select(gene_id, gene_strand, isoform_A, isoform_B, event_vector) %>%
  unnest(event_vector) %>%
  filter(event_type == "SE")

cat(sprintf("Total SE events detected: %d\n", nrow(se_events)))
cat(sprintf("SE events with exon_number: %d\n", sum(!is.na(se_events$exon_number))))
cat(sprintf("SE events with NA exon_number: %d\n", sum(is.na(se_events$exon_number))))

if (sum(is.na(se_events$exon_number)) > 0) {
  cat("\n⚠️  WARNING: Found SE events with NA exon_number\n")
  cat("This should NOT happen - terminal exons should be Alt_TSS/Alt_TES\n")
} else {
  cat("\n✅ VERIFICATION: All SE events have exon numbers (internal exons only)\n")
}

cat("\n")

# ═══════════════════════════════════════════════════════════════════
# SUMMARY STATISTICS
# ═══════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY STATISTICS\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

summary_stats <- results %>%
  group_by(gene_strand) %>%
  summarise(
    n_transitions = n(),
    pct_with_alt_tss = 100 * sum(n_Alt_TSS > 0) / n(),
    pct_with_alt_tes = 100 * sum(n_Alt_TES > 0) / n(),
    pct_with_se = 100 * sum(n_SE > 0) / n(),
    pct_with_a5ss = 100 * sum(n_A5SS > 0) / n(),
    pct_with_a3ss = 100 * sum(n_A3SS > 0) / n(),
    mean_alt_tss = mean(n_Alt_TSS),
    mean_alt_tes = mean(n_Alt_TES),
    .groups = "drop"
  )

print(summary_stats)

cat("\n✓ Examples displayed successfully!\n")
