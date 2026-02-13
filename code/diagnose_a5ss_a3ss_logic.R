#!/usr/bin/env Rscript
# Diagnostic: A5SS/A3SS Detection Logic
#
# Purpose: Investigate why no A5SS/A3SS events are detected
#
# Key Question: Do exons in the same union exon group ever have
#               different coordinates (same start OR same end, but not both)?

library(tidyverse)

cat("\n")
cat("═══ A5SS/A3SS Detection Diagnostic ═══\n\n")

# Configuration
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v4.0_reference_based"
output_dir <- file.path(base_dir, "major_isoforms")

# Load union models
cat("Loading union models...\n")
union_models <- readRDS(file.path(output_dir, "union_exon_models_major.rds"))
cat("  Genes:", length(union_models), "\n\n")

# ============================================================================
# Analyze Union Exon Groups
# ============================================================================

cat("═══ Analyzing Union Exon Groups ═══\n\n")

coord_analysis <- list()

for (gene_id in names(union_models)[1:100]) {  # First 100 genes
  model <- union_models[[gene_id]]

  for (i in seq_along(model$union_exons)) {
    exon_group <- model$union_exons[[i]]

    # Only analyze groups with 2+ isoforms
    if (nrow(exon_group) < 2) next

    # Get unique coordinates in this group
    unique_starts <- unique(exon_group$start)
    unique_ends <- unique(exon_group$end)

    n_unique_starts <- length(unique_starts)
    n_unique_ends <- length(unique_ends)

    # Classify this group
    same_start <- (n_unique_starts == 1)
    same_end <- (n_unique_ends == 1)

    classification <- case_when(
      same_start && same_end ~ "identical",
      same_start && !same_end ~ "A3SS_candidate",  # same start, different ends
      !same_start && same_end ~ "A5SS_candidate",  # different starts, same end
      !same_start && !same_end ~ "both_differ"
    )

    coord_analysis[[length(coord_analysis) + 1]] <- tibble(
      gene_id = gene_id,
      union_exon_num = i,
      n_isoforms = nrow(exon_group),
      exon_type = exon_group$exon_type[1],
      n_unique_starts = n_unique_starts,
      n_unique_ends = n_unique_ends,
      classification = classification,
      starts = paste(unique_starts, collapse = ","),
      ends = paste(unique_ends, collapse = ",")
    )
  }
}

coord_df <- bind_rows(coord_analysis)

cat("Union exon groups with 2+ isoforms:", nrow(coord_df), "\n\n")

# ============================================================================
# Summary by Classification
# ============================================================================

cat("═══ Classification Summary ═══\n\n")

classification_summary <- coord_df %>%
  count(classification) %>%
  arrange(desc(n))

print(classification_summary)
cat("\n")

# ============================================================================
# Summary by Exon Type
# ============================================================================

cat("═══ By Exon Type ═══\n\n")

type_summary <- coord_df %>%
  count(exon_type, classification) %>%
  pivot_wider(names_from = classification, values_from = n, values_fill = 0)

print(type_summary)
cat("\n")

# ============================================================================
# Examples of Each Type
# ============================================================================

cat("═══ Examples ═══\n\n")

for (class_type in c("identical", "A5SS_candidate", "A3SS_candidate", "both_differ")) {
  examples <- coord_df %>%
    filter(classification == class_type) %>%
    head(3)

  if (nrow(examples) > 0) {
    cat(sprintf("\n%s (%d cases):\n", class_type, sum(coord_df$classification == class_type)))
    for (i in 1:nrow(examples)) {
      ex <- examples[i, ]
      cat(sprintf("  Gene: %s, Union Exon: %d, Type: %s\n",
                  ex$gene_id, ex$union_exon_num, ex$exon_type))
      cat(sprintf("    Starts: %s\n", ex$starts))
      cat(sprintf("    Ends: %s\n", ex$ends))
    }
  }
}

# ============================================================================
# Key Finding
# ============================================================================

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   KEY FINDING                                                  ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

n_identical <- sum(coord_df$classification == "identical")
n_a5ss <- sum(coord_df$classification == "A5SS_candidate")
n_a3ss <- sum(coord_df$classification == "A3SS_candidate")
n_both <- sum(coord_df$classification == "both_differ")

cat("Union exon groups by coordinate pattern:\n")
cat(sprintf("  Identical coordinates:  %5d (%.1f%%)\n",
            n_identical, 100 * n_identical / nrow(coord_df)))
cat(sprintf("  A5SS candidates:        %5d (%.1f%%)\n",
            n_a5ss, 100 * n_a5ss / nrow(coord_df)))
cat(sprintf("  A3SS candidates:        %5d (%.1f%%)\n",
            n_a3ss, 100 * n_a3ss / nrow(coord_df)))
cat(sprintf("  Both differ:            %5d (%.1f%%)\n",
            n_both, 100 * n_both / nrow(coord_df)))
cat("\n")

if (n_identical == nrow(coord_df)) {
  cat("PROBLEM IDENTIFIED:\n")
  cat("  All union exon groups have IDENTICAL coordinates.\n")
  cat("  This means the union exon grouping logic is working CORRECTLY:\n")
  cat("  - Exons are grouped by SHARED boundaries\n")
  cat("  - All exons in a group share at least ONE boundary\n")
  cat("  - The A5SS/A3SS detection code expects to find differences\n")
  cat("  - But by definition of the grouping, this won't happen!\n\n")

  cat("SOLUTION:\n")
  cat("  A5SS/A3SS events are NOT detected by comparing exons within\n")
  cat("  the same union exon group. They should be detected by:\n")
  cat("  1. Comparing ADJACENT union exon boundaries between isoforms\n")
  cat("  2. Looking at junction patterns (what connects to what)\n")
  cat("  3. Finding cases where exons have different boundaries\n")
  cat("     but are in the same RELATIVE position\n\n")
} else {
  cat("RESULT:\n")
  cat("  Found cases where union exon groups have different coordinates.\n")
  cat("  The A5SS/A3SS detection logic may be correct.\n")
  cat("  Need to investigate why events aren't being detected.\n\n")
}

cat("═══ Diagnostic Complete ═══\n\n")
