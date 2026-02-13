#!/usr/bin/env Rscript
# Check for Alternative Splice Sites in Raw Data
#
# Purpose: Before grouping, check if genes actually have A5SS/A3SS
#          by examining consecutive exon pairs (junctions)

library(tidyverse)

cat("\n")
cat("═══ Alternative Splice Site Detection (Junction-Based) ═══\n\n")

# Configuration
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v4.0_reference_based"
output_dir <- file.path(base_dir, "major_isoforms")

# Load data
cat("Loading exon structures...\n")
exon_structures <- readRDS(file.path(output_dir, "exon_structures_expanded.rds"))
cat("  Genes:", length(unique(exon_structures$gene_id)), "\n\n")

# ============================================================================
# Junction-Based A5SS/A3SS Detection
# ============================================================================

cat("═══ Analyzing Junctions ═══\n\n")

# Function to extract junctions from an isoform
get_junctions <- function(exons_df) {
  # Sort by transcript_exon_number
  exons_sorted <- exons_df %>% arrange(transcript_exon_number)

  if (nrow(exons_sorted) < 2) return(tibble())

  junctions <- tibble()

  for (i in 1:(nrow(exons_sorted) - 1)) {
    junction <- tibble(
      donor = exons_sorted$end[i],      # End of current exon
      acceptor = exons_sorted$start[i+1],  # Start of next exon
      upstream_exon = i,
      downstream_exon = i + 1
    )
    junctions <- bind_rows(junctions, junction)
  }

  return(junctions)
}

# Analyze first 20 genes with multiple isoforms
genes_to_check <- exon_structures %>%
  group_by(gene_id) %>%
  summarise(n_isoforms = n_distinct(isoform_id), .groups = "drop") %>%
  filter(n_isoforms >= 2) %>%
  slice_head(n = 20) %>%
  pull(gene_id)

cat("Analyzing", length(genes_to_check), "genes with multiple isoforms\n\n")

a5ss_events <- list()
a3ss_events <- list()

for (gene_id in genes_to_check) {
  gene_exons <- exon_structures %>% filter(gene_id == !!gene_id)
  gene_strand <- gene_exons$strand[1]

  isoforms <- unique(gene_exons$isoform_id)

  # Get junctions for each isoform
  isoform_junctions <- list()
  for (iso in isoforms) {
    iso_exons <- gene_exons %>% filter(isoform_id == iso)
    junctions <- get_junctions(iso_exons)
    if (nrow(junctions) > 0) {
      junctions <- junctions %>% mutate(isoform_id = iso)
      isoform_junctions[[iso]] <- junctions
    }
  }

  if (length(isoform_junctions) < 2) next

  all_junctions <- bind_rows(isoform_junctions)

  # Find alternative splice sites
  # Group junctions by shared donor or acceptor
  for (i in 1:(length(isoform_junctions) - 1)) {
    for (j in (i + 1):length(isoform_junctions)) {
      iso_a <- names(isoform_junctions)[i]
      iso_b <- names(isoform_junctions)[j]

      junc_a <- isoform_junctions[[iso_a]]
      junc_b <- isoform_junctions[[iso_b]]

      # Find junctions with same donor but different acceptor (A3SS on + strand)
      for (row_a in 1:nrow(junc_a)) {
        for (row_b in 1:nrow(junc_b)) {
          donor_a <- junc_a$donor[row_a]
          donor_b <- junc_b$donor[row_b]
          acceptor_a <- junc_a$acceptor[row_a]
          acceptor_b <- junc_b$acceptor[row_b]

          # Same donor, different acceptor
          if (donor_a == donor_b && acceptor_a != acceptor_b) {
            event_type <- if (gene_strand == "+") "A3SS" else "A5SS"
            a3ss_events[[length(a3ss_events) + 1]] <- tibble(
              gene_id = gene_id,
              strand = gene_strand,
              isoform_a = iso_a,
              isoform_b = iso_b,
              event_type = event_type,
              shared_donor = donor_a,
              acceptor_a = acceptor_a,
              acceptor_b = acceptor_b
            )
          }

          # Same acceptor, different donor
          if (acceptor_a == acceptor_b && donor_a != donor_b) {
            event_type <- if (gene_strand == "+") "A5SS" else "A3SS"
            a5ss_events[[length(a5ss_events) + 1]] <- tibble(
              gene_id = gene_id,
              strand = gene_strand,
              isoform_a = iso_a,
              isoform_b = iso_b,
              event_type = event_type,
              shared_acceptor = acceptor_a,
              donor_a = donor_a,
              donor_b = donor_b
            )
          }
        }
      }
    }
  }
}

# ============================================================================
# Results
# ============================================================================

cat("═══ Results ═══\n\n")

a5ss_df <- bind_rows(a5ss_events)
a3ss_df <- bind_rows(a3ss_events)

cat("A5SS events found:", nrow(a5ss_df), "\n")
cat("A3SS events found:", nrow(a3ss_df), "\n\n")

if (nrow(a5ss_df) > 0) {
  cat("First 5 A5SS events:\n")
  print(head(a5ss_df, 5))
  cat("\n")
}

if (nrow(a3ss_df) > 0) {
  cat("First 5 A3SS events:\n")
  print(head(a3ss_df, 5))
  cat("\n")
}

# ============================================================================
# Conclusion
# ============================================================================

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   CONCLUSION                                                   ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

if (nrow(a5ss_df) == 0 && nrow(a3ss_df) == 0) {
  cat("NO A5SS/A3SS events found in first 20 multi-isoform genes.\n\n")
  cat("This could mean:\n")
  cat("  1. These events are rare in the major isoforms\n")
  cat("  2. The union model correctly represents the data\n")
  cat("  3. Alternative splicing is primarily exon skipping, not splice sites\n\n")
} else {
  cat("A5SS/A3SS events EXIST in the data!\n\n")
  cat("This confirms:\n")
  cat("  1. The data contains alternative splice sites\n")
  cat("  2. The union model approach cannot detect them\n")
  cat("  3. Junction-based detection is required\n\n")

  cat("RECOMMENDATION:\n")
  cat("  Replace union-exon-based A5SS/A3SS detection with\n")
  cat("  junction-based approach (as demonstrated above)\n\n")
}

cat("═══ Analysis Complete ═══\n\n")
