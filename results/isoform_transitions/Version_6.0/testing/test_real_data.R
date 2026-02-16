#!/usr/bin/env Rscript
library(tidyverse)
library(ggplot2)
library(patchwork)

source("../scripts/event_detection_functions_v2.R")
source("visualization_functions.R")

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   Real Data Event Detection Testing                          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Load sampled data
cat("Loading sampled data...\n")
gtf_df <- parse_gtf("real_data/sampled_isoforms.gtf")
pairs <- read_tsv("real_data/sampled_pairs.tsv", show_col_types = FALSE)

cat(sprintf("  Loaded: %d exons, %d gene pairs\n\n", nrow(gtf_df), nrow(pairs)))

# Build isoform structures
cat("Building isoform structures...\n")
isoform_structures <- gtf_df %>%
  group_by(gene_id, transcript_id) %>%
  arrange(exon_number) %>%
  summarise(
    seqnames = first(seqnames),
    strand = first(strand),
    exons = list(tibble(
      exon_number = exon_number,
      exon_start = start,
      exon_end = end
    )),
    .groups = "drop"
  )

cat(sprintf("  Built structures for %d transcripts\n\n", nrow(isoform_structures)))

# Run event detection
cat("Running event detection...\n")
results <- list()
for (i in seq_len(nrow(pairs))) {
  pair <- pairs[i, ]

  cat(sprintf("[%d/%d] %s\n", i, nrow(pairs), pair$gene_id))
  cat(sprintf("        %s vs %s\n", pair$isoform_A, pair$isoform_B))

  # Get isoform structures
  struct_a <- isoform_structures %>%
    filter(transcript_id == pair$isoform_A) %>%
    pull(exons) %>%
    .[[1]]

  struct_b <- isoform_structures %>%
    filter(transcript_id == pair$isoform_B) %>%
    pull(exons) %>%
    .[[1]]

  strand <- isoform_structures %>%
    filter(transcript_id == pair$isoform_A) %>%
    pull(strand)

  if (nrow(struct_a) == 0 || nrow(struct_b) == 0) {
    cat("        ⚠ Missing data\n\n")
    next
  }

  # Detect events
  events <- detect_splicing_events_v2(struct_a, struct_b, strand)
  detected <- summarize_events_v2(events)

  results[[i]] <- tibble(
    gene_id = pair$gene_id,
    isoform_A = pair$isoform_A,
    isoform_B = pair$isoform_B,
    strand = strand,
    n_exons_A = nrow(struct_a),
    n_exons_B = nrow(struct_b),
    detected_events = detected
  )

  cat(sprintf("        Events: %s\n", detected))
  cat(sprintf("        Exons: A=%d, B=%d\n\n", nrow(struct_a), nrow(struct_b)))
}

results_df <- bind_rows(results)

# Summary statistics
cat("═══════════════════════════════════════════════════════════════════\n")
cat("DETECTION SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat(sprintf("Total pairs analyzed: %d\n\n", nrow(results_df)))

# Event type counts
event_counts <- results_df %>%
  separate_rows(detected_events, sep = ",") %>%
  count(detected_events, name = "count") %>%
  arrange(desc(count))

cat("Event type frequencies:\n")
print(event_counts, n = Inf)

cat("\n")

# Exon count distribution
cat(sprintf("Exon counts: A median=%d (range %d-%d), B median=%d (range %d-%d)\n",
            median(results_df$n_exons_A),
            min(results_df$n_exons_A),
            max(results_df$n_exons_A),
            median(results_df$n_exons_B),
            min(results_df$n_exons_B),
            max(results_df$n_exons_B)))

# Create visualizations
cat("\nCreating visualizations...\n")

# Create a mock expected table for visualization (real data has no expected labels)
expected_mock <- results_df %>%
  mutate(
    isoform_A = isoform_A,
    isoform_B = isoform_B,
    event_type = detected_events,  # Use detected as "expected" for display
    description = sprintf("%s events detected", detected_events),
    expected_union_exons = NA_integer_,
    event_coordinates = "",
    notes = sprintf("Strand: %s | Exons: A=%d, B=%d", strand, n_exons_A, n_exons_B)
  ) %>%
  select(gene_id, isoform_A, isoform_B, event_type, description,
         expected_union_exons, event_coordinates, notes)

# Create plots
plots <- list()
for (i in seq_len(nrow(results_df))) {
  gene_id <- results_df$gene_id[i]
  cat(sprintf("  [%d/%d] Plotting %s...\n", i, nrow(results_df), gene_id))

  p <- plot_gene_structure(gene_id, gtf_df, expected_mock)
  p <- p + ggtitle(sprintf("#%d: %s", i, gene_id))

  plots[[i]] <- p
}

# Combine into single PDF
combined <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    title = "Real Data Event Detection Test (25 sampled pairs)",
    subtitle = sprintf("SQANTI corrected isoforms from lung cell lines | %s\nColor scheme: GREEN = Reference isoform (A) | PURPLE = Comparison isoform (B)",
                       format(Sys.Date(), "%Y-%m-%d")),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray30", lineheight = 1.3)
    )
  )

output_pdf <- "real_data/real_data_validation.pdf"
ggsave(output_pdf, combined, width = 14, height = nrow(results_df) * 3, limitsize = FALSE)
cat(sprintf("\n✓ Saved visualization: %s\n", output_pdf))

# Save results table
output_tsv <- "real_data/real_data_results.tsv"
results_df %>%
  write_tsv(output_tsv)
cat(sprintf("✓ Saved results: %s\n", output_tsv))

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat("REAL DATA TESTING COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════\n")
