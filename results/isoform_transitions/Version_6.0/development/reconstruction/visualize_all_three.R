#!/usr/bin/env Rscript
#
# Visualize All Three Isoforms
#
# Creates visualizations showing:
#   1. Dominant isoform (target)
#   2. Reconstructed isoform (what we built)
#   3. Comparator isoform (what we started from)
#
# Usage:
#   Rscript visualize_all_three.R <original_gtf> <reconstructed_gtf> <comparator_gtf> <events_tsv> <output_pdf>
#

library(tidyverse)
library(ggplot2)
library(patchwork)

# Source visualization functions
if (file.exists("visualization_functions.R")) {
  source("visualization_functions.R")
} else {
  stop("Cannot find visualization_functions.R in current directory")
}

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   Three-Isoform Visualization                                 ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# ==============================================================================
# Command Line Arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  cat("Usage: visualize_all_three.R <original_gtf> <reconstructed_gtf> <comparator_gtf> <events_tsv> <output_pdf>\n\n")
  cat("Arguments:\n")
  cat("  original_gtf      - GTF with original dominant isoforms\n")
  cat("  reconstructed_gtf - GTF with reconstructed isoforms\n")
  cat("  comparator_gtf    - GTF with comparator isoforms\n")
  cat("  events_tsv        - Events file (to get transcript pairs)\n")
  cat("  output_pdf        - Output PDF path\n")
  quit(status = 1)
}

original_gtf_file <- args[1]
reconstructed_gtf_file <- args[2]
comparator_gtf_file <- args[3]
events_file <- args[4]
output_pdf <- args[5]

# ==============================================================================
# Load Data
# ==============================================================================

cat("Loading data...\n")
cat(sprintf("  Original GTF: %s\n", original_gtf_file))
cat(sprintf("  Reconstructed GTF: %s\n", reconstructed_gtf_file))
cat(sprintf("  Comparator GTF: %s\n", comparator_gtf_file))
cat(sprintf("  Events: %s\n\n", events_file))

original_gtf <- parse_gtf(original_gtf_file)
reconstructed_gtf <- parse_gtf(reconstructed_gtf_file)
comparator_gtf <- parse_gtf(comparator_gtf_file)
events <- read_tsv(events_file, show_col_types = FALSE)

# Get unique transcript pairs
transcript_pairs <- events %>%
  select(gene_id, dominant_transcript_id, comparator_transcript_id) %>%
  distinct()

cat(sprintf("  Original: %d exons from %d transcripts\n",
            nrow(original_gtf),
            n_distinct(original_gtf$transcript_id)))
cat(sprintf("  Reconstructed: %d exons from %d transcripts\n",
            nrow(reconstructed_gtf),
            n_distinct(reconstructed_gtf$transcript_id)))
cat(sprintf("  Comparator: %d exons from %d transcripts\n",
            nrow(comparator_gtf),
            n_distinct(comparator_gtf$transcript_id)))
cat(sprintf("  Transcript pairs: %d\n\n", nrow(transcript_pairs)))

# ==============================================================================
# Create Visualizations
# ==============================================================================

cat("Creating visualizations...\n")

plots <- list()
plot_count <- 0

for (i in seq_len(nrow(transcript_pairs))) {
  pair <- transcript_pairs[i, ]
  gene_id <- pair$gene_id
  dominant_id <- pair$dominant_transcript_id
  comparator_id <- pair$comparator_transcript_id

  cat(sprintf("  [%d/%d] %s... ", i, nrow(transcript_pairs), gene_id))

  # Get dominant exons (original target)
  dominant_exons <- original_gtf %>%
    filter(transcript_id == dominant_id)

  if (nrow(dominant_exons) == 0) {
    cat("⚠ No dominant data\n")
    next
  }

  # Get reconstructed exons (with _reconstructed suffix)
  reconstructed_id <- paste0(dominant_id, "_reconstructed")
  reconstructed_exons <- reconstructed_gtf %>%
    filter(transcript_id == reconstructed_id)

  if (nrow(reconstructed_exons) == 0) {
    cat("⚠ No reconstructed data\n")
    next
  }

  # Get comparator exons (what we started from)
  comparator_exons <- comparator_gtf %>%
    filter(transcript_id == comparator_id)

  if (nrow(comparator_exons) == 0) {
    cat("⚠ No comparator data\n")
    next
  }

  # Combine all three for plotting
  combined_exons <- bind_rows(
    dominant_exons %>% mutate(source_type = "Dominant"),
    reconstructed_exons %>% mutate(source_type = "Reconstructed"),
    comparator_exons %>% mutate(source_type = "Comparator")
  )

  strand <- unique(combined_exons$strand)[1]

  # Prepare plot data with y positions
  # Top to bottom: Dominant (3), Reconstructed (2), Comparator (1)
  plot_data <- combined_exons %>%
    mutate(
      y_position = case_when(
        source_type == "Dominant" ~ 3,
        source_type == "Reconstructed" ~ 2,
        source_type == "Comparator" ~ 1
      )
    )

  # Get coordinate range
  x_min <- min(plot_data$start)
  x_max <- max(plot_data$end)
  x_range <- x_max - x_min
  x_buffer <- x_range * 0.1

  # Prepare introns with chevrons
  intron_data <- prepare_intron_chevrons(plot_data, x_range)

  # Create plot
  p <- ggplot() +
    # Draw intron lines
    {if (!is.null(intron_data$intron_lines)) {
      geom_segment(data = intron_data$intron_lines,
                   aes(x = intron_starts, xend = intron_ends,
                       y = y_position, yend = y_position),
                   linewidth = 0.5, color = "gray40")
    }} +
    # Draw chevrons
    {if (!is.null(intron_data$chevrons)) {
      geom_segment(data = intron_data$chevrons,
                   aes(x = x, xend = xend, y = y, yend = yend),
                   arrow = arrow(length = unit(0.24, "cm"), type = "closed"),
                   color = "gray50", linewidth = 0.8)
    }} +
    # Draw exons
    geom_rect(data = plot_data,
              aes(xmin = start, xmax = end,
                  ymin = y_position - 0.3, ymax = y_position + 0.3,
                  fill = source_type),
              color = "black", linewidth = 0.5) +
    scale_fill_manual(
      values = c("Dominant" = "#27AE60",      # Green - target
                 "Reconstructed" = "#E67E22", # Orange - what we built
                 "Comparator" = "#3498DB"),   # Blue - starting point
      name = "Isoform Type"
    )

  # Add exon coordinates where they fit
  min_width_for_label <- x_range * 0.06
  plot_data_labeled <- plot_data %>%
    mutate(
      exon_width = end - start,
      show_label = exon_width >= min_width_for_label
    ) %>%
    filter(show_label)

  if (nrow(plot_data_labeled) > 0) {
    p <- p +
      geom_text(data = plot_data_labeled,
                aes(x = (start + end) / 2, y = y_position,
                    label = sprintf("%d-%d", start, end)),
                size = 2.5, fontface = "bold", color = "white")
  }

  # Add labels on the left
  label_data <- plot_data %>%
    group_by(source_type, y_position) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(label = case_when(
      source_type == "Dominant" ~ sprintf("Dominant (%s)", dominant_id),
      source_type == "Reconstructed" ~ sprintf("Reconstructed (%s)", reconstructed_id),
      source_type == "Comparator" ~ sprintf("Comparator (%s)", comparator_id)
    ))

  p <- p +
    geom_text(data = label_data,
              aes(x = x_min - x_buffer * 0.3, y = y_position, label = label),
              hjust = 1, size = 3, fontface = "bold") +
    scale_y_continuous(breaks = NULL, limits = c(0.5, 3.5)) +
    coord_cartesian(xlim = c(x_min - x_buffer, x_max + x_buffer), clip = "off") +
    labs(x = sprintf("Genomic Position (Strand: %s)", strand),
         y = NULL) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 10),
      plot.margin = margin(5, 5, 5, 60)
    )

  # Add title
  # Check if dominant matches reconstructed
  match_status <- identical(
    sort(paste(dominant_exons$start, dominant_exons$end, sep = "-")),
    sort(paste(reconstructed_exons$start, reconstructed_exons$end, sep = "-"))
  )

  status_icon <- ifelse(match_status, "✓", "✗")
  status_text <- ifelse(match_status, "MATCH", "MISMATCH")

  title_text <- sprintf("#%d: %s [%s %s]", i, gene_id, status_icon, status_text)

  # Count events for this pair
  n_events <- events %>%
    filter(gene_id == !!gene_id,
           dominant_transcript_id == !!dominant_id,
           comparator_transcript_id == !!comparator_id) %>%
    nrow()

  event_types <- events %>%
    filter(gene_id == !!gene_id,
           dominant_transcript_id == !!dominant_id,
           comparator_transcript_id == !!comparator_id) %>%
    pull(event_type) %>%
    paste(collapse = ", ")

  subtitle_text <- sprintf("Events: %d (%s)", n_events, event_types)

  p <- p +
    ggtitle(title_text, subtitle = subtitle_text) +
    theme(
      plot.title = element_text(size = 11, face = "bold",
                               color = ifelse(match_status, "darkgreen", "darkred")),
      plot.subtitle = element_text(size = 9, color = "gray30")
    )

  plot_count <- plot_count + 1
  plots[[plot_count]] <- p
  cat("✓\n")
}

# ==============================================================================
# Save PDF
# ==============================================================================

if (plot_count == 0) {
  cat("\n⚠ No plots created - no visualization output\n")
  quit(status = 1)
}

cat(sprintf("\nCreated %d plots\n", plot_count))
cat("Combining into PDF...\n")

# Count matches
n_match <- sum(sapply(plots, function(p) {
  grepl("MATCH", p$labels$title) && !grepl("MISMATCH", p$labels$title)
}))
n_mismatch <- plot_count - n_match
match_pct <- 100 * n_match / plot_count

combined_plot <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    title = "Three-Isoform Reconstruction Visualization",
    subtitle = sprintf("Dominant → Reconstructed → Comparator | %d MATCH | %d MISMATCH (%.1f%% match) | %s",
                       n_match, n_mismatch, match_pct,
                       format(Sys.Date(), "%Y-%m-%d")),
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray30")
    )
  )

ggsave(output_pdf, combined_plot,
       width = 14, height = plot_count * 4,
       limitsize = FALSE)

cat(sprintf("\n✓ Saved: %s\n", output_pdf))

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat("VISUALIZATION COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════\n")
