#!/usr/bin/env Rscript
#
# Visualize Verification Results
#
# Creates side-by-side visualizations of:
#   - Original dominant isoform
#   - Reconstructed isoform
#
# Highlights mismatches and shows which tests passed/failed
#
# Usage:
#   Rscript visualize_verification.R <original_gtf> <reconstructed_gtf> <verification_tsv> <output_pdf>
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
cat("║   Verification Results Visualization                         ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# ==============================================================================
# Command Line Arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  cat("Usage: visualize_verification.R <original_gtf> <reconstructed_gtf> <verification_tsv> <output_pdf>\n\n")
  cat("Arguments:\n")
  cat("  original_gtf       - GTF with original dominant isoforms\n")
  cat("  reconstructed_gtf  - GTF with reconstructed isoforms\n")
  cat("  verification_tsv   - Verification results file\n")
  cat("  output_pdf         - Output PDF path\n")
  quit(status = 1)
}

original_gtf_file <- args[1]
reconstructed_gtf_file <- args[2]
verification_file <- args[3]
output_pdf <- args[4]

# ==============================================================================
# Load Data
# ==============================================================================

cat("Loading data...\n")
cat(sprintf("  Original GTF: %s\n", original_gtf_file))
cat(sprintf("  Reconstructed GTF: %s\n", reconstructed_gtf_file))
cat(sprintf("  Verification: %s\n\n", verification_file))

original_gtf <- parse_gtf(original_gtf_file)
reconstructed_gtf <- parse_gtf(reconstructed_gtf_file)
verification <- read_tsv(verification_file, show_col_types = FALSE)

cat(sprintf("  Original: %d exons from %d transcripts\n",
            nrow(original_gtf),
            n_distinct(original_gtf$transcript_id)))
cat(sprintf("  Reconstructed: %d exons from %d transcripts\n",
            nrow(reconstructed_gtf),
            n_distinct(reconstructed_gtf$transcript_id)))
cat(sprintf("  Verification: %d test cases\n\n", nrow(verification)))

# ==============================================================================
# Create Visualizations
# ==============================================================================

# Summary statistics
n_pass <- sum(verification$status == "PASS")
n_fail <- sum(verification$status == "FAIL")
n_error <- sum(verification$status == "ERROR")
pass_pct <- 100 * n_pass / nrow(verification)

cat(sprintf("Results: %d PASS, %d FAIL, %d ERROR (%.1f%% pass)\n\n",
            n_pass, n_fail, n_error, pass_pct))

cat("Creating visualizations...\n")

plots <- list()
plot_count <- 0

for (i in seq_len(nrow(verification))) {
  test_case <- verification[i, ]
  gene_id <- test_case$gene_id
  dominant_id <- test_case$dominant_transcript_id
  comparator_id <- test_case$comparator_transcript_id
  status <- test_case$status

  cat(sprintf("  [%d/%d] %s [%s]... ", i, nrow(verification), gene_id, status))

  # Get original dominant exons
  original_exons <- original_gtf %>%
    filter(transcript_id == dominant_id)

  if (nrow(original_exons) == 0) {
    cat("⚠ No original data\n")
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

  # Combine for plotting
  combined_exons <- bind_rows(
    original_exons %>% mutate(source_type = "Original"),
    reconstructed_exons %>% mutate(source_type = "Reconstructed")
  )

  strand <- unique(combined_exons$strand)[1]

  # Prepare plot data
  plot_data <- combined_exons %>%
    mutate(
      y_position = ifelse(source_type == "Original", 2, 1)
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
      values = c("Original" = "#4A90E2", "Reconstructed" = "#E67E22"),
      name = "Source"
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

  # Add labels
  label_data <- plot_data %>%
    group_by(source_type, y_position) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(label = ifelse(source_type == "Original",
                         sprintf("Original (%s)", dominant_id),
                         sprintf("Reconstructed (%s)", reconstructed_id)))

  p <- p +
    geom_text(data = label_data,
              aes(x = x_min - x_buffer * 0.3, y = y_position, label = label),
              hjust = 1, size = 3, fontface = "bold") +
    scale_y_continuous(breaks = NULL, limits = c(0.5, 2.5)) +
    coord_cartesian(xlim = c(x_min - x_buffer, x_max + x_buffer), clip = "off") +
    labs(x = sprintf("Genomic Position (Strand: %s)", strand),
         y = NULL) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 10),
      plot.margin = margin(5, 5, 5, 50)
    )

  # Add title with status
  status_icon <- case_when(
    status == "PASS" ~ "✓",
    status == "FAIL" ~ "✗",
    status == "ERROR" ~ "⚠",
    TRUE ~ "?"
  )

  title_text <- sprintf("#%d: %s [%s %s]", i, gene_id, status_icon, status)

  subtitle_parts <- c(
    sprintf("Comparator: %s", comparator_id)
  )

  if (status != "PASS") {
    subtitle_parts <- c(subtitle_parts, sprintf("Issue: %s", test_case$reason))
  }

  subtitle_text <- paste(subtitle_parts, collapse = " | ")

  p <- p +
    ggtitle(title_text, subtitle = subtitle_text) +
    theme(
      plot.title = element_text(size = 11, face = "bold",
                               color = ifelse(status == "PASS", "darkgreen", "darkred")),
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

combined_plot <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    title = "Reconstruction Verification Results",
    subtitle = sprintf("%d PASS | %d FAIL | %d ERROR (%.1f%% pass) | %s",
                       n_pass, n_fail, n_error, pass_pct,
                       format(Sys.Date(), "%Y-%m-%d")),
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray30")
    )
  )

ggsave(output_pdf, combined_plot,
       width = 14, height = plot_count * 3,
       limitsize = FALSE)

cat(sprintf("\n✓ Saved: %s\n", output_pdf))

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat("VISUALIZATION COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════\n")
