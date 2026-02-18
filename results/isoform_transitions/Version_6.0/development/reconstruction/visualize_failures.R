#!/usr/bin/env Rscript
#
# Visualize Failing Reconstructions
#
# Picks N failures from verification_new.tsv and plots dominant / reconstructed / comparator
# using the same three-track format as visualize_all_three.R
#

library(tidyverse)
library(ggplot2)
library(patchwork)

source("visualization_functions.R")

# ── Arguments ───────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  cat("Usage: visualize_failures.R <original_gtf> <reconstructed_gtf> <comparator_gtf> <verification_tsv> <output_pdf>\n")
  quit(status = 1)
}

original_gtf_file     <- args[1]
reconstructed_gtf_file <- args[2]
comparator_gtf_file   <- args[3]
verification_file     <- args[4]
output_pdf            <- args[5]

N_PLOTS <- 10   # how many failures to show

# ── Load data ────────────────────────────────────────────────────────────────
cat("Loading GTFs...\n")
original_gtf     <- parse_gtf(original_gtf_file)
reconstructed_gtf <- parse_gtf(reconstructed_gtf_file)
comparator_gtf   <- parse_gtf(comparator_gtf_file)

cat("Loading verification...\n")
verification <- read_tsv(verification_file, show_col_types = FALSE)

# Pick N diverse failures (spread across different dominant transcripts)
failures <- verification %>%
  filter(status == "FAIL") %>%
  arrange(gene_id, dominant_transcript_id) %>%
  group_by(dominant_transcript_id) %>%
  slice(1) %>%          # one failure per dominant (different comparators)
  ungroup() %>%
  slice_head(n = N_PLOTS)

cat(sprintf("Selected %d failures to visualise\n\n", nrow(failures)))

# ── Plot each failure ────────────────────────────────────────────────────────
plots <- list()

for (i in seq_len(nrow(failures))) {
  row <- failures[i, ]
  gene_id        <- row$gene_id
  dominant_id    <- row$dominant_transcript_id
  comparator_id  <- row$comparator_transcript_id
  reconstructed_id <- paste0(dominant_id, "::", comparator_id)

  cat(sprintf("[%d/%d] %s  %s -> %s\n", i, nrow(failures), gene_id, comparator_id, dominant_id))

  dominant_exons    <- original_gtf     %>% filter(transcript_id == dominant_id)
  reconstructed_exons <- reconstructed_gtf %>% filter(transcript_id == reconstructed_id)
  comparator_exons  <- comparator_gtf   %>% filter(transcript_id == comparator_id)

  if (nrow(dominant_exons) == 0 || nrow(reconstructed_exons) == 0 || nrow(comparator_exons) == 0) {
    cat("  ⚠ missing exon data, skipping\n")
    next
  }

  combined_exons <- bind_rows(
    dominant_exons    %>% mutate(source_type = "Dominant"),
    reconstructed_exons %>% mutate(source_type = "Reconstructed"),
    comparator_exons  %>% mutate(source_type = "Comparator")
  )

  strand <- unique(combined_exons$strand)[1]

  plot_data <- combined_exons %>%
    mutate(
      y_position = case_when(
        source_type == "Dominant"      ~ 3,
        source_type == "Reconstructed" ~ 2,
        source_type == "Comparator"    ~ 1
      )
    )

  x_min   <- min(plot_data$start)
  x_max   <- max(plot_data$end)
  x_range <- x_max - x_min
  x_buffer <- x_range * 0.1

  intron_data <- prepare_intron_chevrons(plot_data, x_range)

  p <- ggplot() +
    {if (!is.null(intron_data$intron_lines))
      geom_segment(data = intron_data$intron_lines,
                   aes(x = intron_starts, xend = intron_ends,
                       y = y_position, yend = y_position),
                   linewidth = 0.5, color = "gray40")} +
    {if (!is.null(intron_data$chevrons))
      geom_segment(data = intron_data$chevrons,
                   aes(x = x, xend = xend, y = y, yend = yend),
                   arrow = arrow(length = unit(0.24, "cm"), type = "closed"),
                   color = "gray50", linewidth = 0.8)} +
    geom_rect(data = plot_data,
              aes(xmin = start, xmax = end,
                  ymin = y_position - 0.3, ymax = y_position + 0.3,
                  fill = source_type),
              color = "black", linewidth = 0.5) +
    scale_fill_manual(
      values = c("Dominant"      = "#27AE60",
                 "Reconstructed" = "#E67E22",
                 "Comparator"    = "#3498DB"),
      name = "Isoform"
    )

  # Coordinate labels inside wide exons
  min_width_for_label <- x_range * 0.06
  plot_data_labeled <- plot_data %>%
    mutate(exon_width = end - start) %>%
    filter(exon_width >= min_width_for_label)

  if (nrow(plot_data_labeled) > 0) {
    p <- p +
      geom_text(data = plot_data_labeled,
                aes(x = (start + end) / 2, y = y_position,
                    label = sprintf("%d-%d", start, end)),
                size = 2.5, fontface = "bold", color = "white")
  }

  # Track labels on the left
  label_data <- tibble(
    y_position = c(3, 2, 1),
    label = c(
      sprintf("Dominant\n(%s)", dominant_id),
      sprintf("Reconstructed\n(%s)", reconstructed_id),
      sprintf("Comparator\n(%s)", comparator_id)
    )
  )

  p <- p +
    geom_text(data = label_data,
              aes(x = x_min - x_buffer * 0.3, y = y_position, label = label),
              hjust = 1, size = 2.8, fontface = "bold") +
    scale_y_continuous(breaks = NULL, limits = c(0.3, 3.7)) +
    coord_cartesian(xlim = c(x_min - x_buffer, x_max + x_buffer), clip = "off") +
    labs(x = sprintf("Genomic position (%s strand) — %s", strand, row$gene_id),
         y = NULL) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor    = element_blank(),
      legend.position     = "top",
      plot.margin         = margin(5, 5, 5, 80)
    ) +
    ggtitle(
      sprintf("[FAIL] %s  |  dom: %d exons  recon: %d exons  (%+d)",
              dominant_id,
              row$n_exons_original,
              row$n_exons_reconstructed,
              row$n_exons_reconstructed - row$n_exons_original),
      subtitle = sprintf("Comparator: %s  |  %s", comparator_id, row$reason)
    ) +
    theme(
      plot.title    = element_text(size = 10, face = "bold", color = "darkred"),
      plot.subtitle = element_text(size = 8,  color = "gray30")
    )

  plots[[length(plots) + 1]] <- p
}

# ── Save PDF ─────────────────────────────────────────────────────────────────
n <- length(plots)
if (n == 0) { cat("No plots.\n"); quit(status = 1) }

combined <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    title    = "Failing Reconstructions",
    subtitle = sprintf("%d failures shown | %s", n, format(Sys.Date(), "%Y-%m-%d")),
    theme    = theme(plot.title    = element_text(size = 16, face = "bold"),
                     plot.subtitle = element_text(size = 11, color = "gray40"))
  )

ggsave(output_pdf, combined, width = 16, height = n * 5, limitsize = FALSE)
cat(sprintf("\n✓ Saved %d-page PDF: %s\n", n, output_pdf))
