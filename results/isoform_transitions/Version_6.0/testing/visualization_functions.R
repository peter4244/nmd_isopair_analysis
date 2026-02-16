#!/usr/bin/env Rscript
#
# Visualization Functions Library
# Shared functions for creating consistent, publication-ready plots
# of gene structures, exons, introns, and splicing events
#

library(ggplot2)
library(dplyr)

# ==============================================================================
# Gene Structure Visualization
# ==============================================================================

#' Parse GTF file to data frame
#'
#' @param file Path to GTF file
#' @return Tibble with parsed GTF data
parse_gtf <- function(file) {
  lines <- readLines(file)
  lines <- lines[!startsWith(lines, '#') & lines != '']

  results <- list()
  for (line in lines) {
    parts <- strsplit(line, '\t')[[1]]

    # Only parse exon features, skip transcript/gene features
    if (parts[3] != "exon") {
      next
    }

    attrs_str <- parts[9]
    gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', attrs_str)
    transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', attrs_str)
    exon_number <- sub('.*exon_number "([^"]+)".*', '\\1', attrs_str)

    results[[length(results) + 1]] <- tibble(
      seqnames = parts[1],
      type = parts[3],
      start = as.integer(parts[4]),
      end = as.integer(parts[5]),
      strand = parts[7],
      gene_id = gene_id,
      transcript_id = transcript_id,
      exon_number = as.integer(exon_number)
    )
  }
  dplyr::bind_rows(results)
}


#' Create intron segments with UCSC-style directional chevrons
#'
#' @param plot_data Data frame with exon information
#' @param x_range Total genomic range of the plot (for consistent sizing)
#' @param chevron_spacing_pct Spacing between chevrons (as % of plot range, default 0.05 = 5%)
#' @param chevron_size_pct Size of each chevron (as % of plot range, default 0.02 = 2%)
#' @param min_chevrons Minimum chevrons per intron (default 2, switches to dense mode if needed)
#' @return List with intron_lines and chevron_data for plotting
prepare_intron_chevrons <- function(plot_data,
                                   x_range,
                                   chevron_spacing_pct = 0.05,
                                   chevron_size_pct = 0.02,
                                   min_chevrons = 2) {

  # Calculate introns (gaps between exons)
  introns <- plot_data %>%
    group_by(transcript_id) %>%
    arrange(exon_number) %>%
    summarise(
      intron_starts = if (n() > 1) end[1:(n()-1)] else numeric(0),
      intron_ends = if (n() > 1) start[2:n()] else numeric(0),
      y_position = first(y_position),
      strand = first(strand),
      .groups = "drop"
    ) %>%
    filter(length(intron_starts) > 0) %>%
    tidyr::unnest(c(intron_starts, intron_ends))

  if (nrow(introns) == 0) {
    return(list(intron_lines = NULL, chevrons = NULL))
  }

  # STANDARD MODE: chevron spacing and size based on plot range
  standard_spacing_bp <- x_range * chevron_spacing_pct
  standard_size_bp <- x_range * chevron_size_pct

  # DENSE MODE: for short introns (tighter spacing, smaller chevrons)
  dense_spacing_bp <- x_range * (chevron_spacing_pct / 2.5)  # 2.5x tighter
  dense_size_bp <- x_range * (chevron_size_pct / 2)          # 2x smaller

  # Generate chevron positions along each intron
  chevron_list <- list()
  for (i in 1:nrow(introns)) {
    intron <- introns[i, ]
    # Calculate absolute intron length (handles minus strand where start > end)
    intron_length <- abs(intron$intron_ends - intron$intron_starts)
    intron_left <- min(intron$intron_starts, intron$intron_ends)
    intron_right <- max(intron$intron_starts, intron$intron_ends)

    # Determine if we need dense mode
    n_chevrons_standard <- max(1, floor(intron_length / standard_spacing_bp))
    use_dense_mode <- n_chevrons_standard < min_chevrons

    if (use_dense_mode) {
      # DENSE MODE: ensure minimum chevrons with tighter spacing
      spacing_bp <- dense_spacing_bp
      size_bp <- dense_size_bp
      n_chevrons <- max(min_chevrons, floor(intron_length / spacing_bp))
    } else {
      # STANDARD MODE: consistent spacing across plot
      spacing_bp <- standard_spacing_bp
      size_bp <- standard_size_bp
      n_chevrons <- n_chevrons_standard
    }

    # Check if chevron would cover entire intron or overlap
    chevron_total_space <- n_chevrons * size_bp + (n_chevrons - 1) * spacing_bp
    if (chevron_total_space >= intron_length * 0.9) {
      # Chevrons would be too crowded or cover entire intron - use single centered chevron
      n_chevrons <- 1
    }

    # Skip introns where a single chevron would still be too large
    if (n_chevrons == 1 && size_bp >= intron_length * 0.8) {
      next  # Skip this intron entirely
    }

    # Position chevrons evenly along intron (use genomic left-right regardless of strand)
    if (n_chevrons == 1) {
      # Single chevron in middle
      chevron_positions <- (intron_left + intron_right) / 2
    } else {
      # Multiple chevrons with appropriate spacing
      chevron_positions <- seq(intron_left + spacing_bp / 2,
                              intron_right - spacing_bp / 2,
                              length.out = n_chevrons)
    }

    # Create chevron data
    for (pos in chevron_positions) {
      if (intron$strand == "+") {
        # Right-pointing chevron: >
        chevron_list[[length(chevron_list) + 1]] <- tibble(
          x = pos,
          y = intron$y_position,
          xend = pos + size_bp,
          yend = intron$y_position,
          transcript_id = intron$transcript_id,
          strand = intron$strand,
          mode = if (use_dense_mode) "dense" else "standard"
        )
      } else {
        # Left-pointing chevron: <
        chevron_list[[length(chevron_list) + 1]] <- tibble(
          x = pos,
          y = intron$y_position,
          xend = pos - size_bp,
          yend = intron$y_position,
          transcript_id = intron$transcript_id,
          strand = intron$strand,
          mode = if (use_dense_mode) "dense" else "standard"
        )
      }
    }
  }

  chevrons <- if (length(chevron_list) > 0) dplyr::bind_rows(chevron_list) else NULL

  list(intron_lines = introns, chevrons = chevrons)
}


#' Transform genomic coordinates to visual coordinates with compressed introns
#'
#' @param plot_data Data frame with exon coordinates
#' @param intron_compression Compression factor for introns (default 10 = 10x compression)
#' @return Data frame with visual coordinates and mapping information
create_visual_coordinates <- function(plot_data, intron_compression = 10) {

  # Process each transcript separately
  result_list <- list()

  for (tx in unique(plot_data$transcript_id)) {
    tx_data <- plot_data %>%
      filter(transcript_id == tx) %>%
      arrange(exon_number)

    # Create visual coordinate mapping
    visual_coords <- tibble(
      transcript_id = character(),
      exon_number = integer(),
      genomic_start = integer(),
      genomic_end = integer(),
      visual_start = numeric(),
      visual_end = numeric(),
      is_exon = logical()
    )

    visual_pos <- 0  # Current position in visual space

    for (i in 1:nrow(tx_data)) {
      exon <- tx_data[i, ]

      # Add exon at full scale
      exon_width <- exon$end - exon$start
      visual_coords <- bind_rows(visual_coords, tibble(
        transcript_id = exon$transcript_id,
        exon_number = exon$exon_number,
        genomic_start = exon$start,
        genomic_end = exon$end,
        visual_start = visual_pos,
        visual_end = visual_pos + exon_width,
        is_exon = TRUE
      ))
      visual_pos <- visual_pos + exon_width

      # Add compressed intron (if not last exon)
      if (i < nrow(tx_data)) {
        intron_start <- exon$end
        intron_end <- tx_data$start[i + 1]
        intron_width <- intron_end - intron_start
        compressed_width <- intron_width / intron_compression

        visual_coords <- bind_rows(visual_coords, tibble(
          transcript_id = exon$transcript_id,
          exon_number = NA_integer_,
          genomic_start = intron_start,
          genomic_end = intron_end,
          visual_start = visual_pos,
          visual_end = visual_pos + compressed_width,
          is_exon = FALSE
        ))
        visual_pos <- visual_pos + compressed_width
      }
    }

    # Add strand and y_position info
    visual_coords <- visual_coords %>%
      mutate(
        strand = unique(tx_data$strand)[1],
        y_position = unique(tx_data$y_position)[1]
      )

    result_list[[tx]] <- visual_coords
  }

  bind_rows(result_list)
}


#' Plot gene structure for a single test case
#'
#' @param gene_name Gene/test case name
#' @param gtf_data Parsed GTF data
#' @param expected_data Expected events data frame
#' @param show_coordinates Show exon coordinates as text (default TRUE)
#' @param show_chevrons Show UCSC-style directional chevrons on introns (default TRUE)
#' @param compress_introns Use compressed intron visualization (default TRUE)
#' @param intron_compression Compression factor for introns (default 10)
#' @return ggplot object
plot_gene_structure <- function(gene_name,
                               gtf_data,
                               expected_data,
                               show_coordinates = TRUE,
                               show_chevrons = TRUE,
                               compress_introns = TRUE,
                               intron_compression = 10) {

  # Get exons for this gene
  gene_exons <- gtf_data %>%
    filter(gene_id == gene_name) %>%
    arrange(transcript_id, exon_number)

  if (nrow(gene_exons) == 0) return(NULL)

  # Get expected event info
  expected_info <- expected_data %>%
    filter(gene_id == gene_name)

  strand <- unique(gene_exons$strand)[1]

  # Prepare data for plotting
  plot_data <- gene_exons %>%
    mutate(
      transcript_num = as.numeric(factor(transcript_id)),
      y_position = transcript_num
    )

  # Get coordinate range for plot
  x_min <- min(plot_data$start)
  x_max <- max(plot_data$end)
  x_range <- x_max - x_min
  x_buffer <- x_range * 0.1

  # Prepare introns with chevrons (using plot range for consistent sizing)
  intron_data <- prepare_intron_chevrons(plot_data, x_range)

  # Create base plot
  p <- ggplot() +
    # Draw intron lines
    {if (!is.null(intron_data$intron_lines)) {
      geom_segment(data = intron_data$intron_lines,
                   aes(x = intron_starts, xend = intron_ends,
                       y = y_position, yend = y_position),
                   linewidth = 0.5, color = "gray40")
    }} +
    # Draw UCSC-style chevrons on introns
    {if (show_chevrons && !is.null(intron_data$chevrons)) {
      geom_segment(data = intron_data$chevrons,
                   aes(x = x, xend = xend, y = y, yend = yend),
                   arrow = arrow(length = unit(0.24, "cm"), type = "closed"),
                   color = "gray50", linewidth = 0.8)
    }} +
    # Draw exons as rectangles
    geom_rect(data = plot_data,
              aes(xmin = start, xmax = end,
                  ymin = y_position - 0.3, ymax = y_position + 0.3,
                  fill = transcript_id),
              color = "black", linewidth = 0.5)

  # Add exon coordinates if requested (only if they fit)
  if (show_coordinates) {
    # Calculate which exons are wide enough for labels
    # Rough estimate: need ~8-10 characters, each ~0.7% of x_range at size 2.5
    min_width_for_label <- x_range * 0.06  # Need ~6% of range for label to fit

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
                  size = 2.5, fontface = "bold")
    }
  }

  # Add transcript labels
  p <- p +
    geom_text(data = plot_data %>%
                group_by(transcript_id, y_position) %>%
                slice(1),
              aes(x = x_min - x_buffer * 0.5, y = y_position,
                  label = transcript_id),
              hjust = 1, size = 3, fontface = "bold")

  # Add TSS/TES markers
  if (strand == "+") {
    p <- p +
      # TSS at left (5' end)
      geom_segment(data = plot_data %>%
                     group_by(transcript_id) %>%
                     slice_min(start),
                   aes(x = start, xend = start,
                       y = y_position - 0.5, yend = y_position + 0.5),
                   arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
                   color = "darkgreen", linewidth = 1) +
      # TES at right (3' end)
      geom_point(data = plot_data %>%
                   group_by(transcript_id) %>%
                   slice_max(end),
                 aes(x = end, y = y_position),
                 shape = 15, size = 3, color = "darkred")
  } else {
    p <- p +
      # TSS at right (3' end on minus)
      geom_segment(data = plot_data %>%
                     group_by(transcript_id) %>%
                     slice_max(end),
                   aes(x = end, xend = end,
                       y = y_position - 0.5, yend = y_position + 0.5),
                   arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
                   color = "darkgreen", linewidth = 1) +
      # TES at left (5' end on minus)
      geom_point(data = plot_data %>%
                   group_by(transcript_id) %>%
                   slice_min(start),
                 aes(x = start, y = y_position),
                 shape = 15, size = 3, color = "darkred")
  }

  # Styling
  # Generate enough colors for all transcripts
  n_transcripts <- length(unique(plot_data$transcript_id))
  colors <- if (n_transcripts <= 2) {
    c("#7fc97f", "#beaed4")
  } else {
    # Use a color palette that can handle many colors
    scales::hue_pal()(n_transcripts)
  }

  p <- p +
    scale_fill_manual(values = colors) +
    scale_y_continuous(breaks = NULL) +
    coord_cartesian(xlim = c(x_min - x_buffer, x_max + x_buffer)) +
    labs(
      title = sprintf("%s: %s", gene_name, expected_info$description),
      subtitle = sprintf("Strand: %s  |  Expected: %s",
                        strand, expected_info$event_type),
      x = "Genomic coordinate",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "gray30"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )

  # Add legend annotation
  legend_text <- sprintf(
    "Legend: ▲ TSS (green), ■ TES (red), %s = Direction, Boxes = Exons, Lines = Introns\nNotes: %s",
    if (strand == "+") ">>>" else "<<<",
    expected_info$notes
  )

  p <- p +
    labs(caption = legend_text) +
    theme(plot.caption = element_text(size = 7, hjust = 0, color = "gray40"))

  return(p)
}


#' Create text summary for a gene
#'
#' @param gene_name Gene/test case name
#' @param gtf_data Parsed GTF data
#' @param expected_data Expected events data frame
#' @return Character vector with formatted text summary
create_gene_text_summary <- function(gene_name, gtf_data, expected_data) {
  gene_exons <- gtf_data %>% filter(gene_id == gene_name)
  expected_info <- expected_data %>% filter(gene_id == gene_name)
  strand <- unique(gene_exons$strand)[1]

  output <- c(
    sprintf("TEST CASE: %s", gene_name),
    sprintf("Description: %s", expected_info$description),
    sprintf("Strand: %s", strand),
    sprintf("Expected events: %s", expected_info$event_type),
    sprintf("Notes: %s\n", expected_info$notes)
  )

  # Show exon structures
  transcripts <- unique(gene_exons$transcript_id)
  for (tx in transcripts) {
    tx_exons <- gene_exons %>%
      filter(transcript_id == tx) %>%
      arrange(exon_number)

    output <- c(output, sprintf("  %s:", tx))
    for (i in 1:nrow(tx_exons)) {
      exon <- tx_exons[i, ]
      output <- c(output,
                  sprintf("    Exon %d: [%d-%d] (%d bp)",
                          exon$exon_number, exon$start, exon$end,
                          exon$end - exon$start + 1))
    }
    output <- c(output, "")
  }

  output <- c(output,
              "─────────────────────────────────────────────────────────────────\n")

  return(output)
}
