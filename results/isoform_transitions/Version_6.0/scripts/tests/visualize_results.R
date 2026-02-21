#!/usr/bin/env Rscript
#
# Visualize Test Results — Dominant, Comparator, Reconstructed
#
# Creates per-pair plots showing all three isoform structures with mismatch
# regions highlighted. Generates a multi-page PDF.
#
# Usage:
#   Rscript scripts/tests/visualize_results.R <gtf> <pairs> <output.pdf>
#   Rscript scripts/tests/visualize_results.R <gtf> <pairs> <output.pdf> --ue <union_exons.rds>
#   Rscript scripts/tests/visualize_results.R                # defaults: test_data + results.pdf
#   Rscript scripts/tests/visualize_results.R --failures     # only show FAIL/ERROR pairs
#
# Requires: ggplot2, patchwork

library(tidyverse)
library(ggplot2)
library(patchwork)

# ==============================================================================
# Resolve paths
# ==============================================================================

script_dir <- if (interactive()) {
  "scripts/tests"
} else {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args_all[grep("--file=", args_all)])
  if (length(script_path) > 0) dirname(script_path) else "scripts/tests"
}

project_root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)

# Source the test runner (defines parse_gtf, build_union_exons_from_gtf_data,
# verify_transcript, and sources production functions)
source(file.path(project_root, "scripts/core/event_detection_functions.R"))
source(file.path(project_root, "scripts/core/reconstruction_functions.R"))

default_gtf <- file.path(project_root, "scripts/tests/test_data/test_cases.gtf")
default_pairs <- file.path(project_root, "scripts/tests/test_data/pairs.tsv")
default_output <- file.path(project_root, "scripts/tests/test_data/results.pdf")

# ==============================================================================
# GTF Parser (same as run_tests.R)
# ==============================================================================

parse_gtf <- function(gtf_file) {
  lines <- readLines(gtf_file)
  lines <- lines[!startsWith(lines, '#') & lines != '']
  results <- list()
  for (line in lines) {
    parts <- strsplit(line, '\t')[[1]]
    if (parts[3] != "exon") next
    attrs_str <- parts[9]
    gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', attrs_str)
    transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', attrs_str)
    exon_number_str <- sub('.*exon_number "([^"]+)".*', '\\1', attrs_str)
    exon_number <- if (grepl('exon_number', attrs_str)) as.integer(exon_number_str) else NA_integer_
    results[[length(results) + 1]] <- tibble(
      chr = parts[1], exon_start = as.integer(parts[4]),
      exon_end = as.integer(parts[5]), strand = parts[7],
      gene_id = gene_id, transcript_id = transcript_id,
      exon_number = exon_number)
  }
  bind_rows(results)
}

# ==============================================================================
# Union Exon Builder (same as run_tests.R)
# ==============================================================================

build_union_exons_from_gtf_data <- function(gtf_data) {
  genes <- unique(gtf_data$gene_id)
  all_union_exons <- list()
  ue_counter <- 1
  for (gene in genes) {
    gene_exons <- gtf_data %>% filter(gene_id == !!gene)
    chr <- unique(gene_exons$chr)[1]
    strand <- unique(gene_exons$strand)[1]
    all_boundaries <- sort(unique(c(gene_exons$exon_start, gene_exons$exon_end)))
    if (length(all_boundaries) < 2) next
    boundary_types <- sapply(all_boundaries, function(b) {
      is_start <- b %in% gene_exons$exon_start
      is_end <- b %in% gene_exons$exon_end
      if (is_start && is_end) return("BOTH")
      if (is_start) return("START")
      return("END")
    })
    n <- length(all_boundaries)
    seg_starts <- c(); seg_ends <- c()
    current_start <- all_boundaries[1]
    for (i in seq_len(n - 2) + 1) {
      b <- all_boundaries[i]
      if (boundary_types[i] == "START") {
        seg_starts <- c(seg_starts, current_start); seg_ends <- c(seg_ends, b - 1)
        current_start <- b
      } else if (boundary_types[i] == "END") {
        seg_starts <- c(seg_starts, current_start); seg_ends <- c(seg_ends, b)
        current_start <- b + 1
      } else {
        seg_starts <- c(seg_starts, current_start); seg_ends <- c(seg_ends, b)
        current_start <- b + 1
      }
    }
    seg_starts <- c(seg_starts, current_start); seg_ends <- c(seg_ends, all_boundaries[n])
    segments <- tibble(seg_start = seg_starts, seg_end = seg_ends) %>% filter(seg_start <= seg_end)
    covered <- sapply(seq_len(nrow(segments)), function(j) {
      any(gene_exons$exon_start <= segments$seg_start[j] & gene_exons$exon_end >= segments$seg_end[j])
    })
    covered_segments <- segments[covered, ]
    for (i in seq_len(nrow(covered_segments))) {
      all_union_exons[[ue_counter]] <- tibble(chr = chr, start = covered_segments$seg_start[i],
        end = covered_segments$seg_end[i], union_exon_id = sprintf("ue_%d", ue_counter),
        strand = strand, gene_id = gene)
      ue_counter <- ue_counter + 1
    }
  }
  bind_rows(all_union_exons)
}

# ==============================================================================
# Mismatch Detection
# ==============================================================================

#' Compare two sorted exon sets and identify mismatches
#' Returns a tibble with match_status for each exon position
identify_mismatches <- function(dom_exons, recon_exons, strand) {
  dom_sorted <- dom_exons %>% arrange(exon_start)
  recon_sorted <- recon_exons %>% arrange(exon_start)

  n_dom <- nrow(dom_sorted)
  n_recon <- nrow(recon_sorted)

  mismatches <- list()

  if (n_dom != n_recon) {
    # Different exon counts — mark all reconstructed exons
    for (i in seq_len(n_recon)) {
      mismatches[[i]] <- tibble(
        exon_start = recon_sorted$exon_start[i],
        exon_end = recon_sorted$exon_end[i],
        mismatch = TRUE,
        detail = sprintf("Exon count: %d orig vs %d recon", n_dom, n_recon)
      )
    }
    return(bind_rows(mismatches))
  }

  for (i in seq_len(n_dom)) {
    is_first <- (i == 1)
    is_last <- (i == n_dom)

    if (strand == "+") {
      start_tol <- if (is_first) TSS_TOLERANCE else 0L
      end_tol <- if (is_last) TES_TOLERANCE else 0L
    } else {
      start_tol <- if (is_first) TES_TOLERANCE else 0L
      end_tol <- if (is_last) TSS_TOLERANCE else 0L
    }

    start_diff <- abs(dom_sorted$exon_start[i] - recon_sorted$exon_start[i])
    end_diff <- abs(dom_sorted$exon_end[i] - recon_sorted$exon_end[i])

    has_mismatch <- start_diff > start_tol || end_diff > end_tol
    detail <- if (has_mismatch) {
      sprintf("Exon %d: [%d-%d] vs [%d-%d]", i,
              dom_sorted$exon_start[i], dom_sorted$exon_end[i],
              recon_sorted$exon_start[i], recon_sorted$exon_end[i])
    } else ""

    mismatches[[i]] <- tibble(
      exon_start = recon_sorted$exon_start[i],
      exon_end = recon_sorted$exon_end[i],
      mismatch = has_mismatch,
      detail = detail
    )
  }
  bind_rows(mismatches)
}

# ==============================================================================
# Plotting Function
# ==============================================================================

plot_three_isoforms <- function(dom_exons, comp_exons, recon_exons, events,
                                gene_id, dom_id, comp_id, strand,
                                pair_idx, status, reason) {

  # Prepare plot data
  dom_plot <- dom_exons %>%
    mutate(track = "Dominant", y_pos = 3, label = sprintf("Dominant (%s)", dom_id))
  recon_plot <- recon_exons %>%
    mutate(track = "Reconstructed", y_pos = 2, label = "Reconstructed")
  comp_plot <- comp_exons %>%
    mutate(track = "Comparator", y_pos = 1, label = sprintf("Comparator (%s)", comp_id))

  all_exons <- bind_rows(dom_plot, recon_plot, comp_plot)

  x_min <- min(all_exons$exon_start)
  x_max <- max(all_exons$exon_end)
  x_range <- x_max - x_min
  x_buffer <- max(x_range * 0.05, 50)

  # Identify mismatches between dominant and reconstructed
  mismatch_data <- identify_mismatches(dom_exons, recon_exons, strand)
  mismatched_exons <- mismatch_data %>% filter(mismatch)

  # Build base plot
  p <- ggplot() +
    # Intron lines for each track
    geom_segment(data = dom_plot %>% arrange(exon_start) %>%
                   mutate(next_start = lead(exon_start)) %>% filter(!is.na(next_start)),
                 aes(x = exon_end, xend = next_start, y = y_pos, yend = y_pos),
                 color = "gray50", linewidth = 0.4) +
    geom_segment(data = recon_plot %>% arrange(exon_start) %>%
                   mutate(next_start = lead(exon_start)) %>% filter(!is.na(next_start)),
                 aes(x = exon_end, xend = next_start, y = y_pos, yend = y_pos),
                 color = "gray50", linewidth = 0.4) +
    geom_segment(data = comp_plot %>% arrange(exon_start) %>%
                   mutate(next_start = lead(exon_start)) %>% filter(!is.na(next_start)),
                 aes(x = exon_end, xend = next_start, y = y_pos, yend = y_pos),
                 color = "gray50", linewidth = 0.4) +
    # Exon rectangles
    geom_rect(data = all_exons,
              aes(xmin = exon_start, xmax = exon_end,
                  ymin = y_pos - 0.25, ymax = y_pos + 0.25,
                  fill = track),
              color = "black", linewidth = 0.3) +
    scale_fill_manual(
      values = c("Dominant" = "#27AE60",
                 "Reconstructed" = "#E67E22",
                 "Comparator" = "#3498DB"),
      name = NULL
    )

  # Highlight mismatched exons in the reconstructed track
  if (nrow(mismatched_exons) > 0) {
    p <- p +
      geom_rect(data = mismatched_exons,
                aes(xmin = exon_start, xmax = exon_end,
                    ymin = 2 - 0.35, ymax = 2 + 0.35),
                fill = NA, color = "red", linewidth = 1.2, linetype = "solid") +
      # Add red bracket connecting mismatched regions to dominant
      geom_segment(data = mismatched_exons,
                   aes(x = (exon_start + exon_end) / 2,
                       xend = (exon_start + exon_end) / 2,
                       y = 2.35, yend = 2.65),
                   color = "red", linewidth = 0.6, linetype = "dashed")
  }

  # Track labels
  label_df <- tibble(
    y_pos = c(3, 2, 1),
    label = c(sprintf("Dominant\n(%s)", dom_id),
              "Reconstructed",
              sprintf("Comparator\n(%s)", comp_id))
  )

  p <- p +
    geom_text(data = label_df,
              aes(x = x_min - x_buffer * 1.5, y = y_pos, label = label),
              hjust = 1, size = 2.5, fontface = "bold") +
    scale_y_continuous(breaks = NULL, limits = c(0.3, 3.7)) +
    coord_cartesian(xlim = c(x_min - x_buffer * 3, x_max + x_buffer), clip = "off") +
    labs(x = sprintf("Genomic position (%s strand)", strand), y = NULL) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.margin = margin(5, 10, 5, 80)
    )

  # Title and subtitle
  event_str <- if (nrow(events) > 0) {
    paste(events$event_type, collapse = ", ")
  } else {
    "none"
  }

  title_color <- if (status == "PASS") "darkgreen" else "darkred"
  status_symbol <- if (status == "PASS") "PASS" else status

  p <- p +
    ggtitle(
      sprintf("#%d [%s] %s", pair_idx, status_symbol, gene_id),
      subtitle = sprintf("Events: %s | %s", event_str, reason)
    ) +
    theme(
      plot.title = element_text(size = 10, face = "bold", color = title_color),
      plot.subtitle = element_text(size = 8, color = "gray30")
    )

  return(p)
}

# ==============================================================================
# Main
# ==============================================================================

visualize_tests <- function(gtf_file, pairs_file, output_pdf, ue_rds_file = NULL,
                            failures_only = FALSE) {

  cat("\n======================================================================\n")
  cat("  RECONSTRUCTION VISUALIZATION\n")
  cat("======================================================================\n")
  cat(sprintf("  GTF:    %s\n", gtf_file))
  cat(sprintf("  Pairs:  %s\n", pairs_file))
  cat(sprintf("  Output: %s\n", output_pdf))
  cat("======================================================================\n\n")

  # Parse GTF
  gtf_data <- parse_gtf(gtf_file)
  cat(sprintf("Parsed %d exons, %d genes, %d transcripts\n",
              nrow(gtf_data), n_distinct(gtf_data$gene_id), n_distinct(gtf_data$transcript_id)))

  # Build or load UEs
  gtf_union_exons <- build_union_exons_from_gtf_data(gtf_data)
  if (!is.null(ue_rds_file)) {
    ue_raw <- readRDS(ue_rds_file)
    if ("union_exon_start" %in% names(ue_raw)) {
      prod_ue <- ue_raw %>% rename(start = union_exon_start, end = union_exon_end)
      if ("seqnames" %in% names(prod_ue)) prod_ue <- prod_ue %>% rename(chr = seqnames)
    } else {
      prod_ue <- ue_raw
    }
    prod_genes <- unique(prod_ue$gene_id)
    gtf_only <- setdiff(unique(gtf_union_exons$gene_id), prod_genes)
    union_exons <- bind_rows(prod_ue, gtf_union_exons %>% filter(gene_id %in% gtf_only))
  } else {
    union_exons <- gtf_union_exons
  }

  pairs <- read_tsv(pairs_file, show_col_types = FALSE)
  cat(sprintf("Processing %d pairs...\n\n", nrow(pairs)))

  # Process each pair
  plots <- list()
  plot_count <- 0

  for (i in seq_len(nrow(pairs))) {
    pair <- pairs[i, ]
    dom_exons <- gtf_data %>% filter(gene_id == pair$gene_id, transcript_id == pair$isoform_A)
    comp_exons <- gtf_data %>% filter(gene_id == pair$gene_id, transcript_id == pair$isoform_B)

    if (nrow(dom_exons) == 0 || nrow(comp_exons) == 0) {
      if (!failures_only) {
        cat(sprintf("  [%d] %s: SKIP (missing exons)\n", i, pair$gene_id))
      }
      next
    }

    gene_strand <- unique(dom_exons$strand)[1]

    # Detect events
    events <- tryCatch(
      detect_events_for_pair(dom_exons, comp_exons, pair$gene_id,
                             pair$isoform_A, pair$isoform_B, gene_strand),
      error = function(e) tibble()
    )

    # Reconstruct
    status <- "PASS"
    reason <- ""
    recon_exons <- comp_exons  # Default: no events -> comparator = dominant

    if (nrow(events) > 0) {
      gene_ue <- union_exons %>% filter(gene_id == pair$gene_id)
      comp_for_recon <- comp_exons %>% select(chr, exon_start, exon_end, strand, gene_id, transcript_id)
      recon_exons <- tryCatch(
        reconstruct_dominant_v2(comp_for_recon, events, gene_ue),
        error = function(e) {
          status <<- "ERROR"
          reason <<- e$message
          comp_exons
        }
      )
    }

    # Verify
    if (status != "ERROR") {
      if (nrow(recon_exons) == 0) {
        status <- "ERROR"
        reason <- "Reconstruction produced 0 exons"
        recon_exons <- comp_exons
      } else {
        # Inline verification
        dom_sorted <- dom_exons %>% arrange(exon_start)
        recon_sorted <- recon_exons %>% arrange(exon_start)
        if (nrow(dom_sorted) != nrow(recon_sorted)) {
          status <- "FAIL"
          reason <- sprintf("Exon count: %d orig vs %d recon", nrow(dom_sorted), nrow(recon_sorted))
        } else {
          for (j in seq_len(nrow(dom_sorted))) {
            is_first <- (j == 1); is_last <- (j == nrow(dom_sorted))
            if (gene_strand == "+") {
              s_tol <- if (is_first) TSS_TOLERANCE else 0L
              e_tol <- if (is_last) TES_TOLERANCE else 0L
            } else {
              s_tol <- if (is_first) TES_TOLERANCE else 0L
              e_tol <- if (is_last) TSS_TOLERANCE else 0L
            }
            if (abs(dom_sorted$exon_start[j] - recon_sorted$exon_start[j]) > s_tol ||
                abs(dom_sorted$exon_end[j] - recon_sorted$exon_end[j]) > e_tol) {
              status <- "FAIL"
              reason <- sprintf("Exon %d: [%d-%d] vs [%d-%d]", j,
                                dom_sorted$exon_start[j], dom_sorted$exon_end[j],
                                recon_sorted$exon_start[j], recon_sorted$exon_end[j])
              break
            }
          }
          if (status == "PASS") reason <- "Match"
        }
      }
    }

    if (failures_only && status == "PASS") next

    cat(sprintf("  [%d] %s: %s\n", i, pair$gene_id, status))

    p <- tryCatch(
      plot_three_isoforms(dom_exons, comp_exons, recon_exons, events,
                          pair$gene_id, pair$isoform_A, pair$isoform_B,
                          gene_strand, i, status, reason),
      error = function(e) {
        cat(sprintf("    Plot error: %s\n", e$message))
        NULL
      }
    )

    if (!is.null(p)) {
      plot_count <- plot_count + 1
      plots[[plot_count]] <- p
    }
  }

  if (plot_count == 0) {
    cat("\nNo plots to create.\n")
    return(invisible(NULL))
  }

  # Combine and save
  cat(sprintf("\nSaving %d plots to %s...\n", plot_count, output_pdf))

  # Use patchwork to combine
  combined <- wrap_plots(plots, ncol = 1) +
    plot_annotation(
      title = "Reconstruction Test Results",
      subtitle = sprintf("%d plots | %s", plot_count, Sys.Date()),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray40")
      )
    )

  ggsave(output_pdf, combined,
         width = 14, height = plot_count * 3.5,
         limitsize = FALSE)

  cat(sprintf("Saved: %s\n\n", output_pdf))
}

# ==============================================================================
# CLI
# ==============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  # Parse flags
  ue_rds_file <- NULL
  failures_only <- FALSE

  if ("--ue" %in% args) {
    ue_idx <- which(args == "--ue")
    if (ue_idx < length(args)) {
      ue_rds_file <- args[ue_idx + 1]
      args <- args[-c(ue_idx, ue_idx + 1)]
    }
  }
  if ("--failures" %in% args) {
    failures_only <- TRUE
    args <- args[args != "--failures"]
  }

  if (length(args) == 0) {
    gtf_file <- default_gtf
    pairs_file <- default_pairs
    output_pdf <- default_output
  } else if (length(args) == 3) {
    gtf_file <- args[1]
    pairs_file <- args[2]
    output_pdf <- args[3]
  } else {
    cat("Usage: Rscript scripts/tests/visualize_results.R [<gtf> <pairs> <output.pdf>] [--ue <ue.rds>] [--failures]\n")
    quit(status = 1)
  }

  visualize_tests(gtf_file, pairs_file, output_pdf, ue_rds_file, failures_only)
}
