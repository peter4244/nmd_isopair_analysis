#!/usr/bin/env Rscript
#
# Visualization Functions — Isoform Structure and Event Visualization
#
# Core library for plotting isoform pairs: dominant, comparator, reconstructed,
# with detected events annotated. Works with pipeline RDS outputs.
#
# Requires: ggplot2, patchwork, tidyverse
# Depends on: reconstruction_functions.R, event_detection_functions.R

# ==============================================================================
# Mismatch Detection
# ==============================================================================

#' Compare dominant and reconstructed exon sets, identify mismatches
#' @param dom_exons tibble with exon_start, exon_end (sorted)
#' @param recon_exons tibble with exon_start, exon_end (sorted)
#' @param strand "+" or "-"
#' @return tibble with exon_start, exon_end, mismatch (logical), detail
identify_mismatches <- function(dom_exons, recon_exons, strand) {
  dom_sorted <- dom_exons %>% arrange(exon_start)
  recon_sorted <- recon_exons %>% arrange(exon_start)

  n_dom <- nrow(dom_sorted)
  n_recon <- nrow(recon_sorted)

  mismatches <- list()

  if (n_dom != n_recon) {
    for (i in seq_len(n_recon)) {
      mismatches[[i]] <- tibble(
        exon_start = recon_sorted$exon_start[i],
        exon_end = recon_sorted$exon_end[i],
        mismatch = TRUE,
        detail = sprintf("Exon count: %d dom vs %d recon", n_dom, n_recon)
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
# Event Annotation Helpers
# ==============================================================================

#' Extract genomic coordinates from an event row for annotation
#' @param event one-row tibble from detect_events_for_pair()
#' @return tibble with x_start, x_end, label, direction for plotting
event_to_annotation <- function(event) {
  etype <- event$event_type

  # Parse coordinates depending on event type
  if (etype %in% c("IR", "IR_diff_5", "IR_diff_3", "IR_diff_5_3")) {
    x_start <- event$five_prime
    x_end <- event$three_prime
  } else if (etype %in% c("A5SS", "A3SS", "Partial_IR_5", "Partial_IR_3")) {
    x_start <- event$five_prime
    x_end <- event$three_prime
  } else if (etype %in% c("SE", "Missing_Internal")) {
    # Use five_prime/three_prime which bracket the event
    x_start <- event$five_prime
    x_end <- event$three_prime
  } else if (etype == "Alt_TSS") {
    x_start <- event$five_prime
    x_end <- event$three_prime
  } else if (etype == "Alt_TES") {
    x_start <- event$five_prime
    x_end <- event$three_prime
  } else {
    x_start <- event$five_prime
    x_end <- event$three_prime
  }

  tibble(
    x_start = as.numeric(x_start),
    x_end = as.numeric(x_end),
    label = etype,
    direction = as.character(event$direction)
  )
}

# ==============================================================================
# Single-Pair Plot
# ==============================================================================

#' Plot a single isoform comparison: dominant, comparator, reconstructed + events
#'
#' @param dom_exons tibble with exon_start, exon_end for dominant
#' @param comp_exons tibble with exon_start, exon_end for comparator
#' @param recon_exons tibble with exon_start, exon_end for reconstructed (or NULL to skip)
#' @param events tibble from detect_events_for_pair() (can be 0-row)
#' @param gene_id gene identifier string
#' @param dom_id dominant isoform ID
#' @param comp_id comparator isoform ID
#' @param strand "+" or "-"
#' @param status "PASS", "FAIL", or "ERROR"
#' @param reason status detail string
#' @param pair_label optional label (e.g., pair index or comparison name)
#' @return ggplot object
plot_isoform_pair <- function(dom_exons, comp_exons, recon_exons = NULL, events = NULL,
                               gene_id, dom_id, comp_id, strand,
                               status = NA_character_, reason = "",
                               pair_label = "") {

  show_recon <- !is.null(recon_exons) && nrow(recon_exons) > 0

  # Build track data
  dom_plot <- dom_exons %>%
    arrange(exon_start) %>%
    mutate(track = "Dominant", y_pos = if (show_recon) 3 else 2)

  comp_plot <- comp_exons %>%
    arrange(exon_start) %>%
    mutate(track = "Comparator", y_pos = 1)

  if (show_recon) {
    recon_plot <- recon_exons %>%
      arrange(exon_start) %>%
      mutate(track = "Reconstructed", y_pos = 2)
    all_exons <- bind_rows(dom_plot, recon_plot, comp_plot)
  } else {
    all_exons <- bind_rows(dom_plot, comp_plot)
  }

  x_min <- min(all_exons$exon_start)
  x_max <- max(all_exons$exon_end)
  x_range <- x_max - x_min
  x_buffer <- max(x_range * 0.05, 50)

  # Intron connector function
  make_introns <- function(track_data) {
    track_data %>%
      arrange(exon_start) %>%
      mutate(next_start = lead(exon_start)) %>%
      filter(!is.na(next_start))
  }

  dom_introns <- make_introns(dom_plot)
  comp_introns <- make_introns(comp_plot)

  p <- ggplot()

  # Intron lines
  if (nrow(dom_introns) > 0) {
    p <- p + geom_segment(data = dom_introns,
      aes(x = exon_end, xend = next_start, y = y_pos, yend = y_pos),
      color = "gray50", linewidth = 0.4)
  }
  if (nrow(comp_introns) > 0) {
    p <- p + geom_segment(data = comp_introns,
      aes(x = exon_end, xend = next_start, y = y_pos, yend = y_pos),
      color = "gray50", linewidth = 0.4)
  }

  if (show_recon) {
    recon_introns <- make_introns(recon_plot)
    if (nrow(recon_introns) > 0) {
      p <- p + geom_segment(data = recon_introns,
        aes(x = exon_end, xend = next_start, y = y_pos, yend = y_pos),
        color = "gray50", linewidth = 0.4)
    }
  }

  # Exon rectangles
  p <- p +
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

  # Mismatch highlights on reconstructed track
  if (show_recon) {
    mismatch_data <- identify_mismatches(dom_exons, recon_exons, strand)
    mismatched <- mismatch_data %>% filter(mismatch)
    if (nrow(mismatched) > 0) {
      p <- p +
        geom_rect(data = mismatched,
                  aes(xmin = exon_start, xmax = exon_end,
                      ymin = 2 - 0.35, ymax = 2 + 0.35),
                  fill = NA, color = "red", linewidth = 1.2) +
        geom_segment(data = mismatched,
                     aes(x = (exon_start + exon_end) / 2,
                         xend = (exon_start + exon_end) / 2,
                         y = 2.35, yend = 2.65),
                     color = "red", linewidth = 0.6, linetype = "dashed")
    }
  }

  # Event annotations (brackets below the comparator track)
  if (!is.null(events) && nrow(events) > 0) {
    event_annots <- map_dfr(seq_len(nrow(events)), function(i) {
      event_to_annotation(events[i, ])
    })
    # Filter to events with valid coordinates
    event_annots <- event_annots %>% filter(!is.na(x_start), !is.na(x_end))

    if (nrow(event_annots) > 0) {
      # Stagger vertically if events overlap
      event_annots <- event_annots %>%
        arrange(x_start) %>%
        mutate(row = row_number())

      y_base <- 0.4
      y_step <- 0.3

      for (i in seq_len(nrow(event_annots))) {
        ea <- event_annots[i, ]
        y_level <- y_base - (i - 1) * y_step

        dir_symbol <- if (ea$direction == "LOSS") "\u2193" else "\u2191"
        event_label <- sprintf("%s %s", ea$label, dir_symbol)
        bracket_color <- if (ea$direction == "LOSS") "#C0392B" else "#2980B9"

        p <- p +
          # Bracket
          geom_segment(aes(x = !!ea$x_start, xend = !!ea$x_end,
                           y = !!y_level, yend = !!y_level),
                       color = bracket_color, linewidth = 1.2) +
          # Left tick
          geom_segment(aes(x = !!ea$x_start, xend = !!ea$x_start,
                           y = !!y_level, yend = !!(y_level + 0.12)),
                       color = bracket_color, linewidth = 0.8) +
          # Right tick
          geom_segment(aes(x = !!ea$x_end, xend = !!ea$x_end,
                           y = !!y_level, yend = !!(y_level + 0.12)),
                       color = bracket_color, linewidth = 0.8) +
          # Label
          annotate("text", x = (ea$x_start + ea$x_end) / 2, y = y_level - 0.12,
                   label = event_label, size = 2.7, color = bracket_color, fontface = "bold")
      }
    }
  }

  # Track labels
  if (show_recon) {
    label_df <- tibble(
      y_pos = c(3, 2, 1),
      label = c(sprintf("Dominant\n(%s)", dom_id),
                "Reconstructed",
                sprintf("Comparator\n(%s)", comp_id))
    )
    y_limits <- c(-0.5, 3.7)
  } else {
    label_df <- tibble(
      y_pos = c(2, 1),
      label = c(sprintf("Dominant\n(%s)", dom_id),
                sprintf("Comparator\n(%s)", comp_id))
    )
    y_limits <- c(-0.5, 2.7)
  }

  # Adjust y_limits if many events push annotations down
  if (!is.null(events) && nrow(events) > 3) {
    y_limits[1] <- min(y_limits[1], 0.4 - nrow(events) * 0.3 - 0.3)
  }

  p <- p +
    geom_text(data = label_df,
              aes(x = x_min - x_buffer * 1.5, y = y_pos, label = label),
              hjust = 1, size = 3.5, fontface = "bold") +
    scale_y_continuous(breaks = NULL, limits = y_limits) +
    coord_cartesian(xlim = c(x_min - x_buffer * 3, x_max + x_buffer), clip = "off") +
    labs(x = sprintf("Genomic position (%s strand)", strand), y = NULL) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 10),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      plot.margin = margin(5, 10, 5, 100)
    )

  # Title
  event_str <- if (!is.null(events) && nrow(events) > 0) {
    paste(events$event_type, collapse = ", ")
  } else {
    "none"
  }

  title_text <- if (!is.na(status)) {
    sprintf("%s [%s] %s", pair_label, status, gene_id)
  } else {
    sprintf("%s %s", pair_label, gene_id)
  }
  title_color <- if (!is.na(status) && status == "PASS") "darkgreen" else
                 if (!is.na(status) && status %in% c("FAIL", "ERROR")) "darkred" else "black"

  subtitle_text <- sprintf("Events: %s", event_str)
  if (reason != "" && !is.na(reason)) subtitle_text <- paste0(subtitle_text, " | ", reason)

  p <- p +
    ggtitle(title_text, subtitle = subtitle_text) +
    theme(
      plot.title = element_text(size = 13, face = "bold", color = title_color),
      plot.subtitle = element_text(size = 10, color = "gray30")
    )

  p
}

# ==============================================================================
# Batch Visualization from Pipeline Data
# ==============================================================================

#' Visualize a set of pairs from pipeline RDS data
#'
#' Loads isoform structures, runs reconstruction, and plots each pair.
#'
#' @param profiles_rds path to splicing_choice_profiles.rds (with detailed_events)
#' @param structures_rds path to isoform_structures_filtered.rds (compact)
#' @param union_exons_rds path to union_exons_filtered.rds
#' @param output_pdf output PDF path
#' @param n_sample number of pairs to randomly sample (NULL = all)
#' @param seed random seed for reproducible sampling
#' @param failures_only only plot FAIL/ERROR pairs
#' @return invisible list of plot objects
visualize_profiles <- function(profiles_rds, structures_rds, union_exons_rds,
                                output_pdf, n_sample = 20, seed = 42,
                                failures_only = FALSE) {

  cat("\n======================================================================\n")
  cat("  ISOFORM PAIR VISUALIZATION\n")
  cat("======================================================================\n")
  cat(sprintf("  Profiles:  %s\n", profiles_rds))
  cat(sprintf("  Output:    %s\n", output_pdf))
  if (!is.null(n_sample)) cat(sprintf("  Sampling:  %d pairs (seed=%d)\n", n_sample, seed))
  cat("======================================================================\n\n")

  # Load data
  cat("Loading data...\n")
  profiles <- readRDS(profiles_rds)
  structures_compact <- readRDS(structures_rds)
  union_exons_raw <- readRDS(union_exons_rds)

  cat(sprintf("  %d profiles, %d structures\n", nrow(profiles), nrow(structures_compact)))

  # Expand structures
  structures <- structures_compact %>%
    rowwise() %>%
    mutate(exon_data = list(tibble(
      exon_number = seq_along(exon_starts),
      exon_start = exon_starts,
      exon_end = exon_ends
    ))) %>%
    unnest(exon_data) %>%
    select(isoform_id, gene_id, seqnames, strand, exon_number, exon_start, exon_end) %>%
    ungroup()

  # Prepare union exons for reconstruction
  union_exons_recon <- union_exons_raw %>%
    rename(start = union_exon_start, end = union_exon_end, chr = seqnames)

  # Sample pairs
  if (!is.null(n_sample) && n_sample < nrow(profiles)) {
    set.seed(seed)
    idx <- sample(seq_len(nrow(profiles)), n_sample)
    profiles_subset <- profiles[idx, ]
  } else {
    profiles_subset <- profiles
  }

  cat(sprintf("  Processing %d pairs...\n\n", nrow(profiles_subset)))

  # Process each pair
  plots <- list()
  plot_count <- 0

  for (i in seq_len(nrow(profiles_subset))) {
    prof <- profiles_subset[i, ]
    gene <- prof$gene_id
    dom_id <- prof$dominant_isoform_id
    comp_id <- prof$non_dominant_isoform_id

    # Get exon structures
    dom_exons <- structures %>%
      filter(isoform_id == dom_id) %>%
      arrange(exon_start)
    comp_exons <- structures %>%
      filter(isoform_id == comp_id) %>%
      arrange(exon_start)

    if (nrow(dom_exons) == 0 || nrow(comp_exons) == 0) {
      cat(sprintf("  [%d] %s: SKIP (missing structures)\n", i, gene))
      next
    }

    gene_strand <- unique(dom_exons$strand)[1]

    # Extract events
    events <- prof$detailed_events[[1]]

    # Reconstruct
    status <- "PASS"
    reason <- "Match"
    recon_exons <- NULL

    comp_for_recon <- comp_exons %>%
      rename(chr = seqnames, transcript_id = isoform_id) %>%
      select(chr, exon_start, exon_end, strand, gene_id, transcript_id)

    if (is.null(events) || nrow(events) == 0) {
      # No events: comparator should equal dominant
      recon_exons <- comp_for_recon
      vr <- tryCatch(
        verify_transcript(dom_exons, recon_exons, gene_strand),
        error = function(e) list(pass = FALSE, reason = e$message)
      )
      status <- if (vr$pass) "PASS" else "FAIL"
      reason <- vr$reason
    } else {
      gene_ue <- union_exons_recon %>% filter(gene_id == gene)
      recon_exons <- tryCatch({
        r <- reconstruct_dominant_v2(comp_for_recon, events, gene_ue)
        if (nrow(r) == 0) {
          status <<- "ERROR"
          reason <<- "Reconstruction produced 0 exons"
          comp_for_recon
        } else {
          vr <- verify_transcript(dom_exons, r, gene_strand)
          status <<- if (vr$pass) "PASS" else "FAIL"
          reason <<- vr$reason
          r
        }
      }, error = function(e) {
        status <<- "ERROR"
        reason <<- e$message
        comp_for_recon
      })
    }

    if (failures_only && status == "PASS") next

    cat(sprintf("  [%d] %s: %s — %s vs %s (%d events)\n",
                i, status, gene, dom_id, comp_id, nrow(events %||% tibble())))

    p <- tryCatch(
      plot_isoform_pair(
        dom_exons = dom_exons,
        comp_exons = comp_exons,
        recon_exons = recon_exons,
        events = events,
        gene_id = gene,
        dom_id = dom_id,
        comp_id = comp_id,
        strand = gene_strand,
        status = status,
        reason = reason,
        pair_label = sprintf("#%d", i)
      ),
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

  # Save PDF
  cat(sprintf("\nSaving %d plots to %s...\n", plot_count, output_pdf))

  # Ensure output directory exists
  out_dir <- dirname(output_pdf)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  combined <- wrap_plots(plots, ncol = 1) +
    plot_annotation(
      title = "Isoform Pair Visualization",
      subtitle = sprintf("%d pairs | %s", plot_count, Sys.Date()),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12, color = "gray40")
      )
    )

  ggsave(output_pdf, combined,
         width = 14, height = plot_count * 4,
         limitsize = FALSE)

  cat(sprintf("Saved: %s\n\n", output_pdf))
  invisible(plots)
}
