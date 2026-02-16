#!/usr/bin/env Rscript
#
# Event Detection Functions Library
# Shared by Script 07 (production) and validation scripts
#
# This file contains all event detection logic to ensure consistency
# between production analysis and validation testing.
#

# ==============================================================================
# Detection Thresholds (Constants)
# ==============================================================================

# Terminal boundary change detection
TSS_TOLERANCE <- 20  # bp tolerance for TSS change detection
TES_TOLERANCE <- 20  # bp tolerance for TES change detection

# Splice site classification boundary
# Events with coordinate differences <SPLICE_SITE_THRESHOLD are classified as:
#   - A5SS (alternative 5' splice site) if donor differs
#   - A3SS (alternative 3' splice site) if acceptor differs
# Events with differences ≥SPLICE_SITE_THRESHOLD are classified as:
#   - Partial_IR (partial intron retention) if one boundary shared
#   - IR (intron retention) if spanning multiple exons
SPLICE_SITE_THRESHOLD <- 100  # bp threshold: <100 = splice site, ≥100 = retention

# Exon overlap detection
# Two exons are considered "overlapping" (same exonic region) if they overlap
# by at least OVERLAP_THRESHOLD of the shorter exon
OVERLAP_THRESHOLD <- 0.5  # proportion (0.5 = 50% overlap required)

# ==============================================================================
# Biological Exon Ordering (Canonical Function)
# ==============================================================================

#' Order exons in biological transcription order (5' to 3')
#'
#' GROUND TRUTH: Strand and genomic coordinates determine biological order.
#' - Plus strand: Transcription proceeds 5'→3' from LOW to HIGH coordinates
#'   → Exon 1 (TSS) has lowest start, last exon (TES) has highest end
#' - Minus strand: Transcription proceeds 5'→3' from HIGH to LOW coordinates
#'   → Exon 1 (TSS) has highest end, last exon (TES) has lowest start
#'
#' This function establishes the canonical ordering used throughout all
#' event detection and validation code.
#'
#' @param exons Tibble/data.frame with exon_start and exon_end columns
#' @param strand Character: "+" or "-"
#' @return Exons sorted in biological order (TSS → TES) with biological_exon_number
order_exons_biological <- function(exons, strand) {
  if (strand == "+") {
    # Plus strand: sort by start position (ascending)
    # TSS has lowest coordinate
    exons <- exons %>%
      arrange(exon_start) %>%
      mutate(biological_exon_number = row_number())
  } else if (strand == "-") {
    # Minus strand: sort by start position (descending)
    # TSS has highest coordinate
    exons <- exons %>%
      arrange(desc(exon_start)) %>%
      mutate(biological_exon_number = row_number())
  } else {
    stop("Invalid strand: must be '+' or '-'")
  }

  return(exons)
}

# ==============================================================================
# Event Detection Functions
# ==============================================================================

#' Calculate overlap between two exons
#'
#' @param exon1 First exon (tibble row with exon_start, exon_end)
#' @param exon2 Second exon (tibble row with exon_start, exon_end)
#' @return Logical - TRUE if overlap ≥ 50% of shorter exon
calculate_overlap <- function(exon1, exon2) {
  # Calculate overlap region
  overlap_start <- max(exon1$exon_start, exon2$exon_start)
  overlap_end <- min(exon1$exon_end, exon2$exon_end)
  overlap_length <- max(0, overlap_end - overlap_start + 1)

  # Calculate exon lengths
  len1 <- exon1$exon_end - exon1$exon_start + 1
  len2 <- exon2$exon_end - exon2$exon_start + 1

  # Overlap percentage relative to shorter exon
  overlap_pct <- overlap_length / min(len1, len2)

  return(overlap_pct >= OVERLAP_THRESHOLD)
}

#' Detect TSS changes between first exons
#'
#' @param exon_dom First exon from dominant isoform (tibble row)
#' @param exon_non_dom First exon from non-dominant isoform (tibble row)
#' @param strand Gene strand
#' @param tolerance TSS tolerance in bp
#' @return Logical - TRUE if TSS changed
detect_tss_change <- function(exon_dom, exon_non_dom, strand, tolerance = TSS_TOLERANCE) {
  # Determine TSS coordinates based on strand
  if (strand == "+") {
    tss_dom <- exon_dom$exon_start
    tss_non_dom <- exon_non_dom$exon_start
  } else {
    tss_dom <- exon_dom$exon_end
    tss_non_dom <- exon_non_dom$exon_end
  }

  tss_diff <- abs(tss_dom - tss_non_dom)
  return(tss_diff > tolerance)
}

#' Detect TES changes between last exons
#'
#' @param exon_dom Last exon from dominant isoform (tibble row)
#' @param exon_non_dom Last exon from non-dominant isoform (tibble row)
#' @param strand Gene strand
#' @param tolerance TES tolerance in bp
#' @return Logical - TRUE if TES changed
detect_tes_change <- function(exon_dom, exon_non_dom, strand, tolerance = TES_TOLERANCE) {
  # Determine TES coordinates based on strand
  if (strand == "+") {
    tes_dom <- exon_dom$exon_end
    tes_non_dom <- exon_non_dom$exon_end
  } else {
    tes_dom <- exon_dom$exon_start
    tes_non_dom <- exon_non_dom$exon_start
  }

  tes_diff <- abs(tes_dom - tes_non_dom)
  return(tes_diff > tolerance)
}

#' Unified detection for shared-boundary and dual-boundary events
#'
#' Detects splicing events based on boundary comparisons:
#' - ONE boundary shared: A5SS, A3SS, or Partial_IR (by size)
#' - BOTH boundaries differ:
#'   - Terminal exons with ≥50% overlap: Partial_IR if internal boundary ≥100bp
#'   - Internal exons: Dual_boundary (wastebin category)
#'
#' Biological rationale for terminal exons:
#' - If terminal exons overlap ≥50%, spliceosome had access to both sites
#' - Internal boundary differences represent real splicing choices
#' - If overlap <50%, just Alt_TSS/Alt_TES (spliceosome never saw "old" site)
#'
#' @param exon_dom Exon from dominant isoform
#' @param exon_non_dom Exon from non-dominant isoform
#' @param strand Gene strand (+ or -)
#' @param is_first_exon Logical - is this a first exon comparison?
#' @param is_last_exon Logical - is this a last exon comparison?
#' @param terminal_has_overlap Logical - do terminal exons overlap ≥50%? (only relevant for terminal exons)
#' @param flanking_exons_dom Other exons in dominant isoform (for IR check)
#' @param flanking_exons_non_dom Other exons in non-dominant isoform (for IR check)
#' @return List with event_type and bp_diff
detect_shared_boundary_event <- function(exon_dom, exon_non_dom, strand,
                                         is_first_exon = FALSE,
                                         is_last_exon = FALSE,
                                         terminal_has_overlap = FALSE,
                                         flanking_exons_dom = NULL,
                                         flanking_exons_non_dom = NULL) {
  # Determine which boundaries are shared and which differ
  if (strand == "+") {
    # Plus strand: start = acceptor, end = donor
    shares_acceptor <- (exon_dom$exon_start == exon_non_dom$exon_start)
    differs_donor <- (exon_dom$exon_end != exon_non_dom$exon_end)
    donor_diff <- if (differs_donor) abs(exon_dom$exon_end - exon_non_dom$exon_end) else 0

    shares_donor <- (exon_dom$exon_end == exon_non_dom$exon_end)
    differs_acceptor <- (exon_dom$exon_start != exon_non_dom$exon_start)
    acceptor_diff <- if (differs_acceptor) abs(exon_dom$exon_start - exon_non_dom$exon_start) else 0
  } else {
    # Minus strand: start = donor, end = acceptor (coordinates reversed)
    shares_acceptor <- (exon_dom$exon_end == exon_non_dom$exon_end)
    differs_donor <- (exon_dom$exon_start != exon_non_dom$exon_start)
    donor_diff <- if (differs_donor) abs(exon_dom$exon_start - exon_non_dom$exon_start) else 0

    shares_donor <- (exon_dom$exon_start == exon_non_dom$exon_start)
    differs_acceptor <- (exon_dom$exon_end != exon_non_dom$exon_end)
    acceptor_diff <- if (differs_acceptor) abs(exon_dom$exon_end - exon_non_dom$exon_end) else 0
  }

  # Classify based on which boundary is shared and size of difference
  event_type <- "none"
  bp_diff <- 0

  if (shares_acceptor && differs_donor) {
    # Acceptor shared, donor differs
    bp_diff <- donor_diff
    if (bp_diff < SPLICE_SITE_THRESHOLD) {
      event_type <- "A5SS"
    } else {
      # Check if differing donor is a terminal boundary (not an internal splice site)
      # Donor = TES on last exons (both strands)
      donor_is_terminal <- is_last_exon

      if (donor_is_terminal) {
        # Terminal boundaries should not trigger Partial_IR
        event_type <- "none"
      } else {
        # Partial_IR at 5' splice site (donor differs)
        if (!is.null(flanking_exons_dom) && !is.null(flanking_exons_non_dom)) {
          spans_flanking <- check_spans_flanking_exons(
            exon_dom, exon_non_dom,
            flanking_exons_dom, flanking_exons_non_dom
          )
          event_type <- if (spans_flanking) "none" else "Partial_IR_5"  # none = will be IR
        } else {
          event_type <- "Partial_IR_5"
        }
      }
    }
  } else if (shares_donor && differs_acceptor) {
    # Donor shared, acceptor differs
    bp_diff <- acceptor_diff
    if (bp_diff < SPLICE_SITE_THRESHOLD) {
      event_type <- "A3SS"
    } else {
      # Check if differing acceptor is a terminal boundary (not an internal splice site)
      # Acceptor = TSS on first exons (both strands)
      acceptor_is_terminal <- is_first_exon

      if (acceptor_is_terminal) {
        # Terminal boundaries should not trigger Partial_IR
        event_type <- "none"
      } else {
        # Partial_IR at 3' splice site (acceptor differs)
        if (!is.null(flanking_exons_dom) && !is.null(flanking_exons_non_dom)) {
          spans_flanking <- check_spans_flanking_exons(
            exon_dom, exon_non_dom,
            flanking_exons_dom, flanking_exons_non_dom
          )
          event_type <- if (spans_flanking) "none" else "Partial_IR_3"  # none = will be IR
        } else {
          event_type <- "Partial_IR_3"
        }
      }
    }
  } else if (differs_acceptor && differs_donor) {
    # BOTH boundaries differ - dual-shift case
    # Treatment depends on whether this is a terminal exon with overlap

    if ((is_first_exon || is_last_exon) && terminal_has_overlap) {
      # Terminal exon with ≥50% overlap: check internal boundary
      # Biological rationale: spliceosome had access to both sites

      if (is_first_exon) {
        # First exon: internal boundary is the DONOR
        # (TSS is terminal, already handled by Alt_TSS)
        if (donor_diff >= SPLICE_SITE_THRESHOLD) {
          event_type <- "Partial_IR_5"
          bp_diff <- donor_diff
        }
      } else if (is_last_exon) {
        # Last exon: internal boundary is the ACCEPTOR
        # (TES is terminal, already handled by Alt_TES)
        if (acceptor_diff >= SPLICE_SITE_THRESHOLD) {
          event_type <- "Partial_IR_3"
          bp_diff <- acceptor_diff
        }
      }
      # If internal boundary <100bp, return "none" (not reported)

    } else if (!is_first_exon && !is_last_exon) {
      # Internal exon: dual-shift is "Dual_boundary" (wastebin category)
      event_type <- "Dual_boundary"
      bp_diff <- max(donor_diff, acceptor_diff)  # Report largest shift
    }
    # Terminal exons with <50% overlap: return "none" (not biologically meaningful)
  }

  return(list(event_type = event_type, bp_diff = bp_diff))
}

#' Check if Partial_IR extension spans into flanking exons
#'
#' If a Partial_IR event extends into a flanking exon in the other isoform,
#' it should be classified as IR (intron retention) instead.
#'
#' @param exon_dom Exon from dominant isoform
#' @param exon_non_dom Exon from non-dominant isoform
#' @param flanking_exons_dom Other exons in dominant isoform
#' @param flanking_exons_non_dom Other exons in non-dominant isoform
#' @return TRUE if extension overlaps flanking exons, FALSE otherwise
check_spans_flanking_exons <- function(exon_dom, exon_non_dom,
                                       flanking_exons_dom, flanking_exons_non_dom) {
  # Determine which exon is longer and check if extension overlaps flanking exons

  # Check if non_dom extension overlaps with flanking exons in dom
  if (nrow(flanking_exons_dom) > 0) {
    for (i in seq_len(nrow(flanking_exons_dom))) {
      flank <- flanking_exons_dom[i, ]
      # Check if non_dom exon overlaps this flanking exon
      overlaps <- (exon_non_dom$exon_start <= flank$exon_end) &&
                  (exon_non_dom$exon_end >= flank$exon_start)
      if (overlaps) return(TRUE)
    }
  }

  # Check if dom extension overlaps with flanking exons in non_dom
  if (nrow(flanking_exons_non_dom) > 0) {
    for (i in seq_len(nrow(flanking_exons_non_dom))) {
      flank <- flanking_exons_non_dom[i, ]
      # Check if dom exon overlaps this flanking exon
      overlaps <- (exon_dom$exon_start <= flank$exon_end) &&
                  (exon_dom$exon_end >= flank$exon_start)
      if (overlaps) return(TRUE)
    }
  }

  return(FALSE)
}

#' Detect IR (intron retention) by checking if exon spans multiple exons
#' This is a simplified version - full IR detection would need more context
#'
#' @param exon Exon to check
#' @param other_exons Other exons to check against
#' @return Logical - TRUE if exon spans 2+ other exons
detect_ir_simple <- function(exon, other_exons) {
  # Check if this exon spans multiple consecutive exons in the other isoform
  overlapping <- other_exons %>%
    filter(
      (exon_start >= exon$exon_start & exon_start <= exon$exon_end) |
      (exon_end >= exon$exon_start & exon_end <= exon$exon_end) |
      (exon_start <= exon$exon_start & exon_end >= exon$exon_end)
    )

  return(nrow(overlapping) >= 2)
}

#' Detect SE (skipped exon) events
#' SE = exon present in only one isoform, flanked by shared/comparable exons
#'
#' @param comparison Comparison dataframe with union exons and exon_status
#' @return Integer count of SE events
detect_se <- function(comparison) {
  n_se <- 0

  for (i in seq_len(nrow(comparison))) {
    exon_status <- comparison[i, ]$exon_status

    # Skip if this exon is shared or in neither isoform
    if (exon_status %in% c("shared", "neither")) next

    # Check if flanking exons are shared (comparable)
    has_prev <- i > 1
    has_next <- i < nrow(comparison)

    prev_shared <- has_prev && comparison[i-1, ]$exon_status == "shared"
    next_shared <- has_next && comparison[i+1, ]$exon_status == "shared"

    # SE if both flanking exons are shared/comparable
    if (prev_shared && next_shared) {
      n_se <- n_se + 1
    }
  }

  return(n_se)
}
