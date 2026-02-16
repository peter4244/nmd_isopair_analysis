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
#' @param is_first_exon Logical - are BOTH exons first exons? (for shared boundary terminal checks)
#' @param is_last_exon Logical - are BOTH exons last exons? (for shared boundary terminal checks)
#' @param is_first_exon_comp Logical - is comparison exon a first exon? (for overlap-based TSS check)
#' @param is_last_exon_comp Logical - is comparison exon a last exon? (for overlap-based TES check)
#' @param terminal_has_overlap Logical - do terminal exons overlap ≥50%? (only relevant for terminal exons)
#' @param flanking_exons_dom Other exons in dominant isoform (for IR check)
#' @param flanking_exons_non_dom Other exons in non-dominant isoform (for IR check)
#' @return List with event_type, bp_diff, and direction (GAIN/LOSS)
detect_shared_boundary_event <- function(exon_dom, exon_non_dom, strand,
                                         is_first_exon = FALSE,
                                         is_last_exon = FALSE,
                                         is_first_exon_comp = FALSE,
                                         is_last_exon_comp = FALSE,
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
  direction <- NULL  # GAIN/LOSS for overlap-based detection

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

  # ===========================================================================
  # OVERLAP-BASED DETECTION (fallback when no exact shared boundary)
  # ===========================================================================
  # If no event detected via shared boundary, check if exons overlap
  # This allows detection of splice site variations even when boundaries don't match exactly

  if (event_type == "none") {
    # Check if exons overlap by at least 1bp
    exons_overlap <- (exon_dom$exon_start <= exon_non_dom$exon_end) &&
                     (exon_dom$exon_end >= exon_non_dom$exon_start)

    # Skip overlap-based detection for monoexonic vs multi-exonic comparisons
    # (Let IR detection handle those cases)
    has_flanking_dom <- !is.null(flanking_exons_dom) && nrow(flanking_exons_dom) > 0
    has_flanking_non_dom <- !is.null(flanking_exons_non_dom) && nrow(flanking_exons_non_dom) > 0

    # Only apply overlap-based detection if BOTH isoforms are multi-exonic
    # OR if we don't have flanking exon information (can't determine)
    skip_overlap <- (is.null(flanking_exons_dom) || is.null(flanking_exons_non_dom)) ||
                    (!has_flanking_dom && has_flanking_non_dom) ||
                    (has_flanking_dom && !has_flanking_non_dom)

    if (exons_overlap && !skip_overlap) {
      # Get 5' and 3' coordinates (strand-aware)
      if (strand == "+") {
        # Plus strand: 5' = start, 3' = end
        dom_5prime <- exon_dom$exon_start
        dom_3prime <- exon_dom$exon_end
        comp_5prime <- exon_non_dom$exon_start
        comp_3prime <- exon_non_dom$exon_end
      } else {
        # Minus strand: 5' = end, 3' = start
        dom_5prime <- exon_dom$exon_end
        dom_3prime <- exon_dom$exon_start
        comp_5prime <- exon_non_dom$exon_end
        comp_3prime <- exon_non_dom$exon_start
      }

      # Initialize results for overlap-based detection
      events_detected <- list()
      direction <- NULL

      # ===========================================================================
      # Check 5' end (skip if comparison exon is first - its 5' end is TSS, terminal)
      # ===========================================================================
      if (!is_first_exon_comp) {
        result_5prime <- check_boundary_within_exon(
          comp_boundary = comp_5prime,
          dom_boundary = dom_5prime,
          dom_exon_start = if (strand == "+") dom_5prime else dom_3prime,
          dom_exon_end = if (strand == "+") dom_3prime else dom_5prime,
          flanking_exons = if (!is.null(flanking_exons_dom)) {
            # Get upstream flanking exon (5' direction)
            if (strand == "+") {
              # Plus: upstream = lower coordinates
              flanking_exons_dom %>% filter(exon_end < dom_5prime) %>%
                arrange(desc(exon_end)) %>% slice(1)
            } else {
              # Minus: upstream = higher coordinates
              flanking_exons_dom %>% filter(exon_start > dom_5prime) %>%
                arrange(exon_start) %>% slice(1)
            }
          } else NULL,
          strand = strand,
          is_5prime = TRUE
        )

        if (result_5prime$event_type != "none") {
          events_detected[[length(events_detected) + 1]] <- result_5prime
        }
      }

      # ===========================================================================
      # Check 3' end (skip if comparison exon is last - its 3' end is TES, terminal)
      # ===========================================================================
      if (!is_last_exon_comp) {
        result_3prime <- check_boundary_within_exon(
          comp_boundary = comp_3prime,
          dom_boundary = dom_3prime,
          dom_exon_start = if (strand == "+") dom_5prime else dom_3prime,
          dom_exon_end = if (strand == "+") dom_3prime else dom_5prime,
          flanking_exons = if (!is.null(flanking_exons_dom)) {
            # Get downstream flanking exon (3' direction)
            if (strand == "+") {
              # Plus: downstream = higher coordinates
              flanking_exons_dom %>% filter(exon_start > dom_3prime) %>%
                arrange(exon_start) %>% slice(1)
            } else {
              # Minus: downstream = lower coordinates
              flanking_exons_dom %>% filter(exon_end < dom_3prime) %>%
                arrange(desc(exon_end)) %>% slice(1)
            }
          } else NULL,
          strand = strand,
          is_5prime = FALSE
        )

        if (result_3prime$event_type != "none") {
          events_detected[[length(events_detected) + 1]] <- result_3prime
        }
      }

      # Return first detected event (prioritize 5' over 3' if both)
      if (length(events_detected) > 0) {
        event_type <- events_detected[[1]]$event_type
        bp_diff <- events_detected[[1]]$bp_diff
        direction <- events_detected[[1]]$direction
      }
    }
  }

  return(list(event_type = event_type, bp_diff = bp_diff, direction = direction))
}

#' Check if a boundary from comparison exon triggers an event relative to dominant exon
#'
#' Implements the overlap-based detection logic for boundaries that don't match exactly
#'
#' @param comp_boundary The boundary coordinate from comparison exon to check
#' @param dom_boundary The corresponding boundary coordinate from dominant exon
#' @param dom_exon_start Start coordinate of dominant exon (genomic)
#' @param dom_exon_end End coordinate of dominant exon (genomic)
#' @param flanking_exons Flanking exon(s) in dominant isoform (upstream for 5', downstream for 3')
#' @param strand Strand (+ or -)
#' @param is_5prime Logical - checking 5' end (TRUE) or 3' end (FALSE)
#' @return List with event_type, bp_diff, and direction
check_boundary_within_exon <- function(comp_boundary, dom_boundary,
                                       dom_exon_start, dom_exon_end,
                                       flanking_exons, strand, is_5prime) {
  event_type <- "none"
  bp_diff <- 0
  direction <- NULL

  # Determine if comp_boundary is within dom exon
  within_dom <- (comp_boundary >= dom_exon_start) && (comp_boundary <= dom_exon_end)

  if (within_dom) {
    # =========================================================================
    # CASE A: Comparison boundary is WITHIN dominant exon
    # =========================================================================
    # Comparison exon is shorter on this side → LOSS
    distance <- abs(comp_boundary - dom_boundary)

    # If distance is 0, boundaries are identical - no event
    if (distance == 0) {
      return(list(event_type = "none", bp_diff = 0, direction = NULL))
    }

    bp_diff <- distance

    if (distance < SPLICE_SITE_THRESHOLD) {
      # <100bp difference → Alternative splice site
      # Naming based on FUNCTIONAL splice site (donor=5', acceptor=3')
      # Plus: 5' positional = acceptor, 3' positional = donor
      # Minus: 5' positional = acceptor, 3' positional = donor (same mapping!)
      event_type <- if (is_5prime) "A3SS" else "A5SS"
    } else {
      # ≥100bp difference → Partial retention
      event_type <- if (is_5prime) "Partial_IR_3" else "Partial_IR_5"
    }
    direction <- "LOSS"

  } else {
    # =========================================================================
    # CASE B: Comparison boundary extends BEYOND dominant exon
    # =========================================================================
    # Comparison exon is longer on this side → check flanking exon

    if (!is.null(flanking_exons) && nrow(flanking_exons) > 0) {
      flanking <- flanking_exons[1, ]

      # Check if comp_boundary overlaps or extends past the flanking exon
      if (strand == "+") {
        if (is_5prime) {
          # 5' on plus: comp extends to lower coordinates (upstream)
          # Check if it overlaps flanking exon (which has lower coords)
          overlaps_flanking <- comp_boundary <= flanking$exon_end
        } else {
          # 3' on plus: comp extends to higher coordinates (downstream)
          # Check if it overlaps flanking exon (which has higher coords)
          overlaps_flanking <- comp_boundary >= flanking$exon_start
        }
      } else {
        if (is_5prime) {
          # 5' on minus: comp extends to higher coordinates (upstream)
          # Check if it overlaps flanking exon (which has higher coords)
          overlaps_flanking <- comp_boundary >= flanking$exon_start
        } else {
          # 3' on minus: comp extends to lower coordinates (downstream)
          # Check if it overlaps flanking exon (which has lower coords)
          overlaps_flanking <- comp_boundary <= flanking$exon_end
        }
      }

      if (overlaps_flanking) {
        # Sub-case B1: Extends past flanking exon → IR
        event_type <- "IR"
        bp_diff <- abs(comp_boundary - dom_boundary)
        direction <- "GAIN"
      } else {
        # Sub-case B2: Extends beyond dom but not to flanking
        # Calculate distance between boundaries
        distance <- abs(comp_boundary - dom_boundary)
        bp_diff <- distance

        if (distance < SPLICE_SITE_THRESHOLD) {
          # Positional 5' = acceptor (A3SS), Positional 3' = donor (A5SS)
          event_type <- if (is_5prime) "A3SS" else "A5SS"
        } else {
          event_type <- if (is_5prime) "Partial_IR_3" else "Partial_IR_5"
        }
        direction <- "GAIN"
      }
    } else {
      # No flanking exon data - just compare boundaries
      distance <- abs(comp_boundary - dom_boundary)
      bp_diff <- distance

      if (distance < SPLICE_SITE_THRESHOLD) {
        # Positional 5' = acceptor (A3SS), Positional 3' = donor (A5SS)
        event_type <- if (is_5prime) "A3SS" else "A5SS"
      } else {
        event_type <- if (is_5prime) "Partial_IR_3" else "Partial_IR_5"
      }
      direction <- "GAIN"
    }
  }

  return(list(event_type = event_type, bp_diff = bp_diff, direction = direction))
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
