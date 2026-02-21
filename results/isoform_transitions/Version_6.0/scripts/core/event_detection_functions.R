#!/usr/bin/env Rscript
#
# Event Detection Functions Library
# Shared by Script 07 (production) and validation scripts
#
# This file contains all event detection logic to ensure consistency
# between production analysis and validation testing.
#

# ==============================================================================
# CRITICAL: Event Direction Semantics (GAIN vs LOSS)
# ==============================================================================
#
# IMPORTANT: GAIN and LOSS are defined from the COMPARATOR'S PERSPECTIVE
#
# When comparing: Dominant isoform vs Comparator isoform
#
# LOSS = Comparator LOST sequence (dominant has MORE)
#   - Reconstruction: ADD regions to comparator to build dominant
#   - Example: Dominant has extra 5' exon → comparator LOST that exon
#
# GAIN = Comparator GAINED sequence (dominant has LESS)
#   - Reconstruction: REMOVE regions from comparator to build dominant
#   - Example: Comparator has extra 3' exon → comparator GAINED that exon
#
# Alt_TSS/Alt_TES Examples:
#   Plus strand:
#     - dom_tss < comp_tss → dominant starts earlier → comparator LOST 5' → LOSS
#     - dom_tes > comp_tes → dominant ends later → comparator LOST 3' → LOSS
#   Minus strand:
#     - dom_tss > comp_tss → dominant starts earlier → comparator LOST 5' → LOSS
#     - dom_tes < comp_tes → dominant ends later → comparator LOST 3' → LOSS
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
  if (tss_diff > tolerance) return(TRUE)

  # Even if within coordinate tolerance, detect structural change when first
  # exons don't overlap — indicates an extra terminal exon, not just a shift
  first_overlap <- exon_dom$exon_start <= exon_non_dom$exon_end &&
                   exon_dom$exon_end >= exon_non_dom$exon_start
  if (!first_overlap) return(TRUE)

  return(FALSE)
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
  if (tes_diff > tolerance) return(TRUE)

  # Even if within coordinate tolerance, detect structural change when last
  # exons don't overlap — indicates an extra terminal exon, not just a shift
  last_overlap <- exon_dom$exon_start <= exon_non_dom$exon_end &&
                  exon_dom$exon_end >= exon_non_dom$exon_start
  if (!last_overlap) return(TRUE)

  return(FALSE)
}

#' Compute missing terminal exons for Alt_TSS events
#'
#' When dominant has longer 5' end (LOSS direction), computes the exonic regions
#' from dominant that are missing in comparator. Returns only exonic regions
#' (not introns) as comma-separated coordinate ranges.
#'
#' @param dom_exons Dominant exons (ordered TSS → TES)
#' @param comp_exons Comparator exons (ordered TSS → TES)
#' @param strand Gene strand
#' @return String: comma-separated ranges like "1000-1200,1500-1700,2000-2150"
compute_missing_terminal_exons_tss <- function(dom_exons, comp_exons, strand) {
  # Get TSS positions
  if (strand == "+") {
    dom_tss <- dom_exons$exon_start[1]
    comp_tss <- comp_exons$exon_start[1]
    # Dominant extends further 5' (lower coordinates)
    if (dom_tss >= comp_tss) return("")  # No missing region
    cutoff <- comp_tss  # Everything before this is missing
  } else {
    dom_tss <- dom_exons$exon_end[1]
    comp_tss <- comp_exons$exon_end[1]
    # Dominant extends further 5' (higher coordinates)
    if (dom_tss <= comp_tss) return("")  # No missing region
    cutoff <- comp_tss  # Everything after this is missing
  }

  # Walk through dominant exons from TSS, collect missing exonic regions
  missing_ranges <- list()

  for (i in seq_len(nrow(dom_exons))) {
    exon <- dom_exons[i, ]

    if (strand == "+") {
      # Plus: walk from low to high coordinates
      if (exon$exon_end < cutoff) {
        # Entire exon is missing
        missing_ranges[[length(missing_ranges) + 1]] <- sprintf("%d-%d", exon$exon_start, exon$exon_end)
      } else if (exon$exon_start < cutoff && exon$exon_end >= cutoff) {
        # Partial exon missing (5' portion)
        missing_ranges[[length(missing_ranges) + 1]] <- sprintf("%d-%d", exon$exon_start, cutoff - 1)
        break  # Reached comparator TSS
      } else {
        # exon_start >= cutoff, no more missing exons
        break
      }
    } else {
      # Minus: walk from high to low coordinates
      if (exon$exon_start > cutoff) {
        # Entire exon is missing
        missing_ranges[[length(missing_ranges) + 1]] <- sprintf("%d-%d", exon$exon_start, exon$exon_end)
      } else if (exon$exon_start <= cutoff && exon$exon_end > cutoff) {
        # Partial exon missing (5' portion)
        missing_ranges[[length(missing_ranges) + 1]] <- sprintf("%d-%d", cutoff + 1, exon$exon_end)
        break  # Reached comparator TSS
      } else {
        # exon_end <= cutoff, no more missing exons
        break
      }
    }
  }

  # Return comma-separated ranges
  if (length(missing_ranges) > 0) {
    return(paste(missing_ranges, collapse = ","))
  } else {
    return("")
  }
}

#' Compute missing terminal exons for Alt_TES events
#'
#' When dominant has longer 3' end (LOSS direction), computes the exonic regions
#' from dominant that are missing in comparator. Returns only exonic regions
#' (not introns) as comma-separated coordinate ranges.
#'
#' @param dom_exons Dominant exons (ordered TSS → TES)
#' @param comp_exons Comparator exons (ordered TSS → TES)
#' @param strand Gene strand
#' @return String: comma-separated ranges like "2150-2300,2500-2700,3000-3200"
compute_missing_terminal_exons_tes <- function(dom_exons, comp_exons, strand) {
  # Get TES positions
  last_dom_idx <- nrow(dom_exons)
  last_comp_idx <- nrow(comp_exons)

  if (strand == "+") {
    dom_tes <- dom_exons$exon_end[last_dom_idx]
    comp_tes <- comp_exons$exon_end[last_comp_idx]
    # Dominant extends further 3' (higher coordinates)
    if (dom_tes <= comp_tes) return("")  # No missing region
    cutoff <- comp_tes  # Everything after this is missing
  } else {
    dom_tes <- dom_exons$exon_start[last_dom_idx]
    comp_tes <- comp_exons$exon_start[last_comp_idx]
    # Dominant extends further 3' (lower coordinates)
    if (dom_tes >= comp_tes) return("")  # No missing region
    cutoff <- comp_tes  # Everything before this is missing
  }

  # Walk through dominant exons from TES (backwards), collect missing exonic regions
  missing_ranges <- list()

  for (i in last_dom_idx:1) {
    exon <- dom_exons[i, ]

    if (strand == "+") {
      # Plus: walk from high to low coordinates (backwards from TES)
      if (exon$exon_start > cutoff) {
        # Entire exon is missing
        missing_ranges[[length(missing_ranges) + 1]] <- sprintf("%d-%d", exon$exon_start, exon$exon_end)
      } else if (exon$exon_start <= cutoff && exon$exon_end > cutoff) {
        # Partial exon missing (3' portion)
        missing_ranges[[length(missing_ranges) + 1]] <- sprintf("%d-%d", cutoff + 1, exon$exon_end)
        break  # Reached comparator TES
      } else {
        # exon_end <= cutoff, no more missing exons
        break
      }
    } else {
      # Minus: walk from low to high coordinates (backwards from TES)
      if (exon$exon_end < cutoff) {
        # Entire exon is missing
        missing_ranges[[length(missing_ranges) + 1]] <- sprintf("%d-%d", exon$exon_start, exon$exon_end)
      } else if (exon$exon_end >= cutoff && exon$exon_start < cutoff) {
        # Partial exon missing (3' portion)
        missing_ranges[[length(missing_ranges) + 1]] <- sprintf("%d-%d", exon$exon_start, cutoff - 1)
        break  # Reached comparator TES
      } else {
        # exon_start >= cutoff, no more missing exons
        break
      }
    }
  }

  # Reverse the list (we walked backwards) and return comma-separated ranges
  if (length(missing_ranges) > 0) {
    missing_ranges <- rev(missing_ranges)
    return(paste(missing_ranges, collapse = ","))
  } else {
    return("")
  }
}

#' Compute orphaned terminal exons at the TSS end
#'
#' For Alt_TSS events, the "other" isoform may have terminal exons that sit
#' within the reference isoform's span but don't overlap any reference exon.
#' These are "orphaned" exons that must be removed (LOSS) or added (GAIN)
#' during reconstruction.
#'
#' Walk walk_ordered from its biological first exon (TSS) toward TES.
#' Collect each exon that has no overlap with any exon in ref_ordered.
#' Stop as soon as an overlapping exon is found.
#'
#' @param walk_ordered Exons to walk (ordered TSS → TES)
#' @param ref_ordered  Reference exons to test overlap against (ordered TSS → TES)
#' @return String: comma-separated "start-end" ranges, or "" if none
compute_orphan_terminal_exons_tss <- function(walk_ordered, ref_ordered) {
  orphan_ranges <- list()
  for (i in seq_len(nrow(walk_ordered))) {
    exon <- walk_ordered[i, ]
    has_overlap <- any(sapply(seq_len(nrow(ref_ordered)), function(j) {
      r <- ref_ordered[j, ]
      exon$exon_start <= r$exon_end && exon$exon_end >= r$exon_start
    }))
    if (has_overlap) break
    orphan_ranges[[length(orphan_ranges) + 1]] <-
      sprintf("%d-%d", exon$exon_start, exon$exon_end)
  }
  if (length(orphan_ranges) > 0) paste(orphan_ranges, collapse = ",") else ""
}

#' Compute orphaned terminal exons at the TES end
#'
#' For Alt_TES events, the "other" isoform may have terminal exons that sit
#' within the reference isoform's span but don't overlap any reference exon.
#' These are "orphaned" exons that must be removed (LOSS) or added (GAIN)
#' during reconstruction.
#'
#' Walk walk_ordered from its biological last exon (TES) toward TSS.
#' Collect each exon that has no overlap with any exon in ref_ordered.
#' Stop as soon as an overlapping exon is found.
#'
#' @param walk_ordered Exons to walk (ordered TSS → TES)
#' @param ref_ordered  Reference exons to test overlap against (ordered TSS → TES)
#' @return String: comma-separated "start-end" ranges in TSS→TES order, or "" if none
compute_orphan_terminal_exons_tes <- function(walk_ordered, ref_ordered) {
  orphan_ranges <- list()
  for (i in rev(seq_len(nrow(walk_ordered)))) {
    exon <- walk_ordered[i, ]
    has_overlap <- any(sapply(seq_len(nrow(ref_ordered)), function(j) {
      r <- ref_ordered[j, ]
      exon$exon_start <= r$exon_end && exon$exon_end >= r$exon_start
    }))
    if (has_overlap) break
    orphan_ranges[[length(orphan_ranges) + 1]] <-
      sprintf("%d-%d", exon$exon_start, exon$exon_end)
  }
  # Walked backwards; reverse to restore TSS→TES order
  if (length(orphan_ranges) > 0) paste(rev(orphan_ranges), collapse = ",") else ""
}

#' Unified detection for shared-boundary and dual-boundary events
#'
#' Detects splicing events based on boundary comparisons:
#' - ONE boundary shared: A5SS, A3SS, or Partial_IR (by size)
#' - BOTH boundaries differ:
#'   - Terminal exons with ≥50% overlap: Partial_IR if internal boundary ≥100bp
#'   - Internal exons: decomposed into two events (one per boundary)
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
#' @param is_first_exon_dom Logical - is the dominant exon a first exon? (TSS already handled by Alt_TSS)
#' @param is_last_exon_dom Logical - is the dominant exon a last exon? (TES already handled by Alt_TES)
#' @param is_first_exon_comp Logical - is comparison exon a first exon? (for overlap-based TSS check)
#' @param is_last_exon_comp Logical - is comparison exon a last exon? (for overlap-based TES check)
#' @param flanking_exons_dom Other exons in dominant isoform (for IR check)
#' @param flanking_exons_non_dom Other exons in non-dominant isoform (for IR check)
#' @return List with event_type, bp_diff, and direction (GAIN/LOSS)
detect_shared_boundary_event <- function(exon_dom, exon_non_dom, strand,
                                         is_first_exon = FALSE,
                                         is_last_exon = FALSE,
                                         is_first_exon_dom = FALSE,
                                         is_last_exon_dom = FALSE,
                                         is_first_exon_comp = FALSE,
                                         is_last_exon_comp = FALSE,
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
  second_event <- NULL  # Optional second event for asymmetric dual-boundary pairs

  if (shares_acceptor && differs_donor) {
    # Acceptor shared, donor differs
    bp_diff <- donor_diff

    # Compute direction: which isoform is longer?
    # Donor differs means one exon extends further at the donor boundary
    if (strand == "+") {
      # Plus: donor = exon end (higher coordinate = extends further)
      direction <- if (exon_dom$exon_end > exon_non_dom$exon_end) "LOSS" else "GAIN"
    } else {
      # Minus: donor = exon start (lower coordinate = extends further)
      direction <- if (exon_dom$exon_start < exon_non_dom$exon_start) "LOSS" else "GAIN"
    }

    # Check if differing donor is a terminal boundary (not an internal splice site)
    # Donor = TES on last exons (both strands).
    # Suppress ONLY when BOTH exons are last exons (both TES): in that case Alt_TES
    # already captures the boundary difference.
    # If only one exon is last (e.g., comparator's only exon vs dominant's non-last exon),
    # the non-terminal isoform's donor is a real splice site and must be detected.
    donor_is_terminal <- is_last_exon  # TRUE only when both dominant AND comparator are last

    if (donor_is_terminal) {
      # Both TES exons — Alt_TES handles this boundary, don't double-count
      event_type <- "none"
    } else if (bp_diff < SPLICE_SITE_THRESHOLD) {
      event_type <- "A5SS"
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
  } else if (shares_donor && differs_acceptor) {
    # Donor shared, acceptor differs
    bp_diff <- acceptor_diff

    # Compute direction: which isoform is longer?
    # Acceptor differs means one exon extends further at the acceptor boundary
    if (strand == "+") {
      # Plus: acceptor = exon start (lower coordinate = extends further)
      direction <- if (exon_dom$exon_start < exon_non_dom$exon_start) "LOSS" else "GAIN"
    } else {
      # Minus: acceptor = exon end (higher coordinate = extends further)
      direction <- if (exon_dom$exon_end > exon_non_dom$exon_end) "LOSS" else "GAIN"
    }

    # Check if differing acceptor is a terminal boundary (not an internal splice site)
    # Acceptor = TSS on first exons (both strands).
    # Suppress ONLY when BOTH exons are first exons (both TSS): in that case Alt_TSS
    # already captures the boundary difference.
    # If only one exon is first (e.g., comparator's only exon vs dominant's non-first exon),
    # the non-terminal isoform's acceptor is a real splice site and must be detected.
    acceptor_is_terminal <- is_first_exon  # TRUE only when both dominant AND comparator are first

    if (acceptor_is_terminal) {
      # Both TSS exons — Alt_TSS handles this boundary, don't double-count
      event_type <- "none"
    } else if (bp_diff < SPLICE_SITE_THRESHOLD) {
      event_type <- "A3SS"
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
  } else if (differs_acceptor && differs_donor) {
    # BOTH boundaries differ - dual-shift case
    # Treatment depends on whether any involved exon is terminal

    # A terminal exon from either isoform means one boundary is TSS/TES (already
    # captured by Alt_TSS/Alt_TES); only the internal-facing boundary is detected here.
    is_first_any <- is_first_exon || is_first_exon_dom || is_first_exon_comp
    is_last_any  <- is_last_exon  || is_last_exon_dom  || is_last_exon_comp

    if (is_first_any || is_last_any) {
      # A terminal exon is involved: primary = internal-facing boundary.
      # The outer (TSS/TES) boundary is captured by Alt_TSS/Alt_TES.
      # When asymmetric terminal (only ONE exon is terminal), the splice-site-facing
      # boundary is a real splice site — emit it as second_event.

      if (is_first_any) {
        # A first exon is involved: TSS boundary already handled by Alt_TSS.
        # Internal boundary is the DONOR.
        bp_diff <- donor_diff
        if (donor_diff < SPLICE_SITE_THRESHOLD) {
          event_type <- "A5SS"
        } else {
          event_type <- "Partial_IR_5"
        }
        # Compute direction based on donor (which exon extends further)
        if (strand == "+") {
          # Plus: donor = exon end (higher coordinate = extends further)
          direction <- if (exon_dom$exon_end > exon_non_dom$exon_end) "LOSS" else "GAIN"
        } else {
          # Minus: donor = exon start (lower coordinate = extends further)
          direction <- if (exon_dom$exon_start < exon_non_dom$exon_start) "LOSS" else "GAIN"
        }
        # Emit acceptor boundary as second event when asymmetric OR when both
        # first AND last exons are involved simultaneously.
        is_asymmetric_first <- xor(is_first_exon_dom, is_first_exon_comp)
        if (is_asymmetric_first || is_last_any) {
          acc_type <- if (acceptor_diff < SPLICE_SITE_THRESHOLD) "A3SS" else "Partial_IR_3"
          if (strand == "+") {
            acc_dir <- if (exon_dom$exon_start < exon_non_dom$exon_start) "LOSS" else "GAIN"
          } else {
            acc_dir <- if (exon_dom$exon_end > exon_non_dom$exon_end) "LOSS" else "GAIN"
          }
          second_event <- list(event_type = acc_type, bp_diff = acceptor_diff, direction = acc_dir)
        }

      } else if (is_last_any) {
        # A last exon is involved: TES boundary already handled by Alt_TES.
        # Internal boundary is the ACCEPTOR.
        bp_diff <- acceptor_diff
        if (acceptor_diff < SPLICE_SITE_THRESHOLD) {
          event_type <- "A3SS"
        } else {
          event_type <- "Partial_IR_3"
        }
        # Compute direction based on acceptor (which exon extends further)
        if (strand == "+") {
          # Plus: acceptor = exon start (lower coordinate = extends further)
          direction <- if (exon_dom$exon_start < exon_non_dom$exon_start) "LOSS" else "GAIN"
        } else {
          # Minus: acceptor = exon end (higher coordinate = extends further)
          direction <- if (exon_dom$exon_end > exon_non_dom$exon_end) "LOSS" else "GAIN"
        }
        # Emit donor boundary as second event when asymmetric (only one exon is last).
        is_asymmetric_last <- xor(is_last_exon_dom, is_last_exon_comp)
        if (is_asymmetric_last) {
          don_type <- if (donor_diff < SPLICE_SITE_THRESHOLD) "A5SS" else "Partial_IR_5"
          if (strand == "+") {
            don_dir <- if (exon_dom$exon_end > exon_non_dom$exon_end) "LOSS" else "GAIN"
          } else {
            don_dir <- if (exon_dom$exon_start < exon_non_dom$exon_start) "LOSS" else "GAIN"
          }
          second_event <- list(event_type = don_type, bp_diff = donor_diff, direction = don_dir)
        }
      }

    } else {
      # Purely internal exon pair: both boundaries differ.
      # Decompose into two independent splice-site events — one for each boundary.
      # Primary event: donor boundary (A5SS or Partial_IR_5)
      bp_diff <- donor_diff
      if (donor_diff < SPLICE_SITE_THRESHOLD) {
        event_type <- "A5SS"
      } else {
        event_type <- "Partial_IR_5"
      }
      # Direction for donor boundary
      if (strand == "+") {
        direction <- if (exon_dom$exon_end > exon_non_dom$exon_end) "LOSS" else "GAIN"
      } else {
        direction <- if (exon_dom$exon_start < exon_non_dom$exon_start) "LOSS" else "GAIN"
      }

      # Second event: acceptor boundary (A3SS or Partial_IR_3)
      acc_type <- if (acceptor_diff < SPLICE_SITE_THRESHOLD) "A3SS" else "Partial_IR_3"
      if (strand == "+") {
        acc_dir <- if (exon_dom$exon_start < exon_non_dom$exon_start) "LOSS" else "GAIN"
      } else {
        acc_dir <- if (exon_dom$exon_end > exon_non_dom$exon_end) "LOSS" else "GAIN"
      }
      second_event <- list(event_type = acc_type, bp_diff = acceptor_diff, direction = acc_dir)
    }
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
      # Check 5' end (skip if EITHER exon is first - 5' end is TSS, terminal)
      # ===========================================================================
      if (!is_first_exon_comp && !is_first_exon) {
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
      # Check 3' end (skip if EITHER exon is last - 3' end is TES, terminal)
      # ===========================================================================
      if (!is_last_exon_comp && !is_last_exon) {
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

  return(list(event_type = event_type, bp_diff = bp_diff, direction = direction,
              second_event = second_event))
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
        # Sub-case B1: Extends into flanking exon
        # Check distance first - small differences are splice site variations, not IR
        distance <- abs(comp_boundary - dom_boundary)
        bp_diff <- distance

        if (distance < SPLICE_SITE_THRESHOLD) {
          # <100bp: Treat as splice site variation even if it extends into flanking
          # This is an A5SS/A3SS with GAIN, not IR
          # Positional 5' = acceptor (A3SS), Positional 3' = donor (A5SS)
          event_type <- if (is_5prime) "A3SS" else "A5SS"
          direction <- "GAIN"
        } else {
          # ≥100bp: True intron retention
          event_type <- "IR"
          direction <- "GAIN"
        }
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

# ==============================================================================
# Junction Helper Functions
# ==============================================================================

#' Compute splice junctions from ordered exons
#'
#' @param exons_ordered Exons in biological order (TSS → TES)
#' @return Character vector of "left:right" junction strings (genomic, left < right)
compute_junctions <- function(exons_ordered) {
  if (nrow(exons_ordered) < 2) return(character(0))
  vapply(seq_len(nrow(exons_ordered) - 1), function(i) {
    e1 <- exons_ordered[i, ]
    e2 <- exons_ordered[i + 1, ]
    # For any two consecutive non-overlapping exons:
    #   min(end1, end2) = end of the genomically-left exon
    #   max(start1, start2) = start of the genomically-right exon
    # This correctly gives left:right (left < right) for both strands.
    sprintf("%d:%d", min(e1$exon_end, e2$exon_end), max(e1$exon_start, e2$exon_start))
  }, character(1))
}

#' Filter junctions where EITHER endpoint falls within [range_start, range_end]
#'
#' Captures junctions that "touch" the event range — e.g. junctions flanking a skipped exon.
#'
#' @param jxn_vec Character vector of "left:right" junction strings
#' @param range_start Genomic range start
#' @param range_end Genomic range end
#' @return Filtered junction vector
junctions_touching_range <- function(jxn_vec, range_start, range_end) {
  if (length(jxn_vec) == 0) return(character(0))
  jxn_vec[vapply(jxn_vec, function(jxn) {
    cc <- as.integer(strsplit(jxn, ":")[[1]])
    (cc[1] >= range_start && cc[1] <= range_end) ||
    (cc[2] >= range_start && cc[2] <= range_end)
  }, logical(1))]
}

#' Filter junctions where BOTH endpoints fall within [range_start, range_end]
#'
#' Captures junctions fully contained in the event range — e.g. dom junctions within an IR region.
#'
#' @param jxn_vec Character vector of "left:right" junction strings
#' @param range_start Genomic range start
#' @param range_end Genomic range end
#' @return Filtered junction vector
junctions_within_range <- function(jxn_vec, range_start, range_end) {
  if (length(jxn_vec) == 0) return(character(0))
  jxn_vec[vapply(jxn_vec, function(jxn) {
    cc <- as.integer(strsplit(jxn, ":")[[1]])
    cc[1] >= range_start && cc[2] <= range_end
  }, logical(1))]
}

#' Filter junctions where the junction span CONTAINS [range_start, range_end]
#'
#' Captures junctions that "skip over" the event range — e.g. a junction skipping a cassette exon.
#'
#' @param jxn_vec Character vector of "left:right" junction strings
#' @param range_start Genomic range start
#' @param range_end Genomic range end
#' @return Filtered junction vector
junctions_spanning_range <- function(jxn_vec, range_start, range_end) {
  if (length(jxn_vec) == 0) return(character(0))
  jxn_vec[vapply(jxn_vec, function(jxn) {
    cc <- as.integer(strsplit(jxn, ":")[[1]])
    cc[1] <= range_start && cc[2] >= range_end
  }, logical(1))]
}

#' Format junction vector as comma-separated string
#' @param jxn_vec Character vector of junction strings
#' @return Comma-separated string, or "" if empty
format_junctions <- function(jxn_vec) {
  if (length(jxn_vec) == 0) return("")
  paste(jxn_vec, collapse = ",")
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

# ==============================================================================
# Isoform Overlap Check
# ==============================================================================

#' Check if two isoforms have any exon-level overlap
#'
#' Used as a guard between Step 1 (terminal boundary detection) and Step 2
#' (intron retention detection). If no overlap exists, further event detection
#' is not meaningful and evaluation should end for the pair.
#'
#' @param comp_ordered Comparator exons ordered biologically (TSS → TES)
#' @param dom_ordered Dominant exons ordered biologically (TSS → TES)
#' @return List with:
#'   - has_overlap: logical, TRUE if any exon pair overlaps
#'   - dom_exon_coords: character, dominant exon coords as "start-end, ..." string
check_isoform_overlap <- function(comp_ordered, dom_ordered) {
  has_overlap <- FALSE
  for (i in seq_len(nrow(comp_ordered))) {
    for (j in seq_len(nrow(dom_ordered))) {
      if (comp_ordered$exon_start[i] <= dom_ordered$exon_end[j] &&
          comp_ordered$exon_end[i] >= dom_ordered$exon_start[j]) {
        has_overlap <- TRUE
        break
      }
    }
    if (has_overlap) break
  }

  dom_exon_coords <- paste(
    sprintf("%d-%d", dom_ordered$exon_start, dom_ordered$exon_end),
    collapse = ", "
  )

  list(has_overlap = has_overlap, dom_exon_coords = dom_exon_coords)
}

# ==============================================================================
# Beyond-Boundary Event Detection
# ==============================================================================

# ==============================================================================
# Event Detection for Transcript Pair
# ==============================================================================

#' Detect all events between two transcripts
#'
#' Runs the full event detection pipeline for a single dominant/comparator pair.
#' Detection proceeds in this order:
#'   1. Boundary determination (TSS/TES positions, genomic spans)
#'   2. Within-boundary events (only when isoforms overlap):
#'      2a. IR GAIN / IR LOSS  (intron retention)
#'      2b. A5SS / A3SS / Partial_IR  (splice-site changes on overlapping exons)
#'      2c. SE LOSS / SE GAIN  (skipped exons; strict flanking required)
#'   3. Terminal events (always runs, including non-overlapping isoforms):
#'      Alt_TSS / Alt_TES
#'
#' @param dominant_exons Exons from dominant transcript
#' @param comparator_exons Exons from comparator transcript
#' @param gene_id Gene identifier
#' @param dominant_tid Dominant transcript ID
#' @param comparator_tid Comparator transcript ID
#' @param strand Gene strand
#' @return Tibble with detected events (0 rows if none)
detect_events_for_pair <- function(dominant_exons, comparator_exons,
                                   gene_id, dominant_tid, comparator_tid, strand) {

  # Order exons biologically (TSS → TES)
  dom_ordered  <- order_exons_biological(dominant_exons, strand)
  comp_ordered <- order_exons_biological(comparator_exons, strand)

  # Pre-compute splice junctions for dom and comp (used for junction fields on all events)
  dom_junctions_all  <- compute_junctions(dom_ordered)
  comp_junctions_all <- compute_junctions(comp_ordered)

  # ===========================================================================
  # STEP 1: Boundary Determination
  # ===========================================================================

  first_dom  <- dom_ordered[1, ]
  first_comp <- comp_ordered[1, ]
  last_dom   <- dom_ordered[nrow(dom_ordered), ]
  last_comp  <- comp_ordered[nrow(comp_ordered), ]

  comp_genomic_start <- min(comp_ordered$exon_start)
  comp_genomic_end   <- max(comp_ordered$exon_end)
  dom_genomic_start  <- min(dom_ordered$exon_start)
  dom_genomic_end    <- max(dom_ordered$exon_end)

  all_events <- list()

  # Helper: check if a single exon overlaps any exon in an exon set (≥1bp)
  exon_overlaps_any <- function(exon, exon_set) {
    any(sapply(seq_len(nrow(exon_set)), function(k) {
      e <- exon_set[k, ]
      exon$exon_start <= e$exon_end && exon$exon_end >= e$exon_start
    }))
  }

  # ===========================================================================
  # STEP 2: Within-Boundary Events (skipped if isoforms do not overlap)
  # ===========================================================================

  overlap_check <- check_isoform_overlap(comp_ordered, dom_ordered)

  if (overlap_check$has_overlap) {

    # -------------------------------------------------------------------------
    # STEP 2a: IR Detection
    # -------------------------------------------------------------------------
    ir_exon_pairs <- list()

    # IR GAIN: comparator exon spans multiple dominant exons
    for (i in seq_len(nrow(comp_ordered))) {
      if (detect_ir_simple(comp_ordered[i, ], dom_ordered)) {
        comp_exon <- comp_ordered[i, ]

        overlapping_regions     <- list()
        overlapping_dom_indices <- integer()
        first_dom_overlap       <- NULL
        last_dom_overlap        <- NULL

        for (j in seq_len(nrow(dom_ordered))) {
          dom_exon <- dom_ordered[j, ]
          overlaps <- (dom_exon$exon_start <= comp_exon$exon_end) &&
                      (dom_exon$exon_end   >= comp_exon$exon_start)
          if (overlaps) {
            # Use FULL dom exon extents (not intersection) so reconstruction gets
            # the correct dom exon boundaries even when comp IR exon starts/ends inside dom
            overlapping_regions[[length(overlapping_regions) + 1]] <-
              sprintf("%d-%d", dom_exon$exon_start, dom_exon$exon_end)
            overlapping_dom_indices <- c(overlapping_dom_indices, j)
            if (is.null(first_dom_overlap)) first_dom_overlap <- dom_exon
            last_dom_overlap <- dom_exon
          }
        }

        # Classify IR_diff type: compare comp's retained exon boundaries against
        # dom's outermost overlapping exon boundaries.
        # diff_5: comp started after dom's first exon (comp missed dom's 5' extension)
        # diff_3: comp ended before dom's last exon (comp missed dom's 3' extension)
        ir_event_type <- "IR"
        if (!is.null(first_dom_overlap)) {
          if (strand == "+") {
            diff_5 <- comp_exon$exon_start > first_dom_overlap$exon_start
            diff_3 <- comp_exon$exon_end   < last_dom_overlap$exon_end
          } else {
            # Minus: 5' = high coords, 3' = low coords
            diff_5 <- comp_exon$exon_end   < first_dom_overlap$exon_end
            diff_3 <- comp_exon$exon_start > last_dom_overlap$exon_start
          }
          ir_event_type <- if (diff_5 && diff_3) "IR_diff_5_3" else
                           if (diff_5) "IR_diff_5" else
                           if (diff_3) "IR_diff_3" else "IR"
        }

        # Junction fields: dom has split exons (junctions present); comp retains intron (no junctions)
        ir_range_start <- min(sapply(overlapping_dom_indices, function(k) dom_ordered$exon_start[k]))
        ir_range_end   <- max(sapply(overlapping_dom_indices, function(k) dom_ordered$exon_end[k]))
        dom_jxns_ir    <- junctions_within_range(dom_junctions_all, ir_range_start, ir_range_end)

        all_events[[length(all_events) + 1]] <- tibble(
          gene_id                     = gene_id,
          dominant_transcript_id      = dominant_tid,
          comparator_transcript_id    = comparator_tid,
          event_type                  = ir_event_type,
          direction                   = "GAIN",
          chr                         = comp_exon$chr,
          five_prime                  = if (strand == "+") comp_exon$exon_start else comp_exon$exon_end,
          three_prime                 = if (strand == "+") comp_exon$exon_end   else comp_exon$exon_start,
          strand                      = strand,
          bp_diff                     = NA_integer_,
          missing_terminal_exons      = "",
          orphan_terminal_exons       = "",
          ir_split_exons              = paste(overlapping_regions, collapse = ", "),
          dom_junctions               = format_junctions(dom_jxns_ir),
          comp_junctions              = ""
        )

        for (dom_idx in overlapping_dom_indices)
          ir_exon_pairs[[length(ir_exon_pairs) + 1]] <- list(comp = i, dom = dom_idx)
      }
    }

    # IR LOSS: dominant exon spans multiple comparator exons
    for (i in seq_len(nrow(dom_ordered))) {
      if (detect_ir_simple(dom_ordered[i, ], comp_ordered)) {
        dom_exon <- dom_ordered[i, ]

        overlapping_regions      <- list()
        overlapping_comp_indices <- integer()
        first_comp_overlap       <- NULL
        last_comp_overlap        <- NULL

        for (j in seq_len(nrow(comp_ordered))) {
          comp_exon <- comp_ordered[j, ]
          overlaps <- (comp_exon$exon_start <= dom_exon$exon_end) &&
                      (comp_exon$exon_end   >= dom_exon$exon_start)
          if (overlaps) {
            # ir_split_exons for LOSS keeps intersection (used by reconstruction to identify
            # exonic union exons within the dom big exon; the intronic fill step handles the rest)
            overlapping_regions[[length(overlapping_regions) + 1]] <-
              sprintf("%d-%d",
                      max(comp_exon$exon_start, dom_exon$exon_start),
                      min(comp_exon$exon_end,   dom_exon$exon_end))
            overlapping_comp_indices <- c(overlapping_comp_indices, j)
            if (is.null(first_comp_overlap)) first_comp_overlap <- comp_exon
            last_comp_overlap <- comp_exon
          }
        }

        # Classify IR_diff type: compare dom's big exon boundaries against
        # comp's outermost overlapping exon boundaries.
        # diff_5: dom starts before comp's first exon (dom has extra 5' extension)
        # diff_3: dom ends after comp's last exon (dom has extra 3' extension)
        ir_event_type <- "IR"
        if (!is.null(first_comp_overlap)) {
          if (strand == "+") {
            diff_5 <- dom_exon$exon_start < first_comp_overlap$exon_start
            diff_3 <- dom_exon$exon_end   > last_comp_overlap$exon_end
          } else {
            # Minus: 5' = high coords, 3' = low coords
            diff_5 <- dom_exon$exon_end   > first_comp_overlap$exon_end
            diff_3 <- dom_exon$exon_start < last_comp_overlap$exon_start
          }
          ir_event_type <- if (diff_5 && diff_3) "IR_diff_5_3" else
                           if (diff_5) "IR_diff_5" else
                           if (diff_3) "IR_diff_3" else "IR"
        }

        # Junction fields: dom has retention (no junctions); comp has split exons (junctions present)
        comp_jxns_ir <- junctions_within_range(comp_junctions_all,
                                               dom_exon$exon_start, dom_exon$exon_end)

        all_events[[length(all_events) + 1]] <- tibble(
          gene_id                     = gene_id,
          dominant_transcript_id      = dominant_tid,
          comparator_transcript_id    = comparator_tid,
          event_type                  = ir_event_type,
          direction                   = "LOSS",
          chr                         = dom_exon$chr,
          five_prime                  = if (strand == "+") dom_exon$exon_start else dom_exon$exon_end,
          three_prime                 = if (strand == "+") dom_exon$exon_end   else dom_exon$exon_start,
          strand                      = strand,
          bp_diff                     = NA_integer_,
          missing_terminal_exons      = "",
          orphan_terminal_exons       = "",
          ir_split_exons              = paste(overlapping_regions, collapse = ", "),
          dom_junctions               = "",
          comp_junctions              = format_junctions(comp_jxns_ir)
        )

        for (comp_idx in overlapping_comp_indices)
          ir_exon_pairs[[length(ir_exon_pairs) + 1]] <- list(comp = comp_idx, dom = i)
      }
    }

    # -------------------------------------------------------------------------
    # STEP 2b: Boundary Events (A5SS, A3SS, Partial_IR) for non-IR exon pairs
    # -------------------------------------------------------------------------
    for (i in seq_len(nrow(comp_ordered))) {
      comp_exon     <- comp_ordered[i, ]
      is_first_comp <- (i == 1)
      is_last_comp  <- (i == nrow(comp_ordered))

      for (j in seq_len(nrow(dom_ordered))) {
        dom_exon     <- dom_ordered[j, ]
        is_first_dom <- (j == 1)
        is_last_dom  <- (j == nrow(dom_ordered))

        overlaps <- (dom_exon$exon_start <= comp_exon$exon_end) &&
                    (dom_exon$exon_end   >= comp_exon$exon_start)
        if (!overlaps) next

        has_ir <- any(sapply(ir_exon_pairs, function(pair) pair$comp == i && pair$dom == j))
        if (has_ir) next

        event_result <- detect_shared_boundary_event(
          dom_exon, comp_exon, strand,
          is_first_exon      = is_first_dom && is_first_comp,
          is_last_exon       = is_last_dom  && is_last_comp,
          is_first_exon_dom  = is_first_dom,
          is_last_exon_dom   = is_last_dom,
          is_first_exon_comp = is_first_comp,
          is_last_exon_comp  = is_last_comp,
          flanking_exons_dom     = NULL,
          flanking_exons_non_dom = NULL
        )

        if (event_result$event_type != "none") {
          event_type <- event_result$event_type
          direction  <- if (!is.null(event_result$direction)) event_result$direction else "-"

          if (event_type %in% c("A3SS", "Partial_IR_3")) {
            if (strand == "+") {
              if (direction == "LOSS") { five_prime <- dom_exon$exon_start;  three_prime <- comp_exon$exon_start - 1 }
              else                    { five_prime <- comp_exon$exon_start;  three_prime <- dom_exon$exon_start - 1 }
            } else {
              if (direction == "LOSS") { five_prime <- dom_exon$exon_end;    three_prime <- comp_exon$exon_end + 1 }
              else                    { five_prime <- comp_exon$exon_end;    three_prime <- dom_exon$exon_end + 1 }
            }
          } else if (event_type %in% c("A5SS", "Partial_IR_5")) {
            if (strand == "+") {
              if (direction == "LOSS") { five_prime <- comp_exon$exon_end + 1; three_prime <- dom_exon$exon_end }
              else                    { five_prime <- dom_exon$exon_end + 1;   three_prime <- comp_exon$exon_end }
            } else {
              if (direction == "LOSS") { five_prime <- comp_exon$exon_start;   three_prime <- dom_exon$exon_start - 1 }
              else                    { five_prime <- dom_exon$exon_start - 1; three_prime <- comp_exon$exon_start }
            }
          } else {
            if (strand == "+") { five_prime <- min(dom_exon$exon_start, comp_exon$exon_start); three_prime <- max(dom_exon$exon_end, comp_exon$exon_end) }
            else               { five_prime <- max(dom_exon$exon_end,   comp_exon$exon_end);   three_prime <- min(dom_exon$exon_start, comp_exon$exon_start) }
          }

          # Junction fields: both isoforms' junctions at the affected splice site
          evt_start     <- min(five_prime, three_prime)
          evt_end       <- max(five_prime, three_prime)
          dom_jxns_evt  <- junctions_touching_range(dom_junctions_all,  evt_start, evt_end)
          comp_jxns_evt <- junctions_touching_range(comp_junctions_all, evt_start, evt_end)

          all_events[[length(all_events) + 1]] <- tibble(
            gene_id                     = gene_id,
            dominant_transcript_id      = dominant_tid,
            comparator_transcript_id    = comparator_tid,
            event_type                  = event_type,
            direction                   = direction,
            chr                         = dom_exon$chr,
            five_prime                  = five_prime,
            three_prime                 = three_prime,
            strand                      = strand,
            bp_diff                     = event_result$bp_diff,
            missing_terminal_exons      = "",
            orphan_terminal_exons       = "",
            ir_split_exons              = "",
            dom_junctions               = format_junctions(dom_jxns_evt),
            comp_junctions              = format_junctions(comp_jxns_evt)
          )

          # Emit second event for asymmetric dual-boundary pairs (one exon is terminal,
          # other is not), where the splice-site-facing boundary is a real splice site.
          if (!is.null(event_result$second_event)) {
            event_type2 <- event_result$second_event$event_type
            direction2  <- event_result$second_event$direction

            if (event_type2 %in% c("A3SS", "Partial_IR_3")) {
              if (strand == "+") {
                if (direction2 == "LOSS") { five_prime2 <- dom_exon$exon_start;  three_prime2 <- comp_exon$exon_start - 1 }
                else                     { five_prime2 <- comp_exon$exon_start;  three_prime2 <- dom_exon$exon_start - 1 }
              } else {
                if (direction2 == "LOSS") { five_prime2 <- dom_exon$exon_end;    three_prime2 <- comp_exon$exon_end + 1 }
                else                     { five_prime2 <- comp_exon$exon_end;    three_prime2 <- dom_exon$exon_end + 1 }
              }
            } else if (event_type2 %in% c("A5SS", "Partial_IR_5")) {
              if (strand == "+") {
                if (direction2 == "LOSS") { five_prime2 <- comp_exon$exon_end + 1; three_prime2 <- dom_exon$exon_end }
                else                     { five_prime2 <- dom_exon$exon_end + 1;   three_prime2 <- comp_exon$exon_end }
              } else {
                if (direction2 == "LOSS") { five_prime2 <- comp_exon$exon_start;   three_prime2 <- dom_exon$exon_start - 1 }
                else                     { five_prime2 <- dom_exon$exon_start - 1; three_prime2 <- comp_exon$exon_start }
              }
            } else {
              five_prime2  <- min(dom_exon$exon_start, comp_exon$exon_start)
              three_prime2 <- max(dom_exon$exon_end,   comp_exon$exon_end)
            }

            evt_start2     <- min(five_prime2, three_prime2)
            evt_end2       <- max(five_prime2, three_prime2)
            dom_jxns_evt2  <- junctions_touching_range(dom_junctions_all,  evt_start2, evt_end2)
            comp_jxns_evt2 <- junctions_touching_range(comp_junctions_all, evt_start2, evt_end2)

            all_events[[length(all_events) + 1]] <- tibble(
              gene_id                     = gene_id,
              dominant_transcript_id      = dominant_tid,
              comparator_transcript_id    = comparator_tid,
              event_type                  = event_type2,
              direction                   = direction2,
              chr                         = dom_exon$chr,
              five_prime                  = five_prime2,
              three_prime                 = three_prime2,
              strand                      = strand,
              bp_diff                     = event_result$second_event$bp_diff,
              missing_terminal_exons      = "",
              orphan_terminal_exons       = "",
              ir_split_exons              = "",
              dom_junctions               = format_junctions(dom_jxns_evt2),
              comp_junctions              = format_junctions(comp_jxns_evt2)
            )
          }
        }
      }
    }

    # -------------------------------------------------------------------------
    # STEP 2c: SE / Missing_Internal Detection
    #
    # For each exon absent in the other isoform (within span, no overlap):
    #   - If strict flanking met (both neighbors overlap other isoform) → SE
    #   - Otherwise → Missing_Internal
    #
    # Both event types drive reconstruction identically:
    #   LOSS: add dom exon to comp    GAIN: remove comp exon from reconstructed
    #
    # Terminal exons (first/last) are excluded here; they are handled by
    # orphan_terminal_exons on the Alt_TSS/Alt_TES events.
    # -------------------------------------------------------------------------

    # LOSS: dom exon absent in comp (within comp span, no comp overlap)
    for (i in seq_len(nrow(dom_ordered))) {
      dom_exon <- dom_ordered[i, ]

      # Must be fully within comp span
      if (dom_exon$exon_start < comp_genomic_start || dom_exon$exon_end > comp_genomic_end) next

      # Must have no overlap with any comp exon
      if (exon_overlaps_any(dom_exon, comp_ordered)) next

      # Strict flanking check determines event label.
      # Terminal dom exons are included here: they have no neighbor on one side so
      # flanking_ok is FALSE → Missing_Internal LOSS. This handles the NoOverlap GAIN
      # case where dom's terminal exon needs to be added to comp.
      has_prev    <- i > 1
      has_next    <- i < nrow(dom_ordered)
      flanking_ok <- has_prev && has_next &&
                     exon_overlaps_any(dom_ordered[i - 1, ], comp_ordered) &&
                     exon_overlaps_any(dom_ordered[i + 1, ], comp_ordered)
      etype <- if (flanking_ok) "SE" else "Missing_Internal"

      # Junction fields:
      #   dom: junctions flanking the SE (both endpoints touch the exon range)
      #   comp: the skipping junction (spans the exon range)
      dom_jxns_se  <- junctions_touching_range(dom_junctions_all,
                                               dom_exon$exon_start, dom_exon$exon_end)
      comp_jxns_se <- junctions_spanning_range(comp_junctions_all,
                                               dom_exon$exon_start, dom_exon$exon_end)

      all_events[[length(all_events) + 1]] <- tibble(
        gene_id                     = gene_id,
        dominant_transcript_id      = dominant_tid,
        comparator_transcript_id    = comparator_tid,
        event_type                  = etype,
        direction                   = "LOSS",
        chr                         = dom_exon$chr,
        five_prime                  = if (strand == "+") dom_exon$exon_start else dom_exon$exon_end,
        three_prime                 = if (strand == "+") dom_exon$exon_end   else dom_exon$exon_start,
        strand                      = strand,
        bp_diff                     = dom_exon$exon_end - dom_exon$exon_start + 1,
        missing_terminal_exons      = "",
        orphan_terminal_exons       = "",
        ir_split_exons              = "",
        dom_junctions               = format_junctions(dom_jxns_se),
        comp_junctions              = format_junctions(comp_jxns_se)
      )
    }

    # GAIN: comp exon absent in dom (within dom span, no dom overlap)
    for (i in seq_len(nrow(comp_ordered))) {
      comp_exon <- comp_ordered[i, ]

      # Must be fully within dom span
      if (comp_exon$exon_start < dom_genomic_start || comp_exon$exon_end > dom_genomic_end) next

      # Must have no overlap with any dom exon
      if (exon_overlaps_any(comp_exon, dom_ordered)) next

      # Exclude terminal comp exons — handled by orphan_terminal_exons on Alt_TSS/Alt_TES
      if (i == 1 || i == nrow(comp_ordered)) next

      # Strict flanking check determines event label
      flanking_ok <- exon_overlaps_any(comp_ordered[i - 1, ], dom_ordered) &&
                     exon_overlaps_any(comp_ordered[i + 1, ], dom_ordered)
      etype <- if (flanking_ok) "SE" else "Missing_Internal"

      # Junction fields:
      #   dom: the skipping junction (spans the exon range)
      #   comp: junctions flanking the SE (both endpoints touch the exon range)
      dom_jxns_se  <- junctions_spanning_range(dom_junctions_all,
                                               comp_exon$exon_start, comp_exon$exon_end)
      comp_jxns_se <- junctions_touching_range(comp_junctions_all,
                                               comp_exon$exon_start, comp_exon$exon_end)

      all_events[[length(all_events) + 1]] <- tibble(
        gene_id                     = gene_id,
        dominant_transcript_id      = dominant_tid,
        comparator_transcript_id    = comparator_tid,
        event_type                  = etype,
        direction                   = "GAIN",
        chr                         = comp_exon$chr,
        five_prime                  = if (strand == "+") comp_exon$exon_start else comp_exon$exon_end,
        three_prime                 = if (strand == "+") comp_exon$exon_end   else comp_exon$exon_start,
        strand                      = strand,
        bp_diff                     = comp_exon$exon_end - comp_exon$exon_start + 1,
        missing_terminal_exons      = "",
        orphan_terminal_exons       = "",
        ir_split_exons              = "",
        dom_junctions               = format_junctions(dom_jxns_se),
        comp_junctions              = format_junctions(comp_jxns_se)
      )
    }

  } else {
    # Non-overlapping isoforms: the TSS/TES terminal walks (step 3) may miss
    # dominant exons that fall between the two walk cutoffs ("gap zone").
    # Detect them here as Missing_Internal LOSS events.
    for (i in seq_len(nrow(dom_ordered))) {
      dom_exon <- dom_ordered[i, ]

      # Must be fully within comp span
      if (dom_exon$exon_start < comp_genomic_start || dom_exon$exon_end > comp_genomic_end) next

      # Must have no overlap with any comp exon
      if (exon_overlaps_any(dom_exon, comp_ordered)) next

      dom_jxns_gap  <- junctions_touching_range(dom_junctions_all,
                                                dom_exon$exon_start, dom_exon$exon_end)
      comp_jxns_gap <- junctions_spanning_range(comp_junctions_all,
                                                dom_exon$exon_start, dom_exon$exon_end)

      all_events[[length(all_events) + 1]] <- tibble(
        gene_id                     = gene_id,
        dominant_transcript_id      = dominant_tid,
        comparator_transcript_id    = comparator_tid,
        event_type                  = "Missing_Internal",
        direction                   = "LOSS",
        chr                         = dom_exon$chr,
        five_prime                  = if (strand == "+") dom_exon$exon_start else dom_exon$exon_end,
        three_prime                 = if (strand == "+") dom_exon$exon_end   else dom_exon$exon_start,
        strand                      = strand,
        bp_diff                     = dom_exon$exon_end - dom_exon$exon_start + 1,
        missing_terminal_exons      = "",
        orphan_terminal_exons       = "",
        ir_split_exons              = "",
        dom_junctions               = format_junctions(dom_jxns_gap),
        comp_junctions              = format_junctions(comp_jxns_gap)
      )
    }
  }  # end if (overlap_check$has_overlap) / else

  # ===========================================================================
  # STEP 3: Terminal Events (always runs, including non-overlapping isoforms)
  # ===========================================================================
  #
  # Alt_TSS/Alt_TES are detected last.  For non-overlapping isoforms, these
  # two events together capture the complete structural difference.
  # missing_terminal_exons is computed by walking the longer isoform's exons
  # from its terminal end to the shorter isoform's TSS/TES boundary.

  # ---------------------------------------------------------------------------
  # Alt_TSS
  # ---------------------------------------------------------------------------
  if (detect_tss_change(first_dom, first_comp, strand)) {
    if (strand == "+") { dom_tss <- first_dom$exon_start; comp_tss <- first_comp$exon_start }
    else               { dom_tss <- first_dom$exon_end;   comp_tss <- first_comp$exon_end   }

    # LOSS: dominant starts earlier (comparator lost 5' sequence)
    # GAIN: comparator starts earlier (comparator gained 5' sequence)
    direction <- if (
      (strand == "+" && dom_tss < comp_tss) ||
      (strand == "-" && dom_tss > comp_tss)
    ) "LOSS" else "GAIN"

    missing_exons <- if (direction == "LOSS")
      compute_missing_terminal_exons_tss(dom_ordered, comp_ordered, strand)
    else
      compute_missing_terminal_exons_tss(comp_ordered, dom_ordered, strand)

    n_regions <- if (missing_exons == "") 0L else length(strsplit(missing_exons, ",")[[1]])

    # Orphaned terminal exons: comp exons to REMOVE for both directions.
    # LOSS: comp TSS exons that don't overlap dom (NoOverlap case).
    # GAIN: comp extra TSS exons beyond dom's TSS boundary.
    # Dom's missing internal exons are now emitted as standalone Missing_Internal events.
    orphan_exons <- compute_orphan_terminal_exons_tss(comp_ordered, dom_ordered)

    # Junction fields: first junction of each isoform (captures the TSS-proximal splice choice)
    dom_jxns_tss  <- if (length(dom_junctions_all)  > 0) dom_junctions_all[1]  else character(0)
    comp_jxns_tss <- if (length(comp_junctions_all) > 0) comp_junctions_all[1] else character(0)

    all_events[[length(all_events) + 1]] <- tibble(
      gene_id                     = gene_id,
      dominant_transcript_id      = dominant_tid,
      comparator_transcript_id    = comparator_tid,
      event_type                  = "Alt_TSS",
      direction                   = direction,
      chr                         = first_dom$chr,
      five_prime                  = dom_tss,
      three_prime                 = comp_tss,
      strand                      = strand,
      bp_diff                     = abs(dom_tss - comp_tss),
      missing_terminal_exons      = missing_exons,
      orphan_terminal_exons       = orphan_exons,
      ir_split_exons              = "",
      dom_junctions               = format_junctions(dom_jxns_tss),
      comp_junctions              = format_junctions(comp_jxns_tss)
    )
  }

  # ---------------------------------------------------------------------------
  # Alt_TES
  # ---------------------------------------------------------------------------
  if (detect_tes_change(last_dom, last_comp, strand)) {
    if (strand == "+") { dom_tes <- last_dom$exon_end;   comp_tes <- last_comp$exon_end   }
    else               { dom_tes <- last_dom$exon_start; comp_tes <- last_comp$exon_start }

    # LOSS: dominant ends later (comparator lost 3' sequence)
    # GAIN: comparator ends later (comparator gained 3' sequence)
    direction <- if (
      (strand == "+" && dom_tes > comp_tes) ||
      (strand == "-" && dom_tes < comp_tes)
    ) "LOSS" else "GAIN"

    missing_exons <- if (direction == "LOSS")
      compute_missing_terminal_exons_tes(dom_ordered, comp_ordered, strand)
    else
      compute_missing_terminal_exons_tes(comp_ordered, dom_ordered, strand)

    n_regions <- if (missing_exons == "") 0L else length(strsplit(missing_exons, ",")[[1]])

    # Orphaned terminal exons: comp exons to REMOVE for both directions.
    # LOSS: comp TES exons that don't overlap dom (NoOverlap case).
    # GAIN: comp extra TES exons beyond dom's TES boundary.
    # Dom's missing internal exons are now emitted as standalone Missing_Internal events.
    orphan_exons <- compute_orphan_terminal_exons_tes(comp_ordered, dom_ordered)

    # Junction fields: last junction of each isoform (captures the TES-proximal splice choice)
    dom_jxns_tes  <- if (length(dom_junctions_all)  > 0) dom_junctions_all[length(dom_junctions_all)]   else character(0)
    comp_jxns_tes <- if (length(comp_junctions_all) > 0) comp_junctions_all[length(comp_junctions_all)] else character(0)

    all_events[[length(all_events) + 1]] <- tibble(
      gene_id                     = gene_id,
      dominant_transcript_id      = dominant_tid,
      comparator_transcript_id    = comparator_tid,
      event_type                  = "Alt_TES",
      direction                   = direction,
      chr                         = last_dom$chr,
      five_prime                  = dom_tes,
      three_prime                 = comp_tes,
      strand                      = strand,
      bp_diff                     = abs(dom_tes - comp_tes),
      missing_terminal_exons      = missing_exons,
      orphan_terminal_exons       = orphan_exons,
      ir_split_exons              = "",
      dom_junctions               = format_junctions(dom_jxns_tes),
      comp_junctions              = format_junctions(comp_jxns_tes)
    )
  }

  if (length(all_events) > 0) return(bind_rows(all_events))
  return(tibble())
}
