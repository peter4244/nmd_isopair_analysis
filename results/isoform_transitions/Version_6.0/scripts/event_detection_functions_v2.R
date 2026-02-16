#!/usr/bin/env Rscript
#
# Event Detection Functions Library - Version 2
# Intron-centric approach with 50% overlap for "same" exons
#
# Key changes from v1:
# - Unified 50% overlap definition for "same" exons across all events
# - Intron-centric detection (focus on introns, not exact boundaries)
# - IR = ANY overlap with multiple exons (entire intron subsumed)
# - Partial_IR/A5SS/A3SS based on intron coordinate differences
# - One-to-one exon matching with greedy algorithm
#

# ==============================================================================
# Detection Thresholds (Constants)
# ==============================================================================

# Exon matching
OVERLAP_THRESHOLD <- 0.5  # 50% overlap required for exons to be "same"

# Intron boundary classification
SPLICE_SITE_THRESHOLD <- 100  # bp threshold: <100 = splice site, ≥100 = retention

# Terminal boundary tolerance
TSS_TOLERANCE <- 20  # bp tolerance for TSS change detection
TES_TOLERANCE <- 20  # bp tolerance for TES change detection

# Micro-exon definition
MICRO_EXON_THRESHOLD <- 30  # bp threshold for micro-exon flagging

# ==============================================================================
# GTF Validation Functions
# ==============================================================================

#' Validate GTF Exon Ordering
#'
#' Checks that exon numbering follows transcriptional order (5' to 3')
#' regardless of strand. For minus strand, exon 1 should be at the highest
#' (rightmost) coordinates, not the lowest.
#'
#' @param gtf_path Path to GTF file
#' @param full_validation If TRUE, validate entire file. If FALSE, validate first 10000 lines (default)
#' @param verbose If TRUE, print detailed validation info
#' @return List with validation results: valid (logical), issues (data.frame), summary (character)
#' @export
validate_gtf_exon_order <- function(gtf_path, full_validation = FALSE, verbose = TRUE) {

  if (!file.exists(gtf_path)) {
    stop(sprintf("GTF file not found: %s", gtf_path))
  }

  if (verbose) {
    cat(sprintf("\n╔══════════════════════════════════════════════════════════════╗\n"))
    cat(sprintf("║   GTF Validation: %s\n", basename(gtf_path)))
    cat(sprintf("╚══════════════════════════════════════════════════════════════╝\n\n"))
  }

  # Read GTF (full or partial)
  if (full_validation) {
    if (verbose) cat("Reading full GTF file...\n")
    gtf <- read_tsv(gtf_path,
                    col_names = c("seqnames", "source", "feature", "start", "end",
                                 "score", "strand", "frame", "attributes"),
                    comment = "#", show_col_types = FALSE)
  } else {
    if (verbose) cat("Reading first 10000 lines from GTF...\n")
    gtf <- read_tsv(gtf_path,
                    col_names = c("seqnames", "source", "feature", "start", "end",
                                 "score", "strand", "frame", "attributes"),
                    comment = "#", show_col_types = FALSE, n_max = 10000)
  }

  # Filter to exons only
  gtf <- gtf %>% filter(feature == "exon")

  if (nrow(gtf) == 0) {
    if (verbose) cat("✗ No exon features found in GTF\n")
    return(list(valid = FALSE, issues = NULL, summary = "No exons found"))
  }

  # Parse attributes
  gtf <- gtf %>%
    mutate(
      gene_id = str_match(attributes, 'gene_id "([^"]+)"')[,2],
      transcript_id = str_match(attributes, 'transcript_id "([^"]+)"')[,2],
      exon_number = as.integer(str_match(attributes, 'exon_number "([^"]+)"')[,2])
    ) %>%
    filter(!is.na(exon_number))

  if (verbose) {
    cat(sprintf("  Exons: %d\n", nrow(gtf)))
    cat(sprintf("  Transcripts: %d\n", length(unique(gtf$transcript_id))))
    cat(sprintf("  Genes: %d\n\n", length(unique(gtf$gene_id))))
  }

  # Validate exon ordering for each transcript
  validation <- gtf %>%
    group_by(transcript_id, strand) %>%
    arrange(transcript_id, exon_number) %>%
    summarise(
      n_exons = n(),
      # For plus strand: exon coordinates should be ascending
      # For minus strand: exon coordinates should be descending
      coords = list(start),
      expected_order = if(first(strand) == "+") {
        all(diff(start) > 0)  # Ascending for plus
      } else {
        all(diff(start) < 0)  # Descending for minus
      },
      .groups = "drop"
    ) %>%
    filter(n_exons > 1)  # Only check multi-exonic transcripts

  issues <- validation %>% filter(!expected_order)

  if (nrow(issues) > 0) {
    if (verbose) {
      cat(sprintf("✗ VALIDATION FAILED\n"))
      cat(sprintf("  %d/%d multi-exonic transcripts have incorrect exon ordering\n\n",
                  nrow(issues), nrow(validation)))
      cat("Examples of incorrect ordering:\n")
      cat("─────────────────────────────────────────────────────────────\n")
      print(issues %>% head(10) %>% select(transcript_id, strand, n_exons))
      cat("\n")
    }

    summary_msg <- sprintf("FAILED: %d/%d transcripts have incorrect exon order",
                          nrow(issues), nrow(validation))
    return(list(valid = FALSE, issues = issues, summary = summary_msg))
  } else {
    if (verbose) {
      cat(sprintf("✓ VALIDATION PASSED\n"))
      cat(sprintf("  All %d multi-exonic transcripts have correct exon ordering\n",
                  nrow(validation)))
      cat(sprintf("  Plus strand: exon coordinates ascending (5' to 3')\n"))
      cat(sprintf("  Minus strand: exon coordinates descending (5' to 3')\n\n"))
    }

    summary_msg <- sprintf("PASSED: %d transcripts validated", nrow(validation))
    return(list(valid = TRUE, issues = NULL, summary = summary_msg))
  }
}

#' Fix GTF Exon Ordering
#'
#' Corrects exon numbering to follow transcriptional order (5' to 3').
#' For minus strand, reverses exon numbering so exon 1 is at highest coordinates.
#'
#' @param gtf_path Path to GTF file to fix
#' @param output_path Path for corrected GTF (default: overwrites input)
#' @return Invisible TRUE if successful
#' @export
fix_gtf_exon_order <- function(gtf_path, output_path = gtf_path) {

  cat(sprintf("\nFixing exon ordering in: %s\n", gtf_path))

  # Read GTF
  gtf <- read_tsv(gtf_path,
                  col_names = c("seqnames", "source", "feature", "start", "end",
                               "score", "strand", "frame", "attributes"),
                  comment = "#", show_col_types = FALSE)

  # Parse attributes
  gtf <- gtf %>%
    mutate(
      gene_id = str_match(attributes, 'gene_id "([^"]+)"')[,2],
      transcript_id = str_match(attributes, 'transcript_id "([^"]+)"')[,2],
      exon_number = as.integer(str_match(attributes, 'exon_number "([^"]+)"')[,2])
    )

  # Fix exon numbering
  gtf_fixed <- gtf %>%
    filter(!is.na(exon_number)) %>%
    group_by(transcript_id) %>%
    mutate(
      # For minus strand, rank by descending coordinate
      # For plus strand, rank by ascending coordinate
      correct_exon_number = if(first(strand) == "-") {
        rank(desc(start), ties.method = "first")
      } else {
        rank(start, ties.method = "first")
      }
    ) %>%
    ungroup() %>%
    mutate(
      exon_number = correct_exon_number,
      attributes = sprintf('gene_id "%s"; transcript_id "%s"; exon_number "%d";',
                          gene_id, transcript_id, exon_number)
    ) %>%
    select(seqnames, source, feature, start, end, score, strand, frame, attributes)

  # Write fixed GTF
  write_tsv(gtf_fixed, output_path, col_names = FALSE, quote = "none")

  cat(sprintf("✓ Fixed GTF written to: %s\n\n", output_path))

  invisible(TRUE)
}

# ==============================================================================
# Helper Functions - Overlap Calculations
# ==============================================================================

#' Check if two exons have ANY overlap
#'
#' @param exon1 First exon (tibble row with exon_start, exon_end)
#' @param exon2 Second exon (tibble row with exon_start, exon_end)
#' @return Logical - TRUE if any overlap exists
has_any_overlap <- function(exon1, exon2) {
  overlap_start <- max(exon1$exon_start, exon2$exon_start)
  overlap_end <- min(exon1$exon_end, exon2$exon_end)
  overlap_length <- overlap_end - overlap_start + 1
  return(overlap_length > 0)
}

#' Calculate overlap proportion between two exons
#'
#' @param exon1 First exon (tibble row with exon_start, exon_end)
#' @param exon2 Second exon (tibble row with exon_start, exon_end)
#' @return Numeric - overlap as proportion of shorter exon (0 to 1)
calculate_overlap_proportion <- function(exon1, exon2) {
  # Calculate overlap region
  overlap_start <- max(exon1$exon_start, exon2$exon_start)
  overlap_end <- min(exon1$exon_end, exon2$exon_end)
  overlap_length <- max(0, overlap_end - overlap_start + 1)

  # Calculate exon lengths
  len1 <- exon1$exon_end - exon1$exon_start + 1
  len2 <- exon2$exon_end - exon2$exon_start + 1

  # Overlap percentage relative to shorter exon
  if (min(len1, len2) == 0) return(0)
  overlap_pct <- overlap_length / min(len1, len2)

  return(overlap_pct)
}

#' Count how many exons from other isoform this exon overlaps with (ANY overlap)
#'
#' @param exon Exon to check (tibble row with exon_start, exon_end)
#' @param other_exons Exons from other isoform (tibble with exon_start, exon_end)
#' @return Integer - count of overlapping exons
count_overlapping_exons <- function(exon, other_exons) {
  if (nrow(other_exons) == 0) return(0)

  count <- 0
  for (i in seq_len(nrow(other_exons))) {
    if (has_any_overlap(exon, other_exons[i, ])) {
      count <- count + 1
    }
  }
  return(count)
}

# ==============================================================================
# IR Detection (ANY overlap with multiple exons)
# ==============================================================================

#' Detect IR: exon has ANY overlap with multiple exons from other isoform
#'
#' @param exons_a Exons from isoform A (tibble with exon_start, exon_end)
#' @param exons_b Exons from isoform B (tibble with exon_start, exon_end)
#' @return List with ir_in_a (indices) and ir_in_b (indices)
detect_ir_v2 <- function(exons_a, exons_b) {
  ir_in_a <- integer()  # Exons in A that span multiple B exons
  ir_in_b <- integer()  # Exons in B that span multiple A exons

  # Check each exon in A
  for (i in seq_len(nrow(exons_a))) {
    n_overlaps <- count_overlapping_exons(exons_a[i, ], exons_b)
    if (n_overlaps >= 2) {
      ir_in_a <- c(ir_in_a, i)
    }
  }

  # Check each exon in B
  for (j in seq_len(nrow(exons_b))) {
    n_overlaps <- count_overlapping_exons(exons_b[j, ], exons_a)
    if (n_overlaps >= 2) {
      ir_in_b <- c(ir_in_b, j)
    }
  }

  return(list(
    ir_in_a = ir_in_a,
    ir_in_b = ir_in_b
  ))
}

# ==============================================================================
# One-to-One Exon Matching
# ==============================================================================

#' Check if terminal exons match based on splice site containment
#'
#' For first exons: match if donor is contained in both
#' For last exons: match if acceptor is contained in both
#'
#' @param exon_a Exon from isoform A
#' @param exon_b Exon from isoform B
#' @param strand Strand
#' @param is_first Logical - are these first exons?
#' @param is_last Logical - are these last exons?
#' @return Logical - TRUE if they match by terminal exon rule
terminal_exons_match <- function(exon_a, exon_b, strand, is_first = FALSE, is_last = FALSE) {
  if (!is_first && !is_last) return(FALSE)

  if (strand == "+") {
    if (is_first) {
      # First exon on plus: donor is exon_end (internal boundary)
      donor_a <- exon_a$exon_end
      donor_b <- exon_b$exon_end
      # Check if either donor is contained in the other exon
      donor_a_in_b <- (donor_a >= exon_b$exon_start && donor_a <= exon_b$exon_end)
      donor_b_in_a <- (donor_b >= exon_a$exon_start && donor_b <= exon_a$exon_end)
      return(donor_a_in_b || donor_b_in_a)
    } else if (is_last) {
      # Last exon on plus: acceptor is exon_start (internal boundary)
      acceptor_a <- exon_a$exon_start
      acceptor_b <- exon_b$exon_start
      # Check if either acceptor is contained in the other exon
      acceptor_a_in_b <- (acceptor_a >= exon_b$exon_start && acceptor_a <= exon_b$exon_end)
      acceptor_b_in_a <- (acceptor_b >= exon_a$exon_start && acceptor_b <= exon_a$exon_end)
      return(acceptor_a_in_b || acceptor_b_in_a)
    }
  } else {
    # Minus strand
    if (is_first) {
      # First exon on minus: donor is exon_start (internal boundary)
      donor_a <- exon_a$exon_start
      donor_b <- exon_b$exon_start
      donor_a_in_b <- (donor_a >= exon_b$exon_start && donor_a <= exon_b$exon_end)
      donor_b_in_a <- (donor_b >= exon_a$exon_start && donor_b <= exon_a$exon_end)
      return(donor_a_in_b || donor_b_in_a)
    } else if (is_last) {
      # Last exon on minus: acceptor is exon_end (internal boundary)
      acceptor_a <- exon_a$exon_end
      acceptor_b <- exon_b$exon_end
      acceptor_a_in_b <- (acceptor_a >= exon_b$exon_start && acceptor_a <= exon_b$exon_end)
      acceptor_b_in_a <- (acceptor_b >= exon_a$exon_start && acceptor_b <= exon_a$exon_end)
      return(acceptor_a_in_b || acceptor_b_in_a)
    }
  }

  return(FALSE)
}

#' Match exons one-to-one using greedy algorithm
#'
#' Different rules for terminal vs internal exons:
#' - Terminal exons: match if splice site contained in both
#' - Internal exons: match if ≥50% overlap
#'
#' @param exons_a Exons from isoform A (tibble with exon_start, exon_end)
#' @param exons_b Exons from isoform B (tibble with exon_start, exon_end)
#' @param strand Strand
#' @param exclude_a Indices of exons in A to exclude (e.g., IR cases)
#' @param exclude_b Indices of exons in B to exclude (e.g., IR cases)
#' @return Tibble with columns: idx_a, idx_b, overlap_pct
match_exons_one_to_one <- function(exons_a, exons_b, strand,
                                    exclude_a = integer(),
                                    exclude_b = integer()) {
  # Calculate all potential matches
  matches <- tibble::tibble(
    idx_a = integer(),
    idx_b = integer(),
    overlap_pct = numeric()
  )

  for (i in seq_len(nrow(exons_a))) {
    if (i %in% exclude_a) next

    is_first_a <- (i == 1)
    is_last_a <- (i == nrow(exons_a))

    for (j in seq_len(nrow(exons_b))) {
      if (j %in% exclude_b) next

      is_first_b <- (j == 1)
      is_last_b <- (j == nrow(exons_b))

      # Check if these could be terminal exons
      is_both_first <- is_first_a && is_first_b
      is_both_last <- is_last_a && is_last_b

      match_found <- FALSE

      # Terminal exon matching rule
      if (is_both_first || is_both_last) {
        match_found <- terminal_exons_match(
          exons_a[i, ], exons_b[j, ], strand,
          is_first = is_both_first,
          is_last = is_both_last
        )
      }

      # If not terminal or didn't match by terminal rule, try overlap rule
      if (!match_found && !is_both_first && !is_both_last) {
        overlap <- calculate_overlap_proportion(exons_a[i, ], exons_b[j, ])
        match_found <- (overlap >= OVERLAP_THRESHOLD)
      }

      if (match_found) {
        # Use 1.0 for terminal matches to prioritize them
        overlap_score <- if (is_both_first || is_both_last) {
          1.0
        } else {
          calculate_overlap_proportion(exons_a[i, ], exons_b[j, ])
        }

        matches <- matches %>%
          tibble::add_row(idx_a = i, idx_b = j, overlap_pct = overlap_score)
      }
    }
  }

  # Greedy matching: pick best overlaps first, one-to-one only
  matched <- tibble::tibble(
    idx_a = integer(),
    idx_b = integer(),
    overlap_pct = numeric()
  )

  matched_a <- integer()
  matched_b <- integer()

  # Sort by overlap percentage (descending)
  matches <- matches %>% dplyr::arrange(desc(overlap_pct))

  for (i in seq_len(nrow(matches))) {
    idx_a <- matches$idx_a[i]
    idx_b <- matches$idx_b[i]

    # Check if already matched
    if (idx_a %in% matched_a || idx_b %in% matched_b) next

    # Add to matched pairs
    matched <- matched %>%
      tibble::add_row(
        idx_a = idx_a,
        idx_b = idx_b,
        overlap_pct = matches$overlap_pct[i]
      )

    matched_a <- c(matched_a, idx_a)
    matched_b <- c(matched_b, idx_b)
  }

  return(matched)
}

# ==============================================================================
# Intron Comparison and Classification
# ==============================================================================

#' Compare intronic boundaries for all matched exon pairs
#'
#' Key change: Compare ALL matched pairs, not just consecutive ones
#' The matching established they're "the same" exon; now we explain boundary differences
#'
#' @param exons_a Exons from isoform A (tibble)
#' @param exons_b Exons from isoform B (tibble)
#' @param matched Matched exon pairs (tibble with idx_a, idx_b)
#' @param strand Strand ("+" or "-")
#' @return Tibble with boundary comparisons and event classifications
compare_introns_for_matched_exons <- function(exons_a, exons_b, matched, strand) {
  if (nrow(matched) == 0) {
    return(tibble::tibble(
      idx_a = integer(),
      idx_b = integer(),
      donor_a = integer(),
      acceptor_a = integer(),
      donor_b = integer(),
      acceptor_b = integer(),
      donor_diff = integer(),
      acceptor_diff = integer(),
      event = character()
    ))
  }

  results <- list()

  # For each matched exon pair, compare their intronic boundaries
  for (i in seq_len(nrow(matched))) {
    idx_a <- matched$idx_a[i]
    idx_b <- matched$idx_b[i]

    exon_a <- exons_a[idx_a, ]
    exon_b <- exons_b[idx_b, ]

    # Determine exon types
    is_first_a <- (idx_a == 1)
    is_last_a <- (idx_a == nrow(exons_a))
    is_first_b <- (idx_b == 1)
    is_last_b <- (idx_b == nrow(exons_b))
    is_both_first <- is_first_a && is_first_b
    is_both_last <- is_last_a && is_last_b

    # Define boundaries based on strand
    if (strand == "+") {
      acceptor_a <- exon_a$exon_start
      donor_a <- exon_a$exon_end
      acceptor_b <- exon_b$exon_start
      donor_b <- exon_b$exon_end
    } else {
      donor_a <- exon_a$exon_start
      acceptor_a <- exon_a$exon_end
      donor_b <- exon_b$exon_start
      acceptor_b <- exon_b$exon_end
    }

    # ===== TERMINAL EXON LOGIC =====
    if (is_both_first || is_both_last) {
      # Identify internal vs terminal boundary
      if (is_both_first) {
        internal_a <- donor_a
        internal_b <- donor_b
        terminal_a <- acceptor_a  # TSS
        terminal_b <- acceptor_b  # TSS
        internal_type <- "donor"
      } else {  # is_both_last
        internal_a <- acceptor_a
        internal_b <- acceptor_b
        terminal_a <- donor_a  # TES
        terminal_b <- donor_b  # TES
        internal_type <- "acceptor"
      }

      # Check containment: is internal boundary contained in both exons?
      internal_a_in_b <- (internal_a >= exon_b$exon_start && internal_a <= exon_b$exon_end)
      internal_b_in_a <- (internal_b >= exon_a$exon_start && internal_b <= exon_a$exon_end)

      if (!internal_a_in_b && !internal_b_in_a) {
        # Internal boundaries don't overlap → different intronic regions
        next
      }

      # Compare ONLY internal boundary for splice site events
      internal_diff <- abs(internal_a - internal_b)

      event <- NULL

      # Only compare INTERNAL boundary for intronic events
      # Terminal boundary differences are handled separately (Alt_TSS/Alt_TES)
      if (internal_diff > 0) {
        # Internal boundary differs → intronic event
        if (internal_diff < SPLICE_SITE_THRESHOLD) {
          event <- if (internal_type == "donor") "A5SS" else "A3SS"
        } else {
          # Partial_IR with direction: which end is shared?
          # internal_type indicates which boundary differs
          # donor differs → acceptor (5') shared, donor (3') extends → Partial_IR_5
          # acceptor differs → donor (3') shared, acceptor (5') extends → Partial_IR_3
          event <- if (internal_type == "donor") "Partial_IR_5" else "Partial_IR_3"
        }
      }
      # If internal boundary shared (internal_diff == 0) → no intronic event

      if (!is.null(event)) {
        results[[length(results) + 1]] <- tibble::tibble(
          idx_a = idx_a,
          idx_b = idx_b,
          donor_a = donor_a,
          acceptor_a = acceptor_a,
          donor_b = donor_b,
          acceptor_b = acceptor_b,
          donor_diff = if (internal_type == "donor") internal_diff else 0,
          acceptor_diff = if (internal_type == "acceptor") internal_diff else 0,
          event = event
        )
      }

    } else {
      # ===== INTERNAL EXON LOGIC =====
      # Both boundaries are intronic - compare both
      donor_diff <- abs(donor_a - donor_b)
      acceptor_diff <- abs(acceptor_a - acceptor_b)

      if (donor_diff > 0 || acceptor_diff > 0) {
        event <- classify_intron_difference(donor_diff, acceptor_diff, strand)

        if (event != "none") {
          results[[length(results) + 1]] <- tibble::tibble(
            idx_a = idx_a,
            idx_b = idx_b,
            donor_a = donor_a,
            acceptor_a = acceptor_a,
            donor_b = donor_b,
            acceptor_b = acceptor_b,
            donor_diff = donor_diff,
            acceptor_diff = acceptor_diff,
            event = event
          )
        }
      }
    }
  }

  if (length(results) == 0) {
    return(tibble::tibble(
      idx_a = integer(),
      idx_b = integer(),
      donor_a = integer(),
      acceptor_a = integer(),
      donor_b = integer(),
      acceptor_b = integer(),
      donor_diff = integer(),
      acceptor_diff = integer(),
      event = character()
    ))
  }

  dplyr::bind_rows(results)
}

#' Classify intron difference as A5SS, A3SS, Partial_IR_5, Partial_IR_3, or none
#'
#' @param donor_diff Difference in donor coordinates (bp)
#' @param acceptor_diff Difference in acceptor coordinates (bp)
#' @param strand Strand ("+" or "-")
#' @return Character - event type
#' @details
#' Partial_IR direction labels:
#'   - Partial_IR_5: 5' end shared, 3' end extends (acceptor differs)
#'   - Partial_IR_3: 3' end shared, 5' end extends (donor differs)
classify_intron_difference <- function(donor_diff, acceptor_diff, strand) {
  # No difference
  if (donor_diff == 0 && acceptor_diff == 0) {
    return("none")
  }

  # One boundary differs
  if (donor_diff > 0 && acceptor_diff == 0) {
    # Donor differs, acceptor shared → acceptor (5') shared, donor (3') extends → Partial_IR_5
    if (donor_diff < SPLICE_SITE_THRESHOLD) {
      return("A5SS")
    } else {
      return("Partial_IR_5")
    }
  }

  if (acceptor_diff > 0 && donor_diff == 0) {
    # Acceptor differs, donor shared → donor (3') shared, acceptor (5') extends → Partial_IR_3
    if (acceptor_diff < SPLICE_SITE_THRESHOLD) {
      return("A3SS")
    } else {
      return("Partial_IR_3")
    }
  }

  # Both boundaries differ
  if (donor_diff > 0 && acceptor_diff > 0) {
    # Check if either is large enough for Partial_IR
    if (donor_diff >= SPLICE_SITE_THRESHOLD || acceptor_diff >= SPLICE_SITE_THRESHOLD) {
      # Determine primary direction based on which diff is larger
      if (donor_diff > acceptor_diff) {
        return("Partial_IR_5")  # Donor change is larger (donor/3' extends more)
      } else {
        return("Partial_IR_3")  # Acceptor change is larger (acceptor/5' extends more)
      }
    } else {
      # Both are small differences - could be both A5SS and A3SS
      # For now, classify as "Dual_splice" to indicate both boundaries changed
      return("Dual_splice")
    }
  }

  return("none")
}

# ==============================================================================
# TSS/TES Detection (Terminal Boundaries)
# ==============================================================================

#' Detect TSS changes between first exons
#'
#' @param exon_a First exon from isoform A
#' @param exon_b First exon from isoform B
#' @param strand Gene strand
#' @param tolerance TSS tolerance in bp
#' @return Logical - TRUE if TSS changed
detect_tss_change <- function(exon_a, exon_b, strand, tolerance = TSS_TOLERANCE) {
  if (strand == "+") {
    tss_a <- exon_a$exon_start
    tss_b <- exon_b$exon_start
  } else {
    tss_a <- exon_a$exon_end
    tss_b <- exon_b$exon_end
  }

  tss_diff <- abs(tss_a - tss_b)
  return(tss_diff > tolerance)
}

#' Detect TES changes between last exons
#'
#' @param exon_a Last exon from isoform A
#' @param exon_b Last exon from isoform B
#' @param strand Gene strand
#' @param tolerance TES tolerance in bp
#' @return Logical - TRUE if TES changed
detect_tes_change <- function(exon_a, exon_b, strand, tolerance = TES_TOLERANCE) {
  if (strand == "+") {
    tes_a <- exon_a$exon_end
    tes_b <- exon_b$exon_end
  } else {
    tes_a <- exon_a$exon_start
    tes_b <- exon_b$exon_start
  }

  tes_diff <- abs(tes_a - tes_b)
  return(tes_diff > tolerance)
}

# ==============================================================================
# Main Detection Function
# ==============================================================================

#' Detect splicing events between two isoforms (V2 - intron-centric)
#'
#' @param exons_a Exons from isoform A (tibble with exon_start, exon_end, sorted by exon_number)
#' @param exons_b Exons from isoform B (tibble with exon_start, exon_end, sorted by exon_number)
#' @param strand Strand ("+" or "-")
#' @return List with event counts and details
detect_splicing_events_v2 <- function(exons_a, exons_b, strand) {

  # Initialize results
  events <- list(
    alt_tss = FALSE,
    alt_tes = FALSE,
    n_ir = 0,
    n_a5ss = 0,
    n_a3ss = 0,
    n_partial_ir = 0,
    n_partial_ir_5 = 0,
    n_partial_ir_3 = 0,
    n_se = 0,
    ir_details = list(),
    intron_details = tibble::tibble(),
    matched_exons = tibble::tibble()
  )

  # STEP 1: Terminal boundaries (TSS/TES)
  # IMPORTANT: Exons must be in transcriptional order (5' to 3')
  # This means exon 1 (row 1) is TSS exon for BOTH plus and minus strands
  # And last exon (row n) is TES exon for BOTH plus and minus strands
  if (nrow(exons_a) > 0 && nrow(exons_b) > 0) {
    events$alt_tss <- detect_tss_change(exons_a[1, ], exons_b[1, ], strand)
    events$alt_tes <- detect_tes_change(
      exons_a[nrow(exons_a), ],
      exons_b[nrow(exons_b), ],
      strand
    )
  }

  # STEP 2: Identify IR (ANY overlap with multiple exons)
  ir_results <- detect_ir_v2(exons_a, exons_b)
  events$n_ir <- length(ir_results$ir_in_a) + length(ir_results$ir_in_b)
  events$ir_details <- ir_results

  # STEP 3: Match remaining exons one-to-one (exclude IR cases)
  matched <- match_exons_one_to_one(
    exons_a, exons_b, strand,
    exclude_a = ir_results$ir_in_a,
    exclude_b = ir_results$ir_in_b
  )
  events$matched_exons <- matched

  # STEP 4: Compare introns between matched exon pairs
  # CRITICAL: Skip for monoexonic comparisons (no internal splice sites)
  # A monoexonic isoform has boundaries that are TSS/TES, not splice sites
  # So Partial_IR/A5SS/A3SS detection requires both isoforms to have ≥2 exons
  if (nrow(exons_a) >= 2 && nrow(exons_b) >= 2) {
    intron_comparisons <- compare_introns_for_matched_exons(
      exons_a, exons_b, matched, strand
    )
    events$intron_details <- intron_comparisons

    # Count events from intron comparisons
    if (nrow(intron_comparisons) > 0) {
      events$n_a5ss <- sum(intron_comparisons$event == "A5SS")
      events$n_a3ss <- sum(intron_comparisons$event == "A3SS")
      events$n_partial_ir <- sum(grepl("^Partial_IR", intron_comparisons$event))
      events$n_partial_ir_5 <- sum(intron_comparisons$event == "Partial_IR_5")
      events$n_partial_ir_3 <- sum(intron_comparisons$event == "Partial_IR_3")
    }
  } else {
    # Monoexonic comparison: no internal splice sites, so no Partial_IR/A5SS/A3SS
    events$intron_details <- tibble::tibble()
  }

  # STEP 5: Identify skipped exons (SE)
  # Exons that didn't match and aren't involved in IR
  matched_a <- matched$idx_a
  matched_b <- matched$idx_b

  # For IR: need to exclude not just the spanning exon, but also the spanned exons
  # If exon in B spans multiple exons in A, those A exons are involved in IR (not SE)
  # If exon in A spans multiple exons in B, those B exons are involved in IR (not SE)

  ir_involved_a <- ir_results$ir_in_a  # Exons in A that span multiple B exons
  ir_involved_b <- ir_results$ir_in_b  # Exons in B that span multiple A exons

  # Also exclude exons that are spanned by IR exons
  for (ir_idx in ir_results$ir_in_b) {
    # This B exon spans multiple A exons - find which A exons it overlaps
    for (i in seq_len(nrow(exons_a))) {
      if (has_any_overlap(exons_b[ir_idx, ], exons_a[i, ])) {
        ir_involved_a <- c(ir_involved_a, i)
      }
    }
  }

  for (ir_idx in ir_results$ir_in_a) {
    # This A exon spans multiple B exons - find which B exons it overlaps
    for (j in seq_len(nrow(exons_b))) {
      if (has_any_overlap(exons_a[ir_idx, ], exons_b[j, ])) {
        ir_involved_b <- c(ir_involved_b, j)
      }
    }
  }

  ir_involved_a <- unique(ir_involved_a)
  ir_involved_b <- unique(ir_involved_b)

  # Additional check: IR exons may also share boundaries with spanned exons
  # This detects Partial_IR or A5SS/A3SS in addition to IR
  # BUT: skip for monoexonic comparisons (boundaries are TSS/TES, not splice sites)
  # CRITICAL: Both isoforms must have ≥2 exons for Partial_IR/A5SS/A3SS at IR boundaries
  # If either is monoexonic, boundary differences are Alt_TSS/Alt_TES only
  if (nrow(exons_a) >= 2 && nrow(exons_b) >= 2) {
    for (ir_idx in ir_results$ir_in_b) {
    # B exon spans multiple A exons - check if it shares boundaries
    ir_exon <- exons_b[ir_idx, ]

    # Skip if IR exon is monoexonic (first AND last exon)
    # Monoexonic boundaries are TSS/TES, not internal splice sites
    is_monoexonic_b <- (nrow(exons_b) == 1)
    if (is_monoexonic_b) next

    for (i in seq_len(nrow(exons_a))) {
      if (!has_any_overlap(ir_exon, exons_a[i, ])) next

      exon_a <- exons_a[i, ]

      # Check boundary sharing (plus strand: acceptor=start, donor=end)
      if (strand == "+") {
        acceptor_shared <- (ir_exon$exon_start == exon_a$exon_start)
        donor_shared <- (ir_exon$exon_end == exon_a$exon_end)

        if (acceptor_shared && !donor_shared) {
          # Shared acceptor, donor differs
          donor_diff <- abs(ir_exon$exon_end - exon_a$exon_end)
          if (donor_diff < SPLICE_SITE_THRESHOLD) {
            events$n_a5ss <- events$n_a5ss + 1
          } else {
            events$n_partial_ir <- events$n_partial_ir + 1
          }
        } else if (donor_shared && !acceptor_shared) {
          # Shared donor, acceptor differs
          acceptor_diff <- abs(ir_exon$exon_start - exon_a$exon_start)
          if (acceptor_diff < SPLICE_SITE_THRESHOLD) {
            events$n_a3ss <- events$n_a3ss + 1
          } else {
            events$n_partial_ir <- events$n_partial_ir + 1
          }
        }
      } else {  # minus strand: acceptor=end, donor=start
        acceptor_shared <- (ir_exon$exon_end == exon_a$exon_end)
        donor_shared <- (ir_exon$exon_start == exon_a$exon_start)

        if (acceptor_shared && !donor_shared) {
          # Shared acceptor, donor differs
          donor_diff <- abs(ir_exon$exon_start - exon_a$exon_start)
          if (donor_diff < SPLICE_SITE_THRESHOLD) {
            events$n_a5ss <- events$n_a5ss + 1
          } else {
            events$n_partial_ir <- events$n_partial_ir + 1
          }
        } else if (donor_shared && !acceptor_shared) {
          # Shared donor, acceptor differs
          acceptor_diff <- abs(ir_exon$exon_end - exon_a$exon_end)
          if (acceptor_diff < SPLICE_SITE_THRESHOLD) {
            events$n_a3ss <- events$n_a3ss + 1
          } else {
            events$n_partial_ir <- events$n_partial_ir + 1
          }
        }
      }
    }
  }

  # Similarly for A exons that span multiple B exons
  for (ir_idx in ir_results$ir_in_a) {
    # A exon spans multiple B exons - check if it shares boundaries
    ir_exon <- exons_a[ir_idx, ]

    # Skip if IR exon is monoexonic (first AND last exon)
    # Monoexonic boundaries are TSS/TES, not internal splice sites
    is_monoexonic_a <- (nrow(exons_a) == 1)
    if (is_monoexonic_a) next

    for (j in seq_len(nrow(exons_b))) {
      if (!has_any_overlap(ir_exon, exons_b[j, ])) next

      exon_b <- exons_b[j, ]

      # Check boundary sharing (plus strand: acceptor=start, donor=end)
      if (strand == "+") {
        acceptor_shared <- (ir_exon$exon_start == exon_b$exon_start)
        donor_shared <- (ir_exon$exon_end == exon_b$exon_end)

        if (acceptor_shared && !donor_shared) {
          # Shared acceptor, donor differs
          donor_diff <- abs(ir_exon$exon_end - exon_b$exon_end)
          if (donor_diff < SPLICE_SITE_THRESHOLD) {
            events$n_a5ss <- events$n_a5ss + 1
          } else {
            events$n_partial_ir <- events$n_partial_ir + 1
          }
        } else if (donor_shared && !acceptor_shared) {
          # Shared donor, acceptor differs
          acceptor_diff <- abs(ir_exon$exon_start - exon_b$exon_start)
          if (acceptor_diff < SPLICE_SITE_THRESHOLD) {
            events$n_a3ss <- events$n_a3ss + 1
          } else {
            events$n_partial_ir <- events$n_partial_ir + 1
          }
        }
      } else {  # minus strand: acceptor=end, donor=start
        acceptor_shared <- (ir_exon$exon_end == exon_b$exon_end)
        donor_shared <- (ir_exon$exon_start == exon_b$exon_start)

        if (acceptor_shared && !donor_shared) {
          # Shared acceptor, donor differs
          donor_diff <- abs(ir_exon$exon_start - exon_b$exon_start)
          if (donor_diff < SPLICE_SITE_THRESHOLD) {
            events$n_a5ss <- events$n_a5ss + 1
          } else {
            events$n_partial_ir <- events$n_partial_ir + 1
          }
        } else if (donor_shared && !acceptor_shared) {
          # Shared donor, acceptor differs
          acceptor_diff <- abs(ir_exon$exon_end - exon_b$exon_end)
          if (acceptor_diff < SPLICE_SITE_THRESHOLD) {
            events$n_a3ss <- events$n_a3ss + 1
          } else {
            events$n_partial_ir <- events$n_partial_ir + 1
          }
        }
      }
    }
  }
  }  # End of monoexonic check for Partial_IR/A5SS/A3SS at IR boundaries

  unmatched_a <- setdiff(seq_len(nrow(exons_a)),
                         c(matched_a, ir_involved_a))
  unmatched_b <- setdiff(seq_len(nrow(exons_b)),
                         c(matched_b, ir_involved_b))

  # SE count: unmatched exons (excluding first/last which are terminal)
  # CRITICAL: SE requires splicing on both sides, so both isoforms must have ≥2 exons
  # A monoexonic isoform cannot provide flanking context for SE detection
  if (nrow(exons_a) >= 2 && nrow(exons_b) >= 2) {
    se_a <- unmatched_a[unmatched_a > 1 & unmatched_a < nrow(exons_a)]
    se_b <- unmatched_b[unmatched_b > 1 & unmatched_b < nrow(exons_b)]
    events$n_se <- length(se_a) + length(se_b)
  } else {
    # If either isoform is monoexonic, SE cannot occur (no flanking exons)
    events$n_se <- 0
  }

  return(events)
}

# ==============================================================================
# Utility Functions
# ==============================================================================

#' Check if exon is a micro-exon (≤30bp)
#'
#' @param exon Exon (tibble row with exon_start, exon_end)
#' @return Logical - TRUE if micro-exon
is_micro_exon <- function(exon) {
  exon_length <- exon$exon_end - exon$exon_start + 1
  return(exon_length <= MICRO_EXON_THRESHOLD)
}

#' Summarize events as string
#'
#' @param events Event detection results from detect_splicing_events_v2
#' @return Character - comma-separated event list
summarize_events_v2 <- function(events) {
  event_list <- character()

  # Standard order: TSS/TES, SE, A5SS/A3SS, Partial_IR/IR
  if (events$alt_tss) event_list <- c(event_list, "Alt_TSS")
  if (events$alt_tes) event_list <- c(event_list, "Alt_TES")
  if (events$n_se > 0) event_list <- c(event_list, "SE")
  if (events$n_a5ss > 0) event_list <- c(event_list, "A5SS")
  if (events$n_a3ss > 0) event_list <- c(event_list, "A3SS")
  if (events$n_partial_ir_5 > 0) event_list <- c(event_list, "Partial_IR_5")
  if (events$n_partial_ir_3 > 0) event_list <- c(event_list, "Partial_IR_3")
  if (events$n_ir > 0) event_list <- c(event_list, "IR")

  if (length(event_list) == 0) return("none")
  return(paste(event_list, collapse = ","))
}
