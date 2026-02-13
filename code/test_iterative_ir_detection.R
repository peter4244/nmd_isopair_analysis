#!/usr/bin/env Rscript
# Test Iterative IR Detection Algorithm
#
# Purpose: Validate corrected iterative detection on specific test genes
#
# Test cases: Genes with existing A5SS/A3SS events (shared splice sites)
# - ENSG00000113575: A3SS + Alt_TSS (4 union exons)
# - ENSG00000172468: Has A5SS/A3SS (4 union exons)
# - ENSG00000180818: Has A5SS/A3SS (4 union exons)

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   TEST ITERATIVE IR DETECTION                                 ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Configuration
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v4.0_reference_based/major_isoforms"
setwd(base_dir)

TSS_TES_THRESHOLD <- 20

# Test genes (genes with A5SS/A3SS events - shared splice sites)
TEST_GENES <- c(
  "ENSG00000113575",  # A3SS + Alt_TSS (4 union exons)
  "ENSG00000172468",  # Has A5SS/A3SS (4 union exons)
  "ENSG00000180818"   # Has A5SS/A3SS (4 union exons)
)

# ============================================================================
# Helper Functions (from detect_events_major_isoforms.R)
# ============================================================================

get_isoform_exon_details <- function(union_exons, isoform_id) {
  exon_details <- list()

  for (i in seq_along(union_exons)) {
    exon_group <- union_exons[[i]]
    iso_row <- exon_group %>% filter(isoform_id == !!isoform_id)

    if (nrow(iso_row) > 0) {
      exon_details[[length(exon_details) + 1]] <- list(
        union_exon_number = i,
        exon_type = iso_row$exon_type[1],
        start = iso_row$start[1],
        end = iso_row$end[1],
        is_first = iso_row$is_first[1],
        is_last = iso_row$is_last[1],
        transcript_exon_number = iso_row$transcript_exon_number[1]
      )
    }
  }

  return(bind_rows(exon_details))
}

extract_junctions <- function(exon_details, gene_strand) {
  if (nrow(exon_details) < 2) {
    return(tibble(
      donor = integer(),
      acceptor = integer(),
      exon_upstream = integer(),
      exon_downstream = integer()
    ))
  }

  # Sort by transcript order
  exons_sorted <- exon_details %>% arrange(transcript_exon_number)

  junctions <- list()

  for (i in 1:(nrow(exons_sorted) - 1)) {
    exon_current <- exons_sorted[i, ]
    exon_next <- exons_sorted[i + 1, ]

    if (gene_strand == "+") {
      donor <- exon_current$end
      acceptor <- exon_next$start
    } else {
      donor <- exon_next$start
      acceptor <- exon_current$end
    }

    junctions[[length(junctions) + 1]] <- list(
      donor = donor,
      acceptor = acceptor,
      exon_upstream = exon_current$union_exon_number,
      exon_downstream = exon_next$union_exon_number
    )
  }

  return(bind_rows(junctions))
}

# ============================================================================
# New Functions for Iterative Detection
# ============================================================================

#' Generate unique signature for each event
generate_event_signature <- function(event) {
  if (event$event_type == "A5SS") {
    paste0("A5SS:", min(event$donor_A, event$donor_B), ":",
           max(event$donor_A, event$donor_B), ":", event$acceptor)
  } else if (event$event_type == "A3SS") {
    paste0("A3SS:", event$donor, ":", min(event$acceptor_A, event$acceptor_B), ":",
           max(event$acceptor_A, event$acceptor_B))
  } else if (grepl("^IR", event$event_type)) {
    # Canonical IR: always from perspective of retaining isoform
    paste0(event$event_type, ":", event$retained_start, ":", event$retained_end,
           ":", event$retaining_isoform)
  } else {
    paste0(event$event_type, ":", event$union_exon_id)
  }
}

#' Extract splice site coordinates from events
extract_splice_site_coords <- function(splice_events) {
  coords_list <- list()

  for (i in seq_along(splice_events)) {
    event <- splice_events[[i]]
    if (event$event_type == "A5SS") {
      coords_list[[length(coords_list) + 1]] <- tibble(
        type = "A5SS",
        donor_A = event$donor_A,
        donor_B = event$donor_B,
        acceptor_A = event$acceptor,
        acceptor_B = event$acceptor
      )
    } else if (event$event_type == "A3SS") {
      coords_list[[length(coords_list) + 1]] <- tibble(
        type = "A3SS",
        donor_A = event$donor,
        donor_B = event$donor,
        acceptor_A = event$acceptor_A,
        acceptor_B = event$acceptor_B
      )
    }
  }

  if (length(coords_list) > 0) {
    return(bind_rows(coords_list))
  } else {
    return(tibble(type = character(), donor_A = integer(), donor_B = integer(),
                  acceptor_A = integer(), acceptor_B = integer()))
  }
}

#' Detect A5SS and A3SS with iteration tracking
detect_splice_sites_with_flags <- function(junctions_A, junctions_B, exons_A, exons_B, iteration) {
  events <- list()

  if (nrow(junctions_A) == 0 || nrow(junctions_B) == 0) {
    return(events)
  }

  # Compare all junction pairs
  for (i in seq_len(nrow(junctions_A))) {
    junc_A <- junctions_A[i, ]

    for (j in seq_len(nrow(junctions_B))) {
      junc_B <- junctions_B[j, ]

      # Check for alternative splice sites
      if (junc_A$donor == junc_B$donor && junc_A$acceptor != junc_B$acceptor) {
        # A3SS: Same donor, different acceptor
        events[[length(events) + 1]] <- tibble(
          event_type = "A3SS",
          donor = junc_A$donor,
          acceptor_A = junc_A$acceptor,
          acceptor_B = junc_B$acceptor,
          iteration = iteration
        )
      } else if (junc_A$donor != junc_B$donor && junc_A$acceptor == junc_B$acceptor) {
        # A5SS: Different donor, same acceptor
        events[[length(events) + 1]] <- tibble(
          event_type = "A5SS",
          donor_A = junc_A$donor,
          donor_B = junc_B$donor,
          acceptor = junc_A$acceptor,
          iteration = iteration
        )
      }
    }
  }

  return(events)
}

#' Update flags based on detected splice events
update_flags_from_splice_events <- function(splice_events, exons_A, exons_B) {
  for (event in splice_events) {
    if (event$event_type == "A5SS") {
      # Mark exons with these donor coordinates
      exons_A <- exons_A %>% mutate(has_a5ss = if_else(end == event$donor_A, TRUE, has_a5ss))
      exons_B <- exons_B %>% mutate(has_a5ss = if_else(end == event$donor_B, TRUE, has_a5ss))
    } else if (event$event_type == "A3SS") {
      # Mark exons with these acceptor coordinates
      exons_A <- exons_A %>% mutate(has_a3ss = if_else(start == event$acceptor_A, TRUE, has_a3ss))
      exons_B <- exons_B %>% mutate(has_a3ss = if_else(start == event$acceptor_B, TRUE, has_a3ss))
    }
  }

  return(list(exons_A = exons_A, exons_B = exons_B))
}

#' Detect internal IR with flag-informed boundary matching
detect_internal_ir_with_flags <- function(exons_A, exons_B, gene_strand, iteration, known_splice_sites) {
  events <- list()

  # Internal exons only
  internal_A <- exons_A %>% filter(!is_first, !is_last)
  internal_B <- exons_B %>% filter(!is_first, !is_last)

  if (nrow(internal_A) == 0 || nrow(internal_B) < 2) {
    return(bind_rows(events))
  }

  # Check if exon in A spans consecutive exons in B
  for (i in seq_len(nrow(internal_A))) {
    exon_A <- internal_A[i, ]

    for (j in 1:(nrow(internal_B) - 1)) {
      exon_B_first <- internal_B[j, ]
      exon_B_last <- internal_B[j + 1, ]

      # Must be consecutive in B
      if (exon_B_last$transcript_exon_number != exon_B_first$transcript_exon_number + 1) {
        next
      }

      # Check if A encompasses both B exons
      encompasses <- (exon_A$start <= exon_B_first$start && exon_A$end >= exon_B_last$end)
      if (!encompasses) next

      # PRECISE FLAG-INFORMED BOUNDARY MATCHING
      # Start boundary: exact match OR explained by known A3SS
      start_matches <- if (exon_A$start == exon_B_first$start) {
        TRUE
      } else if (exon_A$has_a3ss || exon_B_first$has_a3ss) {
        # Verify this difference is explained by a known A3SS
        if (nrow(known_splice_sites) > 0) {
          any(known_splice_sites$type == "A3SS" &
              ((known_splice_sites$acceptor_A == exon_A$start & known_splice_sites$acceptor_B == exon_B_first$start) |
               (known_splice_sites$acceptor_A == exon_B_first$start & known_splice_sites$acceptor_B == exon_A$start)))
        } else {
          FALSE
        }
      } else {
        FALSE
      }

      # End boundary: exact match OR explained by known A5SS
      end_matches <- if (exon_A$end == exon_B_last$end) {
        TRUE
      } else if (exon_A$has_a5ss || exon_B_last$has_a5ss) {
        # Verify this difference is explained by a known A5SS
        if (nrow(known_splice_sites) > 0) {
          any(known_splice_sites$type == "A5SS" &
              ((known_splice_sites$donor_A == exon_A$end & known_splice_sites$donor_B == exon_B_last$end) |
               (known_splice_sites$donor_A == exon_B_last$end & known_splice_sites$donor_B == exon_A$end)))
        } else {
          FALSE
        }
      } else {
        FALSE
      }

      if (start_matches && end_matches) {
        # Determine compound event type
        has_a3ss_diff <- exon_A$start != exon_B_first$start
        has_a5ss_diff <- exon_A$end != exon_B_last$end

        event_type <- if (has_a3ss_diff && has_a5ss_diff) {
          "IR+A3SS+A5SS"
        } else if (has_a3ss_diff) {
          "IR+A3SS"
        } else if (has_a5ss_diff) {
          "IR+A5SS"
        } else {
          "IR"
        }

        # CANONICAL FORM: Always report from perspective of retaining isoform
        events[[length(events) + 1]] <- tibble(
          event_type = event_type,
          retaining_isoform = "A",
          retained_start = exon_A$start,
          retained_end = exon_A$end,
          split_start = exon_B_first$start,
          split_end = exon_B_last$end,
          iteration = iteration
        )
      }
    }
  }

  # SYMMETRIC CHECK: Now check if exons in B span exons in A
  for (i in seq_len(nrow(internal_B))) {
    exon_B <- internal_B[i, ]

    for (j in 1:(nrow(internal_A) - 1)) {
      exon_A_first <- internal_A[j, ]
      exon_A_last <- internal_A[j + 1, ]

      # Must be consecutive in A
      if (exon_A_last$transcript_exon_number != exon_A_first$transcript_exon_number + 1) {
        next
      }

      encompasses <- (exon_B$start <= exon_A_first$start && exon_B$end >= exon_A_last$end)
      if (!encompasses) next

      # Same precise boundary matching
      start_matches <- if (exon_B$start == exon_A_first$start) {
        TRUE
      } else if (exon_B$has_a3ss || exon_A_first$has_a3ss) {
        if (nrow(known_splice_sites) > 0) {
          any(known_splice_sites$type == "A3SS" &
              ((known_splice_sites$acceptor_A == exon_B$start & known_splice_sites$acceptor_B == exon_A_first$start) |
               (known_splice_sites$acceptor_A == exon_A_first$start & known_splice_sites$acceptor_B == exon_B$start)))
        } else {
          FALSE
        }
      } else {
        FALSE
      }

      end_matches <- if (exon_B$end == exon_A_last$end) {
        TRUE
      } else if (exon_B$has_a5ss || exon_A_last$has_a5ss) {
        if (nrow(known_splice_sites) > 0) {
          any(known_splice_sites$type == "A5SS" &
              ((known_splice_sites$donor_A == exon_B$end & known_splice_sites$donor_B == exon_A_last$end) |
               (known_splice_sites$donor_A == exon_A_last$end & known_splice_sites$donor_B == exon_B$end)))
        } else {
          FALSE
        }
      } else {
        FALSE
      }

      if (start_matches && end_matches) {
        has_a3ss_diff <- exon_B$start != exon_A_first$start
        has_a5ss_diff <- exon_B$end != exon_A_last$end

        event_type <- if (has_a3ss_diff && has_a5ss_diff) {
          "IR+A3SS+A5SS"
        } else if (has_a3ss_diff) {
          "IR+A3SS"
        } else if (has_a5ss_diff) {
          "IR+A5SS"
        } else {
          "IR"
        }

        events[[length(events) + 1]] <- tibble(
          event_type = event_type,
          retaining_isoform = "B",
          retained_start = exon_B$start,
          retained_end = exon_B$end,
          split_start = exon_A_first$start,
          split_end = exon_A_last$end,
          iteration = iteration
        )
      }
    }
  }

  return(bind_rows(events))
}

#' Detect monoexonic IR (one isoform is monoexonic)
detect_monoexonic_ir <- function(exons_A, exons_B, gene_strand, iteration) {
  events <- list()

  # Check if either isoform is monoexonic
  is_mono_A <- nrow(exons_A) == 1
  is_mono_B <- nrow(exons_B) == 1

  if (!is_mono_A && !is_mono_B) {
    return(bind_rows(events))
  }

  if (is_mono_A && !is_mono_B) {
    # A is monoexonic, B is multi-exonic
    # Check if A spans all of B
    mono_exon <- exons_A[1, ]
    min_B <- min(exons_B$start)
    max_B <- max(exons_B$end)

    if (mono_exon$start <= min_B && mono_exon$end >= max_B) {
      events[[length(events) + 1]] <- tibble(
        event_type = "IR_monoexonic",
        retaining_isoform = "A",
        retained_start = mono_exon$start,
        retained_end = mono_exon$end,
        split_start = min_B,
        split_end = max_B,
        iteration = iteration
      )
    }
  } else if (is_mono_B && !is_mono_A) {
    # B is monoexonic, A is multi-exonic
    mono_exon <- exons_B[1, ]
    min_A <- min(exons_A$start)
    max_A <- max(exons_A$end)

    if (mono_exon$start <= min_A && mono_exon$end >= max_A) {
      events[[length(events) + 1]] <- tibble(
        event_type = "IR_monoexonic",
        retaining_isoform = "B",
        retained_start = mono_exon$start,
        retained_end = mono_exon$end,
        split_start = min_A,
        split_end = max_A,
        iteration = iteration
      )
    }
  }

  return(bind_rows(events))
}

#' Main iterative detection function
detect_events_iteratively <- function(isoform_A, isoform_B, union_model, gene_strand) {
  # Initialize exons with flags
  exons_A <- get_isoform_exon_details(union_model$union_exons, isoform_A) %>%
    mutate(has_a5ss = FALSE, has_a3ss = FALSE, has_ir = FALSE)

  exons_B <- get_isoform_exon_details(union_model$union_exons, isoform_B) %>%
    mutate(has_a5ss = FALSE, has_a3ss = FALSE, has_ir = FALSE)

  all_events <- list()
  detected_signatures <- character()
  known_splice_sites <- tibble(type = character(), donor_A = integer(), donor_B = integer(),
                                acceptor_A = integer(), acceptor_B = integer())
  iteration <- 0

  cat("\n  Starting iterative detection for", isoform_A, "vs", isoform_B, "\n")

  # Iterative detection loop
  repeat {
    iteration <- iteration + 1
    iteration_events <- list()

    cat("  Iteration", iteration, "...\n")

    # STEP 1: A5SS and A3SS Detection
    junctions_A <- extract_junctions(exons_A, gene_strand)
    junctions_B <- extract_junctions(exons_B, gene_strand)
    splice_events <- detect_splice_sites_with_flags(junctions_A, junctions_B, exons_A, exons_B, iteration)

    if (length(splice_events) > 0) {
      # Generate signatures
      splice_signatures <- sapply(splice_events, function(e) generate_event_signature(e))
      truly_new <- !splice_signatures %in% detected_signatures

      if (any(truly_new)) {
        new_splice_events <- splice_events[truly_new]
        cat("    Found", sum(truly_new), "new A5SS/A3SS events\n")
        iteration_events <- c(iteration_events, new_splice_events)
        detected_signatures <- c(detected_signatures, splice_signatures[truly_new])

        # Update flags
        flag_updates <- update_flags_from_splice_events(new_splice_events, exons_A, exons_B)
        exons_A <- flag_updates$exons_A
        exons_B <- flag_updates$exons_B
        known_splice_sites <- bind_rows(known_splice_sites,
                                        extract_splice_site_coords(new_splice_events))
      }
    }

    # STEP 2: IR Detection
    ir_events_internal <- detect_internal_ir_with_flags(exons_A, exons_B, gene_strand, iteration, known_splice_sites)
    ir_events_monoexonic <- detect_monoexonic_ir(exons_A, exons_B, gene_strand, iteration)

    ir_events <- bind_rows(ir_events_internal, ir_events_monoexonic)

    if (nrow(ir_events) > 0) {
      # Generate signatures
      ir_signatures <- apply(ir_events, 1, function(e) {
        paste0(e$event_type, ":", e$retained_start, ":", e$retained_end, ":", e$retaining_isoform)
      })
      truly_new <- !ir_signatures %in% detected_signatures

      if (any(truly_new)) {
        new_ir_events <- ir_events[truly_new, ]
        cat("    Found", sum(truly_new), "new IR events\n")
        iteration_events <- c(iteration_events, list(new_ir_events))
        detected_signatures <- c(detected_signatures, ir_signatures[truly_new])
      }
    }

    # STEP 3: Check Convergence
    if (length(iteration_events) == 0) {
      cat("  Converged after", iteration, "iterations\n")
      break
    }

    all_events <- c(all_events, iteration_events)

    # Safety check
    if (iteration > 100) {
      warning(paste("Exceeded 100 iterations - possible infinite loop"))
      break
    }
  }

  return(list(events = all_events, iterations = iteration))
}

# ============================================================================
# Main Test
# ============================================================================

cat("Loading data...\n")
union_models <- readRDS("union_exon_models_major.rds")
events_current <- readRDS("isoform_pairs_events.rds")
exon_structures <- readRDS("exon_structures_expanded.rds")

# Create gene-to-strand mapping
gene_strand_map <- exon_structures %>%
  distinct(gene_id, strand) %>%
  mutate(strand = as.character(strand))

cat("\nTesting on", length(TEST_GENES), "genes with alternative splice sites\n")

for (gene_id in TEST_GENES) {
  cat("\n")
  cat("════════════════════════════════════════════════════════════════\n")
  cat("Testing:", gene_id, "\n")
  cat("════════════════════════════════════════════════════════════════\n")

  # Get union model
  model <- union_models[[gene_id]]

  # Get strand from mapping
  gene_strand <- gene_strand_map %>%
    filter(gene_id == !!gene_id) %>%
    pull(strand) %>%
    first()

  cat("Gene name:", model$gene_name, "\n")
  cat("N isoforms:", model$n_isoforms, "\n")
  cat("N union exons:", model$n_union_exons, "\n")
  cat("Strand:", gene_strand, "\n")

  # Get all isoforms for this gene
  iso_ids <- unique(unlist(lapply(model$union_exons, function(e) e$isoform_id)))

  if (length(iso_ids) < 2) {
    cat("Less than 2 isoforms - skipping\n")
    next
  }

  # Test the first pair
  isoform_A <- iso_ids[1]
  isoform_B <- iso_ids[2]

  cat("\nTesting pair:", isoform_A, "vs", isoform_B, "\n")

  # Show exon structure
  exons_A <- get_isoform_exon_details(model$union_exons, isoform_A)
  exons_B <- get_isoform_exon_details(model$union_exons, isoform_B)

  cat("\nIsoform A structure (", nrow(exons_A), "exons):\n")
  print(exons_A %>% select(union_exon_number, start, end, is_first, is_last, transcript_exon_number))

  cat("\nIsoform B structure (", nrow(exons_B), "exons):\n")
  print(exons_B %>% select(union_exon_number, start, end, is_first, is_last, transcript_exon_number))

  # Show junctions
  junctions_A <- extract_junctions(get_isoform_exon_details(model$union_exons, isoform_A), gene_strand)
  junctions_B <- extract_junctions(get_isoform_exon_details(model$union_exons, isoform_B), gene_strand)

  cat("\nJunctions in A (", nrow(junctions_A), "):\n")
  if (nrow(junctions_A) > 0) print(junctions_A)

  cat("\nJunctions in B (", nrow(junctions_B), "):\n")
  if (nrow(junctions_B) > 0) print(junctions_B)

  # Run iterative detection
  result <- detect_events_iteratively(isoform_A, isoform_B, model, gene_strand)

  cat("\n")
  cat("RESULTS:\n")
  cat("--------\n")
  cat("Total events detected:", length(result$events), "\n")
  cat("Converged in", result$iterations, "iterations\n")

  if (length(result$events) > 0) {
    cat("\nDetected events:\n")
    for (i in seq_along(result$events)) {
      event <- result$events[[i]]
      if (is.data.frame(event)) {
        print(event)
      } else {
        print(as_tibble(event))
      }
    }
  } else {
    cat("No events detected (still NO_DIFF)\n")
  }
}

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("Test complete\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("\n")
