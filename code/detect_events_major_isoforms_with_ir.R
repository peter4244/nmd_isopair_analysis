#!/usr/bin/env Rscript
# Event Detection for Major Isoforms - v5.0 WITH IR DETECTION
#
# Purpose:
#   Detect splicing events including intron retention (IR) using iterative approach
#
# Changes from v4.0:
#   - Added iterative detection for A5SS, A3SS, and IR events
#   - IR detection includes: internal IR, terminal IR, and monoexonic IR
#   - Compound events: IR+A5SS, IR+A3SS, IR+A3SS+A5SS
#   - Event signature tracking to prevent infinite loops
#   - Flag-informed boundary matching for precise IR detection
#
# Inputs:
#   - major_isoforms/union_exon_models_major.rds (3,797 genes)
#   - major_isoforms/exon_structures_expanded.rds (strand info)
#
# Outputs:
#   - major_isoforms/isoform_pairs_events_with_ir.rds
#   - major_isoforms/event_summary_with_ir.tsv

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EVENT DETECTION v5.0 - WITH IR DETECTION                    ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Configuration
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v4.0_reference_based"
output_dir <- file.path(base_dir, "major_isoforms")
BATCH_SIZE <- 100
CHECKPOINT_BATCH <- NA  # Run all genes

TSS_TES_THRESHOLD <- 20

# ============================================================================
# Helper Functions
# ============================================================================

get_isoform_union_exons <- function(union_exons, isoform_id) {
  present_exons <- integer()
  for (i in seq_along(union_exons)) {
    exon_group <- union_exons[[i]]
    if (isoform_id %in% exon_group$isoform_id) {
      present_exons <- c(present_exons, i)
    }
  }
  return(present_exons)
}

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

  exons_sorted <- exon_details %>% arrange(transcript_exon_number)
  junctions <- list()

  for (i in 1:(nrow(exons_sorted) - 1)) {
    exon_current <- exons_sorted[i, ]
    exon_next <- exons_sorted[i + 1, ]

    if (gene_strand == "+" || gene_strand == 1) {
      donor <- exon_current$end
      acceptor <- exon_next$start
    } else {
      donor <- exon_next$start
      acceptor <- exon_current$end
    }

    junctions[[length(junctions) + 1]] <- tibble(
      donor = donor,
      acceptor = acceptor,
      exon_upstream = exon_current$union_exon_number,
      exon_downstream = exon_next$union_exon_number
    )
  }

  return(bind_rows(junctions))
}

# ============================================================================
# Iterative Detection Functions
# ============================================================================

#' Generate unique signature for event deduplication
generate_event_signature <- function(event) {
  if (event$event_type == "A5SS") {
    paste0("A5SS:", min(event$position_A, event$position_B), ":",
           max(event$position_A, event$position_B))
  } else if (event$event_type == "A3SS") {
    paste0("A3SS:", min(event$position_A, event$position_B), ":",
           max(event$position_A, event$position_B))
  } else if (grepl("^IR", event$event_type)) {
    # Use positions and union exon numbers to create unique signature
    paste0(event$event_type, ":", event$position_A, ":", event$position_B, ":",
           event$union_exon_A, ":", event$union_exon_B)
  } else {
    paste0(event$event_type, ":", event$position_A, ":", event$position_B)
  }
}

#' Extract splice site coordinates for known splice sites registry
extract_splice_site_coords <- function(splice_events) {
  coords_list <- list()
  for (i in seq_along(splice_events)) {
    event <- splice_events[[i]]
    if (event$event_type == "A5SS") {
      coords_list[[length(coords_list) + 1]] <- tibble(
        type = "A5SS",
        donor_A = event$position_A,
        donor_B = event$position_B,
        acceptor_A = NA_integer_,
        acceptor_B = NA_integer_
      )
    } else if (event$event_type == "A3SS") {
      coords_list[[length(coords_list) + 1]] <- tibble(
        type = "A3SS",
        donor_A = NA_integer_,
        donor_B = NA_integer_,
        acceptor_A = event$position_A,
        acceptor_B = event$position_B
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

#' Detect A5SS and A3SS events
detect_splice_sites_with_flags <- function(junctions_A, junctions_B, exons_A, exons_B, iteration,
                                           gene_id, gene_strand, isoform_A_id, isoform_B_id) {
  events <- list()

  if (nrow(junctions_A) == 0 || nrow(junctions_B) == 0) {
    return(events)
  }

  for (i in seq_len(nrow(junctions_A))) {
    junc_A <- junctions_A[i, ]
    for (j in seq_len(nrow(junctions_B))) {
      junc_B <- junctions_B[j, ]

      if (junc_A$donor == junc_B$donor && junc_A$acceptor != junc_B$acceptor) {
        events[[length(events) + 1]] <- tibble(
          gene_id = gene_id,
          strand = gene_strand,
          isoform_A = isoform_A_id,
          isoform_B = isoform_B_id,
          event_type = "A3SS",
          union_exon_A = junc_A$exon_downstream,
          union_exon_B = junc_B$exon_downstream,
          position_A = junc_A$acceptor,
          position_B = junc_B$acceptor,
          detail = sprintf("A3SS: acceptor A=%d, B=%d (same donor=%d)",
                          junc_A$acceptor, junc_B$acceptor, junc_A$donor),
          iteration = iteration
        )
      } else if (junc_A$donor != junc_B$donor && junc_A$acceptor == junc_B$acceptor) {
        events[[length(events) + 1]] <- tibble(
          gene_id = gene_id,
          strand = gene_strand,
          isoform_A = isoform_A_id,
          isoform_B = isoform_B_id,
          event_type = "A5SS",
          union_exon_A = junc_A$exon_upstream,
          union_exon_B = junc_B$exon_upstream,
          position_A = junc_A$donor,
          position_B = junc_B$donor,
          detail = sprintf("A5SS: donor A=%d, B=%d (same acceptor=%d)",
                          junc_A$donor, junc_B$donor, junc_A$acceptor),
          iteration = iteration
        )
      }
    }
  }

  return(events)
}

#' Update flags based on detected splice events
update_flags_from_splice_events <- function(splice_events, exons_A, exons_B) {
  # splice_events is a list of tibbles, need to bind them first
  if (length(splice_events) > 0) {
    events_df <- bind_rows(splice_events)

    for (i in 1:nrow(events_df)) {
      event <- events_df[i, ]
      if (event$event_type == "A5SS") {
        exons_A <- exons_A %>% mutate(has_a5ss = if_else(end == event$position_A, TRUE, has_a5ss))
        exons_B <- exons_B %>% mutate(has_a5ss = if_else(end == event$position_B, TRUE, has_a5ss))
      } else if (event$event_type == "A3SS") {
        exons_A <- exons_A %>% mutate(has_a3ss = if_else(start == event$position_A, TRUE, has_a3ss))
        exons_B <- exons_B %>% mutate(has_a3ss = if_else(start == event$position_B, TRUE, has_a3ss))
      }
    }
  }
  return(list(exons_A = exons_A, exons_B = exons_B))
}

#' Detect internal IR with flag-informed boundary matching
detect_internal_ir_with_flags <- function(exons_A, exons_B, gene_strand, iteration, known_splice_sites,
                                          gene_id, isoform_A_id, isoform_B_id) {
  events <- list()

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

      # Check if consecutive (skip if NA)
      if (is.na(exon_B_last$transcript_exon_number) || is.na(exon_B_first$transcript_exon_number) ||
          exon_B_last$transcript_exon_number != exon_B_first$transcript_exon_number + 1) {
        next
      }

      encompasses <- (exon_A$start <= exon_B_first$start && exon_A$end >= exon_B_last$end)
      if (!encompasses) next

      # Precise boundary matching
      start_matches <- if (exon_A$start == exon_B_first$start) {
        TRUE
      } else if (exon_A$has_a3ss || exon_B_first$has_a3ss) {
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

      end_matches <- if (exon_A$end == exon_B_last$end) {
        TRUE
      } else if (exon_A$has_a5ss || exon_B_last$has_a5ss) {
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

        events[[length(events) + 1]] <- tibble(
          gene_id = gene_id,
          strand = gene_strand,
          isoform_A = isoform_A_id,
          isoform_B = isoform_B_id,
          event_type = event_type,
          union_exon_A = exon_A$union_exon_number,
          union_exon_B = NA_integer_,
          position_A = exon_A$start,
          position_B = exon_B_first$start,
          detail = sprintf("%s: A retains intron (%d-%d spans B exons %d-%d)",
                          event_type, exon_A$start, exon_A$end,
                          exon_B_first$start, exon_B_last$end),
          iteration = iteration
        )
      }
    }
  }

  # Symmetric check: B spanning A
  if (nrow(internal_B) > 0 && nrow(internal_A) > 1) {
    for (i in seq_len(nrow(internal_B))) {
      exon_B <- internal_B[i, ]
      for (j in 1:(nrow(internal_A) - 1)) {
        exon_A_first <- internal_A[j, ]
        exon_A_last <- internal_A[j + 1, ]

        # Check if consecutive (skip if NA)
        if (is.na(exon_A_last$transcript_exon_number) || is.na(exon_A_first$transcript_exon_number) ||
            exon_A_last$transcript_exon_number != exon_A_first$transcript_exon_number + 1) {
          next
        }

      encompasses <- (exon_B$start <= exon_A_first$start && exon_B$end >= exon_A_last$end)
      if (!encompasses) next

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
          gene_id = gene_id,
          strand = gene_strand,
          isoform_A = isoform_A_id,
          isoform_B = isoform_B_id,
          event_type = event_type,
          union_exon_A = NA_integer_,
          union_exon_B = exon_B$union_exon_number,
          position_A = exon_A_first$start,
          position_B = exon_B$start,
          detail = sprintf("%s: B retains intron (%d-%d spans A exons %d-%d)",
                          event_type, exon_B$start, exon_B$end,
                          exon_A_first$start, exon_A_last$end),
          iteration = iteration
        )
      }
    }
    }
  }

  return(bind_rows(events))
}

#' Detect monoexonic IR
detect_monoexonic_ir <- function(exons_A, exons_B, gene_strand, iteration,
                                 gene_id, isoform_A_id, isoform_B_id) {
  events <- list()

  # Deduplicate exons by coordinates to handle monoexonic transcripts
  # that appear in multiple union exon groups
  unique_exons_A <- exons_A %>% distinct(start, end, .keep_all = TRUE)
  unique_exons_B <- exons_B %>% distinct(start, end, .keep_all = TRUE)

  is_mono_A <- nrow(unique_exons_A) == 1
  is_mono_B <- nrow(unique_exons_B) == 1

  if (!is_mono_A && !is_mono_B) {
    return(bind_rows(events))
  }

  if (is_mono_A && !is_mono_B) {
    mono_exon <- unique_exons_A[1, ]
    min_B <- min(unique_exons_B$start)
    max_B <- max(unique_exons_B$end)

    if (mono_exon$start <= min_B && mono_exon$end >= max_B) {
      events[[length(events) + 1]] <- tibble(
        gene_id = gene_id,
        strand = gene_strand,
        isoform_A = isoform_A_id,
        isoform_B = isoform_B_id,
        event_type = "IR_monoexonic",
        union_exon_A = mono_exon$union_exon_number,
        union_exon_B = NA_integer_,
        position_A = mono_exon$start,
        position_B = min_B,
        detail = sprintf("IR_monoexonic: A is monoexonic (%d-%d), B has %d exons",
                        mono_exon$start, mono_exon$end, nrow(unique_exons_B)),
        iteration = iteration
      )
    }
  } else if (is_mono_B && !is_mono_A) {
    mono_exon <- unique_exons_B[1, ]
    min_A <- min(unique_exons_A$start)
    max_A <- max(unique_exons_A$end)

    if (mono_exon$start <= min_A && mono_exon$end >= max_A) {
      events[[length(events) + 1]] <- tibble(
        gene_id = gene_id,
        strand = gene_strand,
        isoform_A = isoform_A_id,
        isoform_B = isoform_B_id,
        event_type = "IR_monoexonic",
        union_exon_A = NA_integer_,
        union_exon_B = mono_exon$union_exon_number,
        position_A = min_A,
        position_B = mono_exon$start,
        detail = sprintf("IR_monoexonic: B is monoexonic (%d-%d), A has %d exons",
                        mono_exon$start, mono_exon$end, nrow(unique_exons_A)),
        iteration = iteration
      )
    }
  }

  return(bind_rows(events))
}

#' Main iterative event detection
detect_events_iteratively <- function(union_exons, isoform_A_id, isoform_B_id, gene_id, gene_strand) {
  # Initialize exons with flags
  exons_A <- get_isoform_exon_details(union_exons, isoform_A_id) %>%
    mutate(has_a5ss = FALSE, has_a3ss = FALSE, has_ir = FALSE)

  exons_B <- get_isoform_exon_details(union_exons, isoform_B_id) %>%
    mutate(has_a5ss = FALSE, has_a3ss = FALSE, has_ir = FALSE)

  all_events <- list()
  detected_signatures <- character()
  known_splice_sites <- tibble(type = character(), donor_A = integer(), donor_B = integer(),
                                acceptor_A = integer(), acceptor_B = integer())
  iteration <- 0

  # Iterative detection loop
  repeat {
    iteration <- iteration + 1
    iteration_events <- list()

    # STEP 1: A5SS/A3SS Detection
    junctions_A <- extract_junctions(exons_A, gene_strand)
    junctions_B <- extract_junctions(exons_B, gene_strand)
    splice_events <- detect_splice_sites_with_flags(junctions_A, junctions_B, exons_A, exons_B, iteration,
                                                     gene_id, gene_strand, isoform_A_id, isoform_B_id)

    if (length(splice_events) > 0) {
      splice_signatures <- sapply(splice_events, function(e) generate_event_signature(e))
      truly_new <- !splice_signatures %in% detected_signatures

      if (any(truly_new)) {
        new_splice_events <- splice_events[truly_new]
        iteration_events <- c(iteration_events, new_splice_events)
        detected_signatures <- c(detected_signatures, splice_signatures[truly_new])

        flag_updates <- update_flags_from_splice_events(new_splice_events, exons_A, exons_B)
        exons_A <- flag_updates$exons_A
        exons_B <- flag_updates$exons_B
        known_splice_sites <- bind_rows(known_splice_sites, extract_splice_site_coords(new_splice_events))
      }
    }

    # STEP 2: IR Detection
    ir_events_internal <- detect_internal_ir_with_flags(exons_A, exons_B, gene_strand, iteration, known_splice_sites,
                                                         gene_id, isoform_A_id, isoform_B_id)
    ir_events_monoexonic <- detect_monoexonic_ir(exons_A, exons_B, gene_strand, iteration,
                                                  gene_id, isoform_A_id, isoform_B_id)
    ir_events <- bind_rows(ir_events_internal, ir_events_monoexonic)

    if (nrow(ir_events) > 0) {
      ir_signatures <- sapply(1:nrow(ir_events), function(i) {
        generate_event_signature(ir_events[i, ])
      })
      truly_new <- !ir_signatures %in% detected_signatures

      if (any(truly_new)) {
        new_ir_events <- ir_events[truly_new, ]
        iteration_events <- c(iteration_events, list(new_ir_events))
        detected_signatures <- c(detected_signatures, ir_signatures[truly_new])
      }
    }

    # STEP 3: Check Convergence
    if (length(iteration_events) == 0) {
      break
    }

    all_events <- c(all_events, iteration_events)

    if (iteration > 100) {
      warning(paste("Exceeded 100 iterations for", gene_id, isoform_A_id, "vs", isoform_B_id))
      break
    }
  }

  # STEP 4: Detect other events (Alt_TSS, Alt_TES, SE)
  terminal_events <- detect_terminal_events(exons_A, exons_B, gene_strand, gene_id, isoform_A_id, isoform_B_id)
  se_events <- detect_exon_skipping(exons_A, exons_B, gene_id, isoform_A_id, isoform_B_id)

  # STEP 5: Combine all events
  final_events <- bind_rows(all_events, terminal_events, se_events)

  if (nrow(final_events) > 0) {
    final_events <- final_events %>% mutate(total_iterations = iteration)
  }

  return(final_events)
}

#' Detect terminal events (Alt_TSS/Alt_TES)
detect_terminal_events <- function(exons_A, exons_B, gene_strand, gene_id, isoform_A_id, isoform_B_id) {
  events <- list()

  first_A <- exons_A %>% filter(is_first) %>% slice(1)
  first_B <- exons_B %>% filter(is_first) %>% slice(1)

  if (nrow(first_A) > 0 && nrow(first_B) > 0) {
    if (gene_strand == "+" || gene_strand == 1) {
      pos_A <- first_A$start
      pos_B <- first_B$start
    } else {
      pos_A <- first_A$end
      pos_B <- first_B$end
    }

    if (abs(pos_A - pos_B) > TSS_TES_THRESHOLD) {
      events[[length(events) + 1]] <- tibble(
        gene_id = gene_id,
        strand = gene_strand,
        isoform_A = isoform_A_id,
        isoform_B = isoform_B_id,
        event_type = "Alt_TSS",
        union_exon_A = first_A$union_exon_number,
        union_exon_B = first_B$union_exon_number,
        position_A = pos_A,
        position_B = pos_B,
        detail = sprintf("Different TSS: A=%d, B=%d (diff=%d bp)",
                        pos_A, pos_B, abs(pos_A - pos_B))
      )
    }
  }

  last_A <- exons_A %>% filter(is_last) %>% slice(1)
  last_B <- exons_B %>% filter(is_last) %>% slice(1)

  if (nrow(last_A) > 0 && nrow(last_B) > 0) {
    if (gene_strand == "+" || gene_strand == 1) {
      pos_A <- last_A$end
      pos_B <- last_B$end
    } else {
      pos_A <- last_A$start
      pos_B <- last_B$start
    }

    if (abs(pos_A - pos_B) > TSS_TES_THRESHOLD) {
      events[[length(events) + 1]] <- tibble(
        gene_id = gene_id,
        strand = gene_strand,
        isoform_A = isoform_A_id,
        isoform_B = isoform_B_id,
        event_type = "Alt_TES",
        union_exon_A = last_A$union_exon_number,
        union_exon_B = last_B$union_exon_number,
        position_A = pos_A,
        position_B = pos_B,
        detail = sprintf("Different TES: A=%d, B=%d (diff=%d bp)",
                        pos_A, pos_B, abs(pos_A - pos_B))
      )
    }
  }

  return(bind_rows(events))
}

#' Detect exon skipping (SE)
detect_exon_skipping <- function(exons_A, exons_B, gene_id, isoform_A_id, isoform_B_id) {
  events <- list()

  exons_A_numbers <- exons_A$union_exon_number
  exons_B_numbers <- exons_B$union_exon_number

  only_A <- setdiff(exons_A_numbers, exons_B_numbers)
  only_B <- setdiff(exons_B_numbers, exons_A_numbers)

  internal_only_A <- exons_A %>% filter(union_exon_number %in% only_A, !is_first, !is_last)
  internal_only_B <- exons_B %>% filter(union_exon_number %in% only_B, !is_first, !is_last)

  if (nrow(internal_only_A) > 0) {
    for (i in 1:nrow(internal_only_A)) {
      events[[length(events) + 1]] <- tibble(
        gene_id = gene_id,
        strand = NA_character_,
        isoform_A = isoform_A_id,
        isoform_B = isoform_B_id,
        event_type = "SE",
        union_exon_A = internal_only_A$union_exon_number[i],
        union_exon_B = NA_integer_,
        position_A = internal_only_A$start[i],
        position_B = NA_integer_,
        detail = sprintf("Exon %d skipped in B", internal_only_A$union_exon_number[i])
      )
    }
  }

  if (nrow(internal_only_B) > 0) {
    for (i in 1:nrow(internal_only_B)) {
      events[[length(events) + 1]] <- tibble(
        gene_id = gene_id,
        strand = NA_character_,
        isoform_A = isoform_A_id,
        isoform_B = isoform_B_id,
        event_type = "SE",
        union_exon_A = NA_integer_,
        union_exon_B = internal_only_B$union_exon_number[i],
        position_A = NA_integer_,
        position_B = internal_only_B$start[i],
        detail = sprintf("Exon %d skipped in A", internal_only_B$union_exon_number[i])
      )
    }
  }

  return(bind_rows(events))
}

# ============================================================================
# Load Data
# ============================================================================

cat("═══ Loading Data ═══\n\n")

cat("Loading union exon models...\n")
union_models <- readRDS(file.path(output_dir, "union_exon_models_major.rds"))
cat("  Genes:", length(union_models), "\n\n")

cat("Loading exon structures for strand info...\n")
exon_structures <- readRDS(file.path(output_dir, "exon_structures_expanded.rds"))
gene_strand_map <- exon_structures %>%
  distinct(gene_id, strand) %>%
  mutate(strand = as.character(strand))
cat("  Strand info for", nrow(gene_strand_map), "genes\n\n")

# ============================================================================
# Generate Isoform Pairs
# ============================================================================

cat("═══ Generating Isoform Pairs ═══\n\n")

all_pairs <- list()

for (gene_id in names(union_models)) {
  model <- union_models[[gene_id]]
  isoform_ids <- unique(unlist(map(model$union_exons, ~.x$isoform_id)))
  n_isoforms <- length(isoform_ids)

  if (n_isoforms >= 2) {
    pairs <- combn(isoform_ids, 2, simplify = FALSE)
    for (pair in pairs) {
      all_pairs[[length(all_pairs) + 1]] <- tibble(
        gene_id = gene_id,
        isoform_A = pair[1],
        isoform_B = pair[2],
        n_isoforms_in_gene = n_isoforms,
        n_union_exons = model$n_union_exons
      )
    }
  }
}

isoform_pairs <- bind_rows(all_pairs)
cat("Total pairs:", nrow(isoform_pairs), "\n")
cat("Genes:", length(unique(isoform_pairs$gene_id)), "\n\n")

# ============================================================================
# Detect Events with Iteration
# ============================================================================

cat("═══ Detecting Events with IR ═══\n\n")

pairs_by_gene <- isoform_pairs %>%
  group_by(gene_id) %>%
  nest() %>%
  ungroup()

all_events <- list()
genes_processed <- 0
start_time <- Sys.time()

for (i in 1:nrow(pairs_by_gene)) {
  gene_id <- pairs_by_gene$gene_id[i]
  gene_pairs <- pairs_by_gene$data[[i]]
  model <- union_models[[gene_id]]

  gene_strand <- gene_strand_map %>%
    filter(gene_id == !!gene_id) %>%
    pull(strand) %>%
    first()

  if (is.na(gene_strand)) {
    cat(sprintf("WARNING: No strand for %s\n", gene_id))
    next
  }

  for (j in 1:nrow(gene_pairs)) {
    pair <- gene_pairs[j, ]

    events <- detect_events_iteratively(
      union_exons = model$union_exons,
      isoform_A_id = pair$isoform_A,
      isoform_B_id = pair$isoform_B,
      gene_id = gene_id,
      gene_strand = gene_strand
    )

    if (nrow(events) > 0) {
      events <- events %>%
        mutate(
          n_isoforms_in_gene = pair$n_isoforms_in_gene,
          n_union_exons = pair$n_union_exons,
          pair_id = paste(gene_id, pair$isoform_A, pair$isoform_B, sep = "|")
        )
      all_events[[length(all_events) + 1]] <- events
    } else {
      all_events[[length(all_events) + 1]] <- tibble(
        gene_id = gene_id,
        strand = gene_strand,
        isoform_A = pair$isoform_A,
        isoform_B = pair$isoform_B,
        event_type = "NO_DIFF",
        union_exon_A = NA_integer_,
        union_exon_B = NA_integer_,
        position_A = NA_integer_,
        position_B = NA_integer_,
        detail = "No events detected (differences below thresholds or complex structural variation)",
        total_iterations = 0L,
        n_isoforms_in_gene = pair$n_isoforms_in_gene,
        n_union_exons = pair$n_union_exons,
        pair_id = paste(gene_id, pair$isoform_A, pair$isoform_B, sep = "|")
      )
    }
  }

  genes_processed <- genes_processed + 1

  if (genes_processed %% 100 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    eta <- (elapsed / genes_processed) * (nrow(pairs_by_gene) - genes_processed)
    cat(sprintf("[%d/%d] %s | Elapsed: %.1f min | ETA: %.1f min\n",
                genes_processed, nrow(pairs_by_gene), gene_id, elapsed, eta))
  }

  if (!is.na(CHECKPOINT_BATCH) && genes_processed == CHECKPOINT_BATCH * BATCH_SIZE) {
    cat("\n═══ CHECKPOINT ═══\n\n")
    events_checkpoint <- bind_rows(all_events)
    saveRDS(events_checkpoint, file.path(output_dir, "events_with_ir_checkpoint.rds"))
    cat("Checkpoint saved. Events:", nrow(events_checkpoint), "\n")
    print(events_checkpoint %>% count(event_type) %>% arrange(desc(n)))
    cat("\nSTOPPING AT CHECKPOINT\n\n")
    break
  }
}

cat("\nEvent detection complete\n\n")

# ============================================================================
# Save Results
# ============================================================================

cat("═══ Saving Results ═══\n\n")

all_events_df <- bind_rows(all_events)

output_file <- file.path(output_dir, "isoform_pairs_events_with_ir.rds")
saveRDS(all_events_df, output_file)
cat("Saved:", output_file, "\n")
cat("Total rows:", nrow(all_events_df), "\n\n")

event_counts <- all_events_df %>%
  count(event_type) %>%
  arrange(desc(n))

print(event_counts)

summary_file <- file.path(output_dir, "event_summary_with_ir.tsv")
write_tsv(event_counts, summary_file)
cat("\nSaved summary:", summary_file, "\n\n")

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EVENT DETECTION WITH IR COMPLETE                            ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")
