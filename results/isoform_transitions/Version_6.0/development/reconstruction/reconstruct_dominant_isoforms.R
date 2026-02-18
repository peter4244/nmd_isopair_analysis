#!/usr/bin/env Rscript
#
# Reconstruct Dominant Isoforms
#
# This script reconstructs dominant isoforms from comparator isoforms
# by applying detected events.
#
# Inputs:
#   - Comparator GTF (comparator isoforms only)
#   - Events file (with all event information)
#   - Union exons file (for IR reconstruction)
#
# Outputs:
#   - Reconstructed dominant GTF
#   - Reconstruction log

library(tidyverse)

# Source reconstruction functions
# Try to find the file in common locations
if (file.exists("reconstruction_functions.R")) {
  source("reconstruction_functions.R")
} else if (file.exists("reconstruction/reconstruction_functions.R")) {
  source("reconstruction/reconstruction_functions.R")
} else {
  stop("Cannot find reconstruction_functions.R")
}

# ==============================================================================
# Parse GTF
# ==============================================================================

#' Parse GTF file to extract exon information
#'
#' @param gtf_file Path to GTF file
#' @return Tibble with exon information
parse_gtf <- function(gtf_file) {
  cat("Reading GTF file:", gtf_file, "\n")

  lines <- readLines(gtf_file)
  lines <- lines[!startsWith(lines, '#') & lines != '']

  results <- list()
  for (line in lines) {
    parts <- strsplit(line, '\t')[[1]]

    # Only parse exon features
    if (parts[3] != "exon") {
      next
    }

    attrs_str <- parts[9]
    gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', attrs_str)
    transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', attrs_str)

    # Extract exon_number if present
    exon_number_str <- sub('.*exon_number "([^"]+)".*', '\\1', attrs_str)
    exon_number <- if (grepl('exon_number', attrs_str)) as.integer(exon_number_str) else NA_integer_

    results[[length(results) + 1]] <- tibble(
      chr = parts[1],
      source = parts[2],
      feature = parts[3],
      exon_start = as.integer(parts[4]),
      exon_end = as.integer(parts[5]),
      score = parts[6],
      strand = parts[7],
      frame = parts[8],
      attributes = attrs_str,
      gene_id = gene_id,
      transcript_id = transcript_id,
      exon_number = exon_number
    )
  }

  gtf_data <- bind_rows(results)

  cat(sprintf("  Parsed %d exon records\n", nrow(gtf_data)))
  cat(sprintf("  %d genes, %d transcripts\n",
              n_distinct(gtf_data$gene_id),
              n_distinct(gtf_data$transcript_id)))

  return(gtf_data)
}

# ==============================================================================
# Load Union Exons
# ==============================================================================

#' Load union exons from tabix-indexed file
#'
#' @param union_exon_file Path to union exon file (.tsv.gz)
#' @return Tibble with union exons
load_union_exons <- function(union_exon_file) {
  cat("Reading union exons:", union_exon_file, "\n")

  # Read gzipped file
  union_exons <- read_tsv(union_exon_file,
                          comment = "#",
                          col_names = c("chr", "start", "end", "union_exon_id", "strand", "gene_id"),
                          col_types = "ciiccc",
                          show_col_types = FALSE)

  cat(sprintf("  Loaded %d union exons\n", nrow(union_exons)))

  return(union_exons)
}

# ==============================================================================
# Write GTF
# ==============================================================================

#' Write reconstructed GTF file
#'
#' @param gtf_data Tibble with GTF records
#' @param output_file Output GTF file path
write_gtf <- function(gtf_data, output_file) {
  cat("Writing GTF file:", output_file, "\n")

  # Reconstruct GTF lines
  gtf_lines <- with(gtf_data, {
    paste(chr, source, feature, exon_start, exon_end, score, strand, frame, attributes, sep = "\t")
  })

  writeLines(gtf_lines, output_file)

  cat(sprintf("  Wrote %d exon records\n", nrow(gtf_data)))
  cat(sprintf("  %d genes, %d transcripts\n",
              n_distinct(gtf_data$gene_id),
              n_distinct(gtf_data$transcript_id)))
}

# ==============================================================================
# Reconstruct Attributes String
# ==============================================================================

#' Reconstruct GTF attributes string for reconstructed transcript
#'
#' @param gene_id Gene ID
#' @param transcript_id Original transcript ID (will be modified)
#' @param exon_number Exon number
#' @return Attributes string
reconstruct_attributes <- function(gene_id, transcript_id, exon_number) {
  # Mark as reconstructed by appending "_reconstructed"
  recon_tid <- paste0(transcript_id, "_reconstructed")

  sprintf('gene_id "%s"; transcript_id "%s"; exon_number "%d";',
          gene_id, recon_tid, exon_number)
}

# ==============================================================================
# Main Reconstruction
# ==============================================================================

#' Reconstruct all dominant isoforms
#'
#' @param comparator_gtf Comparator GTF data
#' @param events Events data
#' @param union_exons Union exon data
#' @return Reconstructed GTF data with logs
reconstruct_all <- function(comparator_gtf, events, union_exons) {

  # Get unique transcript pairs from events
  transcript_pairs <- events %>%
    select(gene_id, dominant_transcript_id, comparator_transcript_id) %>%
    distinct()

  cat(sprintf("\nReconstrucing %d transcript pairs...\n", nrow(transcript_pairs)))

  all_reconstructed <- list()
  reconstruction_log <- list()

  for (i in seq_len(nrow(transcript_pairs))) {
    pair <- transcript_pairs[i, ]

    cat(sprintf("\n[%d/%d] %s: %s -> %s (reconstructed)\n",
                i, nrow(transcript_pairs),
                pair$gene_id,
                pair$comparator_transcript_id,
                pair$dominant_transcript_id))

    # Extract comparator exons
    comp_exons <- comparator_gtf %>%
      filter(
        gene_id == pair$gene_id,
        transcript_id == pair$comparator_transcript_id
      ) %>%
      arrange(exon_start)

    if (nrow(comp_exons) == 0) {
      cat("  Warning: No comparator exons found\n")
      reconstruction_log[[i]] <- tibble(
        gene_id = pair$gene_id,
        comparator_transcript_id = pair$comparator_transcript_id,
        dominant_transcript_id = pair$dominant_transcript_id,
        status = "FAILED",
        reason = "No comparator exons",
        n_events = 0,
        n_exons_input = 0,
        n_exons_output = 0
      )
      next
    }

    # Filter events to this transcript pair
    pair_events <- events %>%
      filter(
        gene_id == pair$gene_id,
        comparator_transcript_id == pair$comparator_transcript_id,
        dominant_transcript_id == pair$dominant_transcript_id
      )

    cat(sprintf("  Input: %d exons\n", nrow(comp_exons)))
    cat(sprintf("  Events: %d (%s)\n",
                nrow(pair_events),
                paste(pair_events$event_type, collapse = ", ")))

    # Reconstruct dominant
    reconstructed_exons <- tryCatch({
      reconstruct_dominant_v2(comp_exons, pair_events, union_exons)
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      reconstruction_log[[i]] <- tibble(
        gene_id = pair$gene_id,
        comparator_transcript_id = pair$comparator_transcript_id,
        dominant_transcript_id = pair$dominant_transcript_id,
        status = "FAILED",
        reason = e$message,
        n_events = nrow(pair_events),
        n_exons_input = nrow(comp_exons),
        n_exons_output = 0
      )
      return(NULL)
    })

    if (is.null(reconstructed_exons)) {
      next
    }

    cat(sprintf("  Output: %d exons\n", nrow(reconstructed_exons)))

    # Reconstruct GTF attributes
    reconstructed_gtf <- reconstructed_exons %>%
      mutate(
        source = comp_exons$source[1],
        feature = "exon",
        score = ".",
        frame = ".",
        exon_number = row_number(),
        attributes = reconstruct_attributes(gene_id, pair$dominant_transcript_id, exon_number)
      ) %>%
      select(chr, source, feature, exon_start, exon_end, score, strand, frame, attributes,
             gene_id, transcript_id, exon_number)

    all_reconstructed[[i]] <- reconstructed_gtf

    # Log success
    reconstruction_log[[i]] <- tibble(
      gene_id = pair$gene_id,
      comparator_transcript_id = pair$comparator_transcript_id,
      dominant_transcript_id = pair$dominant_transcript_id,
      status = "SUCCESS",
      reason = "",
      n_events = nrow(pair_events),
      n_exons_input = nrow(comp_exons),
      n_exons_output = nrow(reconstructed_exons)
    )
  }

  # Combine all reconstructed transcripts
  final_gtf <- bind_rows(all_reconstructed)
  final_log <- bind_rows(reconstruction_log)

  return(list(gtf = final_gtf, log = final_log))
}

# ==============================================================================
# Main Pipeline
# ==============================================================================

main <- function(comparator_gtf_file, events_file, union_exons_file,
                output_gtf, output_log) {
  cat("\n╔════════════════════════════════════════════════════════════════╗\n")
  cat("║   Reconstruct Dominant Isoforms Pipeline                      ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")

  # Load inputs
  comparator_gtf <- parse_gtf(comparator_gtf_file)

  cat("\nReading events file:", events_file, "\n")
  events <- read_tsv(events_file, show_col_types = FALSE)
  cat(sprintf("  %d events loaded\n", nrow(events)))

  union_exons <- load_union_exons(union_exons_file)

  # Reconstruct
  results <- reconstruct_all(comparator_gtf, events, union_exons)

  # Write outputs
  cat("\n")
  write_gtf(results$gtf, output_gtf)

  cat("\nWriting reconstruction log:", output_log, "\n")
  write_tsv(results$log, output_log)

  # Summary statistics
  cat("\n═══════════════════════════════════════════════════════════════════\n")
  cat("RECONSTRUCTION COMPLETE\n")
  cat("═══════════════════════════════════════════════════════════════════\n")

  success_count <- sum(results$log$status == "SUCCESS")
  fail_count <- sum(results$log$status == "FAILED")

  cat(sprintf("Success: %d/%d (%.1f%%)\n",
              success_count,
              nrow(results$log),
              100 * success_count / nrow(results$log)))
  cat(sprintf("Failed:  %d/%d\n", fail_count, nrow(results$log)))

  if (fail_count > 0) {
    cat("\nFailed transcripts:\n")
    failed <- results$log %>% filter(status == "FAILED")
    for (i in seq_len(nrow(failed))) {
      cat(sprintf("  - %s: %s\n",
                  failed$comparator_transcript_id[i],
                  failed$reason[i]))
    }
  }

  cat(sprintf("\nOutputs:\n"))
  cat(sprintf("  GTF: %s\n", output_gtf))
  cat(sprintf("  Log: %s\n", output_log))
}

# ==============================================================================
# Command Line Interface
# ==============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 5) {
    cat("Usage: reconstruct_dominant_isoforms.R <comparator_gtf> <events_file> <union_exons_file> <output_gtf> <output_log>\n")
    quit(status = 1)
  }

  main(args[1], args[2], args[3], args[4], args[5])
} else {
  # Interactive testing
  main(
    comparator_gtf_file = "../synthetic/base_events_comparator.gtf",
    events_file = "../synthetic/base_events_events.tsv",
    union_exons_file = "../synthetic/base_events_union_exons.tsv.gz",
    output_gtf = "../synthetic/base_events_reconstructed.gtf",
    output_log = "../synthetic/reconstruction_log.tsv"
  )
}
