#!/usr/bin/env Rscript
#
# Extract Comparator GTF
#
# This script creates a GTF file containing only comparator isoforms
# based on the events file.
#
# Inputs:
#   - Original GTF (all isoforms)
#   - Events file (with comparator_transcript_id column)
#
# Outputs:
#   - Comparator-only GTF

library(tidyverse)

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
      start = as.integer(parts[4]),
      end = as.integer(parts[5]),
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
# Write GTF
# ==============================================================================

#' Write GTF file
#'
#' @param gtf_data Tibble with GTF records
#' @param output_file Output GTF file path
write_gtf <- function(gtf_data, output_file) {
  cat("Writing GTF file:", output_file, "\n")

  # Reconstruct GTF lines
  gtf_lines <- with(gtf_data, {
    paste(chr, source, feature, start, end, score, strand, frame, attributes, sep = "\t")
  })

  writeLines(gtf_lines, output_file)

  cat(sprintf("  Wrote %d exon records\n", nrow(gtf_data)))
  cat(sprintf("  %d genes, %d transcripts\n",
              n_distinct(gtf_data$gene_id),
              n_distinct(gtf_data$transcript_id)))
}

# ==============================================================================
# Main Pipeline
# ==============================================================================

main <- function(original_gtf, events_file, output_gtf) {
  cat("\n╔════════════════════════════════════════════════════════════════╗\n")
  cat("║   Extract Comparator GTF                                      ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")

  # Parse original GTF
  gtf_data <- parse_gtf(original_gtf)

  # Load events file
  cat("\nReading events file:", events_file, "\n")
  events <- read_tsv(events_file, show_col_types = FALSE)
  cat(sprintf("  %d events loaded\n", nrow(events)))

  # Extract unique comparator transcript IDs
  comparator_ids <- unique(events$comparator_transcript_id)
  cat(sprintf("\n%d unique comparator transcripts identified\n", length(comparator_ids)))

  # Filter GTF to comparator transcripts only
  comparator_gtf <- gtf_data %>%
    filter(transcript_id %in% comparator_ids) %>%
    arrange(chr, start, transcript_id)

  cat(sprintf("\nFiltered to %d exons from %d comparator transcripts\n",
              nrow(comparator_gtf),
              n_distinct(comparator_gtf$transcript_id)))

  # Write comparator GTF
  write_gtf(comparator_gtf, output_gtf)

  cat("\n═══════════════════════════════════════════════════════════════════\n")
  cat("COMPARATOR GTF EXTRACTION COMPLETE\n")
  cat("═══════════════════════════════════════════════════════════════════\n")
  cat(sprintf("Output: %s\n", output_gtf))
}

# ==============================================================================
# Command Line Interface
# ==============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 3) {
    cat("Usage: extract_comparator_gtf.R <original_gtf> <events_file> <output_gtf>\n")
    quit(status = 1)
  }

  main(args[1], args[2], args[3])
} else {
  # Interactive testing
  main(
    original_gtf = "../synthetic/TestData/exons/base_events.gtf",
    events_file = "../synthetic/base_events_events.tsv",
    output_gtf = "../synthetic/base_events_comparator.gtf"
  )
}
