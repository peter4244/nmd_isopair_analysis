#!/usr/bin/env Rscript
#
# Create Union Exons and Junctions from GTF
#
# Inputs:
#   - GTF file with all isoforms
#   - Dominant isoform specification (gene_id, dominant_transcript_id)
#
# Outputs:
#   - Union exon file (tabix-indexed)
#   - Junction file (tabix-indexed)
#
# All coordinates are 1-based (GTF convention)

library(tidyverse)

# ==============================================================================
# Parse GTF
# ==============================================================================

#' Parse GTF file to extract exon records
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
      start = as.integer(parts[4]),  # 1-based
      end = as.integer(parts[5]),    # 1-based
      strand = parts[7],
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
# Create Union Exons
# ==============================================================================

#' Create union exons by merging overlapping exons
#'
#' @param exons Tibble with exon coordinates (chr, start, end, strand)
#' @return Tibble with merged union exons
create_union_exons <- function(exons) {
  if (nrow(exons) == 0) {
    return(tibble(chr = character(), start = integer(), end = integer(),
                  strand = character()))
  }

  # Sort by position
  exons <- exons %>%
    arrange(start, end)

  # Merge overlapping exons
  union_exons <- list()
  current_start <- exons$start[1]
  current_end <- exons$end[1]

  for (i in 2:nrow(exons)) {
    if (exons$start[i] <= current_end) {
      # Overlapping, extend current union exon
      current_end <- max(current_end, exons$end[i])
    } else {
      # Non-overlapping, save current and start new
      union_exons[[length(union_exons) + 1]] <- tibble(
        chr = exons$chr[1],
        start = current_start,
        end = current_end,
        strand = exons$strand[1]
      )
      current_start <- exons$start[i]
      current_end <- exons$end[i]
    }
  }

  # Add final union exon
  union_exons[[length(union_exons) + 1]] <- tibble(
    chr = exons$chr[1],
    start = current_start,
    end = current_end,
    strand = exons$strand[1]
  )

  return(bind_rows(union_exons))
}

#' Create union exons for all genes
#'
#' @param gtf_data Parsed GTF data
#' @return Tibble with all union exons and assigned IDs
create_all_union_exons <- function(gtf_data) {
  cat("\nCreating union exons...\n")

  # Group by gene and create union exons
  all_union_exons <- gtf_data %>%
    group_by(gene_id, chr, strand) %>%
    summarise(
      exons = list(tibble(start = start, end = end)),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      union_exons = list(create_union_exons(
        tibble(chr = chr, start = exons$start, end = exons$end, strand = strand)
      ))
    ) %>%
    select(gene_id, union_exons) %>%
    unnest(union_exons)

  # Assign sequential union_exon_ids
  all_union_exons <- all_union_exons %>%
    mutate(union_exon_id = sprintf("ue_%d", row_number())) %>%
    select(chr, start, end, union_exon_id, strand)

  cat(sprintf("  Created %d union exons\n", nrow(all_union_exons)))

  return(all_union_exons)
}

# ==============================================================================
# Create Junctions
# ==============================================================================

#' Extract junctions from a single transcript
#'
#' @param exons Tibble with exon coordinates for one transcript
#' @param strand Strand (+ or -)
#' @return Tibble with junctions
extract_junctions <- function(exons, strand) {
  if (nrow(exons) <= 1) {
    return(tibble(chr = character(), five_prime = integer(),
                  three_prime = integer(), strand = character()))
  }

  # Sort exons by genomic position
  exons <- exons %>% arrange(start)

  junctions <- list()
  for (i in 1:(nrow(exons) - 1)) {
    if (strand == "+") {
      # Plus strand: 5' = donor = end of upstream exon
      #              3' = acceptor = start of downstream exon
      five_prime <- exons$end[i]
      three_prime <- exons$start[i + 1]
    } else {
      # Minus strand: 5' = donor = start of downstream exon
      #               3' = acceptor = end of upstream exon
      five_prime <- exons$start[i + 1]
      three_prime <- exons$end[i]
    }

    junctions[[length(junctions) + 1]] <- tibble(
      chr = exons$chr[i],
      five_prime = five_prime,
      three_prime = three_prime,
      strand = strand
    )
  }

  return(bind_rows(junctions))
}

#' Create junctions for all transcripts
#'
#' @param gtf_data Parsed GTF data
#' @return Tibble with all unique junctions and assigned IDs
create_all_junctions <- function(gtf_data) {
  cat("\nCreating junctions...\n")

  # Extract junctions for each transcript
  all_junctions <- gtf_data %>%
    group_by(transcript_id, chr, strand) %>%
    arrange(start) %>%
    summarise(
      exons = list(tibble(chr = chr, start = start, end = end)),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      junctions = list(extract_junctions(exons, strand))
    ) %>%
    select(junctions) %>%
    unnest(junctions)

  # Get unique junctions
  unique_junctions <- all_junctions %>%
    distinct(chr, five_prime, three_prime, strand)

  # Assign sequential junction_ids and add start/end for tabix
  unique_junctions <- unique_junctions %>%
    mutate(
      junction_id = sprintf("j_%d", row_number()),
      start = pmin(five_prime, three_prime),  # min for tabix
      end = pmax(five_prime, three_prime),    # max for tabix
      type = ""  # Empty for now, can add backsplice detection later
    ) %>%
    select(junction_id, chr, start, end, five_prime, three_prime, strand, type)

  cat(sprintf("  Created %d unique junctions\n", nrow(unique_junctions)))

  return(unique_junctions)
}

# ==============================================================================
# Write and Index Files
# ==============================================================================

#' Write union exon file and index with tabix
#'
#' @param union_exons Tibble with union exons
#' @param output_file Output file path
write_union_exons <- function(union_exons, output_file) {
  cat("\nWriting union exon file:", output_file, "\n")

  # Sort by chr, start
  union_exons_sorted <- union_exons %>%
    arrange(chr, start)

  # Write with header
  write_lines("#chr\tstart\tend\tunion_exon_id\tstrand", output_file)
  write_tsv(union_exons_sorted, output_file, append = TRUE, col_names = FALSE)

  # Compress and index
  system(sprintf("bgzip -f %s", output_file))
  system(sprintf("tabix -f -s 1 -b 2 -e 3 %s.gz", output_file))

  cat(sprintf("  Wrote %d union exons\n", nrow(union_exons_sorted)))
  cat(sprintf("  Indexed: %s.gz\n", output_file))
}

#' Write junction file and index with tabix
#'
#' @param junctions Tibble with junctions
#' @param output_file Output file path
write_junctions <- function(junctions, output_file) {
  cat("\nWriting junction file:", output_file, "\n")

  # Sort by chr, start (for tabix)
  junctions_sorted <- junctions %>%
    arrange(chr, start)

  # Write with header
  write_lines("#junction_id\tchr\tstart\tend\t5_prime\t3_prime\tstrand\ttype", output_file)
  write_tsv(junctions_sorted, output_file, append = TRUE, col_names = FALSE)

  # Compress and index (using start/end columns for tabix)
  system(sprintf("bgzip -f %s", output_file))
  system(sprintf("tabix -f -s 2 -b 3 -e 4 %s.gz", output_file))

  cat(sprintf("  Wrote %d junctions\n", nrow(junctions_sorted)))
  cat(sprintf("  Indexed: %s.gz\n", output_file))
}

# ==============================================================================
# Main Pipeline
# ==============================================================================

main <- function(gtf_file, output_prefix) {
  cat("\n╔════════════════════════════════════════════════════════════════╗\n")
  cat("║   Union Exon and Junction Creation Pipeline                  ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")

  # Parse GTF
  gtf_data <- parse_gtf(gtf_file)

  # Create union exons
  union_exons <- create_all_union_exons(gtf_data)

  # Create junctions
  junctions <- create_all_junctions(gtf_data)

  # Write output files
  union_exon_file <- paste0(output_prefix, "_union_exons.tsv")
  junction_file <- paste0(output_prefix, "_junctions.tsv")

  write_union_exons(union_exons, union_exon_file)
  write_junctions(junctions, junction_file)

  cat("\n═══════════════════════════════════════════════════════════════════\n")
  cat("PIPELINE COMPLETE\n")
  cat("═══════════════════════════════════════════════════════════════════\n")
  cat(sprintf("Union exons: %s.gz\n", union_exon_file))
  cat(sprintf("Junctions:   %s.gz\n", junction_file))
}

# ==============================================================================
# Command Line Interface
# ==============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 2) {
    cat("Usage: create_union_exons_and_junctions.R <gtf_file> <output_prefix>\n")
    cat("\n")
    cat("Arguments:\n")
    cat("  gtf_file       Input GTF file\n")
    cat("  output_prefix  Output file prefix (will create <prefix>_union_exons.tsv.gz and <prefix>_junctions.tsv.gz)\n")
    quit(status = 1)
  }

  main(args[1], args[2])
}
