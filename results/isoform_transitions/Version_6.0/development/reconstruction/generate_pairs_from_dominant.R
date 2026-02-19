#!/usr/bin/env Rscript
#
# Generate Isoform Pairs from Dominant Mapping
#
# Takes a simple dominant isoform mapping (gene_id, dominant_isoform_id) and
# a GTF file, then generates all pairs comparing dominant to non-dominant
# isoforms within each gene.
#
# Inputs:
#   - GTF file (all isoforms)
#   - Dominant mapping file (gene_id, dominant_isoform_id)
#
# Output:
#   - Pairs file (gene_id, isoform_A, isoform_B) where:
#     * isoform_A = dominant (from mapping)
#     * isoform_B = each non-dominant isoform in the gene
#

library(tidyverse)

#' Parse GTF file to get transcript information
#'
#' @param gtf_file Path to GTF file
#' @return Tibble with gene and transcript information
parse_gtf_transcripts <- function(gtf_file) {
  cat("Reading GTF file:", gtf_file, "\n")

  lines <- readLines(gtf_file)
  lines <- lines[!startsWith(lines, '#') & lines != '']

  # Extract exon records (GTF may not have explicit transcript records)
  exon_lines <- grep('\texon\t', lines, value = TRUE)

  results <- list()
  for (line in exon_lines) {
    parts <- strsplit(line, '\t')[[1]]

    attrs_str <- parts[9]
    gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', attrs_str)
    transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', attrs_str)

    results[[length(results) + 1]] <- tibble(
      chr = parts[1],
      start = as.integer(parts[4]),
      end = as.integer(parts[5]),
      strand = parts[7],
      gene_id = gene_id,
      transcript_id = transcript_id
    )
  }

  exons <- bind_rows(results)

  # Get transcript-level info by grouping exons
  transcripts <- exons %>%
    group_by(gene_id, transcript_id, chr, strand) %>%
    summarize(
      start = min(start),
      end = max(end),
      .groups = 'drop'
    )

  cat(sprintf("  Parsed %d transcripts from %d genes\n",
              nrow(transcripts), n_distinct(transcripts$gene_id)))

  return(transcripts)
}

#' Generate pairs from dominant mapping
#'
#' @param gtf_transcripts Tibble from parse_gtf_transcripts
#' @param dominant_mapping Tibble with gene_id and dominant_isoform_id columns
#' @return Tibble with pairs
generate_pairs <- function(gtf_transcripts, dominant_mapping) {
  cat("\nGenerating pairs...\n")

  pairs_list <- list()

  for (i in seq_len(nrow(dominant_mapping))) {
    gene <- dominant_mapping$gene_id[i]
    dominant <- dominant_mapping$dominant_isoform_id[i]

    # Get all transcripts for this gene
    gene_transcripts <- gtf_transcripts %>%
      filter(gene_id == gene)

    if (nrow(gene_transcripts) == 0) {
      warning(sprintf("Gene %s not found in GTF", gene))
      next
    }

    # Check if dominant exists in GTF
    if (!dominant %in% gene_transcripts$transcript_id) {
      warning(sprintf("Dominant isoform %s not found in GTF for gene %s", dominant, gene))
      next
    }

    # Get gene metadata
    chr <- gene_transcripts$chr[1]
    strand <- gene_transcripts$strand[1]

    # Get all non-dominant isoforms
    non_dominant <- gene_transcripts %>%
      filter(transcript_id != dominant) %>%
      pull(transcript_id)

    if (length(non_dominant) == 0) {
      cat(sprintf("  Gene %s: Only 1 isoform (dominant), skipping\n", gene))
      next
    }

    # Create a pair for each non-dominant isoform
    for (comparator in non_dominant) {
      pairs_list[[length(pairs_list) + 1]] <- tibble(
        gene_id = gene,
        isoform_A = dominant,
        isoform_B = comparator,
        chr = chr,
        strand = strand
      )
    }

    cat(sprintf("  Gene %s: %d pairs (1 dominant vs %d non-dominant)\n",
                gene, length(non_dominant), length(non_dominant)))
  }

  if (length(pairs_list) == 0) {
    stop("No pairs generated. Check that gene IDs match between mapping and GTF.")
  }

  pairs <- bind_rows(pairs_list)

  cat(sprintf("\nGenerated %d total pairs from %d genes\n",
              nrow(pairs), n_distinct(pairs$gene_id)))

  return(pairs)
}

#' Main function
#'
#' @param gtf_file Path to GTF file
#' @param dominant_file Path to dominant mapping file
#' @param output_file Path to output pairs file
main <- function(gtf_file, dominant_file, output_file) {
  cat("\n╔════════════════════════════════════════════════════════════════╗\n")
  cat("║   Generate Pairs from Dominant Mapping                       ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")

  # Read dominant mapping
  cat("Reading dominant mapping:", dominant_file, "\n")
  dominant_mapping <- read_tsv(dominant_file,
                               col_types = cols(
                                 gene_id = col_character(),
                                 dominant_isoform_id = col_character()
                               ),
                               show_col_types = FALSE)

  cat(sprintf("  Loaded %d dominant isoform assignments\n\n", nrow(dominant_mapping)))

  # Check required columns
  if (!all(c("gene_id", "dominant_isoform_id") %in% colnames(dominant_mapping))) {
    stop("Dominant mapping file must have columns: gene_id, dominant_isoform_id")
  }

  # Parse GTF
  gtf_transcripts <- parse_gtf_transcripts(gtf_file)

  # Generate pairs
  pairs <- generate_pairs(gtf_transcripts, dominant_mapping)

  # Write output
  cat("\nWriting pairs file:", output_file, "\n")
  write_tsv(pairs, output_file)

  cat(sprintf("  ✓ Wrote %d pairs\n", nrow(pairs)))

  # Summary statistics
  cat("\n═══════════════════════════════════════════════════════════════════\n")
  cat("SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════════\n\n")

  pairs_per_gene <- pairs %>%
    count(gene_id) %>%
    pull(n)

  cat(sprintf("Genes: %d\n", n_distinct(pairs$gene_id)))
  cat(sprintf("Total pairs: %d\n", nrow(pairs)))
  cat(sprintf("Pairs per gene: min=%d, median=%d, max=%d\n",
              min(pairs_per_gene), median(pairs_per_gene), max(pairs_per_gene)))
  cat("\n")

  return(pairs)
}

# Command line interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 3) {
    cat("Usage: generate_pairs_from_dominant.R <gtf_file> <dominant_mapping> <output_pairs>\n")
    cat("\n")
    cat("Arguments:\n")
    cat("  gtf_file         GTF file with all isoforms\n")
    cat("  dominant_mapping TSV file with columns: gene_id, dominant_isoform_id\n")
    cat("  output_pairs     Output TSV file with pairs\n")
    cat("\n")
    quit(status = 1)
  }

  main(args[1], args[2], args[3])
}
