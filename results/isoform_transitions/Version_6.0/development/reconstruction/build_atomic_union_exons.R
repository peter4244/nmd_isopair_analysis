#!/usr/bin/env Rscript
#
# Build Atomic Union Exons
#
# Creates atomic (non-overlapping) union exon segments by splitting at all boundaries.
# Each atomic segment represents the smallest unit that is either fully included
# or fully excluded in any given isoform.

library(tidyverse)

#' Build atomic union exons for a gene
#'
#' @param gtf_file Path to GTF file
#' @param output_file Path to output union exons file
build_atomic_union_exons_from_gtf <- function(gtf_file, output_file) {
  cat("\n╔════════════════════════════════════════════════════════════════╗\n")
  cat("║   Build Atomic Union Exons                                    ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")

  # Parse GTF
  cat("Reading GTF file:", gtf_file, "\n")
  lines <- readLines(gtf_file)
  lines <- lines[!startsWith(lines, '#') & lines != '']

  exons <- list()
  for (line in lines) {
    parts <- strsplit(line, '\t')[[1]]
    if (parts[3] != "exon") next

    attrs_str <- parts[9]
    gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', attrs_str)
    transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', attrs_str)

    exons[[length(exons) + 1]] <- tibble(
      chr = parts[1],
      start = as.integer(parts[4]),
      end = as.integer(parts[5]),
      strand = parts[7],
      gene_id = gene_id,
      transcript_id = transcript_id
    )
  }

  exons_df <- bind_rows(exons)
  cat(sprintf("  Parsed %d exons from %d genes\n",
              nrow(exons_df), n_distinct(exons_df$gene_id)))

  # Build atomic union exons per gene
  cat("\nBuilding atomic union exons per gene...\n")
  genes <- unique(exons_df$gene_id)

  all_union_exons <- list()
  ue_counter <- 1

  for (gene in genes) {
    gene_exons <- exons_df %>% filter(gene_id == !!gene)

    # Get unique chr and strand for this gene
    chr <- unique(gene_exons$chr)[1]
    strand <- unique(gene_exons$strand)[1]

    # Collect all unique boundaries and track their types (START vs END)
    # In GTF: coordinates are 1-based inclusive [start, end]
    # An exon 100-200 contains bases 100, 101, ..., 200
    all_boundaries <- sort(unique(c(gene_exons$start, gene_exons$end)))

    if (length(all_boundaries) < 2) next

    # Determine boundary types: START, END, or BOTH
    boundary_types <- sapply(all_boundaries, function(b) {
      is_start <- b %in% gene_exons$start
      is_end <- b %in% gene_exons$end
      if (is_start && is_end) return("BOTH")
      if (is_start) return("START")
      return("END")
    })

    # Create atomic segments using boundary-type-aware splitting:
    # - START boundary (middle): split before it → [..., b-1] then [b, ...]
    # - END boundary (middle): split after it → [..., b] then [b+1, ...]
    # - First boundary: always segment start
    # - Last boundary: always segment end

    n <- length(all_boundaries)
    seg_starts <- c()
    seg_ends <- c()

    current_start <- all_boundaries[1]

    for (i in 2:(n-1)) {
      b <- all_boundaries[i]
      if (boundary_types[i] == "START") {
        # Split before this START boundary
        seg_starts <- c(seg_starts, current_start)
        seg_ends <- c(seg_ends, b - 1)
        current_start <- b
      } else if (boundary_types[i] == "END") {
        # Split after this END boundary
        seg_starts <- c(seg_starts, current_start)
        seg_ends <- c(seg_ends, b)
        current_start <- b + 1
      } else {
        # BOTH: this is a transition point (end of one exon, start of another)
        # Split at the boundary: end current segment at b, start next at b+1
        seg_starts <- c(seg_starts, current_start)
        seg_ends <- c(seg_ends, b)
        current_start <- b + 1
      }
    }

    # Add final segment
    seg_starts <- c(seg_starts, current_start)
    seg_ends <- c(seg_ends, all_boundaries[n])

    segments <- tibble(
      seg_start = seg_starts,
      seg_end = seg_ends
    )

    # Filter to segments that are actually covered by at least one exon
    covered_segments <- segments %>%
      rowwise() %>%
      filter({
        # Check if this segment is covered by any exon
        any(
          gene_exons$start <= seg_start &
          gene_exons$end >= seg_end
        )
      }) %>%
      ungroup()

    # Create union exons from covered segments
    for (i in seq_len(nrow(covered_segments))) {
      all_union_exons[[ue_counter]] <- tibble(
        chr = chr,
        start = covered_segments$seg_start[i],
        end = covered_segments$seg_end[i],
        union_exon_id = sprintf("ue_%d", ue_counter),
        strand = strand,
        gene_id = gene
      )
      ue_counter <- ue_counter + 1
    }
  }

  union_exons <- bind_rows(all_union_exons)

  cat(sprintf("  Created %d atomic union exons\n", nrow(union_exons)))

  # Write output
  cat("\nWriting output:", output_file, "\n")

  # Sort by chr, start for tabix
  union_exons <- union_exons %>%
    arrange(chr, start)

  # Write with header comment (# prefix so tabix skips it)
  write_lines("#chr\tstart\tend\tunion_exon_id\tstrand\tgene_id", output_file)
  write_tsv(union_exons, output_file, append = TRUE, col_names = FALSE)

  # Compress with bgzip
  cat("Compressing with bgzip...\n")
  bgzip_status <- system(sprintf("bgzip -f %s", output_file))
  if (bgzip_status != 0) {
    stop(sprintf("bgzip failed with status %d", bgzip_status))
  }

  # Index with tabix
  cat("Indexing with tabix...\n")
  tabix_status <- system(sprintf("tabix -f -s 1 -b 2 -e 3 %s.gz", output_file))
  if (tabix_status != 0) {
    stop(sprintf("tabix failed with status %d", tabix_status))
  }

  cat(sprintf("\n✓ Complete: %s.gz (indexed)\n\n", output_file))
  return(union_exons)
}

# Command line interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 2) {
    cat("Usage: build_atomic_union_exons.R <gtf_file> <output_file>\n")
    quit(status = 1)
  }

  build_atomic_union_exons_from_gtf(args[1], args[2])
}
