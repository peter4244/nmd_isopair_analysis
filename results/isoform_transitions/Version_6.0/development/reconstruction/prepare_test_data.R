#!/usr/bin/env Rscript
#
# Prepare Test Data
#
# Step 1: Validate GTF meets specifications
# Step 2: Randomly select 1 dominant isoform per gene
#
# Usage: Rscript prepare_test_data.R <input_gtf> <output_dominant_mapping>
#

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("Usage: prepare_test_data.R <input_gtf> <output_dominant_mapping>\n")
  quit(status = 1)
}

input_gtf <- args[1]
output_mapping <- args[2]

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   Prepare Test Data                                           ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# ==============================================================================
# Step 1: Validate GTF
# ==============================================================================

cat("STEP 1: VALIDATE GTF\n")
cat("════════════════════════════════════════════════════════════════\n")
cat(sprintf("Reading: %s\n", input_gtf))

lines <- readLines(input_gtf)
lines <- lines[!startsWith(lines, '#') & lines != '']

cat(sprintf("  Total lines: %d\n", length(lines)))

# Parse exons
exons <- list()
for (line in lines) {
  parts <- strsplit(line, '\t')[[1]]

  if (length(parts) < 9) {
    cat(sprintf("  ⚠ WARNING: Line with < 9 fields: %s\n", substr(line, 1, 50)))
    next
  }

  # Only parse exon features
  if (parts[3] != "exon") {
    next
  }

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

cat(sprintf("  Exon records: %d\n", nrow(exons_df)))
cat(sprintf("  Genes: %d\n", n_distinct(exons_df$gene_id)))
cat(sprintf("  Transcripts: %d\n", n_distinct(exons_df$transcript_id)))

# Check for required attributes
if (any(is.na(exons_df$gene_id))) {
  stop("ERROR: Some exons missing gene_id attribute")
}
if (any(is.na(exons_df$transcript_id))) {
  stop("ERROR: Some exons missing transcript_id attribute")
}

# Check for valid coordinates
if (any(exons_df$start > exons_df$end)) {
  stop("ERROR: Some exons have start > end")
}

# Check for valid strands
valid_strands <- c("+", "-")
if (!all(exons_df$strand %in% valid_strands)) {
  stop("ERROR: Some exons have invalid strand (must be + or -)")
}

cat("  ✓ GTF validation passed\n\n")

# ==============================================================================
# Step 2: Select Dominant Isoforms
# ==============================================================================

cat("STEP 2: SELECT DOMINANT ISOFORMS\n")
cat("════════════════════════════════════════════════════════════════\n")

# Get genes with multiple isoforms
gene_isoform_counts <- exons_df %>%
  group_by(gene_id) %>%
  summarize(
    n_isoforms = n_distinct(transcript_id),
    chr = first(chr),
    strand = first(strand),
    .groups = "drop"
  ) %>%
  filter(n_isoforms >= 2)  # Only genes with 2+ isoforms

cat(sprintf("  Genes with 2+ isoforms: %d\n", nrow(gene_isoform_counts)))

# For each gene, randomly select one isoform as dominant
set.seed(2026)

dominant_mapping <- list()
for (i in seq_len(nrow(gene_isoform_counts))) {
  gene <- gene_isoform_counts$gene_id[i]

  # Get all isoforms for this gene
  isoforms <- exons_df %>%
    filter(gene_id == gene) %>%
    pull(transcript_id) %>%
    unique()

  # Randomly select one as dominant
  dominant <- sample(isoforms, 1)

  dominant_mapping[[i]] <- tibble(
    gene_id = gene,
    dominant_isoform_id = dominant
  )
}

dominant_df <- bind_rows(dominant_mapping)

cat(sprintf("  Selected %d dominant isoforms\n\n", nrow(dominant_df)))

# Write output
cat(sprintf("Writing: %s\n", output_mapping))
write_tsv(dominant_df, output_mapping)

cat("\n✓ COMPLETE\n")
cat(sprintf("  Output: %s\n", output_mapping))
cat(sprintf("  Genes: %d\n", nrow(dominant_df)))
