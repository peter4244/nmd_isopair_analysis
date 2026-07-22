#!/usr/bin/env Rscript
# Investigate Unclassifiable Genes
# Show concrete examples of overlapping exons with no shared boundaries

library(tidyverse)

# Load data
exon_structures <- readRDS('results/isoform_transitions/v4.0_reference_based/major_isoforms/exon_structures_major_isoforms.rds')
filtered_genes <- read_tsv('results/isoform_transitions/v4.0_reference_based/major_isoforms/filtered_genes_major.tsv', show_col_types = FALSE)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   INVESTIGATING UNCLASSIFIABLE GENES                          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Get unclassifiable genes
unclass_genes <- filtered_genes %>%
  filter(status == "unclassifiable") %>%
  pull(gene_id)

cat(sprintf("Total unclassifiable genes: %d\n\n", length(unclass_genes)))

# Function to find overlapping exons with no shared boundaries
find_unclassifiable_patterns <- function(gene_id) {
  # Unnest exon coordinates
  gene_exons <- exon_structures %>%
    filter(gene_id == !!gene_id) %>%
    select(isoform_id, strand, exon_starts, exon_ends) %>%
    unnest(cols = c(exon_starts, exon_ends)) %>%
    rename(start = exon_starts, end = exon_ends) %>%
    distinct() %>%
    arrange(start, end)

  if (nrow(gene_exons) < 2) return(NULL)

  patterns <- list()

  for (i in 1:(nrow(gene_exons)-1)) {
    for (j in (i+1):nrow(gene_exons)) {
      ex_i <- gene_exons[i,]
      ex_j <- gene_exons[j,]

      # Check overlap
      overlaps <- (ex_i$start < ex_j$end) & (ex_j$start < ex_i$end)

      # Check boundary sharing
      shares_start <- ex_i$start == ex_j$start
      shares_end <- ex_i$end == ex_j$end

      # Unclassifiable pattern
      if (overlaps && !shares_start && !shares_end) {
        patterns[[length(patterns) + 1]] <- tibble(
          exon_1 = sprintf("%s: %d-%d", ex_i$isoform_id, ex_i$start, ex_i$end),
          exon_2 = sprintf("%s: %d-%d", ex_j$isoform_id, ex_j$start, ex_j$end),
          overlap_start = max(ex_i$start, ex_j$start),
          overlap_end = min(ex_i$end, ex_j$end),
          overlap_length = min(ex_i$end, ex_j$end) - max(ex_i$start, ex_j$start)
        )
      }
    }
  }

  if (length(patterns) > 0) {
    bind_rows(patterns) %>% mutate(gene_id = gene_id)
  } else {
    NULL
  }
}

# Examine first 5 unclassifiable genes
cat("═══════════════════════════════════════════════════════════════\n")
cat("EXAMPLES OF UNCLASSIFIABLE PATTERNS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (i in 1:min(5, length(unclass_genes))) {
  gene_id <- unclass_genes[i]

  cat(sprintf("GENE %d: %s\n", i, gene_id))
  cat("───────────────────────────────────────────────────────────────\n")

  # Get gene info
  gene_exons <- exon_structures %>%
    filter(gene_id == !!gene_id) %>%
    select(isoform_id, strand, exon_starts, exon_ends) %>%
    unnest(cols = c(exon_starts, exon_ends)) %>%
    rename(start = exon_starts, end = exon_ends) %>%
    distinct() %>%
    arrange(start, end)

  cat(sprintf("Number of unique exons: %d\n", nrow(gene_exons)))
  cat(sprintf("Isoforms: %d\n\n", n_distinct(gene_exons$isoform_id)))

  # Find unclassifiable patterns
  patterns <- find_unclassifiable_patterns(gene_id)

  if (!is.null(patterns)) {
    cat(sprintf("Found %d overlapping exon pairs with no shared boundaries:\n\n", nrow(patterns)))

    for (j in 1:min(3, nrow(patterns))) {
      p <- patterns[j,]
      cat(sprintf("  Pattern %d:\n", j))
      cat(sprintf("    Exon 1: %s\n", p$exon_1))
      cat(sprintf("    Exon 2: %s\n", p$exon_2))
      cat(sprintf("    Overlap: %d bp (%d-%d)\n", p$overlap_length, p$overlap_start, p$overlap_end))
      cat("\n")
    }

    if (nrow(patterns) > 3) {
      cat(sprintf("  ... and %d more patterns\n\n", nrow(patterns) - 3))
    }
  } else {
    cat("No unclassifiable patterns found (odd - should have been caught)\n\n")
  }

  # Show all exons for this gene
  cat("All exons in gene:\n")
  print(gene_exons %>% head(15))

  if (nrow(gene_exons) > 15) {
    cat(sprintf("... and %d more exons\n", nrow(gene_exons) - 15))
  }

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
}

# Summary statistics
cat("\n")
cat("SUMMARY ANALYSIS OF UNCLASSIFIABLE PATTERNS\n")
cat("───────────────────────────────────────────────────────────────\n")

all_patterns <- map_dfr(head(unclass_genes, 50), find_unclassifiable_patterns)

if (!is.null(all_patterns) && nrow(all_patterns) > 0) {
  cat(sprintf("\nAnalyzed first 50 unclassifiable genes:\n"))
  cat(sprintf("  Total overlapping pairs found: %d\n", nrow(all_patterns)))
  cat(sprintf("  Mean overlap length: %.0f bp\n", mean(all_patterns$overlap_length)))
  cat(sprintf("  Median overlap length: %.0f bp\n", median(all_patterns$overlap_length)))
  cat(sprintf("  Range: %d - %d bp\n", min(all_patterns$overlap_length), max(all_patterns$overlap_length)))

  cat("\n")
  cat("Overlap length distribution:\n")
  hist_breaks <- c(0, 10, 50, 100, 500, 1000, 5000, Inf)
  hist_labels <- c("1-10", "11-50", "51-100", "101-500", "501-1000", "1001-5000", ">5000")

  all_patterns <- all_patterns %>%
    mutate(overlap_bin = cut(overlap_length, breaks = hist_breaks, labels = hist_labels))

  overlap_dist <- all_patterns %>%
    count(overlap_bin) %>%
    mutate(pct = 100 * n / sum(n))

  print(overlap_dist)
}

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   INVESTIGATION COMPLETE                                       ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
