#!/usr/bin/env Rscript
# Quick check of summary statistics from generated files

library(tidyverse)
library(data.table)

# Read all isoform proportion files
files <- list.files("../tmp", pattern = "^isoform_proportions_.*\\.tsv$", full.names = TRUE)

# Read all files, skipping headers
all_data <- map_df(files, ~fread(.x, skip = 19, data.table = FALSE))

# Set column names
colnames(all_data) <- c("cell_type", "treatment", "sample", "gene_id",
                        "transcript_id", "raw_count", "gene_total",
                        "isoform_proportion", "n_isoforms")

# Calculate summary stats (only counting transcripts with expression)
summary_stats <- all_data %>%
  filter(raw_count > 0) %>%  # Only count transcripts with actual expression
  group_by(cell_type, treatment) %>%
  summarize(
    n_samples = n_distinct(sample),
    n_genes = n_distinct(gene_id),
    n_transcripts = n_distinct(transcript_id),
    mean_isoforms_per_gene = n_distinct(transcript_id) / n_distinct(gene_id),
    .groups = "drop"
  ) %>%
  arrange(cell_type, treatment)

print(summary_stats, n = 20)
