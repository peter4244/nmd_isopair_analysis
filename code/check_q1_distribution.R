#!/usr/bin/env Rscript
# Check Q1 delta_isoforms distribution

library(tidyverse)
library(data.table)

# Read Q1 results
q1_data <- fread("../tmp/isoform_interpretation_q1_isoforms_per_gene_2026.1.3.tsv",
                 data.table = FALSE)

cat("Q1 data dimensions:", nrow(q1_data), "rows\n\n")

cat("Summary of delta_isoforms:\n")
print(summary(q1_data$delta_isoforms))

cat("\n\nDistribution of delta_isoforms:\n")
table_delta <- table(q1_data$delta_isoforms)
print(head(sort(table_delta, decreasing = TRUE), 20))

cat("\n\nBreakdown by cell type:\n")
breakdown <- q1_data %>%
  group_by(cell_type) %>%
  summarize(
    n_genes = n(),
    min_delta = min(delta_isoforms),
    max_delta = max(delta_isoforms),
    mean_delta = mean(delta_isoforms),
    median_delta = median(delta_isoforms),
    n_positive = sum(delta_isoforms > 0),
    n_negative = sum(delta_isoforms < 0),
    n_zero = sum(delta_isoforms == 0)
  )
print(breakdown)
