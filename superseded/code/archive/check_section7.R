#!/usr/bin/env Rscript
# Extract Section 7 table: Preview Results

library(tidyverse)
library(data.table)

# Read DD DMSO file
dd_dmso <- fread("../tmp/isoform_proportions_dd_dmso_2026.1.3.tsv",
                 skip = 19, data.table = FALSE)

colnames(dd_dmso) <- c("cell_type", "treatment", "sample", "gene_id",
                        "transcript_id", "raw_count", "gene_total",
                        "isoform_proportion", "n_isoforms")

cat("DD + DMSO results:", nrow(dd_dmso), "rows\n")
cat("Unique samples:", paste(unique(dd_dmso$sample), collapse = ", "), "\n")
cat("Unique genes:", n_distinct(dd_dmso$gene_id), "\n")
cat("Unique transcripts:", n_distinct(dd_dmso$transcript_id), "\n\n")

# Show genes with multiple isoforms
multi_isoform <- dd_dmso %>%
  filter(n_isoforms > 1) %>%
  arrange(desc(n_isoforms), gene_id, sample, desc(isoform_proportion))

cat("Rows with multiple isoforms:", nrow(multi_isoform), "\n\n")

# Show top 50 rows
cat("Top 50 rows:\n")
print(head(multi_isoform, 50), row.names = FALSE)

# Show summary of genes with most isoforms
cat("\n\nGenes with most isoforms:\n")
top_genes <- multi_isoform %>%
  distinct(gene_id, n_isoforms) %>%
  arrange(desc(n_isoforms)) %>%
  head(20)
print(top_genes, row.names = FALSE)
