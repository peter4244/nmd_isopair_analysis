library(tidyverse)

# Load expected labels
expected <- read_tsv("synthetic/TestData/annotations/base_events.tsv",
                     comment = "#", show_col_types = FALSE)

# Load GTF to see which tests have exon definitions
gtf <- read_tsv("synthetic/TestData/exons/base_events.gtf",
                col_names = c("seqnames", "source", "feature", "start", "end",
                             "score", "strand", "frame", "attributes"),
                comment = "#", show_col_types = FALSE)

gtf_genes <- gtf %>%
  mutate(gene_id = str_match(attributes, 'gene_id "([^"]+)"')[,2]) %>%
  pull(gene_id) %>%
  unique()

# Create status table
status <- expected %>%
  mutate(
    has_gtf = gene_id %in% gtf_genes,
    status = if_else(has_gtf, "✓", "✗ MISSING")
  ) %>%
  select(status, gene_id, event_type, description) %>%
  arrange(status, gene_id)

cat("\n╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║  Current Test Cases and Expected Labels                                 ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("Total tests: %d\n", nrow(status)))
cat(sprintf("With GTF exons: %d\n", sum(status$has_gtf)))
cat(sprintf("Missing GTF exons: %d\n\n", sum(!status$has_gtf)))

print(status, n = 100)
