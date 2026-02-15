#!/usr/bin/env Rscript
library(tidyverse)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   Sample Real Data for Event Detection Testing               ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Paths
sqanti_gtf <- "/Users/petecastaldi/claude_projects/nmd/results/sqanti_runs/merged_collapsed/isoforms_corrected.gtf"
output_gtf <- "real_data/sampled_isoforms.gtf"
output_pairs <- "real_data/sampled_pairs.tsv"

# Create output directory
dir.create("real_data", showWarnings = FALSE)

cat("Step 1: Loading GTF file...\n")
cat(sprintf("  Reading: %s\n", sqanti_gtf))

# Read GTF (only exon features)
gtf <- read_tsv(
  sqanti_gtf,
  col_names = c("seqnames", "source", "feature", "start", "end",
                "score", "strand", "frame", "attributes"),
  col_types = "ccciicccc",
  comment = "#"
) %>%
  filter(feature == "exon")

cat(sprintf("  Loaded %s exon records\n\n", format(nrow(gtf), big.mark = ",")))

# Extract gene_id and transcript_id from attributes
cat("Step 2: Parsing gene and transcript IDs...\n")
gtf <- gtf %>%
  mutate(
    gene_id = str_match(attributes, 'gene_id "([^"]+)"')[, 2],
    transcript_id = str_match(attributes, 'transcript_id "([^"]+)"')[, 2]
  ) %>%
  filter(!is.na(gene_id), !is.na(transcript_id))

# Count isoforms per gene
cat("Step 3: Finding genes with multiple isoforms...\n")
gene_isoform_counts <- gtf %>%
  group_by(gene_id) %>%
  summarise(
    n_isoforms = n_distinct(transcript_id),
    strand = first(strand),
    .groups = "drop"
  ) %>%
  filter(n_isoforms >= 2) %>%
  arrange(desc(n_isoforms))

cat(sprintf("  Found %s genes with 2+ isoforms\n", format(nrow(gene_isoform_counts), big.mark = ",")))
cat(sprintf("  Isoforms per gene: median=%d, max=%d\n\n",
            median(gene_isoform_counts$n_isoforms),
            max(gene_isoform_counts$n_isoforms)))

# Sample 20 genes with diverse isoform counts
cat("Step 4: Sampling 20 genes...\n")
set.seed(42)

# Stratify by isoform count to get diversity
sampled_genes <- gene_isoform_counts %>%
  mutate(isoform_bin = cut(n_isoforms,
                            breaks = c(2, 3, 5, 10, Inf),
                            labels = c("2", "3-4", "5-9", "10+"))) %>%
  group_by(isoform_bin) %>%
  slice_sample(n = 5) %>%
  ungroup() %>%
  select(gene_id, n_isoforms, strand)

cat("  Sampled genes:\n")
print(sampled_genes, n = 20)

# For each gene, randomly select 2 isoforms
cat("\nStep 5: Selecting isoform pairs...\n")
pairs <- list()
for (i in seq_len(nrow(sampled_genes))) {
  gene <- sampled_genes$gene_id[i]

  # Get all isoforms for this gene
  gene_isoforms <- gtf %>%
    filter(gene_id == gene) %>%
    pull(transcript_id) %>%
    unique()

  # Sample 2 isoforms
  if (length(gene_isoforms) >= 2) {
    selected <- sample(gene_isoforms, 2)
    pairs[[i]] <- tibble(
      gene_id = gene,
      isoform_A = selected[1],
      isoform_B = selected[2],
      n_isoforms_in_gene = sampled_genes$n_isoforms[i],
      strand = sampled_genes$strand[i]
    )
  }
}
pairs_df <- bind_rows(pairs)

cat(sprintf("  Created %d isoform pairs\n\n", nrow(pairs_df)))

# Extract GTF records for selected ISOFORMS (not all isoforms in the gene)
cat("Step 6: Extracting GTF records for selected isoform pairs...\n")
selected_transcripts <- c(pairs_df$isoform_A, pairs_df$isoform_B)
selected_gtf <- gtf %>%
  filter(transcript_id %in% selected_transcripts)

cat(sprintf("  Extracted %s exon records\n\n", format(nrow(selected_gtf), big.mark = ",")))

# Write outputs
cat("Step 7: Writing output files...\n")

# Write GTF by reading original lines and filtering (preserves exact format)
cat("  Extracting original GTF lines for selected pairs only...\n")
gtf_lines <- readLines(sqanti_gtf)

# Get transcript IDs for the pairs (should be exactly 2 per gene)
keep_transcripts <- selected_transcripts

# Filter lines that contain any of our transcript IDs
selected_lines <- gtf_lines[grepl(paste0("transcript_id \"(",
                                         paste(keep_transcripts, collapse = "|"),
                                         ")\""), gtf_lines)]

writeLines(selected_lines, output_gtf)
cat(sprintf("  ✓ GTF: %s (%d lines)\n", output_gtf, length(selected_lines)))

# Write pairs TSV
pairs_df %>%
  write_tsv(output_pairs)
cat(sprintf("  ✓ Pairs: %s\n", output_pairs))

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat("SAMPLING COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")
cat(sprintf("Sampled: %d genes, %d isoform pairs\n", nrow(sampled_genes), nrow(pairs_df)))
cat(sprintf("GTF size: %.1f MB\n", file.size(output_gtf) / 1e6))
cat("\nNext steps:\n")
cat("  1. Run: Rscript test_real_data.R\n")
cat("  2. Review: real_data/real_data_validation.pdf\n")
