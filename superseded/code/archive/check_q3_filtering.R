#!/usr/bin/env Rscript
# Check Q3 upregulated transcript filtering

library(tidyverse)
library(data.table)

# Load DD DGE results
dd_dge <- read_csv("../longread_dge/nmd_dge_dd_2026.1.2.csv", show_col_types = FALSE)

cat("DD DGE file:\n")
cat("  Total transcripts:", nrow(dd_dge), "\n")

# Filter for upregulated (logFC > 0, adj.P.Val < 0.05)
dd_upregulated <- dd_dge %>%
  filter(logFC > 0, adj.P.Val < 0.05)

cat("  Upregulated (logFC > 0, adj.P < 0.05):", nrow(dd_upregulated), "\n\n")

# Load DD DMSO isoform proportions
dd_dmso <- fread("../tmp/isoform_proportions_dd_dmso_2026.1.3.tsv",
                 skip = 20, data.table = FALSE)

colnames(dd_dmso) <- c("cell_type", "treatment", "sample", "gene_id",
                        "transcript_id", "raw_count", "gene_total",
                        "isoform_proportion", "n_isoforms")

cat("DD DMSO isoform proportions:\n")
cat("  Total rows:", nrow(dd_dmso), "\n")
cat("  Unique transcripts:", n_distinct(dd_dmso$transcript_id), "\n")
cat("  Transcripts with raw_count > 0:",
    n_distinct(dd_dmso$transcript_id[dd_dmso$raw_count > 0]), "\n\n")

# Check overlap
upregulated_in_proportions <- dd_upregulated %>%
  filter(txid %in% dd_dmso$transcript_id)

cat("Upregulated transcripts found in isoform proportion file:\n")
cat("  Count:", nrow(upregulated_in_proportions), "\n\n")

upregulated_with_expression <- dd_upregulated %>%
  inner_join(
    dd_dmso %>% filter(raw_count > 0) %>% distinct(transcript_id),
    by = c("txid" = "transcript_id")
  )

cat("Upregulated transcripts with expression (raw_count > 0) in DMSO:\n")
cat("  Count:", nrow(upregulated_with_expression), "\n\n")

# Check what's missing
missing <- dd_upregulated %>%
  filter(!txid %in% dd_dmso$transcript_id)

cat("Upregulated transcripts NOT in isoform proportion file:\n")
cat("  Count:", nrow(missing), "\n")
if (nrow(missing) > 0) {
  cat("  Top 10 examples:\n")
  print(head(missing %>% select(txid, logFC, adj.P.Val, AveExpr), 10))

  cat("\n  Summary of missing transcripts:\n")
  cat("    Mean logFC:", round(mean(missing$logFC), 2), "\n")
  cat("    Mean AveExpr:", round(mean(missing$AveExpr), 2), "\n")
  cat("    Min AveExpr:", round(min(missing$AveExpr), 2), "\n")
}
