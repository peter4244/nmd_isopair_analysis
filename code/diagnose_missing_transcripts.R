#!/usr/bin/env Rscript
# Diagnose why DD upregulated transcripts are missing

library(tidyverse)
library(data.table)

# Load DD DGE
dd_dge <- read_csv("../longread_dge/nmd_dge_dd_2026.1.2.csv", show_col_types = FALSE)
dd_upregulated <- dd_dge %>%
  filter(logFC > 0, adj.P.Val < 0.05)

cat("DD upregulated transcripts:", nrow(dd_upregulated), "\n\n")

# Load ALL isoform proportion files to check filter across all cell types
files <- list.files("../tmp", pattern = "^isoform_proportions_.*\\.tsv$", full.names = TRUE)
all_props <- map_df(files, ~fread(.x, skip = 20, data.table = FALSE))
colnames(all_props) <- c("cell_type", "treatment", "sample", "gene_id",
                          "transcript_id", "raw_count", "gene_total",
                          "isoform_proportion", "n_isoforms")

cat("Loaded proportions from", length(files), "files\n")
cat("Total unique transcripts in proportion files:", n_distinct(all_props$transcript_id), "\n\n")

# Check which DD upregulated transcripts are in the filtered data
in_filtered <- dd_upregulated %>%
  filter(txid %in% all_props$transcript_id)

not_in_filtered <- dd_upregulated %>%
  filter(!txid %in% all_props$transcript_id)

cat("DD upregulated in filtered data:", nrow(in_filtered), "\n")
cat("DD upregulated NOT in filtered data:", nrow(not_in_filtered), "\n\n")

# Now check the RAW counts file to see if these transcripts exist at all
cat("Loading raw counts file...\n")
raw_counts <- fread("../data/qdf_raw_counts_2026.1.3.csv", data.table = FALSE)

# Get transcript IDs from first column
transcript_ids_raw <- raw_counts[[1]]

cat("Raw counts file has", length(transcript_ids_raw), "transcripts\n")

# Check if missing transcripts are in raw counts
missing_in_raw <- not_in_filtered %>%
  filter(!txid %in% transcript_ids_raw)

cat("Missing transcripts not in raw counts file:", nrow(missing_in_raw), "\n")

# Check filter stats for a sample of missing transcripts
if (nrow(not_in_filtered) > 0 && any(not_in_filtered$txid %in% transcript_ids_raw)) {
  cat("\nChecking filter statistics for missing transcripts...\n")

  # Get sample columns
  sample_cols <- colnames(raw_counts)[3:ncol(raw_counts)]

  # Check first 10 missing transcripts that ARE in raw counts
  missing_in_raw_data <- not_in_filtered %>%
    filter(txid %in% transcript_ids_raw) %>%
    head(10)

  for (i in 1:min(10, nrow(missing_in_raw_data))) {
    tx <- missing_in_raw_data$txid[i]
    tx_data <- raw_counts[transcript_ids_raw == tx, sample_cols]

    cat(sprintf("\n%s (logFC=%.2f, adj.P=%.2e):\n",
                tx, missing_in_raw_data$logFC[i], missing_in_raw_data$adj.P.Val[i]))
    cat("  Max count across all samples:", max(unlist(tx_data)), "\n")
    cat("  Samples with ≥2 counts:", sum(unlist(tx_data) >= 2), "\n")
    cat("  Total counts:", sum(unlist(tx_data)), "\n")
  }
}
