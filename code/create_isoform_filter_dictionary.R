#!/usr/bin/env Rscript

# Create filtered isoform dictionary using 0.25 CPM threshold
# Filter: isoform must have >=0.25 CPM in at least one sample in EITHER DMSO OR Smg1i
# within each cell type

library(tidyverse)
library(edgeR)

cat("=== CREATING FILTERED ISOFORM DICTIONARY ===\n\n")
cat("Filter criteria: >= 0.25 CPM in at least one sample in EITHER DMSO OR Smg1i\n")
cat("Stratified by cell type\n\n")

# Load DGEList
cat("Loading DGEList...\n")
dge_isoform <- readRDS("../results/rds/dge_isoform_2026.1.18.rds")

# Extract normalized CPM
cpm_normalized <- cpm(dge_isoform, normalized = TRUE, log = FALSE)
cat("Loaded", nrow(cpm_normalized), "transcripts x", ncol(cpm_normalized), "samples\n\n")

# Convert to long format
cpm_df <- as.data.frame(cpm_normalized)
cpm_df$transcript_id <- dge_isoform$genes$txid

cpm_long <- cpm_df %>%
  pivot_longer(cols = -transcript_id, names_to = "bamid", values_to = "cpm") %>%
  left_join(dge_isoform$samples %>% select(bamid, sample, treatment, ct),
            by = "bamid") %>%
  rename(cell_type = ct)

# For each cell type and isoform, find max CPM in DMSO and Smg1i
isoform_max_cpm <- cpm_long %>%
  group_by(cell_type, transcript_id) %>%
  summarize(
    max_cpm_dmso = max(cpm[treatment == "DMSO"], na.rm = TRUE),
    max_cpm_smg1i = max(cpm[treatment == "Smg1i"], na.rm = TRUE),
    max_cpm_either = pmax(max_cpm_dmso, max_cpm_smg1i),
    .groups = "drop"
  )

# Apply filter
threshold <- 0.25
filtered_isoforms <- isoform_max_cpm %>%
  filter(max_cpm_either >= threshold) %>%
  select(cell_type, transcript_id, max_cpm_dmso, max_cpm_smg1i, max_cpm_either)

cat("Filtering results:\n\n")
filter_summary <- filtered_isoforms %>%
  group_by(cell_type) %>%
  summarize(
    n_isoforms_pass = n(),
    .groups = "drop"
  ) %>%
  left_join(
    isoform_max_cpm %>% group_by(cell_type) %>% summarize(n_total = n(), .groups = "drop"),
    by = "cell_type"
  ) %>%
  mutate(
    n_excluded = n_total - n_isoforms_pass,
    pct_retained = 100 * n_isoforms_pass / n_total
  ) %>%
  select(cell_type, n_total, n_isoforms_pass, n_excluded, pct_retained)

print(filter_summary)
cat("\n")

# Check against significant isoforms from longread_dge
cat("Comparing to longread_dge significant results...\n\n")

dge_files <- list.files("../longread_dge", pattern = "^nmd_dge_.*\\.csv$", full.names = TRUE)

all_dge <- map_df(dge_files, ~{
  df <- read_csv(.x, show_col_types = FALSE)
  basename <- basename(.x)
  parts <- str_split(basename, "_")[[1]]
  celltype <- parts[3]
  df$cell_type_file <- celltype
  df
})

celltype_map <- c(
  "at2" = "AT",
  "dd" = "DD",
  "ddali" = "DD_ALI",
  "doali" = "DO",
  "fb" = "FB",
  "mv" = "MV"
)

all_dge$cell_type <- celltype_map[all_dge$cell_type_file]
all_dge <- all_dge %>% filter(!is.na(cell_type))

sig_isoforms <- all_dge %>%
  filter(adj.P.Val < 0.05) %>%
  select(cell_type, transcript_id = txid, logFC, adj.P.Val)

# Check coverage
sig_coverage <- sig_isoforms %>%
  left_join(
    filtered_isoforms %>% select(cell_type, transcript_id) %>% mutate(in_dictionary = TRUE),
    by = c("cell_type", "transcript_id")
  ) %>%
  mutate(in_dictionary = replace_na(in_dictionary, FALSE))

coverage_summary <- sig_coverage %>%
  group_by(cell_type) %>%
  summarize(
    n_significant = n(),
    n_in_dictionary = sum(in_dictionary),
    n_excluded = sum(!in_dictionary),
    pct_excluded = 100 * sum(!in_dictionary) / n(),
    .groups = "drop"
  )

cat("Significant isoforms coverage:\n\n")
print(coverage_summary)
cat("\n")

# Characteristics of excluded significant isoforms
excluded_sig <- sig_coverage %>%
  filter(!in_dictionary) %>%
  group_by(cell_type) %>%
  summarize(
    n_excluded = n(),
    mean_logFC = mean(logFC, na.rm = TRUE),
    median_logFC = median(logFC, na.rm = TRUE),
    mean_adj_pval = mean(adj.P.Val, na.rm = TRUE),
    n_upregulated = sum(logFC > 0),
    n_downregulated = sum(logFC < 0),
    .groups = "drop"
  )

cat("Characteristics of excluded significant isoforms:\n\n")
print(excluded_sig)
cat("\n")

# Save the dictionary
cat("Saving filtered isoform dictionary...\n")
write_tsv(filtered_isoforms, "../tmp/filtered_isoform_dictionary_0.25cpm_2026.1.3.tsv")

# Also save as an RDS for faster loading
saveRDS(filtered_isoforms, "../tmp/filtered_isoform_dictionary_0.25cpm_2026.1.3.rds")

# Save summary statistics
write_tsv(filter_summary, "../tmp/isoform_filter_summary_2026.1.3.tsv")
write_tsv(coverage_summary, "../tmp/isoform_filter_sig_coverage_2026.1.3.tsv")
write_tsv(excluded_sig, "../tmp/isoform_filter_excluded_sig_2026.1.3.tsv")

cat("\nSaved:\n")
cat("  - tmp/filtered_isoform_dictionary_0.25cpm_2026.1.3.tsv (TSV format)\n")
cat("  - tmp/filtered_isoform_dictionary_0.25cpm_2026.1.3.rds (RDS format)\n")
cat("  - tmp/isoform_filter_summary_2026.1.3.tsv\n")
cat("  - tmp/isoform_filter_sig_coverage_2026.1.3.tsv\n")
cat("  - tmp/isoform_filter_excluded_sig_2026.1.3.tsv\n\n")

cat("=== SUMMARY ===\n\n")
cat("Total isoforms in dictionary:", nrow(filtered_isoforms), "\n")
cat("Unique isoforms:", n_distinct(filtered_isoforms$transcript_id), "\n")
cat("Significant isoforms retained:", sum(coverage_summary$n_in_dictionary),
    "out of", sum(coverage_summary$n_significant),
    "(", round(100 * sum(coverage_summary$n_in_dictionary) / sum(coverage_summary$n_significant), 2), "%)\n")
