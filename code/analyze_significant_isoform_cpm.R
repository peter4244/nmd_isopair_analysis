#!/usr/bin/env Rscript

# Analyze CPM distribution for significant isoforms from longread_dge
# Stratified by cell type

library(tidyverse)
library(edgeR)

# Load longread DGE results
cat("Loading longread DGE results...\n")
dge_files <- list.files("../longread_dge", pattern = "^nmd_dge_.*\\.csv$", full.names = TRUE)

all_dge <- map_df(dge_files, ~{
  df <- read_csv(.x, show_col_types = FALSE)

  # Extract cell type from filename
  basename <- basename(.x)
  parts <- str_split(basename, "_")[[1]]
  celltype <- parts[3]

  df$cell_type_file <- celltype
  df
})

# Map cell type codes
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

cat("Loaded", nrow(all_dge), "total DGE results\n\n")

# Identify significant isoforms (adj.P.Val < 0.05)
sig_isoforms <- all_dge %>%
  filter(adj.P.Val < 0.05) %>%
  select(cell_type, transcript_id = txid, logFC, adj.P.Val, AveExpr)

cat("Significant isoforms by cell type:\n")
sig_summary <- sig_isoforms %>%
  group_by(cell_type) %>%
  summarize(
    n_significant = n(),
    n_upregulated = sum(logFC > 0),
    n_downregulated = sum(logFC < 0)
  )
print(sig_summary)
cat("\n")

# Load DGEList for normalized CPM values
cat("Loading DGEList...\n")
dge_isoform <- readRDS("../results/rds/dge_isoform_2026.1.3.rds")

# Extract normalized CPM
cpm_normalized <- cpm(dge_isoform, normalized = TRUE, log = FALSE)
cat("Calculated normalized CPM:", nrow(cpm_normalized), "transcripts x", ncol(cpm_normalized), "samples\n\n")

# Convert CPM matrix to long format
cpm_df <- as.data.frame(cpm_normalized)
cpm_df$transcript_id <- dge_isoform$genes$txid

cpm_long <- cpm_df %>%
  pivot_longer(cols = -transcript_id, names_to = "bamid", values_to = "cpm") %>%
  left_join(dge_isoform$samples %>% select(bamid, sample, treatment, ct),
            by = "bamid") %>%
  rename(cell_type = ct)

cat("CPM long format:", nrow(cpm_long), "rows\n\n")

# Match significant isoforms to CPM values
sig_cpm <- sig_isoforms %>%
  left_join(cpm_long, by = c("cell_type", "transcript_id"))

cat("Matched significant isoforms to CPM values\n")
cat("Total rows:", nrow(sig_cpm), "\n\n")

# Calculate distribution statistics by cell type and treatment
cat("=== CPM DISTRIBUTION FOR SIGNIFICANT ISOFORMS ===\n\n")

distribution_stats <- sig_cpm %>%
  group_by(cell_type, treatment) %>%
  summarize(
    n_isoforms = n_distinct(transcript_id),
    n_observations = n(),
    mean_cpm = mean(cpm, na.rm = TRUE),
    median_cpm = median(cpm, na.rm = TRUE),
    sd_cpm = sd(cpm, na.rm = TRUE),
    q25_cpm = quantile(cpm, 0.25, na.rm = TRUE),
    q75_cpm = quantile(cpm, 0.75, na.rm = TRUE),
    min_cpm = min(cpm, na.rm = TRUE),
    max_cpm = max(cpm, na.rm = TRUE),
    n_below_5 = sum(cpm < 5, na.rm = TRUE),
    pct_below_5 = 100 * sum(cpm < 5, na.rm = TRUE) / n(),
    n_zero = sum(cpm == 0, na.rm = TRUE),
    pct_zero = 100 * sum(cpm == 0, na.rm = TRUE) / n(),
    .groups = "drop"
  )

print(distribution_stats)
cat("\n")

# Also calculate per-isoform average CPM (averaging across samples within each treatment/cell type)
isoform_avg_cpm <- sig_cpm %>%
  group_by(cell_type, treatment, transcript_id) %>%
  summarize(
    avg_cpm = mean(cpm, na.rm = TRUE),
    .groups = "drop"
  )

# Distribution of per-isoform averages
cat("=== PER-ISOFORM AVERAGE CPM DISTRIBUTION ===\n\n")

isoform_distribution <- isoform_avg_cpm %>%
  group_by(cell_type, treatment) %>%
  summarize(
    n_isoforms = n(),
    mean_avg_cpm = mean(avg_cpm, na.rm = TRUE),
    median_avg_cpm = median(avg_cpm, na.rm = TRUE),
    q25_avg_cpm = quantile(avg_cpm, 0.25, na.rm = TRUE),
    q75_avg_cpm = quantile(avg_cpm, 0.75, na.rm = TRUE),
    min_avg_cpm = min(avg_cpm, na.rm = TRUE),
    max_avg_cpm = max(avg_cpm, na.rm = TRUE),
    n_avg_below_5 = sum(avg_cpm < 5, na.rm = TRUE),
    pct_avg_below_5 = 100 * sum(avg_cpm < 5, na.rm = TRUE) / n(),
    n_avg_zero = sum(avg_cpm == 0, na.rm = TRUE),
    pct_avg_zero = 100 * sum(avg_cpm == 0, na.rm = TRUE) / n(),
    .groups = "drop"
  )

print(isoform_distribution)
cat("\n")

# Check how many significant isoforms would be excluded by a 5 CPM filter
cat("=== IMPACT OF 5 CPM FILTER (requiring >=5 in EITHER condition) ===\n\n")

# For each isoform, check if it has >=5 CPM in at least one sample in either DMSO or Smg1i
isoform_max_cpm <- sig_cpm %>%
  group_by(cell_type, transcript_id) %>%
  summarize(
    max_cpm_dmso = max(cpm[treatment == "DMSO"], na.rm = TRUE),
    max_cpm_smg1i = max(cpm[treatment == "Smg1i"], na.rm = TRUE),
    max_cpm_either = max(cpm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    passes_filter = max_cpm_either >= 5
  )

filter_impact <- isoform_max_cpm %>%
  group_by(cell_type) %>%
  summarize(
    n_significant = n(),
    n_pass_filter = sum(passes_filter),
    n_excluded = sum(!passes_filter),
    pct_excluded = 100 * sum(!passes_filter) / n(),
    .groups = "drop"
  )

print(filter_impact)
cat("\n")

# Save results
write_tsv(distribution_stats, "../tmp/significant_isoform_cpm_distribution_2026.1.3.tsv")
write_tsv(isoform_distribution, "../tmp/significant_isoform_avg_cpm_distribution_2026.1.3.tsv")
write_tsv(filter_impact, "../tmp/significant_isoform_filter_impact_2026.1.3.tsv")

cat("Results saved to tmp/\n")
