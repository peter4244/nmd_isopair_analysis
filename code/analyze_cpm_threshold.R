#!/usr/bin/env Rscript

# Analyze CPM distribution in Smg1i and test various thresholds
# Focus on EITHER condition logic (>=threshold in DMSO OR Smg1i)

library(tidyverse)
library(edgeR)

# Load longread DGE results
cat("Loading longread DGE results...\n")
dge_files <- list.files("../longread_dge", pattern = "^nmd_dge_.*\\.csv$", full.names = TRUE)

all_dge <- map_df(dge_files, ~{
  df <- read_csv(.x, show_col_types = FALSE)
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

# Identify significant isoforms
sig_isoforms <- all_dge %>%
  filter(adj.P.Val < 0.05) %>%
  select(cell_type, transcript_id = txid, logFC, adj.P.Val)

cat("Loaded", nrow(sig_isoforms), "significant isoform results\n\n")

# Load DGEList
cat("Loading DGEList...\n")
dge_isoform <- readRDS("../results/rds/dge_isoform_2026.1.3.rds")

# Extract normalized CPM
cpm_normalized <- cpm(dge_isoform, normalized = TRUE, log = FALSE)

# Convert to long format
cpm_df <- as.data.frame(cpm_normalized)
cpm_df$transcript_id <- dge_isoform$genes$txid

cpm_long <- cpm_df %>%
  pivot_longer(cols = -transcript_id, names_to = "bamid", values_to = "cpm") %>%
  left_join(dge_isoform$samples %>% select(bamid, sample, treatment, ct),
            by = "bamid") %>%
  rename(cell_type = ct)

cat("CPM data loaded\n\n")

# ========================================================================
# PART 1: Full CPM distribution in Smg1i samples (all isoforms)
# ========================================================================

cat("=== FULL CPM DISTRIBUTION IN SMG1I (ALL ISOFORMS) ===\n\n")

smg1i_cpm_dist <- cpm_long %>%
  filter(treatment == "Smg1i") %>%
  group_by(cell_type) %>%
  summarize(
    n_samples = n_distinct(sample),
    n_observations = n(),
    n_isoforms = n_distinct(transcript_id),
    mean_cpm = mean(cpm, na.rm = TRUE),
    median_cpm = median(cpm, na.rm = TRUE),
    q25_cpm = quantile(cpm, 0.25, na.rm = TRUE),
    q75_cpm = quantile(cpm, 0.75, na.rm = TRUE),
    q90_cpm = quantile(cpm, 0.90, na.rm = TRUE),
    q95_cpm = quantile(cpm, 0.95, na.rm = TRUE),
    min_cpm = min(cpm, na.rm = TRUE),
    max_cpm = max(cpm, na.rm = TRUE),
    n_zero = sum(cpm == 0, na.rm = TRUE),
    pct_zero = 100 * sum(cpm == 0, na.rm = TRUE) / n(),
    n_below_0.25 = sum(cpm < 0.25, na.rm = TRUE),
    pct_below_0.25 = 100 * sum(cpm < 0.25, na.rm = TRUE) / n(),
    n_below_1 = sum(cpm < 1, na.rm = TRUE),
    pct_below_1 = 100 * sum(cpm < 1, na.rm = TRUE) / n(),
    n_below_5 = sum(cpm < 5, na.rm = TRUE),
    pct_below_5 = 100 * sum(cpm < 5, na.rm = TRUE) / n(),
    .groups = "drop"
  )

print(smg1i_cpm_dist)
cat("\n")

# ========================================================================
# PART 2: Test different thresholds with EITHER condition logic
# ========================================================================

cat("=== THRESHOLD TESTING (EITHER DMSO OR SMG1I) ===\n\n")

# For each isoform in each cell type, get max CPM in DMSO and Smg1i
isoform_max_cpm <- cpm_long %>%
  group_by(cell_type, transcript_id) %>%
  summarize(
    max_cpm_dmso = max(cpm[treatment == "DMSO"], na.rm = TRUE),
    max_cpm_smg1i = max(cpm[treatment == "Smg1i"], na.rm = TRUE),
    max_cpm_either = pmax(max_cpm_dmso, max_cpm_smg1i),
    .groups = "drop"
  )

# Add significance status
isoform_max_cpm <- isoform_max_cpm %>%
  left_join(
    sig_isoforms %>% select(cell_type, transcript_id) %>% distinct() %>% mutate(is_significant = TRUE),
    by = c("cell_type", "transcript_id")
  ) %>%
  mutate(is_significant = replace_na(is_significant, FALSE))

# Test multiple thresholds
thresholds <- c(0.1, 0.25, 0.5, 1, 2, 5)

threshold_results <- map_df(thresholds, function(threshold) {

  isoform_max_cpm %>%
    mutate(passes_filter = max_cpm_either >= threshold) %>%
    group_by(cell_type) %>%
    summarize(
      threshold = threshold,
      # All isoforms
      n_total_isoforms = n(),
      n_pass_filter = sum(passes_filter),
      pct_pass_filter = 100 * sum(passes_filter) / n(),
      # Significant isoforms only
      n_significant = sum(is_significant),
      n_sig_pass = sum(is_significant & passes_filter),
      n_sig_excluded = sum(is_significant & !passes_filter),
      pct_sig_excluded = 100 * sum(is_significant & !passes_filter) / sum(is_significant),
      .groups = "drop"
    )
})

cat("Results for each threshold:\n\n")

for (thresh in thresholds) {
  cat("--- Threshold:", thresh, "CPM ---\n")
  result <- threshold_results %>% filter(threshold == thresh)
  print(result %>% select(-threshold))
  cat("\n")
}

# ========================================================================
# PART 3: Detailed breakdown for 0.25 CPM threshold
# ========================================================================

cat("=== DETAILED ANALYSIS: 0.25 CPM THRESHOLD ===\n\n")

threshold_0.25_detail <- isoform_max_cpm %>%
  mutate(passes_0.25 = max_cpm_either >= 0.25) %>%
  group_by(cell_type, is_significant, passes_0.25) %>%
  summarize(
    n_isoforms = n(),
    mean_max_cpm_dmso = mean(max_cpm_dmso, na.rm = TRUE),
    median_max_cpm_dmso = median(max_cpm_dmso, na.rm = TRUE),
    mean_max_cpm_smg1i = mean(max_cpm_smg1i, na.rm = TRUE),
    median_max_cpm_smg1i = median(max_cpm_smg1i, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(cell_type, desc(is_significant), desc(passes_0.25))

print(threshold_0.25_detail)
cat("\n")

# Which condition drives passing the filter for significant isoforms?
sig_filter_source <- isoform_max_cpm %>%
  filter(is_significant) %>%
  mutate(
    passes_0.25 = max_cpm_either >= 0.25,
    passes_dmso_only = max_cpm_dmso >= 0.25 & max_cpm_smg1i < 0.25,
    passes_smg1i_only = max_cpm_smg1i >= 0.25 & max_cpm_dmso < 0.25,
    passes_both = max_cpm_dmso >= 0.25 & max_cpm_smg1i >= 0.25
  ) %>%
  group_by(cell_type) %>%
  summarize(
    n_significant = n(),
    n_pass = sum(passes_0.25),
    n_excluded = sum(!passes_0.25),
    n_dmso_only = sum(passes_dmso_only),
    n_smg1i_only = sum(passes_smg1i_only),
    n_both = sum(passes_both),
    pct_dmso_only = 100 * sum(passes_dmso_only) / sum(passes_0.25),
    pct_smg1i_only = 100 * sum(passes_smg1i_only) / sum(passes_0.25),
    pct_both = 100 * sum(passes_both) / sum(passes_0.25),
    .groups = "drop"
  )

cat("For significant isoforms passing 0.25 CPM filter, which condition drives it?\n\n")
print(sig_filter_source)
cat("\n")

# ========================================================================
# PART 4: Characteristics of excluded significant isoforms
# ========================================================================

cat("=== CHARACTERISTICS OF EXCLUDED SIGNIFICANT ISOFORMS (0.25 CPM threshold) ===\n\n")

excluded_sig <- isoform_max_cpm %>%
  filter(is_significant, max_cpm_either < 0.25) %>%
  left_join(sig_isoforms, by = c("cell_type", "transcript_id"))

excluded_summary <- excluded_sig %>%
  group_by(cell_type) %>%
  summarize(
    n_excluded = n(),
    mean_logFC = mean(logFC, na.rm = TRUE),
    median_logFC = median(logFC, na.rm = TRUE),
    mean_adj_pval = mean(adj.P.Val, na.rm = TRUE),
    median_adj_pval = median(adj.P.Val, na.rm = TRUE),
    mean_max_cpm_either = mean(max_cpm_either, na.rm = TRUE),
    median_max_cpm_either = median(max_cpm_either, na.rm = TRUE),
    n_upregulated = sum(logFC > 0),
    n_downregulated = sum(logFC < 0),
    .groups = "drop"
  )

print(excluded_summary)
cat("\n")

# Save results
write_tsv(smg1i_cpm_dist, "../tmp/smg1i_cpm_full_distribution_2026.1.3.tsv")
write_tsv(threshold_results, "../tmp/threshold_testing_results_2026.1.3.tsv")
write_tsv(sig_filter_source, "../tmp/sig_isoforms_filter_source_2026.1.3.tsv")
write_tsv(excluded_summary, "../tmp/excluded_sig_isoforms_summary_2026.1.3.tsv")

cat("Results saved to tmp/\n")
