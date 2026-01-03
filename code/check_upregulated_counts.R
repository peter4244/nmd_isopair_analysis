#!/usr/bin/env Rscript
library(tidyverse)

# Load DGE results
all_dge <- map_df(list.files("../longread_dge", pattern = "^nmd_dge", full.names = TRUE), ~{
  df <- read_csv(.x, show_col_types = FALSE)
  df$cell_type_file <- str_split(basename(.x), "_")[[1]][3]
  df
})

celltype_map <- c("at2" = "AT", "dd" = "DD", "ddali" = "DD_ALI",
                  "doali" = "DO", "fb" = "FB", "mv" = "MV")
all_dge$cell_type <- celltype_map[all_dge$cell_type_file]
all_dge <- all_dge %>% filter(!is.na(cell_type))

upregulated <- all_dge %>%
  filter(logFC > 0, adj.P.Val < 0.05)

cat("\n=== UPREGULATED TRANSCRIPT COUNTS ===\n\n")
cat("Total upregulated rows (across all cell types):", nrow(upregulated), "\n")
cat("Unique transcripts (may appear in multiple cell types):", n_distinct(upregulated$txid), "\n\n")

cat("By cell type:\n")
by_ct <- upregulated %>%
  group_by(cell_type) %>%
  summarize(
    n_upregulated = n(),
    n_unique_transcripts = n_distinct(txid)
  )
print(by_ct)

cat("\n=== CHECKING BASELINE PROPORTIONS ===\n\n")
# Load baseline proportions
baseline <- read_tsv("../tmp/isoform_interpretation_q3_upregulated_baseline_2026.1.3.tsv",
                     show_col_types = FALSE)

cat("Baseline proportions rows:", nrow(baseline), "\n")
cat("Unique transcripts in baseline:", n_distinct(baseline$transcript_id), "\n\n")

# Check coverage by cell type
coverage <- upregulated %>%
  select(cell_type, transcript_id = txid) %>%
  distinct() %>%
  left_join(
    baseline %>% select(cell_type, transcript_id) %>% distinct() %>% mutate(has_baseline = TRUE),
    by = c("cell_type", "transcript_id")
  ) %>%
  group_by(cell_type) %>%
  summarize(
    n_upregulated = n(),
    n_with_baseline = sum(has_baseline, na.rm = TRUE),
    n_missing_baseline = sum(is.na(has_baseline)),
    pct_with_baseline = 100 * n_with_baseline / n()
  )

cat("Coverage of baseline proportions by cell type:\n")
print(coverage)
