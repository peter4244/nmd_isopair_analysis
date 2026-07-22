#!/usr/bin/env Rscript
# Check diagonal position by cell type

library(tidyverse)
library(edgeR)

# Load DGEList
dge_isoform <- readRDS("../results/rds/dge_isoform_2026.1.3.rds")
cpm_normalized <- cpm(dge_isoform, normalized = TRUE, log = FALSE)

# Load DGE results
all_dge <- map_df(list.files("../longread_dge", pattern = "^nmd_dge", full.names = TRUE), ~{
  df <- read_csv(.x, show_col_types = FALSE)
  celltype_map <- c("at2" = "AT", "dd" = "DD", "ddali" = "DD_ALI",
                    "do" = "DO", "fb" = "FB", "mv" = "MV")
  df$cell_type <- celltype_map[str_extract(basename(.x), "at2|dd|ddali|do|fb|mv")]
  df
})

upregulated <- all_dge %>%
  filter(logFC > 0, adj.P.Val < 0.05) %>%
  select(cell_type, transcript_id = txid, logFC, AveExpr)

# Convert CPM to long format
cpm_df <- as.data.frame(cpm_normalized)
cpm_df$transcript_id <- dge_isoform$genes$txid

cpm_long <- cpm_df %>%
  pivot_longer(cols = -transcript_id, names_to = "bamid", values_to = "cpm") %>%
  left_join(dge_isoform$samples %>% select(bamid, sample, treatment, ct),
            by = "bamid") %>%
  rename(cell_type = ct)

# Calculate mean CPM
mean_cpm <- cpm_long %>%
  inner_join(upregulated, by = c("cell_type", "transcript_id"),
             relationship = "many-to-many") %>%
  group_by(cell_type, treatment, transcript_id) %>%
  summarize(mean_cpm = mean(cpm), .groups = "drop") %>%
  pivot_wider(names_from = treatment,
              values_from = mean_cpm,
              names_prefix = "mean_cpm_") %>%
  filter(mean_cpm_DMSO > 0, mean_cpm_Smg1i > 0) %>%
  mutate(
    above_diagonal = mean_cpm_Smg1i > mean_cpm_DMSO,
    below_diagonal = mean_cpm_Smg1i < mean_cpm_DMSO
  )

# Summary by cell type
by_celltype <- mean_cpm %>%
  group_by(cell_type) %>%
  summarize(
    n_transcripts = n(),
    pct_above = 100 * sum(above_diagonal) / n(),
    pct_below = 100 * sum(below_diagonal) / n(),
    median_cpm_dmso = median(mean_cpm_DMSO),
    median_cpm_smg1i = median(mean_cpm_Smg1i)
  )

cat("Diagonal position by cell type:\n\n")
print(by_celltype, n = 20)

# Check if it's related to expression level
cat("\n\nRelationship to expression level:\n")
mean_cpm_with_dge <- mean_cpm %>%
  left_join(upregulated, by = c("cell_type", "transcript_id"),
            relationship = "many-to-many") %>%
  distinct(cell_type, transcript_id, .keep_all = TRUE)

below_summary <- mean_cpm_with_dge %>%
  group_by(below_diagonal) %>%
  summarize(
    n = n(),
    median_aveexpr = median(AveExpr, na.rm = TRUE),
    median_logfc = median(logFC, na.rm = TRUE),
    median_cpm_dmso = median(mean_cpm_DMSO),
    median_cpm_smg1i = median(mean_cpm_Smg1i)
  )

print(below_summary)
