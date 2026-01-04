#!/usr/bin/env Rscript
# Test Q3 CPM calculation

library(tidyverse)
library(edgeR)
library(data.table)

# Load DGEList
dge_isoform <- readRDS("../results/rds/dge_isoform_2026.1.3.rds")
cpm_normalized <- cpm(dge_isoform, normalized = TRUE, log = FALSE)

cat("CPM dimensions:", dim(cpm_normalized), "\n")

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
  select(cell_type, transcript_id = txid, logFC)

cat("Upregulated transcripts:", nrow(upregulated), "\n")

# Convert CPM to long format
cpm_long <- as.data.frame(cpm_normalized) %>%
  rownames_to_column("transcript_id") %>%
  pivot_longer(cols = -transcript_id, names_to = "bamid", values_to = "cpm")

cat("CPM long format:", nrow(cpm_long), "rows\n")

# Add sample metadata
cpm_with_meta <- cpm_long %>%
  left_join(dge_isoform$samples %>% select(bamid, sample, treatment, ct),
            by = "bamid") %>%
  rename(cell_type = ct)

cat("After adding metadata:", nrow(cpm_with_meta), "rows\n")
cat("Cell types:", paste(unique(cpm_with_meta$cell_type), collapse = ", "), "\n")
cat("Treatments:", paste(unique(cpm_with_meta$treatment), collapse = ", "), "\n")

# Join with upregulated
cpm_upregulated <- cpm_with_meta %>%
  inner_join(upregulated, by = c("cell_type", "transcript_id"))

cat("After join with upregulated:", nrow(cpm_upregulated), "rows\n")

# Calculate mean CPM
mean_cpm <- cpm_upregulated %>%
  group_by(cell_type, treatment, transcript_id) %>%
  summarize(mean_cpm = mean(cpm), .groups = "drop")

cat("Mean CPM calculated:", nrow(mean_cpm), "rows\n")

# Pivot wider
mean_cpm_wide <- mean_cpm %>%
  pivot_wider(names_from = treatment,
              values_from = mean_cpm,
              names_prefix = "mean_cpm_")

cat("After pivot_wider:", nrow(mean_cpm_wide), "rows\n")
cat("Columns:", paste(colnames(mean_cpm_wide), collapse = ", "), "\n")
