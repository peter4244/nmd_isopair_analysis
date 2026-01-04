#!/usr/bin/env Rscript
# Verify Q3 CPM fix - check that join works and points are above diagonal

library(tidyverse)
library(edgeR)
library(data.table)

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
  select(cell_type, transcript_id = txid, logFC)

cat("Upregulated transcripts:", nrow(upregulated), "\n")

# Convert CPM to long format with correct transcript IDs
cpm_df <- as.data.frame(cpm_normalized)
cpm_df$transcript_id <- dge_isoform$genes$txid

cpm_long <- cpm_df %>%
  pivot_longer(cols = -transcript_id, names_to = "bamid", values_to = "cpm") %>%
  left_join(dge_isoform$samples %>% select(bamid, sample, treatment, ct),
            by = "bamid") %>%
  rename(cell_type = ct)

cat("CPM long format rows:", nrow(cpm_long), "\n")

# Join with upregulated
cpm_upregulated <- cpm_long %>%
  inner_join(upregulated, by = c("cell_type", "transcript_id"))

cat("After join with upregulated:", nrow(cpm_upregulated), "rows\n")
cat("Unique transcripts:", n_distinct(cpm_upregulated$transcript_id), "\n\n")

# Calculate mean CPM
mean_cpm <- cpm_upregulated %>%
  group_by(cell_type, treatment, transcript_id) %>%
  summarize(mean_cpm = mean(cpm), .groups = "drop") %>%
  pivot_wider(names_from = treatment,
              values_from = mean_cpm,
              names_prefix = "mean_cpm_")

cat("Mean CPM rows:", nrow(mean_cpm), "\n")
cat("Columns:", paste(colnames(mean_cpm), collapse = ", "), "\n\n")

# Check how many points are above vs below diagonal
mean_cpm_filtered <- mean_cpm %>%
  filter(mean_cpm_DMSO > 0, mean_cpm_Smg1i > 0)

above_diagonal <- sum(mean_cpm_filtered$mean_cpm_Smg1i > mean_cpm_filtered$mean_cpm_DMSO)
below_diagonal <- sum(mean_cpm_filtered$mean_cpm_Smg1i < mean_cpm_filtered$mean_cpm_DMSO)
on_diagonal <- sum(mean_cpm_filtered$mean_cpm_Smg1i == mean_cpm_filtered$mean_cpm_DMSO)

cat("Points above diagonal:", above_diagonal, sprintf("(%.1f%%)", 100*above_diagonal/nrow(mean_cpm_filtered)), "\n")
cat("Points below diagonal:", below_diagonal, sprintf("(%.1f%%)", 100*below_diagonal/nrow(mean_cpm_filtered)), "\n")
cat("Points on diagonal:", on_diagonal, sprintf("(%.1f%%)", 100*on_diagonal/nrow(mean_cpm_filtered)), "\n\n")

# Show examples of points below diagonal (if any)
if (below_diagonal > 0) {
  cat("Examples of points below diagonal:\n")
  below_examples <- mean_cpm_filtered %>%
    filter(mean_cpm_Smg1i < mean_cpm_DMSO) %>%
    left_join(upregulated, by = c("cell_type", "transcript_id")) %>%
    arrange(desc(mean_cpm_DMSO / mean_cpm_Smg1i)) %>%
    head(10)

  print(below_examples %>%
    select(cell_type, transcript_id, mean_cpm_DMSO, mean_cpm_Smg1i, logFC) %>%
    mutate(ratio = mean_cpm_DMSO / mean_cpm_Smg1i))
}
