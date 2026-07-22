# Test filtering criteria
library(tidyverse)
library(data.table)

# Load one file to test
dd_dmso <- fread("../tmp/isoform_proportions_dd_dmso_2026.1.3.tsv",
                 skip = 19, data.table = FALSE)

dd_smg1i <- fread("../tmp/isoform_proportions_dd_smg1i_2026.1.3.tsv",
                  skip = 19, data.table = FALSE)

colnames(dd_dmso) <- c("cell_type", "treatment", "sample", "gene_id",
                        "transcript_id", "raw_count", "gene_total",
                        "isoform_proportion", "n_isoforms")
colnames(dd_smg1i) <- colnames(dd_dmso)

# Combine both treatments
dd_all <- bind_rows(dd_dmso, dd_smg1i)

cat("Starting with:", nrow(dd_all), "rows\n")
cat("Unique transcripts:", n_distinct(dd_all$transcript_id), "\n\n")

# Calculate filter metrics per transcript per treatment
filter_stats <- dd_all %>%
  group_by(cell_type, treatment, transcript_id, gene_id) %>%
  summarize(
    total_counts = sum(raw_count),
    n_samples_with_2plus = sum(raw_count >= 2),
    n_samples = n(),
    .groups = "drop"
  )

cat("Filter stats calculated\n\n")

# Apply filter PER TREATMENT
pass_filter <- filter_stats %>%
  filter(total_counts > 10 & n_samples_with_2plus >= 2)

cat("Transcripts passing filter (per treatment):\n")
cat("  DMSO:", sum(pass_filter$treatment == "DMSO"), "\n")
cat("  Smg1i:", sum(pass_filter$treatment == "Smg1i"), "\n")

# Get unique transcripts that pass in EITHER treatment
transcripts_to_keep <- pass_filter %>%
  distinct(transcript_id) %>%
  pull(transcript_id)

cat("\nUnique transcripts passing in EITHER treatment:",
    length(transcripts_to_keep), "\n")

# Filter original data
dd_filtered <- dd_all %>%
  filter(transcript_id %in% transcripts_to_keep)

cat("\nFiltered rows:", nrow(dd_filtered), "\n")
cat("Reduction:", round(100 * (1 - nrow(dd_filtered)/nrow(dd_all)), 1), "%\n")

# Show some examples of what was filtered
filtered_out <- filter_stats %>%
  anti_join(pass_filter, by = c("transcript_id", "treatment")) %>%
  arrange(desc(total_counts))

cat("\nExamples of filtered transcripts:\n")
print(head(filtered_out, 20))
