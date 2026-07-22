#!/usr/bin/env Rscript
# Quick summary of splicing profiles by comparison type
#
# Usage:
#   Rscript scripts/dev/summarize_comparisons.R                                    # default (nonNMD_0.95)
#   Rscript scripts/dev/summarize_comparisons.R --threshold-dir nonNMD_0.50        # lenient threshold
library(tidyverse)

# Parse CLI arguments
args <- commandArgs(trailingOnly = TRUE)
threshold_dir <- "nonNMD_0.95"
if ("--threshold-dir" %in% args) {
  td_idx <- which(args == "--threshold-dir")
  if (td_idx < length(args)) threshold_dir <- args[td_idx + 1]
}
base_dir <- file.path("comparisons", threshold_dir)
cat(sprintf("Using threshold directory: %s/\n\n", base_dir))

# Load deduplicated profiles
profiles <- readRDS(file.path(base_dir, "deduplicated", "all_splicing_profiles.rds"))

# Load all per-comparison pairs and tag them
pair_files <- list.files(base_dir, pattern = "pairs[.]tsv$",
                          recursive = TRUE, full.names = TRUE)
pair_files <- pair_files[!grepl("deduplicated", pair_files)]

all_pairs <- map_dfr(pair_files, function(f) {
  p <- read_tsv(f, show_col_types = FALSE)
  if (nrow(p) == 0) return(tibble())
  # Path: .../comparisons/nonNMD_X.XX/C{1-4}/{run}/pairs.tsv
  parts <- str_split(f, "/")[[1]]
  comp_idx <- which(parts %in% c("C1", "C2", "C3", "C4"))
  tibble(comparison = parts[comp_idx], run = parts[comp_idx + 1]) %>%
    bind_cols(p)
})

# Create join key
profiles <- profiles %>%
  mutate(pair_key = paste(dominant_isoform_id, non_dominant_isoform_id, sep = "::"))

all_pairs <- all_pairs %>%
  mutate(pair_key = paste(dominant_isoform_id, comparator_isoform_id, sep = "::"))

# Join profiles to comparison membership
tagged <- all_pairs %>%
  left_join(profiles %>% select(pair_key, n_events, n_a5ss, n_a3ss, n_se,
    n_missing_internal, n_ir, n_ir_diff, n_partial_ir, n_alt_tss, n_alt_tes,
    n_exons_dom, n_exons_non_dom, length_dom, length_non_dom,
    tss_changed, tes_changed, n_differences),
    by = "pair_key")

cat("============================================================\n")
cat("  PROFILE SUMMARY BY COMPARISON TYPE\n")
cat("============================================================\n\n")

for (comp in c("C1", "C2", "C3", "C4")) {
  d <- tagged %>% filter(comparison == comp)
  if (nrow(d) == 0) {
    cat(sprintf("--- %s: 0 pairs ---\n\n", comp))
    next
  }
  cat(sprintf("--- %s: %d pairs (%d genes, %d runs) ---\n", comp,
    nrow(d), n_distinct(d$gene_id), n_distinct(d$run)))
  cat(sprintf("  Dominant exons:      %.1f mean (%.0f-%.0f)\n",
    mean(d$n_exons_dom), min(d$n_exons_dom), max(d$n_exons_dom)))
  cat(sprintf("  Comparator exons:    %.1f mean (%.0f-%.0f)\n",
    mean(d$n_exons_non_dom), min(d$n_exons_non_dom), max(d$n_exons_non_dom)))
  cat(sprintf("  Dominant length:     %.0f bp mean\n", mean(d$length_dom)))
  cat(sprintf("  Comparator length:   %.0f bp mean\n", mean(d$length_non_dom)))
  cat(sprintf("  UE differences:      %.1f mean\n", mean(d$n_differences)))
  cat(sprintf("  Events per pair:     %.1f mean (%.0f-%.0f)\n",
    mean(d$n_events), min(d$n_events), max(d$n_events)))
  cat(sprintf("  TSS changed:         %d/%d (%.0f%%)\n",
    sum(d$tss_changed), nrow(d), 100*mean(d$tss_changed)))
  cat(sprintf("  TES changed:         %d/%d (%.0f%%)\n",
    sum(d$tes_changed), nrow(d), 100*mean(d$tes_changed)))
  cat("  Event breakdown (mean per pair):\n")
  cat(sprintf("    Alt_TSS: %.2f  Alt_TES: %.2f\n", mean(d$n_alt_tss), mean(d$n_alt_tes)))
  cat(sprintf("    SE: %.2f  Missing_Internal: %.2f\n", mean(d$n_se), mean(d$n_missing_internal)))
  cat(sprintf("    A5SS: %.2f  A3SS: %.2f\n", mean(d$n_a5ss), mean(d$n_a3ss)))
  cat(sprintf("    IR: %.2f  IR_diff: %.2f  Partial_IR: %.2f\n",
    mean(d$n_ir), mean(d$n_ir_diff), mean(d$n_partial_ir)))
  cat("\n")
}

cat("============================================================\n")
cat("  TOTAL EVENTS BY COMPARISON\n")
cat("============================================================\n\n")

for (comp in c("C1", "C2", "C3", "C4")) {
  d <- tagged %>% filter(comparison == comp)
  if (nrow(d) == 0) next
  cat(sprintf("  %s (%d pairs): %d total events (%.1f/pair)\n",
    comp, nrow(d), sum(d$n_events), mean(d$n_events)))
  cat(sprintf("    TSS:%d TES:%d SE:%d MissInt:%d A5SS:%d A3SS:%d IR:%d IRdiff:%d PartIR:%d\n",
    sum(d$n_alt_tss), sum(d$n_alt_tes), sum(d$n_se), sum(d$n_missing_internal),
    sum(d$n_a5ss), sum(d$n_a3ss), sum(d$n_ir), sum(d$n_ir_diff), sum(d$n_partial_ir)))
}

cat("\n============================================================\n")
cat("  PER-RUN BREAKDOWN\n")
cat("============================================================\n\n")

run_summary <- tagged %>%
  group_by(comparison, run) %>%
  summarize(
    n_pairs = n(),
    n_genes = n_distinct(gene_id),
    mean_events = round(mean(n_events), 1),
    pct_tss = round(100*mean(tss_changed)),
    pct_tes = round(100*mean(tes_changed)),
    mean_diffs = round(mean(n_differences), 1),
    .groups = "drop"
  ) %>%
  filter(n_pairs > 0)

print(as.data.frame(run_summary), row.names = FALSE)

cat("\n============================================================\n")
cat("  C1 vs C3/C4: NMD vs NON-NMD STRUCTURAL CONTRAST\n")
cat("============================================================\n\n")

c1c2 <- tagged %>% filter(comparison %in% c("C1", "C2"))
c3c4 <- tagged %>% filter(comparison %in% c("C3", "C4"))

if (nrow(c1c2) > 0 && nrow(c3c4) > 0) {
  cat(sprintf("  NMD comparisons (C1+C2):     %d pairs, %.1f events/pair, %.1f UE diffs/pair\n",
    nrow(c1c2), mean(c1c2$n_events), mean(c1c2$n_differences)))
  cat(sprintf("  Non-NMD comparisons (C3+C4): %d pairs, %.1f events/pair, %.1f UE diffs/pair\n",
    nrow(c3c4), mean(c3c4$n_events), mean(c3c4$n_differences)))

  cat("\n  Event rates (mean per pair):\n")
  cat(sprintf("                     C1+C2    C3+C4\n"))
  cat(sprintf("    Alt_TSS:         %.2f     %.2f\n", mean(c1c2$n_alt_tss), mean(c3c4$n_alt_tss)))
  cat(sprintf("    Alt_TES:         %.2f     %.2f\n", mean(c1c2$n_alt_tes), mean(c3c4$n_alt_tes)))
  cat(sprintf("    SE:              %.2f     %.2f\n", mean(c1c2$n_se), mean(c3c4$n_se)))
  cat(sprintf("    Missing_Int:     %.2f     %.2f\n", mean(c1c2$n_missing_internal), mean(c3c4$n_missing_internal)))
  cat(sprintf("    A5SS:            %.2f     %.2f\n", mean(c1c2$n_a5ss), mean(c3c4$n_a5ss)))
  cat(sprintf("    A3SS:            %.2f     %.2f\n", mean(c1c2$n_a3ss), mean(c3c4$n_a3ss)))
  cat(sprintf("    IR:              %.2f     %.2f\n", mean(c1c2$n_ir), mean(c3c4$n_ir)))
  cat(sprintf("    IR_diff:         %.2f     %.2f\n", mean(c1c2$n_ir_diff), mean(c3c4$n_ir_diff)))
  cat(sprintf("    Partial_IR:      %.2f     %.2f\n", mean(c1c2$n_partial_ir), mean(c3c4$n_partial_ir)))
} else {
  cat("  Insufficient data for contrast (need pairs in both C1/C2 and C3/C4).\n")
}
cat("\n")
