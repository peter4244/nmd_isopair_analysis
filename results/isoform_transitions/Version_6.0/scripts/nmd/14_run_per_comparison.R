#!/usr/bin/env Rscript
#
# Script: 14_run_per_comparison.R
# Version: 6.0
# Purpose: Filter deduplicated profiles to each comparison×run and run Scripts 09-13
#
# Input:
#   - comparisons/nonNMD_0.50/deduplicated/splicing_choice_profiles.rds (39,647 profiles)
#   - comparisons/nonNMD_0.50/{C}/{run}/pairs.tsv (per comparison×run pair lists)
#
# Output (per comparison×run):
#   - comparisons/nonNMD_0.50/{C}/{run}/splicing_choice_profiles.rds
#   - comparisons/nonNMD_0.50/{C}/{run}/results/ (Scripts 09-13 outputs)
#
# Comparisons:
#   C1: Dominant non-NMD (DMSO) vs dominant NMD-sensitive (Smg1i)
#   C2: Top non-NMD by CPM vs top NMD-sensitive by CPM
#   C4: Dominant non-NMD vs next-best non-NMD by CPM (baseline)
#   C3 is excluded: C4 pairs are a strict subset of C3, so C3 adds pseudo-replication
#   without additional statistical value.
#
# Cell type exclusion — DO:
#   DO is excluded from per-comparison downstream runs because it yields too few
#   NMD pairs for meaningful analysis (C1: 1 pair, C2: 6 pairs at 0.50 threshold).
#   However, DO still contributes to the all_samples isoform classification in
#   Script 07: NMD-sensitive isoforms are defined as the union across ANY cell type
#   (so DO-specific NMD isoforms are included), and non-NMD isoforms require the
#   intersection across ALL cell types (so DO must also pass). Removing DO from
#   Script 07 would slightly expand the non-NMD set (one fewer cell type to agree)
#   and slightly shrink the NMD set (lose DO-only NMD isoforms). We accept this
#   because DO's influence on the all_samples definitions is marginal, and the
#   individual cell-type runs (which exclude DO) are the primary unit of analysis.
#
# Flags:
#   --threshold T          Non-NMD threshold (default: 0.50)
#   --comparisons C1,C2    Comma-separated comparisons to run (default: C1,C2,C4)
#   --runs r1,r2           Comma-separated runs to run (default: all except do)
#   --min-profiles N       Minimum profiles to run downstream (default: 50)
#   --skip-existing        Skip comparison×run if results already exist
#   --dry-run              Show what would be run without executing
#

library(tidyverse)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 14: Per-Comparison Downstream Analysis Runner         ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Paths
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"
setwd(base_dir)

# ═══════════════════════════════════════════════════════════════════
# Parse CLI arguments
# ═══════════════════════════════════════════════════════════════════

args <- commandArgs(trailingOnly = TRUE)

threshold <- "0.50"
if ("--threshold" %in% args) {
  idx <- which(args == "--threshold")
  threshold <- args[idx + 1]
}

comparisons <- c("C1", "C2", "C4")
if ("--comparisons" %in% args) {
  idx <- which(args == "--comparisons")
  comparisons <- str_split(args[idx + 1], ",")[[1]]
}

all_runs <- c("all_samples", "at", "dd", "dd_ali", "fb", "mv")  # exclude do (too few NMD pairs)
runs <- all_runs
if ("--runs" %in% args) {
  idx <- which(args == "--runs")
  runs <- str_split(args[idx + 1], ",")[[1]]
}

min_profiles <- 50
if ("--min-profiles" %in% args) {
  idx <- which(args == "--min-profiles")
  min_profiles <- as.integer(args[idx + 1])
}

skip_existing <- "--skip-existing" %in% args
dry_run <- "--dry-run" %in% args

comp_dir <- file.path("comparisons", paste0("nonNMD_", threshold))

cat(sprintf("  Threshold:     %s\n", threshold))
cat(sprintf("  Comparisons:   %s\n", paste(comparisons, collapse = ", ")))
cat(sprintf("  Runs:          %s\n", paste(runs, collapse = ", ")))
cat(sprintf("  Min profiles:  %d\n", min_profiles))
cat(sprintf("  Skip existing: %s\n", skip_existing))
cat(sprintf("  Dry run:       %s\n\n", dry_run))

# ═══════════════════════════════════════════════════════════════════
# SECTION 1: Load Deduplicated Profiles
# ═══════════════════════════════════════════════════════════════════

cat("Loading deduplicated profiles...\n")

profiles_path <- file.path(comp_dir, "deduplicated", "splicing_choice_profiles.rds")
if (!file.exists(profiles_path)) {
  stop(sprintf("ERROR: Profiles not found at '%s'. Run Script 08 first.", profiles_path))
}

profiles <- readRDS(profiles_path)
cat(sprintf("  Loaded %d profiles covering %d genes\n\n",
            nrow(profiles), n_distinct(profiles$gene_id)))

# ═══════════════════════════════════════════════════════════════════
# SECTION 2: Build Run Matrix
# ═══════════════════════════════════════════════════════════════════

cat("Building comparison × run matrix...\n\n")

run_matrix <- expand_grid(comparison = comparisons, run = runs) %>%
  mutate(
    pairs_file = file.path(comp_dir, comparison, run, "pairs.tsv"),
    output_dir = file.path(comp_dir, comparison, run, "results"),
    profiles_out = file.path(comp_dir, comparison, run, "splicing_choice_profiles.rds"),
    pairs_exist = file.exists(pairs_file)
  )

# Check which pairs files exist
missing_pairs <- run_matrix %>% filter(!pairs_exist)
if (nrow(missing_pairs) > 0) {
  cat("  WARNING: Missing pairs files:\n")
  for (i in seq_len(nrow(missing_pairs))) {
    cat(sprintf("    %s/%s\n", missing_pairs$comparison[i], missing_pairs$run[i]))
  }
  cat("\n")
}

run_matrix <- run_matrix %>% filter(pairs_exist)

# ═══════════════════════════════════════════════════════════════════
# SECTION 3: Filter Profiles Per Comparison×Run
# ═══════════════════════════════════════════════════════════════════

cat("Filtering profiles per comparison × run...\n\n")

# Track results
run_summary <- list()

for (i in seq_len(nrow(run_matrix))) {
  comp <- run_matrix$comparison[i]
  run_name <- run_matrix$run[i]
  pairs_file <- run_matrix$pairs_file[i]
  output_dir <- run_matrix$output_dir[i]
  profiles_out <- run_matrix$profiles_out[i]

  cat(sprintf("─── %s / %s ───\n", comp, run_name))

  # Check if results already exist
  if (skip_existing && file.exists(file.path(output_dir, "complexity_relationship_results.rds"))) {
    cat("  Skipping (results already exist)\n\n")
    run_summary[[length(run_summary) + 1]] <- tibble(
      comparison = comp, run = run_name, n_pairs = NA_integer_,
      n_profiles = NA_integer_, status = "skipped_existing"
    )
    next
  }

  # Load pairs
  pairs <- read_tsv(pairs_file, show_col_types = FALSE)
  cat(sprintf("  Pairs in file: %d\n", nrow(pairs)))

  # Inner join: match pairs to profiles
  # pairs.tsv has: gene_id, dominant_isoform_id, comparator_isoform_id
  # profiles has: gene_id, dominant_isoform_id, non_dominant_isoform_id
  filtered <- profiles %>%
    inner_join(
      pairs %>% select(gene_id, dominant_isoform_id, comparator_isoform_id),
      by = c("gene_id", "dominant_isoform_id",
             "non_dominant_isoform_id" = "comparator_isoform_id")
    )

  cat(sprintf("  Matched profiles: %d (%d genes)\n", nrow(filtered), n_distinct(filtered$gene_id)))

  # Check minimum threshold
  if (nrow(filtered) < min_profiles) {
    cat(sprintf("  SKIPPED: n = %d < %d minimum\n\n", nrow(filtered), min_profiles))
    run_summary[[length(run_summary) + 1]] <- tibble(
      comparison = comp, run = run_name, n_pairs = nrow(pairs),
      n_profiles = nrow(filtered), status = "skipped_low_n"
    )
    next
  }

  if (dry_run) {
    cat("  [DRY RUN] Would run Scripts 09-13\n\n")
    run_summary[[length(run_summary) + 1]] <- tibble(
      comparison = comp, run = run_name, n_pairs = nrow(pairs),
      n_profiles = nrow(filtered), status = "dry_run"
    )
    next
  }

  # Save filtered profiles
  dir.create(dirname(profiles_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(filtered, profiles_out)
  cat(sprintf("  Saved: %s\n", profiles_out))

  # ─── Run Script 09 (Complexity) ───
  cat("  Running Script 09 (Complexity)...\n")
  tryCatch({
    # Script 09 expects --input and --output-dir
    # It reads profiles, computes complexity bins, saves profiles_with_bins
    system2("Rscript",
            c("scripts/nmd/09_analyze_complexity_relationship.R",
              "--input", profiles_out,
              "--output-dir", output_dir),
            stdout = FALSE, stderr = FALSE)
    cat("    Done\n")
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", e$message))
  })

  # ─── Run Script 10 (Co-occurrence) ───
  # Depends on Script 09's profiles_with_bins
  cat("  Running Script 10 (Co-occurrence)...\n")
  profiles_with_bins <- file.path(output_dir, "splicing_choice_profiles_with_bins.rds")
  tryCatch({
    system2("Rscript",
            c("scripts/nmd/10_analyze_cooccurrence.R",
              "--input", profiles_with_bins,
              "--output-dir", output_dir),
            stdout = FALSE, stderr = FALSE)
    cat("    Done\n")
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", e$message))
  })

  # ─── Run Scripts 11, 12, 13 (independent of each other) ───
  cat("  Running Script 11 (Spatial)...\n")
  tryCatch({
    system2("Rscript",
            c("scripts/nmd/11_analyze_spatial_patterns.R",
              "--input", profiles_out,
              "--output-dir", output_dir),
            stdout = FALSE, stderr = FALSE)
    cat("    Done\n")
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", e$message))
  })

  cat("  Running Script 12 (Functional)...\n")
  tryCatch({
    system2("Rscript",
            c("scripts/nmd/12_analyze_functional_context.R",
              "--input", profiles_out,
              "--output-dir", output_dir),
            stdout = FALSE, stderr = FALSE)
    cat("    Done\n")
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", e$message))
  })

  cat("  Running Script 13 (Patterns)...\n")
  tryCatch({
    system2("Rscript",
            c("scripts/nmd/13_analyze_patterns.R",
              "--input", profiles_with_bins,
              "--output-dir", output_dir),
            stdout = FALSE, stderr = FALSE)
    cat("    Done\n")
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", e$message))
  })

  cat("\n")
  run_summary[[length(run_summary) + 1]] <- tibble(
    comparison = comp, run = run_name, n_pairs = nrow(pairs),
    n_profiles = nrow(filtered), status = "completed"
  )
}

# ═══════════════════════════════════════════════════════════════════
# SECTION 4: Summary
# ═══════════════════════════════════════════════════════════════════

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

summary_df <- bind_rows(run_summary)

cat("Run results:\n")
print(summary_df, n = Inf)

cat(sprintf("\n  Completed: %d\n", sum(summary_df$status == "completed")))
cat(sprintf("  Skipped (low n): %d\n", sum(summary_df$status == "skipped_low_n")))
cat(sprintf("  Skipped (existing): %d\n", sum(summary_df$status == "skipped_existing")))

# Save run summary
summary_path <- file.path(comp_dir, "per_comparison_run_summary.tsv")
write_tsv(summary_df, summary_path)
cat(sprintf("\n  Summary saved: %s\n", summary_path))

cat("\n✓ Step 14 complete\n\n")
