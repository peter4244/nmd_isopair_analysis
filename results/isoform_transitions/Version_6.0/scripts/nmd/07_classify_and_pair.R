#!/usr/bin/env Rscript

################################################################################
# Script 07: Classify Isoforms and Generate Comparison Pairs
################################################################################
#
# Purpose:
#   Classify isoforms as NMD-sensitive or non-NMD based on DE results, then
#   generate explicit pairs for 4 comparisons x 7 runs (all_samples + 6 cell
#   types). Deduplicates pairs for efficient event detection in core/08.
#
# Input:
#   - data/expression_data_filtered.rds       (tidy CPM: isoform_id, sample_id, cpm, gene_id)
#   - data/sample_metadata.rds                (sample_id, treatment, ct, ...)
#   - data/isoform_structures_filtered.rds    (compact exon structures)
#   - 6 DE CSVs from longread_dge/            (txid, adj.P.Val, logFC, hgnc_id)
#
# Output:
#   - comparisons/nonNMD_{threshold}/{C1,C2,C3,C4}/{run}/classification.rds
#   - comparisons/nonNMD_{threshold}/{C1,C2,C3,C4}/{run}/pairs.tsv
#   - comparisons/nonNMD_{threshold}/deduplicated/all_pairs.tsv
#
# Usage:
#   Rscript scripts/nmd/07_classify_and_pair.R                           # full run (threshold=0.95)
#   Rscript scripts/nmd/07_classify_and_pair.R --non-nmd-threshold 0.50  # lenient threshold
#   Rscript scripts/nmd/07_classify_and_pair.R --test 50                 # limit to 50 genes
#
################################################################################

library(tidyverse)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args
n_genes_limit <- NULL
if (test_mode) {
  test_idx <- which(args == "--test")
  if (test_idx < length(args)) {
    next_arg <- args[test_idx + 1]
    if (!startsWith(next_arg, "--")) {
      n_genes_limit <- as.integer(next_arg)
    }
  }
}

# Parse --non-nmd-threshold (default: 0.95)
non_nmd_threshold <- 0.95
if ("--non-nmd-threshold" %in% args) {
  thr_idx <- which(args == "--non-nmd-threshold")
  if (thr_idx < length(args)) {
    non_nmd_threshold <- as.numeric(args[thr_idx + 1])
    if (is.na(non_nmd_threshold) || non_nmd_threshold < 0 || non_nmd_threshold > 1) {
      stop("--non-nmd-threshold must be a number between 0 and 1")
    }
  }
}

cat("\n")
cat("==================================================================\n")
cat("   STEP 7: Classify Isoforms and Generate Comparison Pairs\n")
cat(sprintf("   Non-NMD threshold: adj.P.Val > %s\n", format(non_nmd_threshold, nsmall = 2)))
cat(sprintf("   Output directory:  comparisons/nonNMD_%s/\n", format(non_nmd_threshold, nsmall = 2)))
cat("==================================================================\n")
cat("\n")

if (test_mode) {
  if (!is.null(n_genes_limit)) {
    cat(sprintf("*** TEST MODE: Limiting to first %d genes ***\n\n", n_genes_limit))
  } else {
    cat("*** TEST MODE ***\n\n")
  }
}

# ==============================================================================
# 0. Configuration
# ==============================================================================

# DE CSV paths — one per cell type
de_dir <- "/Users/petecastaldi/claude_projects/nmd/longread_dge"
de_date <- "2026.1.18"

# Mapping from sample_metadata ct values to DE file name suffixes
ct_to_de_suffix <- c(
  "AT"     = "at2",
  "DD"     = "dd",
  "DD_ALI" = "ddali",
  "DO"     = "doali",
  "FB"     = "fb",
  "MV"     = "mv"
)

# Thresholds
NMD_ADJ_P_THRESHOLD     <- 0.05
NON_NMD_ADJ_P_THRESHOLD <- non_nmd_threshold
MIN_PROPORTION           <- 0.05
DOMINANCE_THRESHOLD      <- 0.50

# Output base — nested under threshold-named subdirectory
output_base <- file.path("comparisons", sprintf("nonNMD_%s", format(NON_NMD_ADJ_P_THRESHOLD, nsmall = 2)))

# ==============================================================================
# 1. Input Validation
# ==============================================================================

validate_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required input file not found: %s", path))
  }
}

cat("Validating inputs...\n")
validate_file("data/expression_data_filtered.rds")
validate_file("data/sample_metadata.rds")
validate_file("data/isoform_structures_filtered.rds")

# Validate DE files
de_files <- setNames(

  file.path(de_dir, paste0("nmd_dge_", ct_to_de_suffix, "_", de_date, ".csv")),
  names(ct_to_de_suffix)
)
for (ct in names(de_files)) {
  validate_file(de_files[[ct]])
}
cat("  All input files found.\n\n")

# ==============================================================================
# 2. Load Data
# ==============================================================================

cat("Loading data...\n")
expression_data <- readRDS("data/expression_data_filtered.rds")
sample_metadata <- readRDS("data/sample_metadata.rds")
isoform_structures <- readRDS("data/isoform_structures_filtered.rds")

cat(sprintf("  Expression: %d rows (%d isoforms, %d samples)\n",
            nrow(expression_data),
            n_distinct(expression_data$isoform_id),
            n_distinct(expression_data$sample_id)))
cat(sprintf("  Samples: %d (%s)\n",
            nrow(sample_metadata),
            paste(names(table(sample_metadata$treatment)), collapse = "/", sep = "")))
cat(sprintf("  Structures: %d isoforms\n", nrow(isoform_structures)))

# In test mode, limit to first N genes
if (test_mode && !is.null(n_genes_limit)) {
  test_genes <- unique(expression_data$gene_id)[1:n_genes_limit]
  expression_data <- expression_data %>% filter(gene_id %in% test_genes)
  isoform_structures <- isoform_structures %>% filter(gene_id %in% test_genes)
  cat(sprintf("  TEST MODE: Filtered to %d genes (%d expression rows, %d structures)\n\n",
              length(test_genes), nrow(expression_data), nrow(isoform_structures)))
}

# Set of isoforms with structural data (required for event detection)
valid_structure_ids <- unique(isoform_structures$isoform_id)

# Join treatment info onto expression data
expression_data <- expression_data %>%
  left_join(sample_metadata %>% select(sample_id, treatment, ct), by = "sample_id")

if (any(is.na(expression_data$treatment))) {
  n_missing <- sum(is.na(expression_data$treatment))
  stop(sprintf("FATAL: %d expression rows have no matching sample metadata", n_missing))
}

# ==============================================================================
# 3. Load and Merge DE Results
# ==============================================================================

cat("\nLoading DE results...\n")
# Restrict DE to isoforms present in our (possibly subsetted) expression data
expr_isoform_ids <- unique(expression_data$isoform_id)

de_results <- list()
for (ct in names(de_files)) {
  de <- read_csv(de_files[[ct]], show_col_types = FALSE) %>%
    select(txid, adj.P.Val, logFC, hgnc_id) %>%
    rename(isoform_id = txid, gene_id_de = hgnc_id) %>%
    mutate(de_ct = ct) %>%
    filter(isoform_id %in% expr_isoform_ids)
  de_results[[ct]] <- de
  cat(sprintf("  %s: %d isoforms (in expression data)\n", ct, nrow(de)))
}

# Cell types present in data (as.character to avoid factor coercion issues)
cell_types <- sort(unique(as.character(sample_metadata$ct)))
cat(sprintf("\nCell types from metadata: %s\n", paste(cell_types, collapse = ", ")))
cat(sprintf("DE files for: %s\n", paste(names(de_files), collapse = ", ")))

# Verify all cell types have DE files
missing_ct <- setdiff(cell_types, names(de_files))
if (length(missing_ct) > 0) {
  stop(sprintf("No DE file mapping for cell types: %s", paste(missing_ct, collapse = ", ")))
}

# ==============================================================================
# 4. Classify Isoforms (all_samples and per-cell-type)
# ==============================================================================

cat("\nClassifying isoforms...\n")

# --- all_samples classification ---
# NMD-sensitive (all_samples): adj.P.Val < 0.05 AND logFC > 0 in ANY cell type
all_de <- bind_rows(de_results)

nmd_sensitive_all <- all_de %>%
  filter(adj.P.Val < NMD_ADJ_P_THRESHOLD, logFC > 0) %>%
  distinct(isoform_id) %>%
  pull(isoform_id)

# Non-NMD (all_samples): must appear in ALL 6 DE files AND have adj.P.Val > threshold in EVERY one
isoforms_per_de <- all_de %>%
  group_by(isoform_id) %>%
  summarize(
    n_de_files = n_distinct(de_ct),
    min_adj_p = min(adj.P.Val),
    .groups = "drop"
  )

non_nmd_all <- isoforms_per_de %>%
  filter(n_de_files == length(de_files)) %>%
  pull(isoform_id)

# From those present in all files, check each has adj.P.Val > threshold
non_nmd_all <- all_de %>%
  filter(isoform_id %in% non_nmd_all) %>%
  group_by(isoform_id) %>%
  summarize(min_adj_p = min(adj.P.Val), .groups = "drop") %>%
  filter(min_adj_p > NON_NMD_ADJ_P_THRESHOLD) %>%
  pull(isoform_id)

# Ensure no overlap
overlap <- intersect(nmd_sensitive_all, non_nmd_all)
if (length(overlap) > 0) {
  warning(sprintf("%d isoforms classified as both NMD-sensitive and non-NMD (all_samples). Removing from non-NMD.", length(overlap)))
  non_nmd_all <- setdiff(non_nmd_all, overlap)
}

cat(sprintf("  all_samples: %d NMD-sensitive, %d non-NMD\n",
            length(nmd_sensitive_all), length(non_nmd_all)))

# --- Per-cell-type classification ---
nmd_sensitive_ct <- list()
non_nmd_ct <- list()

for (ct in cell_types) {
  de <- de_results[[ct]]

  nmd_ct <- de %>%
    filter(adj.P.Val < NMD_ADJ_P_THRESHOLD, logFC > 0) %>%
    pull(isoform_id) %>%
    unique()

  non_ct <- de %>%
    filter(adj.P.Val > NON_NMD_ADJ_P_THRESHOLD) %>%
    pull(isoform_id) %>%
    unique()

  # Remove overlap
  overlap_ct <- intersect(nmd_ct, non_ct)
  if (length(overlap_ct) > 0) {
    non_ct <- setdiff(non_ct, overlap_ct)
  }

  nmd_sensitive_ct[[ct]] <- nmd_ct
  non_nmd_ct[[ct]] <- non_ct
  cat(sprintf("  %s: %d NMD-sensitive, %d non-NMD\n", ct, length(nmd_ct), length(non_ct)))
}

# ==============================================================================
# 5. Helper Functions
# ==============================================================================

#' Apply 5% expression filter
#' @param expr_df expression data (isoform_id, sample_id, cpm, gene_id)
#' @param sample_ids samples to consider for the filter
#' @return character vector of isoform_ids passing filter
apply_5pct_filter <- function(expr_df, sample_ids) {
  expr_subset <- expr_df %>% filter(sample_id %in% sample_ids)

  # Gene totals per sample
  gene_totals <- expr_subset %>%
    group_by(gene_id, sample_id) %>%
    summarize(gene_cpm = sum(cpm), .groups = "drop")

  # Isoform proportions per sample
  proportions <- expr_subset %>%
    left_join(gene_totals, by = c("gene_id", "sample_id")) %>%
    mutate(proportion = if_else(gene_cpm == 0, 0, cpm / gene_cpm))

  # Keep isoforms where max proportion >= 5% in any sample
  passing <- proportions %>%
    group_by(isoform_id) %>%
    summarize(max_prop = max(proportion), .groups = "drop") %>%
    filter(max_prop >= MIN_PROPORTION) %>%
    pull(isoform_id)

  passing
}

#' Calculate dominance proportions
#' @param expr_df expression data
#' @param sample_ids samples for the dominance calculation
#' @return tibble with isoform_id, gene_id, mean_cpm, dominance_proportion
calc_dominance <- function(expr_df, sample_ids) {
  expr_subset <- expr_df %>% filter(sample_id %in% sample_ids)

  # Gene totals per sample
  gene_totals <- expr_subset %>%
    group_by(gene_id, sample_id) %>%
    summarize(gene_cpm = sum(cpm), .groups = "drop")

  # Per-sample proportions
  proportions <- expr_subset %>%
    left_join(gene_totals, by = c("gene_id", "sample_id")) %>%
    mutate(proportion = if_else(gene_cpm == 0, 0, cpm / gene_cpm))

  # Average across samples
  dominance <- proportions %>%
    group_by(isoform_id, gene_id) %>%
    summarize(
      mean_cpm = mean(cpm),
      dominance_proportion = mean(proportion),
      .groups = "drop"
    )

  dominance
}

#' Generate pairs for one comparison type in one run
#' @param comparison Character: "C1", "C2", "C3", or "C4"
#' @param nmd_ids NMD-sensitive isoform IDs for this run
#' @param non_nmd_ids non-NMD isoform IDs for this run
#' @param expr_df expression data (already filtered to relevant samples)
#' @param all_sample_ids all in-scope sample IDs
#' @param smg1i_ids Smg1i sample IDs
#' @param dmso_ids DMSO sample IDs
#' @param filter_scope "all" for treatment-agnostic, "dmso" for DMSO-only
#' @return list with classification tibble and pairs tibble
generate_pairs <- function(comparison, nmd_ids, non_nmd_ids, expr_df,
                           all_sample_ids, smg1i_ids, dmso_ids, filter_scope) {

  # Step 1: Apply 5% filter
  if (filter_scope == "all") {
    passing_5pct <- apply_5pct_filter(expr_df, all_sample_ids)
  } else {
    passing_5pct <- apply_5pct_filter(expr_df, dmso_ids)
  }

  # Filter to isoforms passing 5% AND having structural data
  nmd_passing <- intersect(nmd_ids, intersect(passing_5pct, valid_structure_ids))
  non_nmd_passing <- intersect(non_nmd_ids, intersect(passing_5pct, valid_structure_ids))

  # Step 2: Calculate dominance for relevant treatment groups
  nmd_dominance <- calc_dominance(expr_df, smg1i_ids) %>%
    filter(isoform_id %in% nmd_passing)

  non_nmd_dominance <- calc_dominance(expr_df, dmso_ids) %>%
    filter(isoform_id %in% non_nmd_passing)

  # Build classification
  classification <- bind_rows(
    nmd_dominance %>%
      mutate(nmd_class = "NMD_sensitive",
             treatment_group = "Smg1i",
             is_dominant = dominance_proportion > DOMINANCE_THRESHOLD),
    non_nmd_dominance %>%
      mutate(nmd_class = "non_NMD",
             treatment_group = "DMSO",
             is_dominant = dominance_proportion > DOMINANCE_THRESHOLD)
  )

  # Step 3: Form pairs based on comparison type
  pairs <- tibble(
    gene_id = character(),
    dominant_isoform_id = character(),
    comparator_isoform_id = character(),
    dominant_class = character(),
    comparator_class = character(),
    dominant_dominance_prop = numeric(),
    comparator_dominance_prop = numeric(),
    dominant_mean_cpm = numeric(),
    comparator_mean_cpm = numeric()
  )

  if (comparison %in% c("C1", "C2")) {
    # NMD vs non-NMD comparisons
    # Dominant = non-NMD dominant in DMSO
    # Comparator = NMD-sensitive (dominant in Smg1i for C1, top by CPM for C2)

    dom_non_nmd <- non_nmd_dominance %>%
      filter(dominance_proportion > DOMINANCE_THRESHOLD)

    if (comparison == "C1") {
      comp_nmd <- nmd_dominance %>%
        filter(dominance_proportion > DOMINANCE_THRESHOLD)
    } else {
      # C2: top NMD-sensitive by CPM (no dominance requirement)
      comp_nmd <- nmd_dominance %>%
        group_by(gene_id) %>%
        slice_max(mean_cpm, n = 1, with_ties = FALSE) %>%
        ungroup()
    }

    # For each gene, match dominant non-NMD with NMD comparator
    shared_genes <- intersect(dom_non_nmd$gene_id, comp_nmd$gene_id)

    if (length(shared_genes) > 0) {
      pairs <- inner_join(
        dom_non_nmd %>%
          filter(gene_id %in% shared_genes) %>%
          group_by(gene_id) %>%
          slice_max(dominance_proportion, n = 1, with_ties = FALSE) %>%
          ungroup() %>%
          select(gene_id,
                 dominant_isoform_id = isoform_id,
                 dominant_dominance_prop = dominance_proportion,
                 dominant_mean_cpm = mean_cpm),
        comp_nmd %>%
          filter(gene_id %in% shared_genes) %>%
          group_by(gene_id) %>%
          {if (comparison == "C1")
            slice_max(., dominance_proportion, n = 1, with_ties = FALSE)
           else
            slice_max(., mean_cpm, n = 1, with_ties = FALSE)} %>%
          ungroup() %>%
          select(gene_id,
                 comparator_isoform_id = isoform_id,
                 comparator_dominance_prop = dominance_proportion,
                 comparator_mean_cpm = mean_cpm),
        by = "gene_id"
      ) %>%
        filter(dominant_isoform_id != comparator_isoform_id) %>%
        mutate(dominant_class = "non_NMD", comparator_class = "NMD_sensitive")
    }

  } else if (comparison == "C3") {
    # Non-NMD dominant vs all other non-NMD
    dom_non_nmd <- non_nmd_dominance %>%
      filter(dominance_proportion > DOMINANCE_THRESHOLD)

    # For each gene with a dominant non-NMD, pair with every other non-NMD isoform
    pairs_list <- list()
    for (i in seq_len(nrow(dom_non_nmd))) {
      g <- dom_non_nmd$gene_id[i]
      d <- dom_non_nmd$isoform_id[i]

      others <- non_nmd_dominance %>%
        filter(gene_id == g, isoform_id != d)

      if (nrow(others) > 0) {
        for (j in seq_len(nrow(others))) {
          pairs_list[[length(pairs_list) + 1]] <- tibble(
            gene_id = g,
            dominant_isoform_id = d,
            comparator_isoform_id = others$isoform_id[j],
            dominant_class = "non_NMD",
            comparator_class = "non_NMD",
            dominant_dominance_prop = dom_non_nmd$dominance_proportion[i],
            comparator_dominance_prop = others$dominance_proportion[j],
            dominant_mean_cpm = dom_non_nmd$mean_cpm[i],
            comparator_mean_cpm = others$mean_cpm[j]
          )
        }
      }
    }
    if (length(pairs_list) > 0) {
      pairs <- bind_rows(pairs_list)
    }

  } else if (comparison == "C4") {
    # Non-NMD dominant vs next-best non-NMD by CPM
    dom_non_nmd <- non_nmd_dominance %>%
      filter(dominance_proportion > DOMINANCE_THRESHOLD)

    pairs_list <- list()
    for (i in seq_len(nrow(dom_non_nmd))) {
      g <- dom_non_nmd$gene_id[i]
      d <- dom_non_nmd$isoform_id[i]

      next_best <- non_nmd_dominance %>%
        filter(gene_id == g, isoform_id != d) %>%
        slice_max(mean_cpm, n = 1, with_ties = FALSE)

      if (nrow(next_best) == 1) {
        pairs_list[[length(pairs_list) + 1]] <- tibble(
          gene_id = g,
          dominant_isoform_id = d,
          comparator_isoform_id = next_best$isoform_id,
          dominant_class = "non_NMD",
          comparator_class = "non_NMD",
          dominant_dominance_prop = dom_non_nmd$dominance_proportion[i],
          comparator_dominance_prop = next_best$dominance_proportion,
          dominant_mean_cpm = dom_non_nmd$mean_cpm[i],
          comparator_mean_cpm = next_best$mean_cpm
        )
      }
    }
    if (length(pairs_list) > 0) {
      pairs <- bind_rows(pairs_list)
    }
  }

  list(classification = classification, pairs = pairs)
}

# ==============================================================================
# 6. Generate All 28 Analysis Sets
# ==============================================================================

cat("\nGenerating pairs for 4 comparisons x 7 runs...\n\n")

comparisons <- c("C1", "C2", "C3", "C4")
runs <- c("all_samples", cell_types)

# Track pair counts
pair_counts <- tibble(
  comparison = character(),
  run = character(),
  n_pairs = integer(),
  n_genes = integer()
)

for (comp in comparisons) {
  for (run in runs) {
    cat(sprintf("  %s / %s: ", comp, run))

    # Determine sample scope
    if (run == "all_samples") {
      run_samples <- sample_metadata$sample_id
      nmd_ids <- nmd_sensitive_all
      non_nmd_ids <- non_nmd_all
    } else {
      run_samples <- sample_metadata %>%
        filter(ct == run) %>%
        pull(sample_id)
      nmd_ids <- nmd_sensitive_ct[[run]]
      non_nmd_ids <- non_nmd_ct[[run]]
    }

    smg1i_samples <- sample_metadata %>%
      filter(sample_id %in% run_samples, treatment == "Smg1i") %>%
      pull(sample_id)
    dmso_samples <- sample_metadata %>%
      filter(sample_id %in% run_samples, treatment == "DMSO") %>%
      pull(sample_id)

    # Filter expression to in-scope samples
    expr_run <- expression_data %>% filter(sample_id %in% run_samples)

    # Cell-type-specific re-filtering: only keep isoforms with non-negligible
    # expression in this cell type (max CPM > 0 across run samples)
    if (run != "all_samples") {
      active_isoforms <- expr_run %>%
        group_by(isoform_id) %>%
        summarize(max_cpm = max(cpm), .groups = "drop") %>%
        filter(max_cpm > 0) %>%
        pull(isoform_id)
      nmd_ids <- intersect(nmd_ids, active_isoforms)
      non_nmd_ids <- intersect(non_nmd_ids, active_isoforms)
    }

    # Determine 5% filter scope
    filter_scope <- if (comp %in% c("C1", "C2")) "all" else "dmso"

    # Generate pairs
    result <- generate_pairs(
      comparison = comp,
      nmd_ids = nmd_ids,
      non_nmd_ids = non_nmd_ids,
      expr_df = expr_run,
      all_sample_ids = run_samples,
      smg1i_ids = smg1i_samples,
      dmso_ids = dmso_samples,
      filter_scope = filter_scope
    )

    # Save outputs
    run_label <- tolower(run)
    out_dir <- file.path(output_base, comp, run_label)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    saveRDS(result$classification, file.path(out_dir, "classification.rds"))

    if (nrow(result$pairs) > 0) {
      write_tsv(result$pairs, file.path(out_dir, "pairs.tsv"))
    } else {
      # Write header-only file so downstream knows the schema
      write_tsv(result$pairs, file.path(out_dir, "pairs.tsv"))
    }

    n_p <- nrow(result$pairs)
    n_g <- n_distinct(result$pairs$gene_id)
    cat(sprintf("%d pairs (%d genes)\n", n_p, n_g))

    pair_counts <- bind_rows(pair_counts, tibble(
      comparison = comp, run = run, n_pairs = n_p, n_genes = n_g
    ))
  }
}

# ==============================================================================
# 7. Summary
# ==============================================================================

cat("\n")
cat("==================================================================\n")
cat("  PAIR COUNT SUMMARY\n")
cat("==================================================================\n")

pair_counts_wide <- pair_counts %>%
  mutate(label = sprintf("%d (%d genes)", n_pairs, n_genes)) %>%
  select(comparison, run, label) %>%
  pivot_wider(names_from = comparison, values_from = label)

print(as.data.frame(pair_counts_wide), row.names = FALSE)

cat(sprintf("\nTotal pairs across all 28 sets: %d\n", sum(pair_counts$n_pairs)))

# ==============================================================================
# 8. Deduplication
# ==============================================================================

cat("\nDeduplicating pairs...\n")

# Collect all pairs across all 28 sets
all_pair_files <- list.files(output_base, pattern = "pairs\\.tsv$",
                              recursive = TRUE, full.names = TRUE)

all_pairs <- map_dfr(all_pair_files, function(f) {
  p <- read_tsv(f, show_col_types = FALSE)
  if (nrow(p) == 0) return(tibble())
  # Extract comparison and run from path
  # Path: .../comparisons/nonNMD_X.XX/C{1-4}/{run}/pairs.tsv
  parts <- str_split(f, "/")[[1]]
  comp_idx <- which(parts %in% c("C1", "C2", "C3", "C4"))
  comp <- parts[comp_idx]
  run  <- parts[comp_idx + 1]
  p %>% mutate(source_comparison = comp, source_run = run)
})

if (nrow(all_pairs) > 0) {
  # Unique pairs by (dominant_isoform_id, comparator_isoform_id)
  dedup_pairs <- all_pairs %>%
    distinct(gene_id, dominant_isoform_id, comparator_isoform_id, .keep_all = TRUE) %>%
    select(gene_id, dominant_isoform_id, comparator_isoform_id)

  dedup_dir <- file.path(output_base, "deduplicated")
  dir.create(dedup_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(dedup_pairs, file.path(dedup_dir, "all_pairs.tsv"))

  cat(sprintf("  Total pairs across 28 sets: %d\n", nrow(all_pairs)))
  cat(sprintf("  Unique (deduplicated) pairs: %d\n", nrow(dedup_pairs)))
  cat(sprintf("  Deduplication ratio: %.1f%%\n",
              100 * (1 - nrow(dedup_pairs) / nrow(all_pairs))))

  # Print mapping: which pairs belong to which comparison x run
  pair_membership <- all_pairs %>%
    group_by(dominant_isoform_id, comparator_isoform_id) %>%
    summarize(
      n_sets = n(),
      sets = paste(paste0(source_comparison, "/", source_run), collapse = ", "),
      .groups = "drop"
    )

  shared_pairs <- pair_membership %>% filter(n_sets > 1)
  cat(sprintf("  Pairs appearing in multiple sets: %d\n", nrow(shared_pairs)))
} else {
  cat("  WARNING: No pairs generated across any comparison x run!\n")
  dedup_dir <- file.path(output_base, "deduplicated")
  dir.create(dedup_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(tibble(gene_id = character(), dominant_isoform_id = character(),
                    comparator_isoform_id = character()),
            file.path(dedup_dir, "all_pairs.tsv"))
}

# ==============================================================================
# 9. Validation
# ==============================================================================

cat("\nValidation checks...\n")

if (nrow(all_pairs) > 0) {
  # Check all isoform IDs exist in structures
  all_iso_ids <- unique(c(dedup_pairs$dominant_isoform_id,
                           dedup_pairs$comparator_isoform_id))
  missing_structs <- setdiff(all_iso_ids, valid_structure_ids)
  if (length(missing_structs) > 0) {
    warning(sprintf("  %d isoform IDs in pairs have no structural data!", length(missing_structs)))
  } else {
    cat(sprintf("  All %d isoform IDs in deduplicated pairs have structural data.\n",
                length(all_iso_ids)))
  }

  # Verify C4 pairs are subset of C3 pairs per run
  for (run in runs) {
    run_label <- tolower(run)
    c3_path <- file.path(output_base, "C3", run_label, "pairs.tsv")
    c4_path <- file.path(output_base, "C4", run_label, "pairs.tsv")
    if (file.exists(c3_path) && file.exists(c4_path)) {
      c3 <- read_tsv(c3_path, show_col_types = FALSE)
      c4 <- read_tsv(c4_path, show_col_types = FALSE)
      if (nrow(c4) > 0) {
        c4_keys <- paste(c4$dominant_isoform_id, c4$comparator_isoform_id)
        c3_keys <- paste(c3$dominant_isoform_id, c3$comparator_isoform_id)
        if (all(c4_keys %in% c3_keys)) {
          cat(sprintf("  C4 subset of C3 for %s: OK (%d/%d)\n",
                      run, nrow(c4), nrow(c3)))
        } else {
          missing_n <- sum(!c4_keys %in% c3_keys)
          warning(sprintf("  C4 NOT subset of C3 for %s: %d pairs missing from C3!",
                          run, missing_n))
        }
      }
    }
  }
}

# ==============================================================================
# 10. Write Run Log
# ==============================================================================

run_log_path <- file.path(output_base, "run_log.txt")
cat(sprintf("\nWriting run log to: %s\n", run_log_path))

log_lines <- c(
  "======================================================================",
  "  RUN LOG — Isoform Classification and Pairing (Script 07)",
  "======================================================================",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Working directory: %s", getwd()),
  "",
  "",
  "----------------------------------------------------------------------",
  "  1. PARAMETERS",
  "----------------------------------------------------------------------",
  "",
  sprintf("Non-NMD adj.P.Val threshold:   > %s", format(NON_NMD_ADJ_P_THRESHOLD, nsmall = 2)),
  sprintf("NMD-sensitive adj.P.Val threshold: < %s", format(NMD_ADJ_P_THRESHOLD, nsmall = 2)),
  sprintf("NMD-sensitive logFC requirement:   > 0"),
  sprintf("Minimum isoform proportion:    >= %s (5%% filter)", format(MIN_PROPORTION, nsmall = 2)),
  sprintf("Dominance threshold:           > %s", format(DOMINANCE_THRESHOLD, nsmall = 2)),
  sprintf("Output directory:              %s/", output_base),
  "",
  "",
  "----------------------------------------------------------------------",
  "  2. CLASSIFICATION DEFINITIONS",
  "----------------------------------------------------------------------",
  "",
  "NMD-sensitive isoform:",
  "  An isoform with adj.P.Val < 0.05 AND logFC > 0 in the Smg1i vs DMSO",
  "  differential expression analysis. These are transcripts upregulated",
  "  when NMD is inhibited, indicating they are normally degraded by the",
  "  NMD pathway.",
  "",
  "Non-NMD isoform:",
  sprintf("  An isoform with adj.P.Val > %s in the Smg1i vs DMSO analysis,",
          format(NON_NMD_ADJ_P_THRESHOLD, nsmall = 2)),
  "  indicating no significant change when NMD is inhibited. These are",
  "  transcripts that are not subject to NMD-mediated degradation.",
  "",
  "all_samples classification:",
  "  - NMD-sensitive: union across cell types (NMD in ANY cell type)",
  "  - Non-NMD: intersection across cell types (non-NMD in ALL cell types,",
  "    must appear in all 6 DE files with adj.P.Val above threshold in each)",
  "",
  "Per-cell-type classification:",
  "  - Each cell type classified independently using its own DE results",
  "  - Overlaps (isoforms meeting both NMD and non-NMD criteria) are",
  "    removed from the non-NMD set",
  "",
  "",
  "----------------------------------------------------------------------",
  "  3. COMPARISON DEFINITIONS",
  "----------------------------------------------------------------------",
  "",
  "C1 — Dominant non-NMD vs dominant NMD-sensitive (1 pair/gene)",
  "  Dominant: non-NMD isoform with highest dominance proportion (>50%)",
  "            in DMSO samples",
  "  Comparator: NMD-sensitive isoform with highest dominance proportion",
  "              (>50%) in Smg1i samples",
  "  5% filter scope: all samples (both treatments)",
  "  Purpose: Compare the main stable transcript to the main NMD target",
  "",
  "C2 — Top non-NMD by CPM vs top NMD-sensitive by CPM (1 pair/gene)",
  "  Dominant: non-NMD isoform with highest dominance proportion (>50%)",
  "            in DMSO samples",
  "  Comparator: NMD-sensitive isoform with highest mean CPM in Smg1i",
  "              samples (no dominance requirement)",
  "  5% filter scope: all samples (both treatments)",
  "  Purpose: Compare the main stable transcript to the most expressed",
  "           NMD target, even if it is not dominant",
  "",
  "C3 — Dominant non-NMD vs all other non-NMD (multiple pairs/gene)",
  "  Dominant: non-NMD isoform with highest dominance proportion (>50%)",
  "            in DMSO samples",
  "  Comparator: every other non-NMD isoform for the same gene",
  "  5% filter scope: DMSO samples only",
  "  Purpose: Characterize structural diversity among stable (non-NMD)",
  "           transcripts within a gene — serves as control for C1/C2",
  "",
  "C4 — Dominant non-NMD vs next-best non-NMD by CPM (1 pair/gene)",
  "  Dominant: non-NMD isoform with highest dominance proportion (>50%)",
  "            in DMSO samples",
  "  Comparator: non-NMD isoform with highest mean CPM after excluding",
  "              the dominant",
  "  5% filter scope: DMSO samples only",
  "  Purpose: Focused control — compare the top two stable transcripts",
  "",
  "Each comparison is run across 7 runs:",
  sprintf("  all_samples, %s", paste(cell_types, collapse = ", ")),
  "",
  "",
  "----------------------------------------------------------------------",
  "  4. INPUT DATA",
  "----------------------------------------------------------------------",
  "",
  sprintf("Expression data:      data/expression_data_filtered.rds"),
  sprintf("  Isoforms: %d", n_distinct(expression_data$isoform_id)),
  sprintf("  Samples:  %d", n_distinct(expression_data$sample_id)),
  sprintf("  Rows:     %d", nrow(expression_data)),
  "",
  sprintf("Sample metadata:      data/sample_metadata.rds"),
  sprintf("  Samples:  %d", nrow(sample_metadata)),
  sprintf("  Cell types: %s", paste(cell_types, collapse = ", ")),
  sprintf("  Treatments: %s", paste(sort(unique(as.character(sample_metadata$treatment))), collapse = ", ")),
  "",
  sprintf("Isoform structures:   data/isoform_structures_filtered.rds"),
  sprintf("  Isoforms: %d", nrow(isoform_structures)),
  "",
  "DE files:",
  paste0("  ", names(de_files), ": ", de_files),
  ""
)

# Classification results
log_lines <- c(log_lines,
  "",
  "----------------------------------------------------------------------",
  "  5. CLASSIFICATION RESULTS",
  "----------------------------------------------------------------------",
  "",
  sprintf("all_samples: %d NMD-sensitive, %d non-NMD",
          length(nmd_sensitive_all), length(non_nmd_all)),
  ""
)
for (ct in cell_types) {
  log_lines <- c(log_lines,
    sprintf("%s: %d NMD-sensitive, %d non-NMD",
            ct, length(nmd_sensitive_ct[[ct]]), length(non_nmd_ct[[ct]])))
}

# Pair counts
log_lines <- c(log_lines,
  "",
  "",
  "----------------------------------------------------------------------",
  "  6. PAIR COUNTS",
  "----------------------------------------------------------------------",
  ""
)
for (i in seq_len(nrow(pair_counts))) {
  log_lines <- c(log_lines,
    sprintf("%-4s / %-12s: %d pairs (%d genes)",
            pair_counts$comparison[i], pair_counts$run[i],
            pair_counts$n_pairs[i], pair_counts$n_genes[i]))
}
log_lines <- c(log_lines,
  "",
  sprintf("Total pairs across all 28 sets: %d", sum(pair_counts$n_pairs))
)

# Deduplication
if (nrow(all_pairs) > 0) {
  log_lines <- c(log_lines,
    "",
    "",
    "----------------------------------------------------------------------",
    "  7. DEDUPLICATION",
    "----------------------------------------------------------------------",
    "",
    sprintf("Total pairs across 28 sets:     %d", nrow(all_pairs)),
    sprintf("Unique (deduplicated) pairs:     %d", nrow(dedup_pairs)),
    sprintf("Deduplication ratio:             %.1f%%",
            100 * (1 - nrow(dedup_pairs) / nrow(all_pairs))),
    sprintf("Pairs in multiple sets:          %d", nrow(shared_pairs)),
    "",
    sprintf("Deduplicated pairs file: %s/deduplicated/all_pairs.tsv", output_base)
  )
} else {
  log_lines <- c(log_lines,
    "",
    "",
    "----------------------------------------------------------------------",
    "  7. DEDUPLICATION",
    "----------------------------------------------------------------------",
    "",
    "WARNING: No pairs generated across any comparison x run."
  )
}

log_lines <- c(log_lines, "", "")

writeLines(log_lines, run_log_path)
cat(sprintf("  Run log written (%d lines).\n", length(log_lines)))

cat("\n")
cat("==================================================================\n")
cat("  Step 7 complete.\n")
cat(sprintf("  Outputs in: %s/\n", output_base))
cat(sprintf("  Deduplicated pairs: %s/deduplicated/all_pairs.tsv\n", output_base))
cat(sprintf("  Run log: %s\n", run_log_path))
cat("==================================================================\n")
cat("\n")
