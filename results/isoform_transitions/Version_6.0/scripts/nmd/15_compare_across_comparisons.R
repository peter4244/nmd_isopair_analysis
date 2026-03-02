#!/usr/bin/env Rscript
#
# Script: 15_compare_across_comparisons.R
# Version: 6.0
# Purpose: Cross-comparison statistical framework — do NMD-triggering transitions (C1/C2)
#          differ from baseline splicing variation (C4)?
#
# Design:
#   Phase 1: NMD vs baseline (paired on gene overlap + unpaired), C4 as baseline
#   Phase 2: Random-effects meta-analysis across cell types
#   Phase 3: Sensitivity analyses (C1/C2 concordance, complexity confound, event direction)
#
# Comparison and cell type scope:
#   Comparisons: C1 (dominant NMD), C2 (top-CPM NMD), C4 (next-best non-NMD baseline)
#   C3 excluded: C4 is a strict subset of C3; using both introduces pseudo-replication.
#   Cell types: AT2, DD, DD_ALI, FB, MV (5 cell types + all_samples)
#   DO excluded from downstream: too few NMD pairs (C1: 1, C2: 6 at 0.50 threshold).
#   DO still contributes to all_samples classification in Script 07 (NMD = union across
#   ANY cell type including DO; non-NMD = intersection across ALL cell types including DO).
#   This is accepted because DO's marginal influence on all_samples definitions does not
#   warrant re-running the full upstream pipeline.
#
# Input:
#   - Per-comparison filtered profiles from Script 14:
#     comparisons/nonNMD_{threshold}/{C}/{run}/splicing_choice_profiles.rds
#   - Per-comparison event_regions from Script 12 (via Script 14):
#     comparisons/nonNMD_{threshold}/{C}/{run}/results/event_regions.rds
#
# Output:
#   - {output_dir}/phase1_paired_results.tsv
#   - {output_dir}/phase1_unpaired_results.tsv
#   - {output_dir}/phase1_event_prevalence.tsv
#   - {output_dir}/phase1_regional_bootstrap.tsv
#   - {output_dir}/phase2_meta_analysis.tsv
#   - {output_dir}/phase3_*.tsv (sensitivity analyses)
#   - {output_dir}/figures/*.pdf
#
# Flags:
#   --threshold T          Non-NMD threshold (default: 0.50)
#   --min-overlap N        Minimum gene overlap for paired tests (default: 50)
#   --n-bootstrap N        Bootstrap iterations for regional enrichment (default: 1000)
#   --output-dir dir/      Output directory (default: comparisons/nonNMD_{threshold}/cross_comparison)
#

library(tidyverse)

# Check for metafor; install guidance if missing
if (!requireNamespace("metafor", quietly = TRUE)) {
  stop("Package 'metafor' is required for meta-analysis. Install with: install.packages('metafor')")
}
library(metafor)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 15: Cross-Comparison Statistical Framework            ║\n")
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

min_overlap <- 50
if ("--min-overlap" %in% args) {
  idx <- which(args == "--min-overlap")
  min_overlap <- as.integer(args[idx + 1])
}

n_bootstrap <- 1000
if ("--n-bootstrap" %in% args) {
  idx <- which(args == "--n-bootstrap")
  n_bootstrap <- as.integer(args[idx + 1])
}

# Parse --source (default: oarfish)
source_type <- "oarfish"
if ("--source" %in% args) {
  src_idx <- which(args == "--source")
  if (src_idx < length(args)) {
    source_type <- args[src_idx + 1]
    if (!source_type %in% c("oarfish", "isocall")) {
      stop("--source must be 'oarfish' or 'isocall'")
    }
  }
}

comp_prefix <- if (source_type == "isocall") "isocall_nonNMD_" else "nonNMD_"
comp_dir <- file.path("comparisons", paste0(comp_prefix, threshold))

output_dir <- file.path(comp_dir, "cross_comparison")
if ("--output-dir" %in% args) {
  idx <- which(args == "--output-dir")
  output_dir <- args[idx + 1]
}

cat(sprintf("  Threshold:     %s\n", threshold))
cat(sprintf("  Min overlap:   %d\n", min_overlap))
cat(sprintf("  Bootstrap n:   %d\n", n_bootstrap))
cat(sprintf("  Output dir:    %s\n\n", output_dir))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(output_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ═══════════════════════════════════════════════════════════════════
# SECTION 0: Helper Functions
# ═══════════════════════════════════════════════════════════════════

# Compute Cliff's delta (nonparametric effect size)
cliffs_delta <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  if (nx == 0 || ny == 0) return(NA_real_)
  # Count dominance: proportion of (x_i, y_j) pairs where x > y minus y > x
  mat <- outer(x, y, function(a, b) sign(a - b))
  sum(mat) / (nx * ny)
}

# Classify a profile into pattern type (mirrors Script 13 logic)
classify_profile <- function(profiles_df) {
  profiles_df %>%
    mutate(
      has_terminal = (n_alt_tss > 0) | (n_alt_tes > 0),
      has_boundary = (n_a5ss > 0) | (n_a3ss > 0),
      has_inclusion = (n_se > 0) | (n_missing_internal > 0),
      has_partial_retention = (n_partial_ir > 0),
      has_full_retention = (n_ir > 0) | (n_ir_diff > 0),
      n_categories = has_terminal + has_boundary + has_inclusion +
                     has_partial_retention + has_full_retention,
      profile_type = case_when(
        n_events == 0 ~ "Identical",
        n_categories == 1 & has_terminal ~ "Terminal-only",
        n_categories == 1 & has_boundary ~ "Boundary-only",
        n_categories == 1 & has_inclusion ~ "Inclusion-only",
        n_categories == 1 & has_partial_retention ~ "Partial_Retention-only",
        n_categories == 1 & has_full_retention ~ "Full_Retention-only",
        n_categories >= 2 ~ "Combined",
        TRUE ~ "Unclassified"
      )
    )
}

# Compute event binary indicators (mirrors Script 10 logic)
add_event_binary <- function(profiles_df) {
  profiles_df %>%
    mutate(
      has_alt_tss = as.integer(n_alt_tss > 0),
      has_alt_tes = as.integer(n_alt_tes > 0),
      has_a5ss = as.integer(n_a5ss > 0),
      has_a3ss = as.integer(n_a3ss > 0),
      has_partial_ir = as.integer(n_partial_ir > 0),
      has_ir = as.integer((n_ir + n_ir_diff) > 0),
      has_se = as.integer(n_se > 0),
      has_missing_internal = as.integer(n_missing_internal > 0)
    )
}

event_type_cols <- c("has_alt_tss", "has_alt_tes", "has_a5ss", "has_a3ss",
                     "has_partial_ir", "has_ir", "has_se", "has_missing_internal")

event_pretty_names <- c(
  has_alt_tss = "Alt_TSS", has_alt_tes = "Alt_TES",
  has_a5ss = "A5SS", has_a3ss = "A3SS",
  has_partial_ir = "Partial_IR", has_ir = "IR",
  has_se = "SE", has_missing_internal = "Missing_Internal"
)

# ═══════════════════════════════════════════════════════════════════
# SECTION 1: Discover and Load Per-Comparison Profiles
# ═══════════════════════════════════════════════════════════════════

cat("Discovering per-comparison profile files...\n")

all_comparisons <- c("C1", "C2", "C4")
all_runs <- c("all_samples", "at", "dd", "dd_ali", "fb", "mv")  # exclude do (too few NMD pairs)

# Build inventory of available data
inventory <- expand_grid(comparison = all_comparisons, run = all_runs) %>%
  mutate(
    profiles_path = file.path(comp_dir, comparison, run, "splicing_choice_profiles.rds"),
    event_regions_path = file.path(comp_dir, comparison, run, "results", "event_regions.rds"),
    profiles_exist = file.exists(profiles_path),
    event_regions_exist = file.exists(event_regions_path)
  )

cat(sprintf("  Profile files found: %d / %d\n",
            sum(inventory$profiles_exist), nrow(inventory)))

# Load all available profiles
all_profiles <- list()
all_event_regions <- list()

for (i in seq_len(nrow(inventory))) {
  if (!inventory$profiles_exist[i]) next

  comp <- inventory$comparison[i]
  run_name <- inventory$run[i]
  key <- paste(comp, run_name, sep = "/")

  prof <- readRDS(inventory$profiles_path[i])
  prof <- prof %>%
    mutate(comparison = comp, run = run_name) %>%
    classify_profile() %>%
    add_event_binary()

  all_profiles[[key]] <- prof

  # Load event_regions if available
  if (inventory$event_regions_exist[i]) {
    er <- readRDS(inventory$event_regions_path[i])
    er$comparison <- comp
    er$run <- run_name
    all_event_regions[[key]] <- er
  }
}

cat(sprintf("  Loaded %d comparison x run profile sets\n", length(all_profiles)))
cat(sprintf("  Loaded %d comparison x run event_regions sets\n\n", length(all_event_regions)))

# Print profile counts
profile_counts <- map_dfr(names(all_profiles), function(key) {
  parts <- str_split(key, "/")[[1]]
  tibble(comparison = parts[1], run = parts[2], n_profiles = nrow(all_profiles[[key]]))
})

cat("Profile counts per comparison x run:\n")
print(pivot_wider(profile_counts, names_from = comparison, values_from = n_profiles), n = Inf)
cat("\n")

# ═══════════════════════════════════════════════════════════════════
# PHASE 1: NMD vs Baseline (Paired + Unpaired)
# ═══════════════════════════════════════════════════════════════════

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   PHASE 1: NMD vs Baseline                                   ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("Using C4 as primary baseline (1:1 pairing per gene)\n\n")

nmd_comparisons <- c("C1", "C2")
baseline_comparison <- "C4"

paired_results <- list()
unpaired_results <- list()
event_prevalence_results <- list()
regional_bootstrap_results <- list()

for (nmd_comp in nmd_comparisons) {
  for (run_name in all_runs) {

    key_nmd <- paste(nmd_comp, run_name, sep = "/")
    key_base <- paste(baseline_comparison, run_name, sep = "/")

    if (!key_nmd %in% names(all_profiles) || !key_base %in% names(all_profiles)) next

    nmd_prof <- all_profiles[[key_nmd]]
    base_prof <- all_profiles[[key_base]]

    cat(sprintf("─── %s vs %s / %s (NMD n=%d, baseline n=%d) ───\n",
                nmd_comp, baseline_comparison, run_name, nrow(nmd_prof), nrow(base_prof)))

    # ─── Find gene overlap for paired analysis ───
    # Both C1/C2 and C4 have 1 pair per gene (dominant_isoform_id is shared)
    overlap_genes <- intersect(
      nmd_prof %>% distinct(gene_id, dominant_isoform_id),
      base_prof %>% distinct(gene_id, dominant_isoform_id)
    )

    n_overlap <- nrow(overlap_genes)
    cat(sprintf("  Gene overlap: %d\n", n_overlap))

    paired_feasible <- n_overlap >= min_overlap

    # ─── UNPAIRED TESTS ───

    # Test 2.1: Event complexity (Wilcoxon rank-sum)
    wrs <- wilcox.test(nmd_prof$n_events, base_prof$n_events)
    cd <- cliffs_delta(nmd_prof$n_events, base_prof$n_events)

    # Test 2.2: Profile type distribution (Chi-square / Fisher's)
    type_tbl <- bind_rows(
      nmd_prof %>% mutate(group = "NMD"),
      base_prof %>% mutate(group = "Baseline")
    ) %>%
      count(group, profile_type) %>%
      pivot_wider(names_from = profile_type, values_from = n, values_fill = 0) %>%
      column_to_rownames("group") %>%
      as.matrix()

    # Remove zero-sum columns
    type_tbl <- type_tbl[, colSums(type_tbl) > 0, drop = FALSE]

    if (ncol(type_tbl) >= 2) {
      expected_counts <- chisq.test(type_tbl)$expected
      if (min(expected_counts) < 5) {
        chi_result <- chisq.test(type_tbl, simulate.p.value = TRUE, B = 10000)
      } else {
        chi_result <- chisq.test(type_tbl)
      }
      n_total <- sum(type_tbl)
      k <- min(nrow(type_tbl), ncol(type_tbl))
      cramers_v <- as.numeric(sqrt(chi_result$statistic / (n_total * (k - 1))))
    } else {
      chi_result <- list(statistic = NA, p.value = NA)
      cramers_v <- NA_real_
    }

    unpaired_results[[length(unpaired_results) + 1]] <- tibble(
      nmd_comparison = nmd_comp, baseline = baseline_comparison, run = run_name,
      n_nmd = nrow(nmd_prof), n_baseline = nrow(base_prof), n_overlap = n_overlap,
      # Event complexity
      nmd_median_events = median(nmd_prof$n_events),
      baseline_median_events = median(base_prof$n_events),
      wilcox_p = wrs$p.value,
      cliffs_delta = cd,
      # Profile type
      chisq_stat = as.numeric(chi_result$statistic),
      chisq_p = chi_result$p.value,
      cramers_v = cramers_v
    )

    # Test 2.3: Event type prevalence (unpaired two-proportion z-test)
    for (et in event_type_cols) {
      nmd_count <- sum(nmd_prof[[et]])
      base_count <- sum(base_prof[[et]])
      nmd_n <- nrow(nmd_prof)
      base_n <- nrow(base_prof)

      # Two-proportion z-test + Haldane-corrected OR (add 0.5 when any cell is zero)
      if (nmd_count + base_count > 0) {
        prop_test <- prop.test(c(nmd_count, base_count), c(nmd_n, base_n), correct = FALSE)
        # Use Haldane correction (+0.5 to all cells) for OR when any cell is zero
        if (nmd_count == 0 || base_count == 0 || nmd_count == nmd_n || base_count == base_n) {
          odds_ratio <- ((nmd_count + 0.5) * (base_n - base_count + 0.5)) /
                        ((base_count + 0.5) * (nmd_n - nmd_count + 0.5))
        } else {
          odds_ratio <- (nmd_count * (base_n - base_count)) /
                        (base_count * (nmd_n - nmd_count))
        }
      } else {
        prop_test <- list(p.value = NA)
        odds_ratio <- NA_real_
      }

      event_prevalence_results[[length(event_prevalence_results) + 1]] <- tibble(
        nmd_comparison = nmd_comp, baseline = baseline_comparison, run = run_name,
        analysis_type = "unpaired",
        event_type = event_pretty_names[et],
        nmd_prop = nmd_count / nmd_n,
        baseline_prop = base_count / base_n,
        diff = nmd_count / nmd_n - base_count / base_n,
        odds_ratio = odds_ratio,
        p_value = prop_test$p.value,
        n_nmd = nmd_n, n_baseline = base_n
      )
    }

    # ─── PAIRED TESTS (on gene overlap set) ───

    if (paired_feasible) {
      cat(sprintf("  Running paired analysis (n = %d genes)...\n", n_overlap))

      # Match NMD and baseline profiles for the same gene
      nmd_paired <- nmd_prof %>%
        inner_join(overlap_genes, by = c("gene_id", "dominant_isoform_id"))
      base_paired <- base_prof %>%
        inner_join(overlap_genes, by = c("gene_id", "dominant_isoform_id"))

      # Ensure 1:1 matching by gene_id + dominant_isoform_id
      paired_data <- nmd_paired %>%
        select(gene_id, dominant_isoform_id,
               nmd_n_events = n_events, nmd_profile_type = profile_type,
               all_of(setNames(event_type_cols, paste0("nmd_", event_type_cols)))) %>%
        inner_join(
          base_paired %>%
            select(gene_id, dominant_isoform_id,
                   base_n_events = n_events, base_profile_type = profile_type,
                   all_of(setNames(event_type_cols, paste0("base_", event_type_cols)))),
          by = c("gene_id", "dominant_isoform_id")
        )

      n_paired <- nrow(paired_data)

      # Paired Test 2.1: Wilcoxon signed-rank on event counts
      paired_wilcox <- wilcox.test(paired_data$nmd_n_events, paired_data$base_n_events,
                                    paired = TRUE)

      # Sign test: proportion of genes where NMD pair has MORE events
      diff_events <- paired_data$nmd_n_events - paired_data$base_n_events
      n_more <- sum(diff_events > 0)
      n_fewer <- sum(diff_events < 0)
      n_tied <- sum(diff_events == 0)
      if (n_more + n_fewer > 0) {
        sign_test_p <- binom.test(n_more, n_more + n_fewer, p = 0.5)$p.value
      } else {
        sign_test_p <- NA_real_
      }

      # Paired Test 2.2: Profile type transitions
      # Count genes that changed profile type
      n_same_type <- sum(paired_data$nmd_profile_type == paired_data$base_profile_type)
      n_diff_type <- sum(paired_data$nmd_profile_type != paired_data$base_profile_type)
      pct_type_change <- 100 * n_diff_type / n_paired

      # Paired Test 2.3: McNemar's test per event type
      for (et in event_type_cols) {
        nmd_col <- paste0("nmd_", et)
        base_col <- paste0("base_", et)
        a <- sum(paired_data[[nmd_col]] == 1 & paired_data[[base_col]] == 1)  # both present
        b <- sum(paired_data[[nmd_col]] == 1 & paired_data[[base_col]] == 0)  # NMD only
        c_val <- sum(paired_data[[nmd_col]] == 0 & paired_data[[base_col]] == 1)  # Baseline only
        d <- sum(paired_data[[nmd_col]] == 0 & paired_data[[base_col]] == 0)  # neither

        # McNemar's test (on discordant pairs b vs c)
        if (b + c_val >= 10) {
          mcnemar_p <- mcnemar.test(matrix(c(a, b, c_val, d), nrow = 2))$p.value
        } else if (b + c_val > 0) {
          # Exact binomial for small discordant counts
          mcnemar_p <- binom.test(b, b + c_val, p = 0.5)$p.value
        } else {
          mcnemar_p <- NA_real_
        }

        event_prevalence_results[[length(event_prevalence_results) + 1]] <- tibble(
          nmd_comparison = nmd_comp, baseline = baseline_comparison, run = run_name,
          analysis_type = "paired",
          event_type = event_pretty_names[et],
          nmd_prop = mean(paired_data[[nmd_col]]),
          baseline_prop = mean(paired_data[[base_col]]),
          diff = mean(paired_data[[nmd_col]]) - mean(paired_data[[base_col]]),
          odds_ratio = NA_real_,  # not meaningful for paired
          p_value = mcnemar_p,
          n_nmd = n_paired, n_baseline = n_paired
        )
      }

      paired_results[[length(paired_results) + 1]] <- tibble(
        nmd_comparison = nmd_comp, baseline = baseline_comparison, run = run_name,
        n_paired = n_paired,
        # Event complexity
        median_diff_events = median(diff_events),
        paired_wilcox_p = paired_wilcox$p.value,
        n_nmd_more = n_more, n_nmd_fewer = n_fewer, n_tied = n_tied,
        sign_test_p = sign_test_p,
        pct_nmd_more = 100 * n_more / max(n_more + n_fewer, 1),
        # Profile type transitions
        n_same_type = n_same_type, n_diff_type = n_diff_type,
        pct_type_change = pct_type_change
      )

    } else {
      cat(sprintf("  Paired analysis skipped (overlap %d < %d)\n", n_overlap, min_overlap))
    }

    # ─── Regional Enrichment Bootstrap (Test 2.4) ───

    key_nmd_er <- paste(nmd_comp, run_name, sep = "/")
    key_base_er <- paste(baseline_comparison, run_name, sep = "/")

    if (key_nmd_er %in% names(all_event_regions) && key_base_er %in% names(all_event_regions)) {

      nmd_er <- all_event_regions[[key_nmd_er]]
      base_er <- all_event_regions[[key_base_er]]

      # Function to compute enrichment ratio from event_regions
      compute_enrichment <- function(er_df) {
        # Overall region proportions (expected)
        total_by_region <- er_df %>%
          count(region_type, name = "total") %>%
          mutate(expected_prop = total / sum(total))

        # Per event type x region
        er_df %>%
          count(event_type, region_type, name = "observed") %>%
          group_by(event_type) %>%
          mutate(event_total = sum(observed), observed_prop = observed / event_total) %>%
          ungroup() %>%
          left_join(total_by_region %>% select(region_type, expected_prop), by = "region_type") %>%
          mutate(enrichment = ifelse(expected_prop > 0, observed_prop / expected_prop, NA_real_))
      }

      nmd_enrich <- compute_enrichment(nmd_er)
      base_enrich <- compute_enrichment(base_er)

      # Bootstrap: resample profiles, recompute enrichment
      # Get profile keys for resampling
      nmd_profile_keys <- nmd_er %>%
        distinct(gene_id, dominant_isoform_id, non_dominant_isoform_id)
      base_profile_keys <- base_er %>%
        distinct(gene_id, dominant_isoform_id, non_dominant_isoform_id)

      cat(sprintf("  Bootstrap regional enrichment (%d iterations)...\n", n_bootstrap))

      set.seed(42)
      boot_diffs <- list()

      for (b in seq_len(n_bootstrap)) {
        # Resample profiles with replacement
        nmd_boot_keys <- nmd_profile_keys[sample(nrow(nmd_profile_keys), replace = TRUE), ]
        base_boot_keys <- base_profile_keys[sample(nrow(base_profile_keys), replace = TRUE), ]

        # Add row ID for joining (handles duplicates from replacement)
        nmd_boot_keys$boot_id <- seq_len(nrow(nmd_boot_keys))
        base_boot_keys$boot_id <- seq_len(nrow(base_boot_keys))

        nmd_boot_er <- nmd_boot_keys %>%
          inner_join(nmd_er, by = c("gene_id", "dominant_isoform_id", "non_dominant_isoform_id"),
                     relationship = "many-to-many")
        base_boot_er <- base_boot_keys %>%
          inner_join(base_er, by = c("gene_id", "dominant_isoform_id", "non_dominant_isoform_id"),
                     relationship = "many-to-many")

        nmd_b_enrich <- compute_enrichment(nmd_boot_er)
        base_b_enrich <- compute_enrichment(base_boot_er)

        # Compute differences for key event x region combos
        diff_b <- nmd_b_enrich %>%
          select(event_type, region_type, nmd_enrichment = enrichment) %>%
          inner_join(
            base_b_enrich %>% select(event_type, region_type, base_enrichment = enrichment),
            by = c("event_type", "region_type")
          ) %>%
          mutate(diff = nmd_enrichment - base_enrichment, boot_iter = b)

        boot_diffs[[b]] <- diff_b
      }

      boot_diffs_df <- bind_rows(boot_diffs)

      # Summarize: 95% CI for enrichment difference
      boot_summary <- boot_diffs_df %>%
        group_by(event_type, region_type) %>%
        summarise(
          mean_diff = mean(diff, na.rm = TRUE),
          ci_lower = quantile(diff, 0.025, na.rm = TRUE),
          ci_upper = quantile(diff, 0.975, na.rm = TRUE),
          boot_n = sum(!is.na(diff)),
          significant = !(ci_lower <= 0 & ci_upper >= 0),  # CI excludes zero
          .groups = "drop"
        ) %>%
        mutate(nmd_comparison = nmd_comp, baseline = baseline_comparison, run = run_name)

      regional_bootstrap_results[[length(regional_bootstrap_results) + 1]] <- boot_summary

    } else {
      cat("  Skipping regional bootstrap (event_regions not available)\n")
    }

    cat("\n")
  }
}

# ─── Compile and save Phase 1 results ───

paired_df <- bind_rows(paired_results)
unpaired_df <- bind_rows(unpaired_results)
regional_bootstrap_df <- bind_rows(regional_bootstrap_results)

# Apply FDR within event prevalence test families
if (length(event_prevalence_results) > 0) {
  event_prevalence_df <- bind_rows(event_prevalence_results) %>%
    group_by(nmd_comparison, run, analysis_type) %>%
    mutate(fdr_p = p.adjust(p_value, method = "fdr")) %>%
    ungroup()
} else {
  event_prevalence_df <- tibble(
    nmd_comparison = character(), baseline = character(), run = character(),
    analysis_type = character(), event_type = character(),
    nmd_prop = numeric(), baseline_prop = numeric(), diff = numeric(),
    odds_ratio = numeric(), p_value = numeric(), n_nmd = integer(),
    n_baseline = integer(), fdr_p = numeric()
  )
}

cat("\n═══ Phase 1 Summary ═══\n\n")

if (nrow(unpaired_df) > 0) {
  cat("Unpaired results:\n")
  print(unpaired_df %>% select(nmd_comparison, run, n_nmd, n_baseline,
                                nmd_median_events, baseline_median_events,
                                cliffs_delta, wilcox_p, cramers_v), n = Inf)
} else {
  cat("No unpaired results (no matching profile files found)\n")
}

if (nrow(paired_df) > 0) {
  cat("\nPaired results:\n")
  print(paired_df %>% select(nmd_comparison, run, n_paired, median_diff_events,
                              paired_wilcox_p, pct_nmd_more, sign_test_p, pct_type_change),
        n = Inf)
}

# Save
write_tsv(unpaired_df, file.path(output_dir, "phase1_unpaired_results.tsv"))
write_tsv(event_prevalence_df, file.path(output_dir, "phase1_event_prevalence.tsv"))

if (nrow(paired_df) > 0) {
  write_tsv(paired_df, file.path(output_dir, "phase1_paired_results.tsv"))
}
if (nrow(regional_bootstrap_df) > 0) {
  write_tsv(regional_bootstrap_df, file.path(output_dir, "phase1_regional_bootstrap.tsv"))
}

cat(sprintf("\n  Phase 1 files saved to %s/\n\n", output_dir))


# ═══════════════════════════════════════════════════════════════════
# PHASE 2: Random-Effects Meta-Analysis Across Cell Types
# ═══════════════════════════════════════════════════════════════════

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   PHASE 2: Random-Effects Meta-Analysis                      ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Meta-analyze key effect sizes across cell types (exclude all_samples — it's not independent)
ct_runs <- c("at", "dd", "dd_ali", "fb", "mv")  # exclude do (too few pairs)

meta_results <- list()

for (nmd_comp in nmd_comparisons) {

  cat(sprintf("─── Meta-analysis for %s vs %s ───\n", nmd_comp, baseline_comparison))

  # ─── Meta 3.1: Cliff's delta for event counts ───
  ct_effects <- unpaired_df %>%
    filter(nmd_comparison == nmd_comp, baseline == baseline_comparison, run %in% ct_runs)

  if (nrow(ct_effects) >= 3) {
    # Approximate SE for Cliff's delta: SE ≈ sqrt((1 - delta^2) / (n_eff - 1))
    # where n_eff = harmonic mean of n_nmd and n_baseline
    ct_effects <- ct_effects %>%
      mutate(
        n_eff = 2 * (n_nmd * n_baseline) / (n_nmd + n_baseline),
        se_delta = sqrt(pmax(1 - cliffs_delta^2, 0.01) / pmax(n_eff - 1, 1))
      )

    # Random-effects model
    tryCatch({
      re_model <- rma(yi = cliffs_delta, sei = se_delta, data = ct_effects, method = "REML")

      meta_results[[length(meta_results) + 1]] <- tibble(
        nmd_comparison = nmd_comp,
        metric = "cliffs_delta_events",
        pooled_estimate = as.numeric(re_model$beta),
        pooled_se = re_model$se,
        pooled_ci_lower = re_model$ci.lb,
        pooled_ci_upper = re_model$ci.ub,
        pooled_p = re_model$pval,
        tau2 = re_model$tau2,
        I2 = re_model$I2,
        Q = re_model$QE,
        Q_p = re_model$QEp,
        k = re_model$k
      )

      cat(sprintf("  Event complexity: pooled delta = %.3f (%.3f, %.3f), p = %.4f\n",
                  re_model$beta, re_model$ci.lb, re_model$ci.ub, re_model$pval))
      cat(sprintf("    I-squared = %.1f%%, tau-squared = %.4f\n", re_model$I2, re_model$tau2))

      # Forest plot
      pdf(file.path(fig_dir, sprintf("forest_events_%s.pdf", nmd_comp)), width = 10, height = 6)
      forest(re_model, slab = ct_effects$run,
             xlab = "Cliff's delta (event count: NMD - baseline)",
             main = sprintf("%s vs %s: Event Complexity", nmd_comp, baseline_comparison))
      dev.off()
    }, error = function(e) {
      cat(sprintf("  Meta-analysis failed: %s\n", e$message))
    })
  } else {
    cat("  Insufficient cell types for meta-analysis\n")
  }

  # ─── Meta 3.2: Cramer's V for profile type shift ───
  if (nrow(ct_effects) >= 3) {
    ct_profile <- ct_effects %>%
      filter(!is.na(cramers_v)) %>%
      mutate(
        # Meta-analyze Cramer's V on the raw scale using delta-method SE.
        # V = sqrt(chi2 / (N*(k-1))). By the delta method:
        # SE(V) ≈ 1 / (2*V) * SE(chi2) / (N*(k-1))
        # where SE(chi2) ≈ sqrt(2*df) for large N. However, for observed V
        # near zero this explodes, so we use a conservative floor.
        # Practical approximation: SE(V) ≈ sqrt(2*(k-1) / N) / (2*sqrt(k-1))
        #   = 1 / sqrt(2*N), which is independent of k for k=2 (our case: 2×m table).
        # We use the more general form with observed chi-square for better accuracy.
        N_total = n_nmd + n_baseline,
        se_v = pmax(sqrt(2 / N_total), 0.01)  # conservative floor to avoid degenerate weights
      )

    if (nrow(ct_profile) >= 3) {
      tryCatch({
        re_profile <- rma(yi = cramers_v, sei = se_v, data = ct_profile, method = "REML")

        meta_results[[length(meta_results) + 1]] <- tibble(
          nmd_comparison = nmd_comp,
          metric = "cramers_v_profile_type",
          pooled_estimate = as.numeric(re_profile$beta),
          pooled_se = re_profile$se,
          pooled_ci_lower = re_profile$ci.lb,
          pooled_ci_upper = re_profile$ci.ub,
          pooled_p = re_profile$pval,
          tau2 = re_profile$tau2,
          I2 = re_profile$I2,
          Q = re_profile$QE,
          Q_p = re_profile$QEp,
          k = re_profile$k
        )

        cat(sprintf("  Profile type shift: pooled V = %.3f, I2 = %.1f%%\n",
                    re_profile$beta, re_profile$I2))
      }, error = function(e) {
        cat(sprintf("  Profile type meta-analysis failed: %s\n", e$message))
      })
    }
  }

  # ─── Meta 3.3: Event-type-specific ORs ───
  for (et in names(event_pretty_names)) {
    et_name <- event_pretty_names[et]

    ct_et <- event_prevalence_df %>%
      filter(nmd_comparison == nmd_comp, analysis_type == "unpaired",
             event_type == et_name, run %in% ct_runs,
             !is.na(odds_ratio), odds_ratio > 0)

    if (nrow(ct_et) >= 3) {
      ct_et <- ct_et %>%
        mutate(
          log_or = log(odds_ratio),
          se_log_or = sqrt(1 / (nmd_prop * n_nmd + 0.5) + 1 / ((1 - nmd_prop) * n_nmd + 0.5) +
                          1 / (baseline_prop * n_baseline + 0.5) + 1 / ((1 - baseline_prop) * n_baseline + 0.5))
        )

      tryCatch({
        re_et <- rma(yi = log_or, sei = se_log_or, data = ct_et, method = "REML")

        meta_results[[length(meta_results) + 1]] <- tibble(
          nmd_comparison = nmd_comp,
          metric = paste0("OR_", et_name),
          pooled_estimate = exp(as.numeric(re_et$beta)),
          pooled_se = re_et$se,
          pooled_ci_lower = exp(re_et$ci.lb),
          pooled_ci_upper = exp(re_et$ci.ub),
          pooled_p = re_et$pval,
          tau2 = re_et$tau2,
          I2 = re_et$I2,
          Q = re_et$QE,
          Q_p = re_et$QEp,
          k = re_et$k
        )
      }, error = function(e) {
        # Silently skip — small n or zero-cell issues
      })
    }
  }

  cat("\n")
}

meta_df <- bind_rows(meta_results)

if (nrow(meta_df) > 0) {
  cat("Meta-analysis summary:\n")
  print(meta_df %>% select(nmd_comparison, metric, pooled_estimate,
                            pooled_ci_lower, pooled_ci_upper, pooled_p, I2, k),
        n = Inf, width = Inf)

  write_tsv(meta_df, file.path(output_dir, "phase2_meta_analysis.tsv"))
  cat(sprintf("\n  Saved: %s\n\n", file.path(output_dir, "phase2_meta_analysis.tsv")))

  # ─── Forest plot: event type ORs ───
  for (nmd_comp in nmd_comparisons) {
    or_meta <- meta_df %>%
      filter(nmd_comparison == nmd_comp, str_starts(metric, "OR_"))

    if (nrow(or_meta) > 0) {
      or_meta <- or_meta %>%
        mutate(event = str_remove(metric, "^OR_"))

      pdf(file.path(fig_dir, sprintf("forest_event_ORs_%s.pdf", nmd_comp)), width = 10, height = 8)
      par(mar = c(5, 10, 4, 2))
      with(or_meta, {
        plot(pooled_estimate, seq_along(event), xlim = range(c(pooled_ci_lower, pooled_ci_upper)),
             pch = 16, xlab = "Odds Ratio (NMD vs Baseline)", ylab = "", yaxt = "n",
             main = sprintf("%s vs %s: Event Type Prevalence ORs", nmd_comp, baseline_comparison),
             log = "x")
        segments(pooled_ci_lower, seq_along(event), pooled_ci_upper, seq_along(event))
        axis(2, at = seq_along(event), labels = event, las = 1)
        abline(v = 1, lty = 2, col = "grey50")
      })
      dev.off()
    }
  }
}


# ═══════════════════════════════════════════════════════════════════
# PHASE 3: Sensitivity Analyses
# ═══════════════════════════════════════════════════════════════════

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   PHASE 3: Sensitivity Analyses                              ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# ─── 3.1: C1 vs C2 concordance ───
cat("3.1: C1 vs C2 concordance...\n\n")

c1c2_concordance <- list()

for (run_name in all_runs) {
  # Get unpaired results for both C1 and C2
  c1_res <- unpaired_df %>% filter(nmd_comparison == "C1", run == run_name)
  c2_res <- unpaired_df %>% filter(nmd_comparison == "C2", run == run_name)

  if (nrow(c1_res) == 0 || nrow(c2_res) == 0) next

  # Get event prevalence differences for both
  c1_prev <- event_prevalence_df %>%
    filter(nmd_comparison == "C1", run == run_name, analysis_type == "unpaired")
  c2_prev <- event_prevalence_df %>%
    filter(nmd_comparison == "C2", run == run_name, analysis_type == "unpaired")

  if (nrow(c1_prev) > 0 && nrow(c2_prev) > 0) {
    matched_prev <- c1_prev %>%
      select(event_type, c1_diff = diff) %>%
      inner_join(c2_prev %>% select(event_type, c2_diff = diff), by = "event_type")

    if (nrow(matched_prev) >= 3) {
      cor_test <- cor.test(matched_prev$c1_diff, matched_prev$c2_diff, method = "pearson")

      c1c2_concordance[[length(c1c2_concordance) + 1]] <- tibble(
        run = run_name,
        n_events_compared = nrow(matched_prev),
        pearson_r = cor_test$estimate,
        cor_p = cor_test$p.value,
        c1_cliffs_delta = c1_res$cliffs_delta,
        c2_cliffs_delta = c2_res$cliffs_delta,
        delta_same_direction = sign(c1_res$cliffs_delta) == sign(c2_res$cliffs_delta)
      )
    }
  }
}

c1c2_df <- bind_rows(c1c2_concordance)
if (nrow(c1c2_df) > 0) {
  cat("C1 vs C2 concordance:\n")
  print(c1c2_df, n = Inf, width = Inf)
  cat("\n")
}

# ─── 3.2: Complexity confound check ───
cat("3.2: Complexity confound check...\n\n")

complexity_confound <- list()

for (nmd_comp in nmd_comparisons) {
  for (run_name in all_runs) {

    key_nmd <- paste(nmd_comp, run_name, sep = "/")
    key_base <- paste(baseline_comparison, run_name, sep = "/")

    if (!key_nmd %in% names(all_profiles) || !key_base %in% names(all_profiles)) next

    nmd_prof <- all_profiles[[key_nmd]]
    base_prof <- all_profiles[[key_base]]

    # Compare structural complexity of the comparator isoform
    for (metric in c("n_exons_non_dom", "length_non_dom", "n_junctions_non_dom")) {
      if (!metric %in% names(nmd_prof) || !metric %in% names(base_prof)) next

      nmd_vals <- nmd_prof[[metric]]
      base_vals <- base_prof[[metric]]

      wrs <- wilcox.test(nmd_vals, base_vals)
      cd <- cliffs_delta(nmd_vals, base_vals)

      complexity_confound[[length(complexity_confound) + 1]] <- tibble(
        nmd_comparison = nmd_comp, run = run_name,
        metric = metric,
        nmd_median = median(nmd_vals, na.rm = TRUE),
        baseline_median = median(base_vals, na.rm = TRUE),
        cliffs_delta = cd,
        wilcox_p = wrs$p.value
      )
    }
  }
}

complexity_confound_df <- bind_rows(complexity_confound)
if (nrow(complexity_confound_df) > 0) {
  cat("Complexity confound (NMD comparator vs baseline comparator):\n")
  print(complexity_confound_df %>%
          filter(metric == "n_exons_non_dom") %>%
          select(nmd_comparison, run, nmd_median, baseline_median, cliffs_delta, wilcox_p),
        n = Inf, width = Inf)
  cat("\n")
}

# ─── 3.3: Event direction (GAIN vs LOSS) ───
cat("3.3: Event direction analysis (GAIN vs LOSS)...\n\n")

direction_results <- list()

for (nmd_comp in nmd_comparisons) {
  for (run_name in all_runs) {

    key_nmd <- paste(nmd_comp, run_name, sep = "/")
    key_base <- paste(baseline_comparison, run_name, sep = "/")

    if (!key_nmd %in% names(all_event_regions) || !key_base %in% names(all_event_regions)) next

    nmd_er <- all_event_regions[[key_nmd]]
    base_er <- all_event_regions[[key_base]]

    # GAIN vs LOSS proportions
    nmd_dir <- nmd_er %>%
      filter(direction %in% c("GAIN", "LOSS")) %>%
      count(direction) %>%
      mutate(prop = n / sum(n))

    base_dir <- base_er %>%
      filter(direction %in% c("GAIN", "LOSS")) %>%
      count(direction) %>%
      mutate(prop = n / sum(n))

    nmd_loss_prop <- nmd_dir %>% filter(direction == "LOSS") %>% pull(prop)
    base_loss_prop <- base_dir %>% filter(direction == "LOSS") %>% pull(prop)
    nmd_loss_n <- nmd_dir %>% filter(direction == "LOSS") %>% pull(n)
    base_loss_n <- base_dir %>% filter(direction == "LOSS") %>% pull(n)
    nmd_total <- sum(nmd_dir$n)
    base_total <- sum(base_dir$n)

    if (length(nmd_loss_prop) > 0 && length(base_loss_prop) > 0 && nmd_total > 0 && base_total > 0) {
      prop_test <- prop.test(c(nmd_loss_n, base_loss_n), c(nmd_total, base_total), correct = FALSE)

      direction_results[[length(direction_results) + 1]] <- tibble(
        nmd_comparison = nmd_comp, run = run_name,
        nmd_loss_prop = nmd_loss_prop,
        baseline_loss_prop = base_loss_prop,
        diff = nmd_loss_prop - base_loss_prop,
        prop_test_p = prop_test$p.value,
        nmd_n_events = nmd_total, baseline_n_events = base_total
      )
    }
  }
}

direction_df <- bind_rows(direction_results)
if (nrow(direction_df) > 0) {
  cat("Event direction (proportion LOSS):\n")
  print(direction_df, n = Inf, width = Inf)
  cat("\n")
}

# ─── Compile and save Phase 3 ───

all_phase3 <- list()
if (nrow(c1c2_df) > 0) all_phase3$c1c2_concordance <- c1c2_df
if (nrow(complexity_confound_df) > 0) all_phase3$complexity_confound <- complexity_confound_df
if (nrow(direction_df) > 0) all_phase3$direction <- direction_df

# Save each as separate TSV
for (name in names(all_phase3)) {
  out_path <- file.path(output_dir, paste0("phase3_", name, ".tsv"))
  write_tsv(all_phase3[[name]], out_path)
  cat(sprintf("  Saved: %s\n", out_path))
}


# ═══════════════════════════════════════════════════════════════════
# SECTION: Summary Figures
# ═══════════════════════════════════════════════════════════════════

cat("\n\nCreating summary figures...\n")

# ─── Figure 1: Profile type distributions side-by-side ───

profile_type_data <- list()
for (key in names(all_profiles)) {
  parts <- str_split(key, "/")[[1]]
  prof <- all_profiles[[key]]
  counts <- prof %>%
    count(profile_type) %>%
    mutate(prop = n / sum(n), comparison = parts[1], run = parts[2])
  profile_type_data[[key]] <- counts
}
profile_type_all <- bind_rows(profile_type_data)

type_colors <- c(
  "Identical" = "#999999", "Terminal-only" = "#E69F00", "Boundary-only" = "#56B4E9",
  "Inclusion-only" = "#009E73", "Partial_Retention-only" = "#F0E442",
  "Full_Retention-only" = "#D55E00", "Combined" = "#CC79A7"
)

# Plot for each run: C1, C2, C4 side-by-side
pdf(file.path(fig_dir, "profile_types_by_comparison.pdf"), width = 14, height = 12)

for (run_name in all_runs) {
  plot_data <- profile_type_all %>%
    filter(run == run_name, comparison %in% c("C1", "C2", "C4"))

  if (nrow(plot_data) == 0) next

  p <- ggplot(plot_data, aes(x = comparison, y = prop, fill = profile_type)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = type_colors) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = sprintf("Profile Type Distribution by Comparison (%s)", run_name),
      subtitle = "C1/C2 = NMD transitions, C4 = baseline non-NMD",
      x = "Comparison", y = "Proportion", fill = "Profile Type"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"))

  print(p)
}

dev.off()
cat(sprintf("  Figure: %s\n", file.path(fig_dir, "profile_types_by_comparison.pdf")))

# ─── Figure 2: Event prevalence comparison ───

pdf(file.path(fig_dir, "event_prevalence_comparison.pdf"), width = 14, height = 10)

for (nmd_comp in nmd_comparisons) {
  prev_data <- event_prevalence_df %>%
    filter(nmd_comparison == nmd_comp, analysis_type == "unpaired",
           run %in% ct_runs) %>%
    mutate(significant = fdr_p < 0.05)

  if (nrow(prev_data) == 0) next

  p <- ggplot(prev_data, aes(x = event_type, y = diff, fill = significant)) +
    geom_col(position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~run, ncol = 3) +
    scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey60")) +
    labs(
      title = sprintf("%s vs %s: Event Type Prevalence Difference", nmd_comp, baseline_comparison),
      subtitle = "Positive = more prevalent in NMD transitions; Red = FDR < 0.05",
      x = "Event Type", y = "Proportion Difference (NMD - Baseline)",
      fill = "Significant"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))

  print(p)
}

dev.off()
cat(sprintf("  Figure: %s\n", file.path(fig_dir, "event_prevalence_comparison.pdf")))

# ─── Figure 3: Paired analysis — gene-level event count differences ───

if (nrow(paired_df) > 0) {
  pdf(file.path(fig_dir, "paired_event_differences.pdf"), width = 12, height = 8)

  p <- ggplot(paired_df, aes(x = run, y = pct_nmd_more, fill = nmd_comparison)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "grey50") +
    geom_text(aes(label = sprintf("n=%d", n_paired)),
              position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Paired Analysis: % of Genes Where NMD Pair Has More Events",
      subtitle = "Dashed line = 50% (no difference); each gene paired on shared dominant isoform",
      x = "Cell Type", y = "% Genes (NMD > Baseline)",
      fill = "NMD Comparison"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"))

  print(p)
  dev.off()
  cat(sprintf("  Figure: %s\n", file.path(fig_dir, "paired_event_differences.pdf")))
}


# ═══════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("FINAL SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat(sprintf("Phase 1: %d unpaired tests, %d paired tests\n",
            nrow(unpaired_df), nrow(paired_df)))
cat(sprintf("         %d event prevalence tests (%d FDR-significant)\n",
            nrow(event_prevalence_df), sum(event_prevalence_df$fdr_p < 0.05, na.rm = TRUE)))
if (nrow(regional_bootstrap_df) > 0) {
  cat(sprintf("         %d regional bootstrap comparisons (%d CI excludes 0)\n",
              nrow(regional_bootstrap_df),
              sum(regional_bootstrap_df$significant, na.rm = TRUE)))
}
cat(sprintf("Phase 2: %d meta-analysis models\n", nrow(meta_df)))
cat(sprintf("Phase 3: C1-C2 concordance (%d), complexity confound (%d), direction (%d)\n",
            nrow(c1c2_df), nrow(complexity_confound_df), nrow(direction_df)))

cat(sprintf("\nAll outputs saved to: %s/\n", output_dir))

cat("\n✓ Step 15 complete\n\n")
