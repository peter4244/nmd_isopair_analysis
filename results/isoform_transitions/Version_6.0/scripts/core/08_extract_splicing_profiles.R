#!/usr/bin/env Rscript

################################################################################
# Script 08: Extract Splicing Choice Profiles
################################################################################
#
# Purpose:
#   Build splicing choice profiles by comparing non-dominant to dominant isoforms.
#   Each profile characterizes the structural differences between isoforms.
#   Performs comprehensive event detection on FILTERED data from Script 06.
#
# Input:
#   - data/isoform_union_exons_annotated_filtered.rds
#   - data/dominant_isoforms_filtered.rds  (only in auto-generation mode; skipped with --pairs-file)
#   - data/union_exons_filtered.rds
#   - data/isoform_structures_filtered.rds
#
# Output:
#   - data/splicing_choice_profiles.rds (default, or --output path)
#   - data/splicing_choice_profiles_intermediate.rds (early output for development)
#
# Usage:
#   Rscript scripts/core/08_extract_splicing_profiles.R                              # full run
#   Rscript scripts/core/08_extract_splicing_profiles.R --test 50                    # limit to 50 genes
#   Rscript scripts/core/08_extract_splicing_profiles.R --reconstruction_check       # with reconstruction verification
#   Rscript scripts/core/08_extract_splicing_profiles.R --pairs-file pairs.tsv       # explicit pairs
#   Rscript scripts/core/08_extract_splicing_profiles.R --output out/profiles.rds    # custom output path
#   Rscript scripts/core/08_extract_splicing_profiles.R --run-log path/run_log.txt        # append to run log
#   Rscript scripts/core/08_extract_splicing_profiles.R --no-store                        # skip profile store
#   Rscript scripts/core/08_extract_splicing_profiles.R --clear-store                     # clear store, recompute all
#   Rscript scripts/core/08_extract_splicing_profiles.R --test 50 --reconstruction_check  # combined
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

reconstruction_check <- "--reconstruction_check" %in% args

pairs_file <- NULL
if ("--pairs-file" %in% args) {
  pf_idx <- which(args == "--pairs-file")
  if (pf_idx < length(args)) {
    pairs_file <- args[pf_idx + 1]
  } else {
    stop("--pairs-file requires a file path argument")
  }
}

output_path <- NULL  # will be set after data_dir is known
if ("--output" %in% args) {
  out_idx <- which(args == "--output")
  if (out_idx < length(args)) {
    output_path <- args[out_idx + 1]
  } else {
    stop("--output requires a file path argument")
  }
}

run_log_path <- NULL
if ("--run-log" %in% args) {
  rl_idx <- which(args == "--run-log")
  if (rl_idx < length(args)) {
    run_log_path <- args[rl_idx + 1]
  } else {
    stop("--run-log requires a file path argument")
  }
}

use_store <- !("--no-store" %in% args)
clear_store <- "--clear-store" %in% args

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

# Set data directory based on source
data_dir <- if (source_type == "isocall") "data/isocall" else "data"

# Default output path if not specified via --output
if (is.null(output_path)) {
  output_path <- file.path(data_dir, "splicing_choice_profiles.rds")
}

# ==============================================================================
# Input Validation Helpers
# ==============================================================================

validate_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required input file not found: %s\nHave you run the prerequisite scripts?", path))
  }
}

validate_columns <- function(df, required_cols, name) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing required columns: %s", name, paste(missing, collapse = ", ")))
  }
}

validate_overlap <- function(ids_a, ids_b, name_a, name_b, min_pct = 10) {
  overlap <- length(intersect(ids_a, ids_b))
  pct <- 100 * overlap / length(ids_a)
  cat(sprintf("  Cross-check: %d/%d %s found in %s (%.1f%%)\n",
              overlap, length(ids_a), name_a, name_b, pct))
  if (pct < min_pct) {
    stop(sprintf("FATAL: Only %.1f%% of %s found in %s. Files may be from different pipeline runs.",
                 pct, name_a, name_b))
  }
}

# ==============================================================================
# DETECTION THRESHOLDS (adjust these for sensitivity analyses)
# ==============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 8: Extract Splicing Choice Profiles                   ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

if (test_mode) {
  if (!is.null(n_genes_limit)) {
    cat(sprintf("*** TEST MODE: Limiting to first %d genes ***\n\n", n_genes_limit))
  } else {
    cat("*** TEST MODE ***\n\n")
  }
}

if (reconstruction_check) {
  cat("*** RECONSTRUCTION CHECK: On-the-fly verification enabled ***\n\n")
}

if (!is.null(pairs_file)) {
  cat(sprintf("*** PAIRS FILE: %s ***\n\n", pairs_file))
}

if (output_path != file.path(data_dir, "splicing_choice_profiles.rds")) {
  cat(sprintf("*** OUTPUT: %s ***\n\n", output_path))
}

if (!use_store) {
  cat("*** PROFILE STORE: Disabled (--no-store) ***\n\n")
} else if (clear_store) {
  cat("*** PROFILE STORE: Will clear existing store (--clear-store) ***\n\n")
}

# ==============================================================================
# 0. Validate Inputs
# ==============================================================================

cat("Validating inputs...\n")
validate_file(file.path(data_dir, "isoform_union_exons_annotated_filtered.rds"))
if (is.null(pairs_file)) {
  validate_file(file.path(data_dir, "dominant_isoforms_filtered.rds"))
}
validate_file(file.path(data_dir, "union_exons_filtered.rds"))
validate_file(file.path(data_dir, "isoform_structures_filtered.rds"))
validate_file("scripts/core/event_detection_functions.R")
cat("  All required input files found.\n\n")

# ==============================================================================
# 0a. Load Event Detection Functions (Shared Library)
# ==============================================================================

cat("Loading event detection functions...\n")
source("scripts/core/event_detection_functions.R")
cat("  Functions loaded\n")

if (reconstruction_check) {
  cat("Loading reconstruction functions (--reconstruction_check active)...\n")
  validate_file("scripts/core/reconstruction_functions.R")
  source("scripts/core/reconstruction_functions.R")
  cat("  reconstruct_dominant_v2() and verify_transcript() loaded\n")
}
cat("\n")

cat("Detection thresholds:\n")
cat(sprintf("  TSS tolerance: %d bp\n", TSS_TOLERANCE))
cat(sprintf("  TES tolerance: %d bp\n", TES_TOLERANCE))
cat(sprintf("  Splice site threshold: %d bp\n\n", SPLICE_SITE_THRESHOLD))

# ==============================================================================
# 0b. Legacy - Function definitions now in event_detection_functions.R
# ==============================================================================
# All detection functions are now sourced from the shared library above

# ==============================================================================
# 1. Load Data
# ==============================================================================

cat("Loading filtered data...\n")
annotated_mapping <- readRDS(file.path(data_dir, "isoform_union_exons_annotated_filtered.rds"))
dominant_isoforms <- NULL
if (is.null(pairs_file)) {
  dominant_isoforms <- readRDS(file.path(data_dir, "dominant_isoforms_filtered.rds"))
}
union_exons <- readRDS(file.path(data_dir, "union_exons_filtered.rds"))
isoform_structures_compact <- readRDS(file.path(data_dir, "isoform_structures_filtered.rds"))

cat(sprintf("  Annotated mappings: %d\n", nrow(annotated_mapping)))
if (!is.null(dominant_isoforms)) {
  cat(sprintf("  Dominant isoforms: %d\n", nrow(dominant_isoforms)))
} else {
  cat("  Dominant isoforms: (from pairs file)\n")
}
cat(sprintf("  Union exons: %d\n", nrow(union_exons)))
cat(sprintf("  Isoform structures (compact): %d\n", nrow(isoform_structures_compact)))

# Column validation
validate_columns(annotated_mapping, c("gene_id", "isoform_id", "union_exon_id", "region_type"),
                 "annotated_mapping")
if (!is.null(dominant_isoforms)) {
  validate_columns(dominant_isoforms, c("gene_id", "dominant_isoform_id"),
                   "dominant_isoforms")
}
validate_columns(union_exons, c("gene_id", "union_exon_id", "union_exon_start", "union_exon_end", "strand"),
                 "union_exons")
validate_columns(isoform_structures_compact,
                 c("isoform_id", "gene_id", "seqnames", "strand", "exon_starts", "exon_ends"),
                 "isoform_structures")

# Cross-file consistency checks
cat("\nCross-file consistency checks:\n")
validate_overlap(unique(union_exons$gene_id),
                 unique(isoform_structures_compact$gene_id),
                 "union_exons gene_ids", "isoform_structures")
if (!is.null(dominant_isoforms)) {
  validate_overlap(unique(dominant_isoforms$gene_id),
                   unique(union_exons$gene_id),
                   "dominant_isoforms gene_ids", "union_exons")
  validate_overlap(unique(dominant_isoforms$dominant_isoform_id),
                   unique(isoform_structures_compact$isoform_id),
                   "dominant isoform_ids", "isoform_structures")
}
cat("")

# Expand isoform_structures to one row per exon (needed for event detection)
cat("\nExpanding isoform structures to one row per exon...\n")
isoform_structures <- isoform_structures_compact %>%
  rowwise() %>%
  mutate(
    exon_data = list(tibble(
      exon_number = seq_along(exon_starts),
      exon_start = exon_starts,
      exon_end = exon_ends
    ))
  ) %>%
  unnest(exon_data) %>%
  select(isoform_id, gene_id, seqnames, strand, exon_number, exon_start, exon_end) %>%
  ungroup()

cat(sprintf("  Expanded to %d exon rows\n", nrow(isoform_structures)))

# ==============================================================================
# 1b. Profile Store: Load, Validate, Clear
# ==============================================================================
# Splicing profiles are purely structural (depend only on exon coordinates and
# union exons, not on classification thresholds or expression). A central store
# at data/splicing_profile_store.rds accumulates profiles across runs. core/08
# checks the store, only computes new pairs, and appends results.
#
# Store key: dominant_isoform_id + non_dominant_isoform_id
# Note: pairs file uses "comparator_isoform_id" which maps to "non_dominant_isoform_id"
# in the profile output. The store uses the profile column names.

store_path <- file.path(data_dir, "splicing_profile_store.rds")
if (clear_store && file.exists(store_path)) {
  file.remove(store_path)
  cat("\n  Cleared existing profile store.\n")
}

profile_store <- NULL
profile_store_keys <- character(0)
n_store_cached <- 0L
n_store_new <- 0L
n_store_reverified <- 0L
cached_profiles_list <- list()

if (use_store && file.exists(store_path)) {
  cat("\nLoading profile store...\n")
  stored <- readRDS(store_path)

  # Validate fingerprint: cached profiles depend on exon coordinates, union exons,
  # annotated mappings, and detection thresholds. If any change, store is stale.
  current_fp <- list(
    isoform_structures_mtime = file.mtime(file.path(data_dir, "isoform_structures_filtered.rds")),
    union_exons_mtime = file.mtime(file.path(data_dir, "union_exons_filtered.rds")),
    annotated_mapping_mtime = file.mtime(file.path(data_dir, "isoform_union_exons_annotated_filtered.rds")),
    tss_tolerance = TSS_TOLERANCE,
    tes_tolerance = TES_TOLERANCE,
    splice_site_threshold = SPLICE_SITE_THRESHOLD
  )

  if (!identical(stored$fingerprint, current_fp)) {
    cat("  WARNING: Store fingerprint mismatch — inputs have changed.\n")
    for (key in names(current_fp)) {
      if (!identical(stored$fingerprint[[key]], current_fp[[key]])) {
        cat(sprintf("    Changed: %s\n", key))
      }
    }
    cat("  Discarding stale store. All pairs will be recomputed.\n")
  } else {
    profile_store <- stored$profiles
    profile_store_keys <- paste(profile_store$dominant_isoform_id,
                                 profile_store$non_dominant_isoform_id, sep = "::")
    cat(sprintf("  Profile store: %d existing profiles (fingerprint OK)\n",
                nrow(profile_store)))
  }
} else if (use_store) {
  cat("\n  No existing profile store found. All pairs will be computed.\n")
}

# Prepare union exons for reconstruction (column names reconstruction expects)
if (reconstruction_check) {
  union_exons_for_recon <- union_exons %>%
    rename(start = union_exon_start, end = union_exon_end, chr = seqnames)
}

# ==============================================================================
# 2. Identify Non-Dominant Isoforms
# ==============================================================================

isoforms_classified <- NULL
if (is.null(pairs_file)) {
  cat("\nIdentifying non-dominant isoforms...\n")

  # Get list of all isoforms per gene
  all_isoforms <- annotated_mapping %>%
    distinct(gene_id, isoform_id)

  # Mark which are dominant
  isoforms_classified <- all_isoforms %>%
    left_join(
      dominant_isoforms %>% select(gene_id, dominant_isoform_id),
      by = "gene_id"
    ) %>%
    mutate(
      is_dominant = (isoform_id == dominant_isoform_id),
      role = if_else(is_dominant, "dominant", "non_dominant")
    )

  cat(sprintf("  Dominant isoforms: %d\n", sum(isoforms_classified$is_dominant)))
  cat(sprintf("  Non-dominant isoforms: %d\n", sum(!isoforms_classified$is_dominant)))
} else {
  cat("\nSkipping non-dominant identification (using --pairs-file).\n")
}

# ==============================================================================
# 2b. Load Explicit Pairs File (if specified)
# ==============================================================================

explicit_pairs <- NULL
valid_pairs <- NULL

if (!is.null(pairs_file)) {
  cat(sprintf("\nLoading explicit pairs from: %s\n", pairs_file))
  explicit_pairs <- read_tsv(pairs_file, show_col_types = FALSE)
  validate_columns(explicit_pairs, c("gene_id", "dominant_isoform_id"),
                   "pairs_file")

  # One-vs-all mode: expand rows with empty/NA comparator_isoform_id
  if (!"comparator_isoform_id" %in% names(explicit_pairs)) {
    explicit_pairs$comparator_isoform_id <- NA_character_
  }
  one_vs_all_rows <- explicit_pairs %>%
    filter(is.na(comparator_isoform_id) | comparator_isoform_id == "")
  explicit_rows <- explicit_pairs %>%
    filter(!is.na(comparator_isoform_id) & comparator_isoform_id != "")

  if (nrow(one_vs_all_rows) > 0) {
    cat(sprintf("  One-vs-all mode: expanding %d rows...\n", nrow(one_vs_all_rows)))
    # For each one-vs-all row, pair dominant against all other isoforms for that gene
    expanded <- list()
    for (i in seq_len(nrow(one_vs_all_rows))) {
      g <- one_vs_all_rows$gene_id[i]
      d <- one_vs_all_rows$dominant_isoform_id[i]
      other_isos <- isoform_structures_compact %>%
        filter(gene_id == g, isoform_id != d) %>%
        pull(isoform_id)
      if (length(other_isos) > 0) {
        expanded[[length(expanded) + 1]] <- tibble(
          gene_id = g,
          dominant_isoform_id = d,
          comparator_isoform_id = other_isos
        )
      }
    }
    if (length(expanded) > 0) {
      expanded_pairs <- bind_rows(expanded)
      cat(sprintf("  Expanded to %d pairs\n", nrow(expanded_pairs)))
      explicit_pairs <- bind_rows(explicit_rows, expanded_pairs)
    } else {
      explicit_pairs <- explicit_rows
    }
  }

  cat(sprintf("  Total explicit pairs: %d\n", nrow(explicit_pairs)))

  # Validate all isoforms exist in structures
  all_pair_isoforms <- unique(c(explicit_pairs$dominant_isoform_id,
                                 explicit_pairs$comparator_isoform_id))
  missing <- setdiff(all_pair_isoforms, unique(isoform_structures_compact$isoform_id))
  if (length(missing) > 0) {
    warning(sprintf("%d isoforms in pairs file not found in structures: %s",
                    length(missing), paste(head(missing, 5), collapse = ", ")))
  }
}

# ==============================================================================
# 3. Build Comparison Profiles with Event Detection
# ==============================================================================

cat("\nBuilding splicing choice profiles with event detection...\n")
cat("  Comparing each non-dominant to dominant isoform...\n")

# For each gene, compare non-dominant to dominant
splicing_profiles <- list()

# Only process genes that have dominant info AND structure data AND union exons
genes_with_structures <- unique(isoform_structures$gene_id)
genes_with_ue <- unique(union_exons$gene_id)

if (!is.null(pairs_file)) {
  # Explicit pairs mode: filter to valid pairs
  valid_pairs <- explicit_pairs %>%
    filter(gene_id %in% genes_with_structures,
           gene_id %in% genes_with_ue,
           dominant_isoform_id %in% unique(isoform_structures$isoform_id),
           comparator_isoform_id %in% unique(isoform_structures$isoform_id))

  cat(sprintf("  Valid pairs (after filtering): %d/%d\n",
              nrow(valid_pairs), nrow(explicit_pairs)))

  genes_with_dominant <- unique(valid_pairs$gene_id)
} else {
  # Auto-generation mode (default)
  genes_with_dominant <- intersect(
    unique(dominant_isoforms$gene_id),
    intersect(genes_with_structures, genes_with_ue)
  )
}

cat(sprintf("  Genes with both dominant and structure data: %d\n", length(genes_with_dominant)))

# Apply gene limit in test mode
if (test_mode && !is.null(n_genes_limit) && n_genes_limit < length(genes_with_dominant)) {
  genes_with_dominant <- genes_with_dominant[1:n_genes_limit]
  cat(sprintf("  TEST MODE: Limited to first %d genes\n", n_genes_limit))
}

# Process in batches with progress reporting
batch_size <- 1000
n_batches <- ceiling(length(genes_with_dominant) / batch_size)

for (batch_idx in 1:n_batches) {
  start_idx <- (batch_idx - 1) * batch_size + 1
  end_idx <- min(batch_idx * batch_size, length(genes_with_dominant))
  batch_genes <- genes_with_dominant[start_idx:end_idx]

  for (gene in batch_genes) {
    # Get gene strand
    gene_strand <- union_exons %>%
      filter(gene_id == gene) %>%
      pull(strand) %>%
      unique()

    if (length(gene_strand) != 1) next

    if (!is.null(pairs_file)) {
      # Explicit pairs mode: iterate over each dominant for this gene
      gene_pairs <- valid_pairs %>% filter(gene_id == gene)
      dom_isos_for_gene <- unique(gene_pairs$dominant_isoform_id)
    } else {
      # Auto-generation mode (default)
      dom_isos_for_gene <- dominant_isoforms %>%
        filter(gene_id == gene) %>%
        pull(dominant_isoform_id)

      if (length(dom_isos_for_gene) != 1) next
    }

    for (dom_iso in dom_isos_for_gene) {

    if (!is.null(pairs_file)) {
      non_dom_isos <- gene_pairs %>%
        filter(dominant_isoform_id == dom_iso) %>%
        pull(comparator_isoform_id)
    } else {
      non_dom_isos <- isoforms_classified %>%
        filter(gene_id == gene, !is_dominant) %>%
        pull(isoform_id)
    }

    if (length(non_dom_isos) == 0) next

    # Get structure for dominant isoform
    dom_structure <- isoform_structures %>%
      filter(isoform_id == dom_iso) %>%
      arrange(exon_number)

    if (nrow(dom_structure) == 0) next

    # Get union exon usage for dominant
    dom_exons <- annotated_mapping %>%
      filter(isoform_id == dom_iso) %>%
      select(union_exon_id, region_type_dom = region_type)

    # Compare each non-dominant to dominant
    for (non_dom_iso in non_dom_isos) {

      # --- Profile store: check for cached profile ---
      if (use_store && length(profile_store_keys) > 0) {
        pair_key <- paste(dom_iso, non_dom_iso, sep = "::")
        if (pair_key %in% profile_store_keys) {
          cached_idx <- match(pair_key, profile_store_keys)
          cached_row <- profile_store[cached_idx, ]

          # Re-verify if --reconstruction_check active and cached status is NA
          if (reconstruction_check && is.na(cached_row$reconstruction_status)) {
            non_dom_structure_rv <- isoform_structures %>%
              filter(isoform_id == non_dom_iso) %>% arrange(exon_number)

            if (nrow(non_dom_structure_rv) > 0) {
              rv_result <- tryCatch({
                events_rv <- cached_row$detailed_events[[1]]
                comp_for_recon <- non_dom_structure_rv %>%
                  rename(chr = seqnames, transcript_id = isoform_id) %>%
                  select(chr, exon_start, exon_end, strand, gene_id, transcript_id)

                if (nrow(events_rv) == 0) {
                  dom_for_verify <- dom_structure %>%
                    rename(chr = seqnames, transcript_id = isoform_id) %>%
                    select(chr, exon_start, exon_end, strand, gene_id, transcript_id)
                  vr <- verify_transcript(dom_for_verify, comp_for_recon, gene_strand)
                  list(status = if (vr$pass) "PASS" else "FAIL", reason = vr$reason)
                } else {
                  gene_ue <- union_exons_for_recon %>% filter(gene_id == gene)
                  reconstructed <- reconstruct_dominant_v2(comp_for_recon, events_rv, gene_ue)
                  if (nrow(reconstructed) == 0) {
                    list(status = "ERROR", reason = "Reconstruction produced 0 exons")
                  } else {
                    vr <- verify_transcript(dom_structure, reconstructed, gene_strand)
                    list(status = if (vr$pass) "PASS" else "FAIL", reason = vr$reason)
                  }
                }
              }, error = function(e) {
                list(status = "ERROR", reason = paste0("Re-verify exception: ", e$message))
              })

              cached_row$reconstruction_status <- rv_result$status
              cached_row$reconstruction_reason <- rv_result$reason
              n_store_reverified <- n_store_reverified + 1L
            }
          }

          cached_profiles_list[[pair_key]] <- cached_row
          n_store_cached <- n_store_cached + 1L
          next
        }
      }

      # Get structure for non-dominant isoform
      non_dom_structure <- isoform_structures %>%
        filter(isoform_id == non_dom_iso) %>%
        arrange(exon_number)

      if (nrow(non_dom_structure) == 0) next

      # Get union exon usage for non-dominant
      non_dom_exons <- annotated_mapping %>%
        filter(isoform_id == non_dom_iso) %>%
        select(union_exon_id, region_type_non_dom = region_type)

      # Full outer join to find differences
      comparison <- union_exons %>%
        filter(gene_id == gene) %>%
        select(gene_id, union_exon_id, union_exon_number, union_exon_start, union_exon_end) %>%
        left_join(dom_exons, by = "union_exon_id") %>%
        left_join(non_dom_exons, by = "union_exon_id") %>%
        mutate(
          in_dominant = !is.na(region_type_dom),
          in_non_dominant = !is.na(region_type_non_dom),
          exon_status = case_when(
            in_dominant & in_non_dominant ~ "shared",
            in_dominant & !in_non_dominant ~ "dominant_only",
            !in_dominant & in_non_dominant ~ "non_dominant_only",
            TRUE ~ "neither"
          )
        )

      # ===== EVENT DETECTION (validated hierarchical algorithm) =====
      # Uses detect_events_for_pair() which performs hierarchical detection:
      # IR → Boundary (A5SS/A3SS/Partial_IR) → SE → Missing_Internal → Terminal

      # Prepare exon structures for detect_events_for_pair()
      dom_exons_for_detection <- dom_structure %>%
        rename(chr = seqnames, transcript_id = isoform_id)

      comp_exons_for_detection <- non_dom_structure %>%
        rename(chr = seqnames, transcript_id = isoform_id)

      # Run hierarchical event detection
      events <- detect_events_for_pair(
        dom_exons_for_detection, comp_exons_for_detection,
        gene, dom_iso, non_dom_iso, gene_strand
      )

      # On-the-fly reconstruction check
      recon_status <- NA_character_
      recon_reason <- NA_character_

      if (reconstruction_check) {
        recon_result <- tryCatch({
          comp_exons_for_recon <- non_dom_structure %>%
            rename(chr = seqnames, transcript_id = isoform_id) %>%
            select(chr, exon_start, exon_end, strand, gene_id, transcript_id)

          if (nrow(events) == 0) {
            # No events: comparator IS the reconstruction — verify directly
            # Use consistent column names (rename dom_structure to match comp format)
            dom_for_verify <- dom_structure %>%
              rename(chr = seqnames, transcript_id = isoform_id) %>%
              select(chr, exon_start, exon_end, strand, gene_id, transcript_id)
            vr <- verify_transcript(dom_for_verify, comp_exons_for_recon, gene_strand)
            list(status = if (vr$pass) "PASS" else "FAIL", reason = vr$reason)
          } else {
            gene_ue <- union_exons_for_recon %>% filter(gene_id == gene)

            reconstructed <- reconstruct_dominant_v2(comp_exons_for_recon, events, gene_ue)

            if (nrow(reconstructed) == 0) {
              list(status = "ERROR", reason = "Reconstruction produced 0 exons")
            } else {
              vr <- verify_transcript(dom_structure, reconstructed, gene_strand)
              list(status = if (vr$pass) "PASS" else "FAIL", reason = vr$reason)
            }
          }
        }, error = function(e) {
          list(status = "ERROR", reason = paste0("Exception: ", e$message))
        })

        recon_status <- recon_result$status
        recon_reason <- recon_result$reason
      }

      # Aggregate event counts (backward compatible with Scripts 08-12)
      if (nrow(events) > 0) {
        n_a5ss <- sum(events$event_type == "A5SS")
        n_a3ss <- sum(events$event_type == "A3SS")
        n_partial_ir <- sum(events$event_type %in% c("Partial_IR_5", "Partial_IR_3"))
        n_ir <- sum(events$event_type == "IR")
        n_se <- sum(events$event_type == "SE")
        n_missing_internal <- sum(events$event_type == "Missing_Internal")
        n_ir_diff <- sum(grepl("^IR_diff", events$event_type))
        n_alt_tss <- sum(events$event_type == "Alt_TSS")
        n_alt_tes <- sum(events$event_type == "Alt_TES")
        tss_changed <- n_alt_tss > 0
        tes_changed <- n_alt_tes > 0
      } else {
        n_a5ss <- 0L; n_a3ss <- 0L; n_partial_ir <- 0L
        n_ir <- 0L; n_se <- 0L; n_missing_internal <- 0L
        n_ir_diff <- 0L; n_alt_tss <- 0L; n_alt_tes <- 0L
        tss_changed <- FALSE; tes_changed <- FALSE
      }
      n_dual_boundary <- 0L  # No longer used; decomposed into two events

      # Get complexity metrics from compact isoform_structures
      dom_complexity <- isoform_structures_compact %>%
        filter(isoform_id == dom_iso) %>%
        select(n_exons, n_junctions, tx_start, tx_end) %>%
        mutate(length = tx_end - tx_start + 1)

      non_dom_complexity <- isoform_structures_compact %>%
        filter(isoform_id == non_dom_iso) %>%
        select(n_exons, n_junctions, tx_start, tx_end) %>%
        mutate(length = tx_end - tx_start + 1)

      # Calculate summary statistics for this comparison
      profile <- tibble(
        gene_id = gene,
        dominant_isoform_id = dom_iso,
        non_dominant_isoform_id = non_dom_iso,

        # Exon counts
        n_union_exons_total = nrow(comparison),
        n_exons_shared = sum(comparison$exon_status == "shared"),
        n_exons_dominant_only = sum(comparison$exon_status == "dominant_only"),
        n_exons_non_dominant_only = sum(comparison$exon_status == "non_dominant_only"),

        # Structural info
        n_exons_in_dominant = sum(comparison$in_dominant),
        n_exons_in_non_dominant = sum(comparison$in_non_dominant),

        # Complexity metrics for dominant isoform
        n_exons_dom = dom_complexity$n_exons,
        n_junctions_dom = dom_complexity$n_junctions,
        length_dom = dom_complexity$length,

        # Complexity metrics for non-dominant isoform (used by Script 08)
        n_exons_non_dom = non_dom_complexity$n_exons,
        n_junctions_non_dom = non_dom_complexity$n_junctions,
        length_non_dom = non_dom_complexity$length,

        # TSS/TES changes
        tss_changed = tss_changed,
        tes_changed = tes_changed,

        # Splicing events (event-level counts)
        n_a5ss = n_a5ss,
        n_a3ss = n_a3ss,
        n_partial_ir = n_partial_ir,
        n_ir = n_ir,
        n_se = n_se,
        n_missing_internal = n_missing_internal,
        n_ir_diff = n_ir_diff,
        n_alt_tss = n_alt_tss,
        n_alt_tes = n_alt_tes,

        # Detailed events from hierarchical detection
        n_events = nrow(events),
        detailed_events = list(events),

        # Total differences (needed for Script 08)
        n_differences = n_exons_dominant_only + n_exons_non_dominant_only,

        # Store full comparison as nested tibble
        comparison_detail = list(comparison),

        # Reconstruction check results (NA when --reconstruction_check not used)
        reconstruction_status = recon_status,
        reconstruction_reason = recon_reason
      )

      splicing_profiles[[paste0(gene, "_", dom_iso, "_", non_dom_iso)]] <- profile
      n_store_new <- n_store_new + 1L
    }
    } # end for (dom_iso ...)
  }

  # Progress report
  cat(sprintf("  Processed batch %d/%d (%d genes, %.1f%% complete)\n",
              batch_idx, n_batches, length(batch_genes),
              100 * end_idx / length(genes_with_dominant)))

  # Save intermediate results after batch 1 for downstream development
  if (batch_idx == 1 && length(splicing_profiles) > 0) {
    cat("\n  ═══════════════════════════════════════════════════════════\n")
    cat("  CHECKPOINT: Saving intermediate results for development\n")
    cat("  ═══════════════════════════════════════════════════════════\n")
    intermediate_profiles <- bind_rows(splicing_profiles)
    saveRDS(intermediate_profiles, file.path(data_dir, "splicing_choice_profiles_intermediate.rds"))
    cat(sprintf("  ✓ Saved %d profiles to:\n", nrow(intermediate_profiles)))
    cat(sprintf("     %s/splicing_choice_profiles_intermediate.rds\n", data_dir))
    cat(sprintf("  ✓ Covering %d genes\n", n_distinct(intermediate_profiles$gene_id)))
    cat("\n  This file can be used for developing Scripts 08-12\n")
    cat("  while the full processing continues...\n")
    cat("  ═══════════════════════════════════════════════════════════\n\n")
  }
}

# Combine newly computed profiles
new_profiles_df <- bind_rows(splicing_profiles)

# Merge cached profiles from store
cached_profiles_df <- if (length(cached_profiles_list) > 0) {
  bind_rows(cached_profiles_list)
} else {
  tibble()
}

splicing_profiles_df <- bind_rows(new_profiles_df, cached_profiles_df)

if (use_store && (n_store_cached > 0 || n_store_new > 0)) {
  cat(sprintf("\n  From profile store: %d cached, %d newly computed", n_store_cached, n_store_new))
  if (n_store_reverified > 0) {
    cat(sprintf(", %d re-verified", n_store_reverified))
  }
  cat("\n")
}

cat(sprintf("\n  Total splicing choice profiles: %d\n", nrow(splicing_profiles_df)))
cat(sprintf("  Covering %d genes\n", n_distinct(splicing_profiles_df$gene_id)))

# ==============================================================================
# 4. Summary Statistics
# ==============================================================================

cat("\nProfile statistics:\n")

cat(sprintf("  Mean union exons per gene: %.1f\n",
            mean(splicing_profiles_df$n_union_exons_total)))
cat(sprintf("  Mean shared exons: %.1f (%.1f%%)\n",
            mean(splicing_profiles_df$n_exons_shared),
            100 * mean(splicing_profiles_df$n_exons_shared / splicing_profiles_df$n_union_exons_total)))
cat(sprintf("  Mean dominant-only exons: %.1f (%.1f%%)\n",
            mean(splicing_profiles_df$n_exons_dominant_only),
            100 * mean(splicing_profiles_df$n_exons_dominant_only / splicing_profiles_df$n_union_exons_total)))
cat(sprintf("  Mean non-dominant-only exons: %.1f (%.1f%%)\n",
            mean(splicing_profiles_df$n_exons_non_dominant_only),
            100 * mean(splicing_profiles_df$n_exons_non_dominant_only / splicing_profiles_df$n_union_exons_total)))

cat("\nComplexity metrics:\n")
cat(sprintf("  Dominant isoforms - mean exons: %.1f, junctions: %.1f, length: %.0f bp\n",
            mean(splicing_profiles_df$n_exons_dom, na.rm = TRUE),
            mean(splicing_profiles_df$n_junctions_dom, na.rm = TRUE),
            mean(splicing_profiles_df$length_dom, na.rm = TRUE)))
cat(sprintf("  Non-dominant isoforms - mean exons: %.1f, junctions: %.1f, length: %.0f bp\n",
            mean(splicing_profiles_df$n_exons_non_dom, na.rm = TRUE),
            mean(splicing_profiles_df$n_junctions_non_dom, na.rm = TRUE),
            mean(splicing_profiles_df$length_non_dom, na.rm = TRUE)))

cat("\nTSS/TES changes:\n")
cat(sprintf("  Profiles with TSS change: %d (%.1f%%)\n",
            sum(splicing_profiles_df$tss_changed, na.rm = TRUE),
            100 * mean(splicing_profiles_df$tss_changed, na.rm = TRUE)))
cat(sprintf("  Profiles with TES change: %d (%.1f%%)\n",
            sum(splicing_profiles_df$tes_changed, na.rm = TRUE),
            100 * mean(splicing_profiles_df$tes_changed, na.rm = TRUE)))

cat("\nSplicing events detected:\n")
cat(sprintf("  Total A5SS events: %d (mean %.2f per profile)\n",
            sum(splicing_profiles_df$n_a5ss), mean(splicing_profiles_df$n_a5ss)))
cat(sprintf("  Total A3SS events: %d (mean %.2f per profile)\n",
            sum(splicing_profiles_df$n_a3ss), mean(splicing_profiles_df$n_a3ss)))
cat(sprintf("  Total Partial_IR events: %d (mean %.2f per profile)\n",
            sum(splicing_profiles_df$n_partial_ir), mean(splicing_profiles_df$n_partial_ir)))
cat(sprintf("  Total IR events: %d (mean %.2f per profile)\n",
            sum(splicing_profiles_df$n_ir), mean(splicing_profiles_df$n_ir)))
cat(sprintf("  Total SE events: %d (mean %.2f per profile)\n",
            sum(splicing_profiles_df$n_se), mean(splicing_profiles_df$n_se)))
cat(sprintf("  Total Missing_Internal events: %d (mean %.2f per profile)\n",
            sum(splicing_profiles_df$n_missing_internal), mean(splicing_profiles_df$n_missing_internal)))

# Distribution of profile types
cat("\nProfile complexity:\n")
complexity_bins <- splicing_profiles_df %>%
  mutate(
    complexity = case_when(
      n_differences == 0 ~ "identical",
      n_differences <= 2 ~ "simple (1-2 diffs)",
      n_differences <= 5 ~ "moderate (3-5 diffs)",
      TRUE ~ "complex (6+ diffs)"
    )
  ) %>%
  count(complexity) %>%
  mutate(percentage = 100 * n / sum(n))

for (i in 1:nrow(complexity_bins)) {
  cat(sprintf("  %s: %d (%.1f%%)\n",
              complexity_bins$complexity[i],
              complexity_bins$n[i],
              complexity_bins$percentage[i]))
}

# ==============================================================================
# 5. Save Output
# ==============================================================================

cat("\nSaving splicing choice profiles...\n")
# Ensure output directory exists
output_dir <- dirname(output_path)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
saveRDS(splicing_profiles_df, output_path)
cat(sprintf("  ✓ %s\n", output_path))

# Update central profile store
if (use_store && nrow(new_profiles_df) > 0) {
  current_fp <- list(
    isoform_structures_mtime = file.mtime(file.path(data_dir, "isoform_structures_filtered.rds")),
    union_exons_mtime = file.mtime(file.path(data_dir, "union_exons_filtered.rds")),
    annotated_mapping_mtime = file.mtime(file.path(data_dir, "isoform_union_exons_annotated_filtered.rds")),
    tss_tolerance = TSS_TOLERANCE,
    tes_tolerance = TES_TOLERANCE,
    splice_site_threshold = SPLICE_SITE_THRESHOLD
  )

  if (!is.null(profile_store)) {
    # Merge: new/re-verified profiles first (take precedence), then remaining store profiles
    updated_store <- bind_rows(new_profiles_df, cached_profiles_df, profile_store) %>%
      distinct(dominant_isoform_id, non_dominant_isoform_id, .keep_all = TRUE)
  } else {
    updated_store <- splicing_profiles_df
  }

  saveRDS(list(profiles = updated_store, fingerprint = current_fp),
          store_path)
  cat(sprintf("  ✓ Store updated: %d total profiles (+%d new)\n",
              nrow(updated_store), nrow(new_profiles_df)))
} else if (use_store && nrow(new_profiles_df) == 0 && n_store_cached > 0) {
  cat("  Store unchanged (all profiles were cached).\n")
}

cat("\n✓ Step 8 complete\n")
cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Splicing profiles created: %d\n", nrow(splicing_profiles_df)))
cat(sprintf("Genes covered: %d\n", n_distinct(splicing_profiles_df$gene_id)))
cat(sprintf("Mean differences per profile: %.1f exons\n",
            mean(splicing_profiles_df$n_differences)))
cat("\n")
cat("Event-level counts (variable, events can co-occur):\n")
cat(sprintf("  A5SS: %d\n", sum(splicing_profiles_df$n_a5ss)))
cat(sprintf("  A3SS: %d\n", sum(splicing_profiles_df$n_a3ss)))
cat(sprintf("  Partial_IR: %d\n", sum(splicing_profiles_df$n_partial_ir)))
cat(sprintf("  IR: %d\n", sum(splicing_profiles_df$n_ir)))
cat(sprintf("  SE: %d\n", sum(splicing_profiles_df$n_se)))
cat(sprintf("  Missing_Internal: %d\n", sum(splicing_profiles_df$n_missing_internal)))
cat(sprintf("  IR_diff: %d\n", sum(splicing_profiles_df$n_ir_diff)))
cat(sprintf("  Alt_TSS: %d\n", sum(splicing_profiles_df$n_alt_tss)))
cat(sprintf("  Alt_TES: %d\n", sum(splicing_profiles_df$n_alt_tes)))
cat("\n")
cat(sprintf("Total events detected: %d\n", sum(splicing_profiles_df$n_events)))
cat(sprintf("  Mean events per profile: %.1f\n",
            mean(splicing_profiles_df$n_events)))

if (reconstruction_check) {
  cat("\n")
  cat("───────────────────────────────────────────────────────────────────\n")
  cat("RECONSTRUCTION CHECK RESULTS\n")
  cat("───────────────────────────────────────────────────────────────────\n")
  rc_pass  <- sum(splicing_profiles_df$reconstruction_status == "PASS", na.rm = TRUE)
  rc_fail  <- sum(splicing_profiles_df$reconstruction_status == "FAIL", na.rm = TRUE)
  rc_error <- sum(splicing_profiles_df$reconstruction_status == "ERROR", na.rm = TRUE)
  rc_total <- rc_pass + rc_fail + rc_error
  cat(sprintf("  PASS:  %d/%d (%.1f%%)\n", rc_pass, rc_total, 100 * rc_pass / rc_total))
  cat(sprintf("  FAIL:  %d/%d (%.1f%%)\n", rc_fail, rc_total, 100 * rc_fail / rc_total))
  cat(sprintf("  ERROR: %d/%d (%.1f%%)\n", rc_error, rc_total, 100 * rc_error / rc_total))

  if (rc_fail > 0) {
    cat("\nFailed pairs:\n")
    failures <- splicing_profiles_df %>% filter(reconstruction_status == "FAIL")
    for (i in seq_len(min(nrow(failures), 20))) {
      f <- failures[i, ]
      cat(sprintf("  [FAIL] %s: %s vs %s — %s\n",
                  f$gene_id, f$dominant_isoform_id, f$non_dominant_isoform_id,
                  f$reconstruction_reason))
    }
    if (nrow(failures) > 20) {
      cat(sprintf("  ... and %d more\n", nrow(failures) - 20))
    }
  }

  if (rc_error > 0) {
    cat("\nError summary:\n")
    errors <- splicing_profiles_df %>%
      filter(reconstruction_status == "ERROR") %>%
      count(reconstruction_reason, sort = TRUE)
    for (i in seq_len(nrow(errors))) {
      cat(sprintf("  %s: %d\n", errors$reconstruction_reason[i], errors$n[i]))
    }
  }

  if (rc_pass == rc_total && rc_total > 0) {
    cat("\n  PERFECT RECONSTRUCTION — All pairs match!\n")
  }
}

cat("═══════════════════════════════════════════════════════════════════\n")

# ==============================================================================
# 6. Append to Run Log (if --run-log specified)
# ==============================================================================

if (!is.null(run_log_path)) {
  cat(sprintf("\nAppending to run log: %s\n", run_log_path))

  log_lines <- c(
    "======================================================================",
    "  RUN LOG — Event Detection and Splicing Profiles (Script 08)",
    "======================================================================",
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("Working directory: %s", getwd()),
    "",
    "",
    "----------------------------------------------------------------------",
    "  8. EVENT DETECTION PARAMETERS",
    "----------------------------------------------------------------------",
    "",
    sprintf("TSS tolerance:          %d bp", TSS_TOLERANCE),
    sprintf("TES tolerance:          %d bp", TES_TOLERANCE),
    sprintf("Splice site threshold:  %d bp", SPLICE_SITE_THRESHOLD),
    "",
    "Event detection hierarchy:",
    "  IR -> Boundary (A5SS/A3SS/Partial_IR) -> SE -> Missing_Internal -> Terminal",
    "",
    "Event types detected:",
    "  Alt_TSS, Alt_TES, A5SS, A3SS, SE, Missing_Internal,",
    "  IR, IR_diff_5, IR_diff_3, IR_diff_5_3, Partial_IR_5, Partial_IR_3",
    ""
  )

  # Input/output
  log_lines <- c(log_lines,
    "",
    "----------------------------------------------------------------------",
    "  9. INPUT / OUTPUT",
    "----------------------------------------------------------------------",
    ""
  )

  if (!is.null(pairs_file)) {
    log_lines <- c(log_lines,
      sprintf("Pairs file:         %s", pairs_file),
      sprintf("  Total pairs:      %d", nrow(explicit_pairs))
    )
    if (!is.null(valid_pairs)) {
      log_lines <- c(log_lines,
        sprintf("  Valid pairs:      %d", nrow(valid_pairs)),
        sprintf("  Genes:            %d", length(genes_with_dominant))
      )
    }
  } else {
    log_lines <- c(log_lines,
      "Pairs source:       auto-generated from dominant_isoforms_filtered.rds",
      sprintf("  Genes:            %d", length(genes_with_dominant))
    )
  }

  log_lines <- c(log_lines,
    "",
    sprintf("Output file:        %s", output_path),
    sprintf("  Profiles:         %d", nrow(splicing_profiles_df)),
    sprintf("  Genes:            %d", n_distinct(splicing_profiles_df$gene_id)),
    ""
  )

  if (use_store) {
    log_lines <- c(log_lines,
      "Profile store:",
      sprintf("  From store:       %d cached, %d newly computed", n_store_cached, n_store_new),
      if (n_store_reverified > 0) sprintf("  Re-verified:      %d", n_store_reverified) else NULL,
      ""
    )
  }

  # Event summary
  log_lines <- c(log_lines,
    "",
    "----------------------------------------------------------------------",
    "  10. EVENT DETECTION RESULTS",
    "----------------------------------------------------------------------",
    "",
    sprintf("Total profiles:             %d", nrow(splicing_profiles_df)),
    sprintf("Total events detected:      %d", sum(splicing_profiles_df$n_events)),
    sprintf("Mean events per profile:    %.1f", mean(splicing_profiles_df$n_events)),
    sprintf("Mean UE differences:        %.1f", mean(splicing_profiles_df$n_differences)),
    "",
    "Event counts:",
    sprintf("  Alt_TSS:          %d", sum(splicing_profiles_df$n_alt_tss)),
    sprintf("  Alt_TES:          %d", sum(splicing_profiles_df$n_alt_tes)),
    sprintf("  A5SS:             %d", sum(splicing_profiles_df$n_a5ss)),
    sprintf("  A3SS:             %d", sum(splicing_profiles_df$n_a3ss)),
    sprintf("  SE:               %d", sum(splicing_profiles_df$n_se)),
    sprintf("  Missing_Internal: %d", sum(splicing_profiles_df$n_missing_internal)),
    sprintf("  IR:               %d", sum(splicing_profiles_df$n_ir)),
    sprintf("  IR_diff:          %d", sum(splicing_profiles_df$n_ir_diff)),
    sprintf("  Partial_IR:       %d", sum(splicing_profiles_df$n_partial_ir)),
    "",
    sprintf("TSS changed: %d/%d (%.1f%%)",
            sum(splicing_profiles_df$tss_changed, na.rm = TRUE),
            nrow(splicing_profiles_df),
            100 * mean(splicing_profiles_df$tss_changed, na.rm = TRUE)),
    sprintf("TES changed: %d/%d (%.1f%%)",
            sum(splicing_profiles_df$tes_changed, na.rm = TRUE),
            nrow(splicing_profiles_df),
            100 * mean(splicing_profiles_df$tes_changed, na.rm = TRUE)),
    ""
  )

  # Complexity breakdown
  complexity_bins_log <- splicing_profiles_df %>%
    mutate(
      complexity = case_when(
        n_differences == 0 ~ "identical",
        n_differences <= 2 ~ "simple (1-2 diffs)",
        n_differences <= 5 ~ "moderate (3-5 diffs)",
        TRUE ~ "complex (6+ diffs)"
      )
    ) %>%
    count(complexity) %>%
    mutate(percentage = 100 * n / sum(n))

  log_lines <- c(log_lines,
    "Complexity distribution:",
    sprintf("  %s: %d (%.1f%%)",
            complexity_bins_log$complexity,
            complexity_bins_log$n,
            complexity_bins_log$percentage),
    ""
  )

  # Reconstruction check
  if (reconstruction_check) {
    rc_pass  <- sum(splicing_profiles_df$reconstruction_status == "PASS", na.rm = TRUE)
    rc_fail  <- sum(splicing_profiles_df$reconstruction_status == "FAIL", na.rm = TRUE)
    rc_error <- sum(splicing_profiles_df$reconstruction_status == "ERROR", na.rm = TRUE)
    rc_total <- rc_pass + rc_fail + rc_error

    log_lines <- c(log_lines,
      "",
      "----------------------------------------------------------------------",
      "  11. RECONSTRUCTION VERIFICATION",
      "----------------------------------------------------------------------",
      "",
      sprintf("PASS:   %d/%d (%.1f%%)", rc_pass, rc_total, 100 * rc_pass / max(rc_total, 1)),
      sprintf("FAIL:   %d/%d (%.1f%%)", rc_fail, rc_total, 100 * rc_fail / max(rc_total, 1)),
      sprintf("ERROR:  %d/%d (%.1f%%)", rc_error, rc_total, 100 * rc_error / max(rc_total, 1)),
      ""
    )

    if (rc_pass == rc_total && rc_total > 0) {
      log_lines <- c(log_lines, "PERFECT RECONSTRUCTION -- All pairs verified.", "")
    }

    if (rc_fail > 0) {
      failures_log <- splicing_profiles_df %>% filter(reconstruction_status == "FAIL")
      log_lines <- c(log_lines, sprintf("Failed pairs (%d):", nrow(failures_log)))
      for (i in seq_len(min(nrow(failures_log), 20))) {
        f <- failures_log[i, ]
        log_lines <- c(log_lines,
          sprintf("  %s: %s vs %s -- %s",
                  f$gene_id, f$dominant_isoform_id, f$non_dominant_isoform_id,
                  f$reconstruction_reason))
      }
      if (nrow(failures_log) > 20) {
        log_lines <- c(log_lines, sprintf("  ... and %d more", nrow(failures_log) - 20))
      }
      log_lines <- c(log_lines, "")
    }
  }

  log_lines <- c(log_lines, "")

  # Append to existing log file
  cat(log_lines, file = run_log_path, sep = "\n", append = TRUE)
  cat(sprintf("  Appended %d lines to run log.\n", length(log_lines)))
}

cat("\n")
