#!/usr/bin/env Rscript
# Union Exon Model Construction - Full Dataset with Incremental Saving
# Processes all genes in batches to avoid memory issues

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   UNION EXON MODEL CONSTRUCTION - BATCHED VERSION             ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Configuration
output_dir <- "results/isoform_transitions/v4.0_reference_based"
TSS_TES_TOLERANCE <- 20
BATCH_SIZE <- 1000  # Save every 1000 genes

# Load helper functions (same as before)
cat("Loading helper functions...\n")

# Helper functions from build_union_exon_model_full.R
calculate_dominant_isoform <- function(gene_isoforms, cpm_data, sample_metadata) {
  isoform_scores <- gene_isoforms %>%
    mutate(
      median_expression = map_dbl(isoform_ref, function(iso) {
        if (!iso %in% rownames(cpm_data)) return(0)
        iso_cpm <- cpm_data[iso, ]
        celltype_means <- tapply(log2(iso_cpm + 1), sample_metadata$cell_type, mean)
        median(celltype_means, na.rm = TRUE)
      })
    ) %>%
    arrange(desc(median_expression))
  return(isoform_scores$isoform_ref[1])
}

extract_all_exons <- function(gene_isoforms, exon_structures) {
  all_exons <- list()
  for (i in seq_len(nrow(gene_isoforms))) {
    isoform_id <- gene_isoforms$isoform_ref[i]
    exon_data <- exon_structures %>% filter(isoform_id == !!isoform_id)
    if (nrow(exon_data) == 0) next
    starts <- exon_data$exon_starts[[1]]
    ends <- exon_data$exon_ends[[1]]
    n_exons <- length(starts)
    if (n_exons == 0) next
    for (j in seq_len(n_exons)) {
      all_exons[[length(all_exons) + 1]] <- tibble(
        isoform_id = isoform_id,
        exon_index = j,
        start = starts[j],
        end = ends[j],
        is_first = (j == 1),
        is_last = (j == n_exons)
      )
    }
  }
  return(bind_rows(all_exons))
}

group_first_exons <- function(first_exons, dominant_iso) {
  # OPTION 1: Include ALL first exons, no coordinate matching required
  # Since we now use explicit Alt_TSS detection (outside union model iteration),
  # we don't need to group first exons by coordinate proximity.
  # Each first exon gets its own group to ensure it's included in the union model.

  if (nrow(first_exons) == 0) return(list())

  first_groups <- list()
  for (i in seq_len(nrow(first_exons))) {
    first_groups[[i]] <- first_exons[i, , drop = FALSE]
  }

  return(first_groups)
}

group_last_exons <- function(last_exons, dominant_iso) {
  # OPTION 1: Include ALL last exons, no coordinate matching required
  # Since we now use explicit Alt_TES detection (outside union model iteration),
  # we don't need to group last exons by coordinate proximity.
  # Each last exon gets its own group to ensure it's included in the union model.

  if (nrow(last_exons) == 0) return(list())

  last_groups <- list()
  for (i in seq_len(nrow(last_exons))) {
    last_groups[[i]] <- last_exons[i, , drop = FALSE]
  }

  return(last_groups)
}

detect_ir_exons <- function(internal_exons) {
  # Optimized O(n²) algorithm
  if (nrow(internal_exons) == 0) {
    return(list(clean_exons = internal_exons, ir_exons = NULL, n_ir_removed = 0))
  }

  exon_starts <- internal_exons$start
  exon_ends <- internal_exons$end
  ir_indices <- integer(0)

  for (i in seq_len(nrow(internal_exons))) {
    big_start <- exon_starts[i]
    big_end <- exon_ends[i]
    left_candidates <- which(exon_starts == big_start & exon_ends < big_end)
    right_candidates <- which(exon_ends == big_end & exon_starts > big_start)

    for (left_idx in left_candidates) {
      left_end <- exon_ends[left_idx]
      for (right_idx in right_candidates) {
        right_start <- exon_starts[right_idx]
        if (left_end < right_start) {
          ir_indices <- c(ir_indices, i)
          break
        }
      }
      if (i %in% ir_indices) break
    }
  }

  ir_indices <- unique(ir_indices)

  if (length(ir_indices) > 0) {
    internal_clean <- internal_exons[-ir_indices, ]
    ir_exons <- internal_exons[ir_indices, ]
  } else {
    internal_clean <- internal_exons
    ir_exons <- NULL
  }

  return(list(clean_exons = internal_clean, ir_exons = ir_exons, n_ir_removed = length(ir_indices)))
}

group_internal_exons <- function(internal_exons_clean) {
  if (nrow(internal_exons_clean) == 0) return(list())
  internal_groups <- list()
  processed <- rep(FALSE, nrow(internal_exons_clean))
  group_id <- 1
  for (i in seq_len(nrow(internal_exons_clean))) {
    if (processed[i]) next
    current <- internal_exons_clean[i, ]
    # FIX: Exclude already-processed rows to prevent duplicates
    same_group <- which(!processed & (internal_exons_clean$start == current$start | internal_exons_clean$end == current$end))
    if (length(same_group) == 0) next  # Skip if no unprocessed matches
    internal_groups[[group_id]] <- internal_exons_clean[same_group, ]
    processed[same_group] <- TRUE
    group_id <- group_id + 1
  }
  return(internal_groups)
}

filter_uncategorizable <- function(internal_groups) {
  if (length(internal_groups) < 2) return(list(keep_gene = TRUE))
  for (i in seq_len(length(internal_groups) - 1)) {
    for (j in (i + 1):length(internal_groups)) {
      variants_i <- internal_groups[[i]]
      variants_j <- internal_groups[[j]]
      for (vi in seq_len(nrow(variants_i))) {
        for (vj in seq_len(nrow(variants_j))) {
          v1 <- variants_i[vi, ]
          v2 <- variants_j[vj, ]
          overlap <- (v1$start < v2$end && v1$end > v2$start)
          share_boundary <- (v1$start == v2$start || v1$end == v2$end)
          if (overlap && !share_boundary) {
            return(list(keep_gene = FALSE, reason = "overlapping_no_shared_boundary"))
          }
        }
      }
    }
  }
  return(list(keep_gene = TRUE))
}

assign_exon_numbers <- function(first_groups, internal_groups, last_groups) {
  union_exons <- list()
  exon_num <- 1
  for (group in first_groups) {
    union_exons[[exon_num]] <- list(exon_number = exon_num, exon_type = "first", variants = group)
    exon_num <- exon_num + 1
  }
  if (length(internal_groups) > 0) {
    internal_positions <- sapply(internal_groups, function(g) min(g$start))
    sorted_internal <- internal_groups[order(internal_positions)]
    for (group in sorted_internal) {
      union_exons[[exon_num]] <- list(exon_number = exon_num, exon_type = "internal", variants = group)
      exon_num <- exon_num + 1
    }
  }
  for (group in last_groups) {
    union_exons[[exon_num]] <- list(exon_number = exon_num, exon_type = "last", variants = group)
    exon_num <- exon_num + 1
  }
  return(union_exons)
}

cat("  Helper functions loaded\n\n")

# Load data
cat("Loading data...\n")
genes_to_process <- readRDS(file.path(output_dir, "genes_for_union_model.rds"))
isoforms_to_process <- readRDS(file.path(output_dir, "isoforms_for_union_model.rds"))
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))
dge <- readRDS("results/rds/dge_isoform_2026.1.20.rds")
cpm_data <- edgeR::cpm(dge)
sample_metadata <- dge$samples
cat("  Genes to process:", nrow(genes_to_process), "\n")
cat("  Isoforms:", nrow(isoforms_to_process), "\n")
cat("  Loaded exon structures:", nrow(exon_structures), "isoforms\n")
cat("  Loaded CPM data\n\n")

# Initialize or load existing results
results_file <- file.path(output_dir, "union_exon_models_full.rds")
filtered_file <- file.path(output_dir, "filtered_genes_full.rds")

if (file.exists(results_file)) {
  cat("Loading existing results...\n")
  union_models <- readRDS(results_file)
  # Always initialize as empty list (don't load old filtered_genes)
  filtered_genes <- list()
  completed_genes <- names(union_models)
  cat("  Previously completed:", length(completed_genes), "genes\n")
  start_idx <- length(completed_genes) + 1
} else {
  union_models <- list()
  filtered_genes <- list()
  start_idx <- 1
}

# Process genes
cat("\nProcessing genes in batches of", BATCH_SIZE, "...\n")
cat("  Progress will be reported every 100 genes\n")
cat("  Starting from gene", start_idx, "\n\n")

start_time <- Sys.time()
genes_in_batch <- 0

for (gene_idx in start_idx:nrow(genes_to_process)) {
  gene <- genes_to_process$gene_id[gene_idx]

  # Progress reporting
  if (gene_idx %% 100 == 0 || gene_idx == start_idx) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    genes_per_min <- (gene_idx - start_idx + 1) / max(elapsed, 0.01)
    remaining_genes <- nrow(genes_to_process) - gene_idx
    eta_mins <- remaining_genes / genes_per_min

    cat(sprintf("[%d/%d] %s | Elapsed: %.1f min | ETA: %.1f min | RAM: %.0f MB\n",
                gene_idx, nrow(genes_to_process), gene,
                elapsed, eta_mins,
                as.numeric(system("ps -o rss= -p $PPID", intern = TRUE)) / 1024))
  }

  # Process gene (same logic as before)
  gene_isoforms <- isoforms_to_process %>%
    filter(gene_id == !!gene) %>%
    rename(isoform_ref = isoform_id)

  if (nrow(gene_isoforms) < 2) next

  dominant <- tryCatch({
    calculate_dominant_isoform(gene_isoforms, cpm_data, sample_metadata)
  }, error = function(e) NULL)

  if (is.null(dominant)) next

  all_exons <- tryCatch({
    extract_all_exons(gene_isoforms, exon_structures)
  }, error = function(e) NULL)

  if (is.null(all_exons) || nrow(all_exons) == 0) next

  first_exons <- all_exons %>% filter(is_first)
  last_exons <- all_exons %>% filter(is_last)
  internal_exons <- all_exons %>% filter(!is_first, !is_last)

  first_groups <- tryCatch({
    group_first_exons(first_exons, dominant)
  }, error = function(e) list())

  last_groups <- tryCatch({
    group_last_exons(last_exons, dominant)
  }, error = function(e) list())

  ir_result <- tryCatch({
    detect_ir_exons(internal_exons)
  }, error = function(e) list(clean_exons = internal_exons, ir_exons = NULL, n_ir_removed = 0))

  internal_groups <- tryCatch({
    group_internal_exons(ir_result$clean_exons)
  }, error = function(e) list())

  filter_result <- tryCatch({
    filter_uncategorizable(internal_groups)
  }, error = function(e) list(keep_gene = FALSE, reason = "error"))

  if (!filter_result$keep_gene) {
    filtered_genes[[gene]] <- list(
      gene_id = gene,
      reason = filter_result$reason,
      n_isoforms = nrow(gene_isoforms)
    )
    next
  }

  union_exons <- tryCatch({
    assign_exon_numbers(first_groups, internal_groups, last_groups)
  }, error = function(e) NULL)

  if (is.null(union_exons)) next

  union_models[[gene]] <- list(
    gene_id = gene,
    dominant_isoform = dominant,
    n_isoforms = nrow(gene_isoforms),
    union_exons = union_exons,
    n_ir_removed = ir_result$n_ir_removed,
    n_union_exons = length(union_exons)
  )

  genes_in_batch <- genes_in_batch + 1

  # Save incrementally every BATCH_SIZE genes
  if (genes_in_batch >= BATCH_SIZE) {
    cat(sprintf("  [CHECKPOINT] Saving results at gene %d...\n", gene_idx))
    saveRDS(union_models, results_file)
    if (length(filtered_genes) > 0) {
      saveRDS(filtered_genes, filtered_file)
    }
    genes_in_batch <- 0
  }
}

# Final save
cat("\nSaving final results...\n")
saveRDS(union_models, results_file)
cat("  Saved: union_exon_models_full.rds\n")

if (length(filtered_genes) > 0) {
  filtered_df <- bind_rows(filtered_genes)
  saveRDS(filtered_df, filtered_file)
  write_tsv(filtered_df, file.path(output_dir, "filtered_genes_full.tsv"))
  cat("  Saved: filtered_genes_full.rds\n")
  cat("  Saved: filtered_genes_full.tsv\n")
}

# Summary statistics
end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   PROCESSING COMPLETE                                         ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("  Total runtime: %.1f minutes\n", total_time))
cat(sprintf("  Genes attempted: %d\n", nrow(genes_to_process)))
cat(sprintf("  Union models built: %d\n", length(union_models)))
cat(sprintf("  Genes filtered: %d\n", length(filtered_genes)))
cat(sprintf("  Success rate: %.1f%%\n\n",
            100 * length(union_models) / nrow(genes_to_process)))

# Save summary statistics
if (length(union_models) > 0) {
  union_stats <- map_dfr(union_models, function(m) {
    tibble(
      gene_id = m$gene_id,
      n_isoforms = m$n_isoforms,
      n_union_exons = m$n_union_exons,
      n_ir_removed = m$n_ir_removed
    )
  })

  write_tsv(union_stats, file.path(output_dir, "union_model_statistics.tsv"))
  cat("  Saved: union_model_statistics.tsv\n\n")

  cat("  Union exon statistics:\n")
  cat(sprintf("    Mean union exons per gene: %.1f\n", mean(union_stats$n_union_exons)))
  cat(sprintf("    Median: %d\n", median(union_stats$n_union_exons)))
  cat(sprintf("    Range: %d - %d\n", min(union_stats$n_union_exons), max(union_stats$n_union_exons)))
  cat(sprintf("    Total IR exons removed: %d\n", sum(union_stats$n_ir_removed)))
  cat(sprintf("    Genes with IR: %d\n\n", sum(union_stats$n_ir_removed > 0)))
}

cat("═══ NEXT STEP ═══\n\n")
cat("  Run event detection on", length(union_models), "genes\n")
cat("  Command: Rscript code/detect_events_from_union_model_full.R\n\n")
