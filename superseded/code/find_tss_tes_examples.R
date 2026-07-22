#!/usr/bin/env Rscript
# Targeted search for Alt_TSS and Alt_TES examples

library(tidyverse)

output_dir <- "results/isoform_transitions/v3.0_reference_based"
union_models <- readRDS(file.path(output_dir, "union_exon_models_full.rds"))
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))

TSS_TES_TOLERANCE <- 20

cat("═══════════════════════════════════════════════════════════════\n")
cat("TARGETED SEARCH FOR Alt_TSS AND Alt_TES EXAMPLES\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Look for genes where first/last union exons have variants with >20bp differences
alt_tss_genes <- list()
alt_tes_genes <- list()

for (gene_id in names(union_models)) {
  model <- union_models[[gene_id]]

  # Skip if too complex
  if (model$n_union_exons > 15 || model$n_isoforms > 4) next

  # Check first exon
  first_exon <- model$union_exons[[1]]
  if (first_exon$exon_type == "first" && nrow(first_exon$variants) >= 2) {
    # Calculate max difference in start positions
    starts <- first_exon$variants$start
    max_diff <- max(starts) - min(starts)

    if (max_diff > TSS_TES_TOLERANCE) {
      alt_tss_genes[[gene_id]] <- list(
        gene_id = gene_id,
        model = model,
        start_diff = max_diff
      )
    }
  }

  # Check last exon
  last_exon <- model$union_exons[[model$n_union_exons]]
  if (last_exon$exon_type == "last" && nrow(last_exon$variants) >= 2) {
    # Calculate max difference in end positions
    ends <- last_exon$variants$end
    max_diff <- max(ends) - min(ends)

    if (max_diff > TSS_TES_TOLERANCE) {
      alt_tes_genes[[gene_id]] <- list(
        gene_id = gene_id,
        model = model,
        end_diff = max_diff
      )
    }
  }

  # Stop if we have enough
  if (length(alt_tss_genes) >= 3 && length(alt_tes_genes) >= 3) break
}

cat("Found", length(alt_tss_genes), "genes with Alt_TSS potential\n")
cat("Found", length(alt_tes_genes), "genes with Alt_TES potential\n\n")

# Load helper functions
source_code <- readLines("code/detect_events_from_union_model_full.R")
helper_start <- which(grepl("# Helper Functions", source_code))[1]
helper_end <- which(grepl("cat.*Helper functions loaded", source_code))[1] - 1
eval(parse(text = source_code[helper_start:helper_end]))

# Get strand helper
get_gene_strand <- function(model, isoforms) {
  strand_lookup <- exon_structures %>%
    filter(isoform_id %in% isoforms) %>%
    select(isoform_id, strand) %>%
    distinct()

  strand_lookup %>%
    filter(isoform_id == model$dominant_isoform) %>%
    pull(strand) %>%
    .[1]
}

get_all_isoforms <- function(model) {
  unique(unlist(lapply(model$union_exons, function(ue) ue$variants$isoform_id)))
}

# Display Alt_TSS examples
if (length(alt_tss_genes) > 0) {
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Alt_TSS EXAMPLES\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  for (i in 1:min(2, length(alt_tss_genes))) {
    ex <- alt_tss_genes[[i]]
    model <- ex$model

    cat("───────────────────────────────────────────────────────────────\n")
    cat("Example", i, ":", ex$gene_id, "(TSS diff:", ex$start_diff, "bp)\n")
    cat("───────────────────────────────────────────────────────────────\n")

    isoforms <- get_all_isoforms(model)
    gene_strand <- get_gene_strand(model, isoforms)

    cat("  Union exons:", model$n_union_exons, "\n")
    cat("  Isoforms:", model$n_isoforms, "\n")
    cat("  Strand:", ifelse(gene_strand == "+", "Plus (+)", "Minus (-)"), "\n\n")

    cat("  First Union Exon Variants:\n")
    print(model$union_exons[[1]]$variants %>% select(isoform_id, start, end))

    # Pick two variants with different starts
    first_exon <- model$union_exons[[1]]$variants
    isoform_A <- first_exon$isoform_id[1]
    isoform_B <- first_exon$isoform_id[nrow(first_exon)]

    cat("\n  Transition:", isoform_A, "→", isoform_B, "\n\n")

    exons_A <- get_isoform_exons(model$union_exons, isoform_A)
    exons_B <- get_isoform_exons(model$union_exons, isoform_B)

    cat("  Isoform A structure:\n")
    print(exons_A %>% select(exon_number, exon_type, start, end, is_first, is_last))

    cat("\n  Isoform B structure:\n")
    print(exons_B %>% select(exon_number, exon_type, start, end, is_first, is_last))

    # Detect events
    events <- detect_events_from_union(model$union_exons, isoform_A, isoform_B, gene_strand)

    cat("\n  Detected Events:\n")
    print(events)

    if ("Alt_TSS" %in% events$event_type) {
      cat("\n  ✓ Alt_TSS detected!\n")
    } else {
      cat("\n  ✗ No Alt_TSS detected (may be combined with other differences)\n")
    }
    cat("\n")
  }
}

# Display Alt_TES examples
if (length(alt_tes_genes) > 0) {
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Alt_TES EXAMPLES\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  for (i in 1:min(2, length(alt_tes_genes))) {
    ex <- alt_tes_genes[[i]]
    model <- ex$model

    cat("───────────────────────────────────────────────────────────────\n")
    cat("Example", i, ":", ex$gene_id, "(TES diff:", ex$end_diff, "bp)\n")
    cat("───────────────────────────────────────────────────────────────\n")

    isoforms <- get_all_isoforms(model)
    gene_strand <- get_gene_strand(model, isoforms)

    cat("  Union exons:", model$n_union_exons, "\n")
    cat("  Isoforms:", model$n_isoforms, "\n")
    cat("  Strand:", ifelse(gene_strand == "+", "Plus (+)", "Minus (-)"), "\n\n")

    cat("  Last Union Exon Variants:\n")
    print(model$union_exons[[model$n_union_exons]]$variants %>% select(isoform_id, start, end))

    # Pick two variants with different ends
    last_exon <- model$union_exons[[model$n_union_exons]]$variants
    isoform_A <- last_exon$isoform_id[1]
    isoform_B <- last_exon$isoform_id[nrow(last_exon)]

    cat("\n  Transition:", isoform_A, "→", isoform_B, "\n\n")

    exons_A <- get_isoform_exons(model$union_exons, isoform_A)
    exons_B <- get_isoform_exons(model$union_exons, isoform_B)

    cat("  Isoform A structure:\n")
    print(exons_A %>% select(exon_number, exon_type, start, end, is_first, is_last))

    cat("\n  Isoform B structure:\n")
    print(exons_B %>% select(exon_number, exon_type, start, end, is_first, is_last))

    # Detect events
    events <- detect_events_from_union(model$union_exons, isoform_A, isoform_B, gene_strand)

    cat("\n  Detected Events:\n")
    print(events)

    if ("Alt_TES" %in% events$event_type) {
      cat("\n  ✓ Alt_TES detected!\n")
    } else {
      cat("\n  ✗ No Alt_TES detected (may be combined with other differences)\n")
    }
    cat("\n")
  }
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("SEARCH COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n")
