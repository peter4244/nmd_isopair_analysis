#!/usr/bin/env Rscript
# Create comprehensive gene and isoform inclusion/exclusion tracking report

library(tidyverse)

output_dir <- "results/isoform_transitions/v3.0_reference_based"

cat("═══════════════════════════════════════════════════════════════\n")
cat("GENE AND ISOFORM TRACKING REPORT\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Load all relevant data
union_models_all <- readRDS(file.path(output_dir, "union_exon_models_full.rds"))
isoforms_for_union <- readRDS(file.path(output_dir, "isoforms_for_union_model.rds"))
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))

# Get all unique genes and isoforms from source data
all_genes <- unique(isoforms_for_union$gene_id)
all_isoforms <- unique(isoforms_for_union$isoform_id)

cat("STARTING UNIVERSE:\n")
cat("  Total genes:", length(all_genes), "\n")
cat("  Total isoforms:", nrow(isoforms_for_union), "\n\n")

# Track exclusions at each stage

# Stage 1: Union model construction exclusions
genes_in_union_models <- names(union_models_all)
genes_no_union_model <- setdiff(all_genes, genes_in_union_models)

cat("STAGE 1: Union model construction\n")
cat("  Genes with union models:", length(genes_in_union_models), "\n")
cat("  Genes excluded (no union model):", length(genes_no_union_model), "\n")
if (length(genes_no_union_model) > 0) {
  cat("  Reason: Single isoform or model construction failure\n")
}
cat("\n")

# Stage 2: Union exon count filter (>1 and <=20)
exon_count_valid <- sapply(union_models_all, function(m) {
  m$n_union_exons > 1 && m$n_union_exons <= 20
})
genes_pass_exon_filter <- names(union_models_all)[exon_count_valid]
genes_fail_exon_filter <- names(union_models_all)[!exon_count_valid]

cat("STAGE 2: Union exon count filter (>1 and ≤20)\n")
cat("  Genes passing:", length(genes_pass_exon_filter), "\n")
cat("  Genes excluded:", length(genes_fail_exon_filter), "\n")

# Count by reason
too_few <- sum(sapply(union_models_all[genes_fail_exon_filter], function(m) m$n_union_exons <= 1))
too_many <- sum(sapply(union_models_all[genes_fail_exon_filter], function(m) m$n_union_exons > 20))
cat(sprintf("    Too few exons (≤1): %d\n", too_few))
cat(sprintf("    Too many exons (>20): %d\n", too_many))
cat("\n")

# Stage 3: Fusion gene filter
fusion_pattern <- "ENSG[0-9]+_ENSG[0-9]+"
is_fusion <- grepl(fusion_pattern, genes_pass_exon_filter)
fusion_genes <- genes_pass_exon_filter[is_fusion]
genes_after_fusion_filter <- genes_pass_exon_filter[!is_fusion]

cat("STAGE 3: Fusion gene exclusion\n")
cat("  Fusion genes excluded:", length(fusion_genes), "\n")
cat("  Genes remaining:", length(genes_after_fusion_filter), "\n")

# Count isoforms in fusion genes
fusion_isoforms <- sum(sapply(fusion_genes, function(g) union_models_all[[g]]$n_isoforms))
cat("  Isoforms in fusion genes:", fusion_isoforms, "\n")
cat("\n")

# Final counts
cat("FINAL INCLUSION FOR EVENT DETECTION:\n")
cat("  Genes included:", length(genes_after_fusion_filter), "\n")
cat("  Genes excluded (total):", length(all_genes) - length(genes_after_fusion_filter), "\n")

total_isoforms_included <- sum(sapply(genes_after_fusion_filter, function(g) {
  union_models_all[[g]]$n_isoforms
}))
cat("  Isoforms included:", total_isoforms_included, "\n")
cat("  Isoforms excluded (total):", nrow(isoforms_for_union) - total_isoforms_included, "\n\n")

# Create detailed exclusion table
exclusions <- tibble(
  gene_id = character(),
  exclusion_stage = character(),
  exclusion_reason = character(),
  n_union_exons = integer(),
  n_isoforms = integer()
)

# Add Stage 1 exclusions
if (length(genes_no_union_model) > 0) {
  exclusions <- bind_rows(exclusions, tibble(
    gene_id = genes_no_union_model,
    exclusion_stage = "1_union_model",
    exclusion_reason = "no_union_model_created",
    n_union_exons = NA_integer_,
    n_isoforms = NA_integer_
  ))
}

# Add Stage 2 exclusions
if (length(genes_fail_exon_filter) > 0) {
  exclusions <- bind_rows(exclusions, tibble(
    gene_id = genes_fail_exon_filter,
    exclusion_stage = "2_exon_count",
    exclusion_reason = ifelse(
      sapply(union_models_all[genes_fail_exon_filter], function(m) m$n_union_exons <= 1),
      "too_few_exons",
      "too_many_exons"
    ),
    n_union_exons = sapply(union_models_all[genes_fail_exon_filter], function(m) m$n_union_exons),
    n_isoforms = sapply(union_models_all[genes_fail_exon_filter], function(m) m$n_isoforms)
  ))
}

# Add Stage 3 exclusions
if (length(fusion_genes) > 0) {
  exclusions <- bind_rows(exclusions, tibble(
    gene_id = fusion_genes,
    exclusion_stage = "3_fusion_gene",
    exclusion_reason = "fusion_gene",
    n_union_exons = sapply(union_models_all[fusion_genes], function(m) m$n_union_exons),
    n_isoforms = sapply(union_models_all[fusion_genes], function(m) m$n_isoforms)
  ))
}

# Create inclusion table
inclusions <- tibble(
  gene_id = genes_after_fusion_filter,
  included = TRUE,
  n_union_exons = sapply(union_models_all[genes_after_fusion_filter], function(m) m$n_union_exons),
  n_isoforms = sapply(union_models_all[genes_after_fusion_filter], function(m) m$n_isoforms)
)

# Save reports
saveRDS(exclusions, file.path(output_dir, "gene_exclusions_detailed.rds"))
write_tsv(exclusions, file.path(output_dir, "gene_exclusions_detailed.tsv"))

saveRDS(inclusions, file.path(output_dir, "gene_inclusions.rds"))
write_tsv(inclusions, file.path(output_dir, "gene_inclusions.tsv"))

cat("REPORTS SAVED:\n")
cat("  gene_exclusions_detailed.tsv/.rds\n")
cat("  gene_inclusions.tsv/.rds\n\n")

# Summary by exclusion reason
cat("EXCLUSION SUMMARY BY REASON:\n")
exclusion_summary <- exclusions %>%
  count(exclusion_reason, sort = TRUE)
print(exclusion_summary)

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("TRACKING REPORT COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n")
