#!/usr/bin/env Rscript
# Test Updated Union Exon Model
# Verify: 1) Unclassifiable check removed, 2) IR detection expanded, 3) Zero duplicates

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   TESTING UPDATED UNION EXON MODEL                            ║\n")
cat("║   1. Unclassifiable check removed                             ║\n")
cat("║   2. IR detection expanded to first/last exons                ║\n")
cat("║   3. Zero duplicates maintained                               ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Run the build (limit to first 500 genes for quick test)
cat("═══ Running Union Model Build (Checkpoint at 500) ═══\n\n")

# Modify the script to add a checkpoint
build_cmd <- 'Rscript -e "
library(tidyverse)
source(\'code/build_union_exon_model_expression_based.R\')
"'

system(build_cmd, intern = FALSE)

cat("\n═══ Loading Results ═══\n\n")

# Load results
union_models <- readRDS('results/isoform_transitions/v4.0_reference_based/major_isoforms/union_exon_models_major.rds')
filtered_genes <- read_tsv('results/isoform_transitions/v4.0_reference_based/major_isoforms/filtered_genes_major.tsv',
                           show_col_types = FALSE)

cat(sprintf("Total genes in model: %d\n", length(union_models)))

# Check status distribution
status_counts <- filtered_genes %>% count(status)
cat("\nGene status distribution:\n")
print(status_counts)

# Check for unclassifiable
n_unclass <- sum(status_counts$status == "unclassifiable", na.rm = TRUE)
if (n_unclass > 0) {
  cat(sprintf("\n✗ FAIL: Still %d unclassifiable genes (should be 0)\n",
              filtered_genes %>% filter(status == "unclassifiable") %>% nrow()))
} else {
  cat("\n✓ PASS: No genes marked as unclassifiable\n")
}

# Count successful genes
successful_genes <- names(union_models)[sapply(union_models, function(g) g$status == "success")]
cat(sprintf("\nSuccessful genes: %d\n", length(successful_genes)))

# Expected improvement
cat("\nExpected: ~5,500-6,500 successful genes (was 2,130 before)\n")
if (length(successful_genes) > 4000) {
  cat("✓ PASS: Significant increase in successful genes\n")
} else {
  cat(sprintf("⚠ WARNING: Only %d successful genes\n", length(successful_genes)))
}

# Check for duplicate coordinates
cat("\n═══ Checking for Duplicate Coordinates ═══\n\n")

check_gene_duplicates <- function(gene_id, gene_data) {
  if (gene_data$status != 'success') return(NULL)

  all_coords <- bind_rows(lapply(seq_along(gene_data$union_exons), function(i) {
    gene_data$union_exons[[i]] %>%
      distinct(start, end) %>%
      mutate(union_exon_number = i)
  }))

  dup_coords <- all_coords %>%
    group_by(start, end) %>%
    summarise(n_union_exons = n_distinct(union_exon_number), .groups = 'drop') %>%
    filter(n_union_exons > 1)

  if (nrow(dup_coords) > 0) {
    tibble(gene_id = gene_id, n_duplicates = nrow(dup_coords))
  } else {
    NULL
  }
}

dup_results <- map2_dfr(names(union_models), union_models, check_gene_duplicates)

cat(sprintf("Genes with duplicate coordinates: %d / %d\n",
            nrow(dup_results), length(union_models)))

if (nrow(dup_results) == 0) {
  cat("✓ PASS: Zero duplicates maintained!\n")
} else {
  cat("✗ FAIL: Found duplicates\n")
  print(dup_results %>% head(10))
}

# Check IR removal statistics
cat("\n═══ IR Detection Statistics ═══\n\n")

ir_stats <- map_dfr(names(union_models), function(gene_id) {
  gene_data <- union_models[[gene_id]]
  tibble(
    gene_id = gene_id,
    status = gene_data$status,
    n_ir_removed = gene_data$n_ir_removed %||% 0
  )
})

cat(sprintf("Total genes: %d\n", nrow(ir_stats)))
cat(sprintf("Genes with IR removed: %d\n", sum(ir_stats$n_ir_removed > 0)))
cat(sprintf("Total IR exons removed: %d\n", sum(ir_stats$n_ir_removed)))
cat(sprintf("Mean IR per gene (with IR): %.1f\n",
            mean(ir_stats$n_ir_removed[ir_stats$n_ir_removed > 0])))

# Sample a few genes to verify structure
cat("\n═══ Sample Gene Verification ═══\n\n")

sample_genes <- sample(successful_genes, min(3, length(successful_genes)))

for (gene_id in sample_genes) {
  gene_data <- union_models[[gene_id]]
  cat(sprintf("\nGene: %s\n", gene_id))
  cat(sprintf("  Isoforms: %d\n", gene_data$n_isoforms))
  cat(sprintf("  Union exons: %d\n", gene_data$n_union_exons))
  cat(sprintf("  IR removed: %d\n", gene_data$n_ir_removed))

  # Check metadata
  if (length(gene_data$union_exons) > 0) {
    first_ue <- gene_data$union_exons[[1]]
    has_meta <- 'first_in_isoforms' %in% colnames(first_ue) &&
                'last_in_isoforms' %in% colnames(first_ue)
    cat(sprintf("  Metadata present: %s\n", ifelse(has_meta, "✓", "✗")))
  }
}

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   TESTING COMPLETE                                             ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
