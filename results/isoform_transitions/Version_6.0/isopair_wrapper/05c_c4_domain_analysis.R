#!/usr/bin/env Rscript
# ==============================================================================
# 05c_c4_domain_analysis.R
#
# Expanded domain analysis for ALL coding C4 pairs (not just productive
# frameshifts). Queries biomaRt for Pfam/InterPro domains by gene, maps
# domain positions to the splice change boundary, classifies domains as
# conserved, disrupted, or in the novel region.
#
# Input:
#   - data_mashr/analysis_cache/fc_c4_*.rds (frame comparison, all CTs + allsamples)
#   - data_mashr/analysis_cache/all_c4_protein_comparison_allsamples.rds
#   - data_mashr/cds.rds
#
# Output:
#   - data_mashr/analysis_cache/c4_domain_analysis.rds
#
# Usage:
#   Rscript 05c_c4_domain_analysis.R
# ==============================================================================

library(dplyr, warn.conflicts = FALSE)

cat("=== C4 Domain Analysis (All Coding Pairs) ===\n")
cat("Started:", format(Sys.time()), "\n\n")

data_dir <- "data_mashr"
cache_dir <- file.path(data_dir, "analysis_cache")
cds <- readRDS(file.path(data_dir, "cds.rds"))

# Load protein comparison (has actual_shared_aa = divergence point per pair)
c4p <- readRDS(file.path(cache_dir, "all_c4_protein_comparison_allsamples.rds"))
comp <- c4p$comparison

# Filter to coding pairs with known frame category (exclude no_shared_cds)
comp <- comp[!comp$frame_category %in% c("no_shared_cds"), ]
cat("Coding C4 pairs (excl no_shared_cds):", nrow(comp), "\n")

# We need gene IDs to query biomaRt. Get from profiles.
profiles_c4 <- readRDS(file.path(data_dir, "profiles_c4_allsamples.rds"))
gene_lookup <- setNames(profiles_c4$gene_id, profiles_c4$comparator_isoform_id)
comp$gene_id <- gene_lookup[comp$comparator_isoform_id]

# Extract Ensembl gene IDs (strip version) for biomaRt query
comp$ensembl_gene <- sub("\\..*", "", comp$gene_id)
unique_genes <- unique(comp$ensembl_gene)
cat("Unique genes to query:", length(unique_genes), "\n\n")

# ==============================================================================
# 1. Query biomaRt for Pfam/InterPro domains
# ==============================================================================

cat("Querying biomaRt for domain annotations...\n")
library(biomaRt)

mart <- tryCatch({
  useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 113)
}, error = function(e) {
  cat("  Ensembl v113 failed, trying current...\n")
  useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
})

# Batch query: all genes at once
# Get Pfam domains with AA positions
cat("  Querying Pfam domains...\n")
pfam_raw <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                  "pfam", "pfam_start", "pfam_end"),
  filters = "ensembl_gene_id",
  values = unique_genes,
  mart = mart
)
pfam_raw <- pfam_raw[!is.na(pfam_raw$pfam) & pfam_raw$pfam != "", ]
cat("  Pfam annotations:", nrow(pfam_raw), "rows for",
    length(unique(pfam_raw$ensembl_gene_id)), "genes\n")

# Get InterPro domains
cat("  Querying InterPro domains...\n")
interpro_raw <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                  "interpro", "interpro_start", "interpro_end",
                  "interpro_description"),
  filters = "ensembl_gene_id",
  values = unique_genes,
  mart = mart
)
interpro_raw <- interpro_raw[!is.na(interpro_raw$interpro) & interpro_raw$interpro != "", ]
cat("  InterPro annotations:", nrow(interpro_raw), "rows for",
    length(unique(interpro_raw$ensembl_gene_id)), "genes\n\n")

# ==============================================================================
# 2. Map domains to splice change positions
# ==============================================================================

cat("Mapping domains to splice boundaries...\n")

# For each pair, the splice change boundary is at actual_shared_aa (AA position
# where ref and comp diverge). Domains before this position are conserved;
# domains spanning it are disrupted; domains after are in the novel region.

# Combine Pfam and InterPro into unified format
# Use gene-level domains (not transcript-specific) — take the union per gene
pfam_domains <- pfam_raw |>
  dplyr::select(ensembl_gene_id, domain_id = pfam,
         domain_start = pfam_start, domain_end = pfam_end) |>
  dplyr::mutate(domain_source = "Pfam") |>
  dplyr::distinct()

interpro_domains <- interpro_raw |>
  dplyr::select(ensembl_gene_id, domain_id = interpro,
         domain_start = interpro_start, domain_end = interpro_end) |>
  dplyr::mutate(domain_source = "InterPro") |>
  dplyr::distinct()

all_domains <- bind_rows(pfam_domains, interpro_domains)

# For each pair, classify each domain
domain_results <- list()
n_with_domains <- 0L
n_disrupted <- 0L

for (i in seq_len(nrow(comp))) {
  if (i %% 500 == 0) cat("  Processed", i, "of", nrow(comp), "\n")

  gene <- comp$ensembl_gene[i]
  shared_aa <- comp$actual_shared_aa[i]
  frame_cat <- comp$frame_category[i]

  gene_doms <- all_domains[all_domains$ensembl_gene_id == gene, ]
  if (nrow(gene_doms) == 0) next

  n_with_domains <- n_with_domains + 1L

  if (is.na(shared_aa) || shared_aa <= 0) {
    # No shared protein — all domains are in novel region
    gene_doms$position_class <- "novel_region"
  } else {
    gene_doms$position_class <- ifelse(
      gene_doms$domain_end <= shared_aa, "conserved",
      ifelse(gene_doms$domain_start > shared_aa, "novel_region", "disrupted"))
  }

  if (any(gene_doms$position_class == "disrupted")) n_disrupted <- n_disrupted + 1L

  gene_doms$comparator_isoform_id <- comp$comparator_isoform_id[i]
  gene_doms$gene_id <- comp$gene_id[i]
  gene_doms$frame_category <- frame_cat
  gene_doms$actual_shared_aa <- shared_aa

  domain_results[[length(domain_results) + 1]] <- gene_doms
}

domain_df <- bind_rows(domain_results)

cat("\nDomain analysis complete:\n")
cat("  Pairs with domain annotations:", n_with_domains, "of", nrow(comp), "\n")
cat("  Pairs with disrupted domains:", n_disrupted, "\n")
cat("  Total domain annotations:", nrow(domain_df), "\n")
cat("  Position classes:", paste(names(table(domain_df$position_class)),
    table(domain_df$position_class), sep = "=", collapse = ", "), "\n")

# ==============================================================================
# 3. Per-cell-type analysis
# ==============================================================================

cat("\nComputing per-cell-type domain disruption rates...\n")

# Load per-CT frame comparisons and protein data
ct_results <- list()
for (ct in c("at", "dd", "fb", "mv")) {
  fc_file <- file.path(cache_dir, sprintf("fc_c4_%s.rds", ct))
  if (!file.exists(fc_file)) next

  fc <- readRDS(fc_file)
  ps <- fc$pair_summary
  coding_ps <- ps[!ps$frame_category %in% c("non_coding", "no_shared_cds"), ]

  # Match to allsamples domain results by comparator_isoform_id
  ct_comp_ids <- coding_ps$comparator_isoform_id
  ct_domains <- domain_df[domain_df$comparator_isoform_id %in% ct_comp_ids, ]

  # Count pairs with any disrupted domain
  pairs_with_disrupted <- length(unique(
    ct_domains$comparator_isoform_id[ct_domains$position_class == "disrupted"]))
  pairs_with_any_domain <- length(unique(ct_domains$comparator_isoform_id))

  ct_results[[ct]] <- data.frame(
    cell_type = toupper(ct),
    n_coding_pairs = nrow(coding_ps),
    n_with_domains = pairs_with_any_domain,
    n_disrupted = pairs_with_disrupted,
    pct_disrupted = round(100 * pairs_with_disrupted / max(pairs_with_any_domain, 1), 1)
  )

  cat(sprintf("  %s: %d coding pairs, %d with domains, %d disrupted (%.1f%%)\n",
    toupper(ct), nrow(coding_ps), pairs_with_any_domain, pairs_with_disrupted,
    100 * pairs_with_disrupted / max(pairs_with_any_domain, 1)))
}

ct_summary <- bind_rows(ct_results)

# ==============================================================================
# 4. Which domains are most frequently disrupted?
# ==============================================================================

cat("\nMost frequently disrupted domains:\n")
disrupted_domains <- domain_df[domain_df$position_class == "disrupted", ]
top_disrupted <- sort(table(disrupted_domains$domain_id), decreasing = TRUE)
cat("  Top 10:\n")
for (i in seq_len(min(10, length(top_disrupted)))) {
  cat(sprintf("    %s: %d pairs\n", names(top_disrupted)[i], top_disrupted[i]))
}

# ==============================================================================
# 5. Save results
# ==============================================================================

output <- list(
  domain_annotations = domain_df,
  per_ct_summary = ct_summary,
  top_disrupted = as.data.frame(top_disrupted),
  metadata = list(
    run_timestamp = Sys.time(),
    n_pairs_analyzed = nrow(comp),
    n_pairs_with_domains = n_with_domains,
    n_pairs_disrupted = n_disrupted,
    n_unique_genes = length(unique_genes)
  )
)

output_path <- file.path(cache_dir, "c4_domain_analysis.rds")
saveRDS(output, output_path)
cat("\nSaved to:", output_path, "\n")
cat("Finished:", format(Sys.time()), "\n")
