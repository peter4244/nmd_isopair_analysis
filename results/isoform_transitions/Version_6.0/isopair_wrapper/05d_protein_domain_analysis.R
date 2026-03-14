#!/usr/bin/env Rscript
# ==============================================================================
# 05d_protein_domain_analysis.R
#
# Comprehensive protein-level analysis of normal alternative splicing.
# Uses ALL non-NMD C4 pairs (not gene-matched), per cell type and all_samples.
#
# For each cell type:
#   1. Load unmatched C4 profiles (all genes with ≥2 non-NMD isoforms)
#   2. Run compareIsoformFrames() to get frame categories
#   3. Compare protein sequences (SQANTI FASTA) to find divergence point
#   4. Query domains (biomaRt Pfam/InterPro) and classify domain effects:
#      - gained_intact: domain entirely in GAINED comparator sequence
#      - disrupted_lost: domain spans divergence or entirely in LOST region
#   5. Test domain over-representation: which domains are preferentially
#      affected by splicing, relative to their prevalence in reference isoforms?
#
# Output: data_mashr/analysis_cache/protein_domain_full_analysis.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(Isopair)
  library(dplyr)
  library(biomaRt)
})

cat("=== Protein & Domain Analysis (All Non-NMD Pairs) ===\n")
cat("Started:", format(Sys.time()), "\n\n")

data_dir  <- "data_mashr"
cache_dir <- file.path(data_dir, "analysis_cache")
FASTA_PATH <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.faa"

cds       <- readRDS(file.path(data_dir, "cds.rds"))
cds_exons <- readRDS(file.path(data_dir, "cds_exons.rds"))

cell_types <- c("allsamples", "at", "dd", "do", "fb", "mv")
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
orf_lookup <- setNames(cds$orf_length, cds$isoform_id)

# ==============================================================================
# 1. Load unmatched C4 profiles and run frame comparison
# ==============================================================================

cat("--- Step 1: Frame comparison for unmatched C4 pairs ---\n")

fc_results <- list()
profiles_all <- list()

for (ct in cell_types) {
  prof_file <- file.path(data_dir, sprintf("profiles_c4_%s.rds", ct))
  if (!file.exists(prof_file)) { cat("  Skipping", ct, "- no profile file\n"); next }

  prof <- readRDS(prof_file)
  profiles_all[[ct]] <- prof
  cat(sprintf("  %s: %d pairs\n", ct, nrow(prof)))

  # Run frame comparison (or load from cache if gene-matched version exists)
  # We need to run fresh for the full unmatched set
  cat(sprintf("  Running compareIsoformFrames for %s...\n", ct))
  fc <- tryCatch(
    compareIsoformFrames(prof, cds_exons, cds),
    error = function(e) { cat("    Error:", e$message, "\n"); NULL })
  if (!is.null(fc)) fc_results[[ct]] <- fc
}

# ==============================================================================
# 2. Extract protein sequences for all coding C4 pairs
# ==============================================================================

cat("\n--- Step 2: Protein sequence extraction ---\n")

# Collect all unique isoform IDs needing protein sequences
all_iso_ids <- character(0)
for (ct in names(fc_results)) {
  ps <- fc_results[[ct]]$pair_summary
  coding_ps <- ps[ps$comparator_isoform_id %in% coding_ids &
                   ps$reference_isoform_id %in% coding_ids &
                   !ps$frame_category %in% c("non_coding", "no_shared_cds"), ]
  all_iso_ids <- union(all_iso_ids, c(coding_ps$reference_isoform_id,
                                       coding_ps$comparator_isoform_id))
}
cat("  Unique isoforms needing protein sequences:", length(all_iso_ids), "\n")

# Extract from FASTA using awk
id_file <- tempfile(fileext = ".txt")
writeLines(all_iso_ids, id_file)

awk_cmd <- sprintf(
  "awk 'BEGIN{while((getline id < \"%s\") > 0) want[id]=1} /^>/{if(seq && found) print hdr \"\\t\" seq; hdr=substr($1,2); found=(hdr in want); seq=\"\"; next} found{seq=seq$0} END{if(seq && found) print hdr \"\\t\" seq}' '%s'",
  id_file, FASTA_PATH)

cat("  Running awk extraction from protein FASTA...\n")
raw_output <- system(awk_cmd, intern = TRUE)
unlink(id_file)

# Parse into named vector
prot_seqs <- character(0)
if (length(raw_output) > 0) {
  split_lines <- strsplit(raw_output, "\t", fixed = TRUE)
  seq_names <- vapply(split_lines, `[`, character(1), 1)
  seq_seqs <- vapply(split_lines, `[`, character(1), 2)
  prot_seqs <- setNames(seq_seqs, seq_names)
}
cat("  Extracted", length(prot_seqs), "protein sequences\n")

# ==============================================================================
# 3. Protein comparison: find divergence point for each pair
# ==============================================================================

cat("\n--- Step 3: Protein comparison (all CTs) ---\n")

compare_proteins <- function(ref_seq, comp_seq) {
  # Character-by-character prefix comparison
  ref_aa <- strsplit(ref_seq, "")[[1]]
  comp_aa <- strsplit(comp_seq, "")[[1]]
  min_len <- min(length(ref_aa), length(comp_aa))

  if (min_len == 0) return(list(shared_aa = 0L, ref_len = length(ref_aa),
                                 comp_len = length(comp_aa), pct_novel = 100))

  # Find first mismatch
  shared <- 0L
  for (k in seq_len(min_len)) {
    if (ref_aa[k] == comp_aa[k]) shared <- shared + 1L else break
  }

  novel_aa <- length(comp_aa) - shared
  pct_novel <- round(100 * novel_aa / max(length(comp_aa), 1), 1)

  list(shared_aa = shared, ref_len = length(ref_aa),
       comp_len = length(comp_aa), pct_novel = pct_novel)
}

protein_results <- list()

for (ct in names(fc_results)) {
  ps <- fc_results[[ct]]$pair_summary
  coding_ps <- ps[ps$comparator_isoform_id %in% coding_ids &
                   ps$reference_isoform_id %in% coding_ids &
                   !ps$frame_category %in% c("non_coding", "no_shared_cds"), ]

  cat(sprintf("  %s: comparing %d coding pairs...\n", ct, nrow(coding_ps)))

  comp_rows <- list()
  for (i in seq_len(nrow(coding_ps))) {
    ref_id <- coding_ps$reference_isoform_id[i]
    comp_id <- coding_ps$comparator_isoform_id[i]

    ref_seq <- prot_seqs[ref_id]
    comp_seq <- prot_seqs[comp_id]

    if (is.na(ref_seq) || is.na(comp_seq)) {
      comp_rows[[i]] <- data.frame(
        reference_isoform_id = ref_id, comparator_isoform_id = comp_id,
        frame_category = coding_ps$frame_category[i],
        actual_shared_aa = NA_integer_, ref_protein_len = NA_integer_,
        comp_protein_len = NA_integer_, pct_novel_protein = NA_real_,
        stringsAsFactors = FALSE)
      next
    }

    res <- compare_proteins(ref_seq, comp_seq)
    comp_rows[[i]] <- data.frame(
      reference_isoform_id = ref_id, comparator_isoform_id = comp_id,
      frame_category = coding_ps$frame_category[i],
      actual_shared_aa = res$shared_aa, ref_protein_len = res$ref_len,
      comp_protein_len = res$comp_len, pct_novel_protein = res$pct_novel,
      stringsAsFactors = FALSE)
  }

  ct_comparison <- bind_rows(comp_rows)

  # Add gene_id from profiles
  gene_lk <- setNames(profiles_all[[ct]]$gene_id, profiles_all[[ct]]$comparator_isoform_id)
  ct_comparison$gene_id <- gene_lk[ct_comparison$comparator_isoform_id]
  ct_comparison$ensembl_gene <- sub("\\..*", "", ct_comparison$gene_id)
  ct_comparison$cell_type <- ct

  protein_results[[ct]] <- ct_comparison
  cat(sprintf("    Compared: %d, missing seq: %d\n",
    sum(!is.na(ct_comparison$actual_shared_aa)), sum(is.na(ct_comparison$actual_shared_aa))))
}

all_protein <- bind_rows(protein_results)
cat("  Total protein comparisons:", nrow(all_protein), "\n")

# ==============================================================================
# 4. Domain queries (biomaRt)
# ==============================================================================

cat("\n--- Step 4: Domain annotations (biomaRt) ---\n")

unique_genes <- unique(all_protein$ensembl_gene[!is.na(all_protein$ensembl_gene)])
cat("  Unique genes:", length(unique_genes), "\n")

mart <- tryCatch({
  useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 113)
}, error = function(e) {
  cat("  Ensembl v113 failed, trying current...\n")
  useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
})

cat("  Querying Pfam domains...\n")
pfam_raw <- getBM(
  attributes = c("ensembl_gene_id", "pfam", "pfam_start", "pfam_end"),
  filters = "ensembl_gene_id", values = unique_genes, mart = mart)
pfam_raw <- pfam_raw[!is.na(pfam_raw$pfam) & pfam_raw$pfam != "", ]
cat("    Pfam:", nrow(pfam_raw), "annotations for", length(unique(pfam_raw$ensembl_gene_id)), "genes\n")

cat("  Querying InterPro domains...\n")
interpro_raw <- getBM(
  attributes = c("ensembl_gene_id", "interpro", "interpro_start", "interpro_end",
                  "interpro_description"),
  filters = "ensembl_gene_id", values = unique_genes, mart = mart)
interpro_raw <- interpro_raw[!is.na(interpro_raw$interpro) & interpro_raw$interpro != "", ]
cat("    InterPro:", nrow(interpro_raw), "annotations for", length(unique(interpro_raw$ensembl_gene_id)), "genes\n")

# Combine into unified format
pfam_doms <- pfam_raw |>
  dplyr::transmute(ensembl_gene_id, domain_id = pfam,
                   domain_start = pfam_start, domain_end = pfam_end,
                   domain_name = pfam, domain_source = "Pfam")

# InterPro has descriptions
interpro_doms <- interpro_raw |>
  dplyr::transmute(ensembl_gene_id, domain_id = interpro,
                   domain_start = interpro_start, domain_end = interpro_end,
                   domain_name = interpro_description, domain_source = "InterPro")

all_domains <- dplyr::bind_rows(pfam_doms, interpro_doms) |> dplyr::distinct()
cat("  Combined unique domain annotations:", nrow(all_domains), "\n")

# ==============================================================================
# 5. Classify domain effects for each pair
# ==============================================================================

cat("\n--- Step 5: Domain effect classification ---\n")
cat("  Categories:\n")
cat("    gained_intact: domain entirely in GAINED comparator sequence (comp longer)\n")
cat("    disrupted_lost: domain spans divergence or in LOST region (non-functional)\n")
cat("    conserved: domain entirely before divergence (intact in both)\n\n")

domain_effect_rows <- list()

for (i in seq_len(nrow(all_protein))) {
  if (i %% 2000 == 0) cat("  Processed", i, "of", nrow(all_protein), "\n")

  gene <- all_protein$ensembl_gene[i]
  shared_aa <- all_protein$actual_shared_aa[i]
  ref_len <- all_protein$ref_protein_len[i]
  comp_len <- all_protein$comp_protein_len[i]
  comp_id <- all_protein$comparator_isoform_id[i]
  ref_id <- all_protein$reference_isoform_id[i]
  ct <- all_protein$cell_type[i]

  if (is.na(gene) || is.na(shared_aa)) next

  gene_doms <- all_domains[all_domains$ensembl_gene_id == gene, ]
  if (nrow(gene_doms) == 0) next

  # Classify each domain
  for (j in seq_len(nrow(gene_doms))) {
    d_start <- gene_doms$domain_start[j]
    d_end <- gene_doms$domain_end[j]

    if (is.na(d_start) || is.na(d_end)) next

    if (d_end <= shared_aa) {
      # Entirely before divergence — conserved in both
      effect <- "conserved"
    } else if (d_start > shared_aa) {
      # Entirely after divergence
      if (!is.na(comp_len) && d_start <= comp_len) {
        # Domain position exists in comparator (comparator is long enough)
        # If comparator is longer than reference in this region, it's gained
        if (!is.na(ref_len) && d_start > ref_len) {
          effect <- "gained_intact"
        } else {
          # Domain is in the divergent region — could be disrupted or novel
          effect <- "disrupted_lost"
        }
      } else {
        # Domain position beyond comparator length — lost
        effect <- "disrupted_lost"
      }
    } else {
      # Domain spans the divergence point
      effect <- "disrupted_lost"
    }

    domain_effect_rows[[length(domain_effect_rows) + 1]] <- data.frame(
      cell_type = ct,
      comparator_isoform_id = comp_id,
      reference_isoform_id = ref_id,
      gene_id = all_protein$gene_id[i],
      ensembl_gene_id = gene,
      domain_id = gene_doms$domain_id[j],
      domain_name = gene_doms$domain_name[j],
      domain_source = gene_doms$domain_source[j],
      domain_start = d_start, domain_end = d_end,
      actual_shared_aa = shared_aa,
      effect = effect,
      stringsAsFactors = FALSE)
  }
}

domain_effects <- bind_rows(domain_effect_rows)
cat("  Total domain effect annotations:", nrow(domain_effects), "\n")
cat("  Effect distribution:\n")
print(table(domain_effects$effect))

# ==============================================================================
# 6. Reference domain prevalence (for over-representation test)
# ==============================================================================

cat("\n--- Step 6: Reference domain prevalence ---\n")

# For each cell type, count how many reference isoforms contain each domain
# This is the baseline: "how common is this domain in the dominant proteome?"

ref_domain_prevalence <- list()

for (ct in names(fc_results)) {
  ps <- fc_results[[ct]]$pair_summary
  coding_ps <- ps[ps$comparator_isoform_id %in% coding_ids &
                   ps$reference_isoform_id %in% coding_ids &
                   !ps$frame_category %in% c("non_coding", "no_shared_cds"), ]

  ref_genes <- unique(sub("\\..*", "",
    profiles_all[[ct]]$gene_id[profiles_all[[ct]]$comparator_isoform_id %in%
                                coding_ps$comparator_isoform_id]))

  # Count domain instances in reference genes
  ref_doms <- all_domains[all_domains$ensembl_gene_id %in% ref_genes, ]
  ref_counts <- ref_doms |>
    dplyr::group_by(domain_id) |>
    dplyr::summarise(n_ref_genes = dplyr::n_distinct(ensembl_gene_id), .groups = "drop")
  ref_counts$cell_type <- ct
  ref_counts$n_total_genes <- length(ref_genes)

  ref_domain_prevalence[[ct]] <- ref_counts
}

ref_prevalence <- bind_rows(ref_domain_prevalence)

# ==============================================================================
# 7. Over-representation test (per CT)
# ==============================================================================

cat("\n--- Step 7: Domain over-representation testing ---\n")

overrep_results <- list()

for (ct in names(fc_results)) {
  ct_effects <- domain_effects[domain_effects$cell_type == ct, ]
  ct_ref <- ref_prevalence[ref_prevalence$cell_type == ct, ]
  n_total_genes <- ct_ref$n_total_genes[1]

  # Count disrupted genes per domain
  disrupted <- ct_effects[ct_effects$effect == "disrupted_lost", ]
  disrupted_counts <- disrupted |>
    dplyr::group_by(domain_id) |>
    dplyr::summarise(
      n_disrupted = dplyr::n_distinct(comparator_isoform_id),
      domain_name = domain_name[1],
      .groups = "drop")

  # Merge with reference prevalence
  merged <- merge(disrupted_counts, ct_ref[, c("domain_id", "n_ref_genes")],
                  by = "domain_id", all.x = TRUE)
  merged$n_ref_genes[is.na(merged$n_ref_genes)] <- 0L

  # Fisher's exact test per domain:
  # Is this domain disrupted more often than expected given its reference prevalence?
  # 2x2: (disrupted/not disrupted) x (has domain in ref / doesn't)
  n_coding_pairs <- nrow(all_protein[all_protein$cell_type == ct &
                                      !is.na(all_protein$actual_shared_aa), ])

  merged$fisher_p <- NA_real_
  merged$fisher_OR <- NA_real_
  for (i in seq_len(nrow(merged))) {
    a <- merged$n_disrupted[i]           # disrupted AND has domain
    b <- merged$n_ref_genes[i] - a       # not disrupted AND has domain
    c <- n_coding_pairs - merged$n_ref_genes[i]  # disrupted AND no domain (approx)
    d <- n_coding_pairs - a - b - c      # not disrupted AND no domain

    # Ensure non-negative
    b <- max(b, 0L)
    c <- max(c, 0L)
    d <- max(d, 0L)

    ft <- tryCatch(fisher.test(matrix(c(a, b, c, d), nrow = 2)),
                   error = function(e) NULL)
    if (!is.null(ft)) {
      merged$fisher_p[i] <- ft$p.value
      merged$fisher_OR[i] <- as.numeric(ft$estimate)
    }
  }

  merged$cell_type <- ct
  merged$n_coding_pairs <- n_coding_pairs
  merged$padj <- p.adjust(merged$fisher_p, method = "BH")
  merged$significant <- merged$padj < 0.05

  overrep_results[[ct]] <- merged

  n_sig <- sum(merged$significant, na.rm = TRUE)
  cat(sprintf("  %s: %d domains tested, %d significant (FDR < 0.05)\n",
    ct, nrow(merged), n_sig))
}

overrep_all <- bind_rows(overrep_results)

# ==============================================================================
# 8. Summary statistics
# ==============================================================================

cat("\n--- Summary ---\n")

# Per-CT protein comparison summary
for (ct in cell_types) {
  ct_prot <- all_protein[all_protein$cell_type == ct & !is.na(all_protein$actual_shared_aa), ]
  if (nrow(ct_prot) == 0) next
  cat(sprintf("  %s: %d coding pairs, median novel protein = %.1f%%\n",
    ct, nrow(ct_prot), median(ct_prot$pct_novel_protein, na.rm = TRUE)))
}

# ==============================================================================
# 9. Save results
# ==============================================================================

output <- list(
  protein_comparisons = all_protein,
  domain_effects = domain_effects,
  ref_prevalence = ref_prevalence,
  overrep_tests = overrep_all,
  all_domains = all_domains,
  metadata = list(
    run_timestamp = Sys.time(),
    cell_types = cell_types,
    fasta_path = FASTA_PATH,
    n_protein_comparisons = nrow(all_protein),
    n_domain_effects = nrow(domain_effects),
    n_unique_genes = length(unique_genes)
  )
)

output_path <- file.path(cache_dir, "protein_domain_full_analysis.rds")
saveRDS(output, output_path)
cat("\nSaved to:", output_path, "\n")
cat("Finished:", format(Sys.time()), "\n")
