#!/usr/bin/env Rscript
# ==============================================================================
# 05e_protein_analysis_final.R
#
# Final protein-level analysis of normal alternative splicing using:
#   - Isoform-specific domain positions from hmmscan (Pfam-A)
#   - Protein sequence comparisons (SQANTI FASTA)
#   - Mass-spec validation against PeptideAtlas (local lookup)
#
# Uses ALL non-NMD C4 pairs (unmatched), deduplicated across cell types.
# Domain effects classified per isoform pair, then mapped back to cell types.
#
# Input:
#   - data_mashr/analysis_cache/hmmscan_domains.domtblout (from hmmscan)
#   - data_mashr/c4_coding_proteins.faa (protein sequences)
#   - data_mashr/profiles_c4_*.rds (unmatched C4 profiles)
#   - PeptideAtlas: resources/peptideatlas/human_peptides_sorted.txt
#
# Output:
#   - data_mashr/analysis_cache/protein_analysis_final.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(Isopair)
  library(dplyr)
})

cat("=== Final Protein & Domain Analysis ===\n")
cat("Started:", format(Sys.time()), "\n\n")

data_dir <- "data_mashr"
cache_dir <- file.path(data_dir, "analysis_cache")
PROTEIN_FASTA <- file.path(data_dir, "c4_coding_proteins.faa")
HMMSCAN_OUT <- file.path(cache_dir, "hmmscan_domains.domtblout")
PEPTIDE_ATLAS <- "/Users/petecastaldi/claude_projects/nmd/resources/peptideatlas/human_peptides_sorted.txt"

cds <- readRDS(file.path(data_dir, "cds.rds"))
cds_exons <- readRDS(file.path(data_dir, "cds_exons.rds"))
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]

cell_types <- c("allsamples", "at", "dd", "fb", "mv")

# ==============================================================================
# 1. Load and deduplicate C4 pairs across cell types
# ==============================================================================

cat("--- Step 1: Load and deduplicate C4 pairs ---\n")

all_pairs <- list()
ct_membership <- list()  # which CTs each pair belongs to

for (ct in cell_types) {
  prof_file <- file.path(data_dir, sprintf("profiles_c4_%s.rds", ct))
  if (!file.exists(prof_file)) next
  prof <- readRDS(prof_file)
  # Filter to coding pairs
  prof_coding <- prof[prof$reference_isoform_id %in% coding_ids &
                       prof$comparator_isoform_id %in% coding_ids, ]
  keys <- paste(prof_coding$reference_isoform_id,
                prof_coding$comparator_isoform_id, sep = "|")
  for (k in keys) {
    if (is.null(ct_membership[[k]])) ct_membership[[k]] <- character(0)
    ct_membership[[k]] <- c(ct_membership[[k]], ct)
  }
  all_pairs[[ct]] <- prof_coding
  cat(sprintf("  %s: %d coding pairs\n", ct, nrow(prof_coding)))
}

# Deduplicate: keep one copy of each unique (ref, comp) pair
all_keys <- names(ct_membership)
dedup_pairs <- do.call(rbind, lapply(all_keys, function(k) {
  parts <- strsplit(k, "\\|")[[1]]
  data.frame(reference_isoform_id = parts[1], comparator_isoform_id = parts[2],
             stringsAsFactors = FALSE)
}))
# Add gene_id from any CT that has this pair
gene_lk <- character(0)
for (ct in names(all_pairs)) {
  g <- setNames(all_pairs[[ct]]$gene_id, paste(all_pairs[[ct]]$reference_isoform_id,
    all_pairs[[ct]]$comparator_isoform_id, sep = "|"))
  gene_lk[names(g)] <- g
}
dedup_pairs$gene_id <- gene_lk[paste(dedup_pairs$reference_isoform_id,
  dedup_pairs$comparator_isoform_id, sep = "|")]

cat(sprintf("  Total unique coding pairs: %d (from %d across CTs)\n",
  nrow(dedup_pairs), length(all_keys)))

# ==============================================================================
# 2. Load protein sequences and compare
# ==============================================================================

cat("\n--- Step 2: Protein sequence comparison ---\n")

# Read protein FASTA
fasta_lines <- readLines(PROTEIN_FASTA)
headers <- grep("^>", fasta_lines)
prot_seqs <- character(0)
for (i in seq_along(headers)) {
  id <- sub("^>", "", fasta_lines[headers[i]])
  end_line <- if (i < length(headers)) headers[i + 1] - 1 else length(fasta_lines)
  seq <- paste(fasta_lines[(headers[i] + 1):end_line], collapse = "")
  prot_seqs[id] <- seq
}
cat(sprintf("  Loaded %d protein sequences\n", length(prot_seqs)))

# Compare each deduplicated pair
cat("  Comparing proteins...\n")
dedup_pairs$actual_shared_aa <- NA_integer_
dedup_pairs$ref_protein_len <- NA_integer_
dedup_pairs$comp_protein_len <- NA_integer_
dedup_pairs$pct_novel_protein <- NA_real_

for (i in seq_len(nrow(dedup_pairs))) {
  if (i %% 5000 == 0) cat("    Processed", i, "of", nrow(dedup_pairs), "\n")
  ref_seq <- prot_seqs[dedup_pairs$reference_isoform_id[i]]
  comp_seq <- prot_seqs[dedup_pairs$comparator_isoform_id[i]]
  if (is.na(ref_seq) || is.na(comp_seq)) next

  ref_aa <- strsplit(ref_seq, "")[[1]]
  comp_aa <- strsplit(comp_seq, "")[[1]]
  min_len <- min(length(ref_aa), length(comp_aa))

  shared <- 0L
  if (min_len > 0) {
    for (k in seq_len(min_len)) {
      if (ref_aa[k] == comp_aa[k]) shared <- shared + 1L else break
    }
  }

  dedup_pairs$actual_shared_aa[i] <- shared
  dedup_pairs$ref_protein_len[i] <- length(ref_aa)
  dedup_pairs$comp_protein_len[i] <- length(comp_aa)
  novel_aa <- length(comp_aa) - shared
  dedup_pairs$pct_novel_protein[i] <- round(100 * novel_aa / max(length(comp_aa), 1), 1)
}

n_compared <- sum(!is.na(dedup_pairs$actual_shared_aa))
cat(sprintf("  Compared: %d pairs\n", n_compared))

# Flag pairs where proteins are identical despite splice event
dedup_pairs$proteins_identical <- dedup_pairs$actual_shared_aa == dedup_pairs$comp_protein_len &
  dedup_pairs$comp_protein_len == dedup_pairs$ref_protein_len
n_identical <- sum(dedup_pairs$proteins_identical, na.rm = TRUE)
cat(sprintf("  Proteins identical despite splicing: %d (%.1f%%)\n",
  n_identical, 100 * n_identical / n_compared))

# ==============================================================================
# 3. Parse hmmscan domain results
# ==============================================================================

cat("\n--- Step 3: Parse hmmscan domains ---\n")

if (!file.exists(HMMSCAN_OUT)) {
  stop("hmmscan output not found: ", HMMSCAN_OUT, "\n  Run hmmscan first.")
}

# Parse domtblout format (space-separated, comment lines start with #)
dom_lines <- readLines(HMMSCAN_OUT)
dom_lines <- dom_lines[!grepl("^#", dom_lines)]

if (length(dom_lines) == 0) {
  stop("No domain hits in hmmscan output")
}

# domtblout columns (first 22 are fixed):
# 1:target_name 2:target_acc 3:tlen 4:query_name 5:query_acc 6:qlen
# 7:full_E 8:full_score 9:full_bias 10:dom_# 11:dom_of
# 12:cE 13:iE 14:dom_score 15:dom_bias
# 16:hmm_from 17:hmm_to 18:ali_from 19:ali_to 20:env_from 21:env_to
# 22:acc 23+:description

parsed_doms <- do.call(rbind, lapply(dom_lines, function(line) {
  fields <- strsplit(trimws(line), "\\s+")[[1]]
  if (length(fields) < 22) return(NULL)
  data.frame(
    domain_name = fields[1],     # Pfam name (e.g., "Tubulin")
    domain_acc = fields[2],      # Pfam accession (e.g., "PF00091.32")
    isoform_id = fields[4],      # query = our isoform ID
    isoform_len = as.integer(fields[6]),
    domain_evalue = as.numeric(fields[13]),  # per-domain i-Evalue
    domain_score = as.numeric(fields[14]),
    ali_from = as.integer(fields[18]),  # alignment start in query (AA position)
    ali_to = as.integer(fields[19]),    # alignment end in query
    env_from = as.integer(fields[20]), # envelope start (slightly wider)
    env_to = as.integer(fields[21]),   # envelope end
    accuracy = as.numeric(fields[22]),
    description = paste(fields[23:length(fields)], collapse = " "),
    stringsAsFactors = FALSE)
}))

# Use envelope coordinates (env_from, env_to) — more conservative than alignment
# Strip Pfam version from accession
parsed_doms$pfam_acc <- sub("\\..*", "", parsed_doms$domain_acc)

cat(sprintf("  Parsed %d domain hits for %d isoforms\n",
  nrow(parsed_doms), length(unique(parsed_doms$isoform_id))))
cat(sprintf("  Unique Pfam domains: %d\n", length(unique(parsed_doms$pfam_acc))))

# ==============================================================================
# 4. Classify domain effects for each pair
# ==============================================================================

cat("\n--- Step 4: Domain effect classification ---\n")

# Build per-isoform domain lookup
iso_domains <- split(parsed_doms, parsed_doms$isoform_id)

domain_effect_rows <- list()
n_processed <- 0L

for (i in seq_len(nrow(dedup_pairs))) {
  if (i %% 5000 == 0) cat("    Processed", i, "of", nrow(dedup_pairs), "\n")

  shared_aa <- dedup_pairs$actual_shared_aa[i]
  if (is.na(shared_aa)) next

  ref_id <- dedup_pairs$reference_isoform_id[i]
  comp_id <- dedup_pairs$comparator_isoform_id[i]
  ref_len <- dedup_pairs$ref_protein_len[i]
  comp_len <- dedup_pairs$comp_protein_len[i]

  # Get domains for BOTH reference and comparator
  ref_doms <- iso_domains[[ref_id]]
  comp_doms <- iso_domains[[comp_id]]

  # Process reference domains: which are conserved vs disrupted/lost in comparator?
  if (!is.null(ref_doms) && nrow(ref_doms) > 0) {
    for (j in seq_len(nrow(ref_doms))) {
      d_start <- ref_doms$env_from[j]
      d_end <- ref_doms$env_to[j]

      if (d_end <= shared_aa) {
        effect <- "conserved"
      } else if (d_start > shared_aa) {
        effect <- "disrupted_lost"  # in reference but after divergence = lost in comp
      } else {
        effect <- "disrupted_lost"  # spans divergence
      }

      domain_effect_rows[[length(domain_effect_rows) + 1]] <- data.frame(
        pair_key = paste(ref_id, comp_id, sep = "|"),
        isoform_source = "reference",
        pfam_acc = ref_doms$pfam_acc[j],
        domain_name = ref_doms$domain_name[j],
        domain_start = d_start, domain_end = d_end,
        actual_shared_aa = shared_aa,
        effect = effect,
        stringsAsFactors = FALSE)
    }
  }

  # Process comparator domains: which are gained (not in reference)?
  if (!is.null(comp_doms) && nrow(comp_doms) > 0) {
    for (j in seq_len(nrow(comp_doms))) {
      d_start <- comp_doms$env_from[j]
      d_end <- comp_doms$env_to[j]

      if (d_end <= shared_aa) {
        # Conserved region — domain in shared prefix (already counted from ref side)
        next
      } else if (d_start > shared_aa) {
        # Entirely after divergence in comparator = gained intact
        effect <- "gained_intact"
      } else {
        # Spans divergence in comparator — partially new
        effect <- "gained_intact"  # domain is present in comp but disrupted from ref perspective
      }

      domain_effect_rows[[length(domain_effect_rows) + 1]] <- data.frame(
        pair_key = paste(ref_id, comp_id, sep = "|"),
        isoform_source = "comparator",
        pfam_acc = comp_doms$pfam_acc[j],
        domain_name = comp_doms$domain_name[j],
        domain_start = d_start, domain_end = d_end,
        actual_shared_aa = shared_aa,
        effect = effect,
        stringsAsFactors = FALSE)
    }
  }
  n_processed <- n_processed + 1L
}

domain_effects <- bind_rows(domain_effect_rows)
cat(sprintf("  Processed %d pairs\n", n_processed))
cat(sprintf("  Total domain effects: %d\n", nrow(domain_effects)))
cat("  Effect distribution:\n")
print(table(domain_effects$effect))

# ==============================================================================
# 5. Over-representation testing
# ==============================================================================

cat("\n--- Step 5: Domain over-representation testing ---\n")

# Reference domain prevalence: how many unique pairs have each domain in the reference?
ref_effects <- domain_effects[domain_effects$isoform_source == "reference", ]
ref_prevalence <- ref_effects |>
  dplyr::group_by(pfam_acc) |>
  dplyr::summarise(n_ref_pairs = dplyr::n_distinct(pair_key), .groups = "drop")

# Disrupted domain counts
disrupted <- ref_effects[ref_effects$effect == "disrupted_lost", ]
disrupted_counts <- disrupted |>
  dplyr::group_by(pfam_acc) |>
  dplyr::summarise(
    n_disrupted = dplyr::n_distinct(pair_key),
    domain_name = domain_name[1],
    .groups = "drop")

# Merge and test
n_total_pairs <- nrow(dedup_pairs[!is.na(dedup_pairs$actual_shared_aa), ])
overrep <- merge(disrupted_counts, ref_prevalence, by = "pfam_acc", all.x = TRUE)

overrep$fisher_p <- NA_real_
overrep$fisher_OR <- NA_real_
for (i in seq_len(nrow(overrep))) {
  a <- overrep$n_disrupted[i]
  b <- overrep$n_ref_pairs[i] - a
  c <- n_total_pairs - overrep$n_ref_pairs[i]
  d <- n_total_pairs - a - b - c
  b <- max(b, 0L); c <- max(c, 0L); d <- max(d, 0L)
  ft <- tryCatch(fisher.test(matrix(c(a, b, c, d), nrow = 2)), error = function(e) NULL)
  if (!is.null(ft)) {
    overrep$fisher_p[i] <- ft$p.value
    overrep$fisher_OR[i] <- as.numeric(ft$estimate)
  }
}
overrep$padj <- p.adjust(overrep$fisher_p, method = "BH")
overrep$significant <- overrep$padj < 0.05
overrep <- overrep[order(overrep$fisher_p), ]

n_sig <- sum(overrep$significant, na.rm = TRUE)
cat(sprintf("  Domains tested: %d, significant (FDR < 0.05): %d\n",
  nrow(overrep), n_sig))

# Per-CT over-representation
ct_overrep <- list()
for (ct in cell_types) {
  ct_keys <- names(ct_membership)[sapply(ct_membership, function(x) ct %in% x)]
  ct_ref <- ref_effects[ref_effects$pair_key %in% ct_keys, ]
  ct_disrupted <- ct_ref[ct_ref$effect == "disrupted_lost", ]
  n_ct_pairs <- sum(dedup_pairs$actual_shared_aa[paste(dedup_pairs$reference_isoform_id,
    dedup_pairs$comparator_isoform_id, sep = "|") %in% ct_keys] > -Inf, na.rm = TRUE)

  ct_dis_counts <- ct_disrupted |>
    dplyr::group_by(pfam_acc) |>
    dplyr::summarise(n_disrupted = dplyr::n_distinct(pair_key), .groups = "drop")
  ct_dis_counts$cell_type <- ct
  ct_dis_counts$n_ct_pairs <- n_ct_pairs
  ct_overrep[[ct]] <- ct_dis_counts
}
ct_overrep_all <- bind_rows(ct_overrep)

# ==============================================================================
# 6. Mass-spec validation (PeptideAtlas)
# ==============================================================================

cat("\n--- Step 6: Mass-spec validation (PeptideAtlas) ---\n")

# Load PeptideAtlas peptides
pa_peptides <- readLines(PEPTIDE_ATLAS)
cat(sprintf("  PeptideAtlas peptides loaded: %d\n", length(pa_peptides)))

# For each pair with novel protein, generate tryptic peptides from the gained region
# Trypsin cuts after K or R (not before P)
tryptic_digest <- function(protein_seq, min_len = 7, max_len = 30) {
  # Split at K or R (not followed by P)
  fragments <- strsplit(protein_seq, "(?<=[KR])(?!P)", perl = TRUE)[[1]]
  fragments <- fragments[nchar(fragments) >= min_len & nchar(fragments) <= max_len]
  unique(fragments)
}

ms_results <- list()
n_with_novel <- 0L
n_with_hits <- 0L

for (i in seq_len(nrow(dedup_pairs))) {
  if (i %% 5000 == 0) cat("    Processed", i, "of", nrow(dedup_pairs), "\n")

  shared_aa <- dedup_pairs$actual_shared_aa[i]
  comp_id <- dedup_pairs$comparator_isoform_id[i]
  ref_id <- dedup_pairs$reference_isoform_id[i]
  comp_seq <- prot_seqs[comp_id]
  ref_seq <- prot_seqs[ref_id]

  if (is.na(shared_aa) || is.na(comp_seq) || is.na(ref_seq)) next
  if (shared_aa >= nchar(comp_seq)) next  # no novel region

  # Extract novel region (after divergence point)
  novel_region <- substr(comp_seq, shared_aa + 1, nchar(comp_seq))
  if (nchar(novel_region) < 7) next  # too short for tryptic peptides

  n_with_novel <- n_with_novel + 1L

  # Generate tryptic peptides from novel region
  novel_peptides <- tryptic_digest(novel_region)
  if (length(novel_peptides) == 0) next

  # Filter to peptides NOT in reference protein
  ref_peptides <- tryptic_digest(ref_seq)
  unique_peptides <- setdiff(novel_peptides, ref_peptides)
  if (length(unique_peptides) == 0) next

  # Check against PeptideAtlas (binary search via %in% on sorted vector)
  hits <- unique_peptides %in% pa_peptides
  n_hits <- sum(hits)

  if (n_hits > 0) n_with_hits <- n_with_hits + 1L

  ms_results[[length(ms_results) + 1]] <- data.frame(
    pair_key = paste(ref_id, comp_id, sep = "|"),
    comparator_isoform_id = comp_id,
    reference_isoform_id = ref_id,
    gene_id = dedup_pairs$gene_id[i],
    n_novel_peptides = length(unique_peptides),
    n_pa_hits = n_hits,
    hit_peptides = if (n_hits > 0) paste(unique_peptides[hits], collapse = ";") else "",
    stringsAsFactors = FALSE)
}

ms_df <- bind_rows(ms_results)

cat(sprintf("  Pairs with novel protein region: %d\n", n_with_novel))
cat(sprintf("  Pairs queried (with unique tryptic peptides): %d\n", nrow(ms_df)))
cat(sprintf("  Pairs with PeptideAtlas hits: %d (%.1f%%)\n",
  n_with_hits, 100 * n_with_hits / max(nrow(ms_df), 1)))
cat(sprintf("  Total unique peptides queried: %d\n", sum(ms_df$n_novel_peptides)))
cat(sprintf("  Total PeptideAtlas hits: %d\n", sum(ms_df$n_pa_hits)))

# ==============================================================================
# 7. Save results
# ==============================================================================

output <- list(
  dedup_pairs = dedup_pairs,
  ct_membership = ct_membership,
  hmmscan_domains = parsed_doms,
  domain_effects = domain_effects,
  overrep_tests = overrep,
  ct_overrep = ct_overrep_all,
  massspec = ms_df,
  metadata = list(
    run_timestamp = Sys.time(),
    n_unique_pairs = nrow(dedup_pairs),
    n_compared = n_compared,
    n_proteins_identical = n_identical,
    n_domain_hits = nrow(parsed_doms),
    n_domain_effects = nrow(domain_effects),
    n_overrep_significant = n_sig,
    n_ms_pairs_with_hits = n_with_hits,
    peptide_atlas_size = length(pa_peptides),
    cell_types = cell_types
  )
)

output_path <- file.path(cache_dir, "protein_analysis_final.rds")
saveRDS(output, output_path)
cat("\nSaved to:", output_path, "\n")
cat("Finished:", format(Sys.time()), "\n")
