#!/usr/bin/env Rscript
# ==============================================================================
# 05f_protein_analysis_v2.R
#
# Protein-level analysis using Isopair::mapSpliceToProtein() for precise
# genomic-to-protein coordinate mapping. Replaces the naive prefix comparison
# approach with per-event protein regions.
#
# Uses ALL non-NMD C4 pairs (unmatched), deduplicated across cell types.
#
# Steps:
#   1. Load & deduplicate C4 pairs
#   2. Run mapSpliceToProtein() → protein-level event regions
#   3. Load domain annotations (hmmscan if available, else biomaRt from 05d)
#   4. Classify domain effects using precise event regions
#   5. Over-representation testing (reference prevalence baseline)
#   6. PeptideAtlas mass-spec validation using event regions
#
# Output: data_mashr/analysis_cache/protein_analysis_v2.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(Isopair)
  library(dplyr)
})

cat("=== Protein Analysis v2 (mapSpliceToProtein) ===\n")
cat("Started:", format(Sys.time()), "\n\n")

data_dir <- "data_mashr"
cache_dir <- file.path(data_dir, "analysis_cache")
HMMSCAN_OUT <- file.path(cache_dir, "hmmscan_domains.domtblout")
PEPTIDE_ATLAS <- "/Users/petecastaldi/claude_projects/nmd/resources/peptideatlas/human_peptides_sorted.txt"
PROTEIN_FASTA <- file.path(data_dir, "c4_coding_proteins.faa")

cds <- readRDS(file.path(data_dir, "cds.rds"))
structures <- readRDS(file.path(data_dir, "structures.rds"))
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]

cell_types <- c("allsamples", "at", "dd", "fb", "mv")

# ==============================================================================
# 1. Load & deduplicate C4 pairs
# ==============================================================================

cat("--- Step 1: Load & deduplicate C4 pairs ---\n")

ct_membership <- list()
all_profiles <- list()

for (ct in cell_types) {
  prof_file <- file.path(data_dir, sprintf("profiles_c4_%s.rds", ct))
  if (!file.exists(prof_file)) next
  prof <- readRDS(prof_file)
  prof_coding <- prof[prof$reference_isoform_id %in% coding_ids &
                       prof$comparator_isoform_id %in% coding_ids, ]
  keys <- paste(prof_coding$reference_isoform_id,
                prof_coding$comparator_isoform_id, sep = "|")
  for (k in keys) {
    if (is.null(ct_membership[[k]])) ct_membership[[k]] <- character(0)
    ct_membership[[k]] <- c(ct_membership[[k]], ct)
  }
  all_profiles[[ct]] <- prof_coding
  cat(sprintf("  %s: %d coding pairs\n", ct, nrow(prof_coding)))
}

# Deduplicate: one copy per unique (ref, comp)
all_keys <- names(ct_membership)
seen <- character(0)
dedup_rows <- list()
for (ct in cell_types) {
  if (is.null(all_profiles[[ct]])) next
  prof <- all_profiles[[ct]]
  keys <- paste(prof$reference_isoform_id, prof$comparator_isoform_id, sep = "|")
  new_mask <- !keys %in% seen
  if (any(new_mask)) {
    dedup_rows[[length(dedup_rows) + 1]] <- prof[new_mask, ]
    seen <- c(seen, keys[new_mask])
  }
}
dedup_profiles <- dplyr::bind_rows(dedup_rows)
cat(sprintf("  Deduplicated: %d unique coding pairs\n\n", nrow(dedup_profiles)))

# ==============================================================================
# 2. Run mapSpliceToProtein()
# ==============================================================================

cat("--- Step 2: mapSpliceToProtein() ---\n")
protein_events <- mapSpliceToProtein(dedup_profiles, cds, structures, buffer_aa = 5L)
cat(sprintf("  CDS-overlapping events: %d\n", nrow(protein_events)))
cat(sprintf("  Pairs with CDS events: %d\n",
  length(unique(paste(protein_events$reference_isoform_id,
                       protein_events$comparator_isoform_id)))))
cat("  Event types:\n")
print(table(protein_events$event_type))
cat(sprintf("  Frame-affecting: %d (%.1f%%)\n",
  sum(protein_events$event_affects_frame),
  100 * mean(protein_events$event_affects_frame)))

# ==============================================================================
# 3. Load domain annotations
# ==============================================================================

cat("\n--- Step 3: Domain annotations ---\n")

# Try hmmscan first, fall back to biomaRt from 05d
use_hmmscan <- FALSE

if (file.exists(HMMSCAN_OUT)) {
  dom_lines <- readLines(HMMSCAN_OUT)
  dom_lines <- dom_lines[!grepl("^#", dom_lines)]

  if (length(dom_lines) > 1000) {
    cat("  Parsing hmmscan domains...\n")
    parsed_doms <- do.call(rbind, lapply(dom_lines, function(line) {
      fields <- strsplit(trimws(line), "\\s+")[[1]]
      if (length(fields) < 22) return(NULL)
      data.frame(
        domain_name = fields[1],
        pfam_acc = sub("\\..*", "", fields[2]),
        isoform_id = fields[4],
        env_from = as.integer(fields[20]),
        env_to = as.integer(fields[21]),
        domain_evalue = as.numeric(fields[13]),
        description = paste(fields[23:length(fields)], collapse = " "),
        stringsAsFactors = FALSE)
    }))
    n_isoforms_with_doms <- length(unique(parsed_doms$isoform_id))
    cat(sprintf("  hmmscan: %d domain hits for %d isoforms\n",
      nrow(parsed_doms), n_isoforms_with_doms))

    # Check if hmmscan is complete enough (>80% of isoforms processed)
    if (n_isoforms_with_doms > 0.5 * length(unique(c(dedup_profiles$reference_isoform_id,
                                                       dedup_profiles$comparator_isoform_id)))) {
      use_hmmscan <- TRUE
      cat("  Using hmmscan (isoform-specific) domains\n")
    } else {
      cat("  hmmscan incomplete — falling back to biomaRt\n")
    }
  }
}

if (!use_hmmscan) {
  cat("  Using biomaRt domains from 05d\n")
  prot_full <- readRDS(file.path(cache_dir, "protein_domain_full_analysis.rds"))
  # biomaRt domains are gene-level — we'll use them as approximation
  biomart_doms <- prot_full$all_domains
  cat(sprintf("  biomaRt: %d domain annotations\n", nrow(biomart_doms)))
}

# ==============================================================================
# 4. Classify domain effects using protein event regions
# ==============================================================================

cat("\n--- Step 4: Domain effect classification ---\n")

domain_effect_rows <- list()

if (use_hmmscan) {
  # Isoform-specific: look up domains on actual reference AND comparator
  iso_doms <- split(parsed_doms, parsed_doms$isoform_id)

  for (i in seq_len(nrow(protein_events))) {
    if (i %% 5000 == 0) cat("  Processed", i, "of", nrow(protein_events), "\n")

    ref_id <- protein_events$reference_isoform_id[i]
    comp_id <- protein_events$comparator_isoform_id[i]
    evt_start <- protein_events$protein_start_aa[i]
    evt_end <- protein_events$protein_end_aa[i]
    evt_type <- protein_events$event_type[i]
    direction <- protein_events$direction[i]

    # Check reference domains: any overlap with this event region?
    ref_doms <- iso_doms[[ref_id]]
    if (!is.null(ref_doms) && nrow(ref_doms) > 0) {
      for (j in seq_len(nrow(ref_doms))) {
        d_start <- ref_doms$env_from[j]
        d_end <- ref_doms$env_to[j]
        # Does domain overlap the event region?
        if (d_end < evt_start || d_start > evt_end) next  # no overlap
        domain_effect_rows[[length(domain_effect_rows) + 1]] <- data.frame(
          reference_isoform_id = ref_id,
          comparator_isoform_id = comp_id,
          event_type = evt_type,
          direction = direction,
          pfam_acc = ref_doms$pfam_acc[j],
          domain_name = ref_doms$domain_name[j],
          domain_start = d_start, domain_end = d_end,
          event_protein_start = evt_start, event_protein_end = evt_end,
          isoform_source = "reference",
          effect = if (direction == "LOSS") "disrupted_lost" else "disrupted_lost",
          stringsAsFactors = FALSE)
      }
    }

    # Check comparator domains: any in GAINED region?
    if (direction == "GAIN") {
      comp_doms <- iso_doms[[comp_id]]
      if (!is.null(comp_doms) && nrow(comp_doms) > 0) {
        for (j in seq_len(nrow(comp_doms))) {
          d_start <- comp_doms$env_from[j]
          d_end <- comp_doms$env_to[j]
          # Domain must overlap the gained event region
          if (d_end < evt_start || d_start > evt_end) next
          domain_effect_rows[[length(domain_effect_rows) + 1]] <- data.frame(
            reference_isoform_id = ref_id,
            comparator_isoform_id = comp_id,
            event_type = evt_type,
            direction = direction,
            pfam_acc = comp_doms$pfam_acc[j],
            domain_name = comp_doms$domain_name[j],
            domain_start = d_start, domain_end = d_end,
            event_protein_start = evt_start, event_protein_end = evt_end,
            isoform_source = "comparator",
            effect = "gained_intact",
            stringsAsFactors = FALSE)
        }
      }
    }
  }
} else {
  # biomaRt gene-level: map domains by gene and AA position
  gene_lk <- setNames(dedup_profiles$gene_id, dedup_profiles$comparator_isoform_id)
  protein_events$ensembl_gene <- sub("\\..*", "",
    gene_lk[protein_events$comparator_isoform_id])

  for (i in seq_len(nrow(protein_events))) {
    if (i %% 5000 == 0) cat("  Processed", i, "of", nrow(protein_events), "\n")

    gene <- protein_events$ensembl_gene[i]
    evt_start <- protein_events$protein_start_aa[i]
    evt_end <- protein_events$protein_end_aa[i]
    evt_type <- protein_events$event_type[i]
    direction <- protein_events$direction[i]
    if (is.na(gene)) next

    gene_doms <- biomart_doms[biomart_doms$ensembl_gene_id == gene, ]
    if (nrow(gene_doms) == 0) next

    for (j in seq_len(nrow(gene_doms))) {
      d_start <- gene_doms$domain_start[j]
      d_end <- gene_doms$domain_end[j]
      if (is.na(d_start) || is.na(d_end)) next
      if (d_end < evt_start || d_start > evt_end) next  # no overlap

      effect <- if (direction == "GAIN") "gained_intact" else "disrupted_lost"

      domain_effect_rows[[length(domain_effect_rows) + 1]] <- data.frame(
        reference_isoform_id = protein_events$reference_isoform_id[i],
        comparator_isoform_id = protein_events$comparator_isoform_id[i],
        event_type = evt_type,
        direction = direction,
        pfam_acc = gene_doms$domain_id[j],
        domain_name = gene_doms$domain_name[j],
        domain_start = d_start, domain_end = d_end,
        event_protein_start = evt_start, event_protein_end = evt_end,
        isoform_source = if (direction == "GAIN") "comparator" else "reference",
        effect = effect,
        stringsAsFactors = FALSE)
    }
  }
}

domain_effects <- if (length(domain_effect_rows) > 0) dplyr::bind_rows(domain_effect_rows) else data.frame()
cat(sprintf("  Domain-event overlaps: %d\n", nrow(domain_effects)))
if (nrow(domain_effects) > 0) {
  cat("  Effect distribution:\n")
  print(table(domain_effects$effect))
}

# ==============================================================================
# 5. Over-representation testing
# ==============================================================================

cat("\n--- Step 5: Over-representation testing ---\n")

if (nrow(domain_effects) > 0) {
  # Reference domain prevalence: how many unique pairs have each domain?
  if (use_hmmscan) {
    # Count reference isoforms with each domain
    ref_ids <- unique(dedup_profiles$reference_isoform_id)
    ref_dom_counts <- parsed_doms[parsed_doms$isoform_id %in% ref_ids, ] |>
      dplyr::group_by(pfam_acc) |>
      dplyr::summarise(n_ref = dplyr::n_distinct(isoform_id), .groups = "drop")
  } else {
    # biomaRt: count genes with each domain
    ref_genes <- unique(sub("\\..*", "", dedup_profiles$gene_id))
    ref_dom_counts <- biomart_doms[biomart_doms$ensembl_gene_id %in% ref_genes, ] |>
      dplyr::group_by(domain_id) |>
      dplyr::summarise(n_ref = dplyr::n_distinct(ensembl_gene_id), .groups = "drop")
    names(ref_dom_counts)[1] <- "pfam_acc"
  }

  # Disrupted counts
  disrupted <- domain_effects[domain_effects$effect == "disrupted_lost", ]
  disrupted_counts <- disrupted |>
    dplyr::group_by(pfam_acc) |>
    dplyr::summarise(
      n_disrupted = dplyr::n_distinct(paste(reference_isoform_id, comparator_isoform_id)),
      domain_name = domain_name[1],
      .groups = "drop")

  n_total <- nrow(dedup_profiles)
  overrep <- merge(disrupted_counts, ref_dom_counts, by = "pfam_acc", all.x = TRUE)
  overrep$n_ref[is.na(overrep$n_ref)] <- 0L

  # Total comparisons with ANY domain disrupted (background rate)
  n_any_disrupted <- length(unique(paste(
    domain_effects$comparator_isoform_id[domain_effects$effect == "disrupted_lost"],
    domain_effects$reference_isoform_id[domain_effects$effect == "disrupted_lost"])))

  overrep$fisher_p <- NA_real_
  overrep$fisher_OR <- NA_real_
  for (i in seq_len(nrow(overrep))) {
    # 2x2: (has domain / no domain) x (disrupted / not disrupted)
    #              Disrupted    Not disrupted
    # Has domain      a              b
    # No domain       c              d
    a <- overrep$n_disrupted[i]                          # disrupted AND has domain
    b <- max(overrep$n_ref[i] - a, 0L)                  # has domain AND not disrupted
    c_val <- max(n_any_disrupted - a, 0L)                # disrupted AND no domain (other domains)
    d_val <- max(n_total - a - b - c_val, 0L)            # no domain AND not disrupted
    ft <- tryCatch(fisher.test(matrix(c(a, c_val, b, d_val), nrow = 2)),
                   error = function(e) NULL)
    if (!is.null(ft)) {
      overrep$fisher_p[i] <- ft$p.value
      overrep$fisher_OR[i] <- as.numeric(ft$estimate)
    }
  }
  overrep$padj <- p.adjust(overrep$fisher_p, method = "BH")
  overrep$significant <- overrep$padj < 0.05
  overrep <- overrep[order(overrep$fisher_p), ]

  cat(sprintf("  Domains tested: %d, significant (FDR < 0.05): %d\n",
    nrow(overrep), sum(overrep$significant, na.rm = TRUE)))
} else {
  overrep <- data.frame()
  cat("  No domain effects to test\n")
}

# ==============================================================================
# 6. PeptideAtlas mass-spec using event regions
# ==============================================================================

cat("\n--- Step 6: PeptideAtlas mass-spec validation ---\n")

# Load protein sequences
fasta_lines <- readLines(PROTEIN_FASTA)
headers <- grep("^>", fasta_lines)
prot_seqs <- character(0)
for (i in seq_along(headers)) {
  id <- strsplit(sub("^>", "", fasta_lines[headers[i]]), "\t")[[1]][1]
  end_line <- if (i < length(headers)) headers[i + 1] - 1 else length(fasta_lines)
  seq <- paste(fasta_lines[(headers[i] + 1):end_line], collapse = "")
  prot_seqs[id] <- seq
}
cat(sprintf("  Loaded %d protein sequences\n", length(prot_seqs)))

# Load PeptideAtlas
pa_peptides <- readLines(PEPTIDE_ATLAS)
pa_set <- new.env(hash = TRUE, size = length(pa_peptides))
for (p in pa_peptides) pa_set[[p]] <- TRUE
cat(sprintf("  PeptideAtlas: %d peptides\n", length(pa_peptides)))

# Tryptic digest
tryptic_digest <- function(protein_seq, min_len = 7, max_len = 30) {
  fragments <- strsplit(protein_seq, "(?<=[KR])(?!P)", perl = TRUE)[[1]]
  fragments <- fragments[nchar(fragments) >= min_len & nchar(fragments) <= max_len]
  unique(fragments)
}

# For each pair, extract the event-affected protein region(s) and generate
# splice-specific peptides
pair_keys <- unique(paste(protein_events$reference_isoform_id,
                           protein_events$comparator_isoform_id, sep = "|"))

ms_rows <- list()
n_with_hits <- 0L
total_peptides <- 0L
total_hits <- 0L

cat(sprintf("  Processing %d pairs with CDS events...\n", length(pair_keys)))

for (pk in pair_keys) {
  parts <- strsplit(pk, "\\|")[[1]]
  ref_id <- parts[1]; comp_id <- parts[2]

  comp_seq <- prot_seqs[comp_id]
  ref_seq <- prot_seqs[ref_id]
  if (is.na(comp_seq) || is.na(ref_seq)) next

  comp_seq_clean <- sub("\\*$", "", comp_seq)
  ref_seq_clean <- sub("\\*$", "", ref_seq)

  # Get all events for this pair
  pair_events <- protein_events[protein_events$reference_isoform_id == ref_id &
                                 protein_events$comparator_isoform_id == comp_id, ]

  # Extract affected regions from comparator protein
  # For GAIN events: the comparator has extra sequence in the event region
  # For LOSS events: the comparator is missing sequence — the region around the
  #   event boundary is where the protein diverges
  novel_regions <- character(0)
  for (j in seq_len(nrow(pair_events))) {
    evt_start_aa <- pair_events$protein_start_aa_buffered[j]
    evt_end_aa <- pair_events$protein_end_aa_buffered[j]
    if (evt_start_aa > nchar(comp_seq_clean)) next
    end_pos <- min(evt_end_aa, nchar(comp_seq_clean))
    region <- substr(comp_seq_clean, evt_start_aa, end_pos)
    if (nchar(region) >= 7) novel_regions <- c(novel_regions, region)
  }

  if (length(novel_regions) == 0) next

  # Generate tryptic peptides from affected regions
  novel_peptides <- unique(unlist(lapply(novel_regions, tryptic_digest)))
  if (length(novel_peptides) == 0) next

  # Remove peptides also in reference
  ref_peptides <- tryptic_digest(ref_seq_clean)
  unique_peptides <- setdiff(novel_peptides, ref_peptides)
  if (length(unique_peptides) == 0) next

  # PeptideAtlas lookup
  hits <- sapply(unique_peptides, function(p) !is.null(pa_set[[p]]))
  n_hits <- sum(hits)
  total_peptides <- total_peptides + length(unique_peptides)
  total_hits <- total_hits + n_hits
  if (n_hits > 0) n_with_hits <- n_with_hits + 1L

  ms_rows[[length(ms_rows) + 1]] <- data.frame(
    pair_key = pk,
    comparator_isoform_id = comp_id,
    reference_isoform_id = ref_id,
    n_events = nrow(pair_events),
    n_novel_peptides = length(unique_peptides),
    n_pa_hits = n_hits,
    stringsAsFactors = FALSE)
}

ms_df <- if (length(ms_rows) > 0) dplyr::bind_rows(ms_rows) else data.frame()

cat(sprintf("\n=== Mass-Spec Results ===\n"))
cat(sprintf("  Pairs with splice-specific peptides: %d\n", nrow(ms_df)))
cat(sprintf("  Pairs with PeptideAtlas hits: %d (%.1f%%)\n",
  n_with_hits, 100 * n_with_hits / max(nrow(ms_df), 1)))
cat(sprintf("  Total splice-specific peptides: %d\n", total_peptides))
cat(sprintf("  Total PeptideAtlas hits: %d (%.1f%%)\n",
  total_hits, 100 * total_hits / max(total_peptides, 1)))

# ==============================================================================
# 7. Save
# ==============================================================================

output <- list(
  protein_events = protein_events,
  domain_effects = domain_effects,
  overrep_tests = overrep,
  parsed_doms = if (use_hmmscan) parsed_doms else NULL,
  massspec = ms_df,
  ct_membership = ct_membership,
  metadata = list(
    run_timestamp = Sys.time(),
    n_unique_pairs = nrow(dedup_profiles),
    n_cds_events = nrow(protein_events),
    n_domain_event_overlaps = nrow(domain_effects),
    domain_source = if (use_hmmscan) "hmmscan" else "biomaRt",
    n_ms_pairs_with_hits = n_with_hits,
    peptide_atlas_size = length(pa_peptides)
  )
)

output_path <- file.path(cache_dir, "protein_analysis_v2.rds")
saveRDS(output, output_path)
cat(sprintf("\nSaved to: %s\n", output_path))
cat("Finished:", format(Sys.time()), "\n")
