#!/usr/bin/env Rscript
# ==============================================================================
# 05r_ref_atg_analysis.R
#
# Reference-ATG-anchored NMD analysis. Backend: Isopair::traceReferenceAtg().
#
# For each gene-matched NMD pair (C2) and Control pair (C4):
#   1. Trace the reference ATG through the comparator (Isopair::traceReferenceAtg)
#   2. For effectively-PTC cases, attribute the responsible splice event
#      using the shared attribute_ptc_events() / attribute_3utr_splice()
#      helpers in analysis_functions.R
#
# Categories produced (from Isopair::traceReferenceAtg):
#   effectively_ptc, truncated_no_ejc, no_downstream_ejc, ref_atg_lost,
#   no_ref_cds, mapping_failed
#
# Inputs:
#   - data_mashr/cds.rds, structures.rds, profiles_c2/c4_allsamples.rds
#   - data_mashr/analysis_cache/ptc_c2_allsamples.rds
#   - data_mashr/analysis_cache/fw_c2_allsamples.rds
#   - SQANTI corrected FASTA
#
# Output:
#   - data_mashr/analysis_cache/ref_atg_analysis.rds
#
# Usage:
#   Rscript 05r_ref_atg_analysis.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Isopair)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

FASTA_PATH <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta"
OUTPUT_PATH <- "data_mashr/analysis_cache/ref_atg_analysis.rds"

cat("=== Reference-ATG-Anchored NMD Analysis (Isopair::traceReferenceAtg) ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Load data and define populations
# ==============================================================================
cat("Loading data...\n")

cds <- readRDS("data_mashr/cds.rds")
structures <- readRDS("data_mashr/structures.rds")
profiles_c2 <- readRDS("data_mashr/profiles_c2_allsamples.rds")
profiles_c4 <- readRDS("data_mashr/profiles_c4_allsamples.rds")
ptc_assoc <- readRDS("data_mashr/analysis_cache/ptc_c2_allsamples.rds")

# Gene-match
shared <- inner_join(
  unique(profiles_c2[, c("gene_id", "reference_isoform_id")]),
  unique(profiles_c4[, c("gene_id", "reference_isoform_id")]),
  by = c("gene_id", "reference_isoform_id"))

gm_c2 <- semi_join(profiles_c2, shared, by = c("gene_id", "reference_isoform_id"))
gm_c4 <- semi_join(profiles_c4, shared, by = c("gene_id", "reference_isoform_id"))

coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
c2_coding <- gm_c2 %>%
  filter(reference_isoform_id %in% coding_ids, comparator_isoform_id %in% coding_ids)
c4_coding <- gm_c4 %>%
  filter(reference_isoform_id %in% coding_ids, comparator_isoform_id %in% coding_ids)

# Add original PTC status
ptc_lookup <- ptc_assoc$pair_summary[, c("comparator_isoform_id", "comp_has_ptc")]
c2_coding <- merge(c2_coding, ptc_lookup, by = "comparator_isoform_id", all.x = TRUE)
c2_coding$comp_has_ptc[is.na(c2_coding$comp_has_ptc)] <- FALSE

cat("  Gene-matched coding pairs: C2 =", nrow(c2_coding), " C4 =", nrow(c4_coding), "\n")
cat("  Original PTC+:", sum(c2_coding$comp_has_ptc),
    " PTC-:", sum(!c2_coding$comp_has_ptc), "\n")

# ==============================================================================
# 2. Extract transcript sequences from FASTA
# ==============================================================================
cat("\nExtracting sequences...\n")

all_isoforms <- unique(c(c2_coding$reference_isoform_id,
                          c2_coding$comparator_isoform_id,
                          c4_coding$reference_isoform_id,
                          c4_coding$comparator_isoform_id))

id_file <- tempfile(fileext = ".txt")
writeLines(all_isoforms, id_file)
awk_cmd <- sprintf(
  "awk 'BEGIN{while((getline id < \"%s\") > 0) want[id]=1} /^>/{if(seq && found) print hdr \"\\t\" seq; hdr=substr($1,2); found=(hdr in want); seq=\"\"; next} found{seq=seq$0} END{if(seq && found) print hdr \"\\t\" seq}' '%s'",
  id_file, FASTA_PATH)
t0 <- proc.time()
raw <- system(awk_cmd, intern = TRUE)
cat("  awk completed in", round((proc.time() - t0)[3], 1), "seconds\n")
unlink(id_file)

split_lines <- strsplit(raw, "\t", fixed = TRUE)
seq_vec <- setNames(
  toupper(vapply(split_lines, `[`, character(1), 2)),
  vapply(split_lines, `[`, character(1), 1))
cat("  Sequences loaded:", length(seq_vec), "\n")

# Structures lookup used by Section 4 (transcript-to-genomic mapping for PTC stops)
structs_lookup <- setNames(
  lapply(seq_len(nrow(structures)), function(i) {
    list(starts = structures$exon_starts[[i]],
         ends = structures$exon_ends[[i]],
         strand = structures$strand[i])
  }), structures$isoform_id)

# ==============================================================================
# 3. Trace ref ATG through comparator (Isopair::traceReferenceAtg)
# ==============================================================================
analyze_pairs <- function(pairs, label) {
  cat(sprintf("\n=== Analyzing %s pairs (n=%d) ===\n", label, nrow(pairs)))

  t0 <- proc.time()
  trace <- Isopair::traceReferenceAtg(
    pairs = pairs[, c("reference_isoform_id", "comparator_isoform_id")],
    structures = structures,
    cds_metadata = cds,
    sequences = seq_vec,
    ejc_threshold = 50L,
    resolve_alt_start = FALSE  # legacy behavior; alt-start handled separately if needed
  )
  cat(sprintf("  traceReferenceAtg completed in %.0f seconds\n",
              (proc.time() - t0)[3]))

  # Augment with wrapper-specific columns expected by Section 4 and downstream consumers
  pair_key <- paste(pairs$reference_isoform_id, pairs$comparator_isoform_id, sep = "|")
  trace_key <- paste(trace$reference_isoform_id, trace$comparator_isoform_id, sep = "|")
  trace$gene_id <- pairs$gene_id[match(trace_key, pair_key)]
  trace$strand <- cds$strand[match(trace$comparator_isoform_id, cds$isoform_id)]

  # Backward-compatible column name used by 05_final_report_mashr.Rmd
  trace$comp_orf_length_from_ref_atg <- trace$comp_orf_length

  # Summary
  cat("\n  Results:\n")
  tab <- sort(table(trace$category), decreasing = TRUE)
  for (nm in names(tab)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", nm, tab[nm], 100 * tab[nm] / nrow(trace)))
  }

  trace
}

# Run on C2 (NMD) pairs
c2_results <- analyze_pairs(c2_coding, "NMD (C2)")
c2_results$original_ptc <- c2_coding$comp_has_ptc[
  match(c2_results$comparator_isoform_id, c2_coding$comparator_isoform_id)]

# Run on C4 (Control) pairs
c4_results <- analyze_pairs(c4_coding, "Control (C4)")

# ==============================================================================
# 4. Splice event attribution for effectively-PTC cases
#    (a) Truncated ORF → frameshift / in-frame stop (attribute_ptc_events)
#    (b) Same/longer ORF with downstream EJCs → 3'UTR splice (attribute_3utr_splice)
# ==============================================================================
cat("\n=== Splice Event Attribution (Shared Function) ===\n")

source("analysis_functions.R")

fw_c2 <- readRDS("data_mashr/analysis_cache/fw_c2_allsamples.rds")

# For effectively-PTC cases in C2 (originally PTC-): map stop to genomic coords
eff_ptc_mask <- c2_results$category == "effectively_ptc" & !c2_results$original_ptc
eff_ptc_idx <- which(eff_ptc_mask)
cat("  Reclassified effectively-PTC (originally PTC-):", length(eff_ptc_idx), "\n")

# Map premature stop from transcript-space to genomic coordinates
# (Isopair::transcriptToGenomic — strand-aware)
comp_stop_genomic <- rep(NA_integer_, nrow(c2_results))
for (k in eff_ptc_idx) {
  comp_id <- c2_results$comparator_isoform_id[k]
  strand <- c2_results$strand[k]
  comp_stop_tx <- c2_results$comp_stop_tx_pos[k]
  comp_s <- structs_lookup[[comp_id]]
  if (!is.null(comp_s) && !is.na(comp_stop_tx)) {
    comp_stop_genomic[k] <- Isopair::transcriptToGenomic(
      comp_stop_tx, comp_s$starts, comp_s$ends, strand)
  }
}
c2_results$comp_stop_genomic <- comp_stop_genomic

# Split reclassified pairs into:
#   (a) Truncated ORF: comp_orf < ref_orf → frameshift/in-frame stop attribution
#   (b) Same/longer ORF with downstream EJCs → 3'UTR splice mechanism
eff_ptc_all <- c2_results[eff_ptc_mask, ]
is_truncated <- !is.na(eff_ptc_all$ref_orf_length) &
  eff_ptc_all$comp_orf_length_from_ref_atg < eff_ptc_all$ref_orf_length
is_truncated[is.na(is_truncated)] <- FALSE  # treat NA as non-truncated

cat("  Truncated ORF (frameshift/in-frame stop):", sum(is_truncated), "\n")
cat("  Same/longer ORF (3'UTR splice):", sum(!is_truncated), "\n")

# --- (a) Truncated ORF attribution via shared function ---
trunc_pairs <- eff_ptc_all[is_truncated, ]
attr_mechanism <- rep(NA_character_, nrow(c2_results))
attr_event <- rep(NA_character_, nrow(c2_results))

if (nrow(trunc_pairs) > 0) {
  ptc_stop_vec <- setNames(trunc_pairs$comp_stop_genomic, trunc_pairs$comparator_isoform_id)
  atg_vec <- setNames(trunc_pairs$ref_atg_genomic, trunc_pairs$comparator_isoform_id)
  strand_vec_t <- setNames(trunc_pairs$strand, trunc_pairs$comparator_isoform_id)

  attr_result <- attribute_ptc_events(
    pairs = trunc_pairs[, c("comparator_isoform_id", "reference_isoform_id")],
    fw_events = fw_c2$events,
    profiles = c2_coding,
    ptc_genomic_pos = ptc_stop_vec,
    atg_genomic_pos = atg_vec,
    strand_vec = strand_vec_t,
    is_frameshift_vec = NULL
  )

  attr_lookup <- setNames(seq_len(nrow(attr_result)), attr_result$comparator_isoform_id)
  trunc_idx <- eff_ptc_idx[is_truncated]
  for (k in seq_along(trunc_idx)) {
    comp_id <- c2_results$comparator_isoform_id[trunc_idx[k]]
    j <- attr_lookup[comp_id]
    if (!is.na(j)) {
      attr_mechanism[trunc_idx[k]] <- attr_result$mechanism[j]
      attr_event[trunc_idx[k]] <- attr_result$ptc_causing_event[j]
    }
  }
}

# --- (b) Same/longer ORF attribution: 3'UTR splice ---
nontrunc_pairs <- eff_ptc_all[!is_truncated, ]
nontrunc_idx <- eff_ptc_idx[!is_truncated]

if (nrow(nontrunc_pairs) > 0) {
  utr3_stop_vec <- setNames(nontrunc_pairs$comp_stop_genomic, nontrunc_pairs$comparator_isoform_id)
  utr3_strand_vec <- setNames(nontrunc_pairs$strand, nontrunc_pairs$comparator_isoform_id)

  utr3_result <- attribute_3utr_splice(
    pairs = nontrunc_pairs[, c("comparator_isoform_id", "reference_isoform_id")],
    profiles = c2_coding,
    stop_genomic_pos = utr3_stop_vec,
    strand_vec = utr3_strand_vec
  )

  utr3_lookup <- setNames(seq_len(nrow(utr3_result)), utr3_result$comparator_isoform_id)
  for (k in seq_along(nontrunc_idx)) {
    comp_id <- c2_results$comparator_isoform_id[nontrunc_idx[k]]
    j <- utr3_lookup[comp_id]
    if (!is.na(j)) {
      attr_mechanism[nontrunc_idx[k]] <- utr3_result$mechanism[j]
      attr_event[nontrunc_idx[k]] <- utr3_result$ptc_causing_event[j]
    }
  }
}

c2_results$attr_mechanism <- attr_mechanism
c2_results$attr_event <- attr_event

# Summary of attribution
reclassified <- c2_results[eff_ptc_mask, ]
cat("\n  PTC-causing mechanism for reclassified pairs:\n")
mech_tab <- sort(table(reclassified$attr_mechanism), decreasing = TRUE)
for (nm in names(mech_tab)) {
  cat(sprintf("    %s: %d (%.1f%%)\n", nm, mech_tab[nm],
              100 * mech_tab[nm] / nrow(reclassified)))
}

cat("\n  Event types within Frameshift mechanism:\n")
fs_reclass <- reclassified[reclassified$attr_mechanism == "Frameshift", ]
if (nrow(fs_reclass) > 0) {
  evt_tab <- sort(table(fs_reclass$attr_event), decreasing = TRUE)
  for (nm in names(evt_tab)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", nm, evt_tab[nm],
                100 * evt_tab[nm] / nrow(fs_reclass)))
  }
}

cat("\n  Event types within In-frame stop mechanism:\n")
ifs_reclass <- reclassified[reclassified$attr_mechanism == "In-frame stop", ]
if (nrow(ifs_reclass) > 0) {
  evt_tab2 <- sort(table(ifs_reclass$attr_event), decreasing = TRUE)
  for (nm in names(evt_tab2)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", nm, evt_tab2[nm],
                100 * evt_tab2[nm] / nrow(ifs_reclass)))
  }
}

cat("\n  Event types within 3'UTR splice mechanism:\n")
utr3_reclass <- reclassified[reclassified$attr_mechanism == "3'UTR splice", ]
if (nrow(utr3_reclass) > 0) {
  evt_tab3 <- sort(table(utr3_reclass$attr_event), decreasing = TRUE)
  for (nm in names(evt_tab3)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", nm, evt_tab3[nm],
                100 * evt_tab3[nm] / nrow(utr3_reclass)))
  }
}

# ==============================================================================
# 5. Summary and save
# ==============================================================================
cat("\n=== FINAL SUMMARY ===\n\n")

# Three groups
group2 <- sum(c2_results$ref_atg_exonic_in_comp == TRUE, na.rm = TRUE)
group3 <- sum(c2_results$ref_atg_exonic_in_comp == FALSE, na.rm = TRUE)
cat("Three groups:\n")
cat("  Group 2 (ref ATG available):", group2, sprintf("(%.1f%%)\n", 100 * group2 / nrow(c2_results)))
cat("  Group 3 (ref ATG lost):", group3, sprintf("(%.1f%%)\n", 100 * group3 / nrow(c2_results)))

# Reclassification
cat("\nAmong originally PTC- pairs (n =", sum(!c2_results$original_ptc), "):\n")
ptc_neg <- c2_results[!c2_results$original_ptc, ]
tab <- sort(table(ptc_neg$category), decreasing = TRUE)
for (nm in names(tab)) {
  cat(sprintf("  %s: %d (%.1f%%)\n", nm, tab[nm], 100 * tab[nm] / nrow(ptc_neg)))
}

# Updated accounting
orig_ptc_pos <- sum(c2_results$original_ptc)
reclass_ptc <- sum(c2_results$category == "effectively_ptc" & !c2_results$original_ptc, na.rm = TRUE)
total_ptc <- orig_ptc_pos + reclass_ptc
cat(sprintf("\nUpdated PTC accounting:\n"))
cat(sprintf("  Original PTC+: %d\n", orig_ptc_pos))
cat(sprintf("  Reclassified effectively-PTC: %d\n", reclass_ptc))
cat(sprintf("  Total PTC-mediated: %d of %d (%.1f%%)\n",
            total_ptc, nrow(c2_results), 100 * total_ptc / nrow(c2_results)))

results <- list(
  c2 = c2_results,
  c4 = c4_results,
  summary = list(
    n_c2 = nrow(c2_results),
    n_c4 = nrow(c4_results),
    group2_n = group2,
    group3_n = group3,
    original_ptc_pos = orig_ptc_pos,
    reclassified_ptc = reclass_ptc,
    total_ptc = total_ptc
  ),
  metadata = list(
    run_timestamp = Sys.time(),
    fasta_path = FASTA_PATH,
    isopair_version = as.character(packageVersion("Isopair")),
    backend = "Isopair::traceReferenceAtg"
  )
)

saveRDS(results, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
