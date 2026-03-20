#!/usr/bin/env Rscript
# ==============================================================================
# 05r_ref_atg_analysis.R
#
# Reference-ATG-anchored NMD analysis.
#
# For each gene-matched NMD pair (C2) and Control pair (C4):
#   1. Check if the reference's CDS start ATG is exonic in the comparator
#   2. If yes, trace the ORF from that ATG through the comparator's transcript
#   3. Compare ORF length to the reference, count downstream EJCs
#   4. Reclassify: is this effectively a PTC?
#   5. Attribute splice events for the effectively-PTC cases
#
# Three groups:
#   Group 1: Reference isoforms (non-NMD baseline)
#   Group 2: NMD+ / ref ATG available (ref CDS start is exonic in comparator)
#   Group 3: NMD+ / ref ATG not available (ref CDS start lost)
#
# Inputs:
#   - data_mashr/cds.rds, structures.rds, profiles_c2/c4_allsamples.rds
#   - data_mashr/analysis_cache/ptc_c2_allsamples.rds
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
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

FASTA_PATH <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta"
OUTPUT_PATH <- "data_mashr/analysis_cache/ref_atg_analysis.rds"

cat("=== Reference-ATG-Anchored NMD Analysis ===\n")
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
# 2. Build exon lookup and helper functions
# ==============================================================================
cat("\nBuilding exon lookups...\n")

structs_lookup <- setNames(
  lapply(seq_len(nrow(structures)), function(i) {
    list(starts = structures$exon_starts[[i]],
         ends = structures$exon_ends[[i]],
         strand = structures$strand[i])
  }), structures$isoform_id)

# Check if a genomic position (and the full 3-nt codon) is exonic
is_codon_exonic <- function(pos, iso_id, strand) {
  ex <- structs_lookup[[iso_id]]
  if (is.null(ex) || is.na(pos)) return(FALSE)
  # For + strand: codon is at pos, pos+1, pos+2
  # For - strand: codon is at pos, pos-1, pos-2
  if (strand == "+") {
    positions <- c(pos, pos + 1L, pos + 2L)
  } else {
    positions <- c(pos, pos - 1L, pos - 2L)
  }
  all(sapply(positions, function(p) any(p >= ex$starts & p <= ex$ends)))
}

# Map genomic position to transcript-space (1-based)
genomic_to_transcript <- function(genomic_pos, exon_starts, exon_ends, strand) {
  if (strand == "-") {
    exon_starts <- rev(exon_starts)
    exon_ends <- rev(exon_ends)
  }
  cum <- 0L
  for (j in seq_along(exon_starts)) {
    es <- exon_starts[j]; ee <- exon_ends[j]
    if (strand == "+") {
      if (genomic_pos >= es && genomic_pos <= ee)
        return(cum + (genomic_pos - es) + 1L)
    } else {
      if (genomic_pos >= es && genomic_pos <= ee)
        return(cum + (ee - genomic_pos) + 1L)
    }
    cum <- cum + (ee - es + 1L)
  }
  NA_integer_
}

# Compute junction positions in transcript space
get_junctions <- function(exon_starts, exon_ends, strand) {
  n <- length(exon_starts)
  if (n <= 1) return(integer(0))
  exon_lengths <- exon_ends - exon_starts + 1L
  if (strand == "-") exon_lengths <- rev(exon_lengths)
  cumsum(exon_lengths)[-n]
}

stop_codons <- c("TAA", "TAG", "TGA")

# ==============================================================================
# 3. Extract transcript sequences
# ==============================================================================
cat("Extracting sequences...\n")

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

# ==============================================================================
# 4. Core analysis: trace ref ATG through comparator
# ==============================================================================
analyze_pairs <- function(pairs, label) {
  cat(sprintf("\n=== Analyzing %s pairs (n=%d) ===\n", label, nrow(pairs)))

  ref_cds_data <- cds[match(pairs$reference_isoform_id, cds$isoform_id), ]
  ref_cds_data$cds_5prime <- ifelse(ref_cds_data$strand == "+",
                                     ref_cds_data$cds_start, ref_cds_data$cds_stop)

  results <- data.frame(
    reference_isoform_id = pairs$reference_isoform_id,
    comparator_isoform_id = pairs$comparator_isoform_id,
    gene_id = pairs$gene_id,
    ref_atg_genomic = ref_cds_data$cds_5prime,
    strand = ref_cds_data$strand,
    ref_atg_exonic_in_comp = NA,
    ref_orf_length = NA_integer_,
    comp_orf_length_from_ref_atg = NA_integer_,
    comp_stop_tx_pos = NA_integer_,
    n_downstream_ejc = NA_integer_,
    category = NA_character_,
    stringsAsFactors = FALSE
  )

  t0 <- proc.time()

  for (i in seq_len(nrow(pairs))) {
    ref_id <- pairs$reference_isoform_id[i]
    comp_id <- pairs$comparator_isoform_id[i]
    ref_atg <- ref_cds_data$cds_5prime[i]
    strand <- ref_cds_data$strand[i]

    if (is.na(ref_atg)) { results$category[i] <- "no_ref_cds"; next }

    # Check if ref ATG codon is exonic in comparator
    atg_exonic <- is_codon_exonic(ref_atg, comp_id, strand)
    results$ref_atg_exonic_in_comp[i] <- atg_exonic

    if (!atg_exonic) { results$category[i] <- "ref_atg_lost"; next }
    if (!ref_id %in% names(seq_vec) || !comp_id %in% names(seq_vec)) {
      results$category[i] <- "no_sequence"; next
    }

    ref_s <- structs_lookup[[ref_id]]
    comp_s <- structs_lookup[[comp_id]]

    # Map ref ATG to transcript positions
    ref_tx_pos <- genomic_to_transcript(ref_atg, ref_s$starts, ref_s$ends, strand)
    comp_tx_pos <- genomic_to_transcript(ref_atg, comp_s$starts, comp_s$ends, strand)

    if (is.na(ref_tx_pos) || is.na(comp_tx_pos)) {
      results$category[i] <- "mapping_failed"; next
    }

    ref_seq <- seq_vec[ref_id]
    comp_seq <- seq_vec[comp_id]

    # Verify ATG in both
    if (substr(ref_seq, ref_tx_pos, ref_tx_pos + 2) != "ATG") {
      results$category[i] <- "not_atg_in_ref"; next
    }
    if (substr(comp_seq, comp_tx_pos, comp_tx_pos + 2) != "ATG") {
      results$category[i] <- "not_atg_in_comp"; next
    }

    # Find ORF from ref ATG in reference
    ref_stop <- NA_integer_
    pos <- ref_tx_pos + 3L
    while (pos + 2L <= nchar(ref_seq)) {
      if (substr(ref_seq, pos, pos + 2L) %in% stop_codons) {
        ref_stop <- pos; break
      }
      pos <- pos + 3L
    }
    if (!is.na(ref_stop)) results$ref_orf_length[i] <- ref_stop - ref_tx_pos

    # Find ORF from ref ATG in comparator
    comp_stop <- NA_integer_
    pos <- comp_tx_pos + 3L
    while (pos + 2L <= nchar(comp_seq)) {
      if (substr(comp_seq, pos, pos + 2L) %in% stop_codons) {
        comp_stop <- pos; break
      }
      pos <- pos + 3L
    }

    if (is.na(comp_stop)) { results$category[i] <- "no_stop_in_comp"; next }

    comp_orf_len <- comp_stop - comp_tx_pos
    results$comp_orf_length_from_ref_atg[i] <- comp_orf_len
    results$comp_stop_tx_pos[i] <- comp_stop

    # Count downstream EJCs
    junctions <- get_junctions(comp_s$starts, comp_s$ends, strand)
    stop_end <- comp_stop + 2L
    n_dejc <- sum(junctions > stop_end)
    results$n_downstream_ejc[i] <- n_dejc

    # Classify
    ref_orf_len <- results$ref_orf_length[i]
    if (!is.na(ref_orf_len) && comp_orf_len < ref_orf_len && n_dejc > 0) {
      results$category[i] <- "effectively_ptc"
    } else if (!is.na(ref_orf_len) && comp_orf_len < ref_orf_len && n_dejc == 0) {
      results$category[i] <- "truncated_no_ejc"
    } else if (n_dejc > 0) {
      results$category[i] <- "same_or_longer_with_ejc"
    } else {
      results$category[i] <- "no_downstream_ejc"
    }

    if (i %% 500 == 0) cat(sprintf("  %d / %d (%.0f sec)\n", i, nrow(pairs),
                                    (proc.time() - t0)[3]))
  }

  cat(sprintf("  Completed in %.0f seconds\n", (proc.time() - t0)[3]))

  # Summary
  cat("\n  Results:\n")
  tab <- sort(table(results$category), decreasing = TRUE)
  for (nm in names(tab)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", nm, tab[nm], 100 * tab[nm] / nrow(pairs)))
  }

  results
}

# Run on C2 (NMD) pairs
c2_results <- analyze_pairs(c2_coding, "NMD (C2)")

# Add original PTC status
c2_results$original_ptc <- c2_coding$comp_has_ptc

# Run on C4 (Control) pairs
c4_results <- analyze_pairs(c4_coding, "Control (C4)")

# ==============================================================================
# 5. Splice event attribution for effectively-PTC cases
# ==============================================================================
cat("\n=== Splice Event Attribution ===\n")

# For effectively-PTC cases in C2: find the splice event responsible
eff_ptc_mask <- c2_results$category == "effectively_ptc" & !c2_results$original_ptc
eff_ptc_idx <- which(eff_ptc_mask)
cat("  Reclassified effectively-PTC (originally PTC-):", length(eff_ptc_idx), "\n")

attr_event <- character(nrow(c2_results))
attr_event[] <- NA_character_

for (k in eff_ptc_idx) {
  ref_id <- c2_results$reference_isoform_id[k]
  comp_id <- c2_results$comparator_isoform_id[k]

  # Find this pair in profiles
  pi <- which(c2_coding$reference_isoform_id == ref_id &
                c2_coding$comparator_isoform_id == comp_id)
  if (length(pi) == 0) next

  de <- c2_coding$detailed_events[[pi[1]]]
  if (is.null(de) || nrow(de) == 0) { attr_event[k] <- "no_events"; next }

  ref_atg <- c2_results$ref_atg_genomic[k]
  # Get comp stop in genomic coords: map back from transcript
  comp_s <- structs_lookup[[comp_id]]
  strand <- c2_results$strand[k]

  # Find splice event between ref ATG and premature stop that caused the truncation
  # Filter events to those whose genomic coordinates fall between ATG and stop
  ref_atg_g <- c2_results$ref_atg_genomic[k]
  # Get the premature stop codon's genomic position
  # We need to map comp_stop_tx_pos back to genomic — use the event genomic coords instead
  # Filter events by: event boundaries within the ATG-to-transcript-end region
  # Since we can't easily map stop back to genomic, use all events and prioritize
  # frameshifts (the most likely cause of ORF truncation)

  # Prioritize: frameshift events first, then in-frame events
  found <- FALSE
  # Sort events by proximity to ref ATG for consistent attribution
  ev_pos <- pmin(de$five_prime, de$three_prime)
  ev_order <- order(abs(ev_pos - ref_atg_g))

  for (j in ev_order) {
    bp <- de$bp_diff[j]
    if (!is.na(bp) && bp %% 3 != 0) {
      attr_event[k] <- de$event_type[j]
      found <- TRUE
      break
    }
  }
  if (!found) {
    for (j in ev_order) {
      bp <- de$bp_diff[j]
      if (!is.na(bp) && bp %% 3 == 0 && bp != 0) {
        attr_event[k] <- paste0(de$event_type[j], " (in-frame)")
        found <- TRUE
        break
      }
    }
  }
  if (!found) attr_event[k] <- "unresolved"
}

c2_results$attr_event <- attr_event

# Summary of attribution
reclassified <- c2_results[eff_ptc_mask, ]
cat("\n  Splice events for reclassified effectively-PTC:\n")
evt_tab <- sort(table(reclassified$attr_event), decreasing = TRUE)
for (nm in names(evt_tab)) {
  cat(sprintf("    %s: %d (%.1f%%)\n", nm, evt_tab[nm], 100 * evt_tab[nm] / nrow(reclassified)))
}

# ==============================================================================
# 6. Summary
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

# ==============================================================================
# 7. Save
# ==============================================================================
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
    fasta_path = FASTA_PATH
  )
)

saveRDS(results, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
