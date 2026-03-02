#!/usr/bin/env Rscript
# 01_compute_ptc_status.R
#
# Compute PTC (Premature Termination Codon) status for all coding isoforms.
#
# PTC definition: stop codon >50nt upstream of the last exon-exon junction
# (the canonical NMD "50-nucleotide rule")
#
# Outputs:
#   ptc_status.rds              — PTC features for all coding isoforms
#   gencode_validation.tsv      — cross-tab: computed PTC vs GENCODE NMD biotype
#
# Run from: results/isoform_transitions/Version_6.0/

suppressPackageStartupMessages(library(tidyverse))

# ─── Parse arguments ─────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
source_type <- "oarfish"  # default
if ("--source" %in% args) {
  idx <- which(args == "--source")
  if (idx < length(args)) source_type <- args[idx + 1]
}
stopifnot(source_type %in% c("oarfish", "isocall"))

cat("\n")
cat("==================================================================\n")
cat(sprintf("   01: Compute PTC Status for All Coding Isoforms (%s)\n", source_type))
cat("==================================================================\n\n")

# ─── 1. Load data ────────────────────────────────────────────────────────────

data_dir <- if (source_type == "isocall") "data/isocall/" else "data/"

cat(sprintf("Loading data from %s...\n", data_dir))
cds_meta <- readRDS(file.path(data_dir, "isoform_cds_metadata.rds"))
structs  <- readRDS(file.path(data_dir, "isoform_structures.rds"))

coding <- cds_meta %>%
  filter(coding_status == "coding") %>%
  inner_join(structs, by = "isoform_id")

cat(sprintf("  %d coding isoforms with exon structures\n", nrow(coding)))

# ─── 2. Compute PTC features ────────────────────────────────────────────────

cat("Computing PTC features...\n")

# Returns: ptc_distance, n_downstream_ejcs, stop_in_last_exon
#
# Coordinate conventions:
#   cds_start/cds_stop: always min/max genomic coordinate
#   + strand: stop codon biological position = cds_stop
#   - strand: stop codon biological position = cds_start
#   Exons ordered 5' to 3' biologically for mRNA-space calculations

compute_ptc_features <- function(cds_start, cds_stop, strand,
                                  exon_starts, exon_ends, n_exons) {
  result <- list(ptc_distance = NA_real_, n_downstream_ejcs = NA_integer_,
                 stop_in_last_exon = NA)
  if (n_exons == 1) return(result)  # No EJCs in single-exon transcripts

  # Biological stop codon position (genomic)
  stop_pos <- if (strand == "+") cds_stop else cds_start

  # Order exons 5' → 3'
  if (strand == "+") {
    ord <- order(exon_starts)
  } else {
    ord <- order(exon_starts, decreasing = TRUE)
  }
  e_starts <- exon_starts[ord]
  e_ends <- exon_ends[ord]

  # Exon lengths and cumulative mRNA positions
  exon_lengths <- e_ends - e_starts
  cum_starts <- c(0, cumsum(exon_lengths[-length(exon_lengths)]))

  # Last EJC in mRNA space = cumulative length through penultimate exon
  last_ejc_mRNA <- sum(exon_lengths[-length(exon_lengths)])

  # Find which exon contains the stop codon
  exon_idx <- which(e_starts <= stop_pos & stop_pos <= e_ends)
  if (length(exon_idx) == 0) return(result)
  exon_idx <- exon_idx[1]

  # Stop codon position in mRNA coordinates
  if (strand == "+") {
    offset_in_exon <- stop_pos - e_starts[exon_idx]
  } else {
    offset_in_exon <- e_ends[exon_idx] - stop_pos
  }
  stop_mRNA <- cum_starts[exon_idx] + offset_in_exon

  result$ptc_distance <- last_ejc_mRNA - stop_mRNA
  result$n_downstream_ejcs <- as.integer(n_exons - exon_idx)
  result$stop_in_last_exon <- (exon_idx == n_exons)
  result
}

ptc_results <- coding %>%
  rowwise() %>%
  mutate(
    features = list(compute_ptc_features(
      cds_start, cds_stop, as.character(strand),
      exon_starts, exon_ends, n_exons
    )),
    ptc_distance = features$ptc_distance,
    n_downstream_ejcs = features$n_downstream_ejcs,
    stop_in_last_exon = features$stop_in_last_exon
  ) %>%
  ungroup() %>%
  select(-features) %>%
  mutate(
    has_ptc = case_when(
      n_exons == 1 ~ FALSE,     # Single-exon: no EJCs possible
      is.na(ptc_distance) ~ NA, # Stop codon not found in exons
      ptc_distance > 50 ~ TRUE, # PTC: >50nt from last EJC
      TRUE ~ FALSE              # Not PTC
    ),
    ptc_distance_bin = case_when(
      n_exons == 1 ~ "single_exon",
      is.na(ptc_distance) ~ "unknown",
      ptc_distance <= 50 ~ "0-50 (no PTC)",
      ptc_distance <= 100 ~ "51-100",
      ptc_distance <= 200 ~ "101-200",
      ptc_distance <= 500 ~ "201-500",
      TRUE ~ ">500"
    ),
    is_gencode = startsWith(isoform_id, "ENST"),
    is_novel = if (source_type == "isocall") grepl("\\.novel\\d+$", isoform_id) else startsWith(isoform_id, "PB.")
  )

# ─── 3. Summary statistics ──────────────────────────────────────────────────

cat(sprintf("\n  Single-exon (no EJCs): %d\n",
            sum(ptc_results$n_exons == 1, na.rm = TRUE)))
cat(sprintf("  Stop not in exon (NA): %d\n",
            sum(is.na(ptc_results$ptc_distance) & ptc_results$n_exons > 1)))
cat(sprintf("  PTC (>50nt): %d\n", sum(ptc_results$has_ptc == TRUE, na.rm = TRUE)))
cat(sprintf("  Not PTC (<=50nt): %d\n",
            sum(ptc_results$has_ptc == FALSE & ptc_results$n_exons > 1, na.rm = TRUE)))
novel_label <- if (source_type == "isocall") "novel" else "PacBio novel"
cat(sprintf("  GENCODE: %d, %s: %d\n",
            sum(ptc_results$is_gencode), novel_label, sum(ptc_results$is_novel)))

# Distance distribution for multi-exon coding isoforms
multi_exon <- ptc_results %>% filter(n_exons > 1, !is.na(ptc_distance))
cat(sprintf("\nPTC distance distribution (n = %d multi-exon coding isoforms):\n",
            nrow(multi_exon)))
cat(sprintf("  Median: %.0f\n", median(multi_exon$ptc_distance)))
cat(sprintf("  Mean: %.0f\n", mean(multi_exon$ptc_distance)))
q <- quantile(multi_exon$ptc_distance, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1))
cat("  Quantiles:\n")
print(q)

# Stratify by source
cat("\nPTC rates by isoform source:\n")
ptc_results %>%
  filter(!is.na(has_ptc)) %>%
  mutate(source = case_when(is_gencode ~ "GENCODE", is_novel ~ novel_label, TRUE ~ "other")) %>%
  group_by(source) %>%
  summarise(
    n = n(),
    n_ptc = sum(has_ptc),
    ptc_rate = round(mean(has_ptc), 3),
    median_distance = median(ptc_distance[!is.na(ptc_distance)]),
    .groups = "drop"
  ) %>%
  print()

# ─── 4. GENCODE NMD biotype validation ──────────────────────────────────────

if (source_type == "isocall") {
  cat("\n  Note: Skipping GENCODE NMD biotype validation for isocall\n")
  cat("  (isocall DE files lack biotype column; algorithm already validated by oarfish run)\n")
  gencode_validation <- NULL
} else {
  cat("\n\n")
  cat("==================================================================\n")
  cat("   GENCODE Validation: Computed PTC vs NMD Biotype\n")
  cat("==================================================================\n\n")

  # Load biotype info from DE files (GENCODE isoforms only)
  de_dir <- "/Users/petecastaldi/claude_projects/nmd/longread_dge"
  de_date <- "2026.1.18"

  # Use any one DE file to get biotype annotations (they're the same across cell types)
  de_sample <- read_csv(file.path(de_dir, paste0("nmd_dge_at2_", de_date, ".csv")),
                         show_col_types = FALSE) %>%
    select(txid, biotype) %>%
    rename(isoform_id = txid) %>%
    distinct()

  validation <- ptc_results %>%
    filter(is_gencode, !is.na(has_ptc)) %>%
    inner_join(de_sample, by = "isoform_id") %>%
    mutate(nmd_biotype = biotype == "nonsense_mediated_decay")

  cat("Cross-tabulation: Computed PTC × GENCODE NMD biotype\n\n")
  tab <- table(Computed_PTC = validation$has_ptc,
               GENCODE_NMD_biotype = validation$nmd_biotype)
  print(tab)

  # Concordance
  concordant <- sum(validation$has_ptc & validation$nmd_biotype) +
    sum(!validation$has_ptc & !validation$nmd_biotype)
  cat(sprintf("\nConcordance: %d / %d (%.1f%%)\n",
              concordant, nrow(validation), 100 * concordant / nrow(validation)))

  # Sensitivity/specificity of our PTC call vs GENCODE annotation
  tp <- sum(validation$has_ptc & validation$nmd_biotype)
  fp <- sum(validation$has_ptc & !validation$nmd_biotype)
  fn <- sum(!validation$has_ptc & validation$nmd_biotype)
  tn <- sum(!validation$has_ptc & !validation$nmd_biotype)
  cat(sprintf("Sensitivity (recall): %.1f%% (%d/%d)\n", 100 * tp / (tp + fn), tp, tp + fn))
  cat(sprintf("Specificity: %.1f%% (%d/%d)\n", 100 * tn / (tn + fp), tn, tn + fp))
  cat(sprintf("Precision (PPV): %.1f%% (%d/%d)\n", 100 * tp / (tp + fp), tp, tp + fp))

  gencode_validation <- tibble(
    computed_ptc = c(TRUE, TRUE, FALSE, FALSE),
    gencode_nmd_biotype = c(TRUE, FALSE, TRUE, FALSE),
    count = c(tp, fp, fn, tn)
  ) %>%
    mutate(
      concordant = computed_ptc == gencode_nmd_biotype,
      pct = round(100 * count / sum(count), 2)
    )
}

# ─── 5. Save outputs ────────────────────────────────────────────────────────

output_dir <- if (source_type == "isocall") "results/ptc/results/isocall" else "results/ptc/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save the full PTC status table
ptc_status <- ptc_results %>%
  select(isoform_id, gene_id, n_exons, orf_length,
         ptc_distance, has_ptc, ptc_distance_bin,
         n_downstream_ejcs, stop_in_last_exon,
         is_gencode, is_novel)

saveRDS(ptc_status, file.path(output_dir, "ptc_status.rds"))
cat(sprintf("\nSaved: %s (%d isoforms)\n",
            file.path(output_dir, "ptc_status.rds"), nrow(ptc_status)))

if (!is.null(gencode_validation)) {
  write_tsv(gencode_validation, file.path(output_dir, "gencode_validation.tsv"))
  cat(sprintf("Saved: %s\n", file.path(output_dir, "gencode_validation.tsv")))
}

cat("\nDone.\n")
