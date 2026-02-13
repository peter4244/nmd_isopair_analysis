#!/usr/bin/env Rscript
# Regenerate exon structures for union model filtered isoforms

library(tidyverse)
library(rtracklayer)

output_dir <- "results/isoform_transitions/v3.0_reference_based"

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   REGENERATE EXON STRUCTURES FOR UNION MODEL                 ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Load filtered isoform list
cat("Loading filtered isoforms...\n")
isoforms_to_process <- readRDS(file.path(output_dir, "isoforms_for_union_model.rds"))

isoform_ids <- unique(isoforms_to_process$isoform_id)
cat("Isoforms to extract:", length(isoform_ids), "\n")

# Identify ENST vs PB IDs
enst_ids <- isoform_ids[grepl("^ENST", isoform_ids)]
pb_ids <- isoform_ids[grepl("^PB\\.", isoform_ids)]

cat("  ENST (annotated):", length(enst_ids), "\n")
cat("  PB (novel):", length(pb_ids), "\n\n")

# === Extract from GENCODE GFF3 (Comprehensive) ===
cat("═══ Loading GENCODE GFF3 (Comprehensive) ═══\n")
gencode_gff <- import("reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz")

gencode_exons <- as.data.frame(gencode_gff) %>%
  filter(type == "exon") %>%
  filter(transcript_id %in% enst_ids) %>%
  select(
    isoform_id = transcript_id,
    gene_id,
    seqnames,
    start,
    end,
    strand
  )

cat("  GENCODE exons extracted:", nrow(gencode_exons), "\n")
cat("  GENCODE isoforms matched:", length(unique(gencode_exons$isoform_id)), "/", length(enst_ids), "\n\n")

# === Extract from PacBio GFF ===
cat("═══ Loading PacBio GFF ═══\n")
pb_gff <- import("isoseq/collapse/merge-collapsed.gff")

pb_exons <- as.data.frame(pb_gff) %>%
  filter(type == "exon") %>%
  filter(transcript_id %in% pb_ids) %>%
  select(
    isoform_id = transcript_id,
    gene_id = any_of(c("gene_id", "gene")),
    seqnames,
    start,
    end,
    strand
  )

cat("  PacBio exons extracted:", nrow(pb_exons), "\n")
cat("  PacBio isoforms matched:", length(unique(pb_exons$isoform_id)), "/", length(pb_ids), "\n\n")

# === Combine both sources ===
cat("═══ Combining exon structures ═══\n")
all_exons <- bind_rows(gencode_exons, pb_exons)

cat("  Total exons:", nrow(all_exons), "\n")
cat("  Total isoforms:", length(unique(all_exons$isoform_id)), "\n\n")

# === Build exon structures ===
cat("═══ Building exon structures ═══\n")
exon_structures <- all_exons %>%
  group_by(isoform_id, gene_id, seqnames, strand) %>%
  arrange(isoform_id, start) %>%
  summarize(
    exon_starts = list(start),
    exon_ends = list(end),
    n_exons = n(),
    .groups = "drop"
  )

cat("  Exon structures built:", nrow(exon_structures), "isoforms\n\n")

# === Infer junctions ===
cat("═══ Inferring junctions ═══\n")
exon_structures <- exon_structures %>%
  mutate(
    junctions = map2(exon_ends, exon_starts, function(ends, starts) {
      if (length(ends) <= 1) {
        return(tibble(junction_start = integer(0), junction_end = integer(0)))
      }
      tibble(
        junction_start = ends[-length(ends)] + 1,
        junction_end = starts[-1] - 1
      )
    }),
    n_junctions = map_int(junctions, nrow)
  )

cat("  Junctions inferred\n")
cat("  Isoforms with ≥1 junction:", sum(exon_structures$n_junctions > 0), "\n\n")

# === Summary statistics ===
cat("═══ Summary Statistics ═══\n")
cat("  Exons per isoform:\n")
cat("    Mean:", sprintf("%.1f", mean(exon_structures$n_exons)), "\n")
cat("    Median:", median(exon_structures$n_exons), "\n")
cat("    Range:", min(exon_structures$n_exons), "-", max(exon_structures$n_exons), "\n\n")

cat("  Junctions per isoform:\n")
cat("    Mean:", sprintf("%.1f", mean(exon_structures$n_junctions)), "\n")
cat("    Median:", median(exon_structures$n_junctions), "\n\n")

# === Coverage check ===
coverage <- length(unique(exon_structures$isoform_id)) / length(isoform_ids)
cat("═══ Coverage Check ═══\n")
cat("  Requested isoforms:", length(isoform_ids), "\n")
cat("  Extracted isoforms:", length(unique(exon_structures$isoform_id)), "\n")
cat("  Coverage:", sprintf("%.2f%%", 100 * coverage), "\n\n")

if (coverage >= 0.9999) {
  cat("✅ 100% coverage achieved!\n\n")
} else {
  missing_count <- length(isoform_ids) - length(unique(exon_structures$isoform_id))
  cat("⚠️  Coverage not yet 100%\n")
  cat("  Missing:", missing_count, "isoforms\n\n")

  # Show some missing isoforms
  missing_ids <- setdiff(isoform_ids, unique(exon_structures$isoform_id))
  cat("  Examples of missing isoforms:\n")
  print(head(missing_ids, 10))
  cat("\n")
}

# === Save ===
cat("═══ Saving Results ═══\n")
saveRDS(exon_structures,
        file.path(output_dir, "exon_structures_by_isoform_full.rds"))

cat("  Saved: exon_structures_by_isoform_full.rds\n\n")

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   EXON STRUCTURE REGENERATION COMPLETE                       ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("Next step: Update build_union_exon_model_full.R to use new exon structures file\n\n")
