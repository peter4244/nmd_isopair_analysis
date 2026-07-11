# SF30 / SF31 data_export.R — canonical Isopair-side derivation
#
# Reproduces the population construction for the two figures:
#   SF30 — Stop-to-last-EJC distance dose response
#   SF31 — NMD effect size by downstream EJC count
#
# Population = Subset 2 (n = 1,166) NMD arm — the ref-AUG-traceable,
# ENST-reference, gene-matched pairwise scope defined by SF28. Consumed
# directly from the canonical Subset 2 export
# (`multipanel/figure5_dl_model/data/subset2_n1166_isoforms.tsv`) so both
# SFs hang off Subset 2 in the SF28 flowchart. Single source of truth.
#
# TD2 bias fix (2026-07-09):
#   The prior version anchored on ptc.rds$ptc_distance (distance from the
#   TD2-called stop to the last EJC). TD2 systematically avoids PTC-containing
#   frames, so for occult-PTC comparators the stop is placed on a non-PTC ORF
#   and ptc_distance is attenuated toward zero. Direction of artifact =
#   attenuated dose response.
#
#   Fix: anchor on the reference-AUG-projected stop
#   (`ref_atg_analysis.rds$c2$comp_stop_tx_pos`), which is unbiased. Compute
#   `stop_to_last_ejc_via_ref_nt = last_ejc_tx_pos - comp_stop_tx_pos`, where
#   `last_ejc_tx_pos` comes from the comparator's exon structure. Downstream
#   EJC count is read directly from `ref_atg_analysis.rds$c2$n_downstream_ejc`
#   (also unbiased — computed against the ref-AUG-projected stop).
#
# Output: two TSVs, one per SF dir.

suppressPackageStartupMessages({
  library(data.table)
})

# ── Provenance paths ────────────────────────────────────────────────────
ISOPAIR_ROOT <- "~/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper"
NMD_ROOT     <- "~/claude_projects/nmd"

REF_ATG_RDS    <- file.path(ISOPAIR_ROOT, "data_mashr", "analysis_cache", "ref_atg_analysis.rds")
STRUCTURES_RDS <- file.path(ISOPAIR_ROOT, "data_mashr", "structures.rds")
MASHR_RDS      <- file.path(NMD_ROOT, "isocall_dge/mashr", "mashr_isoform_model_2026.3.10.rds")
SUBSET2_TSV    <- file.path(NMD_ROOT, "figures/multipanel/figure5_dl_model/data",
                             "subset2_n1166_isoforms.tsv")

SF30_DIR <- file.path("~/claude_projects/nmd/figures/SupplementalFigures",
                       "SF30_PTCDistanceDoseResponse", "data")
SF31_DIR <- file.path("~/claude_projects/nmd/figures/SupplementalFigures",
                       "SF31_NMDEffectByEJCCount", "data")

# ── Subset 2 NMD arm (n = 819 4-CT re-scope; was 1,166) — the canonical population ──
message("[data_export] reading Subset 2 export ...")
subset2 <- fread(SUBSET2_TSV)
subset2_nmd_ids <- subset2[arm == "NMD", comparator_isoform_id]
stopifnot(length(subset2_nmd_ids) == 819L)

# NMD-susceptible-in->=1-CT verification: every SF30 population comparator must be
# classified NMD-susceptible in at least one of the four manuscript cell types
# (the all_samples NMD set is the union over AT/DD/FB/MV) -- confirms no comparator
# is carried on a non-manuscript cell type.
.nmd_class <- readRDS(file.path(ISOPAIR_ROOT, "data_mashr", "nmd_classification.rds"))
.nmd_union_4ct <- unique(unlist(lapply(c("AT", "DD", "FB", "MV"),
                                       function(ct) .nmd_class[[ct]]$nmd)))
stopifnot(all(subset2_nmd_ids %in% .nmd_union_4ct))
message(sprintf("[data_export] NMD-susceptible-in->=1-CT: %d/%d comparators confirmed",
                sum(subset2_nmd_ids %in% .nmd_union_4ct), length(subset2_nmd_ids)))

# ── Load source objects ────────────────────────────────────────────────
message("[data_export] reading ref_atg_analysis.rds ...")
ref_atg    <- readRDS(REF_ATG_RDS)
ref_atg_c2 <- as.data.table(ref_atg$c2)[comparator_isoform_id %in% subset2_nmd_ids]
stopifnot(nrow(ref_atg_c2) == 819L)

message("[data_export] reading structures.rds ...")
structures <- as.data.table(readRDS(STRUCTURES_RDS))

message("[data_export] reading mashr_isoform_model_2026.3.10.rds ...")
mashr_model <- readRDS(MASHR_RDS)
pm <- as.data.frame(mashr_model$result$PosteriorMean)
smg1i_cols <- grep("^Smg1i_in_", colnames(pm), value = TRUE)
# 4-CT guard (m1): the cross-cell-type mashr posterior-mean averaging must be
# over exactly the four manuscript cell types (Smg1i_in_{AT,DD,FB,MV}); assert no
# DD_ALI/DO columns leak into the average.
stopifnot(length(smg1i_cols) == 4,
          !any(grepl("DD_ALI|DO_ALI|_DO_|_DO$", smg1i_cols)))
meta_logfc <- data.table(
  txid       = rownames(pm),
  mean_logFC = rowMeans(pm[, smg1i_cols, drop = FALSE], na.rm = TRUE)
)

# ── Compute last-EJC transcript position per Subset 2 NMD comparator ───
# For a transcript with n_exons ≥ 2, the last EJC sits at the 3′ boundary of
# the second-to-last exon (in transcript order). Equivalently:
#   last_ejc_tx_pos = tx_length - width_of_last_transcript_order_exon
# For + strand, last transcript-order exon = last genomic-order exon (n_exons).
# For − strand, last transcript-order exon = first genomic-order exon (1).
message("[data_export] computing last_ejc_tx_pos from structures.rds ...")
struct_sub <- structures[isoform_id %in% subset2_nmd_ids]

last_ejc_tx_pos <- vapply(seq_len(nrow(struct_sub)), function(i) {
  ne <- struct_sub$n_exons[i]
  if (is.na(ne) || ne < 2L) return(NA_integer_)
  starts <- struct_sub$exon_starts[[i]]
  ends   <- struct_sub$exon_ends[[i]]
  widths <- as.integer(ends - starts + 1L)
  tx_len <- sum(widths)
  last_tx_exon_width <- if (struct_sub$strand[i] == "+") widths[ne] else widths[1L]
  as.integer(tx_len - last_tx_exon_width)
}, integer(1))

tx_len_by_iso <- vapply(seq_len(nrow(struct_sub)), function(i) {
  starts <- struct_sub$exon_starts[[i]]
  ends   <- struct_sub$exon_ends[[i]]
  as.integer(sum(ends - starts + 1L))
}, integer(1))

struct_lookup <- data.table(
  isoform_id       = struct_sub$isoform_id,
  n_exons          = struct_sub$n_exons,
  tx_length_nt     = tx_len_by_iso,
  last_ejc_tx_pos  = last_ejc_tx_pos
)

# ── Assemble SF30 table ────────────────────────────────────────────────
sf_base <- merge(
  ref_atg_c2[, .(comparator_isoform_id, reference_isoform_id, category,
                 comp_stop_tx_pos, n_downstream_ejc, original_ptc,
                 comp_orf_length, ref_orf_length)],
  struct_lookup,
  by.x = "comparator_isoform_id", by.y = "isoform_id", all.x = TRUE
)
sf_base[, stop_to_last_ejc_via_ref_nt := last_ejc_tx_pos - comp_stop_tx_pos]

sf30 <- merge(
  sf_base,
  meta_logfc,
  by.x = "comparator_isoform_id", by.y = "txid",
  all.x = TRUE
)
setnames(sf30, "comparator_isoform_id", "txid")

# Retain all 1,166 for provenance; report defined-metric counts.
n_total    <- nrow(sf30)
n_have_dist <- sum(!is.na(sf30$stop_to_last_ejc_via_ref_nt))
n_have_both <- sum(!is.na(sf30$stop_to_last_ejc_via_ref_nt) & !is.na(sf30$mean_logFC))
message(sprintf("[data_export] Subset 2 NMD arm: %d rows", n_total))
message(sprintf("[data_export]   with defined stop-to-last-EJC distance: %d", n_have_dist))
message(sprintf("[data_export]   with defined distance + mashr logFC (SF30 plot n): %d", n_have_both))
message(sprintf("[data_export]   of which effectively_ptc: %d",
                sum(sf30$category == "effectively_ptc" & !is.na(sf30$stop_to_last_ejc_via_ref_nt))))

# ── SF31 population: comparators with n_downstream_ejc >= 1 ────────────
sf31 <- sf30[!is.na(n_downstream_ejc) & n_downstream_ejc >= 1L
             & !is.na(mean_logFC)]
message(sprintf("[data_export] SF31 population (n_downstream_ejc >= 1 & logFC): n = %d",
                 nrow(sf31)))

# ── Write TSVs ─────────────────────────────────────────────────────────
dir.create(SF30_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SF31_DIR, showWarnings = FALSE, recursive = TRUE)

fwrite(sf30, file.path(SF30_DIR, "sf30_ptc_distance_logfc.tsv"), sep = "\t")
fwrite(sf31, file.path(SF31_DIR, "sf31_ejc_count_logfc.tsv"),    sep = "\t")

message("[data_export] wrote:")
message("  ", file.path(SF30_DIR, "sf30_ptc_distance_logfc.tsv"))
message("  ", file.path(SF31_DIR, "sf31_ejc_count_logfc.tsv"))
