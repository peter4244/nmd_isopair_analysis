# SF30 / SF31 data_export.R — canonical Isopair-side derivation
#
# Reproduces the population construction from
#   nmd:results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_mashr.Rmd
# for the two figures:
#   SF30 — PTC distance dose response   (Rmd fig3: NMD Response vs Stop Codon Distance)
#   SF31 — NMD effect size by EJC count (Rmd fig4: mean_logFC vs downstream EJC count)
#
# Both figures use the same population construction:
#   Start:  gene-matched C2 (NMD) comparator isoforms from the Isopair pipeline
#   Add:    ptc_distance, has_ptc, orf_length, n_downstream_ejcs from ptc.rds
#   Add:    mean_logFC from the mashr meta-analysis (4-condition posterior mean)
#
# Output: two TSVs, one per SF dir. Both scripts read those TSVs directly and
# never re-compute — matches the canonical "figure code reads Rmd exports"
# convention documented in project_nmd_analysis_repo memory.

suppressPackageStartupMessages({
  library(data.table)
})

# ── Provenance paths ────────────────────────────────────────────────────
ISOPAIR_ROOT <- "~/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper"
NMD_ROOT     <- "~/claude_projects/nmd"

PTC_RDS       <- file.path(ISOPAIR_ROOT, "data_mashr", "ptc.rds")
MASHR_RDS     <- file.path(NMD_ROOT, "isocall_dge/mashr", "mashr_isoform_model_2026.3.10.rds")
PROFILES_C2   <- file.path(ISOPAIR_ROOT, "data_mashr", "profiles_c2_allsamples.rds")

SF30_DIR <- file.path("~/claude_projects/nmd/figures/SupplementalFigures",
                       "SF30_PTCDistanceDoseResponse", "data")
SF31_DIR <- file.path("~/claude_projects/nmd/figures/SupplementalFigures",
                       "SF31_NMDEffectByEJCCount", "data")

# ── Load the source objects ─────────────────────────────────────────────
message("[data_export] reading ptc.rds ...")
ptc <- readRDS(PTC_RDS)

message("[data_export] reading profiles_c2_allsamples.rds ...")
profiles_c2 <- readRDS(PROFILES_C2)
comp_ids_as <- profiles_c2$comparator_isoform_id

message("[data_export] reading mashr_isoform_model_2026.3.10.rds ...")
mashr_model <- readRDS(MASHR_RDS)
# Match the Rmd's meta_logfc definition — posterior-mean averaged across the
# 4 SMG1i-vs-DMSO conditions (Smg1i_in_AT, Smg1i_in_DD, Smg1i_in_FB, Smg1i_in_MV).
pm <- as.data.frame(mashr_model$result$PosteriorMean)
smg1i_cols <- grep("^Smg1i_in_", colnames(pm), value = TRUE)
meta_logfc <- data.frame(
  txid       = rownames(pm),
  mean_logFC = rowMeans(pm[, smg1i_cols, drop = FALSE], na.rm = TRUE)
)

# ── SF30 population: nmd_logfc_dist (comparators with ptc_distance AND logFC) ──
ptc_sub <- ptc[ptc$isoform_id %in% comp_ids_as & !is.na(ptc$ptc_distance),
                c("isoform_id", "ptc_distance", "has_ptc",
                  "orf_length", "n_downstream_ejcs")]

nmd_logfc_dist <- merge(
  meta_logfc[meta_logfc$txid %in% comp_ids_as, ],
  ptc_sub, by.x = "txid", by.y = "isoform_id"
)
message(sprintf("[data_export] SF30 population: n = %d comparators with ptc_distance + mashr logFC",
                 nrow(nmd_logfc_dist)))

# ── SF31 population: PTC-positive subset ────────────────────────────────
as_ptc_pos <- nmd_logfc_dist[nmd_logfc_dist$has_ptc == TRUE, ]
message(sprintf("[data_export] SF31 population (PTC+ subset): n = %d",
                 nrow(as_ptc_pos)))

# ── Write TSVs ──────────────────────────────────────────────────────────
dir.create(SF30_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SF31_DIR, showWarnings = FALSE, recursive = TRUE)

fwrite(nmd_logfc_dist,
       file.path(SF30_DIR, "sf30_ptc_distance_logfc.tsv"),
       sep = "\t")
fwrite(as_ptc_pos,
       file.path(SF31_DIR, "sf31_ejc_count_logfc.tsv"),
       sep = "\t")

message("[data_export] wrote:")
message("  ", file.path(SF30_DIR, "sf30_ptc_distance_logfc.tsv"))
message("  ", file.path(SF31_DIR, "sf31_ejc_count_logfc.tsv"))
