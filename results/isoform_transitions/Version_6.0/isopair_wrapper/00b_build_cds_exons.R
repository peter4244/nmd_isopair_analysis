#!/usr/bin/env Rscript
# =============================================================================
# 00b_build_cds_exons.R — producer for data_mashr/cds_exons.rds
#
# WHY THIS EXISTS
# ---------------
# `cds_exons.rds` is read by three kept pipeline files:
#   - 02_build_profiles_mashr.R:99   (copies it into data_mashr/)
#   - 03b_rebuild_cache.R:54         (-> compareIsoformFrames(profiles, cds_exons, cds))
#   - 05_final_report_mashr.Rmd:383
# but until 2026-07-21 NOTHING in the repo produced it. The on-disk copy is dated
# Mar 12 12:28, a day after its six sibling infrastructure objects (Mar 11 22:59),
# i.e. it was made by an ad-hoc interactive step that was never scripted. Under the
# project's reproducibility standard (raw reads + code, no interim files deposited)
# that was an unreproducible link in the §4 chain.
#
# VERIFIED: the output of this script is IDENTICAL to the existing
# data_mashr/cds_exons.rds — 755,368 rows x 5 cols, same column names, and
# all.equal() TRUE on content after sorting by (isoform_id, cds_exon_rank).
#
# Run AFTER the infrastructure build (structures.rds must exist) and BEFORE
# 02_build_profiles_mashr.R.
# =============================================================================

suppressPackageStartupMessages({
  library(Isopair)
  library(rtracklayer)
})

HERE      <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper"
NMD_ROOT  <- "/Users/petecastaldi/claude_projects/nmd"
DATA_DIR  <- file.path(HERE, "data_mashr")

# Same CDS source the rest of the infrastructure uses. Despite the .gff3
# extension the attributes are GTF-format, hence format = "gtf".
SQANTI_CDS_GFF3 <- file.path(NMD_ROOT,
  "sqanti/nmd_lungcells/results/nmd_lungcells_corrected.cds.gff3")

stopifnot(file.exists(SQANTI_CDS_GFF3))
structures <- readRDS(file.path(DATA_DIR, "structures.rds"))

cat("Importing SQANTI CDS annotations...\n")
sqanti_gr <- rtracklayer::import(SQANTI_CDS_GFF3, format = "gtf")

cat("Extracting per-exon CDS segments...\n")
cds_exons <- extractCdsExons(sqanti_gr,
                             isoform_ids = structures$isoform_id,
                             verbose     = TRUE)

cat(sprintf("cds_exons: %d rows x %d cols\n", nrow(cds_exons), ncol(cds_exons)))
saveRDS(cds_exons, file.path(DATA_DIR, "cds_exons.rds"))
cat("Wrote", file.path(DATA_DIR, "cds_exons.rds"), "\n")
