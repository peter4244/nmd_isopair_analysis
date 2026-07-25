#!/usr/bin/env Rscript
# =============================================================================
# add_biotype_to_gmap.R — add transcript_biotype to the deposited transcript map.
#
# INTERNAL packaging step. Not shipped; see DEPOSIT_PROVENANCE_INTERNAL.md.
#
# WHY (2026-07-25)
# ---------------------------------------------------------------------------
# Isoform_Level_Quantification.Rmd resolves transcript biotype as:
#     GTF transcript_type  ->  GTF gene_type  ->  LIVE biomaRt query
# Neither deposited nor legacy GTF carries either attribute (both are isocall/SQANTI
# output: gene_name, gene_id, transcript_id only), so the biomaRt branch always fires.
# That makes §1 fail outright without internet — it died on 2026-07-25 when Ensembl
# returned 504 — and, unpinned, silently tracked whatever Ensembl was current.
#
# The deposit had no offline biotype source: gmap_txlevel_ENSGv115_2025.08.12.rds
# carried only IDs and names. This adds one column so §1 reproduces offline and
# deterministically.
#
# SOURCE: GENCODE v49 (= Ensembl 115, matching the map's own ENSGv115 vintage),
#         local under reference_files/. Both annotation files are used:
#           gencode.v49.primary_assembly.annotation.gtf.gz
#           gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz   (scaffolds/patches)
#         The scaffold set alone yields exactly 533,740 transcripts — the map's own row
#         count — confirming the map was built from the same annotation release.
#
# GUARANTEES (all asserted below):
#   1. row count and row ORDER unchanged
#   2. every pre-existing column byte-identical
#   3. transcript_biotype is 100% populated (no NA)
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

REF  <- path.expand("~/claude_projects/nmd/reference_files")
DEP  <- path.expand("~/claude_projects/nmd_deposit_2026/source_data/annotation")
MAP  <- file.path(DEP, "gmap_txlevel_ENSGv115_2025.08.12.rds")
GTF  <- file.path(REF, "gencode.v49.primary_assembly.annotation.gtf.gz")
GFF3 <- file.path(REF, "gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz")
stopifnot("annotation sources missing" = all(file.exists(c(MAP, GTF, GFF3))))

cat("=== extracting transcript_type from GENCODE v49 ===\n")
tmp <- tempfile()
system2("bash", c("-c", shQuote(sprintf(
  'gzcat %s | awk -F"\\t" \'$3=="transcript"{
      match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16)
      match($9,/transcript_type "[^"]+"/); b=substr($9,RSTART+17,RLENGTH-18)
      if (t!="" && b!="") print t"\\t"b }\' > %s', shQuote(GTF), shQuote(tmp)))))
tmp2 <- tempfile()
system2("bash", c("-c", shQuote(sprintf(
  'gzcat %s | awk -F"\\t" \'$3=="transcript"{
      match($9,/transcript_id=[^;]+/); t=substr($9,RSTART+14,RLENGTH-14)
      match($9,/transcript_type=[^;]+/); b=substr($9,RSTART+16,RLENGTH-16)
      if (t!="" && b!="") print t"\\t"b }\' > %s', shQuote(GFF3), shQuote(tmp2)))))

bt <- unique(rbindlist(list(
  fread(tmp,  header = FALSE, col.names = c("txid", "transcript_biotype")),
  fread(tmp2, header = FALSE, col.names = c("txid", "transcript_biotype")))))
bt <- bt[!duplicated(txid)]
cat(sprintf("  %d transcripts with a biotype\n", nrow(bt)))
unlink(c(tmp, tmp2))

cat("\n=== joining onto the transcript map ===\n")
g <- readRDS(MAP)
orig_cols <- names(g)
orig_ids  <- g$ensembl_transcript_id_version
cat(sprintf("  map: %d rows x %d cols\n", nrow(g), ncol(g)))
stopifnot("transcript_biotype already present" = !"transcript_biotype" %in% orig_cols)

g$transcript_biotype <- bt$transcript_biotype[match(g$ensembl_transcript_id_version, bt$txid)]

n_na <- sum(is.na(g$transcript_biotype))
cat(sprintf("  biotype populated: %d / %d (%.2f%%)   NA: %d\n",
            nrow(g) - n_na, nrow(g), 100 * (1 - n_na / nrow(g)), n_na))

# ---- GUARANTEES --------------------------------------------------------------
cat("\n=== verification ===\n")
stopifnot("row order changed"  = identical(g$ensembl_transcript_id_version, orig_ids))
stopifnot("a pre-existing column changed" =
            identical(g[orig_cols], readRDS(MAP)[orig_cols]))
stopifnot("transcript_biotype is not fully populated" = n_na == 0)
cat("  row order unchanged            : TRUE\n")
cat("  pre-existing columns identical : TRUE\n")
cat("  biotype fully populated        : TRUE\n")
cat(sprintf("  distinct biotypes              : %d (top: %s)\n",
            length(unique(g$transcript_biotype)),
            paste(names(sort(table(g$transcript_biotype), decreasing = TRUE))[1:3], collapse = ", ")))

saveRDS(g, MAP)
cat(sprintf("\nwrote %s (%.1f MB)\n", MAP, file.size(MAP) / 1e6))
cat("verification passed.\n")
