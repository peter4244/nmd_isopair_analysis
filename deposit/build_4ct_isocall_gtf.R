#!/usr/bin/env Rscript
# =============================================================================
# build_4ct_isocall_gtf.R — deposit the isocall-universe GTF for the 4 cell types.
#
# INTERNAL packaging step. Not shipped; see DEPOSIT_PROVENANCE_INTERNAL.md.
#
# WHY (2026-07-25)
# ---------------------------------------------------------------------------
# The §3/§4 Isopair pipeline builds its transcript->gene map from the *isocall*
# GTF, not the SQANTI-filtered one: gene_map.rds has 645,273 rows. The deposit
# previously carried only sqanti/nmd_lungcells_filtered.gtf (614,992
# transcripts), which omits ENST00000441168.6 — so a gene map built from it
# would assign that isoform gene_id = NA and silently drop it from its gene's
# isoform-proportion denominator, defeating nmd_isocall_counts_4ct.csv.
#
# Deposited GZIPPED. The published code already reads this file gzipped
# (`nmd_isocall.isoforms.gtf.gz`), and it costs ~45 MB instead of ~1.1 GB.
#
# INPUT   nmd_isocall.isoforms.gtf.gz   (645,273 transcripts)
# OUTPUT  nmd_isocall_4ct.gtf.gz        (614,993 transcripts)
#
# Trimmed to exactly the isoform set of nmd_isocall_counts_4ct.csv, so the two
# isocall-universe files agree. Two-pass so original line order is preserved and
# a gene record is kept whenever the gene retains >= 1 transcript.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

REPO <- path.expand("~/claude_projects/nmd")
SRC  <- file.path(REPO, "isocall/nmd_lungcells/results/call/nmd_isocall.isoforms.gtf.gz")
DEP  <- path.expand("~/claude_projects/nmd_deposit_2026/source_data")
CNT  <- file.path(DEP, "nmd_isocall_counts_4ct.csv")
EXPR <- file.path(REPO, "results/isoform_transitions/Version_6.0/isopair_wrapper",
                  "data_mashr/expression_data.rds")
OUT  <- file.path(DEP, "nmd_isocall_4ct.gtf.gz")
WORK <- tempfile(); dir.create(WORK)

cat("=== keep-set from nmd_isocall_counts_4ct.csv ===\n")
keep <- fread(CNT, select = 1, showProgress = FALSE)[[1]]
cat(sprintf("  %d isoforms\n", length(keep)))
stopifnot("ENST00000441168.6 must be in the isocall count matrix" =
            "ENST00000441168.6" %in% keep)
kf <- file.path(WORK, "keep.txt"); writeLines(keep, kf)

cat("\n=== pass 1: genes retaining >= 1 transcript ===\n")
gf <- file.path(WORK, "genes.txt")
system2("bash", c("-c", shQuote(sprintf(
  'gzcat %s | awk -v kf=%s \'BEGIN{while((getline l < kf)>0) k[l]}
     match($0,/transcript_id "[^"]+"/){
       tid=substr($0,RSTART+15,RLENGTH-16)
       if (tid in k && match($0,/gene_id "[^"]+"/)) print substr($0,RSTART+9,RLENGTH-10)
     }\' | sort -u > %s', shQuote(SRC), shQuote(kf), shQuote(gf)))))
cat(sprintf("  %d genes\n", length(readLines(gf))))

cat("\n=== pass 2: emit in original order, gzipped ===\n")
system2("bash", c("-c", shQuote(sprintf(
  'gzcat %s | awk -v kf=%s -v gf=%s \'BEGIN{while((getline l < kf)>0) k[l]
       while((getline g < gf)>0) okg[g]}
     { tid=""; gid=""
       if (match($0,/transcript_id "[^"]+"/)) tid=substr($0,RSTART+15,RLENGTH-16)
       if (match($0,/gene_id "[^"]+"/))       gid=substr($0,RSTART+9,RLENGTH-10)
       if (tid!="")      { if (tid in k) print }
       else if (gid!="") { if (gid in okg) print }
       else print }\' | gzip -c > %s', shQuote(SRC), shQuote(kf), shQuote(gf), shQuote(OUT)))))

cat("\n=== verification ===\n")
got <- system2("bash", c("-c", shQuote(sprintf(
  'gzcat %s | grep -o \'transcript_id "[^"]*"\' | sort -u | wc -l', shQuote(OUT)))), stdout = TRUE)
got <- as.integer(trimws(got))
cat(sprintf("  transcripts in output : %d (expected %d)\n", got, length(keep)))
stopifnot("transcript count does not match the isocall count matrix" = got == length(keep))

has <- system2("bash", c("-c", shQuote(sprintf(
  'gzcat %s | grep -c ENST00000441168.6 || true', shQuote(OUT)))), stdout = TRUE)
cat(sprintf("  ENST00000441168.6 records : %s\n", trimws(has)))
stopifnot("ENST00000441168.6 absent — the whole point of this file" =
            as.integer(trimws(has)) > 0)

univ <- rownames(readRDS(EXPR))
cat(sprintf("  published Isopair universe covered : %d / %d\n",
            sum(univ %in% keep), length(univ)))
stopifnot("published Isopair isoforms missing" = all(univ %in% keep))

cat(sprintf("\nwrote %s (%.1f MB)\n", OUT, file.size(OUT) / 1e6))
cat("verification passed.\n")
unlink(WORK, recursive = TRUE)
