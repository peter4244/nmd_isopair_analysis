#!/usr/bin/env Rscript
# =============================================================================
# build_4ct_count_matrix.R — build the 4-cell-type isoform count matrix for the
# Zenodo source-data deposit.
#
# This script is deposited alongside its output so the transformation is
# auditable rather than magic (REPRODUCIBILITY_PLAN §2b).
#
# INPUT   SQANTI3 output: nmd_lungcells_filtered.count_matrix.txt
#         645,272 isoforms x 38 samples, comma-delimited despite the .txt suffix.
#         Column names: Sample<N>_<CT>_<donor>_<treatment>
#
# OUTPUT  nmd_lungcells_counts_4ct.csv
#         645,272 isoforms x 26 samples (id + 26 columns)
#         Column names: <donor>_<treatment>_<CT>, canonical CT names (AT2/LAE)
#
# TWO TRANSFORMATIONS, both deliberate:
#   1. SAMPLES  drop the 12 ALI/submerged samples (DD_ALI, DO) that are not in
#      the manuscript. The paper mentions them zero times.
#   2. NAMES    Sample5_AT_001V_DMSO -> 001V_DMSO_AT2 ; internal codes AT/DD are
#      normalised to the published names AT2/LAE at source, so no downstream
#      script needs a relabelling step.
#
# ROWS ARE NOT TOUCHED. All 645,272 isoforms are retained, including the 30,280
# that are all-zero across the 26 retained samples. This is required: §3 ¶1
# ("% transcriptional output lost") sums over the full isoform universe, and
# claims 1.12/1.13 define cell-type-restricted expression over it too.
#
# The script self-verifies against the current DGEList and fails loudly on any
# mismatch — a faithful transformation must reproduce it exactly.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

REPO    <- path.expand("~/claude_projects/nmd")
SQANTI  <- file.path(REPO, "sqanti/nmd_lungcells/results/nmd_lungcells_filtered.count_matrix.txt")
DGE_REF <- file.path(REPO, "nmd_fig_data/dge_isoform_longread_2026.3.3.rds")
OUTDIR  <- path.expand("~/claude_projects/nmd_deposit_2026/source_data")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
OUT     <- file.path(OUTDIR, "nmd_lungcells_counts_4ct.csv")

KEEP_CT  <- c("AT", "DD", "FB", "MV")                 # internal codes retained
CT_RENAME <- c(AT = "AT2", DD = "LAE", FB = "FB", MV = "MV")

cat("=== reading SQANTI count matrix ===\n")
cm <- fread(SQANTI, showProgress = FALSE)
cat(sprintf("  %d isoforms x %d columns\n", nrow(cm), ncol(cm)))
stopifnot("first column should be the isoform id" = names(cm)[1] == "id")

# ---- parse column names ------------------------------------------------------
raw <- names(cm)[-1]
parsed <- rbindlist(lapply(raw, function(x) {
  p  <- strsplit(sub("^Sample\\d+_", "", x), "_")[[1]]
  ct <- paste(p[1:(length(p) - 2)], collapse = "_")     # everything before donor+treatment
  data.table(raw = x, donor = p[length(p) - 1], treatment = p[length(p)], ct = ct)
}))
cat("\ncell types found in the source matrix:\n"); print(parsed[, .N, by = ct][order(-N)])

keep <- parsed[ct %in% KEEP_CT]
drop <- parsed[!ct %in% KEEP_CT]
cat(sprintf("\nretaining %d samples; dropping %d (%s)\n",
            nrow(keep), nrow(drop), paste(sort(unique(drop$ct)), collapse = ", ")))
stopifnot("expected 26 retained samples" = nrow(keep) == 26)

keep[, new_name := paste(donor, treatment, CT_RENAME[ct], sep = "_")]
stopifnot("duplicate output column names" = !any(duplicated(keep$new_name)))
cat("\nrenaming examples:\n")
print(head(keep[, .(raw, new_name)], 4))

# ---- subset + rename ---------------------------------------------------------
out <- cm[, c("id", keep$raw), with = FALSE]
setnames(out, c("id", keep$new_name))
cat(sprintf("\noutput: %d isoforms x %d samples\n", nrow(out), ncol(out) - 1L))
stopifnot("row count must be unchanged" = nrow(out) == nrow(cm))

# ---- VERIFY against the current DGEList --------------------------------------
cat("\n=== verification against the existing DGEList ===\n")
dge <- readRDS(DGE_REF)
ref <- dge$counts
cat(sprintf("  reference DGEList: %d x %d\n", nrow(ref), ncol(ref)))

# the DGEList keeps internal CT codes in its column names; map to canonical
ref_new <- vapply(colnames(ref), function(x) {
  p <- strsplit(x, "_")[[1]]
  paste(p[1], p[2], CT_RENAME[paste(p[3:length(p)], collapse = "_")], sep = "_")
}, character(1))

mat <- as.matrix(out[, -1]); rownames(mat) <- out$id
stopifnot("isoform sets differ"    = setequal(rownames(mat), rownames(ref)))
stopifnot("sample sets differ"     = setequal(colnames(mat), ref_new))
mat_aligned <- round(mat[rownames(ref), match(ref_new, colnames(mat))])
identical_counts <- all(mat_aligned == round(ref))
cat(sprintf("  isoform sets match : TRUE\n  sample sets match  : TRUE\n  COUNTS IDENTICAL   : %s\n",
            identical_counts))
stopifnot("counts differ from the DGEList — transformation is NOT faithful" = identical_counts)

fwrite(out, OUT)
cat(sprintf("\nwrote %s (%.1f MB)\n", OUT, file.size(OUT) / 1e6))
cat("verification passed: the deposited matrix reproduces the DGEList exactly.\n")
