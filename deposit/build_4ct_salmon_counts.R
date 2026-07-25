#!/usr/bin/env Rscript
# =============================================================================
# build_4ct_salmon_counts.R — build the 4-cell-type short-read gene count matrix
# for the Zenodo source-data deposit.
#
# Deposited alongside its output so the transformation is auditable.
#
# INPUT   nf-core/rnaseq salmon output (SummarizedExperiment),
#         salmon.merged.gene_counts_length_scaled.rds — genes x 38 samples.
#         This is an M2 product, and M2 is the accepted reproducibility gap, so
#         it is a genuine starting file rather than a regenerable intermediate.
#
# OUTPUT  salmon_gene_counts_4ct.csv — genes x 26 samples, columns named
#         <donor>_<treatment>_<CT> with canonical CT names, matching the
#         convention used by nmd_lungcells_counts_4ct.csv.
#
# WHY THIS IS NOT A SIMPLE COLUMN FILTER — three traps in the source names:
#   1. The three DD_ALI donors (001V, 027U, 029T) are labelled plain "DD" here,
#      unlike the long-read matrix where DD_ALI is explicit. Selecting on the
#      "DD" prefix would wrongly retain 6 ALI samples.
#   2. The short-read CT code is "AT2" where the long-read code is "AT".
#   3. Two FB columns carry a "_clean" suffix (FB029T_*_clean).
# To avoid re-deriving that logic (and re-deriving its bugs), the 26 samples are
# taken from the SR DGEList itself, which is what the published pipeline built.
#
# Self-verifies: the deposited counts must equal the DGEList counts exactly.
# =============================================================================
suppressPackageStartupMessages({ library(SummarizedExperiment); library(data.table) })

REPO   <- path.expand("~/claude_projects/nmd")
SALMON <- file.path(REPO, "nmd_fig_data/salmon.merged.gene_counts_length_scaled.rds")
DGE_SR <- file.path(REPO, "nmd_fig_data/dge_shortread_gene_2026.3.2.rds")
OUTDIR <- path.expand("~/claude_projects/nmd_deposit_2026/source_data")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
OUT    <- file.path(OUTDIR, "salmon_gene_counts_4ct.csv")

CT_RENAME <- c(AT = "AT2", AT2 = "AT2", DD = "LAE", FB = "FB", MV = "MV")

se  <- readRDS(SALMON)
dge <- readRDS(DGE_SR)
cat(sprintf("salmon : %d genes x %d samples\n", nrow(se), ncol(se)))
cat(sprintf("SR DGEList (the published 26): %d genes x %d samples\n", nrow(dge), ncol(dge)))

# ---- map each published sample back to its salmon column ---------------------
# DGEList name -> salmon name: AT<donor> is AT2<donor>; FB029T_* gained "_clean".
to_salmon <- function(x) {
  cands <- unique(c(x, sub("^AT", "AT2", x), paste0(x, "_clean"),
                    paste0(sub("^AT", "AT2", x), "_clean")))
  hit <- cands[cands %in% colnames(se)]
  if (length(hit) != 1L) stop("cannot uniquely map sample: ", x)
  hit
}
map <- data.table(dge_name = colnames(dge))
map[, salmon_name := vapply(dge_name, to_salmon, character(1))]

# ---- canonical deposit names: <donor>_<treatment>_<CT> ----------------------
map[, `:=`(
  ct    = CT_RENAME[sub("^([A-Z]+?)[0-9].*$", "\\1", dge_name)],
  donor = sub("^[A-Z]+([0-9]+[A-Z])_.*$", "\\1", dge_name),
  trt   = sub("^.*_([^_]+)$", "\\1", sub("_clean$", "", dge_name))
)]
map[, deposit_name := paste(donor, trt, ct, sep = "_")]
stopifnot("unresolved cell type" = !any(is.na(map$ct)))
stopifnot("duplicate deposit names" = !any(duplicated(map$deposit_name)))
cat("\nper cell type:\n"); print(map[, .N, by = ct][order(ct)])
cat("\nmapping examples:\n"); print(head(map[, .(salmon_name, dge_name, deposit_name)], 4))
cat("\ndropped from salmon (not in the published 26):\n")
print(sort(setdiff(colnames(se), map$salmon_name)))

# ---- build + verify ---------------------------------------------------------
cnt <- assay(se, "counts")[, map$salmon_name, drop = FALSE]
colnames(cnt) <- map$deposit_name

ref <- dge$counts
stopifnot("gene sets differ" = setequal(rownames(cnt), rownames(ref)))
aligned <- cnt[rownames(ref), ]
colnames(aligned) <- map$deposit_name[match(colnames(ref), map$dge_name)]
ok <- all(abs(round(aligned[, map$deposit_name[match(colnames(ref), map$dge_name)]]) - round(ref)) < 1e-6)
cat(sprintf("\ngene sets match : TRUE\nCOUNTS IDENTICAL: %s\n", ok))
stopifnot("counts differ from the published DGEList — transformation is NOT faithful" = ok)

out <- data.table(gene_id = rownames(cnt)); out <- cbind(out, as.data.table(cnt))
fwrite(out, OUT)
cat(sprintf("\nwrote %s (%.1f MB) — %d genes x %d samples\n",
            OUT, file.size(OUT)/1e6, nrow(out), ncol(out)-1L))
