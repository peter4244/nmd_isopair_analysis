#!/usr/bin/env Rscript
# =============================================================================
# build_4ct_isocall_counts.R — build the 4-cell-type isocall-universe isoform
# count matrix for the Zenodo source-data deposit.
#
# "isocall universe" = quantified by isocall WITHOUT the SQANTI structural filter,
# then trimmed to features observed in these 26 samples. It is not a raw dump:
# an isoform with zero counts in all 26 is not part of this dataset and is
# removed, exactly as for every other deposited file.
#
# INTERNAL packaging step. Not shipped; see DEPOSIT_PROVENANCE_INTERNAL.md.
#
# WHY THIS FILE EXISTS (added 2026-07-25)
# ---------------------------------------------------------------------------
# The two published analysis branches did NOT start from the same count matrix:
#
#   §1/§2 gene- and isoform-level mashr -> SQANTI-filtered  (645,272 isoforms)
#   §3/§4 Isopair                       -> raw isocall      (645,273 isoforms)
#
# The two differ by exactly one isoform, ENST00000441168.6 (ANKMY1,
# ENSG00000144504.17), which SQANTI drops — it has a degenerate 1-bp exon. That
# isoform is NOT incidental to §3/§4: it clears both the 5% isoform-proportion
# filter and filterByExpr and sits in the published 95,623-isoform Isopair
# analysis universe (verified in expression_data.rds).
#
# Depositing only the filtered matrix therefore covered 95,622 / 95,623 of the
# published Isopair universe, and substituting it is not inert: dropping
# ENST00000441168.6 also removes it from its gene's isoform-proportion
# denominator, which moves sibling ENSG00000144504.17.novel34 from 0.04348 to
# 0.05000 — across the `>= 0.05` boundary, i.e. a two-way (-1, +1) change.
#
# So we deposit BOTH matrices and let each branch read the one it actually used.
#
# INPUT   isocall output: nmd_isocall.count_matrix.txt
#         645,273 isoforms x 38 samples, comma-delimited despite the .txt suffix.
#         Column names: Sample<N>_<CT>_<donor>_<treatment>
#
# OUTPUT  nmd_isocall_counts_4ct.csv
#         26 samples, all-zero rows removed. Columns named <donor>_<treatment>_<CT> — IDENTICAL column
#         names and order to nmd_lungcells_counts_4ct.csv, so the two deposit
#         matrices are directly comparable.
#         Rows trimmed to features observed in these 26 samples, matching the
#         deposit-wide policy (see trim_deposit_to_observed.sh).
#
# Self-verifies and aborts on mismatch:
#   1. every isoform of the published Isopair universe is present
#   2. ENST00000441168.6 specifically is present and non-zero
#   3. counts for shared isoforms are identical to the filtered deposit matrix
#      (both derive from the same quantification, so any drift is a defect)
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

REPO    <- path.expand("~/claude_projects/nmd")
ISOCALL <- file.path(REPO, "isocall/nmd_lungcells/results/call/nmd_isocall.count_matrix.txt")
EXPR    <- file.path(REPO, "results/isoform_transitions/Version_6.0/isopair_wrapper",
                     "data_mashr/expression_data.rds")
OUTDIR  <- path.expand("~/claude_projects/nmd_deposit_2026/source_data")
FILTERED <- file.path(OUTDIR, "nmd_lungcells_counts_4ct.csv")
OUT      <- file.path(OUTDIR, "nmd_isocall_counts_4ct.csv")

KEEP_CT   <- c("AT", "DD", "FB", "MV")
CT_RENAME <- c(AT = "AT2", DD = "LAE", FB = "FB", MV = "MV")

cat("=== reading raw isocall count matrix ===\n")
cm <- fread(ISOCALL, showProgress = FALSE)
cat(sprintf("  %d isoforms x %d columns\n", nrow(cm), ncol(cm)))
stopifnot("first column should be the isoform id" = names(cm)[1] == "id")

# ---- parse column names (same convention as build_4ct_count_matrix.R) --------
parsed <- rbindlist(lapply(names(cm)[-1], function(x) {
  p  <- strsplit(sub("^Sample\\d+_", "", x), "_")[[1]]
  data.table(raw = x, donor = p[length(p) - 1], treatment = p[length(p)],
             ct = paste(p[1:(length(p) - 2)], collapse = "_"))
}))
keep <- parsed[ct %in% KEEP_CT]
cat(sprintf("retaining %d samples; dropping %d (%s)\n", nrow(keep), nrow(parsed) - nrow(keep),
            paste(sort(unique(parsed[!ct %in% KEEP_CT]$ct)), collapse = ", ")))
stopifnot("expected 26 retained samples" = nrow(keep) == 26)

keep[, new_name := paste(donor, treatment, CT_RENAME[ct], sep = "_")]
stopifnot("duplicate output column names" = !any(duplicated(keep$new_name)))

out <- cm[, c("id", keep$raw), with = FALSE]
setnames(out, c("id", keep$new_name))

# ---- align columns to the filtered deposit matrix ----------------------------
cat("\n=== aligning to the filtered deposit matrix ===\n")
filt_cols <- names(fread(FILTERED, nrows = 0))
stopifnot("column sets differ from the filtered deposit matrix" =
            setequal(names(out), filt_cols))
setcolorder(out, filt_cols)
stopifnot("column order not aligned" = identical(names(out), filt_cols))
cat("  column names and order match nmd_lungcells_counts_4ct.csv\n")

# ---- trim to observed features (deposit-wide policy) -------------------------
mat  <- as.matrix(out[, -1])
keep_rows <- rowSums(mat) > 0
cat(sprintf("\ntrimming to observed: %d -> %d isoforms\n", nrow(out), sum(keep_rows)))
out <- out[keep_rows]

# ---- VERIFY ------------------------------------------------------------------
cat("\n=== verification ===\n")

expr <- readRDS(EXPR)
univ <- rownames(expr)
missing <- setdiff(univ, out$id)
cat(sprintf("  published Isopair universe : %d isoforms\n", length(univ)))
cat(sprintf("  missing from this matrix   : %d\n", length(missing)))
stopifnot("published Isopair isoforms are missing from the deposit matrix" =
            length(missing) == 0)

i <- match("ENST00000441168.6", out$id)
stopifnot("ENST00000441168.6 absent — the whole point of this file" = !is.na(i))
tot <- sum(as.numeric(out[i, -1]))
cat(sprintf("  ENST00000441168.6 present  : TRUE (total counts across 26 samples = %g)\n", tot))
stopifnot("ENST00000441168.6 is all-zero in the 26 manuscript samples" = tot > 0)

cat("  cross-checking shared isoforms against the filtered deposit matrix...\n")
filt <- fread(FILTERED, showProgress = FALSE)
shared <- intersect(out$id, filt$id)
a <- as.matrix(out[match(shared, out$id),   -1])
b <- as.matrix(filt[match(shared, filt$id), -1])
same <- identical(dim(a), dim(b)) && all(a == b)
cat(sprintf("  shared isoforms            : %d\n  COUNTS IDENTICAL           : %s\n",
            length(shared), same))
stopifnot("counts disagree with the filtered deposit matrix — one of them is wrong" = same)

fwrite(out, OUT)
cat(sprintf("\nwrote %s (%.1f MB)\n", OUT, file.size(OUT) / 1e6))
cat("verification passed.\n")
