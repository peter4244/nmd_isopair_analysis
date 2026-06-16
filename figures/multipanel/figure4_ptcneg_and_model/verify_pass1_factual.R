#!/usr/bin/env Rscript
# Pass 1: Factual accuracy verification for Figure 4 (NEW 2x2 layout, 2026-06-15)
# Recomputes n, median, IQR per group from long-format TSVs and compares
# against values asserted in the descriptives TSVs.
# Panels A/B: Section A (all-3-ENST + coding) n=72/118/190
# Panels C/D: Section C (ENST-ref + ref-AUG-traceable, re-intersected) n=1050/113/1166

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/Users/petecastaldi/claude_projects/nmd/figures/multipanel/figure4_ptcneg_and_model"
fails <- list()
pass_ct <- 0L
fail_ct <- 0L

check <- function(label, expected, actual, tol = 1e-6) {
  ok <- isTRUE(all.equal(as.numeric(expected), as.numeric(actual), tolerance = tol))
  if (ok) {
    pass_ct <<- pass_ct + 1L
    cat(sprintf("PASS  %-60s exp=%s act=%s\n", label, expected, actual))
  } else {
    fail_ct <<- fail_ct + 1L
    msg <- sprintf("FAIL  %-60s expected=%s actual=%s", label, expected, actual)
    fails[[length(fails) + 1L]] <<- msg
    cat(msg, "\n")
  }
  invisible(ok)
}

panel_check <- function(panel_letter, slug, expected) {
  long_path <- file.path(ROOT, sprintf("data/panel%s_%s_long.tsv", panel_letter, slug))
  desc_path <- file.path(ROOT, sprintf("data/panel%s_%s_descriptives.tsv", panel_letter, slug))
  pair_path <- file.path(ROOT, sprintf("data/panel%s_%s_pairwise.tsv", panel_letter, slug))
  cat(sprintf("\n=== PANEL %s — %s ===\n", panel_letter, long_path))
  d <- fread(long_path)
  desc <- fread(desc_path)
  numcols <- names(d)[vapply(d, is.numeric, logical(1))]
  vcol <- setdiff(numcols, "group")[1]
  cat("value column:", vcol, "\n")
  d[, (vcol) := as.numeric(get(vcol))]
  recomp <- d[, .(n = .N,
                  median = as.numeric(median(get(vcol), na.rm = TRUE)),
                  q25 = as.numeric(quantile(get(vcol), 0.25, na.rm = TRUE)),
                  q75 = as.numeric(quantile(get(vcol), 0.75, na.rm = TRUE))), by = group]
  cat("Recomputed:\n"); print(recomp)
  cat("Stored descriptives:\n"); print(desc)
  for (g in names(expected)) {
    e <- expected[[g]]
    r <- recomp[group == g]
    if (nrow(r) == 0) {
      fail_ct <<- fail_ct + 1L
      fails[[length(fails) + 1L]] <<- sprintf("FAIL  Panel %s group %s missing", panel_letter, g)
      next
    }
    check(sprintf("Panel %s %s long-n", panel_letter, g), e$n, r$n)
    check(sprintf("Panel %s %s long-median", panel_letter, g), e$median, r$median, tol = 0.51)
    check(sprintf("Panel %s %s long-q25", panel_letter, g), e$q25, r$q25, tol = 0.51)
    check(sprintf("Panel %s %s long-q75", panel_letter, g), e$q75, r$q75, tol = 0.51)
    ds <- desc[group == g]
    check(sprintf("Panel %s %s desc-n", panel_letter, g), e$n, ds$n)
    check(sprintf("Panel %s %s desc-median", panel_letter, g), e$median, ds$median, tol = 0.51)
    check(sprintf("Panel %s %s desc-q25", panel_letter, g), e$q25, ds$q25, tol = 0.51)
    check(sprintf("Panel %s %s desc-q75", panel_letter, g), e$q75, ds$q75, tol = 0.51)
  }
}

# Section A panels (n=72/118/190)
panel_check("A", "5utr_length", list(
  `NMD+/PTC+` = list(n = 72,  median = 84.5,  q25 = 36.5,   q75 = 167.75),
  `NMD+/PTC-` = list(n = 118, median = 426.5, q25 = 264.25, q75 = 730.5),
  `Control`   = list(n = 190, median = 138,   q25 = 60,     q75 = 266)
))
panel_check("B", "longest_5utr_orf", list(
  `NMD+/PTC+` = list(n = 72,  median = 0,     q25 = 0,      q75 = 33.75),
  `NMD+/PTC-` = list(n = 118, median = 223.5, q25 = 150.75, q75 = 367.5),
  `Control`   = list(n = 190, median = 15,    q25 = 0,      q75 = 101.25)
))
# Section C panels (n=1050/113/1166 after 1:1 gene-matched re-intersection)
panel_check("C", "5utr_length", list(
  `NMD+/PTC+` = list(n = 1050, median = 242,   q25 = 134,   q75 = 421.75),
  `NMD+/PTC-` = list(n = 113,  median = 488,   q25 = 245,   q75 = 945),
  `Control`   = list(n = 1166, median = 200.5, q25 = 99,    q75 = 380.5)
))
panel_check("D", "longest_5utr_orf", list(
  `NMD+/PTC+` = list(n = 1050, median = 36,  q25 = 0,    q75 = 174),
  `NMD+/PTC-` = list(n = 113,  median = 198, q25 = 75,   q75 = 489),
  `Control`   = list(n = 1166, median = 24,  q25 = 0,    q75 = 141)
))

cat("\n=================================================\n")
cat(sprintf("PASS 1 SUMMARY: %d passes / %d failures\n", pass_ct, fail_ct))
cat("=================================================\n")
if (fail_ct > 0) {
  cat("\nFailures:\n")
  for (m in fails) cat(" -", m, "\n")
  quit(status = 1)
}
