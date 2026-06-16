#!/usr/bin/env Rscript
# Pass 5: METHODS.md verification — check stale TD2 supplement claims at n=900.

METHODS <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/METHODS.md"
lines <- readLines(METHODS)

stale_markers <- list(
  list("n = 900",       "TD2 supplement claimed scope (current 1166/492)"),
  list("n=900",         "TD2 supplement claimed scope (current 1166/492)"),
  list("Median TD2 ORF = 3,168", "Pre-rebuild ORF length (current TD2 broad median = 2,758, occult 2,934)"),
  list("Median ref AUG PWM = 1.09", "Pre-rebuild Kozak (current broad 0.86 / occult 1.04)"),
  list("8.5", "Pre-rebuild length ratio (current 4.69x broad / 7.73x occult)"),
  list("889 of 900", "Pre-rebuild downstream count (current 487/492 occult)"),
  list("+484 nt", "Pre-rebuild median offset (current +476 nt)"),
  list("3.7", "Pre-rebuild Kozak p-value (current 1.2e-38 broad / 2.6e-38 occult)")
)
cat("=== Searching METHODS.md for stale TD2 supplement markers ===\n")
for (m in stale_markers) {
  idx <- grep(m[[1]], lines, fixed = TRUE)
  if (length(idx) > 0) {
    cat(sprintf("  FOUND  '%s'  (%s)\n    lines: %s\n",
                m[[1]], m[[2]], paste(idx, collapse = ", ")))
  } else {
    cat(sprintf("  OK     '%s' not present\n", m[[1]]))
  }
}

cat("\n=== Searching METHODS.md for CDSand3UTR supplement references ===\n")
cds_refs <- grep("CDSand3UTR|cds_and_3utr", lines)
if (length(cds_refs) == 0) {
  cat("  GAP: CDSand3UTR_GENCODEonly supplement is NOT mentioned anywhere in METHODS.md.\n")
} else {
  cat(sprintf("  FOUND CDS supplement refs at lines: %s\n",
              paste(cds_refs, collapse = ", ")))
}
