#!/usr/bin/env Rscript
# Cross-check verifier — compare new GENCODE-scope Rmd HTML to figure-side
# canonical TSVs to confirm both compute pipelines produce identical numbers.
#
# Per rollback policy (v4 plan §S4): on any FAIL, exit non-zero. Do NOT
# auto-fix the divergence — diagnose root cause first.
#
# Usage:
#   Rscript reproducibility/verify_cross_check_new_rmd_vs_figures.R
#
# Tolerance policy:
#   - counts:      exact
#   - medians:     |obs - exp| <= 0.5
#   - percentages: |obs - exp| <= 0.1
#   - p-values:    same order of magnitude (±1)

suppressPackageStartupMessages({ library(data.table) })

NMD_ROOT <- "/Users/petecastaldi/claude_projects/nmd"
RMD_HTML <- file.path(NMD_ROOT,
  "results/isoform_transitions/Version_6.0/isopair_wrapper",
  "05_final_report_gencode_scope_2026-06-15.html")

FIG3 <- file.path(NMD_ROOT, "figures/multipanel/figure3_isopair_and_ptc/data")
FIG4 <- file.path(NMD_ROOT, "figures/multipanel/figure4_ptcneg_and_model/data")
SUPP_CDS <- file.path(NMD_ROOT,
  "figures/SupplementalFigures/CDSand3UTR_GENCODEonly/data")
SUPP_TD2 <- file.path(NMD_ROOT,
  "figures/SupplementalFigures/TD2BiasEvidence/data")

stopifnot(file.exists(RMD_HTML))

raw  <- readLines(RMD_HTML, warn = FALSE)
text <- gsub("<[^>]+>", " ", paste(raw, collapse = " "))
text <- gsub("[[:space:]]+", " ", text)

# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────
parse_num <- function(s) as.numeric(gsub(",", "", s))

# Escape regex special chars present in group labels (NMD+/PTC+ etc.)
esc <- function(s) {
  s <- gsub("+", "\\+", s, fixed = TRUE)
  s <- gsub("-", "\\-", s, fixed = TRUE)
  s
}

# Extract n + median + q25 + q75 for a group, from the table that follows
# a caption anchor in the HTML stream.
extract_row <- function(anchor, group) {
  pat <- sprintf("%s.*?%s ([0-9,]+) ([0-9.]+) ([0-9.]+) ([0-9.]+)",
                  anchor, esc(group))
  m <- regmatches(text, regexpr(pat, text, perl = TRUE))
  if (length(m) == 0) return(list(n = NA, median = NA, q25 = NA, q75 = NA))
  cap <- regmatches(m, regexec(pat, m, perl = TRUE))[[1]]
  list(n = parse_num(cap[2]), median = parse_num(cap[3]),
       q25 = parse_num(cap[4]), q75 = parse_num(cap[5]))
}

extract_one <- function(regex) {
  m <- regmatches(text, regexpr(regex, text, perl = TRUE))
  if (length(m) == 0) return(NA_real_)
  cap <- regmatches(m, regexec(regex, m, perl = TRUE))[[1]]
  if (length(cap) < 2) return(NA_real_)
  parse_num(cap[2])
}

cmp <- function(obs, exp, type = "median",
                 tol_pct = 0.1, tol_median = 0.5) {
  if (is.na(obs)) return(list(status = "FAIL", msg = "not found in Rmd HTML"))
  if (is.na(exp)) return(list(status = "SKIP", msg = "no figure-side value"))
  pass <- switch(type,
    count  = obs == exp,
    pct    = abs(obs - exp) <= tol_pct,
    median = abs(obs - exp) <= tol_median,
    p      = {
      if (obs == 0 || exp == 0) (obs == 0 && exp == 0)
      else abs(floor(log10(abs(obs))) - floor(log10(abs(exp)))) <= 1
    }
  )
  list(status = if (pass) "PASS" else "FAIL",
       msg = sprintf("rmd=%s fig=%s",
                     format(obs, scientific = (type == "p")),
                     format(exp, scientific = (type == "p"))))
}

results <- list()
record <- function(label, status, msg) {
  cat(sprintf("  [%s] %-58s %s\n", status, label, msg))
  results[[length(results) + 1L]] <<- list(label = label,
                                            status = status, msg = msg)
}

# Generic check across a 3-group descriptive panel.
check_descriptives <- function(anchor, fig_tsv, scope_label,
                                groups = c("NMD+/PTC+", "NMD+/PTC-", "Control")) {
  d <- fread(fig_tsv)
  for (g in groups) {
    exp_row <- d[group == g]
    rmd_row <- extract_row(anchor, g)
    r <- cmp(rmd_row$n, exp_row$n, "count")
    record(sprintf("%s n[%s]", scope_label, g), r$status, r$msg)
    r <- cmp(rmd_row$median, exp_row$median, "median")
    record(sprintf("%s median[%s]", scope_label, g), r$status, r$msg)
  }
}

# ─────────────────────────────────────────────────────────────────────
cat("================================================================\n")
cat("Cross-check: new Rmd vs figure-side TSVs\n")
cat("Rmd :  ", RMD_HTML, "\n", sep = "")
cat("================================================================\n\n")

# 1. Figure 4 Panel A — Section A 5'UTR length
cat("--- Figure 4 Panel A (Section A, 5'UTR length) ---\n")
check_descriptives(
  anchor = "5.UTR length .nt. . Section A descriptives",
  fig_tsv = file.path(FIG4, "panelA_5utr_length_descriptives.tsv"),
  scope_label = "Panel A 5'UTR")

# 2. Figure 4 Panel B — Section A longest 5'UTR ORF
cat("\n--- Figure 4 Panel B (Section A, longest 5'UTR ORF) ---\n")
check_descriptives(
  anchor = "Longest 5.UTR ORF .nt. . Section A descriptives",
  fig_tsv = file.path(FIG4, "panelB_longest_5utr_orf_descriptives.tsv"),
  scope_label = "Panel B 5'UTR-ORF")

# 3. Figure 4 Panel C — Section C 5'UTR length
cat("\n--- Figure 4 Panel C (Section C, 5'UTR length) ---\n")
check_descriptives(
  anchor = "5.UTR length .nt. . Section C descriptives",
  fig_tsv = file.path(FIG4, "panelC_5utr_length_descriptives.tsv"),
  scope_label = "Panel C 5'UTR")

# 4. Figure 4 Panel D — Section C longest 5'UTR ORF
cat("\n--- Figure 4 Panel D (Section C, longest 5'UTR ORF) ---\n")
check_descriptives(
  anchor = "Longest 5.UTR ORF .nt. . Section C descriptives",
  fig_tsv = file.path(FIG4, "panelD_longest_5utr_orf_descriptives.tsv"),
  scope_label = "Panel D 5'UTR-ORF")

# 5. CDSand3UTR supplement — CDS length
cat("\n--- Supplement CDSand3UTR (CDS length) ---\n")
check_descriptives(
  anchor = "CDS length .nt.",
  fig_tsv = file.path(SUPP_CDS, "panelA_cds_length_descriptives.tsv"),
  scope_label = "CDS length")

# 6. CDSand3UTR supplement — 3'UTR translation-based
cat("\n--- Supplement CDSand3UTR (3'UTR translation-based) ---\n")
check_descriptives(
  anchor = "translation-based .own stop . tx end.",
  fig_tsv = file.path(SUPP_CDS, "panelC_utr3_translation_descriptives.tsv"),
  scope_label = "3'UTR translation")

# 7. CDSand3UTR supplement — 3'UTR bias-corrected
cat("\n--- Supplement CDSand3UTR (3'UTR bias-corrected) ---\n")
check_descriptives(
  anchor = "bias-corrected .non-PTC-stop based.",
  fig_tsv = file.path(SUPP_CDS, "panelD_utr3_non_ptc_stop_descriptives.tsv"),
  scope_label = "3'UTR bias-corrected")

# 8. Figure 3 Panel F — mechanism breakdown totals
cat("\n--- Figure 3 Panel F (mechanism breakdown) ---\n")
panelF <- fread(file.path(FIG3, "panelF_mechanism_breakdown.tsv"))
mech_totals_fig <- panelF[, .(n = sum(n)), by = mechanism]
for (m in c("Frameshift", "In-frame stop", "3'UTR splice")) {
  exp_n <- mech_totals_fig[mechanism == m, n]
  obs <- extract_one(sprintf("%s ([0-9]+) ",
                              gsub("'", "[’']", m, fixed = FALSE)))
  r <- cmp(obs, exp_n, "count")
  record(sprintf("Mechanism count[%s]", m), r$status, r$msg)
}

# 9. TD2BiasEvidence supplement headlines
cat("\n--- TD2BiasEvidence supplement summary ---\n")
td2_checks <- list(
  list(label="broad n",            regex="broad = ([0-9,]+) pairs",                                                        exp=1166, type="count"),
  list(label="broad TD2/ref ratio", regex="Broad 1,166 [0-9,.]+ [0-9,]+ ([0-9.]+)×",                                       exp=4.7,  type="median", tol=0.05),
  list(label="broad Kozak n",      regex="Broad [0-9,]+ ([0-9,]+) [0-9,]+ [0-9.]+%",                                      exp=446,  type="count"),
  list(label="broad Kozak %",      regex="Broad [0-9,]+ [0-9,]+ [0-9,]+ ([0-9.]+)% [0-9.]+% [0-9.eE+-]+ Occult",          exp=38.3, type="pct"),
  list(label="broad Kozak p",      regex="Broad [0-9,]+ [0-9,]+ [0-9,]+ [0-9.]+% [0-9.]+% ([0-9.eE+-]+) Occult",          exp=1.21e-38, type="p"),
  list(label="broad % downstream", regex="Broad 1,166 [0-9,]+ [0-9,]+ [0-9,]+ ([0-9.]+)%",                                 exp=42.5, type="pct"),
  list(label="occult n",           regex="Occult-PTC ([0-9,]+) [0-9,.]+ [0-9.]+ [0-9.]+×",                                 exp=492,  type="count"),
  list(label="occult TD2/ref ratio", regex="Occult-PTC [0-9,]+ [0-9,.]+ [0-9.]+ ([0-9.]+)×",                               exp=7.7,  type="median", tol=0.05),
  list(label="occult Kozak n",     regex="Occult-PTC 492 ([0-9,]+) [0-9,]+ [0-9.]+%",                                     exp=384,  type="count"),
  list(label="occult Kozak %",     regex="Occult-PTC 492 384 [0-9]+ ([0-9.]+)%",                                          exp=78.0, type="pct"),
  list(label="occult Kozak p",     regex="Occult-PTC 492 384 [0-9]+ [0-9.]+% [0-9.]+% ([0-9.eE+-]+)",                     exp=2.6e-38, type="p"),
  list(label="occult % downstream", regex="Occult-PTC 492 0 487 5 ([0-9.]+)%",                                             exp=99.0, type="pct")
)
for (row in td2_checks) {
  obs <- extract_one(row$regex)
  r <- cmp(obs, row$exp, row$type,
           tol_median = if (!is.null(row$tol)) row$tol else 0.5)
  record(sprintf("[summary] %s", row$label), r$status, r$msg)
}

# ─────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────
n_pass <- sum(vapply(results, function(r) r$status == "PASS", logical(1)))
n_fail <- sum(vapply(results, function(r) r$status == "FAIL", logical(1)))
n_skip <- sum(vapply(results, function(r) r$status == "SKIP", logical(1)))

cat("\n────────────────────────────────────────────────────────────────\n")
cat(sprintf("Summary: %d PASS / %d FAIL / %d SKIP out of %d checks\n",
            n_pass, n_fail, n_skip, length(results)))
cat("────────────────────────────────────────────────────────────────\n")

if (n_fail > 0L) {
  cat("\nROLLBACK POLICY (v4 plan §S4):\n")
  cat("Do NOT auto-fix the divergence. Diagnose root cause first — typical\n")
  cat("causes: (1) Rmd inline code differs from figure-side scope-construction;\n")
  cat("(2) figure TSV is stale (re-run figure data_export.R); (3) the Rmd was\n")
  cat("knit before the latest Rmd edit. Resolve, then re-run this verifier.\n")
  quit(status = 1L)
}

quit(status = 0L)
