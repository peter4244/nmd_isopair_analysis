#!/usr/bin/env Rscript
# Verification Pass 7 — new GENCODE-scope Rmd against expected-value manifest.
#
# Validates that
#   results/isoform_transitions/Version_6.0/isopair_wrapper/
#     05_final_report_gencode_scope_2026-06-15.html
# renders the canonical numbers established during Phase 3 of the v4 plan
# (paper/rmd_update_plan_v4_2026-06-15.md).
#
# Tolerance policy (per S2 of the v4 plan):
#   - counts:      exact match
#   - percentages: |obs - exp| <= 0.1
#   - medians:     |obs - exp| <= 0.5
#   - p-values:    same order of magnitude
#                  (i.e. floor(log10(obs)) == floor(log10(exp)) within ±1)
#
# Rollback policy (per S4): on any FAIL, exit non-zero. Do NOT modify either
# the Rmd or the figure pipeline to silence a discrepancy — diagnose root
# cause first.
#
# Usage:
#   Rscript reproducibility/verify_pass7_new_rmd.R
#
# Exit codes: 0 if all PASS; 1 if any FAIL.

suppressPackageStartupMessages({ library(data.table) })

NMD_ROOT <- "/Users/petecastaldi/claude_projects/nmd"
RMD_HTML <- file.path(NMD_ROOT,
  "results/isoform_transitions/Version_6.0/isopair_wrapper",
  "05_final_report_gencode_scope_2026-06-15.html")

if (!file.exists(RMD_HTML)) {
  stop(sprintf("Rmd HTML not found at %s\nRun rmarkdown::render() first.", RMD_HTML))
}

raw  <- readLines(RMD_HTML, warn = FALSE)
text <- gsub("<[^>]+>", " ", paste(raw, collapse = " "))
text <- gsub("[[:space:]]+", " ", text)

# ─────────────────────────────────────────────────────────────────────
# Manifest of expected values
# ─────────────────────────────────────────────────────────────────────
# Each row: name, type ("count" / "pct" / "median" / "p"),
#           regex (must capture the numeric value in group 1),
#           expected, label (for the report).
manifest <- rbindlist(list(
  # §1 — pop_BC
  data.table(name="pop_BC_c2_count", type="count",
    regex="pop_BC = ([0-9,]+) shared",
    expected=3009, label="§1 pop_BC shared pairs"),
  data.table(name="pop_BC_c2_n",      type="count",
    regex="([0-9,]+) NMD-comp pairs",
    expected=3009, label="§1 pop_BC NMD pairs"),
  data.table(name="pop_BC_c4_n",      type="count",
    regex="([0-9,]+) Control-comp pairs",
    expected=3009, label="§1 pop_BC Control pairs"),
  data.table(name="tx_len_ref",       type="median",
    regex="Reference = ([0-9,]+)",
    expected=2893, label="§1 transcript-length median Reference"),
  data.table(name="tx_len_nmd",       type="median",
    regex="NMD comparator = ([0-9,]+)",
    expected=3049, label="§1 transcript-length median NMD comp"),
  data.table(name="tx_len_ctrl",      type="median",
    regex="Control comparator = ([0-9,]+)",
    expected=2762, label="§1 transcript-length median Control comp"),
  data.table(name="ref_pct_enst",     type="pct",
    regex="([0-9]+\\.[0-9]+)% ENST,",
    expected=69.7, label="§1 Reference %ENST"),
  data.table(name="nmd_pct_novel",    type="pct",
    regex="predominantly novel \\(([0-9]+\\.[0-9]+)%",
    expected=76.6, label="§1 NMD comp %Novel"),

  # §2 scope
  data.table(name="gencode_all3_c2",  type="count",
    regex="gencode_all3_c2 = ([0-9,]+) NMD pairs",
    expected=190, label="§2 scope gencode_all3_c2"),
  data.table(name="gencode_all3_c4",  type="count",
    regex="gencode_all3_c4 = ([0-9,]+) Control pairs",
    expected=190, label="§2 scope gencode_all3_c4"),

  # §2a PTC determination
  data.table(name="secA_ptc_nmd_n",   type="count",
    regex="NMD: ([0-9]+) / 190",
    expected=72, label="§2a NMD PTC+ count"),
  data.table(name="secA_ptc_ctrl_n",  type="count",
    regex="Control: ([0-9]+) / 190",
    expected=4, label="§2a Control PTC+ count"),
  data.table(name="secA_ptc_nmd_pct", type="pct",
    regex="NMD: 72 / 190 = ([0-9]+\\.[0-9]+)%",
    expected=37.9, label="§2a NMD PTC+ rate"),
  data.table(name="secA_ptc_ctrl_pct", type="pct",
    regex="Control: 4 / 190 = ([0-9]+\\.[0-9]+)%",
    expected=2.1, label="§2a Control PTC+ rate"),
  data.table(name="secA_fold",        type="pct",
    regex="Fold enrichment: ([0-9.]+)",
    expected=18.0, label="§2a fold enrichment"),
  data.table(name="secA_or",          type="median",
    regex="Fisher: OR = ([0-9.]+),",
    expected=28.2, label="§2a Fisher OR"),
  data.table(name="secA_fisher_p",    type="p",
    regex="Fisher: OR = [0-9.]+, p = ([0-9.eE+-]+)",
    expected=1.88e-20, label="§2a Fisher p-value"),

  # §2b attribution
  data.table(name="secA_attr_direct", type="count",
    regex="([0-9]+) direct \\(from",
    expected=66, label="§2b direct attributions"),
  data.table(name="secA_attr_ss",     type="count",
    regex="\\+ ([0-9]+) same-stop",
    expected=8, label="§2b same-stop attributions"),
  data.table(name="secA_attr_unique", type="count",
    regex="= 74 total attribution rows .{0,80}?([0-9]+) unique direct",
    expected=69, label="§2b unique direct attributions"),
  data.table(name="secA_attr_unres",  type="count",
    regex="([0-9]+) unresolved out of 72 PTC. pairs",
    expected=3, label="§2b finally unresolved"),
  data.table(name="secA_attr_cov",    type="pct",
    regex="\\(([0-9]+\\.[0-9]+)% attribution coverage",
    expected=95.8, label="§2b attribution coverage"),

  # §3 scope
  data.table(name="refaug_c2",        type="count",
    regex="refaug_matched_c2 = ([0-9,]+) NMD",
    expected=1166, label="§3 scope refaug_matched_c2"),
  data.table(name="refaug_c4",        type="count",
    regex="refaug_matched_c4 = ([0-9,]+) Control",
    expected=1166, label="§3 scope refaug_matched_c4"),
  data.table(name="secC_pct_ptc",     type="pct",
    regex="PTC. rate at this scope: 1,050 / 1,166 = ([0-9]+\\.[0-9]+)%",
    expected=90.1, label="§3b PTC+ rate at Section C"),

  # §4 TD2 bias — broad scope
  data.table(name="broad_n",          type="count",
    regex="broad = ([0-9,]+) pairs",
    expected=1166, label="§4 broad scope n"),
  data.table(name="td2_broad_kozak_n",type="count",
    regex="Broad [0-9,]+ ([0-9,]+) [0-9,]+ [0-9.]+%",
    expected=446, label="§4b broad ref>TD2 Kozak count"),
  data.table(name="td2_broad_kozak_pct", type="pct",
    regex="Broad [0-9,]+ [0-9,]+ [0-9,]+ ([0-9]+\\.[0-9]+)% [0-9.]+% [0-9.eE+-]+ Occult",
    expected=38.3, label="§4b broad ref>TD2 Kozak % (of all)"),
  data.table(name="td2_broad_kozak_p", type="p",
    regex="Broad [0-9,]+ [0-9,]+ [0-9,]+ [0-9]+\\.[0-9]+% [0-9.]+% ([0-9.eE+-]+) Occult",
    expected=1.21e-38, label="§4b broad Kozak p-value"),
  data.table(name="td2_broad_pct_down", type="pct",
    regex="Broad 1,166 [0-9,]+ [0-9,]+ [0-9,]+ ([0-9]+\\.[0-9]+)%",
    expected=42.5, label="§4c broad % TD2 downstream"),

  # §4 TD2 bias — occult-PTC scope
  data.table(name="occ_n",            type="count",
    regex="Occult-PTC ([0-9,]+) [0-9,.]+ [0-9.]+ [0-9.]+×",
    expected=492, label="§4 occult-PTC scope n"),
  data.table(name="td2_occ_kozak_n",  type="count",
    regex="Occult-PTC 492 ([0-9,]+) [0-9,]+ [0-9.]+%",
    expected=384, label="§4b occult-PTC ref>TD2 Kozak count"),
  data.table(name="td2_occ_kozak_pct",type="pct",
    regex="Occult-PTC 492 384 [0-9]+ ([0-9]+\\.[0-9]+)%",
    expected=78.0, label="§4b occult-PTC ref>TD2 Kozak % (of all)"),
  data.table(name="td2_occ_kozak_p",  type="p",
    regex="Occult-PTC 492 384 [0-9]+ 78\\.[0-9]+% [0-9.]+% ([0-9.eE+-]+)",
    expected=2.60e-38, label="§4b occult-PTC Kozak p-value"),
  data.table(name="td2_occ_pct_down", type="pct",
    regex="Occult-PTC 492 0 487 5 ([0-9]+\\.[0-9]+)%",
    expected=99.0, label="§4c occult-PTC % TD2 downstream"),
  # §4a ORF ratios — parse the ratio column of the table
  data.table(name="td2_broad_ratio",  type="median",
    regex="Broad 1,166 [0-9,.]+ [0-9,]+ ([0-9.]+)×",
    expected=4.69, label="§4a broad TD2/ref ORF ratio"),
  data.table(name="td2_occ_ratio",    type="median",
    regex="Occult-PTC 492 [0-9,.]+ [0-9.]+ ([0-9.]+)×",
    expected=7.73, label="§4a occult-PTC TD2/ref ORF ratio")
))

# ─────────────────────────────────────────────────────────────────────
# Verification
# ─────────────────────────────────────────────────────────────────────
parse_num <- function(s) as.numeric(gsub(",", "", s))

check <- function(row) {
  m <- regmatches(text, regexpr(row$regex, text, perl = TRUE))
  if (length(m) == 0) {
    return(list(status = "FAIL", obs = NA_real_, msg = "value not found in HTML"))
  }
  # Extract the captured group
  cap <- regmatches(m, regexec(row$regex, m, perl = TRUE))[[1]]
  if (length(cap) < 2) {
    return(list(status = "FAIL", obs = NA_real_, msg = "regex matched but no capture"))
  }
  obs <- parse_num(cap[2])
  if (is.na(obs)) {
    return(list(status = "FAIL", obs = NA_real_, msg = sprintf("could not parse '%s'", cap[2])))
  }
  exp <- row$expected
  status <- switch(row$type,
    count  = if (obs == exp) "PASS" else "FAIL",
    pct    = if (abs(obs - exp) <= 0.1) "PASS" else "FAIL",
    median = if (abs(obs - exp) <= 0.5) "PASS" else "FAIL",
    p      = {
      if (obs == 0 || exp == 0) {
        if (obs == 0 && exp == 0) "PASS" else "FAIL"
      } else if (abs(floor(log10(abs(obs))) - floor(log10(abs(exp)))) <= 1) {
        "PASS"
      } else "FAIL"
    }
  )
  list(status = status, obs = obs,
       msg = sprintf("obs=%s exp=%s",
                     format(obs, scientific = (row$type == "p")),
                     format(exp, scientific = (row$type == "p"))))
}

cat("================================================================\n")
cat("PASS 7 — verify new GENCODE-scope Rmd against expected manifest\n")
cat("Source: ", RMD_HTML, "\n", sep = "")
cat("================================================================\n\n")

results <- vector("list", nrow(manifest))
n_pass <- 0L; n_fail <- 0L
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i]
  r <- check(row)
  status_str <- if (r$status == "PASS") "PASS" else "FAIL"
  cat(sprintf("  [%s] %-44s %s\n", status_str, row$label, r$msg))
  results[[i]] <- list(name = row$name, label = row$label,
                        status = r$status, obs = r$obs, exp = row$expected)
  if (r$status == "PASS") n_pass <- n_pass + 1L else n_fail <- n_fail + 1L
}

cat("\n────────────────────────────────────────────────────────────────\n")
cat(sprintf("Summary: %d PASS / %d FAIL out of %d checks\n",
            n_pass, n_fail, nrow(manifest)))
cat("────────────────────────────────────────────────────────────────\n")

if (n_fail > 0L) {
  cat("\nROLLBACK POLICY (per v4 plan §S4):\n")
  cat("Do NOT modify either the Rmd or the figure pipeline to silence the\n")
  cat("discrepancy. Diagnose root cause first. The expected manifest is the\n")
  cat("canonical record from Phase 3 — if you intend to update it, document\n")
  cat("the reason in the commit message and update both this script and the\n")
  cat("Rmd in the same commit.\n")
  quit(status = 1L)
}

quit(status = 0L)
