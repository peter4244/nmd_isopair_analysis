#!/usr/bin/env Rscript
# [LEGACY — DEPRECATED] Verification Pass 6 on the legacy Rmd
# (05_final_report_mashr.html, now legacy-bannered as of 2026-06-15).
#
# This verifier targets the PRIOR-SUBMISSION scope. The current canonical
# scope analyses are in
#   results/isoform_transitions/Version_6.0/isopair_wrapper/
#     05_final_report_gencode_scope_2026-06-15.Rmd
# and are validated by:
#   reproducibility/verify_pass7_new_rmd.R               (manifest, 37 checks)
#   reproducibility/verify_cross_check_new_rmd_vs_figures.R  (cross-check, 57)
#
# Run this script ONLY for sensitivity analyses against the prior submission.
# For routine validation of the current scope, use the two verifiers above.
#
# Moved here from figures/multipanel/figure3_isopair_and_ptc/ on 2026-06-15
# (v4 plan §S3) so a verifier of the Rmd no longer lives in a figure folder.
#
# ─────────────────────────────────────────────────────────────────────
# Original header below
# ─────────────────────────────────────────────────────────────────────
# Verification Pass 6 — 5-pass protocol applied to the Rmd source analyses.
# Focused on the values that Figure 3 depends on.

suppressPackageStartupMessages({
  library(data.table)
})

RMD_DIR <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper"
HTML <- file.path(RMD_DIR, "05_final_report_mashr.html")
TABLES <- file.path(RMD_DIR, "tables")
DM <- file.path(RMD_DIR, "data_mashr")
METHODS <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/METHODS.md"

results <- list()
note <- function(pass, label, status, detail = "") {
  cat(sprintf("  [%s] %s%s\n", status, label,
              if (nzchar(detail)) sprintf(" -- %s", detail) else ""))
  results[[length(results) + 1]] <<- list(pass = pass, label = label,
                                          status = status, detail = detail)
}

# ============================================================================
# Sub-pass 6.1 — Factual accuracy: Rmd-side numbers vs raw data
# ============================================================================
cat("================================================================\n")
cat("PASS 6.1 — Factual accuracy (Rmd source analyses)\n")
cat("================================================================\n")

html <- readLines(HTML, warn = FALSE) |> paste(collapse = "\n")

# Flowchart: '3,099 genes' eligibility node
note(6.1, "Flowchart eligible-genes node renders '3,099'",
     if (grepl("3,099", html)) "PASS" else "FAIL")
note(6.1, "Flowchart C2 NMD pairs render '3,009'",
     if (grepl("3,009", html)) "PASS" else "FAIL")
note(6.1, "Flowchart C4 Control pairs render '8,323'",
     if (grepl("8,323", html)) "PASS" else "FAIL")
note(6.1, "Manuscript-stale '3,674' should not appear in current HTML",
     if (!grepl("3,674", html)) "PASS" else "FAIL")

# New ref-AUG-traceable section: 83.5%, 16.3%, etc.
note(6.1, "New canonical section renders 'Reference-AUG-Traceable'",
     if (grepl("Reference-AUG-Traceable Canonical Analysis", html)) "PASS" else "FAIL")
note(6.1, "New canonical section renders NMD PTC rate '83.5'",
     if (grepl("83\\.5", html)) "PASS" else "FAIL")
note(6.1, "New canonical section renders Control PTC rate '16.3'",
     if (grepl("16\\.3", html)) "PASS" else "FAIL")
note(6.1, "New canonical section renders fold enrichment '5.1'",
     if (grepl("5\\.1[^0-9]", html)) "PASS" else "FAIL")
note(6.1, "New canonical section renders pop_traceable NMD count '2,289'",
     if (grepl("2,289", html)) "PASS" else "FAIL")
note(6.1, "New canonical section renders pop_ptc_plus NMD count '1,912'",
     if (grepl("1,912", html)) "PASS" else "FAIL")

# Tables
t2 <- fread(file.path(TABLES, "table2_ptc_rates_allsamples.csv"))
note(6.1, "table2 NMD denominator = 2,496",
     if (t2$n_coding_pairs_NMD == 2496) "PASS" else "FAIL")
note(6.1, "table2 NMD PTC = 1,068 / 42.8%",
     if (t2$NMD_comp_ptc_n == 1068 && t2$NMD_comp_ptc_pct == 42.8) "PASS" else "FAIL")
note(6.1, "table2 Control PTC = 121 / 4.8%",
     if (t2$Ctrl_comp_ptc_n == 121 && t2$Ctrl_comp_ptc_pct == 4.8) "PASS" else "FAIL")

t10 <- fread(file.path(TABLES, "table10_ref_aug_traceable_summary.csv"))
note(6.1, "table10 NMD traceable=2289 PTC+=1912 (83.5%)",
     if (t10$n_NMD_traceable == 2289 && t10$n_NMD_ptc_plus == 1912 &&
         t10$pct_NMD_ptc == 83.5) "PASS" else "FAIL")
note(6.1, "table10 Control traceable=1763 PTC+=288 (16.3%)",
     if (t10$n_Ctrl_traceable == 1763 && t10$n_Ctrl_ptc_plus == 288 &&
         t10$pct_Ctrl_ptc == 16.3) "PASS" else "FAIL")
note(6.1, "table10 fold enrichment=5.1 OR=25.94",
     if (t10$fold_enrichment == 5.1 && t10$fisher_OR == 25.94) "PASS" else "FAIL")

t5b <- fread(file.path(TABLES, "table5b_ptc_event_enrichment.csv"))
se_row <- t5b[event_type == "SE"]
note(6.1, "table5b SE PTC events=520, ctrl events=698 (TD2-only scope)",
     if (se_row$n_ptc_events == 520 && se_row$n_ctrl_events == 698) "PASS" else "FAIL")

# ============================================================================
# Sub-pass 6.2 — Result correctness: Fisher tests reproduce
# ============================================================================
cat("\n================================================================\n")
cat("PASS 6.2 — Result correctness (Rmd source analyses)\n")
cat("================================================================\n")

# table2 Fisher OR reconstruction
a <- t2$NMD_comp_ptc_n; b <- t2$n_coding_pairs_NMD - a
c_ <- t2$Ctrl_comp_ptc_n; d <- t2$n_coding_pairs_Ctrl - c_
ft <- fisher.test(matrix(c(a, b, c_, d), nrow = 2))
or_re <- round(unname(ft$estimate), 2)
note(6.2, "table2 Fisher OR reconstruction matches",
     if (abs(or_re - t2$fisher_OR) < 0.05) "PASS" else "FAIL",
     sprintf("recompute=%.2f tsv=%.2f", or_re, t2$fisher_OR))

# table10 Fisher OR reconstruction
a <- t10$n_NMD_ptc_plus; b <- t10$n_NMD_traceable - a
c_ <- t10$n_Ctrl_ptc_plus; d <- t10$n_Ctrl_traceable - c_
ft <- fisher.test(matrix(c(a, b, c_, d), nrow = 2))
or_re <- round(unname(ft$estimate), 2)
note(6.2, "table10 Fisher OR reconstruction matches",
     if (abs(or_re - t10$fisher_OR) < 0.05) "PASS" else "FAIL",
     sprintf("recompute=%.2f tsv=%.2f", or_re, t10$fisher_OR))

# Independent calculation: NMD PTC rate vs ref_atg_analysis raw counts
ref_atg <- readRDS(file.path(DM, "analysis_cache/ref_atg_analysis.rds"))
trace <- c("effectively_ptc", "no_downstream_ejc", "truncated_no_ejc")
n_trace_c2 <- sum(ref_atg$c2$category %in% trace)
n_eff_c2 <- sum(ref_atg$c2$category == "effectively_ptc")
note(6.2, "table10 NMD numbers match raw ref_atg data",
     if (n_trace_c2 == 2289 && n_eff_c2 == 1912) "PASS" else "FAIL",
     sprintf("ref_atg trace=%d eff=%d", n_trace_c2, n_eff_c2))

# table5b SE OR reconstruction (uses event-level denominator within table5b)
total_ptc <- sum(t5b$n_ptc_events); total_ctrl <- sum(t5b$n_ctrl_events)
a <- se_row$n_ptc_events; b <- total_ptc - a
c_ <- se_row$n_ctrl_events; d <- total_ctrl - c_
ft_se <- fisher.test(matrix(c(a, b, c_, d), nrow = 2))
note(6.2, "table5b SE direction matches",
     if (se_row$direction == "Enriched in PTC-causing" && ft_se$estimate > 1)
       "PASS" else "FAIL",
     sprintf("OR=%.2f", ft_se$estimate))

# ============================================================================
# Sub-pass 6.3 — Documentation accuracy
# ============================================================================
cat("\n================================================================\n")
cat("PASS 6.3 — Documentation accuracy (Rmd source analyses)\n")
cat("================================================================\n")

# Verify inline R values in new section match table10
# Check that "Headline numbers under the canonical Figure 3 scope" section renders correctly
canon_section <- regmatches(html, regexpr("Headline numbers under the canonical[^h]{0,2000}", html, perl = TRUE))
note(6.3, "New canonical-section headline prose found in HTML",
     if (length(canon_section) > 0 && nchar(canon_section) > 100) "PASS" else "FAIL")
if (length(canon_section) > 0) {
  has_2289 <- grepl("2,289", canon_section)
  has_83 <- grepl("83\\.5", canon_section)
  has_16 <- grepl("16\\.3", canon_section)
  has_51 <- grepl("5\\.1", canon_section)
  note(6.3, "All canonical headline numbers in prose section",
       if (has_2289 && has_83 && has_16 && has_51) "PASS" else "FAIL",
       sprintf("2289=%s 83=%s 16=%s 51=%s",
               has_2289, has_83, has_16, has_51))
}

# Verify Cross-reference to METHODS.md and figures/multipanel
note(6.3, "Rmd new section cross-references figures/multipanel directory",
     if (grepl("figures/multipanel/figure3_isopair_and_ptc", html)) "PASS" else "FAIL")
note(6.3, "Rmd new section cross-references METHODS.md",
     if (grepl("Reference-AUG-Traceable Canonical PTC Analysis", html)) "PASS" else "FAIL")

# ============================================================================
# Sub-pass 6.4 — Reproducibility
# ============================================================================
cat("\n================================================================\n")
cat("PASS 6.4 — Reproducibility (Rmd source analyses)\n")
cat("================================================================\n")

html_info <- file.info(HTML)
fresh_html <- as.POSIXct(html_info$mtime) > as.POSIXct("2026-06-13 14:00")
note(6.4, "HTML render completed successfully (file fresh)",
     if (fresh_html) "PASS" else "FAIL",
     sprintf("mtime=%s", format(html_info$mtime, "%Y-%m-%d %H:%M:%S")))

# Check that critical tables for Figure 3 are fresh
critical_tables <- c("table2_ptc_rates_allsamples.csv",
                     "table5b_ptc_event_enrichment.csv",
                     "table10_ref_aug_traceable_summary.csv")
for (tname in critical_tables) {
  t_info <- file.info(file.path(TABLES, tname))
  fresh <- as.POSIXct(t_info$mtime) > as.POSIXct("2026-06-13 14:00")
  note(6.4, sprintf("Table %s is fresh", tname),
       if (fresh) "PASS" else "FAIL")
}

# Check that new chunk exists in Rmd
rmd <- readLines(file.path(RMD_DIR, "05_final_report_mashr.Rmd"))
chunk_idx <- grep("sec2c-ref-aug-canonical", rmd)
note(6.4, "New sec2c-ref-aug-canonical chunk exists in Rmd",
     if (length(chunk_idx) > 0) "PASS" else "FAIL")

# ============================================================================
# Sub-pass 6.5 — METHODS.md vs Rmd code consistency
# ============================================================================
cat("\n================================================================\n")
cat("PASS 6.5 — METHODS.md vs Rmd consistency\n")
cat("================================================================\n")

methods <- readLines(METHODS)
note(6.5, "METHODS.md has new 'Reference-AUG-Traceable Canonical PTC Analysis' section",
     if (any(grepl("Reference-AUG-Traceable Canonical PTC Analysis", methods)))
       "PASS" else "FAIL")
note(6.5, "METHODS.md has 'Pair Construction and Gene-Matching — Detailed Filter Chain' section",
     if (any(grepl("Detailed Filter Chain", methods))) "PASS" else "FAIL")
note(6.5, "METHODS.md correctly attributes TD2 for novel isoforms (not 'SQANTI ORF')",
     if (!any(grepl("SQANTI ORF predictions \\(novel", methods)) &&
         any(grepl("TransDecoder2 \\(TD2\\)", methods)))
       "PASS" else "FAIL")
note(6.5, "METHODS.md category list includes 'mapping_failed'",
     if (any(grepl("mapping_failed", methods))) "PASS" else "FAIL")
note(6.5, "METHODS.md notes 'same_or_longer_with_ejc' is no longer emitted",
     if (any(grepl("same_or_longer_with_ejc.*does not emit|does not emit.*same_or_longer", methods)))
       "PASS" else "FAIL")
note(6.5, "METHODS.md documents 83.5% / 16.3% headline numbers",
     if (any(grepl("83\\.5", methods)) && any(grepl("16\\.3", methods)))
       "PASS" else "FAIL")
note(6.5, "METHODS.md documents the mixed-source attribution chain",
     if (any(grepl("Mixed-source PTC attribution", methods))) "PASS" else "FAIL")

# Verify all Rmd-cited METHODS.md sections exist (partial check)
rmd_methods_refs <- grep("METHODS\\.md", rmd, value = TRUE)
cited_sections <- unique(gsub('.*METHODS\\.md § "([^"]+)".*', "\\1", rmd_methods_refs))
cited_sections <- cited_sections[!grepl("METHODS\\.md", cited_sections)]
# Strip parenthetical qualifiers and "Publication Report Methods >" prefix for matching
strip <- function(s) {
  s <- gsub("Publication Report Methods > ", "", s, fixed = TRUE)
  s <- gsub(" \\([^)]+\\)", "", s)
  s <- gsub(" > ", " ", s)
  trimws(s)
}
cited_simple <- sapply(cited_sections, strip)
missing <- character(0)
for (s in unique(cited_simple)) {
  if (!any(grepl(s, methods, fixed = TRUE))) missing <- c(missing, s)
}
note(6.5, "All Rmd-cited METHODS.md sections exist (lenient match)",
     if (length(missing) == 0) "PASS" else "WARN",
     if (length(missing) > 0) paste("missing:", paste(missing, collapse = "; ")) else "")

# ============================================================================
# Summary
# ============================================================================
cat("\n================================================================\n")
n_pass <- sum(sapply(results, function(r) r$status == "PASS"))
n_warn <- sum(sapply(results, function(r) r$status == "WARN"))
n_fail <- sum(sapply(results, function(r) r$status == "FAIL"))
cat(sprintf("Pass 6 (5 sub-passes) summary: %d PASS / %d WARN / %d FAIL / %d total\n",
            n_pass, n_warn, n_fail, length(results)))
by_pass <- table(sapply(results, function(r) r$pass),
                 sapply(results, function(r) r$status))
print(by_pass)
if (n_fail > 0) {
  cat("\nFAILURES:\n")
  for (r in results) {
    if (r$status == "FAIL") {
      cat(sprintf("  - [%s] %s -- %s\n", r$pass, r$label, r$detail))
    }
  }
}
cat("================================================================\n")
