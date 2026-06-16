#!/usr/bin/env Rscript
# Verification Pass 5 — METHODS.md verification for Figure 3
#
# Cross-check the canonical METHODS.md against the actual code behavior.
# Confirm the new all-3-ENST + own-GENCODE-stop scope is documented; flag
# documentation drift, missing sections, and stale claims.
#
# Scope refresh 2026-06-15 (v2): all-3-ENST + coding-CDS + own-GENCODE-stop.

suppressPackageStartupMessages({
  library(data.table)
})

results <- list()
note <- function(label, status, detail = "") {
  cat(sprintf("  [%s] %s%s\n", status, label,
              if (nzchar(detail)) sprintf(" -- %s", detail) else ""))
  results[[length(results) + 1]] <<- list(label = label, status = status,
                                          detail = detail)
}

METHODS_PATH <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/METHODS.md"
methods <- readLines(METHODS_PATH)

cat("================================================================\n")
cat("PASS 5: METHODS.md verification for Figure 3 (all-3-ENST scope)\n")
cat("================================================================\n\n")

# ---- 1. CDS-source terminology check ----
cat("[1] CDS-source terminology drift\n")
sqanti_orf_mention <- grep("SQANTI ORF predictions", methods, value = TRUE)
note("METHODS line that says 'SQANTI ORF predictions (novel PacBio isoforms)' should say 'TD2'",
     if (length(sqanti_orf_mention) > 0) "GAP" else "PASS",
     "Pete 2026-06-13: TD2 is the actual CDS caller for novel isoforms; SQANTI3 just wraps it.")

td2_mention <- grep("TransDecoder2|TD2", methods, value = TRUE)
note("METHODS does correctly cite TD2 in some places",
     if (length(td2_mention) > 0) "PASS" else "FAIL",
     sprintf("found %d TD2 mentions", length(td2_mention)))

# ---- 2. Reference-ATG Tracing category list completeness ----
# Each actual category emitted by 05r_ref_atg_analysis.R must be mentioned by
# name somewhere in METHODS.md. The historical category `same_or_longer_with_ejc`
# is documented as deprecated (2026-06-13 note) and no longer emitted.
cat("\n[2] Reference-ATG Tracing categories — METHODS vs actual code output\n")
ref_atg <- readRDS("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr/analysis_cache/ref_atg_analysis.rds")
actual_cats <- sort(unique(c(ref_atg$c2$category, ref_atg$c4$category)))
cat(sprintf("  Actual categories in ref_atg_analysis.rds: %s\n",
            paste(actual_cats, collapse = ", ")))
emitted_but_not_documented <- character(0)
for (cat_name in actual_cats) {
  if (!any(grepl(cat_name, methods, fixed = TRUE))) {
    emitted_but_not_documented <- c(emitted_but_not_documented, cat_name)
  }
}
note("All emitted categories documented in METHODS",
     if (length(emitted_but_not_documented) == 0) "PASS" else "GAP",
     if (length(emitted_but_not_documented) > 0)
       paste("missing from METHODS:", paste(emitted_but_not_documented, collapse = ", "))
     else "")
# Check the historical category is acknowledged as deprecated in METHODS
note("Deprecated `same_or_longer_with_ejc` flagged as no-longer-emitted in METHODS",
     if (any(grepl("same_or_longer_with_ejc", methods))) "PASS" else "GAP",
     "Historical category should be explicitly noted as deprecated to avoid re-introduction.")

# ---- 3. New all-3-ENST scope — is it documented? ----
cat("\n[3] All-3-ENST + own-GENCODE-stop canonical scope (Pete 2026-06-15)\n")
note("METHODS has the 'All-3-ENST + Own-GENCODE-Stop Canonical PTC Analysis' section header",
     if (any(grepl("All-3-ENST \\+ Own-GENCODE-Stop Canonical PTC Analysis", methods)))
       "PASS" else "GAP",
     "This is the canonical Figure 3 Panels D/E/F section header.")

note("METHODS documents the all-3-ENST + coding-CDS scope (Pete's 2026-06-15 policy)",
     if (any(grepl("all-3-ENST|all three.*GENCODE.annotated", methods))) "PASS" else "GAP",
     "Pete's policy: Figure 3 Panels D/E/F restrict to pairs where reference + NMD comparator + Control comparator are all GENCODE-annotated coding isoforms; each comparator's own GENCODE-annotated stop drives PTC determination.")

note("METHODS documents the 37.9% NMD / 2.1% Control PTC rate headline",
     if (any(grepl("37\\.9", methods)) && any(grepl("2\\.1%", methods))) "PASS" else "GAP",
     "These are the new manuscript-headline numbers (Figure 3 Panel D).")

note("METHODS documents the 18x enrichment headline",
     if (any(grepl("18.fold|18×|18-fold", methods))) "PASS" else "GAP",
     "Fold enrichment headline for Figure 3 Panel D.")

note("METHODS documents 190 NMD / 190 Control denominator",
     if (any(grepl("190", methods))) "PASS" else "GAP")

note("METHODS documents PTC+ subset = 72 NMD / 4 Control",
     if (any(grepl("72.*4|72 / 190", methods)) || any(grepl("72 / 4", methods)))
       "PASS" else "GAP")

note("METHODS documents n_ptc_attr = 69 attributed events",
     if (any(grepl("69 ", methods)) || any(grepl("69$", methods))) "PASS" else "GAP")

note("METHODS documents ctrl_total = 447",
     if (any(grepl("447", methods))) "PASS" else "GAP")

# ---- 4. Coding-coding pre-filter ----
cat("\n[4] Coding-CDS pre-filter for Figure 3 D/E/F\n")
note("METHODS notes coding-coding restriction for Figure 3 D/E/F",
     if (any(grepl("coding-coding|both isoforms are coding|coding[- ]?CDS|coding_status", methods)))
       "PASS" else "GAP",
     "Important: explains the 301 -> 190 reduction from ENST gene-matched to all-3-ENST + coding-CDS.")

# ---- 5. Own-GENCODE-stop attribution chain (Panel E/F) ----
cat("\n[5] Own-GENCODE-stop PTC attribution chain (Panel E/F)\n")
note("METHODS documents own-GENCODE-stop attribution",
     if (any(grepl("own.*GENCODE.*stop|own-GENCODE-stop|own GENCODE-annotated stop", methods)))
       "PASS" else "GAP",
     "Under the new scope, attribution uses each comparator's own GENCODE-annotated stop position, NOT ref-AUG projection. attribute_ptc_events(ptc_genomic_pos = own_stop_genomic_lookup).")

note("METHODS documents that TD2 is NOT used in the new Figure 3 scope",
     if (any(grepl("no TD2 dependency|No TD2|TD2-free|without TD2", methods, ignore.case = TRUE)))
       "PASS" else "GAP",
     "The new scope is fully TD2-free; this should be explicit in METHODS.")

note("METHODS documents that ref-AUG projection is NOT used in the new Figure 3 scope",
     if (any(grepl("no ref-AUG projection|No ref-AUG projection", methods)))
       "PASS" else "GAP",
     "The new scope uses each isoform's own GENCODE CDS directly, with no ref-AUG projection.")

# ---- 6. Stale numbers / claims in METHODS ----
cat("\n[6] Stale numbers cross-check\n")
# The old headline (83.5% / 16.3% / 5.1×) was for the prior ENST-reference +
# ref-AUG-projected scope. Under the new policy these aren't the Figure 3
# headline rates. METHODS.md may still mention them in the Section C / Figure
# 4 context (which is legitimate), but they should NOT appear as the Figure 3
# Panel D headline.
in_fig3_section <- function(needle) {
  # Look only inside the "All-3-ENST + Own-GENCODE-Stop Canonical PTC Analysis"
  # section block to see if old numbers leaked into the new section.
  start <- grep("All-3-ENST \\+ Own-GENCODE-Stop Canonical PTC Analysis", methods)
  if (length(start) == 0) return(FALSE)
  # End of block = next "###" or "##" heading after `start`.
  rest <- methods[(start[1] + 1):length(methods)]
  end_offset <- min(c(grep("^##? ", rest), length(rest)))
  block <- rest[seq_len(end_offset)]
  any(grepl(needle, block))
}
note("Figure 3 section in METHODS doesn't quote stale 83.5% as the Figure 3 headline",
     if (!in_fig3_section("83\\.5")) "PASS" else "GAP")
note("Figure 3 section in METHODS doesn't quote stale 5.1x as the Figure 3 headline",
     if (!in_fig3_section("5\\.1.fold|5\\.1-fold|5\\.1×")) "PASS" else "GAP")
note("METHODS doesn't cite the old 44% PTC rate (would be stale)",
     if (!any(grepl("44%.*PTC|n_ptc.*44|43\\.9%", methods))) "PASS" else "GAP")
note("METHODS doesn't cite the old 3,674 gene count",
     if (!any(grepl("3,?674.*genes", methods))) "PASS" else "GAP")

# ---- 7. Compute_ptc_rates_row description ----
cat("\n[7] compute_ptc_rates_row() — denominator chain\n")
note("METHODS documents the Stage 2 + coding-coding + gene_id re-match chain",
     if (any(grepl("Stage 2.*gene-match|three-stage filter|Fisher.*denominator|compute_ptc_rates", methods))) "PASS" else "GAP",
     "compute_ptc_rates_row applies a 3-stage filter that was a source of confusion.")

# ---- 8. Figure 3 source-code map ----
cat("\n[8] Figure 3 source code documented in METHODS\n")
note("METHODS covers the data_export.R / panel_e_compute.R framework",
     if (any(grepl("data_export\\.R|panel_e_compute|figures/multipanel/figure3", methods))) "PASS" else "GAP",
     "Figure 3 panels are rendered via the matplotlib pipeline at figures/multipanel/figure3_isopair_and_ptc/.")

# ---- 9. Cross-reference: Rmd → METHODS.md section pointers ----
cat("\n[9] Rmd cross-references to METHODS.md sections\n")
rmd <- readLines("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_mashr.Rmd")
rmd_methods_refs <- grep("METHODS\\.md", rmd, value = TRUE)
cat(sprintf("  Rmd cites METHODS.md in %d places\n", length(rmd_methods_refs)))
cited_sections <- unique(gsub('.*METHODS\\.md § "([^"]+)".*', "\\1", rmd_methods_refs))
cited_sections <- cited_sections[!grepl("METHODS\\.md", cited_sections)]
# Resolve each cited path: strip optional "Publication Report Methods > " prefix
# and check that every "A > B > C" leaf matches some heading in METHODS.md.
missing_sections <- character(0)
for (sec in cited_sections) {
  # Drop the "Publication Report Methods > " breadcrumb prefix if present.
  sec_norm <- sub("^Publication Report Methods > ", "", sec)
  # The leaf of the breadcrumb is what should appear as a Markdown heading.
  leaf <- sub(".*> ", "", sec_norm)
  if (!any(grepl(leaf, methods, fixed = TRUE))) {
    missing_sections <- c(missing_sections, sec)
  }
}
note("All Rmd-cited METHODS.md sections actually exist in METHODS.md",
     if (length(missing_sections) == 0) "PASS" else "GAP",
     if (length(missing_sections) > 0)
       paste("missing:", paste(missing_sections, collapse = "; "))
     else "")

# ============================================================================
# Summary
# ============================================================================
cat("\n================================================================\n")
n_pass <- sum(sapply(results, function(r) r$status == "PASS"))
n_gap  <- sum(sapply(results, function(r) r$status == "GAP"))
n_fail <- sum(sapply(results, function(r) r$status == "FAIL"))
cat(sprintf("Pass 5 summary: %d PASS / %d GAP / %d FAIL / %d total\n",
            n_pass, n_gap, n_fail, length(results)))

if (n_gap > 0) {
  cat("\nGAPS to address in METHODS.md:\n")
  for (r in results) {
    if (r$status == "GAP") {
      cat(sprintf("  - %s\n    %s\n", r$label, r$detail))
    }
  }
}
cat("================================================================\n")
