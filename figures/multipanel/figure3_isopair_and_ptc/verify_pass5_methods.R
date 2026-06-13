#!/usr/bin/env Rscript
# Verification Pass 5 — METHODS.md verification for Figure 3
#
# Cross-check the canonical METHODS.md against the actual code behavior. Flag
# documentation drift, missing sections, and gaps that Task #23 should close.

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
cat("PASS 5: METHODS.md verification for Figure 3\n")
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
cat("\n[2] Reference-ATG Tracing categories — METHODS vs actual code output\n")
ref_atg <- readRDS("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr/analysis_cache/ref_atg_analysis.rds")
actual_cats <- sort(unique(c(ref_atg$c2$category, ref_atg$c4$category)))
cat(sprintf("  Actual categories in ref_atg_analysis.rds: %s\n",
            paste(actual_cats, collapse = ", ")))

# METHODS lists: effectively_ptc, truncated_no_ejc, ref_atg_lost, no_downstream_ejc, same_or_longer_with_ejc
methods_cats_listed <- c("effectively_ptc", "truncated_no_ejc", "ref_atg_lost",
                        "no_downstream_ejc", "same_or_longer_with_ejc")
documented_but_not_emitted <- setdiff(methods_cats_listed, actual_cats)
emitted_but_not_documented <- setdiff(actual_cats, methods_cats_listed)
note("All emitted categories documented in METHODS",
     if (length(emitted_but_not_documented) == 0) "PASS" else "GAP",
     if (length(emitted_but_not_documented) > 0)
       paste("missing from METHODS:", paste(emitted_but_not_documented, collapse = ", "))
     else "")
note("All METHODS-listed categories actually emitted",
     if (length(documented_but_not_emitted) == 0) "PASS" else "GAP",
     if (length(documented_but_not_emitted) > 0)
       paste("documented but not emitted:", paste(documented_but_not_emitted, collapse = ", "))
     else "")

# ---- 3. New ref-AUG-traceable scope — is it documented? ----
cat("\n[3] New ref-AUG-traceable scope (Pete 2026-06-13)\n")
note("METHODS documents the new ref-AUG-traceable scope (Pete's 2026-06-13 policy)",
     if (any(grepl("ref-AUG-traceable|traceable subset|pop_traceable", methods))) "PASS" else "GAP",
     "Pete's policy: PTC analyses use the ref-AUG-traceable subset (categories effectively_ptc + no_downstream_ejc + truncated_no_ejc) — this is now the canonical population, NOT the old 'TD2 + ref-AUG occult-PTC rescue' framing.")

note("METHODS documents the 83.5% NMD / 16.3% Control PTC rate headline",
     if (any(grepl("83\\.5|83\\.5%|1,912 of 2,289", methods))) "PASS" else "GAP",
     "These are the new manuscript-headline numbers.")

# ---- 4. Coding-coding pre-filter for 05r_ref_atg_analysis.R ----
cat("\n[4] Coding-coding pre-filter for ref-AUG analysis\n")
note("METHODS notes that ref-AUG runs on coding-coding pairs only",
     if (any(grepl("both isoforms are coding|coding-coding|coding pairs", methods))) "PASS" else "GAP",
     "Important: explains why 336 NMD pairs are excluded from ref-AUG tracing.")

# ---- 5. Mixed-source attribution chain (Panel E/F) ----
cat("\n[5] Mixed-source attribution chain (new ref-AUG-traceable Panel E/F)\n")
note("METHODS documents the mixed-source PTC attribution",
     if (any(grepl("attribute_ptc_events.*ref_atg|attr_mechanism.*ref-AUG", methods))) "PASS" else "GAP",
     "Under the new scope, attribution comes from 2 sources: (a) ref_atg_analysis$attr_* for ref-AUG-recovered (original_ptc=FALSE); (b) attribute_ptc_events on TD2-PTC+ subset (original_ptc=TRUE). This methodological mixing should be documented.")

# ---- 6. CDS-source-by-isoform-class clarity ----
cat("\n[6] GENCODE/TD2 split clarity for the manuscript audience\n")
note("METHODS explicitly states what % of NMD comparators are novel vs annotated",
     if (any(grepl("80%.*novel|76.8%", methods))) "PASS" else "GAP",
     "Already partly documented (76.8% noted in occult-PTC section). Could be stated more prominently.")

# ---- 7. Compute_ptc_rates_row description ----
cat("\n[7] compute_ptc_rates_row() — denominator chain\n")
note("METHODS documents the Stage 2 + coding-coding + gene_id re-match chain",
     if (any(grepl("Stage 2.*gene-match.*coding|Fisher.*denominator", methods))) "PASS" else "GAP",
     "compute_ptc_rates_row applies a 3-stage filter (Stage 2 gene-match + coding-coding + gene_id re-match) that wasn't documented separately and was a source of confusion. Should be explicitly described.")

# ---- 8. Stale numbers / claims in METHODS ----
cat("\n[8] Stale numbers cross-check\n")
note("METHODS doesn't cite the old 44% PTC rate (would be stale)",
     if (!any(grepl("44%.*PTC|n_ptc.*44|43\\.9%", methods))) "PASS" else "GAP",
     "If METHODS cites 44%, it needs updating.")
note("METHODS doesn't cite the old 3,674 gene count",
     if (!any(grepl("3,?674.*genes", methods))) "PASS" else "GAP")

# ---- 9. New analyses introduced in this session ----
cat("\n[9] New analyses from this session that need METHODS coverage\n")
note("METHODS covers the data_export.R / panel_e_compute.R framework",
     if (any(grepl("data_export\\.R|panel_e_compute|figures/multipanel", methods))) "PASS" else "GAP",
     "Figure 3 panels are now rendered via a separate matplotlib pipeline in figures/multipanel/figure3_isopair_and_ptc/; this is downstream of the canonical Isopair analysis but is its own reproducibility artifact.")

# ---- 10. Cross-reference: Rmd → METHODS.md section pointers ----
cat("\n[10] Rmd cross-references to METHODS.md sections\n")
rmd <- readLines("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_mashr.Rmd")
rmd_methods_refs <- grep("METHODS\\.md", rmd, value = TRUE)
cat(sprintf("  Rmd cites METHODS.md in %d places\n", length(rmd_methods_refs)))
# Verify each cited section exists
cited_sections <- unique(gsub('.*METHODS\\.md § "([^"]+)".*', "\\1", rmd_methods_refs))
cited_sections <- cited_sections[!grepl("METHODS\\.md", cited_sections)]
missing_sections <- character(0)
for (sec in cited_sections) {
  if (!any(grepl(sec, methods, fixed = TRUE))) {
    missing_sections <- c(missing_sections, sec)
  }
}
note("All Rmd-cited METHODS.md sections actually exist in METHODS.md",
     if (length(missing_sections) == 0) "PASS" else "GAP",
     if (length(missing_sections) > 0)
       paste("missing:", paste(missing_sections, collapse = "; "))
     else "")

# ============================================================================
# Summary + Task #23 gap list
# ============================================================================
cat("\n================================================================\n")
n_pass <- sum(sapply(results, function(r) r$status == "PASS"))
n_gap  <- sum(sapply(results, function(r) r$status == "GAP"))
n_fail <- sum(sapply(results, function(r) r$status == "FAIL"))
cat(sprintf("Pass 5 summary: %d PASS / %d GAP (for Task #23) / %d FAIL / %d total\n",
            n_pass, n_gap, n_fail, length(results)))

if (n_gap > 0) {
  cat("\nGAPS to address in Task #23 (upstream Rmd + METHODS.md updates):\n")
  for (r in results) {
    if (r$status == "GAP") {
      cat(sprintf("  - %s\n    %s\n", r$label, r$detail))
    }
  }
}
cat("================================================================\n")
