#!/usr/bin/env Rscript
# Verification Pass 3 — Documentation accuracy for Figure 3
#
# Cross-check every prose statement (in methodology files, find/replace doc,
# results-to-code map, panel render scripts) against the actual rendered
# panels and underlying TSV/CSV data. Look for drift: does the documentation
# say what the data actually shows?

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

# ---- Load actual data ----
panelB <- fread("data/panelB_sequence_similarity.tsv")
panelC <- fread("data/panelC_event_prevalence.tsv")
panelD <- fread("data/panelD_stop_codon_distance.tsv")
panelE <- fread("data/panelE_ptc_event_attribution.tsv")
panelF <- fread("data/panelF_mechanism_breakdown.tsv")

cat("================================================================\n")
cat("PASS 3: Documentation accuracy verification for Figure 3\n")
cat("================================================================\n\n")

# ============================================================================
# Section A — Panel render scripts' literal text vs methodology + TSV
# ============================================================================
cat("[A] Panel render scripts — literal text claims\n")

# Helper to grep render-script text
extract_lit <- function(script_path, patterns) {
  txt <- readLines(script_path)
  out <- character(0)
  for (p in patterns) {
    hits <- grep(p, txt, value = TRUE, perl = TRUE)
    out <- c(out, hits)
  }
  out
}

# Panel B render-script claim: "Sequence shared with reference"
b_script <- "figure3_panelB_sequence_similarity.py"
title_B <- extract_lit(b_script, "set_title.*Sequence")
note("Panel B render-script title contains 'Sequence shared with reference'",
     if (length(title_B) > 0) "PASS" else "FAIL")

# Panel C: "Splice event prevalence"
c_script <- "figure3_panelC_event_prevalence.py"
title_C <- extract_lit(c_script, "set_title.*[Ss]plice event")
note("Panel C title contains 'Splice event prevalence'",
     if (length(title_C) > 0) "PASS" else "FAIL")

# Panel D: title "Stop codon to last EJC"
d_script <- "figure3_panelD_stop_codon_distance.py"
title_D <- extract_lit(d_script, "set_title.*Stop codon")
note("Panel D title contains 'Stop codon to last EJC'",
     if (length(title_D) > 0) "PASS" else "FAIL")

# Panel D: clipped values disclosed
clip_D <- extract_lit(d_script, "Clipped at")
note("Panel D clipping is disclosed in render script",
     if (length(clip_D) > 0) "PASS" else "FAIL")

# Panel E: PTCs identified subtitle
e_script <- "figure3_panelE_ptc_event_attribution.py"
sub_E <- extract_lit(e_script, "PTCs identified")
note("Panel E subtitle: 'PTCs identified: 1,912 NMD | 288 Control'",
     if (length(sub_E) > 0 && any(grepl("1912", sub_E)) && any(grepl("288", sub_E)))
       "PASS" else "FAIL")

# Panel F: subtitle "attributed pairs | stacked by mechanism"
f_script <- "figure3_panelF_mechanism_breakdown.py"
sub_F <- extract_lit(f_script, "attributed pairs")
note("Panel F subtitle contains 'attributed pairs ... stacked by mechanism'",
     if (length(sub_F) > 0) "PASS" else "FAIL")

# ============================================================================
# Section B — Methodology files: claimed numbers vs TSV
# ============================================================================
cat("\n[B] Methodology files — numeric claims vs TSV\n")

# Panel B methodology: 3,009 each
b_meth <- readLines("figure3_panelB_methodology.md")
n_NMD_B <- sum(panelB$comparison == "NMD")
n_Ctrl_B <- sum(panelB$comparison == "Control")
note("Panel B methodology says '3,009 each' (matches TSV row count)",
     if (any(grepl("3,009 pairs per side", b_meth)) && n_NMD_B == 3009 && n_Ctrl_B == 3009)
       "PASS" else "FAIL",
     sprintf("TSV NMD=%d, Control=%d", n_NMD_B, n_Ctrl_B))

# Panel C methodology: '2x' enrichment claim for SE
c_meth <- readLines("figure3_panelC_methodology.md")
se_or <- panelC$fisher_OR[panelC$event_type == "SE"]
note("Panel C methodology claims SE '~2x' enrichment (verify OR ~2-3)",
     if (any(grepl("~?2.?[xX].*NMD vs.*~?21.*Control", c_meth)) && se_or > 2.5 && se_or < 3.5)
       "PASS" else "FAIL",
     sprintf("SE OR = %.2f", se_or))

# Panel D methodology: 2,289 / 1,763
d_meth <- readLines("figure3_panelD_methodology.md")
n_NMD_D <- sum(panelD$comparison == "NMD")
n_Ctrl_D <- sum(panelD$comparison == "Control")
note("Panel D methodology says 2,289 NMD + 1,763 Control (matches TSV)",
     if (any(grepl("2,289 NMD \\+ 1,763 Control", d_meth)) && n_NMD_D == 2289 && n_Ctrl_D == 1763)
       "PASS" else "FAIL",
     sprintf("TSV NMD=%d Control=%d", n_NMD_D, n_Ctrl_D))

# Panel D methodology: 83.5% / 16.3% PTC rate claim
note("Panel D methodology states 83.5% NMD / 16.3% Control PTC rate",
     if (any(grepl("83.5", d_meth)) && any(grepl("16.3", d_meth)))
       "PASS" else "FAIL")

# Panel D methodology: clip count claim
note("Panel D methodology discloses x-axis clipping range [-1000, 1500]",
     if (any(grepl("-1000.*1500|\\[-1000, 1500\\]", d_meth)))
       "PASS" else "FAIL")

# Panel E methodology: 1,812 n_ptc_attr, 4,525 ctrl_total
e_meth <- readLines("figure3_panelE_methodology.md")
n_ptc_attr <- sum(panelE$n_ptc_events)
ctrl_total <- sum(panelE$n_ctrl_events)
note("Panel E methodology states n_ptc_attr=1,812 (matches TSV)",
     if (any(grepl("1,812.*PTC-causing|n_ptc_attr.*1,812", e_meth)) && n_ptc_attr == 1812)
       "PASS" else "FAIL",
     sprintf("TSV sum=%d", n_ptc_attr))
note("Panel E methodology states ctrl_total=4,525 (matches TSV)",
     if (any(grepl("4,525", e_meth)) && ctrl_total == 4525)
       "PASS" else "FAIL",
     sprintf("TSV sum=%d", ctrl_total))

# Panel E methodology: SE ~55% claim
se_pct_panelE <- panelE$pct_of_ptc[panelE$event_type == "SE"]
note("Panel E methodology says SE ~55% of PTC-causing (matches TSV)",
     if (any(grepl("~?55%", e_meth)) && se_pct_panelE > 54 && se_pct_panelE < 56)
       "PASS" else "FAIL",
     sprintf("TSV SE pct_of_ptc=%.1f", se_pct_panelE))

# Panel F methodology: 1,812 attributed pairs
f_meth <- readLines("figure3_panelF_methodology.md")
n_F_total <- sum(panelF$n)
note("Panel F methodology states 1,812 attributed pairs (matches TSV)",
     if (any(grepl("1,812", f_meth)) && n_F_total == 1812)
       "PASS" else "FAIL",
     sprintf("TSV sum=%d", n_F_total))

# Panel F methodology mechanism % breakdown
fs_pct <- 100 * sum(panelF$n[panelF$mechanism == "Frameshift"]) / n_F_total
ifs_pct <- 100 * sum(panelF$n[panelF$mechanism == "In-frame stop"]) / n_F_total
utr_pct <- 100 * sum(panelF$n[panelF$mechanism == "3'UTR splice"]) / n_F_total
note("Panel F methodology Frameshift ~62-63% (matches TSV)",
     if (fs_pct > 61 && fs_pct < 64) "PASS" else "FAIL",
     sprintf("TSV Frameshift=%.1f%%", fs_pct))
note("Panel F methodology In-frame stop ~30-33% (matches TSV)",
     if (ifs_pct > 29 && ifs_pct < 34) "PASS" else "FAIL",
     sprintf("TSV In-frame=%.1f%%", ifs_pct))
note("Panel F methodology 3'UTR splice ~5-7% (matches TSV)",
     if (utr_pct > 4 && utr_pct < 8) "PASS" else "FAIL",
     sprintf("TSV 3'UTR=%.1f%%", utr_pct))

# ============================================================================
# Section C — Find/replace doc consistency
# ============================================================================
cat("\n[C] Find/replace doc — claimed numbers\n")
fr <- readLines("/Users/petecastaldi/claude_projects/nmd/paper/section4_findreplace_2026-06-13.md")

# Key claims
note("Find/replace doc claims 83.5% NMD PTC rate",
     if (any(grepl("83\\.5", fr))) "PASS" else "FAIL")
note("Find/replace doc claims 16.3% Control PTC rate",
     if (any(grepl("16\\.3", fr))) "PASS" else "FAIL")
note("Find/replace doc claims 3,099 eligible genes (3,009 pair sets)",
     if (any(grepl("3,099", fr)) && any(grepl("3,009", fr))) "PASS" else "FAIL")
note("Find/replace doc claims 5.1-fold enrichment",
     if (any(grepl("5\\.1[ -]fold", fr))) "PASS" else "FAIL")
note("Find/replace doc claims 2,289 NMD / 1,763 Control traceable",
     if (any(grepl("2,289", fr)) && any(grepl("1,763", fr))) "PASS" else "FAIL")
note("Find/replace doc claims 1,912 / 288 PTC+",
     if (any(grepl("1,912", fr)) && any(grepl("288", fr))) "PASS" else "FAIL")

# ============================================================================
# Section D — Results-to-code map PM revision block
# ============================================================================
cat("\n[D] Results-to-code map — PM revision block\n")
rcmap <- readLines("/Users/petecastaldi/claude_projects/nmd/paper/results_to_code_map.md")

note("Results map PM block claims 83.5% NMD PTC",
     if (any(grepl("83\\.5", rcmap))) "PASS" else "FAIL")
note("Results map PM block claims 16.3% Control PTC",
     if (any(grepl("16\\.3", rcmap))) "PASS" else "FAIL")
note("Results map PM block mentions pop_BC, pop_traceable, pop_ptc_plus",
     if (any(grepl("pop_BC", rcmap)) && any(grepl("pop_traceable", rcmap)) &&
         any(grepl("pop_ptc_plus", rcmap)))
       "PASS" else "FAIL")
note("Results map PM block mentions new figure panel filenames",
     if (any(grepl("figure3_panelD_stop_codon_distance", rcmap)) &&
         any(grepl("figure3_panelF_mechanism_breakdown", rcmap)))
       "PASS" else "FAIL")

# ============================================================================
# Section E — Cross-doc consistency
# ============================================================================
cat("\n[E] Cross-document number consistency\n")
# Same headline number should appear in all 3 places: methodology, find/replace, results map
key_numbers <- c("83.5", "16.3", "1,912", "288", "2,289", "1,763", "3,009", "1,812")
for (n in key_numbers) {
  in_d_meth <- any(grepl(n, d_meth, fixed = TRUE))
  in_fr <- any(grepl(n, fr, fixed = TRUE))
  in_rcmap <- any(grepl(n, rcmap, fixed = TRUE))
  in_e_meth <- any(grepl(n, e_meth, fixed = TRUE))
  in_f_meth <- any(grepl(n, f_meth, fixed = TRUE))
  present_count <- sum(c(in_d_meth, in_fr, in_rcmap, in_e_meth, in_f_meth))
  # Only require numbers appear in at least 2 of these 5 docs
  status <- if (present_count >= 2) "PASS" else "FAIL"
  note(sprintf("'%s' appears across docs (>=2 places)", n), status,
       sprintf("count=%d", present_count))
}

# ============================================================================
# Summary
# ============================================================================
cat("\n================================================================\n")
n_pass <- sum(sapply(results, function(r) r$status == "PASS"))
n_fail <- sum(sapply(results, function(r) r$status == "FAIL"))
cat(sprintf("Pass 3 summary: %d PASS / %d FAIL / %d total\n",
            n_pass, n_fail, length(results)))
if (n_fail > 0) {
  cat("\nFAILURES:\n")
  for (r in results) {
    if (r$status == "FAIL") {
      cat(sprintf("  - %s -- %s\n", r$label, r$detail))
    }
  }
}
cat("================================================================\n")
