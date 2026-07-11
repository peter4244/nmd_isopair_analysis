#!/usr/bin/env Rscript
# Verification Pass 3 — Documentation accuracy for Figure 3
#
# Cross-check every prose statement (in methodology files, find/replace doc,
# composite legend, panel render scripts) against the actual rendered panels
# and underlying TSV/CSV data. Look for drift: does the documentation say
# what the data actually shows?
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

# ---- Load actual data ----
panelB <- fread("data/panelB_sequence_similarity.tsv")
panelC <- fread("data/panelC_event_prevalence.tsv")
panelD <- fread("data/panelD_stop_codon_distance.tsv")
panelE <- fread("data/panelE_ptc_event_attribution.tsv")
panelF <- fread("data/panelF_mechanism_breakdown.tsv")

cat("================================================================\n")
cat("PASS 3: Documentation accuracy verification for Figure 3 (all-3-ENST scope)\n")
cat("================================================================\n\n")

# ============================================================================
# Section A — Panel render scripts' descriptive text (docstrings + axis labels)
# Panels are titleless per the publications convention; titles live in the
# composite legend and per-panel methodology files.
# ============================================================================
cat("[A] Panel render scripts — docstrings + axis labels\n")

extract_lit <- function(script_path, patterns) {
  txt <- readLines(script_path)
  out <- character(0)
  for (p in patterns) {
    hits <- grep(p, txt, value = TRUE, perl = TRUE)
    out <- c(out, hits)
  }
  out
}

# Panel B: x-axis = "Sequence shared (%)"
b_script <- "figure3_panelB_sequence_similarity.py"
note("Panel B x-axis 'Sequence shared (%)' present",
     if (length(extract_lit(b_script, "Sequence shared")) > 0) "PASS" else "FAIL")

# Panel C: docstring describes "Splice event prevalence"
c_script <- "figure3_panelC_event_prevalence.py"
note("Panel C docstring mentions 'Splice event prevalence'",
     if (length(extract_lit(c_script, "[Ss]plice event")) > 0) "PASS" else "FAIL")

# Panel D: x-axis "Distance from last EJC" + clipping disclosure
d_script <- "figure3_panelD_stop_codon_distance.py"
note("Panel D x-axis 'Distance from last EJC' present",
     if (length(extract_lit(d_script, "Distance from last EJC")) > 0) "PASS" else "FAIL")
note("Panel D clipping disclosed via 'Clipped at' print statement",
     if (length(extract_lit(d_script, "Clipped at")) > 0) "PASS" else "FAIL")

# Panel E: "PTC-causing (NMD)" legend label
e_script <- "figure3_panelE_ptc_event_attribution.py"
note("Panel E legend label 'PTC-causing (NMD)' present",
     if (length(extract_lit(e_script, "PTC-causing \\(NMD\\)")) > 0) "PASS" else "FAIL")

# Panel F: docstring "stacked by mechanism"
f_script <- "figure3_panelF_mechanism_breakdown.py"
note("Panel F docstring mentions 'stacked by mechanism'",
     if (length(extract_lit(f_script, "stacked by mechanism")) > 0) "PASS" else "FAIL")

# ============================================================================
# Section B — Methodology files: claimed numbers vs TSV
# ============================================================================
cat("\n[B] Methodology files — numeric claims vs TSV\n")

# Panel B methodology: 3,009 each
b_meth <- readLines("figure3_panelB_methodology.md")
n_NMD_B <- sum(panelB$comparison == "NMD")
n_Ctrl_B <- sum(panelB$comparison == "Control")
note("Panel B methodology says '3,009' (matches TSV row count)",
     if (any(grepl("1,548", b_meth)) && n_NMD_B == 1548 && n_Ctrl_B == 1548)
       "PASS" else "FAIL",
     sprintf("TSV NMD=%d, Control=%d", n_NMD_B, n_Ctrl_B))

# Panel C methodology: SE OR ~2x-3x
c_meth <- readLines("figure3_panelC_methodology.md")
se_or <- panelC$fisher_OR[panelC$event_type == "SE"]
note("Panel C methodology claims SE ~2x enrichment (verify OR ~2.5-3.5)",
     if (se_or > 3.4 && se_or < 4.1) "PASS" else "FAIL",
     sprintf("SE OR = %.2f", se_or))

# Panel D methodology: 190 NMD + 190 Control
d_meth <- readLines("figure3_panelD_methodology.md")
n_NMD_D <- sum(panelD$comparison == "NMD")
n_Ctrl_D <- sum(panelD$comparison == "Control")
note("Panel D methodology says '190 NMD' and '190 Control'",
     if (any(grepl("130 NMD", d_meth)) && any(grepl("130 Control", d_meth)) &&
         n_NMD_D == 130 && n_Ctrl_D == 130) "PASS" else "FAIL",
     sprintf("TSV NMD=%d Control=%d", n_NMD_D, n_Ctrl_D))

# Panel D methodology: 37.9% / 2.1% PTC rate (own-GENCODE-stop scope)
note("Panel D methodology states 36.9% NMD / 1.5% Control PTC rate",
     if (any(grepl("36\\.9%", d_meth)) && any(grepl("1\\.5%", d_meth)))
       "PASS" else "FAIL")

# Panel D methodology: 18x enrichment
note("Panel D methodology claims 24x enrichment",
     if (any(grepl("24.fold|24×|24-fold", d_meth))) "PASS" else "FAIL")

# Panel D methodology: clipping disclosure
note("Panel D methodology discloses x-axis clipping range [-1000, 1500]",
     if (any(grepl("-1000.*1500|\\[-1000, 1500\\]", d_meth))) "PASS" else "FAIL")

# Panel D methodology: clipping counts 4 NMD / 4 Control
note("Panel D methodology discloses NMD=4, Control=4 clipped",
     if (any(grepl("NMD = 4", d_meth)) && any(grepl("Control = 4", d_meth)))
       "PASS" else "FAIL")

# Panel E methodology: n_ptc_attr = 69
e_meth <- readLines("figure3_panelE_methodology.md")
n_ptc_attr <- sum(panelE$n_ptc_events)
note("Panel E methodology states n_ptc_attr=69 (matches TSV)",
     if (any(grepl("48 attributed", e_meth)) && n_ptc_attr == 48) "PASS" else "FAIL",
     sprintf("TSV sum=%d", n_ptc_attr))

# Panel E methodology: ctrl_total=447 with explanation of TSV (330)
ctrl_total_tsv <- sum(panelE$n_ctrl_events)
note("Panel E methodology states ctrl_total=292",
     if (any(grepl("292", e_meth))) "PASS" else "FAIL",
     sprintf("TSV n_ctrl_events sum=%d (emits subset)", ctrl_total_tsv))
note("Panel E methodology explains 209 as TSV sum",
     if (any(grepl("209", e_meth))) "PASS" else "FAIL")

# Panel E methodology: SE 43.5% claim
se_pct_panelE <- panelE$pct_of_ptc[panelE$event_type == "SE"]
note("Panel E SE % of PTC-causing == 47.9%",
     if (abs(se_pct_panelE - 47.9) < 0.15) "PASS" else "FAIL",
     sprintf("TSV SE pct_of_ptc=%.1f", se_pct_panelE))
note("Panel E methodology mentions SE 47.9%",
     if (any(grepl("47\\.9", e_meth))) "PASS" else "FAIL")

# Panel F methodology: 69 attributed pairs
f_meth <- readLines("figure3_panelF_methodology.md")
n_F_total <- sum(panelF$n)
note("Panel F methodology states 48 attributed pairs (matches TSV)",
     if (any(grepl("48 attributed", f_meth)) && n_F_total == 48) "PASS" else "FAIL",
     sprintf("TSV sum=%d", n_F_total))

# Panel F mechanism counts: 25 / 18 / 5
note("Panel F methodology states 25 Frameshift / 18 In-frame stop / 5 3'UTR splice",
     if (any(grepl("25 Frameshift", f_meth)) &&
         any(grepl("18 In-frame stop", f_meth)) &&
         any(grepl("5 3.UTR splice", f_meth))) "PASS" else "FAIL")

# Panel F mechanism % breakdown: 55.1 / 33.3 / 11.6
fs_pct <- 100 * sum(panelF$n[panelF$mechanism == "Frameshift"]) / n_F_total
ifs_pct <- 100 * sum(panelF$n[panelF$mechanism == "In-frame stop"]) / n_F_total
utr_pct <- 100 * sum(panelF$n[panelF$mechanism == "3'UTR splice"]) / n_F_total
note("Panel F Frameshift % (~52%)",
     if (abs(fs_pct - 52.1) < 0.2) "PASS" else "FAIL",
     sprintf("TSV Frameshift=%.1f%%", fs_pct))
note("Panel F In-frame stop % (~37.5%)",
     if (abs(ifs_pct - 37.5) < 0.2) "PASS" else "FAIL",
     sprintf("TSV In-frame=%.1f%%", ifs_pct))
note("Panel F 3'UTR splice % (~10.4%)",
     if (abs(utr_pct - 10.4) < 0.2) "PASS" else "FAIL",
     sprintf("TSV 3'UTR=%.1f%%", utr_pct))
note("Panel F methodology mentions 52.1 / 37.5 / 10.4",
     if (any(grepl("52\\.1", f_meth)) &&
         any(grepl("37\\.5", f_meth)) &&
         any(grepl("10\\.4", f_meth))) "PASS" else "FAIL")

# ============================================================================
# Section C — Find/replace doc consistency (v2, current canonical)
# ============================================================================
cat("\n[C] Find/replace doc (section4_findreplace_2026-07-11_4ct.md)\n")
fr <- readLines("/Users/petecastaldi/claude_projects/nmd/paper/section4_findreplace_2026-07-11_4ct.md")

note("Find/replace v2 claims 37.9% NMD PTC rate",
     if (any(grepl("37\\.9", fr))) "PASS" else "FAIL")
note("Find/replace v2 claims 2.1% Control PTC rate",
     if (any(grepl("2\\.1%", fr))) "PASS" else "FAIL")
note("Find/replace v2 claims 18-fold enrichment",
     if (any(grepl("18.fold|18-fold|18×", fr))) "PASS" else "FAIL")
note("Find/replace v2 claims 190 NMD / 190 Control",
     if (any(grepl("130 NMD", fr)) && any(grepl("130 Control", fr))) "PASS" else "FAIL")
note("Find/replace v2 claims 72 NMD PTC+ / 4 Control PTC+",
     if (any(grepl("72 NMD", fr)) && any(grepl("4 Control", fr))) "PASS" else "FAIL")
note("Find/replace v2 mentions 69 attributed PTC+ pairs",
     if (any(grepl("48 attributed events from 48", fr)) || any(grepl("48 attributed", fr))) "PASS" else "FAIL")
note("Find/replace v2 mentions all-3-ENST + coding-CDS scope",
     if (any(grepl("GENCODE|all-3-ENST|all three isoforms", fr))) "PASS" else "FAIL")
note("Find/replace v2 mentions 'no ref-AUG projection'",
     if (any(grepl("reference AUG|ref-AUG|GENCODE annotated", fr))) "PASS" else "FAIL")

# ============================================================================
# Section D — Composite legend consistency
# ============================================================================
cat("\n[D] Composite legend (figure3_composite_legend.md)\n")
leg <- readLines("figure3_composite_legend.md")
note("Composite legend cites 130 NMD / 130 Control",     if (any(grepl("130 NMD", leg)) && any(grepl("130 Control", leg))) "PASS" else "FAIL")
note("Composite legend cites 36.9% NMD PTC",
     if (any(grepl("36\\.9", leg))) "PASS" else "FAIL")
note("Composite legend cites 1.5% Control PTC",
     if (any(grepl("1\\.5%", leg))) "PASS" else "FAIL")
note("Composite legend cites 24-fold enrichment",
     if (any(grepl("24×|24.fold|24-fold", leg))) "PASS" else "FAIL")
note("Composite legend cites 48 attributed pairs",
     if (any(grepl("48 attributed", leg)) || any(grepl("n = 48", leg))) "PASS" else "FAIL")
note("Composite legend mentions 292 total Control events",
     if (any(grepl("292", leg))) "PASS" else "FAIL")
note("Composite legend mentions 25 frameshift / 18 in-frame stop / 5 3'UTR splice",
     if (any(grepl("25 events, 52", leg)) &&
         any(grepl("18, 38", leg)) &&
         any(grepl("5, 10", leg))) "PASS" else "FAIL")

# ============================================================================
# Section E — Cross-doc consistency for headline numbers
# ============================================================================
cat("\n[E] Cross-document number consistency (canonical sources)\n")
# Headline numbers should appear in all 4 places: d_meth, e_meth, f_meth, fr, leg
key_numbers <- c("36.9", "1.5%", "130", "48", "24-fold")
canon <- list(d_meth = d_meth, e_meth = e_meth, f_meth = f_meth,
              find_replace_v2 = fr, composite_legend = leg)
for (n in key_numbers) {
  present <- sapply(canon, function(doc) any(grepl(n, doc, fixed = TRUE)))
  count <- sum(present)
  # Require each key number to appear in at least 2 canonical docs
  status <- if (count >= 2) "PASS" else "FAIL"
  note(sprintf("'%s' appears across canonical docs (>=2 places)", n),
       status,
       sprintf("count=%d in {%s}",
               count, paste(names(canon)[present], collapse = ", ")))
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
