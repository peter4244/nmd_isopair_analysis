#!/usr/bin/env Rscript
# Pass 3: Documentation accuracy (NEW 2x2 layout, 2026-06-15).
# Verify that RATIONALE.md §11, the composite legend, and
# paper/section4_findreplace_2026-07-11_4ct.md describe Figure 4 at its CURRENT
# 2x2 / 1050-113-1166 state. Detect stale text from the prior 3x2 / 1489-163-1286
# layout.

suppressPackageStartupMessages({})

ROOT <- "/Users/petecastaldi/claude_projects/nmd/figures/multipanel/figure4_ptcneg_and_model"
PAPER <- "/Users/petecastaldi/claude_projects/nmd/paper"
LEGEND <- file.path(ROOT, "figure4_composite_legend.md")
RATIONALE <- file.path(ROOT, "RATIONALE.md")
FINDREPLACE <- file.path(PAPER, "section4_findreplace_2026-07-11_4ct.md")

read_doc <- function(p) paste(readLines(p, warn = FALSE), collapse = "\n")
legend_txt <- read_doc(LEGEND)
rat_txt <- read_doc(RATIONALE)
fr_txt <- read_doc(FINDREPLACE)

issues <- list()
ok <- list()

flag <- function(label, found, expected_present) {
  if (found == expected_present) {
    ok[[length(ok) + 1L]] <<- label
    cat(sprintf("PASS  %s\n", label))
  } else {
    issues[[length(issues) + 1L]] <<- label
    cat(sprintf("FAIL  %s\n", label))
  }
}

cat("=== LEGEND (figure4_composite_legend.md) ===\n")
flag("Legend mentions n = 48 (Section A PTC+)",  grepl("[^0-9]48 NMD", legend_txt), TRUE)
flag("Legend mentions n = 82 (Section A PTC-)",  grepl("[^0-9]82 NMD", legend_txt), TRUE)
flag("Legend mentions n = 130 Control",        grepl("130 Control", legend_txt), TRUE)
flag("Legend mentions n = 756 (Section C PTC+)", grepl("756", legend_txt), TRUE)
flag("Legend mentions n = 63 (Section C PTC-)",  grepl("[^0-9]63 NMD", legend_txt), TRUE)
flag("Legend mentions n = 819 (Section C Control)", grepl("819", legend_txt), TRUE)
flag("Legend mentions 2x2 / 4-panel structure (A-D)",
     grepl("Panels A.D|four panel|2.2", legend_txt) | grepl("\\bA\\b.*\\bB\\b.*\\bC\\b.*\\bD\\b", legend_txt), TRUE)
flag("Legend does NOT contain 1,489 (old Section C n)",  grepl("1,489|1489", legend_txt), FALSE)
flag("Legend does NOT contain 163 (old PTC- n)",          grepl("[^0-9]163[^0-9]", legend_txt), FALSE)
flag("Legend does NOT contain 1,286 (old Control n)",     grepl("1,286|1286", legend_txt), FALSE)
flag("Legend mentions the TD2 bias supplement",   grepl("TD2 ?Bias|systematic bias of the standard CDS", legend_txt), TRUE)
flag("Legend describes gene+reference matching", grepl("gene.matched|re-intersect|matched on gene and reference", legend_txt), TRUE)

cat("\n=== RATIONALE.md §11 ===\n")
flag("Rationale §11 mentions 2x2 composite",         grepl("2.2 composite|2.2 layout|2.2.*layout", rat_txt), TRUE)
flag("Rationale §11 mentions 12.x8. landscape",      grepl("12.x8.|12.+x.+8|12.*8.*landscape", rat_txt), TRUE)
flag("Rationale §11 mentions consolidated supplement", grepl("TD2BiasEvidence|consolidated", rat_txt), TRUE)
# CRITICAL: Section C numbers should be the new 1050/113/1166, not 1489/163/1286.
flag("Rationale §11 mentions 756 (Section C PTC+ n)", grepl("756", rat_txt), TRUE)
flag("Rationale §11 mentions 819 (Section C Control n)", grepl("819", rat_txt), TRUE)
# Detect stale numbers
flag("Rationale does NOT use 1,489 for Section C PTC+ tables (line 256-258)",
     grepl("\\| \\*\\*NMD\\+/PTC\\+\\*\\* \\| 1,489", rat_txt), FALSE)
flag("Rationale does NOT use 1,286 for Section C Control (line 258)",
     grepl("\\| \\*\\*Control\\*\\* \\| 1,286", rat_txt), FALSE)
flag("Rationale does NOT use 163 for Section C PTC- (line 257)",
     grepl("\\| \\*\\*NMD\\+/PTC-\\*\\* \\| 163", rat_txt), FALSE)
# Also describe the re-intersection step
flag("Rationale describes 1:1 gene-matched re-intersection", grepl("re-intersect", rat_txt), TRUE)
# Comment in data_export.R header should not say 1,489/163/1,286
de_txt <- paste(readLines(file.path(ROOT, "data_export.R")), collapse = "\n")
flag("data_export.R header comment says new 1050/113/1166 (not 1489/163/1286)",
     grepl("1,489/163/1,286|1489/163/1286", de_txt), FALSE)

cat("\n=== section4_findreplace_2026-07-11_4ct.md ===\n")
# Pair 1 (the actual replacement prose at lines ~76-80) must have the new numbers.
flag("FindReplace has 4-CT Section C PTC+ n=756", grepl("756", fr_txt), TRUE)
flag("FindReplace has 4-CT Section C PTC- n=63", grepl("[^0-9]63 NMD|63 had no", fr_txt), TRUE)
flag("FindReplace has 4-CT Section C Control n=819", grepl("819", fr_txt), TRUE)
# Pair 1 should also include the new Section A numbers
flag("FindReplace has 4-CT Section A PTC+ n=48", grepl("48 NMD./PTC.", fr_txt), TRUE)
flag("FindReplace has 4-CT Section A PTC- n=82", grepl("82 NMD./PTC.", fr_txt), TRUE)
# Pair 2 must describe what's now in the supplement (no longer a 2-panel main fig section)
flag("FindReplace acknowledges TD2 evidence", grepl("TD2|SF35|SF34", fr_txt), TRUE)
# Stale summary text
flag("FindReplace does NOT call Figure 4 a '6-panel'",  grepl("6-panel|six-panel|six panel", fr_txt), FALSE)
flag("FindReplace does NOT reference 3x2 / 12.x12. portrait",
     grepl("3.2 layout|12.+x.+12.+portrait|3.2.*portrait", fr_txt), FALSE)
flag("FindReplace summary table uses NEW Section C numbers (1,050/113/1,166)",
     grepl("756|819", fr_txt), TRUE)
# The summary "What §4 is now saying" block (line 17) says 1,489/163/1,286 — flag
flag("FindReplace summary does NOT still claim Section C n=1,489/163/1,286",
     grepl("1,489/163/1,286|n=1,489", fr_txt), FALSE)
# Section B (TD2) is referenced as separate section in the find/replace doc summary
flag("FindReplace does NOT describe Section B as main fig (now supplement)",
     grepl("\\*\\*Section B \\(n=900 \\+ n=1,032", fr_txt), FALSE)
# The Section C descriptive table in find/replace says median 460 nt for PTC- 5'UTR — old, should be 488
flag("FindReplace Section C PTC- 5'UTR median is 488 (not 460)",
     grepl("488", fr_txt) | !grepl("460 nt \\(231\\.5", fr_txt), TRUE)
flag("FindReplace Section C PTC+ 5'UTR IQR uses 134/421 (new) not 137/422 (old)",
     grepl("134.*421|421\\.75", fr_txt) | !grepl("137.422|137.422\\)", fr_txt), TRUE)
# Files list at line 137-143 references panelC_td2 / panelD_kozak / panelE / panelF — stale
flag("FindReplace file list does NOT reference panelC_td2_vs_refaug_length (stale)",
     grepl("panelC_td2_vs_refaug_length", fr_txt), FALSE)
flag("FindReplace file list does NOT reference panelD_kozak (stale)",
     grepl("panelD_kozak", fr_txt), FALSE)
flag("FindReplace file list does NOT reference panelE_5utr_length / panelF (stale)",
     grepl("panelE_5utr|panelF_longest", fr_txt), FALSE)

cat("\n=== SUMMARY ===\n")
cat(sprintf("PASS: %d   FAIL: %d\n", length(ok), length(issues)))
if (length(issues) > 0) {
  cat("Failures:\n")
  for (m in issues) cat(" -", m, "\n")
}
