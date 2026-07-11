#!/usr/bin/env Rscript
# Pass 5: METHODS.md verification (NEW 2x2 layout, 2026-06-15).
# Compare METHODS.md Figure 4 section against (a) actual code in data_export.R,
# (b) the rendered figure structure, and (c) the re-intersection step in
# Section C.

ROOT <- "/Users/petecastaldi/claude_projects/nmd/figures/multipanel/figure4_ptcneg_and_model"
METHODS <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/METHODS.md"

methods_txt <- paste(readLines(METHODS, warn = FALSE), collapse = "\n")
methods_lines <- readLines(METHODS, warn = FALSE)
# Slice the Figure 4 section
fig4_start <- grep("^### Figure 4 / §4 — NMD\\+/PTC− Mechanism Inference", methods_lines)
fig4_end <- grep("^### ", methods_lines)
fig4_end <- fig4_end[fig4_end > fig4_start][1]
fig4_block <- paste(methods_lines[fig4_start:(fig4_end - 1)], collapse = "\n")
cat("=== METHODS.md Figure 4 section (lines ", fig4_start, "-", fig4_end - 1, ") ===\n", sep = "")
cat(fig4_block, "\n\n")

pass_ct <- 0L; fail_ct <- 0L; fails <- list()
pass <- function(m) { pass_ct <<- pass_ct + 1L; cat("PASS  ", m, "\n", sep = "") }
fail <- function(m) { fail_ct <<- fail_ct + 1L; fails[[length(fails) + 1L]] <<- m; cat("FAIL  ", m, "\n", sep = "") }

cat("\n=== Structural claims ===\n")
if (grepl("2.2 layout|2.2 composite|2.2.*layout|2..2.*landscape", fig4_block)) {
  pass("METHODS Fig 4 mentions 2x2 layout")
} else {
  fail("METHODS Fig 4 does NOT mention 2x2 layout")
}
if (grepl("12.x8.|12. .8.|12. by 8|12.+x.+8.+landscape", fig4_block)) {
  pass("METHODS Fig 4 mentions 12x8 landscape")
} else {
  fail("METHODS Fig 4 does NOT mention 12x8 landscape")
}
if (grepl("Panels A.D|Panels A, B, C, D|four panel|4 panel", fig4_block)) {
  pass("METHODS Fig 4 mentions Panels A-D / four panel structure")
} else {
  fail("METHODS Fig 4 does NOT clearly describe 4-panel structure")
}

cat("\n=== Re-intersection claim ===\n")
if (grepl("re-intersect|gene-matched re-intersection|1:1 gene-matched", fig4_block)) {
  pass("METHODS Fig 4 describes the 1:1 gene-matched re-intersection")
} else {
  fail("METHODS Fig 4 does NOT describe the 1:1 gene-matched re-intersection step on (gene_id, reference_isoform_id) after the TR filter")
}

cat("\n=== Section C numbers ===\n")
if (grepl("756", fig4_block)) {
  pass("METHODS Fig 4 mentions 756 (new Section C PTC+)")
} else {
  fail("METHODS Fig 4 does NOT mention 756 (new Section C PTC+ n)")
}
if (grepl("[^,0-9]63[^,0-9]", fig4_block)) {
  pass("METHODS Fig 4 mentions 63 (new Section C PTC- n)")
} else {
  fail("METHODS Fig 4 does NOT mention 63 (new Section C PTC- n)")
}
if (grepl("819", fig4_block)) {
  pass("METHODS Fig 4 mentions 819 (new Section C Control n)")
} else {
  fail("METHODS Fig 4 does NOT mention 819 (new Section C Control n)")
}

if (grepl("1,489|1489", fig4_block)) {
  fail("METHODS Fig 4 still mentions 1,489 (stale OLD Section C PTC+ n)")
} else {
  pass("METHODS Fig 4 does NOT mention 1,489")
}
if (grepl("n = 1,659|n=1,659|1,659/1,286|n = 1,659/1,286", fig4_block)) {
  fail("METHODS Fig 4 still mentions n=1,659 (stale OLD pre-re-intersection universe)")
} else {
  pass("METHODS Fig 4 does NOT mention 1,659 as denominator")
}
if (grepl("1,286|1286", fig4_block)) {
  fail("METHODS Fig 4 still mentions 1,286 (stale OLD Section C Control n)")
} else {
  pass("METHODS Fig 4 does NOT mention 1,286")
}
if (grepl("/163/", fig4_block) || grepl(" 163 NMD", fig4_block)) {
  fail("METHODS Fig 4 still mentions /163/ stale PTC- n")
} else {
  pass("METHODS Fig 4 does NOT use stale /163/ PTC- n")
}

cat("\n=== Supplemental TD2 figure ===\n")
if (grepl("TD2BiasEvidence|consolidated.+supplement|n = 900 occult", fig4_block)) {
  pass("METHODS Fig 4 references the consolidated TD2BiasEvidence supplement")
} else {
  fail("METHODS Fig 4 does NOT reference the consolidated TD2BiasEvidence supplement")
}

cat("\n=== Cross-check: METHODS code references actually exist ===\n")
expected_files <- c(
  "figures/multipanel/figure4_ptcneg_and_model/data_export.R",
  "figures/multipanel/figure4_ptcneg_and_model/figure4_composite.py",
  "figures/multipanel/figure4_ptcneg_and_model/RATIONALE.md",
  "figures/lib/mechanism_class.R"
)
proj_root <- "/Users/petecastaldi/claude_projects/nmd"
for (f in expected_files) {
  ref_in_methods <- grepl(f, methods_txt, fixed = TRUE)
  exists_on_disk <- file.exists(file.path(proj_root, f))
  if (ref_in_methods && exists_on_disk) {
    pass(sprintf("METHODS references %s (exists)", f))
  } else if (ref_in_methods && !exists_on_disk) {
    fail(sprintf("METHODS references %s but file MISSING", f))
  } else if (!ref_in_methods && exists_on_disk) {
    fail(sprintf("METHODS does NOT reference %s (file present)", f))
  } else {
    fail(sprintf("METHODS does not reference %s and file missing", f))
  }
}

cat("\n=================================================\n")
cat(sprintf("PASS 5 SUMMARY: %d passes / %d failures\n", pass_ct, fail_ct))
cat("=================================================\n")
if (fail_ct > 0) {
  cat("\nFailures:\n")
  for (m in fails) cat(" -", m, "\n")
}
