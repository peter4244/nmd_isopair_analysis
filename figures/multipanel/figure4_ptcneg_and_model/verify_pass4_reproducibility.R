#!/usr/bin/env Rscript
# Pass 4: Reproducibility & completeness (NEW 2x2 layout, 2026-06-15).
# Trace each panel: script → data TSV → data_export.R chunk → upstream RDS.
# Confirm data_export.R Section C performs the re-intersection on
# (gene_id, reference_isoform_id) AFTER the TR category filter, and that the
# re-execution produces n = 756/63/819.

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/Users/petecastaldi/claude_projects/nmd/figures/multipanel/figure4_ptcneg_and_model"
CACHE_DIR <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr"
fails <- list(); pass_ct <- 0L; fail_ct <- 0L
pass <- function(m) { pass_ct <<- pass_ct + 1L; cat("PASS  ", m, "\n", sep = "") }
fail <- function(m) { fail_ct <<- fail_ct + 1L; fails[[length(fails) + 1L]] <<- m; cat("FAIL  ", m, "\n", sep = "") }

cat("=== Files presence ===\n")
required <- c(
  "data_export.R",
  "figure4_panelA_5utr_length_all3enst.py",
  "figure4_panelB_longest_5utr_orf_all3enst.py",
  "figure4_panelC_5utr_length_refaug.py",
  "figure4_panelD_longest_5utr_orf_refaug.py",
  "figure4_composite.py",
  "figure4_composite.png", "figure4_composite.pdf",
  "RATIONALE.md", "figure4_composite_legend.md",
  "data/panelA_5utr_length_long.tsv",
  "data/panelA_5utr_length_descriptives.tsv",
  "data/panelA_5utr_length_pairwise.tsv",
  "data/panelB_longest_5utr_orf_long.tsv",
  "data/panelB_longest_5utr_orf_descriptives.tsv",
  "data/panelB_longest_5utr_orf_pairwise.tsv",
  "data/panelC_5utr_length_long.tsv",
  "data/panelC_5utr_length_descriptives.tsv",
  "data/panelC_5utr_length_pairwise.tsv",
  "data/panelD_longest_5utr_orf_long.tsv",
  "data/panelD_longest_5utr_orf_descriptives.tsv",
  "data/panelD_longest_5utr_orf_pairwise.tsv"
)
for (r in required) {
  if (file.exists(file.path(ROOT, r))) pass(sprintf("file exists: %s", r))
  else fail(sprintf("missing: %s", r))
}

cat("\n=== Upstream RDS presence ===\n")
upstream <- c(
  "analysis_cache/ref_atg_analysis.rds",
  "analysis_cache/utr5_features_refaug.rds",
  "analysis_cache/utr5_features_all.rds",
  "cds.rds", "structures.rds",
  "profiles_c2_allsamples.rds", "profiles_c4_allsamples.rds"
)
for (u in upstream) {
  if (file.exists(file.path(CACHE_DIR, u))) pass(sprintf("upstream: %s", u))
  else fail(sprintf("upstream MISSING: %s", u))
}

cat("\n=== mechanism_class helper sourced ===\n")
mc_path <- "/Users/petecastaldi/claude_projects/nmd/figures/lib/mechanism_class.R"
if (file.exists(mc_path)) pass("mechanism_class.R exists") else fail("mechanism_class.R missing")

cat("\n=== data_export.R Section C re-intersection check ===\n")
de_lines <- readLines(file.path(ROOT, "data_export.R"))
de_text <- paste(de_lines, collapse = "\n")
# Look for the re-intersection step on (gene_id, reference_isoform_id) AFTER TR filter
if (grepl("c2_tr <- c2_full\\[category %in% TR\\]", de_text) &&
    grepl("c4_tr <- c4_full\\[category %in% TR\\]", de_text) &&
    grepl("shared_keys <- merge\\(\\s*unique\\(c2_tr\\[, \\.\\(gene_id, reference_isoform_id\\)\\]\\)", de_text)) {
  pass("data_export.R applies TR filter to each side first, then re-intersects on (gene_id, reference_isoform_id)")
} else {
  fail("data_export.R re-intersection step on (gene_id, reference_isoform_id) AFTER TR filter is not found in expected form")
}

# Check the TR vector definition
if (grepl('TR <- c\\("effectively_ptc", "no_downstream_ejc", "truncated_no_ejc"\\)', de_text)) {
  pass("TR category vector includes the 3 ref-AUG-traceable categories")
} else {
  fail("TR category vector definition not found in expected form")
}

cat("\n=== Re-execute Section C universe construction independently ===\n")
ref_atg <- readRDS(file.path(CACHE_DIR, "analysis_cache/ref_atg_analysis.rds"))
profiles_c2 <- as.data.table(readRDS(file.path(CACHE_DIR, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(CACHE_DIR, "profiles_c4_allsamples.rds")))
source(mc_path)

is_enst <- function(x) grepl("^ENST", x)
make_key <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")
shared <- intersect(unique(make_key(profiles_c2)), unique(make_key(profiles_c4)))
pop_BC_c2 <- profiles_c2[make_key(profiles_c2) %in% shared]
pop_BC_c4 <- profiles_c4[make_key(profiles_c4) %in% shared]
sec_C_c2_base <- pop_BC_c2[is_enst(reference_isoform_id)]
sec_C_c4_base <- pop_BC_c4[is_enst(reference_isoform_id)]
c2_full <- as.data.table(ref_atg$c2)[comparator_isoform_id %in% sec_C_c2_base$comparator_isoform_id]
c4_full <- as.data.table(ref_atg$c4)[comparator_isoform_id %in% sec_C_c4_base$comparator_isoform_id]
TR <- c("effectively_ptc", "no_downstream_ejc", "truncated_no_ejc")
c2_tr <- c2_full[category %in% TR]
c4_tr <- c4_full[category %in% TR]
shared_keys <- merge(
  unique(c2_tr[, .(gene_id, reference_isoform_id)]),
  unique(c4_tr[, .(gene_id, reference_isoform_id)]),
  by = c("gene_id", "reference_isoform_id"))
c2_re <- merge(c2_tr, shared_keys, by = c("gene_id", "reference_isoform_id"))
c4_re <- merge(c4_tr, shared_keys, by = c("gene_id", "reference_isoform_id"))
c2_re[, mech5 := mechanism_class(category, comp_orf_length, ref_orf_length)]
c2_re[, mech4 := mechanism_class_4(mech5)]
c2_re[, group3 := fcase(
  mech4 == "NMD+/PTC+", "NMD+/PTC+",
  mech4 %in% c("NMD+/PTC- ORF match", "NMD+/PTC- ORF diff"), "NMD+/PTC-",
  default = "Other"
)]
n_ptc_pos <- nrow(c2_re[group3 == "NMD+/PTC+"])
n_ptc_neg <- nrow(c2_re[group3 == "NMD+/PTC-"])
n_ctrl    <- nrow(c4_re)
cat(sprintf("Recomputed Section C: NMD+/PTC+=%d  NMD+/PTC-=%d  Control=%d\n",
            n_ptc_pos, n_ptc_neg, n_ctrl))
if (n_ptc_pos == 756L) pass("Section C re-intersection: NMD+/PTC+ = 756") else fail(sprintf("Section C NMD+/PTC+ = %d (expected 756)", n_ptc_pos))
if (n_ptc_neg == 63L)  pass("Section C re-intersection: NMD+/PTC- = 63")  else fail(sprintf("Section C NMD+/PTC- = %d (expected 63)", n_ptc_neg))
if (n_ctrl    == 819L) pass("Section C re-intersection: Control = 819")  else fail(sprintf("Section C Control = %d (expected 819)", n_ctrl))

# Also verify NMD and Control share the exact same (gene_id, reference_isoform_id) universe
nmd_keys <- unique(c2_re[, .(gene_id, reference_isoform_id)])
ctrl_keys <- unique(c4_re[, .(gene_id, reference_isoform_id)])
if (nrow(nmd_keys) == nrow(ctrl_keys) && nrow(nmd_keys) == nrow(shared_keys)) {
  pass(sprintf("Section C: NMD and Control share same %d (gene, ref) keys", nrow(shared_keys)))
} else {
  fail(sprintf("Section C key mismatch: NMD=%d Ctrl=%d shared=%d", nrow(nmd_keys), nrow(ctrl_keys), nrow(shared_keys)))
}

cat("\n=== Each panel script reads from the correct TSV ===\n")
panels <- list(
  "figure4_panelA_5utr_length_all3enst.py"      = "data/panelA_5utr_length",
  "figure4_panelB_longest_5utr_orf_all3enst.py" = "data/panelB_longest_5utr_orf",
  "figure4_panelC_5utr_length_refaug.py"        = "data/panelC_5utr_length",
  "figure4_panelD_longest_5utr_orf_refaug.py"   = "data/panelD_longest_5utr_orf"
)
for (script in names(panels)) {
  src <- paste(readLines(file.path(ROOT, script)), collapse = "\n")
  slug <- panels[[script]]
  if (grepl(paste0(basename(slug), "_long\\.tsv"), src) &&
      grepl(paste0(basename(slug), "_pairwise\\.tsv"), src) &&
      grepl(paste0(basename(slug), "_descriptives\\.tsv"), src)) {
    pass(sprintf("%s reads %s_{long,descriptives,pairwise}.tsv", script, basename(slug)))
  } else {
    fail(sprintf("%s does not consistently read %s_*", script, basename(slug)))
  }
}

cat("\n=================================================\n")
cat(sprintf("PASS 4 SUMMARY: %d passes / %d failures\n", pass_ct, fail_ct))
cat("=================================================\n")
if (fail_ct > 0) {
  cat("\nFailures:\n")
  for (m in fails) cat(" -", m, "\n")
  quit(status = 1)
}
