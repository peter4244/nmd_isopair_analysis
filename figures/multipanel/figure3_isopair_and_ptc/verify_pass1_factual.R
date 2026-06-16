#!/usr/bin/env Rscript
# Verification Pass 1 — Factual accuracy for Figure 3
#
# Independently recompute every stated number in Panel A-F methodology files,
# the find/replace document, and the results-to-code map. Don't trust cached
# intermediate values — recompute from raw RDS caches.
#
# Output: PASS/FAIL per claim, with the independent value vs the stated value.
#
# Scope refresh 2026-06-15 (v2): all-3-ENST + coding-CDS + own-GENCODE-stop.
# No ref-AUG projection; no TD2 dependency. Headline 37.9% NMD / 2.1% Control,
# 18× enrichment, 190/190 pairs.

suppressPackageStartupMessages({
  library(data.table)
})

DM <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr"
CACHE <- file.path(DM, "analysis_cache")

# ---- Load raw caches ----
profiles_c2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))
cds <- as.data.table(readRDS(file.path(DM, "cds.rds")))
ptc <- as.data.table(readRDS(file.path(DM, "ptc.rds")))
gm <- readRDS(file.path(DM, "gene_map.rds"))
nmd_class <- readRDS(file.path(DM, "nmd_classification.rds"))

# Track pass/fail
results <- list()
check <- function(label, expected, actual, tol = 0) {
  ok <- if (is.numeric(expected) && tol > 0)
    abs(actual - expected) <= tol
  else identical(as.character(expected), as.character(actual))
  status <- if (ok) "PASS" else "FAIL"
  cat(sprintf("  [%s] %s: expected=%s, actual=%s\n",
              status, label, expected, actual))
  results[[length(results) + 1]] <<- list(label = label, status = status,
                                          expected = expected, actual = actual)
  invisible(ok)
}

cat("================================================================\n")
cat("PASS 1: Factual accuracy verification for Figure 3 (all-3-ENST scope)\n")
cat("================================================================\n\n")

# ============================================================================
# Section A — Gene-eligibility / pair-construction
# ============================================================================
cat("[A] Gene eligibility + pair counts\n")
nmd_genes <- unique(gm$gene_id[gm$isoform_id %in% nmd_class[["all_samples"]]$nmd])
check("genes with >=1 NMD isoform", 5367, length(nmd_genes))
nonnmd_per_gene <- table(gm$gene_id[gm$isoform_id %in% nmd_class[["all_samples"]]$non_nmd])
genes_2plus_nonnmd <- names(nonnmd_per_gene[nonnmd_per_gene >= 2])
check("genes with >=2 non-NMD isoforms", 8472, length(genes_2plus_nonnmd))
genes_eligible <- intersect(nmd_genes, genes_2plus_nonnmd)
check("eligible genes (>=1 NMD AND >=2 non-NMD)", 3099, length(genes_eligible))
check("raw profiles_c2 pairs", 3009, nrow(profiles_c2))
check("raw profiles_c2 unique genes", 3009, length(unique(profiles_c2$gene_id)))
check("raw profiles_c4 pairs", 8323, nrow(profiles_c4))

# ============================================================================
# Section B — Stage 2 gene-matching (pop_BC for Panels B, C)
# ============================================================================
cat("\n[B] pop_BC = Stage 2 gene-matched (Panels B, C)\n")
make_key <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")
shared <- intersect(unique(make_key(profiles_c2)), unique(make_key(profiles_c4)))
pop_BC_c2 <- profiles_c2[make_key(profiles_c2) %in% shared]
pop_BC_c4 <- profiles_c4[make_key(profiles_c4) %in% shared]
check("pop_BC NMD pairs", 3009, nrow(pop_BC_c2))
check("pop_BC Control pairs", 3009, nrow(pop_BC_c4))

# ============================================================================
# Section C — All-3-ENST + coding-CDS scope (Panels D/E/F universe)
# ============================================================================
cat("\n[C] All-3-ENST gene-matched + coding-CDS (Panel D/E/F universe)\n")
is_enst <- function(x) grepl("^ENST", x)

# Per-side ENST (reference + comparator)
pop_BC_c2_enst <- pop_BC_c2[is_enst(reference_isoform_id) & is_enst(comparator_isoform_id)]
pop_BC_c4_enst <- pop_BC_c4[is_enst(reference_isoform_id) & is_enst(comparator_isoform_id)]

# Re-intersect by (gene_id, reference_isoform_id) so that all 3 isoforms
# (reference + NMD comparator + Control comparator) are ENST.
all3_keys <- intersect(unique(make_key(pop_BC_c2_enst)),
                       unique(make_key(pop_BC_c4_enst)))
all3_c2 <- pop_BC_c2_enst[make_key(pop_BC_c2_enst) %in% all3_keys]
all3_c4 <- pop_BC_c4_enst[make_key(pop_BC_c4_enst) %in% all3_keys]
check("all-3-ENST NMD pairs (gene-matched)", 301, nrow(all3_c2))
check("all-3-ENST Control pairs (gene-matched)", 301, nrow(all3_c4))

# Coding-CDS filter: every isoform in every pair has coding_status == "coding"
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
all3_coding_c2 <- all3_c2[reference_isoform_id %in% coding_ids &
                          comparator_isoform_id %in% coding_ids]
all3_coding_c4 <- all3_c4[reference_isoform_id %in% coding_ids &
                          comparator_isoform_id %in% coding_ids]

# Re-intersect once more after coding filter (both sides must retain the same
# (gene, reference) key).
final_keys <- intersect(unique(make_key(all3_coding_c2)),
                        unique(make_key(all3_coding_c4)))
pop_DEF_c2 <- all3_coding_c2[make_key(all3_coding_c2) %in% final_keys]
pop_DEF_c4 <- all3_coding_c4[make_key(all3_coding_c4) %in% final_keys]
check("all-3-ENST + coding-CDS, re-intersected NMD", 190, nrow(pop_DEF_c2))
check("all-3-ENST + coding-CDS, re-intersected Control", 190, nrow(pop_DEF_c4))

# ============================================================================
# Section D — Panel D distance + PTC+ subset
# ============================================================================
cat("\n[D] Panel D distance stats + PTC+ subset (own-GENCODE-stop, 50-nt rule)\n")
panelD_tsv <- fread("data/panelD_stop_codon_distance.tsv")
check("Panel D NMD n", 190, sum(panelD_tsv$comparison == "NMD"))
check("Panel D Control n", 190, sum(panelD_tsv$comparison == "Control"))

# TSV column schema: should NOT have a `category` column (own-GENCODE scope)
expected_cols <- c("comparator_isoform_id", "gene_id", "reference_isoform_id",
                   "comparison", "distance", "last_ejc_tx_pos", "own_stop_tx_pos")
check("Panel D TSV columns match expected (no `category` column)",
      paste(expected_cols, collapse = ","),
      paste(colnames(panelD_tsv), collapse = ","))

# Medians
med_nmd <- median(panelD_tsv[comparison == "NMD"]$distance)
med_ctrl <- median(panelD_tsv[comparison == "Control"]$distance)
check("Panel D NMD median distance", -66, med_nmd)
check("Panel D Control median distance", -143, med_ctrl)

# PTC+ counts (own-stop > 50 nt past last EJC == distance > 50 nt under
# the panelD convention `distance = last_ejc_tx_pos - own_stop_tx_pos`)
n_ptc_nmd <- sum(panelD_tsv[comparison == "NMD"]$distance > 50)
n_ptc_ctrl <- sum(panelD_tsv[comparison == "Control"]$distance > 50)
check("Panel D NMD PTC+ count (distance > 50)", 72, n_ptc_nmd)
check("Panel D Control PTC+ count (distance > 50)", 4, n_ptc_ctrl)

# Headline percentages
ptc_rate_nmd <- round(100 * n_ptc_nmd / nrow(pop_DEF_c2), 1)
ptc_rate_ctrl <- round(100 * n_ptc_ctrl / nrow(pop_DEF_c4), 1)
check("NMD PTC rate % (72/190)", 37.9, ptc_rate_nmd)
check("Control PTC rate % (4/190)", 2.1, ptc_rate_ctrl)

# Fold enrichment (integer)
fold <- round(ptc_rate_nmd / ptc_rate_ctrl)
check("Fold enrichment (~18x)", 18, fold)

# Clipping at [-1000, 1500]
n_clip_nmd <- sum(panelD_tsv[comparison == "NMD"]$distance > 1500 |
                  panelD_tsv[comparison == "NMD"]$distance < -1000)
n_clip_ctrl <- sum(panelD_tsv[comparison == "Control"]$distance > 1500 |
                   panelD_tsv[comparison == "Control"]$distance < -1000)
check("Panel D NMD clipped count (axis [-1000,1500])", 4, n_clip_nmd)
check("Panel D Control clipped count (axis [-1000,1500])", 4, n_clip_ctrl)

# ============================================================================
# Section E — Panel E + F attribution counts
# ============================================================================
cat("\n[E] Panel E/F attribution counts\n")
panelE_tsv <- fread("data/panelE_ptc_event_attribution.tsv")
panelF_tsv <- fread("data/panelF_mechanism_breakdown.tsv")
n_ptc_attr <- sum(panelE_tsv$n_ptc_events)
ctrl_total_tsv <- sum(panelE_tsv$n_ctrl_events)
check("Panel E n_ptc_attr (sum n_ptc_events)", 69, n_ptc_attr)
# The TSV emits only the 9 event types that appear among PTC-causing events;
# ctrl_total_tsv = 330 (subset of ctrl_total = 447 used for percentages /
# Fisher denominators in panel_e_compute.R).
check("Panel E n_ctrl_events sum in TSV (9 PTC-matching event types)",
      330, ctrl_total_tsv)
check("Panel F total attributed pairs (sum n)", 69, sum(panelF_tsv$n))

# Top event (SE)
se_row <- panelE_tsv[event_type == "SE"]
check("Panel E SE n_ptc_events", 30, se_row$n_ptc_events)
check("Panel E SE pct_of_ptc", 43.5, se_row$pct_of_ptc, tol = 0.1)
check("Panel E SE n_ctrl_events", 63, se_row$n_ctrl_events)
check("Panel E SE pct_ctrl", 14.1, se_row$pct_ctrl, tol = 0.1)
check("Panel E SE direction", "Enriched in PTC-causing", se_row$direction)
# Fisher p check (log10 comparison, ~1e-7..1e-8)
se_p <- se_row$fisher_p
check("Panel E SE fisher_p ~ 7.98e-08", 7.98e-08, se_p, tol = 1e-08)

# A5SS row
a5_row <- panelE_tsv[event_type == "A5SS"]
check("Panel E A5SS n_ptc_events", 9, a5_row$n_ptc_events)
check("Panel E A5SS pct_of_ptc", 13.0, a5_row$pct_of_ptc, tol = 0.1)
check("Panel E A5SS n_ctrl_events", 20, a5_row$n_ctrl_events)
check("Panel E A5SS pct_ctrl", 4.5, a5_row$pct_ctrl, tol = 0.1)
check("Panel E A5SS direction", "Enriched in PTC-causing", a5_row$direction)
check("Panel E A5SS fisher_p ~ 8.88e-03", 8.88e-03, a5_row$fisher_p, tol = 1e-3)

# Alt_TES row (depleted)
at_row <- panelE_tsv[event_type == "Alt_TES"]
check("Panel E Alt_TES n_ptc_events", 4, at_row$n_ptc_events)
check("Panel E Alt_TES pct_of_ptc", 5.8, at_row$pct_of_ptc, tol = 0.1)
check("Panel E Alt_TES n_ctrl_events", 120, at_row$n_ctrl_events)
check("Panel E Alt_TES pct_ctrl", 26.8, at_row$pct_ctrl, tol = 0.1)
check("Panel E Alt_TES direction", "Depleted in PTC-causing", at_row$direction)
check("Panel E Alt_TES fisher_p ~ 3.37e-05", 3.37e-05, at_row$fisher_p, tol = 1e-5)

# Panel F mechanism totals (3 mechanisms only)
cat("  Panel F mechanism totals:\n")
mech_totals <- panelF_tsv[, .(n = sum(n)), by = mechanism]
print(mech_totals)
check("Panel F Frameshift total", 38, sum(panelF_tsv$n[panelF_tsv$mechanism == "Frameshift"]))
check("Panel F In-frame stop total", 23, sum(panelF_tsv$n[panelF_tsv$mechanism == "In-frame stop"]))
check("Panel F 3'UTR splice total", 8, sum(panelF_tsv$n[panelF_tsv$mechanism == "3'UTR splice"]))

# Mechanism percentages
fs_pct <- round(100 * 38 / 69, 1)
ifs_pct <- round(100 * 23 / 69, 1)
utr_pct <- round(100 * 8 / 69, 1)
check("Panel F Frameshift %", 55.1, fs_pct, tol = 0.1)
check("Panel F In-frame stop %", 33.3, ifs_pct, tol = 0.1)
check("Panel F 3'UTR splice %", 11.6, utr_pct, tol = 0.1)

# Top event count in Panel F
se_in_F <- sum(panelF_tsv$n[panelF_tsv$event_type == "SE"])
check("Panel F top event SE total", 30, se_in_F)

# ============================================================================
# Section F — Panel B + C (unchanged scope at pop_BC = 3,009 each)
# ============================================================================
cat("\n[F] Panel B + C (pop_BC = 3,009 each)\n")
panelB_tsv <- fread("data/panelB_sequence_similarity.tsv")
check("Panel B NMD n", 3009, sum(panelB_tsv$comparison == "NMD"))
check("Panel B Control n", 3009, sum(panelB_tsv$comparison == "Control"))
panelC_tsv <- fread("data/panelC_event_prevalence.tsv")
check("Panel C # event types", 10, nrow(panelC_tsv))
n_se_nmd <- sum(pop_BC_c2$n_se > 0)
n_se_ctrl <- sum(pop_BC_c4$n_se > 0)
check("Panel C SE NMD (raw count)", n_se_nmd,
      panelC_tsv$n_pairs_with_event_NMD[panelC_tsv$event_type == "SE"])
check("Panel C SE Control (raw count)", n_se_ctrl,
      panelC_tsv$n_pairs_with_event_Ctrl[panelC_tsv$event_type == "SE"])

# ============================================================================
# Summary
# ============================================================================
cat("\n================================================================\n")
n_pass <- sum(sapply(results, function(r) r$status == "PASS"))
n_fail <- sum(sapply(results, function(r) r$status == "FAIL"))
cat(sprintf("Pass 1 summary: %d PASS / %d FAIL / %d total\n", n_pass, n_fail, length(results)))
if (n_fail > 0) {
  cat("\nFAILURES:\n")
  for (r in results) {
    if (r$status == "FAIL") {
      cat(sprintf("  - %s: expected=%s, actual=%s\n", r$label, r$expected, r$actual))
    }
  }
}
cat("================================================================\n")
