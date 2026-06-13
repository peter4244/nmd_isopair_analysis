#!/usr/bin/env Rscript
# Verification Pass 1 — Factual accuracy for Figure 3
#
# Independently recompute every stated number in Panel A-F methodology files,
# the find/replace document, and the results-to-code map. Don't trust cached
# intermediate values — recompute from raw RDS caches.
#
# Output: PASS/FAIL per claim, with the independent value vs the stated value.

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
ref_atg <- readRDS(file.path(CACHE, "ref_atg_analysis.rds"))
div_c2 <- readRDS(file.path(CACHE, "div_c2_allsamples.rds"))
div_c4 <- readRDS(file.path(CACHE, "div_c4_allsamples.rds"))

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
cat("PASS 1: Factual accuracy verification for Figure 3\n")
cat("================================================================\n\n")

# ============================================================================
# Section A — Gene-eligibility / pair-construction (find/replace doc + results map)
# ============================================================================
cat("[A] Gene eligibility + pair counts (from find/replace doc + Rmd flowchart)\n")
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
# Section C — Ref-AUG-traceable subset (pop_traceable for Panel D)
# ============================================================================
cat("\n[C] pop_traceable = ref-AUG tracing succeeded (Panel D)\n")
TRACEABLE <- c("effectively_ptc", "no_downstream_ejc", "truncated_no_ejc")
ra_c2 <- as.data.table(ref_atg$c2)
ra_c4 <- as.data.table(ref_atg$c4)
trace_c2_ids <- ra_c2$comparator_isoform_id[ra_c2$category %in% TRACEABLE]
trace_c4_ids <- ra_c4$comparator_isoform_id[ra_c4$category %in% TRACEABLE]
pop_trace_c2 <- pop_BC_c2[comparator_isoform_id %in% trace_c2_ids]
pop_trace_c4 <- pop_BC_c4[comparator_isoform_id %in% trace_c4_ids]
check("pop_traceable NMD", 2289, nrow(pop_trace_c2))
check("pop_traceable Control", 1763, nrow(pop_trace_c4))

# Category breakdown for NMD ref_atg
cat("  Ref-AUG category counts:\n")
cat_breakdown <- ra_c2[, .N, by = category]
print(cat_breakdown)
check("NMD effectively_ptc", 1912, sum(ra_c2$category == "effectively_ptc"))
check("NMD no_downstream_ejc", 301, sum(ra_c2$category == "no_downstream_ejc"))
check("NMD truncated_no_ejc", 76, sum(ra_c2$category == "truncated_no_ejc"))
check("NMD ref_atg_lost", 351, sum(ra_c2$category == "ref_atg_lost"))
check("NMD mapping_failed", 33, sum(ra_c2$category == "mapping_failed"))
check("Control effectively_ptc", 288, sum(ra_c4$category == "effectively_ptc"))
check("Control no_downstream_ejc", 825, sum(ra_c4$category == "no_downstream_ejc"))
check("Control truncated_no_ejc", 650, sum(ra_c4$category == "truncated_no_ejc"))
check("Control ref_atg_lost", 770, sum(ra_c4$category == "ref_atg_lost"))
check("Control mapping_failed", 46, sum(ra_c4$category == "mapping_failed"))

# Headline percentages
ptc_rate_nmd <- round(100 * sum(ra_c2$category == "effectively_ptc") / nrow(pop_trace_c2), 1)
ptc_rate_ctrl <- round(100 * sum(ra_c4$category == "effectively_ptc") / nrow(pop_trace_c4), 1)
check("NMD PTC rate (1912/2289)", 83.5, ptc_rate_nmd)
check("Control PTC rate (288/1763)", 16.3, ptc_rate_ctrl)
# Fold enrichment
fold <- round(ptc_rate_nmd / ptc_rate_ctrl, 1)
check("Fold enrichment (NMD/Control PTC rate)", 5.1, fold)

# ============================================================================
# Section D — pop_ptc_plus (Panels E, F)
# ============================================================================
cat("\n[D] pop_ptc_plus = effectively_ptc subset (Panels E, F)\n")
pop_ptc_c2 <- pop_BC_c2[comparator_isoform_id %in%
                         ra_c2$comparator_isoform_id[ra_c2$category == "effectively_ptc"]]
pop_ptc_c4 <- pop_BC_c4[comparator_isoform_id %in%
                         ra_c4$comparator_isoform_id[ra_c4$category == "effectively_ptc"]]
check("pop_ptc_plus NMD", 1912, nrow(pop_ptc_c2))
check("pop_ptc_plus Control", 288, nrow(pop_ptc_c4))

# Split by original_ptc
n_orig_T <- sum(ra_c2$category == "effectively_ptc" & ra_c2$original_ptc == TRUE)
n_orig_F <- sum(ra_c2$category == "effectively_ptc" & ra_c2$original_ptc == FALSE)
check("effectively_ptc x original_ptc=TRUE", 1012, n_orig_T)
check("effectively_ptc x original_ptc=FALSE", 900, n_orig_F)

# ============================================================================
# Section E — Panel D distance computation sanity
# ============================================================================
cat("\n[E] Panel D distance stats\n")
panelD_tsv <- fread("data/panelD_stop_codon_distance.tsv")
check("Panel D NMD n", 2289, sum(panelD_tsv$comparison == "NMD"))
check("Panel D Control n", 1763, sum(panelD_tsv$comparison == "Control"))
# Median distance (should be positive for NMD, very negative for Control)
med_nmd <- median(panelD_tsv[comparison == "NMD"]$distance)
med_ctrl <- median(panelD_tsv[comparison == "Control"]$distance)
cat(sprintf("  Panel D NMD median distance: %.0f (expect strongly positive)\n", med_nmd))
cat(sprintf("  Panel D Control median distance: %.0f (expect negative)\n", med_ctrl))
n_pos_nmd <- sum(panelD_tsv[comparison == "NMD"]$distance > 50)
pct_pos_nmd <- round(100 * n_pos_nmd / nrow(pop_trace_c2), 1)
cat(sprintf("  Panel D NMD with distance > 50: %d (%.1f%%) -- should match PTC rate 83.5%%\n",
            n_pos_nmd, pct_pos_nmd))
check("Panel D NMD pct_pos (== effectively_ptc rate)", 83.5, pct_pos_nmd, tol = 1)

# ============================================================================
# Section F — Panel E + F attribution counts
# ============================================================================
cat("\n[F] Panel E/F attribution counts\n")
panelE_tsv <- fread("data/panelE_ptc_event_attribution.tsv")
panelF_tsv <- fread("data/panelF_mechanism_breakdown.tsv")
n_ptc_attr <- sum(panelE_tsv$n_ptc_events)
ctrl_total <- sum(panelE_tsv$n_ctrl_events)
check("Panel E n_ptc_attr (NMD PTC-causing events)", 1812, n_ptc_attr)
check("Panel E ctrl_total (Control all-events baseline)", 4525, ctrl_total)
check("Panel F total stacked count", 1812, sum(panelF_tsv$n))

# SE dominance check (from manuscript prose: "skipped exon... ~50%")
se_pct <- panelE_tsv$pct_of_ptc[panelE_tsv$event_type == "SE"]
cat(sprintf("  Panel E SE %% of PTC-causing: %s (manuscript: 'majority of PTCs caused by exon skipping')\n", se_pct))

# Panel F mechanism totals
mech_totals <- panelF_tsv[, .(n = sum(n)), by = mechanism]
print(mech_totals)
check("Panel F Frameshift total in F", 1117, sum(panelF_tsv$n[panelF_tsv$mechanism == "Frameshift"]))
check("Panel F In-frame stop total in F", 587, sum(panelF_tsv$n[panelF_tsv$mechanism == "In-frame stop"]))
check("Panel F 3'UTR splice total in F", 108, sum(panelF_tsv$n[panelF_tsv$mechanism == "3'UTR splice"]))

# ============================================================================
# Section G — Panel B + C
# ============================================================================
cat("\n[G] Panel B + C sanity\n")
panelB_tsv <- fread("data/panelB_sequence_similarity.tsv")
check("Panel B NMD n", 3009, sum(panelB_tsv$comparison == "NMD"))
check("Panel B Control n", 3009, sum(panelB_tsv$comparison == "Control"))
panelC_tsv <- fread("data/panelC_event_prevalence.tsv")
check("Panel C # event types", 10, nrow(panelC_tsv))
# SE event count from raw data, sanity-check Panel C TSV
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
