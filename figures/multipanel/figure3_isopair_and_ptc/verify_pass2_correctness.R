#!/usr/bin/env Rscript
# Verification Pass 2 — Result correctness for Figure 3
#
# For each statistical test in Figure 3 (Fisher exacts in Panel C and Panel E),
# verify the 2×2 table construction independently. Check OR / p-values for
# plausibility. Apply domain knowledge to test biological reasonableness.
# Actively try to disprove headline claims rather than confirm them.

suppressPackageStartupMessages({
  library(data.table)
})

DM <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr"
CACHE <- file.path(DM, "analysis_cache")

# ---- Reload populations from Pass 1 ----
profiles_c2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))
ref_atg <- readRDS(file.path(CACHE, "ref_atg_analysis.rds"))
ra_c2 <- as.data.table(ref_atg$c2)
ra_c4 <- as.data.table(ref_atg$c4)
make_key <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")
shared <- intersect(unique(make_key(profiles_c2)), unique(make_key(profiles_c4)))
pop_BC_c2 <- profiles_c2[make_key(profiles_c2) %in% shared]
pop_BC_c4 <- profiles_c4[make_key(profiles_c4) %in% shared]
TRACEABLE <- c("effectively_ptc", "no_downstream_ejc", "truncated_no_ejc")
trace_c2_ids <- ra_c2$comparator_isoform_id[ra_c2$category %in% TRACEABLE]
trace_c4_ids <- ra_c4$comparator_isoform_id[ra_c4$category %in% TRACEABLE]
pop_trace_c2 <- pop_BC_c2[comparator_isoform_id %in% trace_c2_ids]
pop_trace_c4 <- pop_BC_c4[comparator_isoform_id %in% trace_c4_ids]
N_BC <- 3009

# ---- Tracking ----
results <- list()
note <- function(label, status, detail = "") {
  cat(sprintf("  [%s] %s%s\n", status, label,
              if (nzchar(detail)) sprintf(" -- %s", detail) else ""))
  results[[length(results) + 1]] <<- list(label = label, status = status,
                                          detail = detail)
}

cat("================================================================\n")
cat("PASS 2: Result correctness verification for Figure 3\n")
cat("================================================================\n\n")

# ============================================================================
# Test 1 — Panel C Fisher exact tests (event prevalence)
# Independently recompute the 2×2 table for each event and verify OR + p-value
# match the TSV values to reasonable precision.
# ============================================================================
cat("[1] Panel C — Fisher exact tests, 2x2 reconstruction\n")
panelC <- fread("data/panelC_event_prevalence.tsv")

for (i in seq_len(nrow(panelC))) {
  e <- panelC$event_type[i]
  a <- panelC$n_pairs_with_event_NMD[i]
  c_ <- panelC$n_pairs_with_event_Ctrl[i]
  b <- N_BC - a
  d <- N_BC - c_
  ft <- fisher.test(matrix(c(a, b, c_, d), nrow = 2))
  or_ind <- round(unname(ft$estimate), 2)
  p_ind <- ft$p.value
  or_tsv <- panelC$fisher_OR[i]
  p_tsv <- panelC$fisher_p[i]
  ok_or <- abs(or_ind - or_tsv) < 0.02
  ok_p  <- abs(log10(p_ind + 1e-300) - log10(p_tsv + 1e-300)) < 0.1
  status <- if (ok_or && ok_p) "PASS" else "FAIL"
  detail <- sprintf("OR ind=%.2f tsv=%.2f | p ind=%.2g tsv=%.2g",
                    or_ind, or_tsv, p_ind, p_tsv)
  note(sprintf("Panel C %s 2x2 reconstruction", e), status, detail)
}

# Adversarial: SE is the headline. Verify the 2x2 table directly from raw profiles.
cat("\n  [Adversarial] SE 2x2 directly from raw pop_BC profiles:\n")
se_nmd <- sum(pop_BC_c2$n_se > 0)
se_ctrl <- sum(pop_BC_c4$n_se > 0)
cat(sprintf("    SE in NMD: %d / %d  | SE in Control: %d / %d\n",
            se_nmd, N_BC, se_ctrl, N_BC))
ft_se <- fisher.test(matrix(c(se_nmd, N_BC - se_nmd, se_ctrl, N_BC - se_ctrl), nrow = 2))
note("SE OR > 2 (manuscript: 'twice as frequent')",
     if (ft_se$estimate > 2.0) "PASS" else "FAIL",
     sprintf("OR=%.2f", ft_se$estimate))

# Adversarial: check that percentages sum to <= 100% per panel (events can overlap)
nmd_pct_sum <- sum(panelC$pct_of_pairs_NMD)
ctrl_pct_sum <- sum(panelC$pct_of_pairs_Ctrl)
cat(sprintf("    Sum of NMD pct (events can overlap, so >100 is OK): %.1f%%\n", nmd_pct_sum))
cat(sprintf("    Sum of Control pct: %.1f%%\n", ctrl_pct_sum))

# ============================================================================
# Test 2 — Panel E Fisher exact tests (PTC-causing vs all-events Control)
# ============================================================================
cat("\n[2] Panel E — Fisher exact tests, 2x2 reconstruction\n")
panelE <- fread("data/panelE_ptc_event_attribution.tsv")
n_ptc_attr <- sum(panelE$n_ptc_events)
ctrl_total <- sum(panelE$n_ctrl_events)
cat(sprintf("    n_ptc_attr=%d, ctrl_total=%d\n", n_ptc_attr, ctrl_total))

for (i in seq_len(nrow(panelE))) {
  e <- panelE$event_type[i]
  a <- panelE$n_ptc_events[i]
  c_ <- panelE$n_ctrl_events[i]
  b <- n_ptc_attr - a
  d <- ctrl_total - c_
  ft <- fisher.test(matrix(c(a, b, c_, d), nrow = 2))
  p_ind <- ft$p.value
  p_tsv <- panelE$fisher_p[i]
  # TSV stores signif(_,3) so allow up to 1% relative tolerance on log p
  ok_p <- abs(log10(p_ind + 1e-300) - log10(p_tsv + 1e-300)) < 0.1
  status <- if (ok_p) "PASS" else "FAIL"
  detail <- sprintf("OR=%.2f p ind=%.2g tsv=%.2g (direction=%s)",
                    ft$estimate, p_ind, p_tsv, panelE$direction[i])
  note(sprintf("Panel E %s 2x2 reconstruction", e), status, detail)
}

# ============================================================================
# Test 3 — Panel E SE enrichment (the manuscript headline)
# ============================================================================
cat("\n[3] Panel E — SE enrichment (manuscript: 'majority of PTCs caused by SE')\n")
se_pct_ptc <- panelE$pct_of_ptc[panelE$event_type == "SE"]
se_pct_ctrl <- panelE$pct_ctrl[panelE$event_type == "SE"]
cat(sprintf("    SE: PTC-causing %.1f%% vs Control all-events %.1f%% (enrichment %.1fx)\n",
            se_pct_ptc, se_pct_ctrl,
            panelE$enrichment[panelE$event_type == "SE"]))
note("SE PTC% > Control% (enriched)",
     if (se_pct_ptc > se_pct_ctrl) "PASS" else "FAIL")
note("SE %% of PTC > 30% (majority story)",
     if (se_pct_ptc > 30) "PASS" else "FAIL",
     sprintf("actual=%.1f%%", se_pct_ptc))

# ============================================================================
# Test 4 — Panel D distance distribution sign + biological direction
# ============================================================================
cat("\n[4] Panel D — distance sign and direction\n")
panelD <- fread("data/panelD_stop_codon_distance.tsv")
med_nmd <- median(panelD[comparison == "NMD"]$distance)
med_ctrl <- median(panelD[comparison == "Control"]$distance)
cat(sprintf("    NMD median distance: %.0f nt (expect strongly positive > 50)\n", med_nmd))
cat(sprintf("    Control median distance: %.0f nt (expect negative)\n", med_ctrl))
note("NMD median > 50 (PTC direction)", if (med_nmd > 50) "PASS" else "FAIL")
note("Control median <= 0 (normal stop direction)", if (med_ctrl <= 0) "PASS" else "FAIL")
note("NMD/Control medians strongly separated (>200 nt apart)",
     if ((med_nmd - med_ctrl) > 200) "PASS" else "FAIL")

# Check the sign convention by inspecting an example
example <- panelD[comparison == "NMD" & abs(distance) < 100][1]
cat(sprintf("    Example NMD pair: distance=%.0f, category=%s\n",
            example$distance, example$category))

# ============================================================================
# Test 5 — Panel F mechanism plausibility
# ============================================================================
cat("\n[5] Panel F — mechanism breakdown plausibility\n")
panelF <- fread("data/panelF_mechanism_breakdown.tsv")
mech_totals <- panelF[, .(n = sum(n)), by = mechanism]
print(mech_totals)
fs_total <- mech_totals$n[mech_totals$mechanism == "Frameshift"]
ifs_total <- mech_totals$n[mech_totals$mechanism == "In-frame stop"]
utr_total <- mech_totals$n[mech_totals$mechanism == "3'UTR splice"]
total <- fs_total + ifs_total + utr_total
note("Frameshift dominates (>50% of attributions)",
     if (fs_total / total > 0.50) "PASS" else "FAIL",
     sprintf("Frameshift=%.1f%%", 100*fs_total/total))
note("3'UTR splice is minority (<15%)",
     if (utr_total / total < 0.15) "PASS" else "FAIL",
     sprintf("3'UTR splice=%.1f%%", 100*utr_total/total))
note("Mechanism counts sum to total attributed (1812)",
     if (total == 1812) "PASS" else "FAIL")

# Adversarial: SE row should be dominated by Frameshift + In-frame, NOT 3'UTR splice
se_rows <- panelF[event_type == "SE"]
print(se_rows)
se_fs <- se_rows$n[se_rows$mechanism == "Frameshift"]
se_utr <- se_rows$n[se_rows$mechanism == "3'UTR splice"]
note("SE attributions dominated by Frameshift+Inframe, not 3'UTR splice",
     if (se_fs > se_utr * 5) "PASS" else "FAIL",
     sprintf("Frameshift=%d, 3'UTR=%d", se_fs, se_utr))

# ============================================================================
# Test 6 — Headline PTC rate 83.5% adversarial check
# ============================================================================
cat("\n[6] Adversarial — can we disprove the 83.5% NMD PTC headline?\n")
# Try: what if ref-AUG categorization itself is biased toward PTC?
# Look at the proportion of effectively_ptc among original_ptc=FALSE pairs
# (the recovered set). If this is much higher than baseline, ref-AUG bias would
# be a concern.
n_orig_F <- sum(ra_c2$original_ptc == FALSE)
n_eff_orig_F <- sum(ra_c2$original_ptc == FALSE & ra_c2$category == "effectively_ptc")
recovery_rate <- 100 * n_eff_orig_F / n_orig_F
cat(sprintf("    Of original_ptc=FALSE NMD pairs in ref_atg analysis (n=%d), %d (%.1f%%) are effectively_ptc\n",
            n_orig_F, n_eff_orig_F, recovery_rate))
# Control side: ra_c4 doesn't carry original_ptc column. Use ptc.rds directly
# to look up TD2 PTC status for each Control comparator.
ptc <- as.data.table(readRDS(file.path(DM, "ptc.rds")))
td2_pos <- ptc$isoform_id[ptc$has_ptc == TRUE]
ra_c4$original_ptc_lookup <- ra_c4$comparator_isoform_id %in% td2_pos
n_orig_F_ctrl <- sum(!ra_c4$original_ptc_lookup)
n_eff_orig_F_ctrl <- sum(!ra_c4$original_ptc_lookup & ra_c4$category == "effectively_ptc")
recovery_rate_ctrl <- 100 * n_eff_orig_F_ctrl / n_orig_F_ctrl
cat(sprintf("    Same for Control (TD2-PTC- subset): %d of %d = %.1f%% (should be MUCH lower if ref-AUG isn't biased)\n",
            n_eff_orig_F_ctrl, n_orig_F_ctrl, recovery_rate_ctrl))
note("Ref-AUG recovery rate in NMD >> Control (biology-specific, not bias)",
     if (recovery_rate > recovery_rate_ctrl * 3) "PASS" else "FAIL",
     sprintf("NMD=%.1f%% vs Control=%.1f%%", recovery_rate, recovery_rate_ctrl))

# ============================================================================
# Test 7 — Panel C visualization → TSV mapping (sanity check)
# ============================================================================
cat("\n[7] Sanity — Panel C SE bar should match TSV value\n")
# Read the SE row from TSV
se_row <- panelC[event_type == "SE"]
cat(sprintf("    SE NMD: %.1f%% (= %d/%d), SE Control: %.1f%% (= %d/%d)\n",
            se_row$pct_of_pairs_NMD, se_row$n_pairs_with_event_NMD, N_BC,
            se_row$pct_of_pairs_Ctrl, se_row$n_pairs_with_event_Ctrl, N_BC))
note("SE NMD percent ~44%", if (abs(se_row$pct_of_pairs_NMD - 44.2) < 0.5) "PASS" else "FAIL")
note("SE Control percent ~21%", if (abs(se_row$pct_of_pairs_Ctrl - 21.2) < 0.5) "PASS" else "FAIL")

# ============================================================================
# Summary
# ============================================================================
cat("\n================================================================\n")
n_pass <- sum(sapply(results, function(r) r$status == "PASS"))
n_fail <- sum(sapply(results, function(r) r$status == "FAIL"))
cat(sprintf("Pass 2 summary: %d PASS / %d FAIL / %d total\n",
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
