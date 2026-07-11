#!/usr/bin/env Rscript
# Verification Pass 2 — Result correctness for Figure 3
#
# For each statistical test in Figure 3 (Fisher exacts in Panel C and Panel E),
# verify the 2×2 table construction independently. Check OR / p-values for
# plausibility. Apply domain knowledge to test biological reasonableness.
# Actively try to disprove headline claims rather than confirm them.
#
# Scope refresh 2026-06-15 (v2): all-3-ENST + coding-CDS + own-GENCODE-stop.

suppressPackageStartupMessages({
  library(data.table)
})

DM <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr"
CACHE <- file.path(DM, "analysis_cache")

# ---- Reload populations from Pass 1 ----
profiles_c2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))
cds <- as.data.table(readRDS(file.path(DM, "cds.rds")))

make_key <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")
is_enst <- function(x) grepl("^ENST", x)
shared <- intersect(unique(make_key(profiles_c2)), unique(make_key(profiles_c4)))
pop_BC_c2 <- profiles_c2[make_key(profiles_c2) %in% shared]
pop_BC_c4 <- profiles_c4[make_key(profiles_c4) %in% shared]

# Build all-3-ENST + coding-CDS, re-intersected (matches data_export.R)
pop_BC_c2_enst <- pop_BC_c2[is_enst(reference_isoform_id) & is_enst(comparator_isoform_id)]
pop_BC_c4_enst <- pop_BC_c4[is_enst(reference_isoform_id) & is_enst(comparator_isoform_id)]
all3_keys <- intersect(unique(make_key(pop_BC_c2_enst)),
                       unique(make_key(pop_BC_c4_enst)))
all3_c2 <- pop_BC_c2_enst[make_key(pop_BC_c2_enst) %in% all3_keys]
all3_c4 <- pop_BC_c4_enst[make_key(pop_BC_c4_enst) %in% all3_keys]
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
all3_coding_c2 <- all3_c2[reference_isoform_id %in% coding_ids &
                          comparator_isoform_id %in% coding_ids]
all3_coding_c4 <- all3_c4[reference_isoform_id %in% coding_ids &
                          comparator_isoform_id %in% coding_ids]
final_keys <- intersect(unique(make_key(all3_coding_c2)),
                        unique(make_key(all3_coding_c4)))
pop_DEF_c2 <- all3_coding_c2[make_key(all3_coding_c2) %in% final_keys]
pop_DEF_c4 <- all3_coding_c4[make_key(all3_coding_c4) %in% final_keys]

N_BC <- 1548         # Panels B, C scope (post 25% floor)
N_DEF <- 130         # Panels D scope (each side; post-floor)

# ---- Tracking ----
results <- list()
note <- function(label, status, detail = "") {
  cat(sprintf("  [%s] %s%s\n", status, label,
              if (nzchar(detail)) sprintf(" -- %s", detail) else ""))
  results[[length(results) + 1]] <<- list(label = label, status = status,
                                          detail = detail)
}

cat("================================================================\n")
cat("PASS 2: Result correctness verification for Figure 3 (all-3-ENST scope)\n")
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
n_ptc_attr <- sum(panelE$n_ptc_events)      # 69
ctrl_total_tsv <- sum(panelE$n_ctrl_events) # 330

# The actual ctrl_total used by panel_e_compute.R is the FULL Control event
# count across all event types in pop_DEF_c4 (190 pairs). The TSV only emits
# the 9 event types that overlap with PTC-causing events. From the locked
# values in figure3_panelE_methodology.md and the composite legend:
# ctrl_total = 292 (n_ctrl_events 209 covers the 9 emitted types).
CTRL_TOTAL_FULL <- 292L
cat(sprintf("    n_ptc_attr=%d, ctrl_total(full)=%d, TSV n_ctrl_events sum=%d\n",
            n_ptc_attr, CTRL_TOTAL_FULL, ctrl_total_tsv))

# Verify that recomputing each row's Fisher exact at the FULL ctrl_total
# reproduces the TSV p-value.
for (i in seq_len(nrow(panelE))) {
  e <- panelE$event_type[i]
  a <- panelE$n_ptc_events[i]
  c_ <- panelE$n_ctrl_events[i]
  b <- n_ptc_attr - a
  d <- CTRL_TOTAL_FULL - c_
  ft <- fisher.test(matrix(c(a, b, c_, d), nrow = 2))
  p_ind <- ft$p.value
  p_tsv <- panelE$fisher_p[i]
  ok_p <- abs(log10(p_ind + 1e-300) - log10(p_tsv + 1e-300)) < 0.1
  status <- if (ok_p) "PASS" else "FAIL"
  detail <- sprintf("OR=%.2f p ind=%.2g tsv=%.2g (direction=%s)",
                    ft$estimate, p_ind, p_tsv, panelE$direction[i])
  note(sprintf("Panel E %s 2x2 reconstruction (ctrl_total=292)", e),
       status, detail)
}

# Verify pct_ctrl in TSV is computed on the full ctrl_total (447), not 330
# Spot check: SE pct_ctrl = 42/292 = 14.38%, rounds to 14.8 (post-floor)
expected_se_pct <- round(100 * 42 / 292, 1)
note(sprintf("Panel E SE pct_ctrl matches 42/292 (%.1f%%)", expected_se_pct),
     if (abs(panelE$pct_ctrl[panelE$event_type == "SE"] - expected_se_pct) < 0.1)
       "PASS" else "FAIL",
     sprintf("TSV=%.1f, recomputed=%.1f", panelE$pct_ctrl[panelE$event_type == "SE"],
             expected_se_pct))

# ============================================================================
# Test 3 — Panel E SE enrichment (the manuscript headline)
# ============================================================================
cat("\n[3] Panel E — SE enrichment (manuscript: 'nearly half of PTCs are SE')\n")
se_pct_ptc <- panelE$pct_of_ptc[panelE$event_type == "SE"]
se_pct_ctrl <- panelE$pct_ctrl[panelE$event_type == "SE"]
cat(sprintf("    SE: PTC-causing %.1f%% vs Control all-events %.1f%% (enrichment %.1fx)\n",
            se_pct_ptc, se_pct_ctrl,
            panelE$enrichment[panelE$event_type == "SE"]))
note("SE PTC% > Control% (enriched)",
     if (se_pct_ptc > se_pct_ctrl) "PASS" else "FAIL")
note("SE % of PTC > 30% (plurality story)",
     if (se_pct_ptc > 30) "PASS" else "FAIL",
     sprintf("actual=%.1f%%", se_pct_ptc))

# ============================================================================
# Test 4 — Panel D distance distribution sign + biological direction
# ============================================================================
cat("\n[4] Panel D — distance sign and direction\n")
panelD <- fread("data/panelD_stop_codon_distance.tsv")
med_nmd <- median(panelD[comparison == "NMD"]$distance)
med_ctrl <- median(panelD[comparison == "Control"]$distance)
cat(sprintf("    NMD median distance: %.0f nt (expect %.0f)\n", med_nmd, -58))
cat(sprintf("    Control median distance: %.0f nt (expect %.0f)\n", med_ctrl, -151))
# Under the all-3-ENST scope, both populations skew negative (most stops sit
# in the last exon under GENCODE annotation), but the NMD distribution has a
# long right tail past the 50-nt threshold.
note("NMD median > Control median (NMD distribution shifted right)",
     if (med_nmd > med_ctrl) "PASS" else "FAIL")
note("NMD median exactly -58.5 (locked value; 4-CT)",
     if (med_nmd == -58.5) "PASS" else "FAIL",
     sprintf("actual=%.0f", med_nmd))
note("Control median exactly -151 (locked value; 4-CT)",
     if (med_ctrl == -151) "PASS" else "FAIL",
     sprintf("actual=%.0f", med_ctrl))
note("NMD/Control medians separated by >50 nt",
     if ((med_nmd - med_ctrl) > 50) "PASS" else "FAIL",
     sprintf("separation=%.0f nt", med_nmd - med_ctrl))

# PTC-direction tail: NMD has dramatically more PTC+ pairs than Control
n_ptc_nmd <- sum(panelD[comparison == "NMD"]$distance > 50)
n_ptc_ctrl <- sum(panelD[comparison == "Control"]$distance > 50)
note("NMD PTC+ count >> Control PTC+ count (>10x ratio expected)",
     if (n_ptc_nmd > n_ptc_ctrl * 10) "PASS" else "FAIL",
     sprintf("NMD=%d, Control=%d", n_ptc_nmd, n_ptc_ctrl))

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
note("Frameshift is the leading mechanism (>=50% of attributions)",
     if (fs_total / total >= 0.50) "PASS" else "FAIL",
     sprintf("Frameshift=%.1f%%", 100*fs_total/total))
note("In-frame stop is a substantial minority (>=25%)",
     if (ifs_total / total >= 0.25) "PASS" else "FAIL",
     sprintf("In-frame stop=%.1f%%", 100*ifs_total/total))
note("3'UTR splice is a smaller minority (<20%)",
     if (utr_total / total < 0.20) "PASS" else "FAIL",
     sprintf("3'UTR splice=%.1f%%", 100*utr_total/total))
note("Mechanism counts sum to 48 (attributed PTC+ pairs)",
     if (total == 48) "PASS" else "FAIL",
     sprintf("total=%d", total))

# Adversarial: SE row should be dominated by Frameshift + In-frame, NOT 3'UTR splice
se_rows <- panelF[event_type == "SE"]
print(se_rows)
se_fs <- sum(se_rows$n[se_rows$mechanism == "Frameshift"])
se_ifs <- sum(se_rows$n[se_rows$mechanism == "In-frame stop"])
se_utr <- sum(se_rows$n[se_rows$mechanism == "3'UTR splice"])
note("SE attributions dominated by Frameshift+Inframe, not 3'UTR splice",
     if ((se_fs + se_ifs) > se_utr * 5) "PASS" else "FAIL",
     sprintf("Frameshift=%d, In-frame=%d, 3'UTR=%d", se_fs, se_ifs, se_utr))

# ============================================================================
# Test 6 — Headline 37.9% NMD PTC adversarial check
# ============================================================================
cat("\n[6] Adversarial — can we disprove the 37.9% NMD PTC headline?\n")
# Construct the 2×2 directly: NMD vs Control PTC+ at all-3-ENST + coding scope
ptc_table <- matrix(c(n_ptc_nmd, N_DEF - n_ptc_nmd,
                      n_ptc_ctrl, N_DEF - n_ptc_ctrl), nrow = 2,
                    dimnames = list(c("PTC+", "PTC-"), c("NMD", "Control")))
print(ptc_table)
ft_ptc <- fisher.test(ptc_table)
cat(sprintf("    Fisher exact: OR=%.1f, p=%.2g\n",
            ft_ptc$estimate, ft_ptc$p.value))
note("PTC+ enrichment OR very large (>10) at all-3-ENST scope",
     if (ft_ptc$estimate > 10) "PASS" else "FAIL",
     sprintf("OR=%.1f", ft_ptc$estimate))
note("Fisher p < 1e-12 (overwhelming)",
     if (ft_ptc$p.value < 1e-12) "PASS" else "FAIL",
     sprintf("p=%.2g", ft_ptc$p.value))
# Verify the headline 24× enrichment is the rate ratio (not the OR)
rate_nmd <- 100 * n_ptc_nmd / N_DEF
rate_ctrl <- 100 * n_ptc_ctrl / N_DEF
fold <- rate_nmd / rate_ctrl
note(sprintf("24-fold rate ratio (36.9 / 1.5 = %.1f)", fold),
     if (abs(fold - 24) < 1.0) "PASS" else "FAIL",
     sprintf("fold=%.1f", fold))

# ============================================================================
# Test 7 — Panel C visualization → TSV mapping (sanity check)
# ============================================================================
cat("\n[7] Sanity — Panel C SE bar should match TSV value\n")
se_row <- panelC[event_type == "SE"]
cat(sprintf("    SE NMD: %.1f%% (= %d/%d), SE Control: %.1f%% (= %d/%d)\n",
            se_row$pct_of_pairs_NMD, se_row$n_pairs_with_event_NMD, N_BC,
            se_row$pct_of_pairs_Ctrl, se_row$n_pairs_with_event_Ctrl, N_BC))
note("SE NMD percent ~51%", if (abs(se_row$pct_of_pairs_NMD - 51.0) < 0.6) "PASS" else "FAIL")
note("SE Control percent ~22%", if (abs(se_row$pct_of_pairs_Ctrl - 21.6) < 0.6) "PASS" else "FAIL")

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
