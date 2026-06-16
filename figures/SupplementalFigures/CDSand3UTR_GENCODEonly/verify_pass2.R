#!/usr/bin/env Rscript
# Pass 2: Result correctness. Recompute pairwise Wilcoxon HL shift / p-value
# from the long TSV and compare to stored *_pairwise.tsv.

suppressPackageStartupMessages({ library(data.table) })

HERE <- "/Users/petecastaldi/claude_projects/nmd/figures/SupplementalFigures/CDSand3UTR_GENCODEonly"
long <- fread(file.path(HERE, "data", "features_190_subset_long.tsv"))

panels <- list(
  list(value_col = "cds_nt",                   slug = "panelA_cds_length"),
  list(value_col = "utr3_to_tx_end_nt",        slug = "panelC_utr3_translation"),
  list(value_col = "utr3_via_non_ptc_stop_nt", slug = "panelD_utr3_non_ptc_stop")
)
pairs <- list(c("NMD+/PTC-", "NMD+/PTC+"),
              c("NMD+/PTC-", "Control"),
              c("NMD+/PTC+", "Control"))

all_ok <- TRUE
for (p in panels) {
  cat(sprintf("\n--- %s ---\n", p$slug))
  stored <- fread(file.path(HERE, "data", paste0(p$slug, "_pairwise.tsv")))
  for (pr in pairs) {
    a <- long[group == pr[1], get(p$value_col)]; a <- a[!is.na(a)]
    b <- long[group == pr[2], get(p$value_col)]; b <- b[!is.na(b)]
    w <- suppressWarnings(wilcox.test(a, b, conf.int = TRUE))
    hl_r <- as.numeric(w$estimate)
    pv_r <- w$p.value
    row <- stored[group_x == pr[1] & group_y == pr[2]]
    hl_s <- as.numeric(row$hl_shift); pv_s <- as.numeric(row$wilcox_p)
    ok_hl <- isTRUE(all.equal(hl_s, hl_r, tolerance = 1e-4))
    ok_pv <- isTRUE(all.equal(pv_s, pv_r, tolerance = 1e-6))
    cat(sprintf("  %-25s HL=%.3f (recomp %.3f) p=%.4g (recomp %.4g) %s\n",
                paste(pr, collapse = " vs "), hl_s, hl_r, pv_s, pv_r,
                if (ok_hl && ok_pv) "OK" else "FAIL"))
    if (!(ok_hl && ok_pv)) all_ok <- FALSE
  }
}

# Check that legend-stated p-values + HL shifts match.
cat("\n=== Legend prose number checks ===\n")
# Legend Panel A: medians 760 vs 834 vs 968; PTC+ vs Control HL shift = -186 nt, p = 5e-3
A <- fread(file.path(HERE, "data", "panelA_cds_length_pairwise.tsv"))
r <- A[group_x == "NMD+/PTC+" & group_y == "Control"]
cat(sprintf("  Panel A PTC+ vs Control: HL=%.0f (legend -186), p=%.2e (legend 5e-3) %s\n",
            r$hl_shift, r$wilcox_p,
            if (round(r$hl_shift) == -186 && abs(r$wilcox_p - 0.005) < 0.001) "OK" else "FAIL"))
# Legend Panel B: medians 869 vs 517; p = 1e-4
C <- fread(file.path(HERE, "data", "panelC_utr3_translation_pairwise.tsv"))
r <- C[group_x == "NMD+/PTC+" & group_y == "Control"]
cat(sprintf("  Panel B PTC+ vs Control: medians %.0f vs %.0f, p=%.2e (legend 1e-4) %s\n",
            r$median_x, r$median_y, r$wilcox_p,
            if (abs(r$wilcox_p - 1e-4) < 1e-4) "OK (~1e-4 OK)" else "CHECK"))
# Legend Panel C: medians 524 vs 517; p = 0.90; PTC- vs Control p = 0.03
D <- fread(file.path(HERE, "data", "panelD_utr3_non_ptc_stop_pairwise.tsv"))
r1 <- D[group_x == "NMD+/PTC+" & group_y == "Control"]
r2 <- D[group_x == "NMD+/PTC-" & group_y == "Control"]
cat(sprintf("  Panel C PTC+ vs Control: medians %.0f vs %.0f, p=%.2f (legend 0.90) %s\n",
            r1$median_x, r1$median_y, r1$wilcox_p,
            if (abs(r1$wilcox_p - 0.90) < 0.02) "OK" else "FAIL"))
cat(sprintf("  Panel C PTC- vs Control: medians %.0f vs %.0f, p=%.2f (legend 0.03) %s\n",
            r2$median_x, r2$median_y, r2$wilcox_p,
            if (abs(r2$wilcox_p - 0.03) < 0.01) "OK" else "FAIL"))

cat(sprintf("\nOverall Pass 2 (recomp): %s\n", if (all_ok) "PASS" else "FAIL"))
