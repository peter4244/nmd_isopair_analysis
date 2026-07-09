#!/usr/bin/env Rscript
# Pass 1: Factual accuracy. Recompute n / median / IQR independently from
# the long-format TSV and compare to the *_descriptives.tsv files.

suppressPackageStartupMessages({ library(data.table) })

HERE <- "/Users/petecastaldi/claude_projects/nmd/figures/SupplementalFigures/CDSand3UTR_GENCODEonly"
long <- fread(file.path(HERE, "data", "features_190_subset_long.tsv"))

cat(sprintf("Total rows: %d\n", nrow(long)))
cat("Group counts:\n"); print(long[, .N, by = group])

panels <- list(
  list(value_col = "cds_nt",                   slug = "panelA_cds_length",
       expected_n = c(`NMD+/PTC+` = 72, `NMD+/PTC-` = 118, `Control` = 190)),
  list(value_col = "utr3_to_tx_end_nt",        slug = "panelC_utr3_translation",
       expected_n = c(`NMD+/PTC+` = 72, `NMD+/PTC-` = 118, `Control` = 190)),
  list(value_col = "utr3_via_non_ptc_stop_nt", slug = "panelD_utr3_non_ptc_stop",
       expected_n = c(`NMD+/PTC+` = 66, `NMD+/PTC-` = 118, `Control` = 190))
)

cat("\n=== Recomputed descriptives vs stored TSVs ===\n")
all_ok <- TRUE
for (p in panels) {
  cat(sprintf("\n--- %s (value=%s) ---\n", p$slug, p$value_col))
  stored <- fread(file.path(HERE, "data", paste0(p$slug, "_descriptives.tsv")))
  recomp <- long[!is.na(get(p$value_col)),
                 .(n = .N,
                   median = median(get(p$value_col)),
                   q25 = quantile(get(p$value_col), 0.25),
                   q75 = quantile(get(p$value_col), 0.75)),
                 by = group][order(match(group, c("NMD+/PTC+","NMD+/PTC-","Control")))]
  for (g in c("NMD+/PTC+","NMD+/PTC-","Control")) {
    s <- stored[group == g]; r <- recomp[group == g]
    ok_n   <- s$n      == r$n
    ok_med <- isTRUE(all.equal(as.numeric(s$median), as.numeric(r$median), tolerance=1e-6))
    ok_q25 <- isTRUE(all.equal(as.numeric(s$q25),    as.numeric(r$q25),    tolerance=1e-6))
    ok_q75 <- isTRUE(all.equal(as.numeric(s$q75),    as.numeric(r$q75),    tolerance=1e-6))
    ok_expected_n <- s$n == p$expected_n[[g]]
    cat(sprintf("  %-10s n=%d (recomp %d, expected %d) median=%.2f (recomp %.2f) q25=%.2f (recomp %.2f) q75=%.2f (recomp %.2f) %s\n",
                g, s$n, r$n, p$expected_n[[g]],
                s$median, r$median, s$q25, r$q25, s$q75, r$q75,
                if (ok_n && ok_med && ok_q25 && ok_q75 && ok_expected_n) "OK" else "FAIL"))
    if (!(ok_n && ok_med && ok_q25 && ok_q75 && ok_expected_n)) all_ok <- FALSE
  }
}

cat(sprintf("\nOverall Pass 1: %s\n", if (all_ok) "PASS" else "FAIL"))
