#!/usr/bin/env Rscript
# Pass 1: Factual accuracy. Recompute n / median / IQR for each panel TSV
# independently and check against headline numbers in README, legend, and
# summary.tsv.

suppressPackageStartupMessages({ library(data.table) })

HERE <- "/Users/petecastaldi/claude_projects/nmd/figures/SupplementalFigures/TD2BiasEvidence"

cat("=== Panel A (broad, length) ===\n")
A <- fread(file.path(HERE, "data", "panelA_td2_vs_refaug_length.tsv"))
print(A[, .(n = .N, median_nt = median(length_nt),
            q25 = quantile(length_nt, 0.25),
            q75 = quantile(length_nt, 0.75)),
        by = ORF_source])
n_A_ref <- nrow(A[ORF_source == "Reference ATG ORF"])
n_A_td2 <- nrow(A[ORF_source == "TD2-selected CDS"])
cat(sprintf("  n_ref=%d  n_td2=%d  (expected 1166 each)\n", n_A_ref, n_A_td2))

cat("\n=== Panel D (occult-PTC, length) ===\n")
D <- fread(file.path(HERE, "data", "panelD_td2_vs_refaug_length_occult.tsv"))
print(D[, .(n = .N, median_nt = median(length_nt),
            q25 = quantile(length_nt, 0.25),
            q75 = quantile(length_nt, 0.75)),
        by = ORF_source])
n_D_ref <- nrow(D[ORF_source == "Reference ATG ORF"])
n_D_td2 <- nrow(D[ORF_source == "TD2-selected CDS"])
cat(sprintf("  n_ref=%d  n_td2=%d  (expected 492 each)\n", n_D_ref, n_D_td2))

cat("\n=== Panel B (broad, Kozak) ===\n")
B <- fread(file.path(HERE, "data", "panelB_kozak.tsv"))
print(B[, .(n = .N, median = median(score),
            q25 = quantile(score, 0.25),
            q75 = quantile(score, 0.75)),
        by = ATG_source])

cat("\n=== Panel E (occult-PTC, Kozak) ===\n")
E <- fread(file.path(HERE, "data", "panelE_kozak_occult.tsv"))
print(E[, .(n = .N, median = median(score),
            q25 = quantile(score, 0.25),
            q75 = quantile(score, 0.75)),
        by = ATG_source])

cat("\n=== Panel C (broad, position) ===\n")
C <- fread(file.path(HERE, "data", "panelC_td2_position.tsv"))
nC <- nrow(C)
nC_eq   <- nrow(C[distance_nt == 0])
nC_down <- nrow(C[distance_nt > 0])
nC_up   <- nrow(C[distance_nt < 0])
med_down <- median(C[distance_nt > 0, distance_nt])
cat(sprintf("  n=%d eq=%d (%.1f%%) down=%d (%.1f%%) up=%d (%.1f%%) median_downstream=%.0f nt\n",
            nC, nC_eq, 100*nC_eq/nC, nC_down, 100*nC_down/nC,
            nC_up, 100*nC_up/nC, med_down))

cat("\n=== Panel F (occult-PTC, position) ===\n")
F_ <- fread(file.path(HERE, "data", "panelF_td2_position_occult.tsv"))
nF <- nrow(F_)
nF_eq   <- nrow(F_[distance_nt == 0])
nF_down <- nrow(F_[distance_nt > 0])
nF_up   <- nrow(F_[distance_nt < 0])
med_down_F <- median(F_[distance_nt > 0, distance_nt])
cat(sprintf("  n=%d eq=%d (%.1f%%) down=%d (%.1f%%) up=%d (%.1f%%) median_downstream=%.0f nt\n",
            nF, nF_eq, 100*nF_eq/nF, nF_down, 100*nF_down/nF,
            nF_up, 100*nF_up/nF, med_down_F))

cat("\n=== Headline-number comparison vs README/legend ===\n")
checks <- list(
  list("Panel A TD2 length median", median(A[ORF_source == "TD2-selected CDS",  length_nt]), 2758),
  list("Panel A ref length median", median(A[ORF_source == "Reference ATG ORF", length_nt]), 588),
  list("Panel D TD2 length median", median(D[ORF_source == "TD2-selected CDS",  length_nt]), 2934),
  list("Panel D ref length median", median(D[ORF_source == "Reference ATG ORF", length_nt]), 380),
  list("Panel B ref Kozak median",  median(B[ATG_source == "Reference ATG",     score]),    0.86),
  list("Panel B TD2 Kozak median",  median(B[ATG_source == "TD2-predicted ATG", score]),    0.29),
  list("Panel E ref Kozak median",  median(E[ATG_source == "Reference ATG",     score]),    1.04),
  list("Panel E TD2 Kozak median",  median(E[ATG_source == "TD2-predicted ATG", score]),   -0.35),
  list("Panel C n_eq",         nC_eq,    578),
  list("Panel C n_down",       nC_down,  496),
  list("Panel C n_up",         nC_up,     92),
  list("Panel C median down",  round(med_down), 471),
  list("Panel F n_down",       nF_down,  487),
  list("Panel F n_up",         nF_up,      5),
  list("Panel F median down (legend says 488)", round(med_down_F), 488)
)
for (chk in checks) {
  actual <- chk[[2]]; expected <- chk[[3]]
  ok <- isTRUE(all.equal(actual, expected, tolerance = 0.02))
  cat(sprintf("  %-45s actual=%.3f  expected=%.3f  %s\n",
              chk[[1]], actual, expected, if (ok) "OK" else "FAIL"))
}
