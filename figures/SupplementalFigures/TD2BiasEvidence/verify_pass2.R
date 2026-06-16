#!/usr/bin/env Rscript
# Pass 2: Result correctness.
#   - Paired Wilcoxon p-values for Panels B / E (Kozak)
#   - Directional counts: ref Kozak > TD2 Kozak in 446/1166 (38.3%) / 384/492 (78.0%)
#   - TD2-downstream-of-ref-AUG fractions (Panel C / F)

suppressPackageStartupMessages({ library(data.table) })

HERE <- "/Users/petecastaldi/claude_projects/nmd/figures/SupplementalFigures/TD2BiasEvidence"

# Panel B — Kozak paired Wilcoxon (broad scope, n=1166)
B <- fread(file.path(HERE, "data", "panelB_kozak.tsv"))
ref <- B[ATG_source == "Reference ATG",     .(comparator_isoform_id, score)]
td2 <- B[ATG_source == "TD2-predicted ATG", .(comparator_isoform_id, score)]
paired_B <- merge(ref, td2, by = "comparator_isoform_id",
                  suffixes = c("_ref", "_td2"))
nB <- nrow(paired_B)
n_ref_higher_B <- sum(paired_B$score_ref > paired_B$score_td2)
n_td2_higher_B <- sum(paired_B$score_ref < paired_B$score_td2)
n_equal_B <- sum(paired_B$score_ref == paired_B$score_td2)
wB <- wilcox.test(paired_B$score_ref, paired_B$score_td2, paired = TRUE)

# Panel E — Kozak paired Wilcoxon (occult-PTC, n=492)
E <- fread(file.path(HERE, "data", "panelE_kozak_occult.tsv"))
ref <- E[ATG_source == "Reference ATG",     .(comparator_isoform_id, score)]
td2 <- E[ATG_source == "TD2-predicted ATG", .(comparator_isoform_id, score)]
paired_E <- merge(ref, td2, by = "comparator_isoform_id",
                  suffixes = c("_ref", "_td2"))
nE <- nrow(paired_E)
n_ref_higher_E <- sum(paired_E$score_ref > paired_E$score_td2)
n_td2_higher_E <- sum(paired_E$score_ref < paired_E$score_td2)
n_equal_E <- sum(paired_E$score_ref == paired_E$score_td2)
wE <- wilcox.test(paired_E$score_ref, paired_E$score_td2, paired = TRUE)

cat("=== Panel B (broad, n=1166) paired Kozak Wilcoxon ===\n")
cat(sprintf("  n_pairs=%d  p=%.4e  (legend 1.2e-38)\n", nB, wB$p.value))
cat(sprintf("  ref > td2: %d / %d (%.1f%%)  (legend 446 / 1166 / 38.3%%)\n",
            n_ref_higher_B, nB, 100*n_ref_higher_B/nB))
cat(sprintf("  ref < td2: %d  ref == td2: %d  (50%% equal means TD2 same start)\n",
            n_td2_higher_B, n_equal_B))

cat("\n=== Panel E (occult, n=492) paired Kozak Wilcoxon ===\n")
cat(sprintf("  n_pairs=%d  p=%.4e  (legend 2.6e-38)\n", nE, wE$p.value))
cat(sprintf("  ref > td2: %d / %d (%.1f%%)  (legend 384 / 492 / 78.0%%)\n",
            n_ref_higher_E, nE, 100*n_ref_higher_E/nE))
cat(sprintf("  ref < td2: %d  ref == td2: %d\n", n_td2_higher_E, n_equal_E))

# Panels C / F — TD2 downstream fractions
C <- fread(file.path(HERE, "data", "panelC_td2_position.tsv"))
F_ <- fread(file.path(HERE, "data", "panelF_td2_position_occult.tsv"))
cat("\n=== Panel C — TD2 vs ref AUG position ===\n")
cat(sprintf("  n=%d  eq=%d (%.1f%%)  down=%d (%.1f%%)  up=%d (%.1f%%)\n",
            nrow(C), sum(C$distance_nt == 0), 100*mean(C$distance_nt == 0),
            sum(C$distance_nt > 0),  100*mean(C$distance_nt > 0),
            sum(C$distance_nt < 0),  100*mean(C$distance_nt < 0)))
cat(sprintf("  median downstream offset: %.0f nt (README +471)\n",
            median(C[distance_nt > 0, distance_nt])))

cat("\n=== Panel F — TD2 vs ref AUG position (occult) ===\n")
cat(sprintf("  n=%d  eq=%d (%.1f%%)  down=%d (%.1f%%)  up=%d (%.1f%%)\n",
            nrow(F_), sum(F_$distance_nt == 0), 100*mean(F_$distance_nt == 0),
            sum(F_$distance_nt > 0),  100*mean(F_$distance_nt > 0),
            sum(F_$distance_nt < 0),  100*mean(F_$distance_nt < 0)))
cat(sprintf("  median downstream offset: %.0f nt (README +488 / summary +476)\n",
            median(F_[distance_nt > 0, distance_nt])))

# Headline length ratios
A <- fread(file.path(HERE, "data", "panelA_td2_vs_refaug_length.tsv"))
D <- fread(file.path(HERE, "data", "panelD_td2_vs_refaug_length_occult.tsv"))
ratio_A <- median(A[ORF_source == "TD2-selected CDS",  length_nt]) /
           median(A[ORF_source == "Reference ATG ORF", length_nt])
ratio_D <- median(D[ORF_source == "TD2-selected CDS",  length_nt]) /
           median(D[ORF_source == "Reference ATG ORF", length_nt])
cat(sprintf("\nPanel A length ratio: %.2fx (legend 4.69x)\n", ratio_A))
cat(sprintf("Panel D length ratio: %.2fx (legend 7.73x)\n", ratio_D))

# Headline % checks
cat("\n=== Headline number cross-check ===\n")
cat(sprintf("  Panel A TD2 same as ref AUG: %d/%d=%.0f%% (README 50%%)\n",
            sum(C$distance_nt == 0), nrow(C), 100*mean(C$distance_nt == 0)))
cat(sprintf("  Panel A TD2 != ref AUG: %d/%d=%.0f%%   (README 50%%)\n",
            sum(C$distance_nt != 0), nrow(C), 100*mean(C$distance_nt != 0)))
