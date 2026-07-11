# Handoff: Isopair reference-selection correctness issue (from SF26 investigation)

## TL;DR
The Isopair **reference isoform** is, for **~31% of genes (934/3,009)**, NOT the gene's
dominant transcript — it's a minor (often near-zero, often novel) isoform. Root cause is a
**classification-set mismatch**: reference eligibility requires *strict* non-NMD, but the gene's
dominant transcript frequently falls in a **"neither" gap** (not NMD-responsive, but also fails
the strict non-NMD threshold). This surfaced as an artifactual low-share spike in **SF26**, but
the reference is the **baseline for the entire pair analysis**, so it likely affects
**Figure 3/4** and any per-pair comparison, not just SF26.

## How reference selection works (confirmed)
File: `results/isoform_transitions/Version_6.0/isopair_wrapper/02_build_profiles_mashr.R`
- Lines **147–154**: C4 pairs built by
  `generatePairsExpression(expr_non_nmd, gene_map_non_nmd, dmso_samp[[ct]], method = "top_two")`
  → "top two non-NMD by CPM in DMSO", **per cell type**. Reference = the top; Control comparator = 2nd.
- `expr_non_nmd` candidate pool = isoforms in the **strict** non-NMD set only.
- C2 (line 167+) reuses the same reference, pairs it with the top NMD isoform in SMG1i.
- `generatePairsExpression` / `identifyDominantIsoforms` live in `isopair_wrapper/analysis_functions.R`
  (or the Isopair package) — NOT yet read in detail; verify the exact candidate-pool filter there.

## The two "non-NMD" definitions that conflict
- **Reference eligibility (strict non-NMD):** `nmd_classification.rds$all_samples$non_nmd`
  = isoforms with `adj.P.Val > 0.30 in ALL 4 cell types` (per NMD_orf_model_v5_4ct CLAUDE.md).
- **SF26 `gene_total` denominator:** everything **NOT** in `$all_samples$nmd` — the looser
  complement, which ALSO includes the **"neither" gap** (fails strict non-NMD in ≥1 CT but not NMD).
- Dominant "neither" isoforms are **excluded from reference candidacy** but **included in the
  denominator** → reference share collapses.

## Worked example (fully verified)
Gene `ENSG00000169288.19`, per-CT DMSO mean expression:
| isoform | mean DMSO | is_nmd | strict non_nmd | eligible ref? |
|---|---|---|---|---|
| ENST00000315567.13 | **85.8** | F | F | no — "neither" (true dominant) |
| novel2 | 2.49 | F | F | no — "neither" |
| novel3/novel1/novel4 | ≤0.68 | **T** | F | no — NMD |
| **novel7** | 0.14 | F | **T** | **← chosen reference** (highest eligible) |
| novel15 | 0.10 | F | T | 2nd eligible |
Selection is "correct" given eligibility (novel7 is the top strict-non-NMD isoform), but the
gene is dominated by two "neither" isoforms (85.8 + 2.49) that can't be references yet ARE in
the denominator → share = 0.14 / ~104 ≈ 0.1%.

## Quantified impact
- 934/3,009 genes (31%) violate `ref_fraction >= 1/n_nonNMD` (impossible if ref were the max).
- All 856 genes with <10% reference share are in the violating set.
- Median reference share: **30.8% (all) vs 54.3% (non-violators)**.
- References that are novel isoforms: **59.6% among violators vs 17.1% among non-violators**.

## Data / code map for the next session
- Reference-share export logic (SF26): recover from git — it was in the now-deleted
  `figures/SupplementalFigures/PairSetDescriptives/data_export.R`
  (`git show f6d96bc^:figures/SupplementalFigures/PairSetDescriptives/data_export.R`),
  lines ~108–140 (nmd_iso_set filter, gene_total_dmso sum, ref_fraction_of_gene).
- `results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr/`:
  `expression_data.rds` (123790×36, isoform×sample; colnames `SampleN_{CT}_{donor}_{treatment}`),
  `gene_map.rds` (isoform_id, gene_id), `nmd_classification.rds` (list: AT/DD/FB/MV/all_samples,
  each `$nmd` and `$non_nmd` character vectors).
- SF26 figure + its committed data: `figures/SupplementalFigures/SF26_ReferenceShare/`
  (`data/ref_expression_fraction.tsv`, `data/descriptives_summary.tsv`).
- Manuscript claim to reconcile: methods text says "high-expressing non-NMD reference isoform"
  (`paper/methods_updates_2026-06-18.md`); SF26 legend says "highest-expressed non-NMD isoform".

## Decisions to make (do NOT change anything without Pete)
1. **Design/scientific:** should "neither"-gap dominant isoforms be eligible references? The
   reference is meant to be the gene's normal baseline; using a 0.14-expression novel isoform as
   the baseline when an 85-expression transcript exists may distort pair comparisons.
2. **Scope:** quantify how many pair-analysis genes (Figure 3/4, SF27) have a "neither"-gap
   dominant isoform displacing the reference — size the exposure before deciding a fix.
3. **SF26 presentation (safe, pending #1):** make numerator/denominator use the same set
   (strict non-NMD only) → median jumps to ~54%, legend becomes accurate. But this only masks
   the deeper question if the reference itself is wrong for the pair analysis.

## Status of everything else (DONE this session, pushed to both remotes @ 6ea9971)
All Supplemental Figure corrections complete: 13-item list + follow-ons (ORF abbreviation sweep,
nt definitions, NMD-abbrev removal, reference-AUG hyphens, SF29 Gain, SF40 threshold line,
SF37 native re-render replacing the stale placeholder, PairSetDescriptives deleted, SF30/SF35
ref-AUG fixes). Final QA-verified docx: `paper/nmd_supplemental_figures_sf24_sf42.docx`.
Deferred (manuscript-side, Pete owns): update methods_updates PairSetDescriptives reference;
"Main ORF selection" Methods cross-ref.
