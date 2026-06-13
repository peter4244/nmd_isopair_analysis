# Figure 3 Panel B — methodology

**Title:** Sequence shared with reference
**Render script:** `figure3_panelB_sequence_similarity.py`
**Data TSV:** `data/panelB_sequence_similarity.tsv`
**Data export:** `data_export.R` ([Panel B] section)
**Scope:** All data (no test-set filter)

## Headline claim

NMD susceptible comparator isoforms share more exonic sequence with the reference (dominant non-NMD isoform of the same gene) than Control comparator isoforms do. The NMD-pair `pct_shared` distribution peaks near 85–90%; the Control-pair distribution is broader and centered lower. Consistent with NMD-causing splicing being a targeted (single-event) perturbation of an otherwise-similar isoform.

## Source data

| File | What it provides |
|---|---|
| `data_mashr/analysis_cache/div_c2_allsamples.rds` | Per-pair exonic sequence similarity (`pct_shared`) for NMD pairs, computed by Isopair `quantifyPairDivergence()`. Row-aligned to `profiles_c2_allsamples.rds`. |
| `data_mashr/analysis_cache/div_c4_allsamples.rds` | Same for Control pairs. |
| `data_mashr/profiles_c{2,4}_allsamples.rds` | Pair-level metadata used to apply Stage 2 gene-matching. |

## Population

**Denominator: 3,009 pairs per side** (Stage 2 gene-matched, `pop_BC` in `data_export.R`).

Sequence-similarity analysis does not require CDS information — `pct_shared` is computed at the exonic level. So this panel uses the broadest meaningful population (Stage 2 gene-matched), not the ref-AUG-traceable subset that Panels D/E/F use.

## Computation

Per-pair `pct_shared` already computed by Isopair's `quantifyPairDivergence()`. The TSV joins the Stage 2 gene-matched index to the divergence cache row positions (both row-aligned). Density plot uses `scipy.stats.gaussian_kde` with Scott's bandwidth on each comparison's pct_shared values.

## Output TSV columns

| Column | Type | Notes |
|---|---|---|
| `comparison` | "NMD" \| "Control" | |
| `pct_shared` | float | Percent of exonic sequence shared between comparator and reference (0–100) |

## Caveats / limitations

1. **Different denominator from D/E/F.** Panel B uses 3,009 (Stage 2 gene-matched); Panels D/E/F use the ref-AUG-traceable subset (2,289 NMD / 1,763 Control or smaller). This deviation is justified — sequence similarity is computed at the exonic level and does not depend on CDS identification, so restricting to ref-AUG-traceable pairs would arbitrarily exclude pairs that have meaningful sequence-similarity data.
2. `gaussian_kde` smooths the empirical distribution. Sharp boundaries (0% and 100%) appear as upturns near the edges; this is a KDE artifact, not biology.
3. Both populations show some right-mode density near 100%, representing pairs that differ from the reference by only TSS / TES choice (no exonic sequence change).
4. Test-only sensitivity check not generated; primary analysis uses all data per project policy.

## Cross-references

- Rmd source: `05_final_report_mashr.Rmd` chunk `goal1-fig1-pct-shared` (line ~1302)
- Isopair package: `R/compare-sets.R` (the `quantifyPairDivergence()` function and `pct_shared` definition)
- `feedback_nmd_analysis_scope_test_vs_all` (all data per project default)
- `feedback_figure_sample_size_consistency` (denominator difference vs Panel D/E/F is justified by analytical need; documented above)

## Composite cross-reference

This panel appears in the Figure 3 composite as the labeled cell — see `figure3_composite_methodology.md` for the layout and per-cell mapping. The composite embeds the panel's pre-rendered PNG (`figure3_panelB.png`); regenerating the composite does NOT re-run this panel. Re-render this panel script first if its data or rendering has changed.
