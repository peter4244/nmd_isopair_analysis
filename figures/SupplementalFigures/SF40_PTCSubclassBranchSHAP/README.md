# Supplemental Figure — Branch SHAP × PTC subclass at full n = 819

Companion to the §5 manuscript sentence:

> "We also examined the Shapley scores for the START and STOP windows
> stratified by isoform PTC subclass, and we observed that the START
> window information is roughly three times more important in the
> NMD+/PTC− than in the NMD+/PTC+ isoforms (SFx)."

## What it shows

Branch-level KernelSHAP at the fusion layer of the deep-learning model
(same source as **Figure 5 Panel C**), stratified by PTC subclass
within the full n = 819 ref-AUG-traceable subset (same scope as Fig 4
C/D and Fig 5 Panel G).

**Style is consistent with Figure 5 Panel C:** branch order ranked by
importance (Structural / Stop / AUG), same colour palette
(`#d95f02` / `#FF6B6B` / `#4ECDC4`), percent above each bar, mean
|SHAP| inside each bar, no in-panel title. Three side-by-side mini-
panels — one per subgroup — share a common y-axis so heights compare
directly.

## Scope: full n = 819, not test-set-restricted

Branch attribution uses the **full-cohort** KernelSHAP table
(`kernel_shap_branch_atg500_stop500_all.tsv`, 39,938 isoforms) — not
the test-set-only file. After intersection with the n = 819 isoform
list:

| Subgroup | n |
|---|---|
| NMD+/PTC+ retained | 735 |
| NMD+/PTC− retained | 54 |
| Control | 781 |
| **Total** | **1,570 / 1,638 (95.8%)** |

The 54 NMD+/PTC− and 735 NMD+/PTC+ counts match the n quoted in
**Figure 5 Panel G**'s legend at the same scope. The ~4% drop comes
from a handful of n = 819 comparator isoforms that fall outside the
KernelSHAP table (mostly low-expression / dropped during model
training).

The sibling **SF41_PTCSubclassPerformance/** SF is still test-set-only
(AUC/AUPRC are the only test-only analyses per
[feedback_nmd_analysis_scope_test_vs_all](~/.claude/projects/-Users-petecastaldi/memory/feedback_nmd_analysis_scope_test_vs_all.md)).

## Headline numbers (full cohort)

Within-subgroup branch share (sum of subgroup-mean |SHAP|):

| Subgroup | n | Structural | Stop | AUG |
|---|---:|---:|---:|---:|
| NMD+/PTC+ | 1,016 | 61.9% (2.28) | 29.1% (1.07) | 9.0% (0.33) |
| NMD+/PTC− retained | 95 | 52.6% (0.99) | 29.2% (0.55) | **18.1% (0.34)** |
| Control | 1,107 | 62.3% (1.14) | 27.1% (0.50) | 10.6% (0.19) |

**ATG share ratio PTC− vs PTC+ = 2.01×** (within-subgroup pooled means)
**ATG share ratio PTC− vs PTC+ = 2.21×** (mean of per-isoform shares)

Both formulations are around 2×; neither is 3×.

## Discrepancy with the manuscript text — find/replace pair

The manuscript says "roughly three times more important". The actual
share-of-total ratio is **~2×**.

Apply this find/replace in the manuscript Google Doc (per repo
convention — manuscript prose changes go to the Doc, not the .md
export):

```
find:    is roughly three times more important in the NMD+/PTC-
replace: is roughly two-fold higher in the NMD+/PTC-
```

Surrounding context (no replacement needed; included for unambiguous
targeting):

> "We also examined the Shapley scores for the START and STOP windows
> stratified by isoform PTC subclass, and we observed that the START
> window information **is roughly three times more important in the
> NMD+/PTC-** than in the NMD+/PTC+ isoforms (SFx)."

## Files

| File | What |
|---|---|
| `data_export.R` | Joins n=819 list with **full-cohort** branch SHAP and (separately) test-set predictions; emits four TSVs in `data/`. The predictions TSV is also dropped into `../SF41_PTCSubclassPerformance/data/` for SF41 reuse |
| `figure_s_branch_shap_by_subclass.py` | Three Panel-C-style mini-panels (validator-clean) |
| `figure_s_branch_shap_by_subclass_legend.md` | Manuscript-style legend |
| `rmd_patch_section_8_5_refaug_subclass.Rmd` | Draft new §8.5 to splice into `NMD_orf_model_v5_4ct/orf_model_report_v5.Rmd` when cluster is back |
| `data/branch_shap_by_subclass_refaug.tsv` | Long-form per-isoform × branch (full cohort) |
| `data/branch_shap_by_subclass_refaug_descriptives.tsv` | Per-group × per-branch means |
| `data/branch_shap_by_subclass_refaug_pairwise.tsv` | Wilcoxon contrasts |
| `data/predprob_by_subclass_refaug.tsv` | Per-isoform predicted prob (test set only; shared with SF41_PTCSubclassPerformance) |

## Regenerating

```bash
Rscript data_export.R
python3 figure_s_branch_shap_by_subclass.py
```

## Cross-references

- **Figure 5 Panel C** — same branch decomposition pooled across the
  whole NMD-positive test set. This SF reproduces Panel C's style
  three times, one per PTC subclass at the n = 819 scope. Panel C's
  pooled NMD numbers (Structural 60.7% / Stop 28.8% / AUG 10.5%) sit
  between this SF's NMD+/PTC+ and NMD+/PTC− columns, weighted toward
  PTC+ by sample size.
- **`SF41_PTCSubclassPerformance/`** — companion SF for the second half of
  the same paragraph (predictive performance lower in NMD+/PTC−).
  That SF stays test-set scoped.
- **Figure 5 Panel G** — same n = 1,016 / 95 / 1,107 subgroup counts;
  that panel uses attention-aggregator output, this SF uses fusion-
  layer KernelSHAP. Two complementary interpretability angles on the
  same scope.
