# Supplemental Figure — DL discrimination by PTC subclass at n = 1,166

Companion to the §5 manuscript sentence:

> "However, overall predictive performance was substantially lower in
> NMD+/PTC− isoforms (SFx)."

## What it shows

Held-out (chr-1/3/5/7 paralog-free) discrimination of the deep-learning
NMD predictor, restricted to the n = 1,166 ref-AUG-traceable subset
(same scope as Figure 4 C/D and Figure 5 Panel G) and stratified into
NMD+/PTC+, NMD+/PTC− retained, and Control.

Two panels:

- **(A)** ROC curves for the two NMD-vs-Control contrasts.
- **(B)** Per-isoform predicted NMD probability per subgroup.

## Headline numbers

| Contrast | AUC | AUPRC | n_pos vs n_neg |
|---|---|---|---|
| NMD+/PTC+ vs Control | 0.96 | 0.94 | 255 vs 276 |
| NMD+/PTC− vs Control | 0.74 | 0.26 | 30 vs 276 |

Mean predicted NMD probability per subclass:
- Control: 0.19
- NMD+/PTC+: 0.85
- NMD+/PTC−: 0.39 (straddles the 0.5 decision threshold)

The manuscript's "substantially lower" claim is well-supported by the
0.96 → 0.74 AUC drop.

## Scope note — test-set only

Predictions used here are out-of-fold predictions on the chr-1/3/5/7
paralog-free test split. The NMD+/PTC− subgroup has n = 30 at the test
split. This is per the test-set-only policy for performance metrics
([feedback_nmd_analysis_scope_test_vs_all](~/.claude/projects/-Users-petecastaldi/memory/feedback_nmd_analysis_scope_test_vs_all.md)).

## Files

| File | What |
|---|---|
| `figure_s_performance_by_subclass.py` | Two-panel matplotlib figure (validator-clean) |
| `figure_s_performance_by_subclass_legend.md` | Manuscript-style legend |
| `data/predprob_by_subclass_refaug.tsv` | Per-isoform predicted probability + group (shared with `PTCSubclassBranchSHAP/`; built by `PTCSubclassBranchSHAP/data_export.R`) |

## Regenerating

The data file is produced by the sibling SF's data export:

```bash
cd ../PTCSubclassBranchSHAP && Rscript data_export.R
cd ../PTCSubclassPerformance && python3 figure_s_performance_by_subclass.py
```

## Cross-references

- **Figure 5 Panel B** — pooled (no subclass) ROC + PR on the same test
  split: AUC = 0.93, AUPRC = 0.83. The pooled number sits between the
  per-subclass numbers here, weighted by the larger PTC+ subgroup.
- **`PTCSubclassBranchSHAP/`** — companion SF for the same paragraph;
  explains *why* PTC− discrimination is lower (loss of stop and
  structural branch signal; the model leans more on the AUG branch but
  not enough to fully compensate).
