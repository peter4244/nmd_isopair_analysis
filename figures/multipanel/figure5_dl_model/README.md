# Figure 5 — Deep-learning model for NMD prediction

Section 5 multipanel composite. The model is trained in
[`peter4244/NMD_orf_model_v5_4ct`](https://github.com/peter4244/NMD_orf_model_v5_4ct)
on the chr-1/3/5/7 paralog-free held-out test set; this folder turns its
exports into Figure 5.

## Panels

| Panel | What | Source |
|---|---|---|
| **A** | Model architecture schematic | `figure5_panelA_architecture.R` (sources `make_architecture_figure.R` from the model repo) |
| **B** | ROC + PR inset (AUC = 0.93, AUPRC = 0.83) | `figure5_panelB_performance.py` |
| **C** | KernelSHAP branch decomposition (Structural / Stop / ATG) | `figure5_panelC_branch_importance.py` |
| **D** | Mean \|SHAP\| per structural feature | `figure5_panelD_structural_features.py` |
| **E** | Per-nucleotide signed SHAP × input around AUG, NMD class | `figure5_panelE_atg_logo.py` (reads `data/motif_logo_atg_joint_atg500_stop500.tsv`) |
| **F** | Same for stop codon | `figure5_panelF_stop_logo.py` (reads `data/motif_logo_stop_joint_atg500_stop500.tsv`) |
| **G** | Path B-strict uORF attention at n=1,166 Subset 2 | `figure5_panelG_uorf_attention.py` + `data_export_n1166.R` |

## Composite

`figure5_composite.py` lays out A in the left column (preserving its
native 0.86 aspect) and B–G in a 3 × 2 grid on the right (cells 6 × 4 in,
matching the 1.5:1 native aspect of B/C/D/G; E and F at native 2:1 sit
inside their cells with small top/bottom whitespace). Output:
`figure5_composite.{pdf,png}` at 22.29 × 12.00 in. Layout is checked by
`validate_multipanel_layout` (passes clean).

Legend prose: [`figure5_composite_legend.md`](figure5_composite_legend.md).

## Status

All seven panels (A–G) are final. Data and renders match the model
exports at `~/claude_projects/NMD_orf_model_v5_4ct/results_4ct/`. Panels
E and F are rendered natively with `logomaker` against the pooled
5-run DeepSHAP motif-logo TSVs (`data/motif_logo_{atg,stop}_joint_atg500_stop500.tsv`);
they share the same coordinate convention (first nucleotide of the codon
= +1, no position 0) so the two windows read in parallel.

Style: all matplotlib panels use the shared ggplot-mimic shim
(`figures/lib/ggplot_style.py`) and `assert_text_within_canvas` before
`savefig`. Layout is checked by `validate_multipanel_layout`.

## Regenerating

```bash
# from this directory
Rscript data_export_n1166.R                       # G inputs
Rscript figure5_panelA_architecture.R             # A
python3 figure5_panelB_performance.py             # B
python3 figure5_panelC_branch_importance.py       # C
python3 figure5_panelD_structural_features.py     # D
python3 figure5_panelE_atg_logo.py                # E (placeholder)
python3 figure5_panelF_stop_logo.py               # F (placeholder)
python3 figure5_panelG_uorf_attention.py          # G
python3 figure5_composite.py                      # composite
```

## Cross-references

- Manuscript §5 (Google Doc): see `paper/results_to_code_map.md` §5 entry.
- Model architecture, training, METHODS: `NMD_orf_model_v5_4ct/METHODS.md`.
- TD2 bias context for Panel G's strict-uORF rule:
  `figures/SupplementalFigures/SF34-SF35_TD2BiasEvidence/README.md`,
  `TD2_BIAS_AUDIT.md` (repo root).
