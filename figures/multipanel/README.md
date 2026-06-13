# `figures/multipanel/` — final manuscript multipanel figures

This folder holds the **manuscript multipanel figures** for the NMD long-read paper. The composition workflow matches the grant work in `~/claude_projects/grants/nmd_2026/` (matplotlib + the shared `figures/lib/` helpers + validator-clean output bar).

## Folder structure

```
multipanel/
├── README.md                    # This file
├── figureN_<short-name>/
│   ├── figureN_composite.py     # Composition script (loads panel data → matplotlib subplot grid)
│   ├── figureN_composite.pdf
│   ├── figureN_composite.png
│   ├── figureN_panelA_<descr>.py    # Individual panel render scripts
│   ├── figureN_panelA_<descr>.pdf
│   ├── figureN_panelA_<descr>.png
│   ├── figureN_panelB_<descr>.py
│   ├── figureN_panelB_<descr>.pdf
│   ├── figureN_panelB_<descr>.png
│   └── ...
└── figureM_<short-name>/
    └── ...
```

**Rule of thumb:** one folder per manuscript multipanel figure. Each folder is self-contained — composite script + per-panel scripts + their PDF/PNG outputs. Per-panel outputs are kept so individual panels can be reused in supplements, talks, or other contexts without re-running the composite.

## Shared infrastructure (`figures/lib/`)

- `figure_primitives.py` — common matplotlib primitives.
- `figure_geometry.py` — `horizontal_layout()` + `geometry_summary()` for sizing.
- `validate_figure_layout.py` — layout validator. **A figure is not "done" until validator returns zero warnings.**
- `principles.md` — language-agnostic figure conventions (two-font-only, structural symmetry, etc.).
- `matplotlib.md` — Python addendum (Arial+font setup, `bbox_inches='tight'` quirk, etc.).
- `ggplot2_patchwork.md` — R-side addendum, for figures still in R.

## Per-figure naming convention

- Composite filename: `figureN_composite.{py,pdf,png}` (no panel suffix).
- Panel filename: `figureN_panel<L>_<short_descriptor>.{py,pdf,png}` where `<L>` is the panel letter (A, B, C, …) and `<short_descriptor>` is snake_case.
- Example for a hypothetical Figure 3 with panels A/B/C: `figure3_panelA_ptc_distance.py`, `figure3_panelB_event_prevalence.py`, `figure3_panelC_attribution.py`, plus `figure3_composite.py`.

## Current figures

| Folder | Manuscript figure | Status |
|---|---|---|
| `figure3_isopair_and_ptc/` | Combined Isopair + PTC figure (= old Fig 3 + old Fig 4) | **Built (2026-06-13).** 6 panels (A–F) at ref-AUG-traceable scope; composite 3×2; 6-pass verification protocol passed; published headline 83.5% NMD PTC vs 16.3% Control. |
| `figure4_ptcneg_and_model/` | New Fig 4 = old Fig 5 + old Fig 6 merge | **Placeholder.** Per Pete's plan: old 5A + 5B dropped; 5F + 6A → supplement; new Fig 4 ≈ 8 panels (3 from old 5 + 5 from old 6). Awaiting kickoff. |
| `figure5_attention_and_subgroup/` | New Fig 5 = old Fig 7 + old Fig 8 merge | **Placeholder.** Old Fig 7 (attention distribution) + old Fig 8 (mechanistic subgroup analyses) — both still needed captions in the manuscript when last checked. Awaiting kickoff. |
| (TBD) | Figure 1 (Isoform Landscape) | TBD |
| (TBD) | Figure 2 (Output Lost + PCI) | TBD |

Net effect of Pete's plan: 8 original manuscript figures → 5 in the new layout (1, 2, new 3, new 4, new 5). See `paper/results_to_code_map.md` for the per-claim provenance map across the new layout.

## Validation workflow (matches grant pattern)

For every composite or panel render:

1. `bash`: `python figureN_composite.py` (or panel script).
2. `bash`: `python ~/claude_projects/nmd/figures/lib/validate_figure_layout.py figureN_composite.pdf`.
3. Confirm zero warnings before declaring "done."
4. Snapshot intermediate states (overview + detail + comparison) at close-to-final milestones for Pete's review.

For data integrity: every number in the figure must be reproducible from a tracked script — no manual annotations baked into the figure file.
