# `figures/lib/` — shared figure tooling

**Before editing anything in this directory or drafting a figure that imports from here, open `~/.claude/memory/figures_workflow_INDEX.md`.** That's the universal entry point — decision gate, doc map, known failure modes.

## What lives here

| File | Purpose |
|---|---|
| `ggplot_style.py` | Style SSOT — `BODY_FS`, `HEADER_FS`, palette (`NMD_COLOR`, `CONTROL_COLOR`, …), `docx_body_fs()`, `apply_ggplot_rcparams()`, `style_axes_ggplot()`, `panel_letter()`, `assert_text_within_canvas()`. Every matplotlib figure imports from here. |
| `validate_figure_layout.py` | Six validators — `validate_figure_layout`, `validate_multipanel_layout`, `assert_style_symmetric`, `assert_docx_readable`, `assert_data_free`, plus `find_open_regions` for query-driven placement. |
| `validate_figure_layout.R` | R port of `validate_figure_layout` for ggplot2 panels. |
| `validate_composite_layout.R` | Patchwork composite validator. |
| `validate_flowchart_dot.R` | Static DOT / Graphviz flowchart validator (SF29-class). |
| `figure_primitives.py` | Reusable primitives for grant-style diagrams (aim cards, arrows, section frames). |
| `figure_geometry.py` | Named-coord constants + `geometry_summary` for structural symmetry. |
| `mechanism_class.R` | PTC-mechanism classification helper for `nmd_responsive` isoforms. |
| `p_to_stars.{py,R}` | p-value → significance-stars formatter (matches Fig 3–5 convention). |
| `principles.md`, `matplotlib.md`, `ggplot2_patchwork.md` | Project-local addenda that supplement the universal principles doc. |

## Universal docs (read these before writing panel code)

1. `~/.claude/memory/figures_workflow_INDEX.md` — start here (decision gate + doc map).
2. `~/.claude/memory/figures_workflow_publications.md` — 9-phase sequenced runbook.
3. `~/.claude/memory/figures_principles.md` — universal principles.
4. `~/.claude/memory/figures_style_publications.md` — publication-specific conventions.
5. `~/.claude/memory/figures_matplotlib.md` — Python/matplotlib tool addendum.
6. `~/.claude/memory/figures_validators.md` — validator catalog.

## Canonical templates

Copy from one of these when starting a new figure:

- **2×2 with centered bottom panel (matplotlib):** `../SupplementalFigures/SF34_TD2Bias_broad/figure_sf34.py`
- **1×2 boxplot + KDE (matplotlib):** `../SupplementalFigures/SF39_AttentionDistribution/figure_s_attention_distribution.py`
- **12-panel schematic composite (R/ggplot + patchwork):** `../SupplementalFigures/SF25_SpliceEventCategories/render_sf25_canonical.R`
- **DOT / Graphviz flowchart (R):** `../SupplementalFigures/SF29_PairAnalysisFlowchart/build_flowchart.R`

## Project-specific overrides

- Palette: `~/.claude/projects/-Users-petecastaldi/memory/project_nmd_figure_palette.md`
- SF numbering + supplement structure: `~/.claude/projects/-Users-petecastaldi/memory/reference_nmd_supplemental_figures_workflow.md`
- Docx build pipeline: `paper/build_supplemental_figures_docx.js` (invoked via `cd paper && node build_supplemental_figures_docx.js`)
- Results-to-code trace: `paper/results_to_code_map.md`
