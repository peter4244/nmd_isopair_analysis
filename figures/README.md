# figures/ — shared figure tooling and style guide

Shared infrastructure for publication-quality figures in the NMD long-read paper. Pete and Yul both pull from here; updates apply to both sides.

## What's in `lib/`

| File | Purpose |
|---|---|
| [`principles.md`](lib/principles.md) | **Read first.** Language-agnostic figure principles — validator-as-bar, two-font hierarchy, structural symmetry, snapshot triplets, reference-first, data-integrity rules. |
| [`matplotlib.md`](lib/matplotlib.md) | Python / matplotlib addendum — Arial font setup, `figure_primitives` / `figure_geometry` usage, `bbox_inches='tight'` quirk. |
| [`ggplot2_patchwork.md`](lib/ggplot2_patchwork.md) | R / ggplot2 + patchwork addendum — validators, patchwork gotchas, inline-R-values workflow. |
| [`figure_primitives.py`](lib/figure_primitives.py) | Reusable drawing primitives: `aim_card`, `ds_box`, `block_arrow_right/left`, `thin_arrow`, `render_bullets`, `exon_pill`, `intron_segment`, airway biology (`draw_ali_section`, `draw_cilia`, `smoke_puff`), etc. |
| [`figure_geometry.py`](lib/figure_geometry.py) | `horizontal_layout` / `vertical_layout` for derived coordinates; `geometry_summary` for pre-render balance checks. |
| [`validate_figure_layout.py`](lib/validate_figure_layout.py) | Python layout validator — font-size escapes, shape overlaps, off-axis text. Run after every render. |
| [`validate_figure_layout.R`](lib/validate_figure_layout.R) | R layout validator for single annotation-based ggplots (architecture diagrams, flowcharts). |
| [`validate_composite_layout.R`](lib/validate_composite_layout.R) | R validator for patchwork composite figures (tag-title collisions, tag overlap, tag font consistency). |

## Quickstart

### Python / matplotlib

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path("figures/lib").resolve()))  # adjust if running from a deeper directory

import matplotlib.pyplot as plt
plt.rcParams['font.family']     = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica Neue', 'Helvetica', 'DejaVu Sans']
plt.rcParams['pdf.fonttype']    = 42

from figure_primitives import aim_card, block_arrow_right
from figure_geometry import horizontal_layout, geometry_summary
from validate_figure_layout import validate_figure_layout

# ... build figure ...

result = validate_figure_layout(fig, ax, verbose=False)
assert result['summary']['n_errors'] == 0
plt.savefig("figure.pdf", bbox_inches='tight')
```

See [`lib/matplotlib.md`](lib/matplotlib.md) for the full workflow.

### R / ggplot2 + patchwork

```r
# Single annotation-based figure:
source("figures/lib/validate_figure_layout.R")
validate_figure_layout(p, fig_width = 16, fig_height = 24)

# Patchwork composite:
source("figures/lib/validate_composite_layout.R")
validate_composite_layout(combined, fig_width = 16, fig_height = 18)
```

See [`lib/ggplot2_patchwork.md`](lib/ggplot2_patchwork.md) for the full workflow.

## Workflow path

The figure scripts in this directory and the canonical analysis Rmd at [`results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_gencode_scope_2026-06-15.Rmd`](../results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_gencode_scope_2026-06-15.Rmd) are **parallel renderers** from the same upstream RDS caches in `results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr/`:

```
data_mashr/ (RDS caches: profiles_c{2,4}, cds, structures, ref_atg_analysis,
             utr5_features_{all,refaug}, expression_data, nmd_classification, ...)
                  │
        ┌─────────┴──────────┐
        ▼                    ▼
  figure scripts        new Rmd
  (this directory)      (canonical analysis report)
        │                    │
        ▼                    ▼
  panel PNGs + TSVs    HTML with inline-R values + embedded panels
        │                    │
        └────────┬───────────┘
                 ▼
   reproducibility/verify_cross_check_new_rmd_vs_figures.R
   (confirms both pipelines agree to the digit; 57 checks)
```

- The **Rmd is the canonical analysis report** — its inline-R values are the primary record of what was computed at which scope, with a full prose narrative anchoring each finding.
- The **figure scripts** in `figures/multipanel/...` and `figures/SupplementalFigures/...` are downstream renderers. Each panel script reads from the same RDS caches the Rmd uses, computes its own descriptives, and writes a small TSV alongside the PNG.
- **Drift detection** between the two sides is automated: [`reproducibility/verify_cross_check_new_rmd_vs_figures.R`](../reproducibility/verify_cross_check_new_rmd_vs_figures.R) parses the Rmd HTML and the figure-side TSVs, then asserts they agree at every shared number (counts exact, percentages ±0.1, medians ±0.5, p-values same order of magnitude). On any FAIL the verifier exits non-zero — diagnose root cause first per v4 plan §S4, never paper over a discrepancy.
- The **legacy report** ([`05_final_report_mashr.Rmd`](../results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_mashr.Rmd)) remains in the tree, legacy-bannered, for sensitivity analyses against the prior submission. Its companion verifier is [`reproducibility/verify_legacy_rmd_reproducibility.R`](../reproducibility/verify_legacy_rmd_reproducibility.R) (deprecated; routine use is the two new verifiers).

## How updates flow

This directory mirrors Pete's personal figure-tooling helpers (originally at `~/.claude/utils/` and `~/.claude/memory/figures_*.md`). The version here is the **shared, repo-tracked** copy — when a primitive or validator is improved by either contributor, the canonical version is the one here.

If you find a bug or add a primitive useful for the manuscript figures, edit the file in this directory, commit, push. The change becomes visible to all collaborators on next `git pull`.

## The figure-discipline summary

The full principles doc lives in [`lib/principles.md`](lib/principles.md). The short version:

- **Validator clean** before showing anyone — every render gets a validator pass and every error gets resolved.
- **Two font sizes only.** Headers bold; body never bold.
- **Structural symmetry** via named coordinate constants and `horizontal_layout`/`geometry_summary` (Python) or `validate_composite_layout()` (R).
- **Define population and denominator BEFORE building.**
- **Show all data points.** Truncated axes must be disclosed.
- **Reference image first** for biology panels.
- **Snapshot triplets** at "close-but-not-final" milestones.
- **Inline computed values** (`` `r ...` `` in R, `f"{x:.2f}"` in Python) — never hardcode results.
- **Empty space → compress, don't pad.** Shrink panels, trim margins, remove titles; don't add filler.
- **Mockup / brainstorming stage:** rigor is relaxed — fast iteration matters more than validator-clean output. Full rigor returns when the figure is heading into the submitted draft.
