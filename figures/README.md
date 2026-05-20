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
