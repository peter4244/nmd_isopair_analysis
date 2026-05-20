---
name: Matplotlib figure conventions
description: Python/matplotlib-specific addendum — read with figures_principles.md
type: reference
---
# Matplotlib figures — Python addendum

**Read first:** [`principles.md`](principles.md) — language-agnostic principles (validator-as-bar, two-font hierarchy, structural symmetry, snapshot triplets, reference-first, deliberate symmetry pass, data-integrity rules). This file only covers matplotlib-specific tooling and quirks.

R sibling: [`ggplot2_patchwork.md`](ggplot2_patchwork.md).

## Tooling at `figures/lib/`

- `validate_figure_layout.py` — layout validator. Checks font-size escapes, shape overlaps (rect-rect, rect-polygon, rect-circle), off-axis text. Run after every render.
- `figure_primitives.py` — reusable drawing primitives.
- `figure_geometry.py` — `horizontal_layout` / `vertical_layout` for derived coords; `geometry_summary` for balance checks.

## Matplotlib-specific quirks

### Use Arial, not Avenir Next

On macOS, Avenir Next exposes only its 700 (Bold) face to matplotlib's font manager, so body text rendered with `fontweight='normal'` falls back to Bold and silently looks bold. Arial has both Regular (400) and Bold (700) faces visible to matplotlib, so the weight hierarchy actually shows. Set early:
```python
import matplotlib.pyplot as plt
plt.rcParams['font.family']     = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica Neue', 'Helvetica', 'DejaVu Sans']
plt.rcParams['pdf.fonttype']    = 42
plt.rcParams['ps.fonttype']     = 42
```

### `bbox_inches='tight'` does not trim everything

It includes `ax.patch`, so empty regions inside the axes data limits stay in the saved image. To remove dead space at the bottom or top, change `ax.set_ylim()` rather than `figsize`.

## Reusable primitives

```python
import sys
from pathlib import Path
# Adjust if running from a deeper directory:
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures" / "lib"))
from figure_primitives import (
    rbox, thin_arrow, block_arrow_right, block_arrow_left,
    aim_card, ds_box, render_bullets, render_insights_timeline,
    exon_pill, intron_segment, smoke_puff, draw_cilia, draw_ali_section,
)
```

What's there:
- **Cards:** `aim_card` (white body + colored stripe + header band), `ds_box` (data-source box: title + italic subtitle).
- **Arrows:** `block_arrow_right`, `block_arrow_left` (filled block arrows), `thin_arrow` (line annotation arrow).
- **Content:** `render_bullets` (labeled bullets), `render_insights_timeline` (numbered colored dots + connector + bullet text).
- **Gene model:** `exon_pill`, `intron_segment`.
- **Airway biology:** `draw_ali_section` — pseudostratified HBEC ALI cross-section with basal/club/ciliated/goblet cells, three stages (`basal`, `inter`, `full`); `draw_cilia`; `smoke_puff`.

Add new generic primitives to the module, not the figure script.

## Symmetry — Python tooling

The principles doc explains *why* symmetry must be structural; this is *how* in matplotlib:

```python
from figure_geometry import horizontal_layout, geometry_summary
SIDE_W, MID_W, ARROW_W, GAP, PAD = 4.00, 4.20, 0.70, 0.05, 0.30
coords, _ = horizontal_layout(PAD, [
    ('SAEC', SIDE_W), ('gap_l1', GAP), ('arrow_l', ARROW_W), ('gap_l2', GAP),
    ('INS',  MID_W),  ('gap_r1', GAP), ('arrow_r', ARROW_W), ('gap_r2', GAP),
    ('LTRC', SIDE_W),
])
geometry_summary(
    page=(FIG_W, FIG_H),
    columns=[('SAEC', *coords['SAEC']), ('INS', *coords['INS']), ('LTRC', *coords['LTRC'])],
    arrows=[('L', *coords['arrow_l']), ('R', *coords['arrow_r'])],
    rows=[('A1', A1_Y0, A1_Y1), ('A2', A2_Y0, A2_Y1), ('A3', A3_Y0, A3_Y1)],
    require_equal_widths=[['SAEC', 'LTRC']],
    require_equal_heights=[['A1', 'A2', 'A3']],
    require_arrows_equal=True,
    require_symmetric_gaps=True,
)
```

`geometry_summary` prints column widths, gaps, arrow widths, and row heights, then verifies declared equalities. Run before plotting; a glance tells you whether the layout balances without rendering.

## Suggested figure-script skeleton

```python
import matplotlib.pyplot as plt
from figure_primitives import aim_card, block_arrow_right, ...
from figure_geometry import horizontal_layout, geometry_summary
from validate_figure_layout import validate_figure_layout

# 1. Font setup (Arial, two sizes)
plt.rcParams['font.sans-serif'] = ['Arial', ...]
HEADER_FS, BODY_FS = 18, 14

# 2. Geometry from primitives
SIDE_W, MID_W, ARROW_W, GAP, PAD = ...
coords, _ = horizontal_layout(PAD, [...])

# 3. Pre-render geometry check
geometry_summary(page=..., columns=..., arrows=..., rows=...,
                 require_equal_widths=..., require_arrows_equal=True)

# 4. Build figure with primitives
fig, ax = plt.subplots(...)
aim_card(ax, ...); block_arrow_right(ax, ...); ...

# 5. Validate, then save
result = validate_figure_layout(fig, ax, verbose=False)
assert result['summary']['n_errors'] == 0
plt.savefig(..., bbox_inches='tight')
```
