# Scientific Figure Guide for R / ggplot2 + patchwork

R/ggplot-specific addendum for publication-quality figures.

**Read first:** [`principles.md`](principles.md) — language-agnostic principles (validator-as-bar, two-font hierarchy, structural symmetry, snapshot triplets, reference-first, deliberate symmetry pass, data-integrity rules). This file only covers R/ggplot/patchwork-specific tooling and quirks.

Python sibling: [`matplotlib.md`](matplotlib.md).

## Validators

Both run after every render; fix all errors before presenting.

**Single annotation-based figures** (architecture diagrams, flowcharts, schematics):
```r
source("figures/lib/validate_figure_layout.R")
validate_figure_layout(p, fig_width = 16, fig_height = 24)
```
7 checks: text-segment collision, text-text overlap, proximity crowding, text-rect centering, horizontal alignment, arrow minimum length, text-rect padding/overflow.

**Patchwork multipanel composites:**
```r
source("figures/lib/validate_composite_layout.R")
validate_composite_layout(combined, fig_width = 16, fig_height = 18)
```
4 checks: tag-title overlap (especially inside `wrap_elements`), tag-tag overlap, missing/duplicate tags, tag font consistency.

The validators catch spatial problems but not aesthetic ones — human review is still needed.

## Patchwork-specific gotchas

- **Don't use `plot_annotation(title = ...)` inside `wrap_elements()`.** It collides with patchwork tags. Use `labs(title = ...)` on the inner plot instead. This is the most common source of `validate_composite_layout()` errors.
- **Tag font size and face must be consistent across panels.** Set theme tag styling once at the composite level, not per-plot.
- **Shared legends are factored out where possible** (`patchwork::plot_layout(guides = "collect")`).

## ggplot conventions

- **Publication theme applied consistently** across all panels (single `theme_*()` definition).
- **Axis labels are informative,** not raw column names.
- **Statistical test where appropriate** — Wilcoxon for distributions, Fisher's exact for proportions. Annotate p-values on the panel; full test details in the figure legend.
- **Inline R values** for any computed quantity in prose: `` `r ...` `` referencing objects from earlier chunks. Never hardcode computed values; never present placeholder numbers as real results.
- **If truncating/clipping axes:** state the threshold and number of clipped values in the legend or surrounding prose.

## Multipanel build checklist

### Before building
- [ ] Population and denominator defined (see principles doc)
- [ ] Plan panel layout — which comparisons go side-by-side, what shares an axis

### Each panel
- [ ] Statistical test applied where appropriate
- [ ] Publication theme consistent
- [ ] Axis labels informative
- [ ] All data points shown (no silent subsetting)

### Composition
- [ ] `validate_composite_layout()` clean
- [ ] Panel tags A, B, C, ... sequential, no gaps
- [ ] Tag font size + face consistent
- [ ] Shared legends collected
- [ ] No `plot_annotation(title=)` inside `wrap_elements()`

### Prose / captions
- [ ] Population description appears BEFORE the figure
- [ ] All numbers via inline R; no hardcoded values
- [ ] Biological interpretation verified computationally; null expectation computed first

## Common mistakes specific to R workflows

- `plot_annotation(title = ...)` inside `wrap_elements()` — tag-title collision. Use `labs(title = ...)` on the inner plot.
- Stray ``` fences from editing — verify chunk boundaries after edits.
- Loading objects in one section that are used in another — load shared objects in a setup chunk.
- Hardcoding computed values into prose — always use inline R.
