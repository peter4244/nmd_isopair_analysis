---
name: Figure-making principles (language-agnostic)
description: Shared principles for publication / grant figures — apply to both R/ggplot and Python/matplotlib workflows
type: reference
---
# Figure-making principles (language-agnostic)

Applies to publication / grant figures regardless of tooling. Language-specific addenda:
- **R / ggplot2 + patchwork:** [`ggplot2_patchwork.md`](ggplot2_patchwork.md)
- **Python / matplotlib:** [`matplotlib.md`](matplotlib.md)

## The bar before showing the user

- **Validator clean.** Run the layout validator after every render and resolve every error before presenting.
  - R: `validate_figure_layout()` and `validate_composite_layout()` from `figures/lib/`.
  - Python: `validate_figure_layout()` from `figures/lib/validate_figure_layout.py`.
  - When new failure modes appear, **extend the validator**; don't work around it.
- **Deliberate symmetry pass.** Before sharing, ask: does the left mirror the right? Do siblings (rows, panels, bullet blocks) match each other in height, in right-edge of wrapped text, in arrow position, in vertical center lines? Treat this as its own pass — not blended with other review.

## Typography

- **Two font sizes only.** One header size, one body size — no third. Mixing more sizes weakens the hierarchy. Typical: HEADER ≈ 18, BODY ≈ 14.
- **Headers bold; body never bold.** Mixing bold body with bold headers kills the visual hierarchy.
- **Tag/label fonts consistent across panels.** Same size and face for A/B/C tags everywhere.
- **No fontsize literals in panel code.** Everything — xticks, yticks, xlabel, ylabel, annotations, medians, brackets — imports `BODY_FS` (or `HEADER_FS`) from `figures/lib/ggplot_style.py`. Local overrides like `fontsize=11`, `LABEL_FS = 12`, or `BODY_FS - 2` create silent cross-panel drift and defeat the single-source-of-truth the module provides. If a size feels wrong, bump the shared constant for this panel; don't sprinkle a one-off.
- **Native font size scales with figure width.** The docx pipeline scales every figure to a 6.5 " content width. A 14 pt native BODY_FS on an 11 " render reads at only 8.3 pt in the delivered docx. Every panel that will end up in the paper's supplemental doc sets its own `BODY_FS = docx_body_fs(NATIVE_W)` at the top so the effective docx size hits the 10 pt target regardless of native width. Floor: 9 pt effective. Ceiling: ~14 pt effective. `assert_docx_readable(fig)` catches violations at write time.
- **Font-size symmetry is a validator target, not a review target.** Call `assert_style_symmetric(fig)` alongside `assert_text_within_canvas(fig)` before every `savefig`. It fails if sibling axes disagree on the same role (xtick vs xtick, ylabel vs ylabel).

## Layout & symmetry

- **Symmetry must be structural, not coincidental.** Don't declare each coord as a free parameter — derive coords from a small set of primitives (page padding, column widths, arrow width, gap) and compute the rest. Left and right arrows then share the same `ARROW_W` and cannot drift apart from a typo.
- **Named coordinate constants up front.** Columns, rows, shared sizes (header height, stripe width). Y-centers like `A1_YC = (A1_Y0 + A1_Y1) / 2` keep arrows aligned across panels.
- **Symmetric arrow convergence.** When multiple arrows converge or fan out, make the pattern symmetric about the center; the middle element should go straight down. Spread arrival points across the target's edge rather than collapsing into one blob.
- **Uniform visual weight.** Don't progressively fade outer elements; uniform weight reads cleaner.
- **Dashed enclosing boxes** communicate "this unit applies repeatedly" without drawing multiple copies.

## Color

- **Color used semantically.** Related elements share a color family. Generic processing stages use neutral grey. Don't assign arbitrary colors to unrelated elements.
- **Color must encode meaning that survives B/W printing or re-coloring.** Don't rely on color alone for category distinctions in critical contrasts; pair with shape, label, or position.

## Data integrity

- **Define population and denominator before building.** State observations, filtering, and the denominator explicitly. Population description appears BEFORE the figure in surrounding prose.
- **State specific counts.** "1,331 isoforms," not "a subset of isoforms." Use inline values from the analysis, not narrative approximations.
- **Show all data points.** Never silently drop subsets. If some data doesn't fit a panel, show it in a companion panel or table.
- **Truncated/clipped axes must be disclosed.** State the threshold and the count of clipped values in the legend or surrounding prose.
- **Mixed denominators in one panel** = split the panel or make the subsetting explicit.

## Workflow

- **Code first, render, examine, then write prose.** Never write interpretive narrative before computing the result. Verify every biological interpretation computationally; if the observed pattern matches random expectation, say so.
- **Reference image before iterating on a biology panel.** Find / ask for a reference (BioRender-style cell-type panel, micrograph) and match colors, proportions, and morphology faithfully. Trying to guess wastes iterations.
- **Snapshot triplets at "close-but-not-final" milestones.** When the user says "save this version, it's close" — copy the script to a named variant and re-export the rendered files alongside before the next iteration. Cheap insurance.
- **Never present placeholder values as real.** If a result isn't computed, say so explicitly.

## Common failure modes

- Mixed denominators in one panel.
- Truncated axes without disclosure.
- Interpretive prose written before the underlying computation.
- Free-floating coordinates that drift apart across edits (cause: not derived from primitives).
- Multiple font sizes accumulating across iterations.
- Body text rendering bold because of font-fallback quirks (verify with a test render).
- "Looks symmetric" without programmatic verification — use the language-specific geometry/balance tools.
