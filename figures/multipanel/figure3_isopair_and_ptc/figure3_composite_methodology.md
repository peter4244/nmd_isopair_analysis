# Figure 3 composite — methodology

**Render script:** `figure3_composite.py`
**Outputs:** `figure3_composite.{pdf,png}`
**Layout:** 2 columns × 3 rows (portrait)
**Scope:** Composition only; each panel's data + computation are documented in its own `figure3_panel{A,B,C,D,E,F}_methodology.md`

**Canonical analysis report:** the numbers in Panels A–F are also computed inline (independently) by [`results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_gencode_scope_2026-06-15.Rmd`](../../../results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_gencode_scope_2026-06-15.Rmd) — Panels A/B/C in §1, Panel D in §2a, Panels E + F in §2b. Both the figure scripts here and the Rmd consume the same upstream RDS caches; the cross-check verifier at [`reproducibility/verify_cross_check_new_rmd_vs_figures.R`](../../../reproducibility/verify_cross_check_new_rmd_vs_figures.R) confirms they agree to the digit.

## What this script does

Reads the 6 pre-rendered panel PNGs (`figure3_panel{A,B,C,D,E,F}_*.png`) and composes them into a single multipanel figure. The composite does no data work — it's a layout-only deliverable.

| Cell | Panel |
|---|---|
| Row 1, Col 1 | A — Isoform pair construction (schematic) |
| Row 1, Col 2 | B — Sequence shared with reference |
| Row 2, Col 1 | C — Splice event prevalence |
| Row 2, Col 2 | D — Stop codon to last EJC distance |
| Row 3, Col 1 | E — PTC-causing event attribution |
| Row 3, Col 2 | F — PTC-causing splice events (mechanism) |

## Layout decisions

The composite is built around a single design rule: **every panel renders at the same aspect ratio, so the grid has no per-row scaling and no distortion.** This is enforced upstream in each panel's render script — not papered over in the composite.

- **Cell geometry: 6.0" × 4.0" (1.5:1).** Each panel script sets `figsize=(6, 4)` and saves *without* `bbox_inches="tight"` (which would crop the saved PNG to a content-dependent aspect ratio and was the source of distortion in the prior 2×3 build). The output PNG is therefore guaranteed to be exactly 1800 × 1200 px at 300 dpi.
- **Composite layout:** plain `GridSpec(3, 2)` with equal cells. `aspect="auto"` on imshow stretches the source PNG to fill the cell — but since cell and source share the same 1.5:1, the stretch ratio is 1:1 and no distortion occurs.
- **Composite figure size:** `FIG_W = 2 × 6.0 = 12.0"`, `FIG_H = 3 × 4.0 = 12.0"` — square. PNG = 3,600 × 3,600 px at 300 dpi.
- **Panel labels `A`–`F`:** 22 pt bold, top-left corner of each cell (transform=ax.transAxes). Matches grant `figure_nmd_prelim_composite_3by1.py` convention.
- **Margins compressed** (`left=0.005, right=0.995, top=0.985, bottom=0.005`, `hspace=0.04`, `wspace=0.02`) so the panel content fills the canvas with breathing room only for the letter labels.
- **No internal frame around cells** (spines disabled, ticks hidden). The composite reads as 6 distinct panels by virtue of their content and label placement.

## Per-panel adjustments needed to hit 1.5:1

The data panels (B/C/D/E/F) compress cleanly from their old larger figsizes to 6×4 — `subplots_adjust()` controls axes-within-canvas placement so axis labels, tick labels, titles, and legends fit at the smaller size. The values are per-panel because each plot has different label real-estate needs:

| Panel | `subplots_adjust(left, right, bottom, top)` | Rationale |
|---|---|---|
| B (KDE) | 0.13, 0.97, 0.16, 0.87 | minimal — axis labels + title |
| C (vertical bar, 10 rotated labels) | 0.12, 0.97, **0.38**, 0.87 | bottom generous, 35°-rotated category labels |
| D (KDE) | 0.13, 0.97, 0.16, 0.87 | same as B |
| E (vertical bar, 11 rotated labels + subtitle) | 0.12, 0.97, **0.38**, **0.80** | top lower for title + subtitle stack; bottom for rotated labels |
| F (horizontal stacked bar + bottom legend) | **0.22**, 0.97, **0.30**, 0.80 | left for y-axis category labels; bottom for legend below x-label (legend `bbox_to_anchor=(0.5, -0.42)`) |

Panel A (schematic) required more invasive changes — at 6×4 the relative font size shrinks the available horizontal label-space, so:
- Font sizes reduced (`HEADER_FS=15, BODY_FS=11` vs 18/14 elsewhere) so left-of-track labels fit.
- Bracket-pair gap widened (`BRACKET_CTRL_X = 3050` vs 2780) so the "NMD pair" label between the two pair brackets clears the right bracket's vertical segment — caught by the layout validator on the first render attempt.
- `xlim` extended to (−1200, 3600) for the wider bracket span and the left-of-track labels.

## Reproducibility

This script depends on the source panel PNGs existing. To regenerate from scratch:

```bash
cd figures/multipanel/figure3_isopair_and_ptc/
Rscript data_export.R                                     # regenerate TSV inputs
~/miniforge3/bin/python figure3_panelA_pair_concept.py     # render panel A
~/miniforge3/bin/python figure3_panelB_sequence_similarity.py
~/miniforge3/bin/python figure3_panelC_event_prevalence.py
~/miniforge3/bin/python figure3_panelD_stop_codon_distance.py
~/miniforge3/bin/python figure3_panelE_ptc_event_attribution.py
~/miniforge3/bin/python figure3_panelF_mechanism_breakdown.py
~/miniforge3/bin/python figure3_composite.py               # compose
```

The whole chain is byte-deterministic on the current cached RDS inputs (verified in Pass 4 / Pass 6.4 of the 5-pass verification protocol).

## Caveats / known issues

1. **Source PNGs not PDFs.** The composite embeds the rendered PNG of each panel, not its vector PDF. At high zoom (or for print), text in the embedded panels is slightly softer than it would be from vector embedding. To switch to PDF-based embedding would require either using `matplotlib.backends.backend_pgf` with vector handling, or composing in a vector-aware tool (Illustrator, Inkscape) instead. Current PNG embedding is the standard grant-figure pattern.

2. **Per-panel internal titles overlap with composite letter labels.** Each panel renders its own title (e.g., "Sequence shared with reference") above its plot, and the composite adds the letter label (A, B, …) just above the cell. Two title-like elements stack at the top of each cell. If you want a cleaner look, drop the in-panel `set_title()` calls and let the composite letter carry the labeling role. Pete to confirm before changing.

3. **Panel A uses smaller fonts** (15/11 pt vs 18/14 pt in B–F) so the left-of-track descriptive labels fit at 6" wide. If Panel A is rendered standalone (not composite), the smaller fonts look mildly disproportionate. Acceptable because Panel A's primary use-case is composite-embedded.

## Cross-references

- Per-panel methodology: `figure3_panel{A,B,C,D,E,F}_methodology.md`
- Source: `figure3_composite.py`
- Grant convention reference: `~/claude_projects/grants/nmd_2026/figure_nmd_prelim_composite_3by1.py`
- Figure principles: `~/.claude/skills/figures/references/figures_principles.md` and `figures/lib/principles.md`
- Verification: covered by Pass 4 (reproducibility) of the 5-pass scientific-report verification protocol — `verify_pass4_reproducibility.sh` re-runs the entire panel + composite chain end-to-end and confirms byte-identical regeneration.

## Layout validator coverage

`figures/lib/validate_figure_layout.py` is a matplotlib annotation-primitive validator — it walks `Text`, `Line`, `Rectangle`, and `Polygon` objects on a live `(fig, ax)` to flag text-segment collisions, text-text overlap, centering, crowding, and shape-shape partial overlap. It has no CLI; it must be called from Python with `(fig, ax)` before `savefig`.

**Every Figure 3 panel's `main()` now runs the validator and prints its report before saving.** Panel A raises on errors (annotation-heavy schematic — errors are unambiguous); B–F print results but don't raise (data plots produce some validator false positives, see below — visual review is still required).

| Panel | Validator status |
|---|---|
| A | **Clean** — 0 errors, 0 warnings. (8 texts, 64 segments, 17 rects.) Raises on any error. |
| B | **Clean** — 0 errors, 0 warnings. (No callouts to collide.) |
| C | **Clean** — 0 errors, 0 warnings. (Significance stars + bars, no overlap.) |
| D | 2 errors, both **false positive** — `text_segment_collision` between the "50-nt PTC threshold" callout (at data-y ≈ 0.0016) and the x-axis spine at y=0. The y-axis range here is ~0.0017 in absolute density units, which the validator's text-bbox estimator collapses near zero. Visually the callout sits at the top of the plot, well clear of the spine. |
| E | 4 errors, all **false positive** — `text_rect_overflow` for the title and subtitle (positioned in axes-coords above the data area) being compared against the leftmost bar's data-coord rect. The mixed-coordinate comparison gives a spurious overlap. Visually the title/subtitle sit cleanly above the bars. |
| F | 6 errors. **2 are false positive** (same axes-vs-data title overflow as Panel E, against the bottom "Other" bar). **4 are over-conservative** — `text_rect_overflow` reporting that the in-segment value labels (259, 109, 165, 737) vertically overflow their bar segments. At INNER_FS=9pt in a 4-inch canvas, each label is ~38 px tall and each bar is ~63 px tall — the labels visually center inside the bars, but the validator's text-height estimator pads the bbox enough to trip the check. Originally Panel F produced **47 errors** before tuning (oversized value labels, label-stack overlap across rows, legend clipping); the remaining 6 are residual conservatism, not real layout problems. |

**Iteration history**: Panel F's original render had `LABEL_MIN=10` (label any segment ≥10) and `BODY_FS=14` for in-segment labels — at 10×7 figsize that fit, but at 6×4 it produced the "mess" described in the validator output (in-segment labels overlapping each other across rows, intruding into the y-axis label column, and the "Frameshift" legend entry truncated at the right edge). Tuned to `LABEL_MIN=100`, `INNER_FS=9`, `LEGEND_FS=11` with `handlelength=1.2 / handletextpad=0.5 / columnspacing=1.0` for the legend. Title bumped to y=1.20 / subtitle to y=1.04 so their bboxes don't collide.

- **The composite itself** embeds pre-rendered PNGs via `imshow` plus 6 letter labels. Running the validator over each composite axis reports 0 errors / 0 warnings, but this is mostly trivial — the load-bearing signal is whether each contributing panel was visually validated.
