# Figure 3 Panel F — methodology

**Title:** PTC-causing splice events (stacked by mechanism)
**Render script:** `figure3_panelF_mechanism_breakdown.py`
**Data TSV:** `data/panelF_mechanism_breakdown.tsv`
**Data export:** `data_export.R` → `panel_e_compute.R` (Panel F block)
**Scope:** All data (no test-set filter)
**Mirrors:** Original Figure 4 Panel D ("PTC-Causing Splice Events (PTC+ Pairs)" / Rmd render `fig_sankey_ptcpos`) — same horizontal stacked bar structure with event types on Y and mechanism color-stack.

## Headline claim

For each splice event Isopair attributes as the cause of a PTC in an NMD comparator, the count of pairs is broken down by **mechanism**:

- **Frameshift** (coral) — splicing changes the reading frame; the new frame hits a stop codon
- **In-frame stop** (blue) — splicing keeps the original frame but introduces a premature stop within an event
- **3'UTR splice** (teal) — comparator and reference share the same stop, but 3'UTR splicing positions a new exon-exon junction >50 nt downstream of the stop

Skipped exon (SE) dominates with **~1,000 attributed pairs** — 737 Frameshift + 259 In-frame stop + a small 3'UTR-splice contribution. A3SS, A5SS, Alt TES, and IR are the next most common. The split across all events is approximately Frameshift (~62%), In-frame stop (~33%), 3'UTR splice (~5%) — confirming the manuscript's "roughly evenly split between frameshift and in-frame splicing events" needs an update to reflect the new attribution chain (which now leans more frameshift after including the 900 ref-AUG-recovered pairs, ~80% of which are frameshift-mechanism).

## Source data

| File | What it provides |
|---|---|
| `all_attr_new` (computed in `panel_e_compute.R`) | Per-pair PTC attribution combining ref_atg_analysis$c2$attr_event/attr_mechanism (for 900 ref-AUG-recovered) + Rmd's `attribute_ptc_events()` output (for 1,012 TD2-PTC+) |
| `analysis_functions.R::shorten_event_labels` | Collapses rare event types into more readable labels (`event_short`) |

## Population

**Denominator: 1,812 attributed PTC+ pairs** in NMD, derived from the new ref-AUG-traceable scope:

1. Stage 2 gene-matched NMD (3,009)
2. Ref-AUG-traceable subset (`pop_traceable`, 2,289)
3. Ref-AUG PTC+ subset (`pop_ptc_plus`, 1,912; = effectively_ptc category)
4. Attribution chain (mixed sources, see Panel E methodology): 1,812 attributed PTC+ pairs

The 100 PTC+ pairs not in `all_attr_new` are unresolved at the attribution step (PTC+ by ref-AUG classification but no specific splice event could be attributed). They're absent from Panel F by definition.

## Computation

In `panel_e_compute.R`:

```text
sankey_agg = group_by(all_attr_new, event_short, mechanism) → count
events ordered by total count descending
mechanism order: Frameshift, In-frame stop, 3'UTR splice
```

Rendered as a horizontal stacked bar in matplotlib with `barh()` stacked by mechanism (stack order left-to-right: 3'UTR splice → In-frame stop → Frameshift). In-segment value labels for stacks ≥ 10. Legend at bottom.

## Output TSV columns

| Column | Type | Notes |
|---|---|---|
| `event_type` | string | Shortened event label (e.g., `SE`, `Alt TES`, `A3SS`, `IR`) |
| `mechanism` | "Frameshift" \| "In-frame stop" \| "3'UTR splice" | |
| `n` | int | Pairs attributed to this (event_type, mechanism) cell |

## Caveats / limitations

1. **Denominator difference from the original Rmd render:** the original `fig_sankey_ptcpos` had 1,074 attributed pairs (on the Stage 2 + coding-coding population but missing ref-AUG-recovered pairs). The new 1,812 includes both the original TD2-PTC+ attributions and the 900 ref-AUG-recovered. The relative event ordering is broadly preserved.
2. **Mechanism shift:** because the 900 ref-AUG-recovered pairs are ~80% frameshift-mechanism (per `ref_atg_analysis$c2$attr_mechanism` distribution), the overall mechanism balance shifts toward frameshift relative to the original Rmd render (which was ~45% frameshift / 47% in-frame / 8% 3'UTR splice).
3. **The 100 unresolved PTC+ pairs** (1,912 - 1,812) are PTC+ by ref-AUG classification but no specific splice event could be attributed. They are NOT counted in Panel F's bars but DO contribute to Panel D's distance density (since they have a defined stop position).
4. **Test-only sensitivity check** not generated; primary analysis uses all data per project policy.

## Cross-references

- `figures/lib/principles.md` — figure-making principles
- `feedback_figure_sample_size_consistency` — denominator matches Panel E (1,812 attributed events) with the larger Panel D / E population well-documented
- `feedback_default_match_original_figure` — Panel F's horizontal stacked bar structure matches the original `fig_sankey_ptcpos` from the Rmd
- Rmd source: `05_final_report_mashr.Rmd` chunks `goal2-ptc-mechanisms` (line ~1923), `goal2-fig5-sankey` (line ~2002)
- `analysis_functions.R`: `attribute_ptc_events()`, `attribute_3utr_splice()`, `shorten_event_labels()`
- Original render preserved: `nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/figures/fig_sankey_ptcpos.{pdf,png}`

## Composite cross-reference

This panel appears in the Figure 3 composite as the labeled cell — see `figure3_composite_methodology.md` for the layout and per-cell mapping. The composite embeds the panel's pre-rendered PNG (`figure3_panelF.png`); regenerating the composite does NOT re-run this panel. Re-render this panel script first if its data or rendering has changed.
