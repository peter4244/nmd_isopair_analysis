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

Skipped exon (SE) dominates with 23 attributed pairs; A5SS (6) and A3SS (5) follow. Total = **48 attributed PTC+ pairs**. The mechanism breakdown across all events is **25 Frameshift / 18 In-frame stop / 5 3'UTR splice = 52.1% / 37.5% / 10.4%** — frameshift is the leading mechanism but in-frame stop events contribute a substantial minority.

## Source data

| File | What it provides |
|---|---|
| `all_attr_new` (computed in `panel_e_compute.R`) | Per-pair PTC attribution combining ref_atg_analysis$c2$attr_event/attr_mechanism (for 900 ref-AUG-recovered) + Rmd's `attribute_ptc_events()` output (for 1,012 TD2-PTC+) |
| `analysis_functions.R::shorten_event_labels` | Collapses rare event types into more readable labels (`event_short`) |

## Population

**Denominator: 69 attributed PTC+ pairs** in NMD, derived from the all-3-ENST + coding-CDS scope:

1. Stage 2 gene-matched NMD (1,548)
2. All 3 isoforms ENST gene-matched (301 NMD / 301 Control)
3. All 3 isoforms ENST + coding-CDS, re-intersected (**190 NMD / 190 Control**)
4. PTC+ subset (own GENCODE stop > 50 nt past last EJC): **72 NMD / 4 Control**
5. Attribution chain (own GENCODE stop input to `attribute_ptc_events` / `attribute_3utr_splice`): **69 directly attributed**

The 3 PTC+ pairs not in `all_attr_new` (72 − 69) are unresolved at the attribution step. They're absent from Panel F by definition.

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

1. **Denominator under all-3-ENST + coding-CDS scope (2026-06-15):** the original Rmd `fig_sankey_ptcpos` had 1,074 attributed pairs (Stage 2 coding-coding, TD2-PTC+ only). Iterations went through 1,812 (mixed source) and 1,434 (§4 ENST-reference + ref-AUG). The current Panel F at the all-3-ENST + coding-CDS scope yields **48 attributed pairs** — fully GENCODE-anchored, no TD2 anywhere, no ref-AUG projection.
2. **Mechanism balance at this scope:** **25 Frameshift / 18 In-frame stop / 5 3'UTR splice** (= 52.1% / 37.5% / 10.4%). Frameshift is the most common but the breakdown is less extreme than at broader scopes (where Frameshift dominated 76%).
3. **The 3 unresolved PTC+ pairs** (72 − 69) are PTC+ by GENCODE-stop classification but no specific splice event could be attributed. They are NOT counted in Panel F's bars but DO contribute to Panel D's distance density.
4. **Test-only sensitivity check** not generated; primary analysis uses all data per project policy.

## Cross-references

- `figures/lib/principles.md` — figure-making principles
- `feedback_figure_sample_size_consistency` — denominator matches Panel E (48 attributed events at all-3-ENST + coding-CDS scope) with the larger Panel D population (190 NMD / 190 Control) well-documented
- `feedback_default_match_original_figure` — Panel F's horizontal stacked bar structure matches the original `fig_sankey_ptcpos` from the Rmd
- Rmd source: `05_final_report_mashr.Rmd` chunks `goal2-ptc-mechanisms` (line ~1923), `goal2-fig5-sankey` (line ~2002)
- `analysis_functions.R`: `attribute_ptc_events()`, `attribute_3utr_splice()`, `shorten_event_labels()`
- Original render preserved: `nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/figures/fig_sankey_ptcpos.{pdf,png}`

## Composite cross-reference

This panel appears in the Figure 3 composite as the labeled cell — see `figure3_composite_methodology.md` for the layout and per-cell mapping. The composite embeds the panel's pre-rendered PNG (`figure3_panelF.png`); regenerating the composite does NOT re-run this panel. Re-render this panel script first if its data or rendering has changed.
