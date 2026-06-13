# Figure 3 Panel E — methodology

**Title:** PTC-causing event attribution
**Render script:** `figure3_panelE_ptc_event_attribution.py`
**Data TSV:** `data/panelE_ptc_event_attribution.tsv`
**Data export:** `data_export.R` ([Panel E+F] section) → sources `panel_e_compute.R`
**Scope:** All data (no test-set filter)

## Headline claim

Among splice events Isopair attributes as the *cause* of a PTC in NMD comparators (under ref-AUG-traced stop positions), **skipped exon (SE) accounts for ~55% of attributed PTC-causing events**, vs. ~9% among all Control events. A3SS (~9%) and A5SS (~7%) are also significantly enriched. Terminal events (Alt_TSS, Alt_TES) are *depleted* among PTC-causing events because they don't typically introduce PTCs in the coding region.

## Source data

| File | What it provides |
|---|---|
| `data_mashr/profiles_c{2,4}_allsamples.rds` | NMD and Control pair profiles |
| `data_mashr/cds.rds` | Coding status + CDS start/stop + strand lookups |
| `data_mashr/ptc.rds` | Per-isoform PTC call (`has_ptc`) — used to split effectively_ptc into TD2-PTC+ and recovered subsets |
| `data_mashr/analysis_cache/fc_c2_allsamples.rds` | Frame comparison (`pair_summary` with `frame_category`) |
| `data_mashr/analysis_cache/fw_c2_allsamples.rds` | Frame walk (`events`) |
| `data_mashr/analysis_cache/ref_atg_analysis.rds` | Ref-AUG-traced ORF metadata (`category`, `attr_event`, `attr_mechanism`) |
| `analysis_functions.R` | `attribute_ptc_events()`, `attribute_3utr_splice()`, `build_cds_lookups()`, `shorten_event_labels()` |

## Population

**Denominator (event-level): n_ptc_attr = 1,812 attributed PTC-causing NMD events; ctrl_total = 4,525 total events in pop_trace_c4 (1,763 Control pairs).**

**Pair-level population: pop_ptc_plus = 1,912 NMD effectively_ptc pairs.** All NMD pairs in this panel are PTC+ under the ref-AUG-traced ORF definition.

## Attribution chain (the new scope's mixed-source attribution)

Pairs in `pop_ptc_plus` (1,912 NMD effectively_ptc) split into two attribution sources based on original (TD2) PTC status:

**(a) Original_ptc = FALSE (900 pairs — ref-AUG-recovered):**
Attribution stored directly in `ref_atg_analysis$c2$attr_event` / `attr_mechanism` columns (computed by `05r_ref_atg_analysis.R` at ref-AUG runtime, using the ref-AUG-traced stop position).

- 845 directly attributed
- 55 unresolved (PTC+ but no specific event could be assigned)

**(b) Original_ptc = TRUE (1,012 pairs — TD2 already called PTC+):**
Re-runs `attribute_ptc_events()` from `analysis_functions.R` using **TD2's stop position** (since these pairs have an original TD2 CDS call). Plus `attribute_3utr_splice()` for same-stop subset.

- 928 diff-stop directly attributed + 79 same-stop 3'UTR-splice attributions
- Some overlap removed in deduplication

**Combined: `all_attr_new`:**
- 1,812 directly attributed PTC+ pairs (= 845 + 928 + 79 − overlap)
- Attribution rate: 1,812 / 1,912 = 94.8%

**Mixing of stop-position sources is intentional but worth noting:** the 900 ref-AUG-recovered pairs use ref-AUG stops for attribution; the 1,012 TD2-PTC+ pairs use TD2 stops (the original Rmd chain). Methodologically cleaner would be to ALSO re-run attribute_ptc_events with ref-AUG stops for the TD2-PTC+ subset, but the existing all_attr chain is preserved here for consistency with the published Rmd analysis pipeline. Task #23 covers the upstream cleanup.

**Control side:** no per-pair mechanism attribution. The Control baseline is the flattened event distribution across `pop_trace_c4`'s 1,763 pair detailed_events (4,525 events total).

## Computation

| Number | Computation |
|---|---|
| `n_ptc_events` | count of PTC-causing events of this type within `all_attr_new` |
| `pct_of_ptc` | `100 * n_ptc_events / 1812`, 1 decimal |
| `n_ctrl_events` | count of events of this type in flattened pop_trace_c4 detailed_events |
| `pct_ctrl` | `100 * n_ctrl_events / 4525`, 1 decimal |
| `enrichment` | `pct_of_ptc / max(pct_ctrl, 0.1)`, 1 decimal |
| `fisher_OR`, `fisher_p` | Fisher's exact on 2×2 (this-event-in-PTC, n_ptc_attr-this, this-event-in-Ctrl, ctrl_total-this) |

11 event types total (10 base + IR_diff_5/3 split + "Other" if shorten_event_labels collapsed any).

## Output TSV columns

| Column | Type | Notes |
|---|---|---|
| `event_type` | string | Shortened event label |
| `n_ptc_events`, `pct_of_ptc` | int / float | NMD PTC-causing event counts and % (of 1,812) |
| `n_ctrl_events`, `pct_ctrl` | int / float | Control all-event counts and % (of 4,525) |
| `enrichment` | float | NMD/Control proportion ratio |
| `fisher_p` | float | Fisher exact p-value |
| `direction` | string | "Enriched in PTC-causing" / "Depleted in PTC-causing" / "Non-significant" |

## Caveats / limitations

1. **Larger denominators than original Rmd:** original `fig5b_ptc_event_proportions` had `n_ptc_attr` = 985 (under the broader Stage 2 + coding-coding scope but missing ref-AUG-recovered pairs). The new 1,812 includes both the original TD2 PTC+ attributions and the 900 ref-AUG-recovered with their ref_atg_analysis-attributed events. The relative event proportions are roughly preserved but absolute counts increase by ~85%.
2. **Mixed stop-position sources** (see Attribution chain section above) — TD2 stops for the TD2-PTC+ subset, ref-AUG stops for the recovered subset. Methodologically suboptimal but matches the existing Rmd pipeline. Task #23 covers the canonicalization.
3. **`mechanism != "3'UTR splice"` is NOT enforced here.** Unlike the prior Panel E iteration, this version includes 3'UTR-splice attributions in the event-type breakdown.
4. **Control baseline includes all events**, not just events of type that could plausibly cause a PTC. This matches the original Rmd convention.
5. **Test-only sensitivity check** not generated; primary analysis uses all data per project policy.

## Cross-references

- Rmd source: `05_final_report_mashr.Rmd` chunks `goal2-ptc-mechanisms` (line ~1923), `goal2-fig5-sankey` (line ~2002), `goal2-table4b-baseline` (line ~2179)
- `analysis_functions.R`: `attribute_ptc_events()` (line 41), `attribute_3utr_splice()` (line 177), `build_cds_lookups()`, `shorten_event_labels()`
- `feedback_figure_sample_size_consistency` — denominator differs from B/C/D because it's event-level, and the underlying pair set differs because Panel E/F restrict to PTC+

## Composite cross-reference

This panel appears in the Figure 3 composite as the labeled cell — see `figure3_composite_methodology.md` for the layout and per-cell mapping. The composite embeds the panel's pre-rendered PNG (`figure3_panelE.png`); regenerating the composite does NOT re-run this panel. Re-render this panel script first if its data or rendering has changed.
