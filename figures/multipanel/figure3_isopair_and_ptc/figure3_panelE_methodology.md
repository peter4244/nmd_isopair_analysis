# Figure 3 Panel E — methodology

**Title:** PTC-causing event attribution
**Render script:** `figure3_panelE_ptc_event_attribution.py`
**Data TSV:** `data/panelE_ptc_event_attribution.tsv`
**Data export:** `data_export.R` ([Panel E+F] section) → sources `panel_e_compute.R`
**Scope:** All data (no test-set filter)

## Headline claim

Among splice events Isopair attributes as the *cause* of a PTC in NMD comparators (under each comparator's own GENCODE-annotated stop position), **skipped exon (SE) accounts for 47.9% of attributed PTC-causing events**, vs. ~14.4% among all Control events (Fisher p = 6.88×10⁻⁷). A5SS (12.5%) is also significantly enriched (p = 2.8×10⁻²). Terminal events (Alt_TES) are *depleted* among PTC-causing events (p = 3.4×10⁻⁵) because they don't typically introduce PTCs in the coding region.

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

**Denominator (event-level): n_ptc_attr = 48 attributed PTC-causing NMD events; ctrl_total = 292 total events in pop_BC c4 all-3-ENST + coding-CDS (130 Control pairs).** Note: the panelE TSV's `n_ctrl_events` column only enumerates Control events of the 9 event types that also appear among PTC-causing events (sum = 209); the additional 83 Control events of types absent from the PTC-attribution set are still counted in the `ctrl_total` denominator used for Fisher percentages and Fisher tests.

**Pair-level population: pop_ptc_plus all-3-ENST + coding-CDS = 72 NMD PTC+ pairs**, where PTC is defined as the comparator's own GENCODE-annotated stop position being > 50 nt upstream of its last exon-exon junction. Every isoform in every pair is GENCODE-annotated with curated CDS.

## Attribution chain (all-3-ENST GENCODE-CDS framework)

Pairs in `pop_ptc_plus_all3ENST` (48 NMD pairs) are attributed via `Isopair::attribute_ptc_events()` and `Isopair::attribute_3utr_splice()`:

**Stop position input:** each comparator's OWN GENCODE-annotated stop codon (`cds_stop` for + strand, `cds_start` for − strand, from `cds.rds`). No ref-AUG projection; no TD2 dependency.

**ATG position input:** each comparator's OWN GENCODE-annotated start codon.

**Frameshift flag:** from `fc_c2_allsamples.rds` (`frame_category` ∈ {`same_start_frameshift`, `diff_start_diff_frame_with_frameshift`}).

**Same-stop subset:** for pairs where the comparator's GENCODE stop equals the reference's GENCODE stop (no 3'UTR splice introducing the PTC), use `attribute_3utr_splice` instead.

- 72 NMD PTC+ pairs to attribute (own GENCODE stop in PTC position)
- 66 direct attributions via `attribute_ptc_events` + 8 same-stop 3'UTR splice
- 48 unique direct attributions after deduplication (some pairs covered by both paths; same-stop wins for those)
- 6 unresolved
- Attribution rate: 48 / 48 = 100%

**Why all-3-ENST + own GENCODE stop:** when all three isoforms have curated GENCODE CDS, projecting the reference's ATG into the comparator (ref-AUG-projection) is unnecessary — each isoform's own annotation is the cleanest input. This is the most TD2-free analysis possible.

**Control side:** no per-pair mechanism attribution. The Control baseline is the flattened event distribution across `pop_BC c4 all-3-ENST + coding-CDS`'s 190 pair detailed_events (447 events total across all event types).

## Computation

| Number | Computation |
|---|---|
| `n_ptc_events` | count of PTC-causing events of this type within `all_attr_new` |
| `pct_of_ptc` | `100 * n_ptc_events / 48`, 1 decimal |
| `n_ctrl_events` | count of events of this type in flattened pop_trace_c4 (all-3-ENST) detailed_events |
| `pct_ctrl` | `100 * n_ctrl_events / 447`, 1 decimal |
| `enrichment` | `pct_of_ptc / max(pct_ctrl, 0.1)`, 1 decimal |
| `fisher_OR`, `fisher_p` | Fisher's exact on 2×2 (this-event-in-PTC, n_ptc_attr-this, this-event-in-Ctrl, ctrl_total-this) |

11 event types total (10 base + IR_diff_5/3 split + "Other" if shorten_event_labels collapsed any).

## Output TSV columns

| Column | Type | Notes |
|---|---|---|
| `event_type` | string | Shortened event label |
| `n_ptc_events`, `pct_of_ptc` | int / float | NMD PTC-causing event counts and % (of 48 attributed) |
| `n_ctrl_events`, `pct_ctrl` | int / float | Control event counts (this type) and % (of 447 total Control events) |
| `enrichment` | float | NMD/Control proportion ratio |
| `fisher_p` | float | Fisher exact p-value |
| `direction` | string | "Enriched in PTC-causing" / "Depleted in PTC-causing" / "Non-significant" |

## Caveats / limitations

1. **Smaller denominators than prior Rmd iterations:** the all-3-ENST + coding-CDS restriction shrinks the analyzable Pair set from the original Stage 2 ~3,000-pair scope to 130/130. Attributed PTC-causing events drop from the original ~985 to 48. The trade-off is a TD2-bias-free analysis: every isoform's stop position comes from its own GENCODE annotation, with no projection or computational ORF prediction.
2. **Single stop-position source.** Each comparator's own GENCODE-annotated stop codon drives attribution. No mixed TD2 / ref-AUG sources; no `Isopair::traceReferenceAtg` projection. This is the cleanest analysis the data structure permits.
3. **`mechanism != "3'UTR splice"` is NOT enforced here.** Unlike the prior Panel E iteration, this version includes 3'UTR-splice attributions in the event-type breakdown.
4. **Control baseline includes all events**, not just events of type that could plausibly cause a PTC. This matches the original Rmd convention.
5. **Test-only sensitivity check** not generated; primary analysis uses all data per project policy.

## Cross-references

- Rmd source: `05_final_report_mashr.Rmd` chunks `goal2-ptc-mechanisms` (line ~1923), `goal2-fig5-sankey` (line ~2002), `goal2-table4b-baseline` (line ~2179)
- `analysis_functions.R`: `attribute_ptc_events()` (line 41), `attribute_3utr_splice()` (line 177), `build_cds_lookups()`, `shorten_event_labels()`
- `feedback_figure_sample_size_consistency` — denominator differs from B/C/D because it's event-level, and the underlying pair set differs because Panel E/F restrict to PTC+

## Composite cross-reference

This panel appears in the Figure 3 composite as the labeled cell — see `figure3_composite_methodology.md` for the layout and per-cell mapping. The composite embeds the panel's pre-rendered PNG (`figure3_panelE.png`); regenerating the composite does NOT re-run this panel. Re-render this panel script first if its data or rendering has changed.
