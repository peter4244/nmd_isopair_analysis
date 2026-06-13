# Figure 3 Panel C — methodology

**Title:** Splice event prevalence
**Render script:** `figure3_panelC_event_prevalence.py`
**Data TSV:** `data/panelC_event_prevalence.tsv`
**Data export:** `data_export.R` ([Panel C] section)
**Scope:** All data (no test-set filter)

## Headline claim

Across 10 splice-event categories detected by Isopair, **skipped exon (SE)** is enriched in NMD pairs vs. Control pairs by roughly 2× (~44% NMD vs. ~21% Control, Fisher p ≈ 10⁻⁸¹). A5SS, A3SS, and Missing_Internal are also significantly enriched in NMD. Intron retention is *depleted* in NMD relative to Control. Terminal events (Alt_TSS, Alt_TES) are highly prevalent in both groups.

## Source data

| File | What it provides |
|---|---|
| `data_mashr/profiles_c{2,4}_allsamples.rds` | Per-pair splice-event counts and `detailed_events` list column |

## Population

**Denominator: 3,009 pairs per side** (Stage 2 gene-matched, `pop_BC` in `data_export.R`).

Splice-event analysis does not require CDS information — events apply to all pairs regardless of coding status. So this panel uses the broadest meaningful population.

## Computation

Re-implementation of the Rmd's `compute_event_freq()` function on `pop_BC` (3,009 each), plus Fisher's exact tests per event type with 2×2 contingency:

```text
For each event type e:
  a = NMD pairs with event e
  b = 3,009 - a
  c = Control pairs with event e
  d = 3,009 - c
  Fisher's exact on matrix(c(a,b,c,d), nrow=2) → OR, p-value
```

Event types (10): `Alt_TSS`, `Alt_TES`, `A5SS`, `A3SS`, `SE`, `Missing_Internal`, `IR`, `IR_diff`, `Partial_IR_5`, `Partial_IR_3`.

## Output TSV columns

| Column | Type | Notes |
|---|---|---|
| `event_type` | string | Isopair canonical event-type label |
| `n_pairs_with_event_NMD` | int | NMD pairs with ≥1 event of this type |
| `pct_of_pairs_NMD` | float | `n / 3009 * 100`, 1 decimal |
| `n_pairs_with_event_Ctrl` | int | Same for Control |
| `pct_of_pairs_Ctrl` | float | Same |
| `fisher_OR` | float | Odds ratio, 2 decimals |
| `fisher_p` | float | Fisher exact p-value |
| `direction` | "Enriched in NMD" \| "Enriched in Control" \| "Non-significant" | At p < 0.05 |

## Caveats / limitations

1. **Different denominator from D/E/F.** Same justification as Panel B — event prevalence does not depend on CDS.
2. Pairs can contribute to multiple event-type counts (a pair with both SE and A3SS counts toward both).
3. The 5'/3' split of Partial_IR comes from the `detailed_events` list column. `IR_diff` is aggregated in this panel (no 5'/3' split).
4. Test-only sensitivity check not generated; primary analysis uses all data per project policy.

## Cross-references

- Rmd source: `05_final_report_mashr.Rmd` chunk `goal1-table1b-event-enrichment` (line ~1338), Rmd helper `compute_event_freq()` (line ~85)
- Isopair package: `R/event-detection.R` (upstream event-detection logic)
- `feedback_figure_sample_size_consistency` (denominator difference vs D/E/F justified)
