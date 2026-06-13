# Figure 3 Panel D — methodology

**Title:** Stop codon to last EJC
**Render script:** `figure3_panelD_stop_codon_distance.py`
**Data TSV:** `data/panelD_stop_codon_distance.tsv`
**Data export:** `data_export.R` ([Panel D] section)
**Scope:** All data (no test-set filter)

## Headline claim

Distance from each comparator isoform's stop codon to its last exon-exon junction. **Positive = upstream of last EJC = PTC direction; negative = stop in or after last exon.** NMD comparators concentrate to the right of the 50-nt PTC threshold (median far upstream, indicating PTCs); Control comparators concentrate to the left (peak around -100 nt, normal stops in the last exon).

## Source data

| File | What it provides |
|---|---|
| `data_mashr/profiles_c{2,4}_allsamples.rds` | Pair-level metadata for Stage 2 gene-matching |
| `data_mashr/structures.rds` | Per-isoform exon coordinates (`exon_starts`, `exon_ends`, `strand`, `n_exons`) — used to compute last-junction transcript position |
| `data_mashr/analysis_cache/ref_atg_analysis.rds` | Reference-AUG-traced stop position (`comp_stop_tx_pos`) per comparator |

## Population

**Denominator: 2,289 NMD + 1,763 Control** = `pop_traceable` in `data_export.R`. Pairs where ref-AUG tracing produced a valid ORF stop position — categories: `effectively_ptc` + `no_downstream_ejc` + `truncated_no_ejc`.

**Why this scope (not pop_BC at 3,009):** Per Pete's policy 2026-06-13: "TD2 CDS annotations are unreliable, so for all analyses that depend on identifying the stop codon, we limit to isoform pairs where reference AUG tracing can be performed." Distance to last EJC requires a defined stop position; we restrict to pairs where ref-AUG tracing gives one.

**Pairs excluded** (NMD: 720 / Control: 1,246):
- `ref_atg_lost`: ref AUG not exonic in comparator (351 NMD / 770 Control)
- `mapping_failed`: tracing pipeline errored (33 NMD / 46 Control)
- Pairs never traced (336 NMD / 430 Control): the 05r_ref_atg_analysis.R pipeline was pre-filtered to coding-coding pairs, so these never went through tracing at all. See methodology Caveat #4 for the consequence.

## Computation

```text
For each comparator in pop_traceable:
  stop position    = ref_atg_analysis::comp_stop_tx_pos   (ref-AUG-traced)
  last EJC pos     = sum(exon_lengths[1:(n_exons-1)]) in transcript order
  distance         = last_ejc_tx_pos - comp_stop_tx_pos
                     # positive = stop UPSTREAM of last EJC = PTC direction
                     # convention matches ptc.rds$ptc_distance
```

Last-EJC transcript position from `structures.rds`:
- + strand: `sum(exon_lengths[1:(n_exons-1)])` in genomic order
- − strand: `sum(rev(exon_lengths)[1:(n_exons-1)])` (transcript runs from the genomic end)

KDE density plot via `scipy.stats.gaussian_kde` with Scott's bandwidth. 50-nt PTC threshold marker drawn at x = 50.

## Output TSV columns

| Column | Type | Notes |
|---|---|---|
| `comparator_isoform_id` | string | |
| `comparison` | "NMD" \| "Control" | |
| `distance` | float | Nucleotides; positive = upstream of last EJC (PTC direction) |
| `category` | string | One of `effectively_ptc` / `no_downstream_ejc` / `truncated_no_ejc` |

## Headline numbers under this scope

- NMD PTC rate (fraction with distance > 50 nt): **83.5%** (= 1,912 of 2,289)
- Control PTC rate: **16.3%** (= 288 of 1,763)
- Aligns with manuscript prose claim "≥85% explained by PTCs" (after ref-AUG tracing)

## Caveats / limitations

1. **Axis clipping**: X-axis is `[-1000, 1500]` nt for legibility. Clipped at this render: NMD = 227 (10%), Control = 89 (5%). Clipped tail values are at very long distances where density is near zero.
2. **KDE smoothing**: `scipy.stats.gaussian_kde` with Scott's bandwidth. Some smoothing across the 50-nt threshold is unavoidable.
3. **The 336 NMD and 430 Control pairs never subjected to ref-AUG tracing** are excluded because the 05r pipeline pre-filtered to coding-coding pairs. To include them would require extending `05r_ref_atg_analysis.R` to run on ALL gene-matched pairs (including those where the comparator was `coding_status == "unknown"` in the initial SQANTI3+TD2 analysis). Task #23 covers this upstream extension.
4. **Test-only sensitivity check** not generated; primary analysis uses all data per project policy.

## Cross-references

- `figures/lib/principles.md` — figure-making principles
- `feedback_figure_sample_size_consistency` — denominator differs from B/C (3,009) and from E/F (1,912/288), but each deviation is justified by the underlying analytical need (Panel D needs stop positions; Panel E/F further restrict to PTC+)
- `feedback_sqanti_cds_ptc_bias` project memory — TD2's anti-PTC bias (novel isoforms) motivates ref-AUG tracing as the canonical CDS analysis
- `feedback_default_match_original_figure` — Panel D's distance-density concept matches the original Rmd's `goal1-fig1-stop-dist` plot at line 1717
- Rmd source: `05_final_report_mashr.Rmd` chunks `goal1-fig1-stop-dist` (line ~1717), `sec2c-ref-atg-load` (line 3268); driver `05r_ref_atg_analysis.R`
- Replaces the previous (AM revision) Panel D, which was a bar chart of "TD2 → ref-AUG combined PTC rate" — that figure is now moved to a supplement (`figure3_supp_ptc_rate_by_cds.*`)

## Composite cross-reference

This panel appears in the Figure 3 composite as the labeled cell — see `figure3_composite_methodology.md` for the layout and per-cell mapping. The composite embeds the panel's pre-rendered PNG (`figure3_panelD.png`); regenerating the composite does NOT re-run this panel. Re-render this panel script first if its data or rendering has changed.
