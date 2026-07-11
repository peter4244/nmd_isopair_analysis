# Figure 3 Panel D — methodology

**Title:** Stop codon to last EJC
**Render script:** `figure3_panelD_stop_codon_distance.py`
**Data TSV:** `data/panelD_stop_codon_distance.tsv`
**Data export:** `data_export.R` ([Panel D] section)
**Scope:** All data (no test-set filter)

## Headline claim

Distance from each comparator isoform's stop codon to its last exon-exon junction. **Positive = upstream of last EJC = PTC direction; negative = stop in or after last exon.** NMD comparators have a long right-tail extending past the 50-nt PTC threshold (36.9% PTC+, median −66 nt); Control comparators concentrate near the last-exon stop position (median −143 nt, 1.5% PTC+) — an 24-fold enrichment of PTCs in NMD substrates.

## Source data

| File | What it provides |
|---|---|
| `data_mashr/profiles_c{2,4}_allsamples.rds` | Pair-level metadata for Stage 2 gene-matching |
| `data_mashr/structures.rds` | Per-isoform exon coordinates (`exon_starts`, `exon_ends`, `strand`, `n_exons`) — used to compute last-junction transcript position |
| `data_mashr/analysis_cache/ref_atg_analysis.rds` | Reference-AUG-traced stop position (`comp_stop_tx_pos`) per comparator |

## Population

**Denominator: 130 NMD + 130 Control** = pop_BC ∩ all-3-ENST ∩ all-3-coding-CDS in `data_export.R`. Gene-matched pairs where **all three isoforms** (reference, NMD comparator, Control comparator) are GENCODE-annotated (ENST-prefixed) AND have CDS annotation (`coding_status == "coding"` in `cds.rds`).

**Stop position used:** each comparator's OWN GENCODE-annotated stop codon (computed from `cds.rds` `cds_start` / `cds_stop` mapped to transcript coordinates via `Isopair::genomicToTranscript`). No ref-AUG projection; no TD2; no `Isopair::traceReferenceAtg` — every isoform's PTC status is determined from its own curated annotation.

**Why this scope:** Pete's clarified policy 2026-06-15 — when all three isoforms in a pair are GENCODE-annotated and have CDS, the cleanest analysis uses each isoform's own GENCODE-annotated CDS directly. No projection of the reference's start codon onto the comparator is needed (or appropriate). This produces a fully TD2-free, GENCODE-anchored analysis.

**Filter cascade** (relative to pop_BC at 1,548 each):
- Both isoforms ENST per side: NMD 525 / Control 993 (single-side)
- All 3 ENST gene-matched: 301 NMD / 301 Control / 301 unique gene-reference combos
- All 3 ENST + coding-CDS, re-intersected: **130 NMD / 130 Control**

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
| `gene_id` | string | |
| `reference_isoform_id` | string | |
| `comparison` | "NMD" \| "Control" | |
| `distance` | float | Nucleotides; `last_ejc_tx_pos − own_stop_tx_pos`; positive = upstream of last EJC (PTC direction) |
| `last_ejc_tx_pos` | integer | Transcript-coordinate position of the last exon–exon junction |
| `own_stop_tx_pos` | integer | Transcript-coordinate position of the comparator's own GENCODE-annotated stop codon |

## Headline numbers under this scope

- NMD PTC rate (fraction with distance > 50 nt under own GENCODE stop): **36.9%** (= 72 of 190)
- Control PTC rate: **1.5%** (= 4 of 190)
- Fold enrichment: **24×** (36.9 / 1.5)
- Numerically distinct from the ref-AUG-projected PTC rate at the ENST-reference scope (89.8% / 16.1%, 5.6×) because (a) the all-3-ENST scope is more restrictive (1,171 → 130 NMD pairs) and (b) the GENCODE-stop measure doesn't try to project the reference's ATG into the comparator; each isoform's own annotated stop drives the PTC determination.

## Caveats / limitations

1. **Axis clipping**: X-axis is `[-1000, 1500]` nt for legibility. Clipped at this render: **NMD = 4 (1.5%), Control = 4 (1.5%)**. Clipped tail values are at very long distances where density is near zero.
2. **KDE smoothing**: `scipy.stats.gaussian_kde` with Scott's bandwidth. Some smoothing across the 50-nt threshold is unavoidable.
3. **Scope construction**: the all-3-ENST + coding-CDS filter cascade reduces 1,548 pop_BC pairs to 130 per side. Pairs where the reference, NMD comparator, or Control comparator is novel (non-ENST) or non-coding are excluded — they are characterized at the broader ref-AUG-traceable scope in Figure 4 Section C.
4. **Test-only sensitivity check** not generated; primary analysis uses all data per project policy.

## Cross-references

- `figures/lib/principles.md` — figure-making principles
- `feedback_figure_sample_size_consistency` — denominator differs from B/C (1,548) and from E/F (48 PTC+ NMD, 2 PTC+ Control), but each deviation is justified by the analytical need (Panel D needs stop positions; Panel E/F further restrict to PTC+)
- `feedback_sqanti_cds_ptc_bias` project memory — TD2's anti-PTC bias (novel isoforms) motivates ref-AUG tracing as the canonical CDS analysis
- `feedback_default_match_original_figure` — Panel D's distance-density concept matches the original Rmd's `goal1-fig1-stop-dist` plot at line 1717
- Rmd source: `05_final_report_mashr.Rmd` chunks `goal1-fig1-stop-dist` (line ~1717), `sec2c-ref-atg-load` (line 3268); driver `05r_ref_atg_analysis.R`
- Replaces the previous (AM revision) Panel D, which was a bar chart of "TD2 → ref-AUG combined PTC rate" — that figure is now moved to a supplement (`figure3_supp_ptc_rate_by_cds.*`)

## Composite cross-reference

This panel appears in the Figure 3 composite as the labeled cell — see `figure3_composite_methodology.md` for the layout and per-cell mapping. The composite embeds the panel's pre-rendered PNG (`figure3_panelD.png`); regenerating the composite does NOT re-run this panel. Re-render this panel script first if its data or rendering has changed.
