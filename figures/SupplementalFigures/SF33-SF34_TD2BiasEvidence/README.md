# Supplemental Figure — TD2 ORF-call bias on the §4 universe (1:1 gene-matched, ENST-reference + ref-AUG-traceable)

Consolidated supplemental figure presenting three orthogonal observations of TransDecoder2 (TD2)'s ORF-call behavior at two complementary scopes within the **same isoform universe as Figure 4 Panels C/D** (1:1 gene-matched NMD and Control comparator pairs, re-intersected on `(gene_id, reference_isoform_id)` after the ref-AUG-traceable category filter). 2×3 layout:

- **Row 1 (Panels A–C; n = 1,166 broad scope)** — Figure 4 Section C universe: ENST-reference NMD comparator pairs in the ref-AUG-traceable category set, with TD2 CDS available and reference AUG exonic in the comparator, restricted to gene/reference keys that also have a matched Control comparator after the same category filter. TD2 agreed with the reference AUG in 50% of pairs (578) and picked a different ATG in 50% (588).
- **Row 2 (Panels D–F; n = 492 occult-PTC subset)** — the subset of Row 1 in which reference-AUG tracing revealed a downstream-EJC stop that TD2's own CDS call missed (`effectively_ptc ∩ original_ptc == FALSE`, within the re-intersected universe). TD2 ATG ≠ ref AUG by construction; the bias is most extreme.

## Files

- `data_export.R` — recomputes all six panels' input data from upstream RDS caches (`ref_atg_analysis.rds`, `cds.rds`, `structures.rds`, `profiles_c{2,4}_allsamples.rds`, and the comparator FASTA). Run with `/usr/local/bin/Rscript data_export.R`. Writes:
  - `data/panelA_td2_vs_refaug_length.tsv`, `panelB_kozak.tsv`, `panelC_td2_position.tsv` — broad scope (n=1,166)
  - `data/panelD_td2_vs_refaug_length_occult.tsv`, `panelE_kozak_occult.tsv`, `panelF_td2_position_occult.tsv` — occult-PTC subset (n=492)
  - `data/summary.tsv` — one-line summary of headline numbers per panel
- `figure_s_td2_bias.py` — renders the 2×3 composite (18″×10″ landscape).
- `figure_s_td2_bias.{pdf,png}` — rendered output.
- `figure_s_td2_bias_legend.md` — composite legend.
- `legacy_main_panelC.py`, `legacy_main_panelD.py` — preserved copies of the prior main-figure panel scripts. Reference only.

## Scope (top vs bottom row)

| | Broad (Row 1) | Occult-PTC (Row 2) |
|---|---|---|
| Filter | pop_BC NMD c2 ∩ ENST-ref ∩ ref-AUG-traceable category ∩ TD2 CDS ∩ ref-AUG ORF computable ∩ ref AUG exonic ∩ re-intersected with Control side | Row 1 ∩ `category == "effectively_ptc"` ∩ `original_ptc == FALSE` |
| **n** | **1,166** | **492** |
| TD2 == ref AUG | 578 (50%) | 0 (0%) |
| TD2 ≠ ref AUG | 588 (50%) | 492 (100%) |

This scope **matches Figure 4 Panels C/D exactly** (1,166 NMD comparator pairs paired 1:1 with 1,166 Control comparator pairs at the same gene/reference keys).

## Panel layout

| | Length | Kozak | Position |
|---|---|---|---|
| **Broad (n=1,166)** | A: KDE | B: paired violin | C: histogram |
| **Occult-PTC (n=492)** | D: KDE | E: paired violin | F: histogram |

## Headline numbers per panel

| Panel | n | Metric | Result |
|---|---:|---|---|
| A | 1,166 | TD2 / ref-AUG ORF length ratio | 4.69× (median 2,758 / 588 nt) |
| B | 1,166 | Paired Kozak (ref vs TD2) | p = 1.2×10⁻³⁸; medians 0.86 vs 0.29; ref Kozak > TD2 Kozak in 446/1,166 (38.3%) |
| C | 1,166 | TD2 position vs ref AUG | 50% same; 43% downstream (+471 nt median); 8% upstream |
| D | 492 | TD2 / ref-AUG ORF length ratio | **7.73×** (median 2,934 / 380 nt) |
| E | 492 | Paired Kozak (ref vs TD2) | **p = 2.6×10⁻³⁸**; medians 1.04 vs −0.35; ref Kozak > TD2 Kozak in **384/492 (78.0%)** |
| F | 492 | TD2 position vs ref AUG | **99.0% downstream** (median +476 nt); 1.0% upstream; 0% same |

The cross-scope comparison shows the bias scales with the bias-relevant population: in the occult-PTC subset (where TD2 had to avoid a real PTC), the ORF-length ratio nearly doubles, the directional Kozak agreement rises from 38% to 78%, and the position spike at 0 disappears entirely.

## Reproducibility

```bash
cd .../figures/SupplementalFigures/TD2BiasEvidence
/usr/local/bin/Rscript data_export.R                              # writes data/ TSVs (broad + occult-PTC)
/Users/petecastaldi/miniforge3/bin/python figure_s_td2_bias.py    # → .pdf, .png
```
