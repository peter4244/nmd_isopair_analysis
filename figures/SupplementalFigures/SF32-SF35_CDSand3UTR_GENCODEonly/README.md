# Supplemental Figure — CDS length and 3′UTR length at two cohort scopes

A 2×3 supplement supporting two §4-related claims at **two cohort scopes**:
(i) CDS length is similar across NMD+/PTC+, NMD+/PTC−, and Control isoforms,
ruling out a CDS-length confound for the 5′UTR / 3′UTR signals; and
(ii) the apparent 3′UTR-length difference between NMD+/PTC+ and Control is
a measurement artifact of the translation-based 3′UTR measure — after bias
correction (walking past the PTC to the next non-PTC stop), PTC+ moves much
closer to Control. The same pattern is shown both at the small GENCODE-only
scope (Row 1) and at the ~6× larger ref-AUG-traceable scope (Row 2).

> Folder name retained as `CDSand3UTR_GENCODEonly` for backward compatibility
> with Rmd / METHODS.md cross-references; the row 2 panels at n=1,166 are
> NOT GENCODE-only (most comparators are novel transcripts at that scope).

## Files

- `data_export.R` — recomputes features from upstream RDS caches (`profiles_c{2,4}_allsamples.rds`, `cds.rds`, `structures.rds`, `ref_atg_analysis.rds`) + sources `figures/lib/mechanism_class.R` for the Section C 3-group classification. Run with `/usr/local/bin/Rscript data_export.R`. Writes:
  - **Row 1 (Section A, n=190)**
    - `data/features_190_subset_long.tsv` — per-isoform CDS / 5′UTR / 3′UTR (both forms)
    - `data/panel{A,B,C,D}_*_descriptives.tsv` — per-group n / median / IQR
    - `data/panel{A,B,C,D}_*_pairwise.tsv` — Wilcoxon HL shift + CI + Cliff's δ + p-value
    - Naming retained from the prior 4-quantity export; the figure consumes the CDS, 3′UTR-translation, and 3′UTR-non-PTC-stop TSVs and ignores the 5′UTR ones because 5′UTR length is shown in main Figure 4 Panel A.
  - **Row 2 (Section C, n=1,166)**
    - `data/features_1166_subset_long.tsv` — per-isoform ref-AUG-projected CDS / 3′UTR
    - `data/panel{D,E,F}_*_refaug_descriptives.tsv` + `data/panel{D,E,F}_*_refaug_pairwise.tsv`
- `figure_s_cds_and_3utr.py` — renders the 2×3 composite (18″×8.5″ landscape).
- `figure_s_cds_and_3utr.{pdf,png}` — rendered output.
- `figure_s_cds_and_3utr_legend.md` — composite legend.

## Panel layout

| Row | Left | Middle | Right |
|---|---|---|---|
| **Row 1 — Section A** (n=190, own-GENCODE-stop classifier) | **A** CDS length | **B** 3′UTR length (translation-based) | **C** 3′UTR length (non-PTC-stop based) |
| **Row 2 — Section C** (n=1,166, ref-AUG-projected classifier) | **D** CDS length (ref-AUG ORF) | **E** 3′UTR length (translation-based) | **F** 3′UTR length (non-PTC-stop based) |

## Scope

- pop_BC (Stage-2 gene-matched): 3,009 / 3,009 NMD / Control pairs
- **Row 1 — Section A** (all-3-ENST gene-matched + coding-CDS + re-intersected on (gene_id, reference_isoform_id)): **190 / 190**; PTC classification = each comparator's own GENCODE-annotated stop, 50-nt rule. Groups: **72 NMD+/PTC+**, **118 NMD+/PTC−**, **190 Control**.
- **Row 2 — Section C** (ENST-reference + ref-AUG-traceable categories + re-intersected on (gene_id, reference_isoform_id) AFTER category filter): **1,166 / 1,166**; classification = `mechanism_class_4()` on ref-AUG-projected categories + ORF-length comparison. Groups: **1,050 NMD+/PTC+**, **116 NMD+/PTC−**, **1,166 Control**. Matches Figure 4 Section C (Panels C/D) scope.

## Methodological note

The two 3′UTR measures are pre-registered in `../../multipanel/figure4_ptcneg_and_model/RATIONALE.md` §4.2:

- **Translation-based** (Panel B): `tx_length(comparator) − own_GENCODE_stop_tx_position`. For NMD+/PTC+ pairs this includes the [PTC → natural stop] coding sequence — upward-biased as a 3′UTR length.
- **Non-PTC-stop based** (Panel C): walk in-frame downstream from the comparator's own GENCODE ATG, skipping any stop with a downstream EJC > 50 nt (i.e., itself a PTC), until reaching the first non-PTC stop; then `tx_length − that_position`. For NMD+/PTC− and Control isoforms (own GENCODE stop is non-PTC by classification) this equals the translation-based measure. For NMD+/PTC+ isoforms it uses the pre-computed `non_ptc_stop_tx_pos` field from `ref_atg_analysis.rds`.

The PTC− vs Control difference under (C) (p = 0.03) replicates — at the much smaller complete-GENCODE-only scope — the broader GENCODE-only 3′UTR finding from the manuscript line 203 paragraph (n = 1,904 NMD vs 22,335 non-NMD outside the Isopair gene-matched set).

## Reproducibility

```bash
cd .../figures/SupplementalFigures/CDSand3UTR_GENCODEonly
/usr/local/bin/Rscript data_export.R                              # writes data/ TSVs
/Users/petecastaldi/miniforge3/bin/python figure_s_cds_and_3utr.py  # → .pdf, .png
```
