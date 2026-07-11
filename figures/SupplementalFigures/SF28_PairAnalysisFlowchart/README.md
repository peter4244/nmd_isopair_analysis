# Supplemental Figure — Pair-analysis cohort flowchart

A Graphviz/DiagrammeR flowchart of the cohort cascade for the manuscript's
pair-based analyses. Backs Figures 3, 4, and their supplements by making the
filter chain from the 4-CT manuscript-scope cohort through pop_BC to each
panel's published n explicit and audit-able in one image.

## Files

- `build_flowchart.R` — standalone R build script. Loads upstream RDS caches
  to pull live counts (filtered isoforms, NMD classification, pair
  construction); inlines the canonical scope numbers (130 / 819 / 348);
  writes both the raw DOT source and a self-contained HTML render.
  Run: `Rscript build_flowchart.R`.
- `figure_s_pair_analysis_flowchart.dot` — canonical Graphviz DOT source.
  This is the reproducible record of the flowchart structure.
- `figure_s_pair_analysis_flowchart.html` — standalone HTML render. Opens
  in any browser; uses Viz.js to render the DOT spec client-side. Use
  Chrome → File → Print → "Save as PDF" to produce a PDF.
- `figure_s_pair_analysis_flowchart_legend.md` — figure legend.

## Layout

```
                 [4-CT manuscript-scope cohort: 26 samples, 123,790 isoforms]
                                       ↓
                            [mashr NMD classification]
                                       ↓
                              [Pair construction]
                                       ↓
                       [pop_BC — 1,548 / 1,548 (HUB)]
                                       ↓
                  ┌────────────────────┴────────────────────┐
                  ↓                                         ↓
        [SECTION A scope, n=130]                  [SECTION C scope, n=819]
        (all-3-ENST + coding-CDS)                 (ENST-ref + ref-AUG-traceable)
        own-GENCODE-stop classifier               mechanism_class_4 classifier
        48 / 82 / 130                              756 / 63 / 819
                  ↓                                         ↓
        ┌─────────┼─────────┐                  ┌────────────┼──────────────┐
        ↓         ↓         ↓                  ↓            ↓              ↓
      Fig 3     Fig 4    CDS/3'UTR          Fig 4        CDS/3'UTR    [TD2BiasEvidence
      D/E/F     A/B      supp Row 1         C/D          supp Row 2    broad n=819]
                                                                             ↓
                                                                       + occult-PTC
                                                                             ↓
                                                                  [TD2BiasEvidence
                                                                   occult-PTC n=348]
```

## Scope

- **Upstream pipeline** (handled in `01_prepare_data_mashr.R`, not in this
  flowchart): raw isocall matrix 645,273 isoforms × 38 samples → outlier
  removal → TMM + 5% filter + `edgeR::filterByExpr` → 123,790 isoforms
  × 36 samples (6 CTs in the matrix).
- **Flowchart starts at**: the **4-CT manuscript-scope cohort = 26 samples**
  (AT 6 · DD 8 · FB 6 · MV 6), the 123,790 filtered isoforms.
- **pop_BC**: Stage-2 gene-matched, 1,548 NMD pairs / 1,548 Control pairs.
- **Section A**: all-3-ENST + coding-CDS + re-intersect → 130 / 130; classifier
  = own-GENCODE-stop + 50-nt rule; groups 72 NMD+/PTC+ · 118 NMD+/PTC− ·
  130 Control. Fisher OR = 37.06, p = 1.58e-14, 24× enrichment.
- **Section C**: ENST-ref + ref-AUG-traceable category + re-intersect →
  819 / 819; classifier = `mechanism_class_4()` on ref-AUG-projected
  ORF; groups 756 NMD+/PTC+ · 63 NMD+/PTC− · 819 Control.
- **TD2BiasEvidence broad**: Section C ∩ TD2 CDS available ∩ ref-AUG ORF
  computable ∩ ref AUG exonic → n = 819.
- **TD2BiasEvidence occult-PTC**: broad ∩ `effectively_ptc` ∩
  `original_ptc == FALSE` → n = 348. TD2/ref ORF 7.67×, ref > TD2 Kozak
  286/348 (82 %, paired p = 9.3e−36), TD2 downstream of ref AUG 99.1 %.

## Number provenance

The flowchart is co-canonical with the inline §1 chunk in
[`05_final_report_gencode_scope_2026-06-15.Rmd`](../../../results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_gencode_scope_2026-06-15.Rmd) §1 "Cohort cascade —
pop_BC through Section A / Section C / TD2 branches". Both consume the same
upstream RDS caches in
`results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr/`. The
canonical scope numbers (130, 819, 348, etc.) are guarded against drift
across the Rmd, this flowchart, and the figure-side TSVs by the verifiers:

- [`reproducibility/verify_pass7_new_rmd.R`](../../../reproducibility/verify_pass7_new_rmd.R) — 37-check expected-value manifest
- [`reproducibility/verify_cross_check_new_rmd_vs_figures.R`](../../../reproducibility/verify_cross_check_new_rmd_vs_figures.R) — 57 cross-checks between Rmd HTML and figure-side TSVs

## Reproducibility

```bash
cd .../figures/SupplementalFigures/PairAnalysisFlowchart
Rscript build_flowchart.R
open figure_s_pair_analysis_flowchart.html   # view in browser
```

For PNG / PDF export there are three options (in increasing fidelity):
- **Quickest**: open the HTML in Chrome and File → Print → "Save as PDF".
- **Online**: paste the `.dot` contents into <https://dreampuf.github.io/GraphvizOnline/> and download as PNG/PDF/SVG.
- **Full reproducibility**: install Graphviz (`brew install graphviz`) +
  `DiagrammeRsvg` + `rsvg` R packages, then
  `Rscript -e "DiagrammeR::export_graph(DiagrammeR::grViz('figure_s_pair_analysis_flowchart.dot'), 'figure_s_pair_analysis_flowchart.png', file_type='png')"`.
