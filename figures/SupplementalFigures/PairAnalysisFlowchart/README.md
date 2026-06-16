# Supplemental Figure — Pair-analysis cohort flowchart

A Graphviz/DiagrammeR flowchart of the cohort cascade for the manuscript's
pair-based analyses. Backs Figures 3, 4, and their supplements by making the
filter chain from the 4-CT manuscript-scope cohort through pop_BC to each
panel's published n explicit and audit-able in one image.

## Files

- `build_flowchart.R` — standalone R build script. Loads upstream RDS caches
  to pull live counts (filtered isoforms, NMD classification, pair
  construction); inlines the canonical scope numbers (190 / 1,166 / 492);
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
                       [pop_BC — 3,009 / 3,009 (HUB)]
                                       ↓
                  ┌────────────────────┴────────────────────┐
                  ↓                                         ↓
        [SECTION A scope, n=190]                  [SECTION C scope, n=1,166]
        (all-3-ENST + coding-CDS)                 (ENST-ref + ref-AUG-traceable)
        own-GENCODE-stop classifier               mechanism_class_4 classifier
        72 / 118 / 190                            1,050 / 116 / 1,166
                  ↓                                         ↓
        ┌─────────┼─────────┐                  ┌────────────┼──────────────┐
        ↓         ↓         ↓                  ↓            ↓              ↓
      Fig 3     Fig 4    CDS/3'UTR          Fig 4        CDS/3'UTR    [TD2BiasEvidence
      D/E/F     A/B      supp Row 1         C/D          supp Row 2    broad n=1,166]
                                                                             ↓
                                                                       + occult-PTC
                                                                             ↓
                                                                  [TD2BiasEvidence
                                                                   occult-PTC n=492]
```

## Scope

- **Upstream pipeline** (handled in `01_prepare_data_mashr.R`, not in this
  flowchart): raw isocall matrix 645,273 isoforms × 38 samples → outlier
  removal → TMM + 5% filter + `edgeR::filterByExpr` → 123,790 isoforms
  × 36 samples (6 CTs in the matrix).
- **Flowchart starts at**: the **4-CT manuscript-scope cohort = 26 samples**
  (AT 6 · DD 8 · FB 6 · MV 6), the 123,790 filtered isoforms.
- **pop_BC**: Stage-2 gene-matched, 3,009 NMD pairs / 3,009 Control pairs.
- **Section A**: all-3-ENST + coding-CDS + re-intersect → 190 / 190; classifier
  = own-GENCODE-stop + 50-nt rule; groups 72 NMD+/PTC+ · 118 NMD+/PTC− ·
  190 Control. Fisher OR = 28.2, p = 1.88e−20, 18× enrichment.
- **Section C**: ENST-ref + ref-AUG-traceable category + re-intersect →
  1,166 / 1,166; classifier = `mechanism_class_4()` on ref-AUG-projected
  ORF; groups 1,050 NMD+/PTC+ · 116 NMD+/PTC− · 1,166 Control.
- **TD2BiasEvidence broad**: Section C ∩ TD2 CDS available ∩ ref-AUG ORF
  computable ∩ ref AUG exonic → n = 1,166.
- **TD2BiasEvidence occult-PTC**: broad ∩ `effectively_ptc` ∩
  `original_ptc == FALSE` → n = 492. TD2/ref ORF 7.73×, ref > TD2 Kozak
  384/492 (78 %, paired p = 2.6e−38), TD2 downstream of ref AUG 99.0 %.

## Number provenance

The flowchart is co-canonical with the inline §1 chunk in
[`05_final_report_gencode_scope_2026-06-15.Rmd`](../../../results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_gencode_scope_2026-06-15.Rmd) §1 "Cohort cascade —
pop_BC through Section A / Section C / TD2 branches". Both consume the same
upstream RDS caches in
`results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr/`. The
canonical scope numbers (190, 1,166, 492, etc.) are guarded against drift
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
