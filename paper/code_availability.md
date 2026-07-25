<!--
Code Availability — draft for the NMD long-read manuscript. Prepared 2026-07-24.
DOIs are Zenodo CONCEPT (all-versions) DOIs (Pete's preference: resolve to latest);
version used in this study noted in brackets. Provenance (verified via Zenodo API):

  repo                      concept DOI                     version DOI (used)
  nmd_isopair_analysis      10.5281/zenodo.21539734         10.5281/zenodo.21539735  [v1.0.0]
  Isopair                   10.5281/zenodo.21536494         10.5281/zenodo.21536495  [v1.0.0]
  Isocall_v1                10.5281/zenodo.21536485         10.5281/zenodo.21536486  [v1.0.0]
  NMD_orf_model_v5_4ct      10.5281/zenodo.21536501         10.5281/zenodo.21539601  [v2.0.0]
                                                            (v1.0.0 = ...21536502, history)
  data                      GEO GSE329233

Before submission: (1) fix minimal Zenodo author metadata on each record (add full
author list + ORCIDs, with Yul); (2) confirm main-repo public name is nmd_isopair_analysis;
(3) in the Google Doc, make "Figure 5" a live cross-reference field.
-->

> ## ⚠️ DO NOT PASTE YET — the code half is superseded
>
> The **data** sentence is final: source data DOI **10.5281/zenodo.21544337** (reserved) and
> GEO **GSE329233**.
>
> The **code** text below still describes four separate repositories. That structure is
> superseded — §1–§5 code is being consolidated into a single repo
> (`nmd_lung_longread_2026`), which will get its own DOI. When that DOI exists, this paragraph
> should cite **3 code DOIs + 1 data DOI**: the consolidated analysis repo, Isopair, Isocall,
> and the Zenodo source data — not the four listed below. `NMD_orf_model_v5_4ct` and
> `nmd_isopair_analysis` become historical records.
>
> See `REPRODUCIBILITY_PLAN.md` decisions 9–11.

**Code availability**

All code required to reproduce the analyses is publicly available on GitHub and permanently archived on Zenodo; the DOIs below are concept (all-versions) DOIs that resolve to the latest release, with the version used in this study noted in brackets. Primary data processing and figure generation are implemented in `nmd_isopair_analysis` (https://github.com/peter4244/nmd_isopair_analysis; DOI: 10.5281/zenodo.21539734 [v1.0.0]), which draws on three supporting repositories: **Isopair**, the gene-matched isoform-pair analysis package for splice-event detection and premature-termination-codon attribution (https://github.com/peter4244/Isopair; DOI: 10.5281/zenodo.21536494 [v1.0.0]); **Isocall**, the long-read (PacBio Iso-Seq) processing pipeline (https://github.com/peter4244/Isocall_v1; DOI: 10.5281/zenodo.21536485 [v1.0.0]); and **NMD_orf_model_v5_4ct**, the deep-learning model predicting NMD from ORF sequence context, comprising model training, evaluation, and attention- and SHAP-based interpretability (https://github.com/peter4244/NMD_orf_model_v5_4ct; DOI: 10.5281/zenodo.21536501 [v2.0.0]). Reproducing Figure 5 requires both NMD_orf_model_v5_4ct, which trains the model and computes the interpretability outputs, and nmd_isopair_analysis, which renders the figure panels and performs the NMD-predictor benchmark. Sequencing data generated in this study are available from the NCBI Gene Expression Omnibus under accession GSE329233, and the source data required to reproduce the analyses are archived at Zenodo (DOI: 10.5281/zenodo.21544337).
