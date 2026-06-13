# Manuscript results → code map

**Manuscript:** *Long Read RNA Sequencing in Primary Lung Cell Types Reveals Principles of Nonsense-Mediated Decay* (Leshem et al., in prep)

**Purpose.** Identify the code that produced every quantitative claim and figure panel in the manuscript. This is **step 1** of the verification workflow — the map is the trace; verification (Pete's 5-step scientific-report protocol from `~/.claude/CLAUDE.md`) follows once the trace is in place.

**Sources used to build this map.** Live Google Doc manuscript (fetched 2026-06-13), the code-map Google Doc shared by Pete (with corrections noted below), `nmd/ONBOARDING.md`, `nmd/CLAUDE.md`, the local `nmd/code/` mirror, and the manuscript Methods text.

**Status of this document.** First-pass complete. Methods + Sections 1–5 all mapped at the per-claim granularity. 144 numbered claims total. Next step: triage which claims to verify first against the 5-step scientific-report protocol.

**Verifiable-locally summary** (claims that can be reproduced end-to-end from `~/claude_projects/` clones on this laptop without cluster or Yul-side access):

| Section | Verifiable / total | Primary blockers |
|---|---|---|
| §1 (Isoform discovery) | ~1 / 15 | Yul-side `Isoform_Landscape.Rmd`, `correlation_analysis.Rmd` |
| §2 (NMD response) | ~10 / 26 | Sharing/specificity computed by Yul; Tan reanalysis Rmd path unknown |
| §3 (Output Lost + PCI) | ~1 / 26 | **All three primary §3 Rmds are Yul-side** (`transcriptional_output.Rmd`, `comparison_analysis.Rmd`, `productive_compensation.Rmd`) — §3 is the most Yul-blocked section |
| §4 (Isopair / splice events / PTC) | ~44 / 46 | Essentially all local — `isopair:` + `nmd:results/.../isopair_wrapper/` |
| §5 (DL model) | ~30 / 31 | Essentially all local — `model:` (trained weights cached) |

About half of the manuscript's quantitative claims are directly verifiable from this laptop, concentrated heavily in §4 and §5.

---

## Conventions

### Repo / canonical-source legend

| Tag | Means |
|---|---|
| `nmd:` | This repo. Canonical at GitHub `peter4244/nmd_isopair_analysis` ≡ Channing GitLab `repjc/nmd_lungcells_2026` ≡ local `~/claude_projects/nmd/`. Pete-side. |
| `yul:` | Yul's `/udd/reyle/nmd_lungcells_2026/code/final/` on Channing. Presumed mirrored to changit at `reyle/nmd_lungcells_2026` (need member access). Not directly readable from this laptop yet — cited by path only. |
| `isopair:` | Isopair package. Canonical at GitHub `peter4244/Isopair` ≡ local `~/claude_projects/Isopair/`. |
| `model:` | DL model. Canonical at GitHub `peter4244/NMD_orf_model_v5_4ct` ≡ local `~/claude_projects/NMD_orf_model_v5_4ct/`. Trained on Northeastern Discovery. |
| `randell:` | Long-read isoform-discovery pipeline at Channing `/proj/regeps/regep00/studies/ExternalCellLines/data/longread/mrna/Randell_Lung_Cells_2025/`. Possibly also a changit repo (TBD — Pete to check with Yul). |

### DIE / DGE provenance

Manuscript-quoted DIE/DGE numbers are **mashr posterior estimates** (lfsr, posterior mean logFC), not raw limma adj.P.Val. The "NMD susceptible" definition is `lfsr < 0.05 AND posterior_mean > 0`, per ONBOARDING §6.

**LR DIE provenance (single Yul-side Rmd doing both limma and mashr):**

```
limma + mashr in one Rmd (Yul-side, canonical)         →  canonical CSV (in this repo)
yul:Isoform-Level_Quantification.Rmd                      nmd:isocall_dge/mashr/nmd_mashr_die_{at,dd,fb,mv}_2026.3.10.csv
```

Pete also has a parallel limma implementation at `nmd:code/isocall_limma_dge_fullmodel_2026.3.1.Rmd`, but per Pete (2026-06-13) Yul did her own limma analysis — so the Pete-side Rmd is **not** the canonical input to the mashr CSVs. Treat Pete's Rmd as a QC / sanity-check parallel implementation. The canonical limma contrasts that feed the manuscript come from Yul's Rmd. *Flag for Yul to confirm.*

**SR DGE provenance (single Yul-side Rmd doing both limma and mashr):**

```
limma + mashr in one Rmd (Yul-side)                    →  canonical CSV (in this repo)
yul:NMD_shortread_dge_fullmodel_2026.5.5.Rmd              nmd:shortread_dge/mashr/nmd_mashr_dge_{at,dd,doali,ddali,fb,mv}_2026.3.10.csv
                                                          nmd:shortread_dge/mashr/mashr_{lfsr,posterior_means,...}_2026.3.10.csv   (4-CT shared)
```

### Outdated material in the code-map Google Doc

Pete flagged "outdated material, especially at the beginning of the document." Concrete corrections applied here:

1. **Long-read quantification used the PacBio Isocall pipeline, not oarfish.** The oarfish minimap2/oarfish_quant pipeline and the downstream `oarfish_gencode49_merged_collapsed_limma_dge_fullmodel_*.Rmd` series are **deprecated**. The current canonical limma step is `nmd:code/isocall_limma_dge_fullmodel_2026.3.1.Rmd`.
2. SR DGE Rmd path: Pete's code-map doc cites Yul's `2026.5.5` vintage. The mashr CSVs are dated `2026.3.10`. **Open question:** did the `2026.5.5` Rmd re-emit the same `2026.3.10` CSVs (unlikely), or is there an intermediate `~2026.3.10` vintage that wrote the CSVs? Flagged for Yul to confirm.

### Code-map gaps

These manuscript analyses are not in Pete's code-map Google Doc and need to be added:

- Isopair package + wrapper pipeline (Section 4)
- DL model training + interpretability (Section 5)
- Tan et al. (2025) mashr reanalysis (Section 2)
- DMSO-only one-vs-rest cell-type-marker mashr analysis (Methods)
- 3′UTR length comparison restricted to GENCODE-annotated CDS (Section 4 final paragraph)

---

## Methods → code

### M1. Cell culture, treatment, RNA isolation

Wet-lab protocols; not in scope for code mapping. Quoted in Methods sections "Isolation and Culturing of Primary Lung Cell Types" and "RNA Isolation."

### M2. Short-read RNA-seq processing (nf-core/rnaseq)

- **Pipeline:** nf-core/rnaseq v3.14.0 + Nextflow 24.04.4 (per Methods).
- **Code:** wherever the nf-core launcher script lives on Channing (path not yet identified — flag for Yul).
- **Outputs feeding downstream analysis:** Salmon-quantified gene/isoform counts via tximport (used only at gene level here; Isocall is canonical for isoform).
- **Where the gene-level counts land:** referenced by `yul:NMD_shortread_dge_fullmodel_2026.5.5.Rmd` (input file path TBD — need to read the Rmd header).

### M3. Long-read RNA-seq processing (PacBio Kinnex → Isocall)

- **Library + sequencing:** PacBio protocol 103-238-700 REV07 (per Methods).
- **Alignment + quantification:** PacBio Isocall pipeline (`github.com/PacificBiosciences/isocall`). NOT oarfish (deprecated per Pete 2026-06-13).
- **Where Isocall outputs land:** Channing `/proj/regeps/.../Randell_Lung_Cells_2025/results/` (specific subfolder TBD — flag for Yul).
- **SQANTI3 classification:** config file at `randell:results/sqanti_runs/merged_collapsed/sqanti3_config_cts_2subj_5reads_2025.12.20_merged-collapsed.yaml`.
- **Filtered count matrix consumed by limma:** `nmd:sqanti/nmd_lungcells/results/nmd_lungcells_filtered.count_matrix.txt` (per `isocall_limma_dge_fullmodel_2026.3.1.Rmd` setup chunk).

### M4. Isoform-landscape characterization (DMSO-only)

- **Code (best guess):** `yul:Isoform_Landscape.Rmd` (per code-map doc).
- **Specifics in Methods:** restricted to DMSO samples; 162,800 SQANTI3-filtered isoforms; per-CT detection threshold = summed count ≥1; Jaccard on binary detection; Spearman/Pearson on mean log₂CPM (edgeR::cpm log=TRUE, prior.count=1).
- **Status:** cite-only (cluster verify with Yul).

### M5. DMSO-only one-vs-rest cell-type-marker mashr

- **Methods text:** "To identify isoforms whose baseline expression distinguishes the four cell types..." — describes restricting DGEList to DMSO, re-applying filterByExpr (min.count=5, min.total.count=10), TMM, voom, design `~ 0 + cell_type`, donor as duplicateCorrelation block, four one-vs-rest contrasts, then mashr.
- **Code:** **gap** — not listed in the code-map Google Doc. Flag for Yul.
- **Output:** the LAE/AT2/MV/FB cell-type-marker counts quoted in Section 1 (28,930 LAE, 22,131 AT2, 21,429 MV, 18,422 FB at 5% FDR).

### M6. Differential expression — SR gene level, LR gene level, LR isoform level

- **SR gene level (limma + mashr in one Rmd):** `yul:NMD_shortread_dge_fullmodel_2026.5.5.Rmd` → CSVs in `nmd:shortread_dge/mashr/`.
- **LR gene level:** code path not in Pete's code-map doc. The manuscript reports "long-read gene-level analysis tested 19,056 genes" — flag for Yul to identify the Rmd that produces this (it may be folded into another script).
- **LR isoform level (limma + mashr in one Rmd, canonical):** `yul:Isoform-Level_Quantification.Rmd` → CSVs in `nmd:isocall_dge/mashr/`. Per Pete (2026-06-13), Yul did her own limma version, so this Rmd is the canonical end-to-end implementation for the manuscript numbers.
- **LR isoform level — parallel Pete-side limma (QC / sanity, NOT canonical for manuscript):** `nmd:code/isocall_limma_dge_fullmodel_2026.3.1.Rmd`. Reads SQANTI-filtered Isocall count matrix; writes limma contrast tables to `nmd:isocall_dge/`. Useful for cross-checking Yul's results but not the source of the mashr CSVs.
- **Convention:** design = `~ cell_type + treatment + cell_type:treatment`, reference = LAE, donor block via duplicateCorrelation.

### M7. mashr (the shared infrastructure)

- **Method:** estimate_null_correlation_simple on 20k random feature subset → mash_1by1 at lfsr<0.05 for strong signals → cov_canonical + data-driven (PCA, ≤5 PCs) → mash fit.
- **Code:** the mashr step happens inside the Rmds above (M5, M6 SR, M6 LR-DIE). There may also be a shared mashr helper function file — flag for Yul.

### M8. GSEA + pathway enrichment (gene-level only)

- **Code:** `nmd:code/gsea_mashr_2026.3.10.R` (per ONBOARDING §9).
- **Output:** `nmd:tmp/gsea_mashr_gene_2026.3.10_run2026-05-18.tsv` (per ONBOARDING).
- **Method:** fgsea against MSigDB (Hallmark, KEGG, Reactome, GO-BP), genes ranked by signed posterior mean logFC, min=15, max=500, FDR<0.05.

### M9. Tan et al. (2025) reanalysis

- **Code:** **gap** — not in code-map doc. Methods describe converting EBSeq PostFC + PPEE to bhat/shat and fitting mashr jointly across the 8 conditions (3 UPF2 + 1 UPF3B × hESC + NPC).
- **Inputs:** Tan et al. Supplementary Tables S1, S2, S4, S6.
- **Status:** flag for Yul / search Yul's `final/` dir for a `tan*` or `wilkinson*` Rmd.

### M10. Transcriptional output lost + PCI

- **Code:** `yul:transcriptional_output.Rmd` (% output lost) and `yul:productive_compensation.Rmd` (PCI). Per code-map doc.
- **Status:** cite-only.

### M11. Isopair pairs analysis

- **Package:** `isopair:` (canonical at `peter4244/Isopair`, vignette `Isopair/vignettes/NMD-attribution.Rmd` is the canonical methods source per ONBOARDING).
- **Pipeline wrapper (force-tracked in this repo):** `nmd:results/isoform_transitions/Version_6.0/isopair_wrapper/` — `01_..05_*.R` + `05_final_report_mashr.Rmd`.
- **3′UTR length analysis on GENCODE-annotated CDS subset (Section 4 final paragraph, n=1,904 NMD vs 22,335 non-NMD):** likely in `05_final_report_mashr.Rmd` ("3'UTR analyses" was added in commit `38c00a5` per `git log`). To confirm by reading the Rmd.

### M12. Deep learning model

- **Code:** `model:` — `03_train.py`, `config.yaml`, `model.py`. Trained on Northeastern Discovery (`/home/p.castaldi/cc/nmd_orf_model_v5_4ct/`).
- **Methods text source:** `model:METHODS.md` (per ONBOARDING §5).
- **Architecture:** two shared-weight CNN branches (AUG window + stop window, 9 channels), linear branch for 5 structural features, fused per-ORF embedding × 5 ORFs, attention aggregation, classification head. BCEWithLogitsLoss w/ pos_weight, Adam (lr=1e-3, differential weight decay), batch=256, fp16 AMP, ReduceLROnPlateau, early stop on val AUC, seed=42.
- **Sweep:** AUG ∈ {100, 500, 1000} × stop ∈ {100, 500, 1000, 2000}; chose AUG=500/stop=500 by held-out AUC.

### M13. Interpretability (DeepSHAP + KernelSHAP + attention)

- **Code:** `model:` — interpretability scripts (filenames TBD; `model:METHODS.md` should enumerate).
- **Joint DeepSHAP:** ORF-0 (AUG + stop + 5 structural) varied while ORFs 1–4 held at observed values; 5 seeds × 500 background transcripts; per-position/channel attributions + per-feature structural attributions.
- **KernelSHAP at embedding level:** decomposes prediction across the 3 sub-encoder embeddings (AUG / stop / structural).
- **Subgroup attribution:** stratifies pooled DeepSHAP by NMD mechanistic subgroups (PTC+, PTC− ref-AUG retained, PTC− ref-AUG lost).
- **Occult-PTC ORF enumeration:** `Isopair::enumerateOrfs()` (in the Isopair package, not the model repo).

---

## Section 1 (Isoform Discovery via Long-read RNA-sequencing) → code

### Paragraph 1 — sequencing summary + SQANTI categories + isoform-count distribution

| # | Claim (verbatim or near-verbatim) | Source code | Notes / status |
|---|---|---|---|
| 1.1 | "26 samples to an average depth of 13 million aligned FLNC reads, 341,638,920 total reads" | Channing Isoseq pipeline output stats (`randell:results/long-read-isoseq/20251202-early-ALL/`). Aggregation script TBD. | Flag for Yul/Pete to identify the aggregation script. |
| 1.2 | "162,800 expressed isoforms across 18,270 genes" | Output of SQANTI3 + filterByExpr. Concrete filterByExpr application is in `nmd:code/isocall_limma_dge_fullmodel_2026.3.1.Rmd` (setup chunk → DGEList → filterByExpr). | Verifiable locally. The "18,270 genes" rollup needs an Rmd that aggregates isoform→gene — flag whether this is in Yul's `Isoform_Landscape.Rmd`. |
| 1.3 | "Isoform length distributions were unimodal and centered near 2.5 kb" (SF: Isoform Length By Cell Type And Treatment) | `yul:Isoform_Landscape.Rmd` (best guess). | Cite-only. |
| 1.4 | "53.7% FSM/ISM (55,770 FSM + 31,667 ISM); 44.4% novel (44,193 NIC + 28,045 NNC); 1.9% other" (SF: SQANTI3 Structural Categories by Cell Type) | SQANTI3 `classification.txt` from `randell:results/sqanti_runs/merged_collapsed/`. Aggregation/percentage script TBD. | Numbers should also be reproducible from the SQANTI classification.txt via a one-off script. |
| 1.5 | "median gene expressed five isoforms (mean 8.9, IQR 2–13, max 97), 13,941 (76.3%) genes ≥2 isoforms, 5,963 (32.6%) ≥10" | Isoform-per-gene aggregation. Likely `yul:Isoform_Landscape.Rmd`. | Cite-only. |

### Paragraph 2 — short-read vs long-read concordance (Figure 1A + SFx)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 1.6 | "mean Pearson r = 0.90 (range 0.83–0.91), mean Spearman ρ = 0.890 (range 0.849–0.901), 26 matched samples, 29,185 common genes" | `yul:correlation_analysis.Rmd` (per code-map doc). | Cite-only. |
| 1.7 | "Per cell type Pearson means 0.880 (FB) — 0.901 (AT2)" | Same. | Cite-only. |
| 1.8 | **Figure 1A** (sample-wise SR↔LR correlation plot) | Render code TBD — likely embedded in `yul:correlation_analysis.Rmd` or in a downstream `Gene-Level_DGE_Summary_mashR.Rmd` plot block. | Cite-only. **Figure-source-script location is the main `figures/` gap.** |

### Paragraph 3 — cell-type-specific expression (Figure 1B + SFx)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 1.9 | "105,938 isoforms expressed at testable levels under DMSO" | DMSO-only filterByExpr re-application step in the one-vs-rest mashr analysis (M5). Code path TBD — flag for Yul. | Cite-only. |
| 1.10 | "28,930 LAE, 22,131 AT2, 21,429 MV, 18,422 FB significant at 5% FDR" | One-vs-rest mashr output (M5). | Cite-only. |
| 1.11 | "Spearman ρ = 0.79 (FB↔MV); 0.43–0.56 between AT2/LAE and others" (SF: Pairwise Expression Similarity) | `yul:Isoform_Landscape.Rmd` (best guess). | Cite-only. |
| 1.12 | "LAE 8,087; AT2 2,267; MV 2,015; FB 1,131 cell-type-restricted isoforms; 0.36–2.24% of each CT's expressed isoforms" | `yul:Isoform_Landscape.Rmd`. | Cite-only. |
| 1.13 | "27.8% (LAE) to 51.7% (MV) of CT-restricted isoforms passed the additional expression filter" | Same. | Cite-only. |
| 1.14 | "Pairwise Jaccard indices 0.86–0.92" | Same. | Cite-only. |
| 1.15 | **Figure 1B** (expressed + restricted isoform counts per CT) | Render code TBD; likely embedded in `yul:Isoform_Landscape.Rmd`. | Cite-only. |

---

## Section 2 (NMD Response) → code

### Paragraph 1 — per-CT DGE/DIE counts + cross-CT sharing of NMD-susceptible features (Figure 1C, 1D)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 2.1 | "25,955 genes per CT" tested at SR level | `yul:NMD_shortread_dge_fullmodel_2026.5.5.Rmd` filterByExpr step → `nmd:shortread_dge/mashr/nmd_mashr_dge_*_2026.3.10.csv` (row count). | Cite-only; CSV row count verifiable locally. |
| 2.2 | "3,122–6,753 significant genes per CT (49–67% NMD susceptible)" — see also Short-Read mashr table (AT2 3,638 / 62.56%; LAE 6,753 / 49.25%; FB 3,122 / 66.18%; MV 3,428 / 67.39%) | Same Rmd; subsetting on `adj.P.Val < 0.05` for "significant" and `lfsr < 0.05 & posterior_mean > 0` for "NMD susceptible". | Verifiable locally by reading the 4 per-CT mashr CSVs. **Watch:** "Significant" here is the *limma adj.P.Val* denominator, not mashr lfsr — the mixed-denominator phrasing is fine but verification needs to use the right column. |
| 2.3 | "162,800 isoforms per CT" tested at LR level | `yul:Isoform-Level_Quantification.Rmd` (filterByExpr) → `nmd:isocall_dge/mashr/nmd_mashr_die_*_2026.3.10.csv`. | Cite-only; CSV row count verifiable locally. |
| 2.4 | "24,803–35,336 significant isoforms per CT (90.9–92.0% NMD susceptible)" — see DIE mashr table (AT2 24,847 / 91.96%; LAE 35,336 / 91.24%; FB 24,803 / 90.92%; MV 27,834 / 90.87%; logFC 1.60–2.97) | Same Rmd; same dual-denominator pattern. | Verifiable locally from CSVs. |
| 2.5 | "34,387 unique NMD susceptible isoforms across CTs" | Likely `yul:interpret_isoform_patterns_mashr_2026.3.10.Rmd` — union over per-CT NMD-susceptible sets. | Cite-only. Verifiable locally by union over the 4 mashr CSVs. |
| 2.6 | "19,803 (57.6%) core targets shared across all CTs; 9,161 (26.6%) CT-specific (LAE 83.4% of those)" | Same Rmd — intersect/setdiff over per-CT NMD-susceptible sets. | Cite-only; locally verifiable from CSVs. |
| 2.7 | **Figure 1C** — volcano plots per CT for LR DIE (logFC vs −log10 lfsr, color by NMD-susceptible threshold) | Figure-render script TBD. Inputs are `nmd:isocall_dge/mashr/nmd_mashr_die_{ct}_2026.3.10.csv`. | Cite-only. **Render script location is an open question.** |
| 2.8 | **Figure 1D** — distribution of posterior mean logFC for NMD-susceptible features at SR-gene vs LR-isoform | Figure-render script TBD. Inputs are both sets of mashr CSVs. | Cite-only. |
| 2.9 | "mean posterior logFC substantially larger at LR isoform than SR gene level" | Same Rmd as 2.8. | Cite-only. |

### Paragraph 2 — NMD-susceptible isoform proportion + expression distribution (SFx panels)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 2.10 | "Isoforms concentrated at low isoform proportions — majority < 50% of parent gene expression, substantial fraction < 10%, span > 4 orders of magnitude of absolute expression" (SF: All Proportion vs Expression: DMSO; SF: NMD Only Proportion vs Expression: DMSO) | Likely `yul:interpret_isoform_patterns_mashr_2026.3.10.Rmd` or `yul:Isoform_Landscape.Rmd`. Inputs: DMSO-only CPM × NMD-susceptible flag from mashr CSVs. | Cite-only. |
| 2.11 | "Clear upward shift in proportion and absolute expression between DMSO and SMG1i" (SF: NMD-Responsive Isoforms: DMSO vs SMG1i) | Same Rmd. | Cite-only. |

### Paragraph 3 — quantitative CT specificity vs pairwise sharing (Figure 1E, SF: Sharing)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 2.12 | Per-CT specificity at SR-gene level: LAE 34.2%, AT2/FB/MV 1.9–5.3% | `yul:interpret_isoform_patterns_mashr_2026.3.10.Rmd` (or sibling) — fraction of NMD-susceptible significant in only that CT. | Cite-only. |
| 2.13 | Per-CT specificity at LR-isoform level: LAE 23.7%, AT2/FB/MV 0.8–3.4% | Same Rmd. | Cite-only. |
| 2.14 | Pairwise sharing at SR-gene level: LAE-involving 51–70%; AT2/FB/MV trio 83–90% | Same Rmd — fraction with concordant direction + within 2-fold magnitude on mashr posteriors. | Cite-only. |
| 2.15 | Pairwise sharing at LR-isoform level: LAE-involving 35–72%; AT2/FB/MV trio 68–83% | Same Rmd. | Cite-only. |
| 2.16 | **Figure 1E** — fraction of NMD-susceptible genes + isoforms classified as CT-specific (bar / dot plot) | Figure-render script TBD. | Cite-only. |

### Paragraph 4 — Tan et al. (2025) reanalysis

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 2.17 | Tan-reported overlap: 16–25% UPF2-dependent NMD-target transcripts shared between hESC and NPC; "176 shared of 709 hESC and 1,088 NPC" | Pulled from Tan et al. 2025 Supplementary Tables S1, S2, S4, S6. Source: their paper, not our code. | External — verifiable from their supp tables. |
| 2.18 | "transcript-level overlap of UPF2-dependent NMD targets between hESCs and NPCs increased to 41–50% (3,069 shared transcripts; 49.8% of hESC and 40.6% of NPC targets)" — concordant direction across both CTs for large majority | Methods M9 describes the procedure (EBSeq PostFC → bhat, PPEE → shat, mashr fit jointly across 8 conditions, lfsr < 0.05 + posterior mean > 0 in same direction). **Code path: gap.** Likely in Yul's `final/` — file probably named `tan*.Rmd` or `wilkinson*.Rmd`. | **Open. Flag for Yul.** |

### Paragraph 5 — GSEA + pathway enrichments (Figure 1F)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 2.19 | "5 pathways significant in 3 of 4 CTs (LAE, FB, MV; FDR 0.002–0.048); AT2 trending (NES 1.65–1.83 below FDR)" — Table X — Pathways Significant | `nmd:code/gsea_mashr_2026.3.10.R` → output TSV `nmd:tmp/gsea_mashr_gene_2026.3.10_run2026-05-18.tsv` (per ONBOARDING §9). Input: signed mashr posterior-mean logFC ranking per CT from `nmd:shortread_dge/mashr/nmd_mashr_dge_*_2026.3.10.csv`. | **Verifiable locally** — both code and output are in the repo. |
| 2.20 | Pathways named: cellular response to topologically incorrect protein, cellular response to unfolded protein, response to topologically incorrect protein, intrinsic apoptotic signaling in response to ER stress, Reactome unfolded protein response | Same GSEA output TSV. | Verifiable locally. |
| 2.21 | "Leading-edge genes converging on ATF4, DDIT3, PPP1R15A, ATF3, CHAC1" | Same GSEA output (leading-edge column from fgsea). | Verifiable locally — need to confirm fgsea was called with leading-edge return. |
| 2.22 | **Figure 1F** — KEGG Protein Processing in ER pathway enrichment for LAE; corresponding panels for other CTs in SFx | Figure-render script TBD. Could be in `yul:Gene-Level_DGE_Summary_mashR.Rmd` (per code-map doc) or in a sibling `*_p2.Rmd`. | Cite-only. |

### Paragraph 6 — CT-specific pathway enrichments (Table X — Top Pathways)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 2.23 | LAE-specific: cilium movement, cilium/flagellum-dependent cell motility | Same GSEA output TSV; LAE-only enriched pathways. | Verifiable locally. |
| 2.24 | AT2-specific: translation-regulatory programs, xenobiotic metabolism | Same. | Verifiable locally. |
| 2.25 | AT2 + MV: regulation of RNA splicing | Same. | Verifiable locally. |
| 2.26 | FB: no CT-specific pathways | Same — null result, recompute. | Verifiable locally. |

### Section 2 figure render-script gap

Figure 1 has 6 panels (A–F). Panel A (SR↔LR sample correlation) was tied to `yul:correlation_analysis.Rmd` in §1. Panels B–F:

- **1B** (CT-restricted isoform bars) — likely `yul:Isoform_Landscape.Rmd`
- **1C** (volcano plots) — render script TBD
- **1D** (logFC distribution) — render script TBD
- **1E** (CT-specific fraction) — likely in `yul:interpret_isoform_patterns_mashr_2026.3.10.Rmd`
- **1F** (KEGG ER pathway, LAE) — render script TBD; likely `yul:Gene-Level_DGE_Summary_mashR.Rmd` or downstream of `nmd:code/gsea_mashr_2026.3.10.R`

This recurring "render script TBD" pattern is itself the open question #8 from §1 — figure rendering isn't isolated to any single script in either the code-map doc or `nmd:code/`.

### Section 2 verifiable-locally summary

These claims can be verified end-to-end on this laptop without cluster or Yul-side access:

- 2.2 (SR sig + NMD-susceptible counts) — directly from `shortread_dge/mashr/*_2026.3.10.csv`
- 2.4 (LR sig + NMD-susceptible counts) — directly from `isocall_dge/mashr/*_2026.3.10.csv`
- 2.5–2.6 (unique / shared / CT-specific isoforms) — same CSVs, recompute the union/intersection
- 2.19–2.21, 2.23–2.26 (GSEA results) — `code/gsea_mashr_2026.3.10.R` + `tmp/gsea_mashr_gene_2026.3.10_run2026-05-18.tsv`

These are the most efficient first targets when we move from mapping into actual verification.

---

## Section 3 (Transcriptional Output Lost + PCI) → code

### Paragraph 1 — % transcriptional output lost (Figure 2A, SFx)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 3.1 | Gene-level % output lost: 11.1% (AT2) — 12.5% (MV) (SF: % Transcriptional Output Lost, Gene Level) | `yul:transcriptional_output.Rmd` — computes Σmax(CPM_SMG1i − CPM_DMSO, 0) / ΣCPM_SMG1i over NMD-susceptible features. Inputs: per-CT mashr CSVs + DMSO/SMG1i CPM tables. | Cite-only. Locally verifiable from the mashr CSVs + CPM matrices in `nmd:shortread_dge/` if Yul wrote out the CPM table. |
| 3.2 | Isoform-level % output lost: 14.4% (AT2) — 18.2% (LAE) | Same Rmd, isoform-level computation. Inputs: per-CT DIE mashr CSVs + DMSO/SMG1i CPM matrices. | Cite-only. |
| 3.3 | Per-isoform percent lost: median 60–67%, IQR 30–100% (SF: Per-Isoform Percent Lost) | Same Rmd — per-isoform fraction = (SMG1i CPM − DMSO CPM)/SMG1i CPM for NMD-susceptible with positive delta. | Cite-only. |
| 3.4 | **Figure 2A** — % output lost at isoform level across CTs | Render script TBD; likely embedded in `yul:transcriptional_output.Rmd`. | Cite-only. |

### Paragraph 2 — isoform vs gene discrepancy + PCI definition (in-text 63–90%)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 3.5 | Isoform-to-gene ratio of % output lost: 1.23–1.52 across CTs | `yul:transcriptional_output.Rmd` — ratio of isoform-level to gene-level totals. | Cite-only. |
| 3.6 | PCI definition (significant unproductive↑ + significant productive↓; LFSR<0.05 in both, signs opposed) | Methods M10. Implementation in `yul:productive_compensation.Rmd` (per code-map doc). | Cite-only. The full pseudo-feature aggregation framework (protein-coding-only, ambiguous-multi-locus exclusion, productive/unproductive split, separate limma + mashr fit) is documented in Methods. |
| 3.7 | "63–90% of genes with significant NMD-isoform accumulation showed concurrent significant productive-isoform decrease" | `yul:productive_compensation.Rmd` — count over genes meeting both PCI criteria, divided by count of genes with significant unproductive accumulation. | Cite-only. **Verification target:** this is a key headline number; pin down the exact denominator and PCI-class threshold from the Rmd when accessible. |

### Paragraph 3 — SR↔LR concordance of NMD calls (Figure 2B)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 3.8 | "23–30% of genes with NMD susceptible isoforms (LR) classified NMD susceptible at gene level (SR)" | `yul:comparison_analysis.Rmd` (per code-map doc; "Short-read vs Long-read Comparison Analysis (logFC, etc.)"). Roll up LR isoform NMD-susceptible flags to gene, intersect with SR gene-level NMD-susceptible set. | Cite-only. |
| 3.9 | "60–68% of SR NMD-susceptible genes had at least one detectable LR NMD-susceptible isoform" | Same Rmd, reverse intersection. | Cite-only. |
| 3.10 | "Approximately 40% of SR NMD-susceptible genes lacked sufficient LR isoform coverage" | Same Rmd — depth threshold + filterByExpr drop-out analysis. **Methods don't explicitly define "sufficient coverage"** — pin down when accessible. | Cite-only. **Watch:** quantitative threshold for "sufficient coverage" is not specified in Methods text; verify against Rmd. |
| 3.11 | **Figure 2B** — SR↔LR NMD-call concordance (Venn / scatter / bar) | Render script TBD; likely embedded in `yul:comparison_analysis.Rmd`. | Cite-only. |

### Paragraph 4 — PCI control analysis on non-NMD genes (SFx)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 3.12 | Productive isoforms in NMD-isoform-containing genes: median logFC −0.41 to −0.67 (across the 4 CTs) | `yul:productive_compensation.Rmd` — control analysis per Methods M10 (parallel limma + mashr on genes lacking NMD-targeted isoforms, then per-CT logFC distribution comparison). | Cite-only. |
| 3.13 | Productive isoforms in non-NMD-isoform genes: median logFC −0.006 to −0.015 | Same. | Cite-only. |
| 3.14 | Wilcoxon p < 10⁻⁵ in all CTs (one-sided rank-sum) | Same; explicit one-sided Wilcoxon between the two productive-logFC distributions. | Cite-only. |

### Paragraph 5 — GPR180 worked example (Figure 2D, 2E, 2F)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 3.15 | GPR180 in FB: 3 novel unproductive isoforms 5.2 → 15.7 CPM under SMG1i | `yul:productive_compensation.Rmd` — gene-specific drill-down on GPR180 in FB. Inputs: per-isoform CPM matrix from FB. | Cite-only. **Verifiable:** if Yul's Rmd outputs a gene-detail TSV for GPR180 we can confirm directly. |
| 3.16 | Productive isoforms 40.5 → 30.2 CPM | Same. | Cite-only. |
| 3.17 | Productive fraction 88.6% → 65.7% | Computed from above. | Locally derivable. |
| 3.18 | **Figure 2D** — GPR180 productive/unproductive isoform expression DMSO vs SMG1i | Render script TBD. Likely in `yul:productive_compensation.Rmd`. | Cite-only. |
| 3.19 | **Figure 2E** — GPR180 isoform structures grouped by productive classification, "rendered using Isopair" | `isopair:` plotting function (specific function TBD — check `isopair:R/plot_*.R` or `isopair:vignettes/`). Driver call in Yul's Rmd. | Cite-only. The fact that the figure caption explicitly attributes rendering to Isopair means the function must be in the Isopair package — locally readable. |
| 3.20 | **Figure 2F** — GPR180 isoform-level logFC, "rendered using Isopair" | Same — Isopair plotting function. | Cite-only. |

### Paragraph 6 — PCI gene baseline properties + logistic regression (Figure 2C)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 3.21 | PCI genes more highly expressed at baseline: 109–164 CPM (PCI) vs 26–48 CPM (non-PCI NMD susceptible) — across the 4 CTs | `yul:productive_compensation.Rmd` — per-CT median (or mean) baseline total gene CPM stratified by PCI status. Inputs: gene-level DMSO CPM + PCI flag. | Cite-only. |
| 3.22 | Baseline productive fraction: 0.98 (PCI) vs 0.82–0.94 (non-PCI NMD susceptible) | Same; per-CT productive fraction = productive pseudo-feature CPM / total gene CPM under DMSO. | Cite-only. |
| 3.23 | Logistic regression: total gene CPM OR 2.1–3.0 per SD (p<10⁻⁵ in all CTs); dominant isoform fraction contributes independently in 3 of 4 CTs | Same Rmd. Methods M10 specifies per-CT logistic regression with standardized predictors (log₁₀ total gene CPM, dominant isoform fraction, productive isoform count). | Cite-only. |
| 3.24 | **Figure 2C** — violin plots: baseline expression and dominant productive isoform fraction, PCI vs non-PCI | Render script TBD; likely in `yul:productive_compensation.Rmd`. | Cite-only. |

### Paragraph 7 — PCI pathway enrichment (SF: compensation pathways)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 3.25 | "PCI in all 4 CTs" set: enriched for oxidative phosphorylation / mito energy metabolism — specifically: mitochondrial electron transport NADH→ubiquinone (p=1.4×10⁻⁶), respiratory chain complex I assembly/function (p=8.3×10⁻⁶), proton motive force / ATP synthesis (p=2.0×10⁻⁵), ATP metabolic process (p=3.3×10⁻²) | `yul:productive_compensation.Rmd` (Methods M10 specifies topGO with elim algorithm + Fisher's exact, background = all pseudo-feature-passable genes per CT, threshold elim Fisher p<0.05, "shared = sig in ≥2 CTs"). | Cite-only. **Watch:** the p-values quoted aren't labeled as cross-CT vs single-CT — verify whether they're per-CT or pooled. Methods M10 says "terms significant in at least two cell types were reported as cross-cell-type shared enrichments," so these are presumably the cross-CT shared list. |
| 3.26 | Additional enrichments: mRNA splicing, ribosomal subunits, proteasomal components | Same Rmd / topGO output. | Cite-only. |

### Section 3 figure render-script gap

Figure 2 has 6 panels (A–F):

- **2A** (% output lost isoform-level) — likely `yul:transcriptional_output.Rmd`
- **2B** (SR↔LR concordance) — likely `yul:comparison_analysis.Rmd`
- **2C** (PCI vs non-PCI violins) — likely `yul:productive_compensation.Rmd`
- **2D** (GPR180 expression bars/lines) — likely `yul:productive_compensation.Rmd`
- **2E** (GPR180 isoform structures) — **Isopair plotting function**, locally readable in `isopair:` package
- **2F** (GPR180 isoform-level logFC) — **Isopair plotting function**, locally readable in `isopair:`

The 2E + 2F Isopair-rendered panels are the first concrete figure-script targets I can identify locally. The other 4 are still TBD on Yul's side.

### Section 3 verifiable-locally summary

- **3.17** (GPR180 productive fraction shift 88.6% → 65.7%) — arithmetic from 3.15/3.16, locally derivable once Yul's per-isoform CPM table for GPR180 is accessible
- **Figure 2E + 2F rendering logic** — Isopair source is `isopair:` (local). The driver calls (with GPR180 as input) live in Yul's Rmd but the actual plotting functions are locally readable.

Most of §3 still requires Yul-side access. Compared to §2 (high local-verifiability), §3 is heavily concentrated in `yul:transcriptional_output.Rmd`, `yul:productive_compensation.Rmd`, and `yul:comparison_analysis.Rmd` — making this section a top priority for getting Yul's repo access set up.

---

## Section 4 (Isopair: attributing NMD to splicing events) → code

> **⚠ POST-2026-06-13 PM REVISION — supersedes the per-claim § 4 entries below for any claim related to Figure 3 (new combined figure of old Fig 3 + Fig 4).**
>
> The figure structure, manuscript numbers, and methodology framework for §4 have been substantially reframed under Pete's clarified policy: "TD2 CDS annotations are unreliable, so for all analyses that depend on identifying the stop codon, we limit to isoform pairs where reference AUG tracing can be performed."
>
> **New population structure** (`figures/multipanel/figure3_isopair_and_ptc/data_export.R`):
>
> | Layer | NMD n | Control n | Used by |
> |---|---|---|---|
> | **pop_BC** = Stage 2 gene-matched | 3,009 | 3,009 | Panels B, C |
> | **pop_traceable** = ref-AUG-traceable subset | 2,289 | 1,763 | Panel D |
> | **pop_ptc_plus** = ref-AUG PTC+ (effectively_ptc) | 1,912 | 288 | Panels E, F |
>
> **New headline numbers**:
> - NMD PTC rate (ref-AUG-defined): **83.5%** (1,912 / 2,289) — replaces the old "44%" claim
> - Control PTC rate: **16.3%** (288 / 1,763) — replaces the old "5%" claim
> - Fold enrichment: **~5.1×** (replaces "15-fold")
> - "85% combined PTC explanation" → now directly **83.5%** (single source, not a TD2 + ref-AUG composite)
> - "3,674 genes" → **3,099 eligible genes** (3,009 yielding pair sets) — vintage drift from older mashr classification
>
> **Old per-claim entries below (4.1–4.46) are SUPERSEDED** for any claim that touches the new Fig 3. The new mapping is:
>
> | New Fig 3 panel | Mirrors old | Source code |
> |---|---|---|
> | New Fig 3 Panel A | (new — schematic) | `figures/multipanel/figure3_isopair_and_ptc/figure3_panelA_pair_concept.py` (ported from `make_pair_concept_figure.R`) |
> | New Fig 3 Panel B (sequence similarity) | old Fig 3B | `figures/multipanel/figure3_isopair_and_ptc/figure3_panelB_sequence_similarity.py` |
> | New Fig 3 Panel C (event prevalence) | old Fig 3C | `figures/multipanel/figure3_isopair_and_ptc/figure3_panelC_event_prevalence.py` |
> | New Fig 3 Panel D (stop-to-EJC distance) | old Fig 4A | `figures/multipanel/figure3_isopair_and_ptc/figure3_panelD_stop_codon_distance.py` |
> | New Fig 3 Panel E (PTC-causing events) | old Fig 4C | `figures/multipanel/figure3_isopair_and_ptc/figure3_panelE_ptc_event_attribution.py` |
> | New Fig 3 Panel F (mechanism breakdown) | old Fig 4D | `figures/multipanel/figure3_isopair_and_ptc/figure3_panelF_mechanism_breakdown.py` |
>
> **Old 3D (gain/loss direction)** and **old 4B (NMD logFC × distance)** moved to supplement per Pete's plan.
>
> **§4 prose updates**: see `paper/section4_findreplace_2026-06-13.md` for Google Doc find/replace pairs.
>
> **Methodology files** (in `figures/multipanel/figure3_isopair_and_ptc/`):
> - `figure3_panelA_methodology.md`, `figure3_panelB_methodology.md`, `figure3_panelC_methodology.md`, `figure3_panelD_methodology.md`, `figure3_panelE_methodology.md`, `figure3_panelF_methodology.md`
>
> **Upstream Rmd updates pending** (Task #23): `compute_ptc_rates_row()` and `goal2-ptc-mechanisms` chunk in `05_final_report_mashr.Rmd` should be updated to use the new ref-AUG-traceable scope as the canonical PTC analysis population. This will let future Rmd renders match Figure 3 directly without per-panel filter chains in `data_export.R`.

**Primary code home for §4:** `nmd:results/isoform_transitions/Version_6.0/isopair_wrapper/` (force-tracked in this repo) — wrapper pipeline (`01_..06_*.R/Rmd`), with `05_final_report_mashr.Rmd` as the main report and named sub-scripts (`05r_ref_atg_analysis.R`, `05k_utr5_all_isoforms.R`, `05l_unified_model.R`, etc.) for specific analyses. Backing package functions live in `isopair:` (`peter4244/Isopair`). Both are locally readable — **§4 is the most directly verifiable-from-this-laptop section in the manuscript.**

**Key package-to-script bridge** (from `05r_ref_atg_analysis.R` header):

```
05r_ref_atg_analysis.R  →  Isopair::traceReferenceAtg()  →  categorical NMD+/PTC- subgroup labels
                            (effectively_ptc, truncated_no_ejc, no_downstream_ejc,
                             ref_atg_lost, no_ref_cds, mapping_failed)
```

### Paragraph 1 — pair-set construction (3,674 genes)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 4.1 | "12 categories of splicing events" enumerated by Isopair | `isopair:R/event-detection.R` — defines the 12 categories. Per Methods M11 + Fig 4D legend: Alt TSS, Alt TES, exon skipping, intron retention, A5SS, A3SS, MXE, alt first/last exon, terminal exon extension/truncation, combination events. Caption defines abbreviations: SE, IR, A3SS, A5SS, Partial IR 5'/3', IR diff 5'/3', Alt_TSS, Missing_Int. | **Verifiable locally.** The 12-category list in the Rmd output should match the enumeration in `event-detection.R`. |
| 4.2 | Correctness verified by ability to reconstruct reference exon structure | `isopair:R/reconstruction.R` — reconstruction function. Per ONBOARDING the validator is built into the package. | **Verifiable locally.** |
| 4.3 | "3,674 genes from which the isoform pairs were drawn" (≥3 expressed isoforms, ≥1 NMD-susceptible) | `nmd:results/.../isopair_wrapper/01_prepare_data_mashr.R` + `02_build_profiles_mashr.R` — pair-set filtering. Per Methods: "≥5% of overall gene expression in either DMSO or SMG1i and ≥5 reads in ≥1 sample"; non-NMD definition was lowered from adj.P.Val > 0.50 to > 0.30 per ONBOARDING §6. | **Verifiable locally.** Run the 01–02 scripts; row count of the resulting profile object should be 3,674. |
| 4.4 | Median 7 isoforms per gene (SF: Isoform Count in Isoform Pair Sets) | `nmd:results/.../05_final_report_mashr.Rmd` — summary statistic chunk. | Verifiable locally. |
| 4.5 | Reference isoform: median 70% of parent gene expression; 75% of references >50% of parent gene expression | Same Rmd — DMSO baseline CPM per gene + per reference isoform. | Verifiable locally. |
| 4.6 | NMD + non-NMD comparators expressed at similar levels in SMG1i samples (SF: Expression Levels of Isoform Pair Sets) | Same Rmd. | Verifiable locally. |
| 4.7 | Transcript length: median 3,067 (reference) / 3,107 (NMD comparator) / 2,906 (Control comparator), p<0.001 for Control vs reference and Control vs NMD (SF: Transcript Length Comparison) | Same Rmd — Wilcoxon test on transcript-length distributions. | Verifiable locally. |
| 4.8 | Flowchart for the 3,674-gene pair construction (SF: NMD pairs flowchart) | Render script TBD. Either embedded in the Rmd or a separate script (could be a manual figure). | Cite-only. |
| 4.9 | **Figure 3A** (schematic illustrating reference / NMD-comparator / Control-comparator pair construction) | Manual illustration likely — no automated render needed (schematic, not data plot). | Cite-only / manual figure. |

### Paragraph 2 — Figure 3 panels B, C, D (sequence similarity, event prevalence, gain/loss direction)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 4.10 | "NMD isoforms shared more sequence content with reference than Control isoforms" — sequence similarity comparison | `isopair:R/compare-sets.R` — sequence-similarity computation. Driver in `05_final_report_mashr.Rmd`. | **Verifiable locally.** |
| 4.11 | **Figure 3B** (sequence similarity of NMD vs Control pairs) | Render code in `05_final_report_mashr.Rmd` (per commit `38c00a5`: "drop standalone figure scripts" — figures are now in-Rmd). | Verifiable locally. |
| 4.12 | "Skipped exon events twice as frequent in NMD pairs vs Controls" — splice-event prevalence | `isopair:R/event-detection.R` + driver in `05_final_report_mashr.Rmd` — per-category event counts in NMD vs Control pairs. | **Verifiable locally.** |
| 4.13 | "Intron retention more common in Control pairs" | Same. | Verifiable locally. |
| 4.14 | **Figure 3C** (prevalence of splice event categories in NMD vs Control) | Render in `05_final_report_mashr.Rmd`. | Verifiable locally. |
| 4.15 | "Terminal events (Alt TSS / TES) and skipped exons more frequently led to sequence GAIN in NMD; intron retention led to sequence gain in Controls" | `isopair:R/event-detection.R` (GAIN/LOSS semantics — clarified in commit `d242798`) + driver in `05_final_report_mashr.Rmd`. | Verifiable locally. |
| 4.16 | **Figure 3D** (gain/loss direction by event type) | Render in `05_final_report_mashr.Rmd`. | Verifiable locally. |

### Paragraph 3 — PTC enrichment + Figure 4

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 4.17 | "44% PTC rate in NMD susceptible isoforms vs 5% in Controls, 15-fold enrichment, p<0.001" | `isopair:R/ptc-attribution.R` + `ptc.R` — PTC detection via 50-nt rule. Driver in `05_final_report_mashr.Rmd`. | **Verifiable locally.** |
| 4.18 | "Large enrichment of stop codons far upstream of the terminal exon junction in NMD isoforms" | `isopair:R/ptc.R` + `spatial.R` — stop-codon-to-last-EJC distance computation. Driver in `05_final_report_mashr.Rmd`. | Verifiable locally. |
| 4.19 | **Figure 4A** (distribution of stop-codon distances from last EJC, NMD vs Control) | Render in `05_final_report_mashr.Rmd`. | Verifiable locally. |
| 4.20 | "Dose-response: NMD logFC increases with distance from stop codon to last EJC, peak ~500 nt" | `05_final_report_mashr.Rmd` — joins mashr posterior logFC with PTC distance. | Verifiable locally. |
| 4.21 | **Figure 4B** (NMD log₂FC × PTC-to-last-EJC distance) | Render in `05_final_report_mashr.Rmd`. | Verifiable locally. |
| 4.22 | "NMD response peaks at 4–5 downstream EJCs from PTC" (SF: NMD Effect Size by Number of EJCs) | Same Rmd — joins mashr posterior logFC with downstream EJC count. | Verifiable locally. |
| 4.23 | "Roughly evenly split between frameshift and in-frame splicing events, smaller contribution from 3′UTR splicing" | `isopair:R/ptc-attribution.R` — PTC introduction-mechanism classifier. Driver in `05_final_report_mashr.Rmd`. | Verifiable locally. |
| 4.24 | "Majority of PTCs caused by exon skipping; significant enrichment for A5SS and A3SS" | Same — PTC-introducing event-type breakdown. | Verifiable locally. |
| 4.25 | **Figure 4C** (PTC-introducing event-type prevalence) | Render in `05_final_report_mashr.Rmd`. | Verifiable locally. |
| 4.26 | **Figure 4D** (frameshift / in-frame / 3′UTR PTC-introduction distribution) | Render in `05_final_report_mashr.Rmd`. | Verifiable locally. |

### Paragraph 4 — NMD+/PTC− analysis + Figure 5 (after panel-drop renumbering)

The Figure 5 panel-drop change (current task #7 from §6 prep): A and B dropped, C→A, D→B, E→C, F→D. **Mapping below uses the NEW panel letters.**

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 4.27 | "56% of NMD susceptible coding isoform pairs did not appear to contain PTCs in initial analysis" (= NMD+/PTC− set, n ≈ 3,674 × 0.56 ≈ 2,057; actual gene-pair count TBD) | `05_final_report_mashr.Rmd` — initial PTC call before ref-AUG tracing, based on `isopair:R/ptc.R` using TD2 CDS. | **Verifiable locally.** |
| 4.28 | "88% of NMD+/PTC− isoforms contain the reference AUG in their transcript sequence" | `05r_ref_atg_analysis.R` (via `Isopair::traceReferenceAtg()`). Categories: `effectively_ptc`, `truncated_no_ejc`, `no_downstream_ejc`, `ref_atg_lost`, `no_ref_cds`, `mapping_failed`. "Has reference AUG" = not `ref_atg_lost` and not `no_ref_cds`. | **Verifiable locally.** |
| 4.29 | "72% of those ORFs contained very early PTCs" — Figure 5A (post-renumber) | `05r_ref_atg_analysis.R` — `effectively_ptc` category. The 72% denominator is the AUG-intact set (88% of NMD+/PTC−). | **Verifiable locally.** |
| 4.30 | **Figure 5A** (most NMD+/PTC− have a PTC when tracing the reference AUG) | Render in `05_final_report_mashr.Rmd` or sibling Rmd. | Verifiable locally. |
| 4.31 | "In NMD+/PTC− with ref-AUG-defined PTC, TD2 called CDS is always longer (length-based bias)" — Figure 5B | `05_final_report_mashr.Rmd` — length comparison of ref-AUG-defined CDS vs TD2-called CDS over the `effectively_ptc` subset. Likely also references `05v_model_comparison.R` (CDS-caller comparison). | **Verifiable locally.** |
| 4.32 | **Figure 5B** (TD2 CDS length vs ref-AUG CDS length) | Render in `05_final_report_mashr.Rmd`. | Verifiable locally. |
| 4.33 | "Kozak strength of reference AUG vs TD2-called AUG: reference is generally more attractive" — Figure 5C | `05_final_report_mashr.Rmd` + Kozak-scoring function (likely `isopair:` or `ORFik` per Methods). | Verifiable locally. |
| 4.34 | **Figure 5C** (Kozak strength comparison) | Render in `05_final_report_mashr.Rmd`. | Verifiable locally. |
| 4.35 | "Splice-event profile of NMD+/PTC− (post-ref-AUG-tracing) indistinguishable from NMD+/PTC+" — Figure 5D | `05_final_report_mashr.Rmd` — re-runs the event-prevalence comparison on the ref-AUG-reclassified subset. | Verifiable locally. |
| 4.36 | **Figure 5D** (splice-event profile, ref-AUG-reclassified PTC+ subset vs original PTC+) | Render in `05_final_report_mashr.Rmd`. | Verifiable locally. |
| 4.37 | "Overall 85% of NMD+ isoforms in our isoform pairs can be explained by the presence of PTCs (original 44% + occult-PTC subset under ref-AUG tracing)" | Computation across the original and ref-AUG-reclassified PTC sets. Logic in `05_final_report_mashr.Rmd`. | Verifiable locally. Headline number worth pinning the exact arithmetic for during verification. |

### Paragraph 5 — Remaining ~15% / three subgroups (5′UTR analysis)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 4.38 | "498 isoforms = remaining ~15% (NMD+ but no PTC even after ref-AUG tracing)" | `05r_ref_atg_analysis.R` — sum of `truncated_no_ejc + no_downstream_ejc + ref_atg_lost + no_ref_cds + mapping_failed`. | **Verifiable locally.** |
| 4.39 | "Subset of 71 isoforms had markedly truncated ORFs despite no downstream EJCs" | `05r_ref_atg_analysis.R` — `truncated_no_ejc` category count. | **Verifiable locally.** |
| 4.40 | "Three subgroups of NMD+/PTC−" | `05l_unified_model.R` — subgroup unification (per commit `38c00a5`: "subgroup unification"). Three subgroups likely = `truncated_no_ejc` + `no_downstream_ejc` + `ref_atg_lost` (the AUG-intact-but-not-effectively-PTC trio). | **Verifiable locally.** |
| 4.41 | "NMD effect size in these two PTC− subgroups was lower than PTC+" (SFx) | `05_final_report_mashr.Rmd` — joins subgroup classification with mashr posterior logFC. | Verifiable locally. |
| 4.42 | "3′UTR mean length shorter than controls in each subgroup — rules out 3′UTR-length as the trigger" | `05_final_report_mashr.Rmd` — 3′UTR length distribution per subgroup vs Control comparator. | Verifiable locally. |
| 4.43 | "5′UTR lengths markedly longer in two of the groups, accounting for 86% of the isoforms" | `05k_utr5_all_isoforms.R` + driver in `05_final_report_mashr.Rmd`. | **Verifiable locally.** |
| 4.44 | "Markedly longer uORFs in these two groups" | `isopair:R/uorf.R` (uORF detection function) + driver in `05_final_report_mashr.Rmd`. | Verifiable locally. |

**Open question (§4 figure-naming):** the manuscript text in this paragraph cites "Fig X, panel A/B/C/D" for the three-subgroup + 3′UTR + 5′UTR + uORF panels, but after the Figure 5 panel-drop (panels A/B = original 5′UTR / TSS-schematic), these panels don't match the current Figure 5 (which now goes A = PTC recovery, B = TD2 length, C = Kozak, D = splice events). **These cited "Fig X" panels probably belong to an unnamed figure that may need to become Figure 6 (or a supplementary).** Worth flagging when reviewing the manuscript text alongside the panel-drop edits.

### Paragraph 6 — 3′UTR length on GENCODE-annotated CDS subset

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 4.45 | "n=1,904 NMD isoforms and 22,335 non-NMD" — analysis restricted to isoforms with GENCODE-annotated CDS to avoid computational-CDS bias | `05_final_report_mashr.Rmd` — added in commit `38c00a5` ("3'UTR analyses"). Filter: GENCODE CDS annotation present. | **Verifiable locally.** |
| 4.46 | "3′UTRs were shorter in NMD susceptible (605 [299–1,119] nt vs 919 [360–1,980] in non-NMD, p<0.001)" | Same. Wilcoxon test on 3′UTR length distributions. | **Verifiable locally.** |

### Section 4 figure render-script summary

In sharp contrast to §1–3, **the figure render scripts for §4 are essentially all local**. Per commit `38c00a5` ("drop standalone figure scripts"), Figure 3, Figure 4, and Figure 5 are rendered inside `05_final_report_mashr.Rmd` itself. The Isopair package provides:

- `isopair:R/visualization.R::plotIsoformPair()` — the function used for Figure 2E and any subsequent isoform-structure renders (also used for Figure 5 schematic panels if any).

This makes Figures 3, 4, 5 the first complete-and-verifiable figure set for the manuscript on this laptop.

### Section 4 verifiable-locally summary

**Essentially all of §4 is locally verifiable.** The full pipeline (Isopair package + wrapper scripts + main Rmd) is in `~/claude_projects/Isopair/` and `~/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/`. Concrete first verification targets:

- **4.3** (3,674 genes) — run 01_prepare_data_mashr.R + 02_build_profiles_mashr.R, count rows in the profile output
- **4.17** (44% vs 5% PTC enrichment) — `05_final_report_mashr.Rmd` chunk
- **4.27–4.29** (NMD+/PTC− → ref-AUG tracing → 88%/72%) — `05r_ref_atg_analysis.R` direct run; categorical breakdown of `traceReferenceAtg()` output
- **4.37** ("85% explained by PTCs") — composite arithmetic; pin down the exact union arithmetic from the Rmd
- **4.45–4.46** (GENCODE-CDS 3′UTR length comparison) — `05_final_report_mashr.Rmd` chunk

§4 will be the most efficient section to verify end-to-end. Worth doing first when we move from mapping to verification, before tackling the Yul-side-blocked §3.

---

## Section 5 (Interpretable Predictive Model for NMD) → code

**Primary code home for §5:** `model:` = `peter4244/NMD_orf_model_v5_4ct` = `~/claude_projects/NMD_orf_model_v5_4ct/`. Trained on Northeastern Discovery (cluster path `/home/p.castaldi/cc/nmd_orf_model_v5_4ct/`). Methods text source of truth: `model:METHODS.md`. **Like §4, §5 is essentially fully verifiable on this laptop** (training itself requires GPU + the HDF5 input; everything downstream — evaluation, attribution, figure rendering — is local).

**Pipeline overview (local files):**

| Stage | Script | What it does |
|---|---|---|
| Data prep | `model:data_prep.py` | Builds HDF5 from transcript sequence + ORF features + labels |
| Training | `model:03_train.py` + `model:config.yaml` + `model:model.py` | Trains the CNN+attention model |
| Eval (AUC/AUPRC) | `model:evaluate.py` | Held-out evaluation |
| Attention attribution | `model:04_interpret_attention.py` | Per-ORF attention weights for each transcript |
| Structural-feature attribution | `model:05_interpret_structural.py` | Per-feature importance |
| Joint DeepSHAP (ORF-0 sequence + structural) | `model:06_export_deepshap_tsv.py` + `model:deepshap.py` + SLURM wrappers `slurm_deepshap_joint.sh`, `slurm_deepshap_joint_orf1_4.sh` | Per-position/channel sequence attributions + per-feature structural attributions; 5 seeds × 500 background |
| KernelSHAP at embedding level (3-branch decomposition) | `model:11_kernel_shap_branches.py` | Exact additive Shapley over the AUG / STOP / structural sub-encoder embeddings |
| Motif analysis (Kozak, UGA, +4) | `model:07_motif_analysis.py` + `scripts/export_joint_motif_logos.py` | Aggregated motif-level attributions for sequence reporting |
| Subgroup-stratified DeepSHAP | `model:08_export_subgroup_deepshap_tsv.py` + `model:09b_export_subgroup_profiles.py` | Per-mechanistic-subgroup attribution comparison |
| GC content + auxiliary channels | `model:09_export_gc_content.py`, `09_export_junction_ordinal.py`, `09_export_polya.py` | Channel-specific exports for SFx panels |
| uORF attention specific analysis | `model:infer_uorf_attention.py` + `audit_uorf_attention.R` + `compute_uorf_attention_metrics.R` | uORF attention follow-up |
| Figure 6A (architecture) | `model:make_architecture_figure.R` | Architecture schematic |
| Figure 6 SHAP panels | `model:make_shap_interpretation_figure.R` | DeepSHAP-based interpretation panels |
| Analysis report | `model:orf_model_report_v5.Rmd` | End-to-end analysis report consuming the TSV exports above |

### Paragraph 1 — model architecture + input encoding + train/test split

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 5.1 | "Five candidate ORFs per transcript, priority: (1) reference CDS ORF when dominant non-NMD isoform's start is present; (2) TD2 CDS when different from (1); (3) remaining slots filled by top Kozak-scoring ORFs from ORFik scan" | `model:data_prep.py` + `model:METHODS.md`. Implementation in data_prep.py — ORF ranking logic. | **Verifiable locally.** |
| 5.2 | "9-channel encoding per position: A/C/G/T (4) + EEJ position + 50-bp rolling GC + 3 reading-frame channels relative to AUG" | `model:data_prep.py` + `model:METHODS.md`. | Verifiable locally. |
| 5.3 | "5 per-ORF structural features: frac_start, frac_stop, is_ref_cds, is_sqanti_cds, n_downstream_ejc" | `model:data_prep.py` + `model:METHODS.md`. | Verifiable locally. |
| 5.4 | "Two shared-weight CNN branches (AUG window, stop window) + linear branch (structural) → fused per-ORF embedding × 5 ORFs → attention aggregation → classification head" | `model:model.py`. | Verifiable locally. |
| 5.5 | "Chr 2, 4 for val; Chr 1, 3, 5, 7 for test; paralogs excluded from test set" | `model:data_prep.py` + paralog list (likely `model:results_4ct/` or a TSV in the repo). | **Verifiable locally.** The paralog exclusion ties back to `nmd:results/.../isopair_wrapper/05u_paralog_annotation.R`. |
| 5.6 | "BCEWithLogitsLoss with pos_weight = n_neg/n_pos; Adam (lr=1e-3, differential weight decay 1e-3 CNN / 1e-4 elsewhere); batch=256; fp16 AMP; ReduceLROnPlateau factor=0.5 patience=5; early-stop patience=10 on val AUC; max 100 epochs; seed=42" | `model:03_train.py` + `model:config.yaml`. | Verifiable locally — open `config.yaml` and confirm values match Methods text. |
| 5.7 | "Sweep over AUG ∈ {100, 500, 1000} × stop ∈ {100, 500, 1000, 2000}; selected AUG=500 / stop=500 by held-out AUC" | Sweep driver script (TBD — `03_train.py` may run individual configs while a sweep orchestrator like SLURM array job picks combinations). The selected config is in `model:config.yaml`. | Cite-only for the sweep results table (config-by-config AUC); the chosen config is verifiable locally. |
| 5.8 | **Figure 6A** — architecture schematic | `model:make_architecture_figure.R`. | **Verifiable locally.** |

### Paragraph 2 — performance + block-level importance (Figure 6B, C, D)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 5.9 | "AUC=0.93, AUPRC=0.83 on held-out test set" | `model:evaluate.py` — produces test-set AUC / AUPRC. | **Verifiable locally** (assuming trained-model weights are in `model:results_4ct/`). |
| 5.10 | **Figure 6B** — ROC / PR curve at AUC=0.93 / AUPRC=0.83 | `model:make_shap_interpretation_figure.R` (or `model:evaluate.py` if rendered there). | Verifiable locally. |
| 5.11 | "Roughly ⅔ of predictive information from ORF structural data, ⅓ from START+STOP sequence" — Shapley block-level decomposition | `model:11_kernel_shap_branches.py` — KernelSHAP at embedding level, exact additive Shapley across the 3 sub-encoders (AUG / STOP / structural). | **Verifiable locally.** |
| 5.12 | "STOP sequence ~3× as important as START sequence" | Same `11_kernel_shap_branches.py` — ratio of STOP to START sub-encoder attributions. | Verifiable locally. |
| 5.13 | **Figure 6C** — block-level Shapley decomposition (3-branch) | `model:make_shap_interpretation_figure.R`. | Verifiable locally. |
| 5.14 | "Of individual ORF structural features, EJC count was by far the most important" | `model:05_interpret_structural.py` + `model:06_export_deepshap_tsv.py` (per-feature DeepSHAP). | Verifiable locally. |
| 5.15 | **Figure 6D** — ranked per-feature structural importance | `model:make_shap_interpretation_figure.R`. | Verifiable locally. |

### Paragraph 3 — nucleotide-level attributions (Figure 6E, F, SFx)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 5.16 | "Importance concentrated within 50 nt preceding the start or stop site" (SFx — Shapley value profile across windows) | `model:06_export_deepshap_tsv.py` — per-position attribution. Aggregation in `model:make_shap_interpretation_figure.R`. | Verifiable locally. |
| 5.17 | "Both windows: more importance for NMD than Control isoforms" | Same — stratify per-position attributions by NMD label. | Verifiable locally. |
| 5.18 | "START window: model learned Kozak sequence; no weight on AUG itself because invariant across training ORFs" | `model:07_motif_analysis.py` + `scripts/export_joint_motif_logos.py`. | Verifiable locally. |
| 5.19 | **Figure 6E** — nucleotide-level attribution around AUG (Kozak motif logo) | `model:make_shap_interpretation_figure.R`. | Verifiable locally. |
| 5.20 | "STOP window: UGA shifts predictions toward NMD; U at +4 has largest importance for that position (matches readthrough biology)" | `model:07_motif_analysis.py` + driver. **This is the headline biological finding for §5.** | **Verifiable locally — high priority verification target given the lit-review §6 dependence on this claim.** |
| 5.21 | **Figure 6F** — nucleotide-level attribution around stop codon (UGA + U+4 highlighted) | `model:make_shap_interpretation_figure.R`. | Verifiable locally. |

### Paragraph 4 — attention distribution (cited as "Fig X panel A and B")

**Open question (§5 figure-naming):** Paragraph 4 references "Fig X, panel A and B" for the attention-distribution comparison (NMD vs Control). This overlaps with Figure 6 panel A (architecture) and panel B (AUC/AUPRC) cited in Paragraphs 1 and 2. The attention panels are therefore probably part of a **separate figure** (Figure 7?) or supplementary, not Figure 6. Worth pinning down with Pete when fixing manuscript figure numbering.

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 5.22 | "ORFs numbered by CDS likelihood; ORF0 = most likely CDS; most attention on ORF0" | `model:04_interpret_attention.py` — exports per-transcript attention vectors over 5 ORFs. | Verifiable locally. |
| 5.23 | "Attention more broadly distributed in NMD than in Control isoforms" | Same — entropy or top-1-share of attention vector, stratified by NMD label. | Verifiable locally. |
| 5.24 | **"Fig X panel A and B"** (attention distribution NMD vs Control) — likely Figure 7 or supplementary, not Figure 6 | Render likely embedded in `model:make_shap_interpretation_figure.R` or in `model:orf_model_report_v5.Rmd`. | Verifiable locally once panel identity is pinned down. |

### Paragraph 5 — GC content + mechanistic subgroups (SFx)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 5.25 | "Strong GC content difference between NMD and Control isoforms, most prominent after the stop codon — reflects PTC turning exonic sequence into 3′UTR-like" | `model:09_export_gc_content.py` — extracts the GC channel attribution; stratify by NMD label. | Verifiable locally. |
| 5.26 | "GC content channel discriminates NMD vs Control in STOP but not START window" (SFx) | Same — START vs STOP window GC profile comparison. | Verifiable locally. |
| 5.27 | "Mechanistic subgroups (PTC+, PTC− ref-AUG retained, PTC− ref-AUG lost) extended from gene-paired set to full test data" | `model:08_export_subgroup_deepshap_tsv.py` + `09b_export_subgroup_profiles.py`. Subgroup labels traced back to `nmd:results/.../isopair_wrapper/05r_ref_atg_analysis.R` outputs. | **Verifiable locally** (cross-repo dependency between model and Isopair wrapper). |
| 5.28 | "PTC+ NMD relies heavily on EJC count for predictions" — subgroup-stratified DeepSHAP | `model:08_export_subgroup_deepshap_tsv.py`. | Verifiable locally. |
| 5.29 | "PTC− ref-AUG-retained + PTC− ref-AUG-lost subgroups rely primarily on 5′UTR fraction + CDS source" | Same. | Verifiable locally. |
| 5.30 | "START window block-level importance ~3× higher in those two subgroups than in PTC+ group" | Same — KernelSHAP block-level (5.11/5.12 logic) stratified by subgroup. | Verifiable locally. |
| 5.31 | "Overall predictive performance substantially lower in those two subgroups" | `model:evaluate.py` stratified by subgroup. | Verifiable locally. |

### Section 5 figure render-script summary

Local. Figure 6A is `model:make_architecture_figure.R`; the remaining Figure 6 panels (B / C / D / E / F) are produced by `model:make_shap_interpretation_figure.R`. The "Fig X panel A and B" attention panels in Paragraph 4 are TBD on figure numbering but render-locally as well.

### Section 5 verifiable-locally summary

**Essentially all of §5 is locally verifiable**, with these caveats:

- **Training itself requires GPU + HDF5 input.** Pete already has the trained-model weights cached in `model:results_4ct/`; verification of AUC/AUPRC is from running `evaluate.py` against those weights, not retraining.
- **KernelSHAP and DeepSHAP runs already executed.** The TSV exports under `model:results_4ct/` (per the `06_export_deepshap_tsv.py` and `11_kernel_shap_branches.py` outputs) are what feed the figure-render scripts. Verification is reproducing the chain from TSV → figure.
- **Cross-repo subgroup dependency:** the PTC+ / PTC− ref-AUG retained / lost labels come from the Isopair wrapper (`05r_ref_atg_analysis.R`). Verification of subgroup attributions requires that those labels match between the two repos.

Headline biological-claim verification targets:
- **5.9** AUC=0.93 / AUPRC=0.83
- **5.11/5.12** block-level ⅔ structural / ⅓ sequence + STOP 3× START
- **5.14** EJC count #1 structural feature
- **5.20** UGA + U+4 attribution at stop codon (the lit-review §6 dependency)
- **5.23** Attention more broadly distributed for NMD
- **5.27–5.30** Subgroup-stratified shift in which features carry predictive weight

---

## Open questions for Pete / Yul

### Repo / access

1. **Yul's changit repo URL** — is `reyle/nmd_lungcells_2026` correct? Pete's `repjc` account currently gets 404 (security-through-obscurity for non-members) when probing. Either Yul adds Pete as a member, or Pete `scp`s the `final/` dir from Channing.
2. **Long-read pipeline repo (Randell_Lung_Cells_2025)** — is there a changit project, or are the scripts cluster-filesystem-only?

### Pete↔Yul Rmd reconciliation

3. **Intermediate SR DGE Rmd vintage** — what produced the `_2026.3.10.csv` mashr outputs in `shortread_dge/mashr/`? Yul's `2026.5.5` Rmd post-dates them.
4. **LR gene-level DE Rmd** — Section 2 quotes "long-read gene-level analysis tested 19,056 genes," but no LR gene-level Rmd is in the code-map doc. Where does it live?
5. **Confirm Pete's `isocall_limma_dge_fullmodel_2026.3.1.Rmd` is QC-only, not canonical** — current map treats it as a parallel/sanity-check implementation per Pete (2026-06-13). Worth a one-time check that the limma contrasts it produces agree with Yul's.

### Unmapped Rmd / script gaps

6. **DMSO-only one-vs-rest cell-type-marker Rmd** (Methods M5; Section 1 numbers 1.9–1.14) — code path TBD.
7. **Tan et al. mashr reanalysis Rmd** (Methods M9; Section 2 claim 2.18) — code path TBD; likely in Yul's `final/`.
8. **Sequencing-summary aggregation script** (Section 1 claim 1.1: "26 samples / 13M FLNC reads / 341,638,920 total reads") — code path TBD.
9. **Section 5 sweep results table** (claim 5.7: AUG × stop window-size sweep results table) — orchestrator script and recorded results TBD.

### Figure-render-script gaps

10. **Figure 1 panels** (A, C, D, F) — render scripts TBD; nominally Yul-side.
11. **Figure 2 panels** (A, B, C, D) — render scripts TBD; nominally Yul-side.
    Panels 2E and 2F (GPR180 isoform structure + logFC) are rendered via `isopair:R/visualization.R::plotIsoformPair()` — locally available.

### Methods-text / claim verification flags

12. **SR↔LR coverage threshold** (Section 3 claim 3.10: "~40% lacked sufficient LR isoform coverage") — Methods don't define "sufficient coverage" quantitatively. Pin down the threshold from `yul:comparison_analysis.Rmd`.
13. **PCI 63–90% denominator** (Section 3 claim 3.7) — exact denominator (genes with significant unproductive accumulation vs all NMD-isoform genes) and PCI-class threshold worth pinning from `yul:productive_compensation.Rmd`.
14. **"85% explained by PTCs"** (Section 4 claim 4.37) — composite arithmetic; pin the exact union arithmetic from `05_final_report_mashr.Rmd`.

### Figure-numbering inconsistencies

15. **Figure 5 panel-drop** (already on task list as #7) — find/replace pairs drafted; Pete to apply in the Google Doc + decide whether to update the figure source code.
16. **"Fig X panel A/B/C/D" orphan in §4 paragraph 5** — these references (three subgroups + 3′UTR + 5′UTR + uORF panels) don't match the post-renumbering Figure 5 content. Likely need to become Figure 6 (or supplementary). To resolve alongside the broader Fig 3+ rework Pete flagged.
17. **§5 attention panels overlap** (claim 5.24) — paragraph 4 cites "Fig X panel A and B" for attention distribution, but those letters are already used for Figure 6 architecture + AUC panels. Likely Figure 7 or supplementary. To resolve alongside Fig 3+ rework.

---

## Next steps (after this map)

In rough priority order:

1. **Pete to triage** which claims most need fast verification vs which can wait for Yul-repo access. The verifiable-locally claims in §2 (GSEA, NMD-susceptible CSV row counts) and almost all of §4 + §5 can begin verification immediately.
2. **Yul-repo access setup** (option A from earlier: Yul adds `repjc` as a member of her changit project), unblocks §1, most of §2, and all of §3.
3. **Figure rework for Fig 3 onwards** (Pete flagged this — distinct workstream).
4. **Then proceed to formal 5-step verification** per `~/.claude/CLAUDE.md` scientific-report protocol, starting with the §4 + §5 verifiable-locally headline claims (44%/5% PTC enrichment; 88%/72% ref-AUG tracing; AUC=0.93; UGA + U+4 stop-codon attribution).

---

*Draft v1.0 — 2026-06-13. First pass through Methods + Sections 1–5 complete (144 numbered claims). Ready for Pete review.*
