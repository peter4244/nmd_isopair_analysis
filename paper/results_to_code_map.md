# Manuscript results → code map

**Manuscript:** *Long Read RNA Sequencing in Primary Lung Cell Types Reveals Principles of Nonsense-Mediated Decay* (Leshem et al., in prep)

**Purpose.** Identify the code that produced every quantitative claim and figure panel in the manuscript. This is **step 1** of the verification workflow — the map is the trace; verification (Pete's 5-step scientific-report protocol from `~/.claude/CLAUDE.md`) follows once the trace is in place.

**Sources used to build this map.** Live Google Doc manuscript (fetched 2026-06-13; §4 audited and corrected 2026-06-16), the code-map Google Doc shared by Pete, `nmd/ONBOARDING.md`, `nmd/CLAUDE.md`, the local `nmd/code/` mirror, the new canonical Rmd `05_final_report_gencode_scope_2026-07-11.Rmd`, and the two verifiers under `reproducibility/` (37-check `verify_pass7_new_rmd.R` manifest + 57-check `verify_cross_check_new_rmd_vs_figures.R`).

**Status of this document.** Revised **2026-07-11** for the 25% reference-share floor (Family A): §4 scopes and all pinned numbers updated to the floored state (pop_BC 3,009→1,548, n=190→130, n=1,166→819, occult 492→348; PTC enrichment 18×→24×). Prior revision 2026-06-16 against the post-cleanup canonical state. The §4 entries reference the new canonical Rmd (`05_final_report_gencode_scope_2026-07-11.Rmd`); the legacy `05_final_report_mashr.Rmd` is deprecated (banners in place) and retained for reference only. Methods + Sections 1–5 mapped at per-claim granularity; 144 numbered claims. §4 numerical claims pinned to the verifier outputs (both pass).

**Verifiable-locally summary** (claims that can be reproduced end-to-end from `~/claude_projects/` clones on this laptop without cluster or Yul-side access):

| Section | **Verified / total** | Status (updated 2026-07-24) |
|---|---|---|
| §1 (Isoform discovery) | **9 + 2 partial / 15** | Verified 2026-07-21: 1.2, 1.4, 1.5, 1.6–1.8, 1.9, 1.10, 1.12 (counts), 1.14. Partial: 1.3, **1.13** (filter narrowed to per-CT 2026-07-24 — see F-14). 1.11 Yul correcting. Open: 1.1 (Channing-side), 1.15 (Fig 1B render). |
| §2 (NMD response) | **18 / 26** | Verified: 2.1–2.6, 2.9, **2.10–2.11 (unblocked 2026-07-24)**, 2.12–2.15, 2.18, 2.19; 2.21/2.23/2.26 robust. **2.24/2.25 RETIRED — no longer in the manuscript (F-13).** 2.17 is a citation of Tan et al., not our computation. Open: figure claims 2.7, 2.8, 2.16, 2.22 and 2.20 (pathway names). |
| §3 (Output Lost + Productive Response) | **24 / 24** | **REWRITTEN 2026-07-24 — see §3.** ¶1–¶5 all verified (3.1–3.23) by executing Yul's own code. Every computation reproduces; open issues are framing/attribution. 🔴 F-3, F-9, F-11. **All 24 claims now addressed.** All three §3 Rmds **and** their inputs are now LOCAL (`code/upstream/`, `nmd_fig_data/`) — no longer Yul-blocked. The PCI/GPR180 framework (old 3.6–3.26) is retired; canonical source is now `productive_response.Rmd`. ⚠️ finding **F-3** open. |
| §4 (Isopair / splice events / PTC) | 46 / 46 | All local. Verified by the pass-7 (37 checks) + cross-check (57 checks) verifiers under `reproducibility/`. Both PASS. |
| §5 (DL model) | ~30 / 31 | Essentially all local — `model:` (trained weights cached) |

About half of the manuscript's quantitative claims are directly verifiable from this laptop, concentrated heavily in §4 and §5. **§4 is now fully verified** against the verifier suite.

---

## Conventions

### Repo / canonical-source legend

| Tag | Means |
|---|---|
| `nmd:` | This repo. Canonical at GitHub `peter4244/nmd_isopair_analysis` ≡ Channing GitLab `repjc/nmd_lungcells_2026` ≡ local `~/claude_projects/nmd/`. Pete-side. |
| `yul:` | Yul's upstream analysis code. **IMPORTED into this repo 2026-07-20 (Phase 3): `yul:<path>` now resolves to `code/upstream/<path>`.** Source of record = Yul's canonical GitHub `YLeshem18/nmd_lungcells_2026` (per Pete, O-1). The prefix is retained as a provenance marker (Yul owns §1–3 + DGE/DIE/mashr; see `code/upstream/IMPORT_PROVENANCE.md`). Verified: all `yul:` citations resolve under `code/upstream/` **except** `yul:productive_compensation.Rmd`, which is a known-wrong filename (actual: `productive_response.Rmd`) left in place because §3 is parked. |
| `isopair:` | Isopair package. Canonical at GitHub `peter4244/Isopair` ≡ local `~/claude_projects/Isopair/`. |
| `model:` | DL model. Canonical at GitHub `peter4244/NMD_orf_model_v5_4ct` ≡ local `~/claude_projects/NMD_orf_model_v5_4ct/`. Trained on Northeastern Discovery. |
| `randell:` | Long-read isoform-discovery pipeline at Channing `/proj/regeps/regep00/studies/ExternalCellLines/data/longread/mrna/Randell_Lung_Cells_2025/`. Possibly also a changit repo (TBD — Pete to check with Yul). |

### DIE / DGE provenance

Manuscript-quoted DIE/DGE numbers are **mashr posterior estimates** (lfsr, posterior mean logFC), not raw limma adj.P.Val. The "NMD susceptible" definition is `lfsr < 0.05 AND posterior_mean > 0`, per ONBOARDING §6.

**LR DIE provenance (single Yul-side Rmd doing both limma and mashr):**

```
limma + mashr in one Rmd (Yul-side, canonical)         →  canonical CSV (in this repo)
yul:Isoform_Level_Quantification.Rmd                      nmd:isocall_dge/mashr/nmd_mashr_die_{at,dd,fb,mv}_2026.3.10.csv
```

Pete also has a parallel limma implementation at `nmd:code/isocall_limma_dge_fullmodel_2026.3.1.Rmd` and a deleted-2026-06-16 Pete-side run at `nmd:isocall_dge/limma/pete/` (see `isocall_dge/limma/README.md`). Per Pete's confirmed ownership model: **Yul owns limma + mashr DGE; Pete owns Isopair + the deep-learning model.** Yul's Rmd is canonical for the manuscript mashr CSVs; the Pete-side Rmd is a QC/sanity-check parallel implementation that should not be cited as canonical.

**SR DGE provenance (single Yul-side Rmd doing both limma and mashr):**

```
limma + mashr in one Rmd (Yul-side)                    →  canonical CSV (in this repo)
yul:NMD_shortread_dge_fullmodel_2026.5.5.Rmd              nmd:shortread_dge/mashr/nmd_mashr_dge_{at,dd,doali,ddali,fb,mv}_2026.3.10.csv
                                                          nmd:shortread_dge/mashr/mashr_{lfsr,posterior_means,...}_2026.3.10.csv   (4-CT shared)
```

### Outdated material in the code-map Google Doc

Pete flagged "outdated material, especially at the beginning of the document." Concrete corrections applied here:

1. **Long-read quantification used the PacBio Isocall pipeline, not oarfish.** The oarfish minimap2/oarfish_quant pipeline and the downstream `oarfish_gencode49_merged_collapsed_limma_dge_fullmodel_*.Rmd` series are **deprecated**. The Pete-side QC limma is `nmd:code/isocall_limma_dge_fullmodel_2026.3.1.Rmd`; the canonical limma is Yul's (M6 below).
2. SR DGE Rmd path: Pete's code-map doc cites Yul's `2026.5.5` vintage. The mashr CSVs are dated `2026.3.10`. **Open question:** did the `2026.5.5` Rmd re-emit the same `2026.3.10` CSVs (unlikely), or is there an intermediate `~2026.3.10` vintage that wrote the CSVs? Flagged for Yul to confirm.

### Code-map gaps

These manuscript analyses are now covered by added code references; the remaining gaps are Yul-side:

- ✓ Isopair package + wrapper pipeline (§4) — now mapped to `isopair:` package + `05_final_report_gencode_scope_2026-07-11.Rmd` + per-panel scripts (§4 section)
- ✓ DL model training + interpretability (§5) — mapped to `model:` (§5 section)
- ✓ 3′UTR length restricted to GENCODE-annotated CDS (§4 final paragraph) — covered by the `CDSand3UTR_GENCODEonly` supplement (`figures/SupplementalFigures/SF33-SF36_CDSand3UTR_GENCODEonly/`) and Rmd §2d / §3a
- **Remaining gap:** Tan et al. (2025) mashr reanalysis (§2) — Yul-side; flag for Yul to identify the Rmd in `final/`
- **Remaining gap:** DMSO-only one-vs-rest cell-type-marker mashr analysis (Methods M5) — Yul-side; flag for Yul

---

## Methods → code

### M1. Cell culture, treatment, RNA isolation

Wet-lab protocols; not in scope for code mapping. Quoted in Methods sections "Isolation and Culturing of Primary Lung Cell Types" and "RNA Isolation."

### M2. Short-read RNA-seq processing (nf-core/rnaseq)

- **Pipeline:** nf-core/rnaseq v3.14.0 + Nextflow 24.04.4 (per Methods).
- **Code:** run by **John Ziniti** (bioinformatician; co-author) — not Yul, not Pete.
- **Status: ACCEPTED GAP — closed as a decision, not a task (Pete, 2026-07-20).** The
  launcher script/config will **not** be obtained or deposited. Do not re-open this as an
  action item.
- **Why this is defensible:** nf-core/rnaseq is a standardized, versioned community
  pipeline, and the Methods pin **v3.14.0 + Nextflow 24.04.4**. The launcher is
  parameterization, not novel logic, so a reader can reconstruct the run from the Methods
  description + the pinned versions. This differs from bespoke code, where a missing
  script means a missing method.
- **The limitation, stated honestly:** the chain is therefore not *literally* one-command
  runnable from GEO to salmon gene counts — that step must be re-created from the Methods.
  The Methods description carries this load, so it must stay accurate and complete.
- **Outputs feeding downstream analysis:** Salmon-quantified gene/isoform counts via tximport (used only at gene level here; Isocall is canonical for isoform).
- **Where the gene-level counts land:** referenced by `yul:NMD_shortread_dge_fullmodel_2026.5.5.Rmd` (input file path TBD — need to read the Rmd header).

### M3. Long-read RNA-seq processing (PacBio Kinnex → Isocall)

- **Library + sequencing:** PacBio protocol 103-238-700 REV07 (per Methods).
- **Alignment + quantification:** **run by Pete**, using his own Nextflow wrapper —
  **`peter4244/Isocall_v1`** (verified reachable on GitHub 2026-07-20; Channing mirror
  `changit:repjc/Isocall_v1`). Ships `main.nf`, `Dockerfile`, `environment.yml`, `conf/`,
  `modules/`. **M3 code gap CLOSED.** Wraps the PacBio Isocall pipeline; NOT oarfish
  (deprecated per Pete 2026-06-13).
- **⚠ Repo hygiene:** a working copy sits at `nmd:isocall_pipeline/` but is
  **gitignored** (`.gitignore:38`) with **no `.gitmodules`** — an untracked nested repo,
  invisible to this repo and to any snapshot. For the citable deposit it must be **cited
  by URL/DOI** (like Isopair) or made a real submodule; it will not travel on its own.
  A second copy, `isocall_pipeline.channing-backup/` (7.5 M), is also nested + untracked.
- **Where Isocall outputs land:** Channing `/proj/regeps/.../Randell_Lung_Cells_2025/results/`
  (interim products — under D-5 these are **not** deposited).
- **SQANTI3 classification:** config file at `randell:results/sqanti_runs/merged_collapsed/sqanti3_config_cts_2subj_5reads_2025.12.20_merged-collapsed.yaml`.
- **Filtered count matrix consumed by limma:** `nmd:sqanti/nmd_lungcells/results/nmd_lungcells_filtered.count_matrix.txt` (per `isocall_limma_dge_fullmodel_2026.3.1.Rmd` setup chunk).

### M3a. Platform correlation (SR ↔ LR)

Added 2026-07-20. Has a standalone Supplemental Methods subsection ("Platform
Correlation") but previously had no M-section — it was covered only at claim level
(rows 1.6–1.8).

- **Code:** `yul:correlation_analysis.Rmd` (analysis). Panel render is
  `yul:Figures/make_panels.R` → panel 1A (see the §2 figure table) — the analysis Rmd
  does *not* render the figure.
- **Inputs:** `dge_gene_unfiltered_2026.1.2.rds`, `salmon.merged.gene_counts_length_scaled.FULL.rds`.
- **⚠ DO NOT "MODERNIZE" `dge_gene_unfiltered_2026.1.2.rds` — verified 2026-07-21.** It
  looks deprecated by every heuristic (oldest date stamp, loose in the project root while
  everything else is filed under `long_read/`/`short_read/`, and a newer
  `dge_shortread_gene_2026.3.2.rds` exists). **It is required.** In
  `render_sr_lr_correlation.R` it is not the expression source (salmon is) — it supplies
  the **metadata that identifies DD_ALI samples**, which appear in salmon column names
  mislabeled as plain `_DD`. The script rebuilds canonical IDs from this DGEList and drops
  any `_DD` column whose `_DD_ALI` counterpart exists. That works **only because this
  object predates the 4-CT restriction and still contains DD_ALI**.
  **Measured on the real objects:** Jan (6-CT, 38 samples) → **6 DD_ALI columns dropped**;
  Mar (4-CT, 26 samples) → **0 dropped**. Substituting the newer file silently
  contaminates Figure 1A with 6 DD_ALI samples — no error, just wrong numbers.
  Claims 1.6–1.8 depend on this.
- **Claims:** 1.6 (mean Pearson r = 0.90, mean Spearman ρ = 0.890, 26 matched samples,
  29,185 common genes), 1.7, 1.8 (Figure 1A).
- **Status: ✅ VERIFIED LOCALLY 2026-07-21** — re-derived end-to-end from the data bundle,
  independent of Yul's rendered output. Every published value reproduces:

  | Quantity | Published | Re-derived |
  |---|---|---|
  | Common genes (post expressed-filter) | 29,185 | **29,185** |
  | Matched samples | 26 | **26** |
  | Mean Pearson r | 0.90 (0.83–0.91) | **0.896** (0.834–0.906) |
  | Mean Spearman ρ | 0.890 (0.849–0.901) | **0.890** (0.849–0.901) |

  First §1–3 claim verified on this laptop. Also confirms (a) the
  `dge_gene_unfiltered_2026.1.2.rds` input is correct — see the do-not-modernize warning
  above — and (b) `salmon...FULL.rds` is the same object as
  `salmon.merged.gene_counts_length_scaled.rds`, renamed when Yul staged the bundle
  (no missing file).
- Also the basis for excluding DD_ALI / DO_ALI from the 4-CT manuscript scope (low
  SR-vs-LR effect-size correlation).

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
- **LR isoform level (limma + mashr in one Rmd, canonical):** `yul:Isoform_Level_Quantification.Rmd` → CSVs in `nmd:isocall_dge/mashr/`. Per Pete's ownership model, Yul owns this — canonical end-to-end implementation for the manuscript numbers.
- **LR isoform level — parallel Pete-side limma (QC / sanity, NOT canonical for manuscript):** `nmd:code/isocall_limma_dge_fullmodel_2026.3.1.Rmd`. Reads SQANTI-filtered Isocall count matrix; writes limma contrast tables to `nmd:isocall_dge/`. Useful for cross-checking Yul's results but not the source of the mashr CSVs. The previously-tracked `nmd:isocall_dge/limma/pete/` snapshot was deleted on 2026-06-16 (see `isocall_dge/limma/README.md`).
- **Convention:** design = `~ cell_type + treatment + cell_type:treatment`, reference = LAE, donor block via duplicateCorrelation.

### M7. mashr (the shared infrastructure)

- **Method:** estimate_null_correlation_simple on 20k random feature subset → mash_1by1 at lfsr<0.05 for strong signals → cov_canonical + data-driven (PCA, ≤5 PCs) → mash fit.
- **Code:** the mashr step happens inside the Rmds above (M5, M6 SR, M6 LR-DIE). There may also be a shared mashr helper function file — flag for Yul.

### M8. GSEA + pathway enrichment (gene-level only)

- **Code:** `nmd:code/gsea_mashr_2026.3.10.R` (per ONBOARDING §9).
- **Output:** `nmd:tmp/gsea_mashr_gene_2026.3.10_run2026-05-18.tsv` (per ONBOARDING).
- **Method:** fgsea against MSigDB (Hallmark, KEGG, Reactome, GO-BP), genes ranked by signed posterior mean logFC, min=15, max=500, FDR<0.05.

### §1–3 supplemental-figure producers — POINTER (added 2026-07-20)

Yul's 2026-07-20 push (`yul/main` commit `78f1a8c`) added per-figure render scripts for
the §1–3 supplements that previously had no in-repo producer, all wired into
`yul:Figures/make_supplemental_figures.Rmd`. **NOT yet threaded into the SF claim rows
below** — that is Phase 2 of the consolidation plan, done after the Yul import + the
deliverables ledger. Recorded here so the producers are not re-discovered:

| Supplement(s) | Producer (`yul:Figures/`) |
|---|---|
| SR↔LR correlation (Fig 1A + SF) | `render_sr_lr_correlation.R` |
| Isoform length | `render_isoform_length.R` |
| Pairwise expression | `render_pairwise_expression.R` |
| SQANTI3 categories | `render_sqanti_categories.R` |
| mashr sharing | `render_mashr_sharing.R` |
| Proportion vs expression | `render_proportion_vs_expression.R` |
| SF20 % output lost, gene level | `render_output_lost_gene.R` |
| SF21 per-isoform output lost | `render_output_lost_per_isoform.R` |
| SF22 output lost by CPM threshold | `render_output_lost_threshold.R` |

RBP rosters now committed at `yul:data/` (`encode_rbp_roster_vannostrand2020.csv`,
`gerstberger_2014_rbp_census.csv`). Tan tables committed at `yul:tan_reanalysis/data/`
**via Git LFS** — see consolidation plan W-1 (LFS + snapshot) and W-2 (redistribution).

### M8a. RNA-binding protein + SR protein enrichment among NMD targets

Added 2026-07-20. This was the map's only true hole — a full Supplemental Methods
subsection with **no M-section and no claim rows**. Code landed 2026-07-20 (commit
`1fb49d8`). Distinct from M8, which is gene-level GSEA/pathway enrichment.

- **Code:** `yul:nmd_rbp_enrichment.Rmd` (RBP-category enrichment: LR DIE mashr NMD
  calls intersected with the Gerstberger 2014 RBP census) and `yul:rbp_sr.Rmd`
  (SR/hnRNP-focused analysis + the ENCODE eCLIP roster). Figure render is
  `yul:Figures/make_sr_isopair.R` (per-gene SR-protein isoform structure +
  donor-paired logFC), orchestrated by `yul:Figures/make_supplemental_figures.Rmd`.
- **Inputs:** `dge_isoform_longread_filtered_2026.3.3.rds`,
  `mashr_isoform_{posterior_means,lfsr}_2026.3.10.csv`,
  `gerstberger_2014_rbp_census.csv`, `encode_rbp_roster_vannostrand2020.csv`.
- **Method:** NMD target = mashr lfsr < 0.05 & posterior mean > 0; background = all
  genes with an isoform tested in the DIE/mashr set; one-sided Fisher per RBP category.
- **Status:** code resolved. **Two open data items:**
  (a) `gerstberger_2014_rbp_census.csv` is a published census — external download, fine;
  (b) `encode_rbp_roster_vannostrand2020.csv` lives in Yul's `New_NMD_Files/` and
  appears **curated, not a straight download** — confirm with Yul whether it needs
  depositing rather than citing.
- **Path caveat (D3):** both Rmds hardcode `/udd/reyle/nmd_lungcells_2026` and set
  `.libPaths("/udd/reyle/Rlibs")` — needs relative-path rewrite before the citable repo.

### M9. Tan et al. (2025) reanalysis

- **Code:** `yul:tan_reanalysis/tan_transcript_reanalysis.R` (+ `tan_reanalysis/README.md`). Gap CLOSED 2026-07-20 (commit `de08b94`). Converts EBSeq PostFC + PPEE to bhat/shat and fits mashr jointly across the 8 conditions (3 UPF2 + 1 UPF3B × hESC + NPC), then a Tan-style binary NMD-target overlap between hESC and NPC.
- **Inputs:** Tan et al. Supplementary Tables S1, S2, S4, S6 (`.xlsx`). **Not redistributed** — published author data, downloaded from the paper's supplementary material. The per-condition file→sheet manifest is in the script's README.
- **Status:** code resolved; inputs external-but-public (acceptable for the citable repo with a download note).

### M10. Transcriptional output lost + productive isoform response

- **Code:** `yul:transcriptional_output.Rmd` (% output lost) and `yul:productive_response.Rmd` (custom NMD-anchored normalization + Kitagawa decomposition). **Both local** under `code/upstream/`.
- ⚠️ **The PCI framework (`productive_compensation.Rmd`) is RETIRED** — that analysis is no longer in the manuscript (2026-07-24 check: `GPR180`/`compensat*` appear 0 times). Do not cite it as the §3 source.
- **Status:** ¶1 verified exact 2026-07-24 (`code/upstream/verify_section3_p1.R`); ¶2–¶5 recomputable. ⚠️ See §3 finding **F-3** — the "% output lost" metric is provably a total variation distance (symmetric; 21–55% NMD-attributable). Interpretation needs Yul's sign-off before submission.

### M11. Isopair pairs analysis (canonical Rmd 2026-07-11; 25% reference-share floor)

- **Reference-share floor (2026-07-11, Family A):** after `generatePairsExpression` selects the reference (rank-1 by DMSO mean CPM in the strict non-NMD pool), a gene is RETAINED only if `mean_DMSO(reference) / Σ_{i ∈ structures∩expressed isoforms of gene} mean_DMSO(i) ≥ 0.25`, DMSO basis = all_samples (13 DMSO libraries). Injected in `02_build_profiles_mashr.R` (`REF_SHARE_FLOOR <- 0.25`). Effect: pop_BC 3,009→1,548, n=190→130, n=1,166→819, occult 492→348. Rationale + verification: `isopair_wrapper/REFERENCE_FLOOR_{RATIONALE,VERIFICATION,SF_INVENTORY,NUMBERS_DELTA}.md`. The floor is NOT re-anchoring — it drops genes whose truly-non-NMD reference is a minor isoform; it does not change which isoform is the reference for retained genes.
- **Package:** `isopair:` (canonical at `peter4244/Isopair`, vignette `Isopair/vignettes/NMD-attribution.Rmd` is the canonical methods source per ONBOARDING). ⚠ The vignette + manuscript Methods do NOT yet describe the floor — add in Phase 7.
- **Pipeline wrapper (force-tracked in this repo):** `nmd:results/isoform_transitions/Version_6.0/isopair_wrapper/` — `01`–`03b` core + per-feature analysis sub-scripts (`05k_utr5_all_isoforms.R`, `05k_b_utr5_refaug.R`, `05r_ref_atg_analysis.R`). (Legacy-report feeders `04`/`04b`/`05l`/`05s`/`05m`/`05t`/`05u`/`05v` removed 2026-07-11; the legacy reports are kept as frozen `.html`.)
- **Canonical report Rmd (NEW, 2026-06-15):** `05_final_report_gencode_scope_2026-07-11.Rmd`. Structure: §1 (pop_BC + Fig 3 A/B/C), §1a (transcript length), §1b (annotation status), §1c (Fig 3 A/B/C embed), §2 (n=130 strict scope), §2a (PTC determination + Fig 3 D), §2b (PTC-causing event attribution + Fig 3 E/F), §2c (5'UTR length + longest 5'UTR ORF + Fig 4 A/B), §2d (CDS/3'UTR sanity + CDSand3UTR supplement row 1), §3 (n=819 broad ref-AUG-traceable scope), §3a (5'UTR + longest ORF + Fig 4 C/D), §3b (PTC rate at Section C), §4 (TD2 bias: §4a ORF length, §4b paired Kozak PWM, §4c TD2 ATG position, TD2BiasEvidence supplement), §5 (cumulative accounting: 756/819).
- **Legacy Rmd (DEPRECATED, retained for reference):** `05_final_report_mashr.Rmd` — has deprecation banners; numbers in this Rmd are no longer canonical.
- **CDS / 3'UTR analysis on GENCODE-annotated CDS subset (manuscript §4 final paragraph):** Rmd §2d (n=130 strict scope) + §3a (n=819 broad scope). The `CDSand3UTR_GENCODEonly` supplement is a 2×3 layout reproducing the bias-correction pattern across both scopes.
- **Shared helpers:** `nmd:figures/lib/mechanism_class.R` (derived helper; replaces the old cached `mechanism_class` column); `nmd:figures/lib/validate_flowchart_dot.R` (DOT static validator); `nmd:figures/lib/validate_figure_layout.py` (scale-aware per-axis tolerances).
- **Verifiers:** `nmd:reproducibility/verify_pass7_new_rmd.R` (37-check manifest pinning every Rmd-rendered figure number to its canonical source) + `nmd:reproducibility/verify_cross_check_new_rmd_vs_figures.R` (57-check cross-check binding Rmd HTML claims to figure-side TSVs). Both PASS.

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

## Tables → code

Added 2026-07-21. The map previously had **no tables section at all**, which made any
"every table traces to a producer" gate vacuous. Manuscript citation counts: ST4 ×6,
ST3 ×2, Tables 1–4 ×2 each.

**Caption convention (do not mis-read):** in the manuscript the caption sits **BELOW**
its table. Line 187's table (header "Number of *Genes* Tested") is captioned Table 1 at
line 194; line 196's table (header "Number of *Transcripts* Tested") is captioned Table 2
at line 203. Table 3's caption notes "FB only reached significance in 3 pathways" and the
table above it has exactly 3 FB rows. All four pairings are correct.

| Table | Content | Producer | Status |
|---|---|---|---|
| **Table 1** | SR gene DGE per CT: tested / significant / NMD-susceptible / % / mean logFC | `nmd:shortread_dge/mashr/nmd_mashr_dge_{at,dd,fb,mv}_2026.3.10.csv` (M6 → M7) | ✅ **VERIFIED 2026-07-21** — all 4 CTs × 5 columns exact (25,955 tested; AT2 3,638/2,276/62.56%/0.80; LAE 6,753/3,326/49.25%/0.94; FB 3,122/2,066/66.18%/0.80; MV 3,428/2,310/67.39%/0.98) |
| **Table 2** | LR isoform DIE per CT, same columns | `nmd:isocall_dge/mashr/nmd_mashr_die_{at,dd,fb,mv}_2026.3.10.csv` (M6 → M7) | ✅ **VERIFIED 2026-07-21** — all 4 CTs × 5 columns exact (162,800 tested; AT2 24,847/22,850/91.96%/2.09; LAE 35,336/32,239/91.24%/2.97; FB 24,803/22,552/90.92%/1.60; MV 27,834/25,294/90.87%/2.68) |
| **Table 3** | Top-5 **gene-level** significant pathways per CT (NES, p.adj, size) | **`code/upstream/Gene-Level_DGE_Summary_mashR.Rmd`** (msigdbr **26.1.0**, seed 42; fgsea `nPermSimple=100000`, size 15–500; ranked by posterior logFC). Verifier: `code/upstream/verify_table3_gsea.R`. | ✅ **VERIFIED 2026-07-21** against the **updated** manuscript after installing msigdbr 26.1.0. AT2 top-5 reproduce exactly: RNA splicing NES 1.48 / padj 0.0071 / size 462; tRNA modification 1.78 / 0.020 / 96; Vitamin D 1.86 / 0.033 / 18; RNA modification 1.59 / 0.037 / 178; tRNA metabolic 1.54 / 0.048 / 207. NES + sizes exact; padj to 3 dp (permutation jitter). The size that was **449 under 25.1.1 is now 462 under 26.1.0** — the version was the whole cause. |
| **Table 4** | Top-5 **isoform-level** significant pathways per CT | **`code/upstream/die_mashr_enrichment_part2_2026.3.10.Rmd`** (msigdbr 26.1.0, seed 42; per-gene ranking = largest-|logFC| isoform, signed logFC; novel `^ENSG` symbols dropped; fgsea **default** `nPermSimple`). Verifier: `code/upstream/verify_table4_gsea.R`. | ✅ **VERIFIED 2026-07-21** — AT2 top-5 exact on rank order, NES, and size: RNA splicing 1.33/438, via-transesterification 1.34/316, capped-intron pre-mRNA 1.33/291, ribosome biogenesis 1.28/325, translation 1.26/364. padj differs (2.1e-17 vs 9.0e-17 etc.) — expected: default `nPermSimple` gives coarse permutation p-values (fgsea warned); NES + size are the deterministic quantities and match. |
| **ST3** | CT-associated isoform counts (28,930 / 22,131 / 21,429 / 18,422); Gerstberger OR 4.0 / 6.9 | Counts = one-vs-rest mashr (M5). OR = RBP enrichment (**M8a**, `yul:nmd_rbp_enrichment.Rmd`). | Partially mapped; not re-derived. |
| **ST4** | Productive-response direction calls, level-vs-composition split, SR-protein enrichment | `yul:productive_response.Rmd` (M10 area) | **Blocked** — §3 is parked (see M10 note); ST4 is the most-cited table in the paper. |

**Workbook artifact:** `nmd:PJC_Tables/NMD_Tables.xlsx` — 4 sheets (`DGE summary`,
`DIE_summary`, `DGEtoIso`, `IsotoDGE`). Hand-maintained, **no generating script**; the
`DGE summary` / `DIE_summary` sheets correspond to Tables 1–2, whose numbers are now
verified against the mashr CSVs directly. Note
`isopair_wrapper/build_paper_tables.R` builds Isopair's **internal** `table1`–`table13`,
a different numbering — it does **not** produce manuscript Tables 1–4.
Housekeeping: `PJC_Tables/~$NMD_Tables.xlsx` is a tracked Excel lock file → untrack.

## Section 1 (Isoform Discovery via Long-read RNA-sequencing) → code

### Paragraph 1 — sequencing summary + SQANTI categories + isoform-count distribution

| # | Claim (verbatim or near-verbatim) | Source code | Notes / status |
|---|---|---|---|
| 1.1 | "26 samples to an average depth of 13 million aligned FLNC reads, 341,638,920 total reads" | Channing Isoseq pipeline output stats (`randell:results/long-read-isoseq/20251202-early-ALL/`). Aggregation script TBD. | Flag for Yul/Pete to identify the aggregation script. |
| 1.2 | ✅ **VERIFIED 2026-07-21** — 162,800 = `nrow(dge_isoform_longread_filtered_2026.3.3.rds)` exactly; 18,270 = distinct SQANTI `associated_gene` over that filtered set (NOT the DGEList's mapped `gene_id`, which gives 17,382 — both correct, different entities; record this so the 17,382 is not later mistaken for a discrepancy). "162,800 expressed isoforms across 18,270 genes" | Output of SQANTI3 + filterByExpr. Concrete filterByExpr application is in `nmd:code/isocall_limma_dge_fullmodel_2026.3.1.Rmd` (setup chunk → DGEList → filterByExpr). | Verifiable locally. The "18,270 genes" rollup needs an Rmd that aggregates isoform→gene — flag whether this is in Yul's `Isoform_Landscape.Rmd`. |
| 1.3 | ⚠️ **PARTLY VERIFIED 2026-07-21** — "centered near 2.5 kb" confirmed (median 2,576 bp; mean 2,896; IQR 1,716–3,758; density peak ≈2,670 bp). **"Unimodal" NOT tested** — a naive local-maxima count on 162,800 points is bandwidth-noise, not evidence; this is a shape claim to judge from the figure. "Isoform length distributions were unimodal and centered near 2.5 kb" (SF: Isoform Length By Cell Type And Treatment) | `yul:Isoform_Landscape.Rmd` (best guess). | Cite-only. |
| 1.4 | ✅ **VERIFIED 2026-07-21** — all categories exact from `nmd_lungcells_classification.txt` over the filtered set: FSM 55,770 · ISM 31,667 · NIC 44,193 · NNC 28,045 · other = fusion 2,719 + genic 406 = 3,125 (1.9%). "53.7% FSM/ISM (55,770 FSM + 31,667 ISM); 44.4% novel (44,193 NIC + 28,045 NNC); 1.9% other" (SF: SQANTI3 Structural Categories by Cell Type) | SQANTI3 `classification.txt` from `randell:results/sqanti_runs/merged_collapsed/`. Aggregation/percentage script TBD. | Numbers should also be reproducible from the SQANTI classification.txt via a one-off script. |
| 1.5 | ✅ **VERIFIED 2026-07-21** — every statistic exact: median 5 · mean 8.9 · IQR 2–13 · max 97 · ≥2: 13,941 (76.3%) · ≥10: 5,963 (32.6%). Grouped by SQANTI `associated_gene` (same key as 1.2). "median gene expressed five isoforms (mean 8.9, IQR 2–13, max 97), 13,941 (76.3%) genes ≥2 isoforms, 5,963 (32.6%) ≥10" | Isoform-per-gene aggregation. Likely `yul:Isoform_Landscape.Rmd`. | Cite-only. |

### Paragraph 2 — short-read vs long-read concordance (Figure 1A + SFx)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 1.6 | "mean Pearson r = 0.90 (range 0.83–0.91), mean Spearman ρ = 0.890 (range 0.849–0.901), 26 matched samples, 29,185 common genes" | `yul:correlation_analysis.Rmd` (per code-map doc). | Cite-only. |
| 1.7 | "Per cell type Pearson means 0.880 (FB) — 0.901 (AT2)" | Same. | Cite-only. |
| 1.8 | **Figure 1A** (sample-wise SR↔LR correlation plot) | Analysis `yul:correlation_analysis.Rmd`; render `yul:Figures/make_panels.R` → `fig_panels/figure_composite.py` (Methods M3a). | Cite-only. Render script RESOLVED 2026-07-20. |

### Paragraph 3 — cell-type-specific expression (Figure 1B + SFx)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 1.9 | ✅ **VERIFIED 2026-07-21 — EXACT (105,938).** Reproduced by replicating the filter in `code/upstream/ct_de.Rmd:172`: subset to DMSO + the 4 kept CTs (13 of 26 samples), then `filterByExpr(y, design = model.matrix(~0+ct), min.count = 5, min.total.count = 10)` → **105,938 of 162,800**. The non-default arguments are the whole story: my earlier guesses (bare `filterByExpr` → 13,263; `group=ct` → 68,933) all missed because they omitted `design=` and the count thresholds. The claim was correct throughout. Producer confirmed by Yul. Its inputs — `dge_isoform_longread_filtered_2026.3.3.rds` and `dge_shortread_gene_filtered_2026.3.2.rds` — are **both available locally**, so this claim is now **re-runnable rather than blocked**. Under standard D-5 (code + raw reads, no interim files) the chain is complete: the code is in-repo even though its outputs are not. *(Superseded diagnosis below.)* ⛔ M5 output not in repo: `105,938` appears in the manuscript, supplement, and this map but in **no data file** anywhere locally; the DMSO-only one-vs-rest mashr (M5) is Yul-side. My own `filterByExpr` on the DMSO subset gives 68,933 (group=ct) or 13,263 (no group) — neither is the published figure, and without the M5 code path the exact filter parameters are unknown, so no discrepancy should be inferred. "105,938 isoforms expressed at testable levels under DMSO" | DMSO-only filterByExpr re-application step in the one-vs-rest mashr analysis (M5). Code path TBD — flag for Yul. | Cite-only. |
| 1.10 | ✅ **VERIFIED 2026-07-21 — EXACT (all four cell types).** The counts are the **mashr marker rule from the Supplemental Methods: `lfsr < 0.05 AND posterior mean log2FC > 1`**, applied to the one-vs-rest mashr fit (seed 42, `estimate_null_correlation_simple` on a 20k subset, `cov_pca` + `cov_canonical`). Reproduces AT2 22,131 / LAE 28,930 / FB 18,422 / MV 21,429 — **exact**. Full chain also exact end-to-end: DMSO filter → 105,938, donor correlation → 0.026, then this mashr criterion. **Correction to my earlier reasoning:** the manuscript prose says "significantly associated", which I mapped to the limma rule (BH p<0.05 & coef>0 → 25,438, ~13% low). That was wrong — the number is the mashr `pm>1` marker count (the map correctly called 1.10 a mashr output). I had tested `pm>0` (36,433, too high) but not the documented `pm>1`. Verifier: `code/upstream/verify_1_10_onevsrest.R`. **Criterion resolved from Supplemental Methods** ("Differential Gene and Transcript Expression"/"Cell Type Specificity"): "significantly associated with a given cell type" = **BH-adjusted p < 0.05 AND positive contrast coefficient** (the mashr `lfsr<0.05 & pm>1` rule is a *different* quantity — "cell type markers"). **Two intermediate checkpoints reproduce EXACTLY:** the DMSO filter → **105,938** (the Methods' own stated number) and the within-donor correlation → **0.026** (also stated). Applying the exact criterion to that fit gives **AT2 17,270 / LAE 25,438 / FB 15,086 / MV 16,738** vs claimed **22,131 / 28,930 / 18,422 / 21,429** — order exact, magnitudes a uniform ~13% low. Tested and ruled out: single vs iterative voom (25,438→25,459, negligible); both-directions counting (46,604, too high); mashr lfsr<0.05 (60k) and lfsr<0.05&pm>0 (36,433). **Because every documented parameter AND two exact intermediates match, this is the one claim this session where matching the Methods does not reproduce the number** — unlike the six earlier "discrepancies" that all resolved to a wrong input/parameter on my side. Most likely an undocumented step or the counts predating a pipeline change. Verifier: `code/upstream/verify_1_10_onevsrest.R`. **→ Genuine question for Yul, now sharply bounded.** Replicated `ct_de.Rmd` end-to-end and **every checkable step matches**: filter → **105,938 exactly** (= claim 1.9); donor column = `id` (confirmed correct — `get_donor()` at line 97 prefers `id`); design, one-vs-rest contrasts, and the mashr recipe (seed 42, `n_random=20000`, `estimate_null_correlation_simple`, `cov_pca` + `cov_canonical`) all as written. **But no criterion in the file reproduces the counts,** and the two *directional* variants bracket the claim:

| Criterion | AT2 | LAE | FB | MV | Order? |
|---|---|---|---|---|---|
| **Claim** | **22,131** | **28,930** | **18,422** | **21,429** | — |
| limma `adj.P<0.05` | 38,647 | 46,604 | 40,907 | 35,437 | ✗ |
| limma `adj.P<0.05 & logFC>1` | 16,871 | 24,611 | 14,498 | 16,187 | ✓ (~24% low) |
| limma `adj.P<0.05 & \|logFC\|>1` | 37,560 | 45,028 | 39,677 | 34,708 | ✗ |
| mashr `lfsr<0.05` *(the script's own reporting, line 376)* | 60,677 | 64,631 | 61,993 | 58,177 | ✗ |
| mashr `lfsr<0.05 & pm>0` | 28,067 | 36,433 | 23,323 | 27,501 | ✓ (~26% high) |

Only the two **directional** rules reproduce the claim's rank order (LAE>AT2>MV>FB), so the real criterion is one-sided — but it sits *between* them. **Deliberately not searching for a threshold that fits**: with the target known, that would be fitting the answer rather than verifying it. **→ Ask Yul: what exact rule produced 22,131 / 28,930 / 18,422 / 21,429?** Producer confirmed (`ct_de.Rmd`); only the counting rule is open. *(Detail of the first pass below.)* Replicated `ct_de.Rmd`'s pipeline end-to-end (DMSO + 4 CTs → `filterByExpr(design, 5, 10)` → voom → `duplicateCorrelation(block=donor)`, consensus ρ = 0.026 → one-vs-rest contrasts). The filter step reproduces **105,938 exactly** (claim 1.9), so the subset and model are right. But counting `adj.P.Val < 0.05` gives **AT2 38,647 / LAE 46,604 / FB 40,907 / MV 35,437** vs the claimed AT2 22,131 / LAE 28,930 / FB 18,422 / MV 21,429 — higher, **and in a different rank order** (mine LAE>FB>AT2>MV; claim LAE>AT2>MV>FB). A changed ranking means a different *criterion*, not a scale factor. `ct_de.Rmd` offers two others: **mashr `lfsr < 0.05`** (`LFSR_CUT`, and the map itself calls 1.10 a mashr output) and the marker rule **`adj.P.Val < 0.05 & |logFC| > 1`** (`LOGFC_MARKER = 1`). Testing those next; the mashr variant requires the joint fit. Producer confirmed by Yul. **1.12 and 1.13** (cell-type-restricted counts) come from `code/upstream/Isoform_Landscape.Rmd`. Outputs are not in the repo, but both scripts are, and 1.9/1.10's inputs are local → **re-runnable**. "28,930 LAE, 22,131 AT2, 21,429 MV, 18,422 FB significant at 5% FDR" | One-vs-rest mashr output (M5). | Cite-only. |
| 1.11 | 🔧 **CONFIRMED ISSUE — YUL IS CORRECTING (2026-07-21).** She replied: "I will look into it, rerun the analysis, and update the text with correct numbers." So the AT2/LAE range in the manuscript is being revised; **re-verify against the updated text**. This is the **second manuscript correction surfaced by verification** (after 9,161→9,154). Detail of the original finding: ⚠️ **PARTIALLY VERIFIED 2026-07-21** — using the exact definition from `code/upstream/Figures/render_pairwise_expression.R` (mean **log**-CPM, `prior.count=1`, DMSO only, 4 CTs): **FB↔MV = 0.796 ≈ 0.79 ✓**. But the AT2/LAE-involving pairs come out **0.508–0.599** (AT2–FB 0.563, AT2–MV 0.599, LAE–FB 0.508, LAE–MV 0.521) vs the claimed **0.43–0.56** — the claimed lower bound is below anything the 4-CT set produces. Hypothesis (untested): the SF included DD_ALI/DO_ALI among "others", which would extend the range down; no 6-CT isoform object exists locally to check. **Do not treat as a discrepancy until the CT set behind the figure is confirmed.** "Spearman ρ = 0.79 (FB↔MV); 0.43–0.56 between AT2/LAE and others" (SF: Pairwise Expression Similarity) | `yul:Isoform_Landscape.Rmd` (best guess). | Cite-only. |
| 1.12 | ✅ **COUNTS VERIFIED EXACT 2026-07-21** — all four reproduce: **LAE 8,087 · AT2 2,267 · MV 2,015 · FB 1,131**, using the Methods definition (expressed in a CT = summed count ≥1 across that CT's DMSO samples; restricted = detected in exactly one of the four). ⚠️ **Percentages not reproduced:** against each CT's expressed-isoform total I get **0.88–5.48%** (FB 0.88, MV 1.50, AT2 1.70, LAE 5.48) vs the claimed **0.36–2.24%** — a different denominator, not yet identified (tested: per-CT expressed, and all 162,800). Counts are the substantive claim and are exact. "LAE 8,087; AT2 2,267; MV 2,015; FB 1,131 cell-type-restricted isoforms; 0.36–2.24% of each CT's expressed isoforms" | `yul:Isoform_Landscape.Rmd`. | Cite-only. |
| 1.13 | ⚠️ **FILTER NARROWED 2026-07-24 (`verify_section1_claim113.R`) — it is PER-CELL-TYPE.** Of nine candidates, per-CT `filterByExpr` is the only one reproducing the claim's ordering (LAE lowest 25.6, MV highest 45.0 vs claimed 27.8 / 51.7); every global filter gets the ordering wrong. Values remain 2–7 points low, so the exact definition is still needed — see **F-14**. 🔴 **Worse: the sentence's two halves are mutually inconsistent.** "Among those that did all reached significance" caps LAE significance at 27.8% (2,248 of 8,087), but **4,499 (55.6%)** of LAE-restricted isoforms are `nmd_responsive` in the LAE DIE — twice the cap. All 8,087 are also present in the DIE table, so the "filter required for differential testing" did not exclude them from the mashr DIE. *(Superseded diagnosis 2026-07-21:* — taking "the additional expression filter" as the DMSO `filterByExpr` set (105,938) gives **LAE 19.5% → MV 31.4%** vs claimed **27.8% → 51.7%**. The *structure* reproduces (LAE lowest, MV highest, same ordering), so the restricted sets are right (1.12 counts are exact) but the second filter is something else. Not a discrepancy — the filter simply is not specified in the Methods. "27.8% (LAE) to 51.7% (MV) of CT-restricted isoforms passed the additional expression filter" | Same. | Cite-only. |
| 1.14 | ✅ **VERIFIED 2026-07-21** — 0.863–0.923 (FB↔MV 0.923 high, AT2↔FB 0.863 low) vs claimed 0.86–0.92. Definition taken from `code/upstream/Isoform_Landscape.Rmd:148` (`get_expressed` = **≥1 count**, not a CPM cutoff), over the filtered 162,800 set using **all samples, both treatments**. Note DMSO-only gives 0.776–0.852 and does *not* match — the claim is all-sample. "Pairwise Jaccard indices 0.86–0.92" | Same. | Cite-only. |
| 1.15 | **Figure 1B** (expressed + restricted isoform counts per CT) | `yul:Figures/make_panels.R` → `fig_panels/figure_composite.py`; input `dge_isoform_longread_2026.3.3.rds`. NOT `yul:Isoform_Landscape.Rmd` (that Rmd holds the §1 landscape analysis, not the panel render). | Cite-only. Render script RESOLVED 2026-07-20. |

---

## Section 2 (NMD Response) → code

### Paragraph 1 — per-CT DGE/DIE counts + cross-CT sharing of NMD-susceptible features (Figure 1C, 1D)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 2.1 | ✅ **VERIFIED 2026-07-21** — 25,955 exact, all 4 CTs. "25,955 genes per CT" tested at SR level | `yul:NMD_shortread_dge_fullmodel_2026.5.5.Rmd` filterByExpr step → `nmd:shortread_dge/mashr/nmd_mashr_dge_*_2026.3.10.csv` (row count). | Cite-only; CSV row count verifiable locally. |
| 2.2 | ✅ **VERIFIED 2026-07-21** — exact — FB 3,122 · AT2 3,638 · MV 3,428 · LAE 6,753; %NMD-of-sig 49.25–67.39% (AT2 62.56%, LAE 49.25% both exact). "3,122–6,753 significant genes per CT (49–67% NMD susceptible)" — see also Short-Read mashr table (AT2 3,638 / 62.56%; LAE 6,753 / 49.25%; FB 3,122 / 66.18%; MV 3,428 / 67.39%) | Same Rmd; subsetting on `adj.P.Val < 0.05` for "significant" and `lfsr < 0.05 & posterior_mean > 0` for "NMD susceptible". | Verifiable locally by reading the 4 per-CT mashr CSVs. **Watch:** "Significant" here is the *limma adj.P.Val* denominator, not mashr lfsr — the mixed-denominator phrasing is fine but verification needs to use the right column. |
| 2.3 | ✅ **VERIFIED 2026-07-21** — 162,800 exact, all 4 CTs (= the filtered set of claim 1.2). "162,800 isoforms per CT" tested at LR level | `yul:Isoform_Level_Quantification.Rmd` (filterByExpr) → `nmd:isocall_dge/mashr/nmd_mashr_die_*_2026.3.10.csv`. | Cite-only; CSV row count verifiable locally. |
| 2.4 | ✅ **VERIFIED 2026-07-21** — exact — FB 24,803 · AT2 24,847 · MV 27,834 · LAE 35,336; %NMD-of-sig 90.87–91.96%. "24,803–35,336 significant isoforms per CT (90.9–92.0% NMD susceptible)" — see DIE mashr table (AT2 24,847 / 91.96%; LAE 35,336 / 91.24%; FB 24,803 / 90.92%; MV 27,834 / 90.87%; logFC 1.60–2.97) | Same Rmd; same dual-denominator pattern. | Verifiable locally from CSVs. |
| 2.5 | ✅ **VERIFIED 2026-07-21** — 34,387 exact (union of per-CT nmd_responsive sets). "34,387 unique NMD susceptible isoforms across CTs" | Likely `yul:interpret_isoform_patterns_mashr_2026.3.10.Rmd` — union over per-CT NMD-susceptible sets. | Cite-only. Verifiable locally by union over the 4 mashr CSVs. |
| 2.6 | ⚠️ **VERIFIED with one small correction 2026-07-21** — 19,803 shared-across-all-4 exact (57.6%); LAE = 83.4% of CT-specific exact. **BUT CT-specific is 9,154, not 9,161.** The re-derived partition is complete and sums exactly to the union: 9,154 + 1,721 (2 CT) + 3,709 (3 CT) + 19,803 (4 CT) = 34,387 = claim 2.5. The published 9,161 is arithmetically inconsistent with the paper's own 34,387 and 19,803 (it would total 34,394). Percentage unaffected — both round to 26.6%. **Manuscript fix: 9,161 → 9,154.** "19,803 (57.6%) core targets shared across all CTs; 9,161 (26.6%) CT-specific (LAE 83.4% of those)" | Same Rmd — intersect/setdiff over per-CT NMD-susceptible sets. | Cite-only; locally verifiable from CSVs. |
| 2.7 | **Figure 1C** — volcano plots per CT for LR DIE (logFC vs −log10 lfsr, color by NMD-susceptible threshold) | `yul:Figures/make_panels.R` → `fig_panels/figure_composite.py`. Inputs are `nmd:isocall_dge/mashr/nmd_mashr_die_{ct}_2026.3.10.csv`. | Cite-only. Render script RESOLVED 2026-07-20. |
| 2.8 | **Figure 1D** — distribution of posterior mean logFC for NMD-susceptible features at SR-gene vs LR-isoform | `yul:Figures/make_panels.R` → `fig_panels/figure_composite.py`. Inputs are both sets of mashr CSVs. | Cite-only. |
| 2.9 | ✅ **VERIFIED 2026-07-21** — confirmed quantitatively: SR gene 0.80–0.98 vs LR isoform 1.60–2.97 → **2.0–3.2× larger** (AT2 0.80/2.09, LAE 0.94/2.97, FB 0.80/1.60, MV 0.98/2.68). "mean posterior logFC substantially larger at LR isoform than SR gene level" | Same Rmd as 2.8. | Cite-only. |

### Paragraph 2 — NMD-susceptible isoform proportion + expression distribution (SFx panels)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 2.10 | ✅ **VERIFIED 2026-07-24** (`code/upstream/verify_section2_claims_210_211.R`, running Yul's `interpret_isoform_patterns_mashr_2026.3.10.Rmd`). The 2026-07-21 block is **resolved**: that attempt joined a stale oarfish-era table covering <17%; this run covers **92.1%** of NMD isoform-CT pairs. "Concentrated at low isoform proportions" — median proportion **0.1%**, **99.0–99.2%** below 50%, **96.0–96.8%** below 10%. "More than four orders of magnitude" — log10 span **4.21 (AT2) / 4.53 (FB) / 4.94 (LAE) / 4.48 (MV)**, all >4. ⚠️ minor: the central 98% spans only ~2.5–2.8 orders; the >4 comes from min–max extremes. | `yul:interpret_isoform_patterns_mashr_2026.3.10.Rmd` chunks `q3-baseline`, `q3-all-iso`, `q3-nmd-dmso` | ✅ |
| 2.11 | ✅ **VERIFIED 2026-07-24** — "clear upward shift in both proportion and absolute expression". Proportion: median **0.1% → 1.0–1.6%**, increased in **96.8–97.9%** of isoforms. Expression: median CPM **0.06–0.08 → 0.52–0.85**, increased in **99.7–100%**. Paired Wilcoxon p ≈ 0 for both, all 4 CTs. | Same Rmd, chunk `q3-dmso-vs-smg1i` | ✅ |

### Paragraph 3 — quantitative CT specificity vs pairwise sharing (Figure 1E, SF: Sharing)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 2.12 | ✅ **VERIFIED 2026-07-21** — exact — LAE 34.2% (1,138); AT2 1.9% (43), FB 2.1% (44), MV 5.3% (122). Per-CT specificity at SR-gene level: LAE 34.2%, AT2/FB/MV 1.9–5.3% | `yul:interpret_isoform_patterns_mashr_2026.3.10.Rmd` (or sibling) — fraction of NMD-susceptible significant in only that CT. | Cite-only. |
| 2.13 | ✅ **VERIFIED 2026-07-21** — exact — LAE 23.7% (7,637); FB 0.8% (179), AT2 2.1% (488), MV 3.4% (850). Per-CT specificity at LR-isoform level: LAE 23.7%, AT2/FB/MV 0.8–3.4% | Same Rmd. | Cite-only. |
| 2.14 | ✅ **VERIFIED 2026-07-21** — exact from `nmd:shortread_dge/mashr/mashr_sharing_2fold_2026.3.10.csv` — LAE-involving 51.5 / 66.8 / 69.8% (= 51–70%); trio 83.2 / 87.3 / 90.5% (= 83–90%; 90.497 rounds to 90). Pairwise sharing at SR-gene level: LAE-involving 51–70%; AT2/FB/MV trio 83–90% | Same Rmd — fraction with concordant direction + within 2-fold magnitude on mashr posteriors. | Cite-only. |
| 2.15 | ✅ **VERIFIED 2026-07-21** — exact from `nmd:isocall_dge/mashr/mashr_isoform_sharing_2fold_2026.3.10.csv` — LAE-involving 34.6 / 62.3 / 72.4% (= 35–72%); trio 68.1 / 80.3 / 83.0% (= 68–83%). Pairwise sharing at LR-isoform level: LAE-involving 35–72%; AT2/FB/MV trio 68–83% | Same Rmd. | Cite-only. |
| 2.16 | **Figure 1E** — fraction of NMD-susceptible genes + isoforms classified as CT-specific (bar / dot plot) | `yul:Figures/make_panels.R` → `fig_panels/figure_composite.py`. | Cite-only. |

### Paragraph 4 — Tan et al. (2025) reanalysis

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 2.17 | ℹ️ **CITATION OF TAN ET AL., NOT OUR COMPUTATION (clarified 2026-07-21).** The "16–25%" is just the two ratios of Tan's own figures: 176/709 = 24.8% and 176/1,088 = 16.2%. These are **not reproducible from the supplementary tables via our pipeline** and should not be treated as a gap in our code: applying Tan's stated cutoffs (PPDE>0.95, PostFC>2) to their transcript tables gives far larger sets — UPF2_KD ESC 1,900 / NPC 2,683 (163 shared), UPF2_iKO 2,729 / 5,207 (322), UPF2_iKO_KD 8,667 / 11,540 (1,694) — so their 709/1,088 are presumably gene-level or carry additional filters. **Verifiable only against the Tan publication.** Tan-reported overlap: 16–25% UPF2-dependent NMD-target transcripts shared between hESC and NPC; "176 shared of 709 hESC and 1,088 NPC" | Pulled from Tan et al. 2025 Supplementary Tables S1, S2, S4, S6. Source: their paper, not our code. | External — verifiable from their supp tables. |
| 2.18 | ✅ **VERIFIED END-TO-END 2026-07-21** — reproduced **exactly**: shared **3,069**, **49.8%** of hESC, **40.6%** of NPC (hESC UPF2 union 6,168; NPC 7,556; Jaccard 0.288). Run independently on this laptop from the raw Tan supplementary `.xlsx` through `code/upstream/tan_reanalysis/tan_transcript_reanalysis.R`: 8-condition transcript-level mashr over 41,595 transcripts common to all conditions (data-driven PCA + canonical covariances, null correlation from a seeded 20k subset), NMD-target call = posterior mean > 0 & lfsr < 0.05, UPF2 union over KD / iKO / iKO+KD per lineage. The strongest verification in the map: external raw inputs → full model refit → exact integer match. "transcript-level overlap of UPF2-dependent NMD targets between hESCs and NPCs increased to 41–50% (3,069 shared transcripts; 49.8% of hESC and 40.6% of NPC targets)" — concordant direction across both CTs for large majority | Methods M9 describes the procedure (EBSeq PostFC → bhat, PPEE → shat, mashr fit jointly across 8 conditions, lfsr < 0.05 + posterior mean > 0 in same direction). **Code path: gap.** Likely in Yul's `final/` — file probably named `tan*.Rmd` or `wilkinson*.Rmd`. | **Open. Flag for Yul.** |

### Paragraph 5 — GSEA + pathway enrichments (Figure 1F)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 2.19 | ✅ **PROVENANCE RESOLVED 2026-07-21 (Yul).** Tables 3–4 come from **her Rmds**, not from Pete's `nmd:code/gsea_mashr_2026.3.10.R`: Table 3 ← `code/upstream/Gene-Level_DGE_Summary_mashR.Rmd`, Table 4 ← `code/upstream/die_mashr_enrichment_part2_2026.3.10.Rmd`, both with **msigdbr 26.1.0, seed 42**. **This explains the mismatch I found:** local msigdbr is **25.1.1**, and a different MSigDB release changes gene-set membership — precisely the 449 vs 462 size difference. ⚠️ **Yul re-ran with the seed and UPDATED Table 3 in the manuscript**, so the published values have changed; re-verify against the current Doc. Reproducing her GSEA locally requires installing **msigdbr 26.1.0** (env gap — see `ENVIRONMENT.md`). *(Superseded note: the earlier BLOCKED status here was correct to withhold judgement — the input really was wrong.)* "5 pathways significant in 3 of 4 CTs (LAE, FB, MV; FDR 0.002–0.048); AT2 trending (NES 1.65–1.83 below FDR)" — Table X — Pathways Significant | `nmd:code/gsea_mashr_2026.3.10.R` → output TSV `nmd:tmp/gsea_mashr_gene_2026.3.10_run2026-05-18.tsv` (per ONBOARDING §9). Input: signed mashr posterior-mean logFC ranking per CT from `nmd:shortread_dge/mashr/nmd_mashr_dge_*_2026.3.10.csv`. | **Verifiable locally** — both code and output are in the repo. |
| 2.20 | Pathways named: cellular response to topologically incorrect protein, cellular response to unfolded protein, response to topologically incorrect protein, intrinsic apoptotic signaling in response to ER stress, Reactome unfolded protein response | Same GSEA output TSV. | Verifiable locally. |
| 2.21 | ✅ **ROBUST 2026-07-21** — all five present in the UPR/topologically-incorrect leading edge (95 genes total): ATF4, DDIT3, PPP1R15A, ATF3, CHAC1. *(Checked against `tmp/gsea_mashr_gene_2026.3.10_run2026-05-18.tsv`, which is NOT the run behind Tables 3–4 (see 2.19). Treat as a robustness signal across runs, not a definitive verification.)* "Leading-edge genes converging on ATF4, DDIT3, PPP1R15A, ATF3, CHAC1" | Same GSEA output (leading-edge column from fgsea). | Verifiable locally — need to confirm fgsea was called with leading-edge return. |
| 2.22 | **Figure 1F** — KEGG Protein Processing in ER pathway enrichment for LAE; corresponding panels for other CTs in SFx | `yul:Figures/make_pathway_allct.R` (renders hsa04141 per CT) → `make_panels.R` → `fig_panels/figure_composite.py`. | Cite-only. |

### Paragraph 6 — CT-specific pathway enrichments (Table X — Top Pathways)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 2.23 | ✅ **ROBUST 2026-07-21** — both recovered as LAE-only significant: `GOBP_CILIUM_MOVEMENT`, `GOBP_CILIUM_OR_FLAGELLUM_DEPENDENT_CELL_MOTILITY` (LAE has 59 CT-specific pathways). *(Checked against `tmp/gsea_mashr_gene_2026.3.10_run2026-05-18.tsv`, which is NOT the run behind Tables 3–4 (see 2.19). Treat as a robustness signal across runs, not a definitive verification.)* LAE-specific: cilium movement, cilium/flagellum-dependent cell motility | Same GSEA output TSV; LAE-only enriched pathways. | Verifiable locally. |
| 2.24 | ⚪ **RETIRED 2026-07-24 — not in the manuscript.** "xenobiotic" / "translation-regulatory" appear 0 times in both the 2026.7.17 and 2026.2.5 versions; the CT-specific pathway specifics moved to "(supplemental results)", a document not on disk. Do NOT re-run GSEA for this. See F-13. | — | retired |
| 2.25 | ⚪ **RETIRED 2026-07-24 — not in the manuscript.** The one "regulation of RNA splicing" hit is a general Discussion statement, not the AT2+MV cell-type-specific finding. Same disposition as 2.24. See F-13. | — | retired |
| 2.26 | ✅ **ROBUST 2026-07-21** — confirmed null: CT-specific counts are LAE 59, AT2 12, MV 5, **FB 0**. *(Checked against `tmp/gsea_mashr_gene_2026.3.10_run2026-05-18.tsv`, which is NOT the run behind Tables 3–4 (see 2.19). Treat as a robustness signal across runs, not a definitive verification.)* FB: no CT-specific pathways | Same — null result, recompute. | Verifiable locally. |

### Section 2 figure render scripts — CLOSED (2026-07-20)

Resolved by `yul:Figures/` (pushed 2026-07-20, commit `de08b94`). Figure 1 is built
in two stages: per-panel R scripts write PNGs into `fig_panels/`, then a matplotlib
stitcher composes them. Orchestrated by `yul:Figures/make_main_figures.Rmd`, which
carries the authoritative panel→script→input table.

| Panel | Content | Render script | Data inputs (`nmd_fig_data/`) |
|---|---|---|---|
| 1A | SR↔LR gene-expression correlation | `yul:Figures/make_panels.R` | `dge_gene_unfiltered_2026.1.2.rds`, `salmon.merged.gene_counts_length_scaled.rds` |
| 1B | Expressed + CT-restricted isoforms (DMSO) | `yul:Figures/make_panels.R` | `dge_isoform_longread_2026.3.3.rds` |
| 1C | LR isoform DIE volcanoes, 4 CTs | `yul:Figures/make_panels.R` | `nmd_mashr_die_{at2,dd,fb,mv}_2026.3.10.csv` |
| 1D | Posterior-mean logFC distribution, SR gene vs LR isoform | `yul:Figures/make_panels.R` | `nmd_mashr_dge_*`, `nmd_mashr_die_*` |
| 1E | Pairwise sharing of NMD effects, SR genes vs LR isoforms | `yul:Figures/make_panels.R` | mashr CSVs |
| 1F | KEGG hsa04141 "protein processing in ER", LAE | `yul:Figures/make_pathway_allct.R` → `make_panels.R` | `nmd_mashr_dge_dd_2026.3.10.csv` |
| — | Composition of A–F into Figure 1 | `yul:Figures/fig_panels/figure_composite.py` | panel PNGs from the above |

**Corrections to the prior guesses** (all four were wrong — do not reinstate): 1A's
*render* is `make_panels.R`, not `correlation_analysis.Rmd` (that Rmd holds the §1
correlation *analysis*, not the panel); 1B is `make_panels.R`, not
`Isoform_Landscape.Rmd`; 1E is `make_panels.R`, not
`interpret_isoform_patterns_mashr_2026.3.10.Rmd`; 1F is `make_pathway_allct.R`, not
`Gene-Level_DGE_Summary_mashR.Rmd` / `gsea_mashr_2026.3.10.R`.

**Open (data, not code):** all six panels read a `nmd_fig_data/` bundle that is not in
either repo (see `yul:Figures/README.md`) — this is deposit workstream D2.
`make_pathway_allct.R` additionally fetches the KEGG `hsa04141` template at run time
if it is not cached in `fig_panels/`.

### Section 2 verifiable-locally summary

These claims can be verified end-to-end on this laptop without cluster or Yul-side access:

- 2.2 (SR sig + NMD-susceptible counts) — directly from `shortread_dge/mashr/*_2026.3.10.csv`
- 2.4 (LR sig + NMD-susceptible counts) — directly from `isocall_dge/mashr/*_2026.3.10.csv`
- 2.5–2.6 (unique / shared / CT-specific isoforms) — same CSVs, recompute the union/intersection
- 2.19–2.21, 2.23–2.26 (GSEA results) — `code/gsea_mashr_2026.3.10.R` + `tmp/gsea_mashr_gene_2026.3.10_run2026-05-18.tsv`

These are the most efficient first targets when we move from mapping into actual verification.

---

## Section 3 (Transcriptional Output Lost + Productive Isoform Response) → code

> **⚠ REWRITTEN 2026-07-24.** This section previously mapped a **PCI (productive
> compensation index)** framework with a **GPR180** worked example (old claim IDs 3.6–3.26).
> **That analysis is no longer in the manuscript.** Checked against
> `paper/NMD manuscript 2026.7.17.md`: `GPR180` and `compensat*` appear **0 times**, and
> `PCI` appears only inside base64 image blobs. §3 was replaced by a **custom NMD-anchored
> normalization + Kitagawa decomposition** analysis whose canonical source is
> `yul:productive_response.Rmd` — **not** the `productive_compensation.Rmd` cited
> previously. Old IDs 3.6–3.26 are retired; the tables below renumber against the current
> text. Do not re-verify the retired claims.
>
> **Access status changed:** all three §3 Rmds are now local under `code/upstream/`
> (`transcriptional_output.Rmd`, `productive_response.Rmd`, `comparison_analysis.Rmd`) and
> every data input they read is on disk (`nmd_fig_data/`). **§3 is recomputable — no longer
> cite-only.**

### Paragraph 1 — % transcriptional output lost (Fig 2A, SF20, SF21) — ✅ VERIFIED

**Verified 2026-07-24** by independent recomputation from the long-read DGEList:
`code/upstream/verify_section3_p1.R` (reads no cached result). Every published number
reproduces exactly.

| # | Claim | Source code | Status |
|---|---|---|---|
| 3.1 | Gene-level % output lost: 11.1% (AT2) — 12.5% (MV) **(SF20)** | `yul:transcriptional_output.Rmd` chunks `gene-level-output` + `gene-level-global-metric`; LR isoform CPM aggregated to gene via `rowsum` | ✅ **EXACT** — recomputed AT2 11.06, FB 11.83, LAE 12.03, MV 12.50; argmin/argmax CTs match. |
| 3.2 | Isoform-level % output lost: 14.4% (AT2) — 18.2% (LAE) **(Fig 2A)** | Same Rmd, chunk `global-metric` | ✅ **EXACT** — recomputed AT2 14.40, FB 14.56, MV 15.53, LAE 18.23; argmin/argmax match. |
| 3.3 | Per-isoform % lost: median 60–67%, IQR 30–100% **(SF21)** | Same Rmd, chunk `boxplot-isoform` (filtered isoform set, `frac_lost > 0`) | ⚠️ **median EXACT** (60.2 MV – 67.1 LAE); **quartile range imprecise** — see F-1. |
| 3.4 | **Figure 2A** — isoform-level % output lost across CTs | `yul:transcriptional_output.Rmd` chunk `barplot-global` (writes `pct_output_lost_barplot.png`) | Producer identified (previously "TBD"). |

**Definition as implemented** (per cell type; paired within donor, matched pairs only —
AT2 3, FB 3, LAE 4, MV 3):

```
% output lost = 100 × Σ_i Σ_d max(CPM_smg1i − CPM_dmso, 0) / Σ_i Σ_d CPM_smg1i
```

Normalization is plain CPM with `normalized.lib.sizes = FALSE` (TMM deliberately avoided —
it would scrub the compositional shift NMD inhibition creates).

### §3 ¶1 findings

**F-1 (documentation, minor).** SF21's stated "IQR 30–100%" does not hold for MV: observed
quartiles are Q1 25.7 / Q3 96.2 (AT2 31.4/100, FB 30.2/100, LAE 31.6/100). Either widen the
stated range to ~26–100% or scope it to AT2/FB/LAE. Note also that the manuscript uses "IQR"
to mean the Q1–Q3 *interval*, not the IQR width (which is 68.4–70.5).

**F-2 (map correction).** The previous map entry described this metric as computed "over
NMD-susceptible features." **It is not** — `global-metric` sums over **all** isoforms in
`results_all` (645,272 isoforms × donor pairs) with no susceptibility filter.

**F-3 (result correctness — RAISE WITH YUL BEFORE SUBMISSION).** The statistic is provably
the **total variation distance** between the DMSO and Smg1i isoform-composition vectors, not
a directional measure of degraded output. Diagnostics:
`code/upstream/verify_section3_p1_metric_diagnostics.R`.

- **Closed composition.** Every CPM column sums to exactly 1e6, so `sum(Δ) = 0` and therefore
  `Σmax(Δ,0) ≡ Σmax(−Δ,0) ≡ ½Σ|Δ|`.
- **D2 symmetry.** Recomputed in the reverse direction (DMSO − Smg1i) the value is
  **identical** in all four CTs (14.40 / 14.56 / 18.23 / 15.53). The statistic cannot by
  itself distinguish output "lost" from output "gained" — the totals are equal by
  construction (the contributing isoform *sets* do differ).
- **D1 attribution.** Only **20.7–54.5%** of the signal comes from isoforms that are
  statistically NMD susceptible (mashr lfsr<0.05 & posterior mean>0): AT2 4.58 of 14.40,
  FB 3.01 of 14.56, LAE 9.93 of 18.23, MV 6.03 of 15.53 percentage points. The remainder is
  sub-threshold NMD and/or noise — `max(Δ,0)` has positive expectation even when true Δ=0.
- **D3 no-treatment baseline.** The same statistic between two **DMSO** samples from
  different donors gives 13.4–33.3%, exceeding the treatment value in 3 of 4 CTs.
  **Caveat:** cross-donor contrasts carry genuine donor biology, so this is an upper bound on
  a no-treatment baseline, not a matched technical null.

**Why it matters.** This metric backs an **Abstract** claim ("NMD degraded 11–18% of
transcriptional output"), **Fig 2A**, **SF20/SF21**, and the **Discussion** comparison to
Fair et al. (15%). The numbers reproduce exactly (step 1 ✅); what is unresolved is whether
the statistic supports the *interpretation* (step 2). Options: restrict the numerator to
NMD-susceptible isoforms (yields 3.0–9.9%), report it explicitly as a compositional-shift /
TVD measure, or add a permutation baseline showing the treatment effect exceeds donor-level
variability.

### Paragraph 2 — custom NMD-anchored normalization + productive DE (Fig 2B, ST4) — ✅ VERIFIED (3.5–3.10)

**Verified 2026-07-24** by `code/upstream/verify_section3_p2.R`, which purls
`productive_response.Rmd`, repoints its `/udd/reyle` paths at local copies of the identical
inputs, and executes **Yul's own code** through the `mashr` chunk. Every quoted number
reproduces. Two **disclosed deviations**, both lookup-only (no data/threshold/computation
changed): the repo's mashr CSVs predate the April-2026 rename, so `resolve_mash_col()` and
`resolve_sr()` were widened to accept `Smg1i_in_AT`/`_in_DD` as aliases for AT2/LAE.

| # | Claim | Source code | Status |
|---|---|---|---|
| 3.5 | Correction preserves gene rankings: Spearman ≥ 0.99 (adjusted vs uncorrected logFC) | `productive_response.Rmd`, `cor(custom, Kitagawa)` | ✅ **VERIFIED** — 0.9980 (LAE) – 0.9997 (FB/AT2); NMD-only 0.9984–0.9997 |
| 3.6 | Non-NMD genes closely matched DMSO vs adjusted-Smg1i: Spearman ≥ 0.96 **(SF23)** | Same, chunk `R5-control-corr` | ✅ **VERIFIED** — 0.9643 (LAE) – 0.9705 (FB); min ≥ 0.96 holds |
| 3.7 | FB treated as provisional due to **donor confounding** | Methods + same Rmd | ⚠️ **MISMATCH — see F-7.** `productive_response.Rmd:61` gives a different reason: "smaller effect sizes; fails part of the reproducibility battery". |
| 3.8 | Median adjusted log2FC near zero: AT2 +0.14, LAE −0.04, FB +0.10, MV +0.05 | Same, chunk `base-build` | ✅ **EXACT** — AT2 0.141, LAE −0.040, FB 0.098, MV 0.051 |
| 3.9 | No significant change in ~83% of genes **(Fig 2B; ST4)** | Same, chunk `mashr` | ⚠️ **numbers reproduce, framing misleading — see F-6.** Per-CT: AT2 82.9, FB 93.8, LAE **59.8**, MV 83.9. "83%" = median across CTs (83.4), not pooled (79.7 all / 75.2 excl FB). |
| 3.10 | Downregulated 88 (AT2) – 1337 (LAE); upregulated 635 (MV) – 888 (LAE) **(ST4)** | Same | ⚠️ **all four numbers EXACT, range framing wrong — see F-5.** FB (down 11, up 306) falls below both stated minima. |
| 3.11 | Up-genes enriched for UPR (leading edge *HSPA5*/BiP, *XBP1*, *EIF2AK3*) and NF-κB/TNFα **(ST4)** | `leading-edge-check` (GO:0006986) + `step16-2-up-NFkB` (`HALLMARK_TNFA_SIGNALING_VIA_NFKB`) | ⚠️ **PARTIAL** (`verify_section3_p2_claim311.R`). **NF-κB/TNFα ✅** — productive pm shifted UP in NMD genes in all 4 CTs (p 4.3e-05…0.008). **UPR ⚠️** — all three named genes appear, but no single CT has all three (AT2: HSPA5+EIF2AK3; LAE/MV: HSPA5+XBP1; FB: HSPA5 only), and a direct Fisher test of UPR membership among up-genes is significant in only 2 of 4 CTs (LAE p=0.004, MV p=0.034; FB p=0.053, **AT2 p=0.44**). See F-12. |
| 3.12 | **Figure 2B** — productive mean logFC distribution, NMD vs non-NMD genes | ✅ **TRACED** — `productive_response.Rmd:505` (chunk `dist-others`): `plot_dist(base, "custom", "Adjusted productive log2FC", xlim=c(-2,2), vline=0)`; `plot_dist()` facets by cell type and fills by `gene_category`, matching the legend | Producer identified |

**Population note.** The manuscript's DE population is genes with ≥1 NMD-susceptible isoform
(`gene_category == "NMD"`): AT2 5,104 · FB 5,101 · LAE 5,535 · MV 5,301 after the expression
floor. An `"unsure"` category (borderline `0.05 ≤ lfsr < 0.20`) is held to the supplement.

### Paragraph 3 — Kitagawa decomposition (ST4) — ✅ VERIFIED (3.13–3.15)

**Verified 2026-07-24** by `code/upstream/verify_section3_p3.R`, calling Yul's own
`kit_nmd_table()` / `kitagawa_decomp()` via the shared harness
`code/upstream/_run_productive_response.R`. **All three percentages reproduce exactly.**

Population (per the code, not restated in the manuscript prose): significant NMD genes
(`dir ∈ {up,down}`, `gene_category=="NMD"`), **FB excluded as provisional** →
AT2 871 · LAE 2,225 · MV 853 = **3,949 genes**. Decomposition identity holds
(max |dP − (level+composition)| = 9.1e−13).

| # | Claim | Recomputed | Status |
|---|---|---|---|
| 3.13 | Level dominant 86.2%, composition 13.5% | level 3,403 = **86.2%**; composition 534 = **13.5%** | ✅ **EXACT** |
| 3.14 | Productive-up: 99.2% level, 0.3% composition | level 2,288 = **99.2%**; composition 7 = **0.3%** | ✅ **EXACT** |
| 3.15 | Productive-down: 32.1% composition, 67.9% transcription | composition 527 = **32.1%**; level 1,115 = **67.9%** | ✅ **EXACT** |
| 3.16 | Conclusion holds under the custom normalization with SR induction as an external signal | — | Not yet run (S6b `down-arm-mechanism`, past the harness stop point) |

**Residual resolved.** The sums (99.7% / 99.5% / 99.9%) are **not** an interaction term and
**not** ties (exact `|composition| == |level|` ties = **0 of 3,949**). They are **NA rows** —
12 overall, 11 in the up-arm, 1 in the down-arm — genes present in `prod_dir_mashr` but
dropped from `base` by the expression floor, so the `left_join` in `kit_nmd_table()` yields
no CPM. NAs sit in the denominator, so the published percentages very slightly understate the
split among genes with a defined driver (86.2% → 86.4% if NA rows are dropped). Cosmetic, but
it explains the arithmetic.

**Sensitivity.** Including FB shifts the headline to 87.1% level / 12.5% composition.
⚠️ The FB exclusion is **not restated in the ¶3 prose** — see F-8.

### Paragraph 4 — two mechanisms partitioning productive-down genes (ST4) — ⚠️ MOSTLY VERIFIED

**Verified 2026-07-24** by `code/upstream/verify_section3_p4.R`, replicating Yul's S6b chunk
`down-arm-mechanism` with her own `classify_mech()` / `theme_enrichment()` / `THEME_SETS`.
Population: `dir=="down"`, FB excluded → transcriptional 1,631 · composition 430 ·
universe 9,196.

| # | Claim | Recomputed | Status |
|---|---|---|---|
| 3.16 | Conclusion holds under the custom normalization + SR induction | 3-way split: transcriptional **70.1%**, composition 16.5%, mixed 13.4% — agrees with Kitagawa's 67.9% level-dominated | ✅ **SUPPORTED** |
| 3.17 | Composition genes × RNA splicing, OR 2.5, p<10⁻³ | OR **2.52**, p = 9.68e-07 | ✅ **EXACT** |
| 3.18 | Composition genes × SR/hnRNP factors, OR 6.8, p<10⁻³ | OR **6.76**, p = 2.61e-04 | ✅ **EXACT** |
| 3.19 | Of 12 canonical SR proteins, every significant productive change was a **decrease** (8 LAE, 3 AT2, 2 MV) | LAE **8**, AT2 **3**, MV **2**, all `down`; **zero** significant increases | ✅ **EXACT** — but only **11 of 12** SRSF genes exist in the data (SRSF12 absent); see F-10 |
| 3.20 | Transcriptional genes **instead** enriched for DNA replication and repair, OR 2.1, p<10⁻³ | DNA **replication** trans OR **2.12** ✅; DNA **repair** trans OR **1.57** ✗ — and neither is specific to the transcriptional class | ⚠️ **PARTIALLY CONTRADICTED — see F-9** |

**Full theme-enrichment table (Yul's `theme_enrichment()` output):**

| theme | trans_OR | trans_p | comp_OR | comp_p |
|---|---|---|---|---|
| Mitotic cell cycle | 0.98 | 5.7e-01 | 1.28 | 9.9e-02 |
| DNA replication | **2.12** | 9.3e-07 | **2.03** | 4.5e-03 |
| DNA repair | **1.57** | 3.4e-05 | **1.88** | 4.2e-04 |
| RNA splicing | 1.19 | 1.0e-01 | **2.52** | 9.7e-07 |
| SR/hnRNP factors | 0.76 | 7.6e-01 | **6.76** | 2.6e-04 |
| ISR/UPR | 0.70 | 9.3e-01 | 0.51 | 9.3e-01 |

The **composition-side** contrast is clean and strongly specific (RNA splicing 1.19→2.52;
SR/hnRNP 0.76→6.76). The **transcriptional-side** contrast is not (see F-9).

### Paragraph 5 — three worked examples (Fig 2C, 2D, 2E) — ⚠️ 2 of 3 VERIFIED

**Verified 2026-07-24** by `code/upstream/verify_section3_p5.R` (per-gene lookups from Yul's
`base` / `sr_long` / `prod_dir_mashr`). All four cell types were pulled for each gene so the
values themselves identify the intended CT.

| # | Claim | Recomputed | Status |
|---|---|---|---|
| 3.21 | *SHMT2* (**MV**): NMD susceptible, **significantly** increased; productive 86.5 → 111.1 CPM, adjusted logFC +0.42, SR logFC +0.56 **(Fig 2C)** | Every value matches **AT2**, not MV (86.5 / 111.1 / 0.421 / 0.561). MV is 44.1 / 54.8 / 0.391 / 1.161. And SHMT2 is **`ns` in all four CTs** | 🔴 **CONTRADICTED — see F-11** |
| 3.22 | *SRSF2*: productive 84 → 50 CPM, adjusted logFC −0.73, fraction 91% → 26% **(Fig 2D)** | **MV**: 84.1 → 50.4, custom **−0.728**, fraction **91.0% → 26.4%** | ✅ **EXACT** — and this **resolves F-4**: the numbers confirm MV, so the Fig 2D legend is right; only the prose omits the CT |
| 3.23 | *PCNA* (**LAE**): productive 251 → 130 CPM, adjusted logFC −0.95, fraction 100% → 93%, SR logFC −0.81 **(Fig 2E)** | **LAE**: 251.5 → 129.9, custom **−0.947**, fraction **99.5% → 92.5%**, sr_logFC **−0.811** | ✅ **EXACT** |
| 3.24 | Fig 2C/2D/2E rendered using **Isopair** | ✅ **TRACED** — `Isopair::plotIsoformPair()`, defined `Isopair/R/visualization.R:47`, exported at `NAMESPACE:41` | Confirmed |

**Quoting-basis note.** SRSF2 and PCNA both quote the **adjusted** SMG1i productive CPM
(50.4→"50"; 129.9→"130"). SHMT2's "111.1" is the **raw** value (`prodCPM_SMG_raw`); its
adjusted value is 116.2 — yet the companion "+0.42" logFC is computed from the adjusted
figure. So the three examples are not quoted on a consistent basis (part of F-11).

### Section 3 status summary (2026-07-24)

| | |
|---|---|
| ¶1 (3.1–3.4) | ✅ verified exact; findings F-1 (minor), F-2 (map), F-3 (substantive) |
| ¶2 (3.5–3.10) | ✅ verified — every number reproduces via Yul's own code; findings F-5, F-6, F-7 |
| ¶2 (3.11–3.12) | ✅ 3.12 traced to `dist-others`; ⚠️ 3.11 partial — NF-κB verified, UPR enrichment holds in 2 of 4 CTs (F-12) |
| ¶3 (3.13–3.15) | ✅ verified exact via Yul's `kit_nmd_table()`; residual = NA joins, not ties; finding F-8 |
| ¶4 (3.16–3.20) | ✅ 3.16/3.17/3.18/3.19 verified exact; ⚠️ **3.20 partially contradicted** (F-9); F-10 minor |
| ¶5 (3.21–3.24) | ✅ 3.22/3.23 exact (SRSF2=MV, PCNA=LAE); 🔴 **3.21 SHMT2 contradicted** (F-11); 3.24 traced to `Isopair::plotIsoformPair()` |
| Blocking access | **none** — code and data both local |
| Open for Yul | 🔴 **F-3** metric interpretation · 🔴 **F-9** DNA repl/repair "instead" · 🔴 **F-11** SHMT2 wrong CT + "significantly" · **F-6** the "83%" · **F-8** LAE drives the 13.5% · **F-5** FB in ranges · **F-7** FB provisional reason · F-1 SF21 quartiles · F-4 SRSF2 cell type · F-10 SR-protein denominator |
| Verification scripts | `verify_section3_p1.R`, `verify_section3_p1_metric_diagnostics.R`, `verify_section3_p2.R`, `verify_section3_p3.R`, `verify_section3_p4.R`, `verify_section3_p5.R`, harness `_run_productive_response.R` (all `code/upstream/`) |
| Running issues list | `paper/VERIFICATION_ISSUES.md` |
---

## Section 4 (Isopair: attributing NMD to splicing events) → code

**Scope framework (current canonical, 2026-07-11; post 25% reference-share floor):** §4 manuscript prose uses three nested analytic scopes, each with explicit gene-matched NMD vs Control pairs. All n's below are POST-floor (pop_BC 1,548; see M11 for the floor definition).

| Scope | Definition | NMD n | Control n | Used for |
|---|---|---|---|---|
| **Subset 1 / `gencode_all3` (strict)** | all 3 transcripts in the triplet are GENCODE-annotated coding (reference + NMD comparator + Control comparator) | 130 | 130 | Headline PTC enrichment (24-fold, 36.9% vs 1.5%); 5'UTR length / longest 5'UTR ORF (Fig 4 A/B); CDS/3'UTR sanity (CDSand3UTR supplement row 1) |
| **Subset 2 / `mechanism_class_4` (broad)** | ref-AUG-traceable; `mechanism_class ∈ {effectively_ptc, truncated_no_ejc, no_downstream_ejc, ref_atg_lost}` | 819 | 819 | Mechanism breakdown; 5'UTR analysis at expanded scope (Fig 4 C/D); cumulative PTC accounting (756/819 = 92%); 3'UTR length comparison (1324/722/922) |
| **Occult-PTC** | within Subset 2 NMD+ where ref AUG defines a PTC but TD2 misses it | 348 | n/a | TD2 bias evidence (PTC rate, ref-AUG vs TD2 Kozak strength, TD2 ATG position) — TD2BiasEvidence supplement |
| **Upstream cell-type cohort** | AT (n=6) + DD (n=8) + FB (n=6) + MV (n=6) = 26 samples; 1,548 genes pass pop_BC filtering | — | — | All upstream pair construction |

**Primary code home for §4:** `nmd:results/isoform_transitions/Version_6.0/isopair_wrapper/` (force-tracked in this repo).

- **Canonical report Rmd:** `05_final_report_gencode_scope_2026-07-11.Rmd` (NEW 2026-06-15, branch-don't-rewrite from the legacy). Sections §1 (pop_BC + Fig 3 A/B/C) → §2 (n=130 strict) → §3 (n=819 broad) → §4 (TD2 bias / n=348 occult-PTC) → §5 (cumulative accounting). Each manuscript-§4 paragraph maps to specific Rmd chunks (table below).
- **Wrapper scripts:** `01_prepare_data_mashr.R`, `02_build_profiles_mashr.R`, `03b_rebuild_cache.R`, plus sub-scripts `05k_utr5_all_isoforms.R`, `05k_b_utr5_refaug.R`, `05r_ref_atg_analysis.R`. (Legacy-report-only feeders `04`/`04b`/`05l`/`05s`/`05m`/`05t`/`05u`/`05v` removed 2026-07-11.)
- **Legacy Rmd (DEPRECATED):** `05_final_report_mashr.Rmd` retained with deprecation banners. **Do not cite for any new analysis or claim**; numbers in it are no longer canonical.
- **Backing package:** `isopair:` (`peter4244/Isopair`).

**Supplements (NEW 2026-06-15, in `figures/SupplementalFigures/`):**

| Supplement | Scope | What it shows |
|---|---|---|
| `CDSand3UTR_GENCODEonly/` | n=130 (row 1) + n=819 (row 2) | 2×3 layout: CDS length, 3'UTR translation-based, 3'UTR non-PTC-stop bias-corrected — replicates the bias-correction pattern across both scopes |
| `TD2BiasEvidence/` | n=819 broad + n=348 occult-PTC | 2×3: TD2-vs-ref-AUG ORF length, paired Kozak PWM, TD2 ATG position |
| `PairAnalysisFlowchart/` | full cohort | DiagrammeR/grViz cohort cascade: 4-CT cohort (26 samples) → NMD classification → Subset 1 (strict, n=130) + Subset 2 (broad, n=819) |

**Manuscript §4 paragraph → Rmd chunk → Figure 3/4/supplement panel map:**

| Manuscript §4 paragraph | Topic | Rmd section / chunk | Figure / panel |
|---|---|---|---|
| ¶1 (pair construction, descriptives) | pop_BC = 1,548 genes; transcript length 2896/3018/2742 nt; 7 isoforms/gene | §1, §1a (`sec1-pop-bc`, `sec1-tx-length`, `sec1-annotation`) | PairAnalysisFlowchart supplement |
| ¶2 (sequence similarity + splice events) | NMD shares more with reference; SE 2× more in NMD; IR more common in Control | §1c (`sec1-fig3abc`) | Fig 3 B, C |
| ¶3 (PTC enrichment + mechanism) | 24-fold (36.9% vs 1.5%); SE 48%/14%, A5SS 13%/4%; frameshift/in-frame/3'UTR = 52%/38%/10% | §2a (`sec2a-ptc`, `sec2a-fig3-panelD`), §2b (`sec2b-attribution`, `sec2b-fig3-panelE-F`) | Fig 3 D, E, F |
| ¶4 (5'UTR + uORF at n=130 and n=819; TD2 vs ref AUG) | 5'UTR longer for NMD+/PTC−; 52% TD2 == ref AUG; 348 occult-PTC; 99% TD2 downstream; 82% ref AUG Kozak stronger | §2c (`sec2c-fig4AB`), §3a (`sec3a-fig4CD`), §4 (`sec4-compute`, `sec4-supp`) | Fig 4 A, B, C, D; TD2BiasEvidence supplement |
| ¶5 (3'UTR + cumulative PTC accounting) | No 3'UTR diff at n=130; expanded set median 1324/722/922; 90% PTC at n=819 | §2d (`sec2d-utr3np`, `sec2d-supp`), §3b, §5 (`sec5-accounting`) | CDSand3UTR_GENCODEonly supplement; cumulative accounting table |

**Verifiers:** `nmd:reproducibility/verify_pass7_new_rmd.R` (37-check manifest) + `nmd:reproducibility/verify_cross_check_new_rmd_vs_figures.R` (57-check cross-check). Both PASS as of 2026-06-15.

**Key package-to-script bridge** (from `05r_ref_atg_analysis.R` header):

```
05r_ref_atg_analysis.R  →  Isopair::traceReferenceAtg()  →  categorical NMD+/PTC- subgroup labels
                            (effectively_ptc, truncated_no_ejc, no_downstream_ejc,
                             ref_atg_lost, no_ref_cds, mapping_failed)
```

**Per-claim entries (4.1–4.46) below have been updated to the current canonical Rmd + scope framing.**

### Paragraph 1 — pair-set construction (1,548 genes; descriptives)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 4.1 | "12 categories of splicing events" enumerated by Isopair | `isopair:R/event-detection.R` — defines the 12 categories. Per Methods M11 + Fig 3F legend: Alt TSS, Alt TES, exon skipping, intron retention, A5SS, A3SS, MXE, alt first/last exon, terminal exon extension/truncation, combination events. Caption defines abbreviations: SE, IR, A3SS, A5SS, Partial IR 5'/3', IR diff 5'/3', Alt_TSS, Missing_Int. | **Verifiable locally.** The 12-category list in the Rmd output should match the enumeration in `event-detection.R`. |
| 4.2 | Correctness verified by ability to reconstruct reference exon structure | `isopair:R/reconstruction.R` — reconstruction function. Per ONBOARDING the validator is built into the package. | **Verifiable locally.** |
| 4.3 | "1,548 genes from which the isoform pairs were drawn" (≥3 expressed isoforms, ≥1 NMD-susceptible) | `nmd:results/.../isopair_wrapper/01_prepare_data_mashr.R` + `02_build_profiles_mashr.R` — pair-set filtering. Per Methods: "≥5% of overall gene expression in either DMSO or SMG1i and ≥5 reads in ≥1 sample"; non-NMD definition was lowered from adj.P.Val > 0.50 to > 0.30 per ONBOARDING §6. Final pop_BC = 1,548 genes. | **Verifiable locally.** New Rmd `sec1-pop-bc` chunk. Pinned by verifier `verify_pass7_new_rmd.R` (1,548). |
| 4.4 | Median 7 isoforms per gene (SF: Isoform Count in Isoform Pair Sets) | New Rmd §1 — pop_BC summary statistic. | Verifiable locally. **§1b descriptive not currently in the pass-7 manifest** — Pete deferred recomputation 2026-06-16; numbers carried over from the legacy Rmd. |
| 4.5 | Reference isoform share (FLOORED): median **67%** of parent gene expression; **71%** of references >50% (all retained refs ≥25% by the floor) | New Rmd §1 `sec1-ref-share` (complement-of-NMD denom, all_samples DMSO). | **DRIFT RESOLVED.** The legacy "70%/75%" measured the gene DOMINANT isoform's share (`dominant_pct = max/total`), mislabeled as the reference share; post-floor the reference *is* the dominant for retained genes, so the honest reference-share stat (67%) matches SF28 (all-iso 66.7%). Renders in the Rmd; not in the pass7 37-check manifest. |
| 4.6 | NMD + non-NMD comparators expressed at similar levels in SMG1i samples (SF: Reference Isoform Share of Gene Expression) | New Rmd §1. | Verifiable locally. |
| 4.7 | Transcript length: median 2,958 (reference) / 2,991 (NMD comparator) / 2,808 (Control comparator), paired Wilcoxon p<0.001 (SF: Transcript Length Comparison) | New Rmd `sec1-tx-length` chunk — paired Wilcoxon test on bias-corrected pop_BC-restricted c4 transcript lengths. | **Verifiable locally.** Pinned by `verify_pass7_new_rmd.R`. *(Legacy 3067/3107/2906 numbers were pre-recomputation; Pete corrected the manuscript prose on 2026-06-16.)* |
| 4.8 | Flowchart for the 1,548-gene pair construction (Supplemental Figure: PairAnalysisFlowchart) | NEW supplement: `figures/SupplementalFigures/SF26_PairAnalysisFlowchart/build_flowchart.R` (DiagrammeR/grViz). Renders to `.dot` + `.html`. Validated by `figures/lib/validate_flowchart_dot.R`. New Rmd embeds the flowchart in §1 (`sec1-flowchart` chunk). | **Verifiable locally.** |
| 4.9 | **Figure 3A** (schematic illustrating reference / NMD-comparator / Control-comparator pair construction) | `figures/multipanel/figure3_isopair_and_ptc/figure3_panelA_pair_concept.py` (matplotlib port of the original `make_pair_concept_figure.R`). | **Verifiable locally.** Embedded by new Rmd `sec1-fig3abc`. |

### Paragraph 2 — Figure 3 panels B, C (sequence similarity, event prevalence)

> Note: previous "Figure 3D = gain/loss direction" has been retired. The current Figure 3 has six panels A–F, with D = stop-to-EJC distance (Paragraph 3 below), E = PTC-causing event attribution, F = mechanism breakdown.

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 4.10 | "NMD isoforms shared more sequence content with reference than Control isoforms" — sequence similarity comparison | `isopair:R/compare-sets.R` — sequence-similarity computation. Driver in new Rmd `sec1-fig3abc` chunk, panel `figure3_panelB_sequence_similarity.py`. | **Verifiable locally.** |
| 4.11 | **Figure 3B** (sequence similarity of NMD vs Control pairs) | `figures/multipanel/figure3_isopair_and_ptc/figure3_panelB_sequence_similarity.py` + `data_export.R`. New Rmd embeds via `sec1-fig3abc`. Scope: pop_BC (n=1,548). | **Verifiable locally.** |
| 4.12 | "Skipped exon events twice as frequent in NMD pairs vs Controls" — splice-event prevalence | `isopair:R/event-detection.R` + per-category event counts via new Rmd §1c. | **Verifiable locally.** |
| 4.13 | "Intron retention more common in Control pairs" | Same. | Verifiable locally. |
| 4.14 | **Figure 3C** (prevalence of splice event categories in NMD vs Control) | `figures/multipanel/figure3_isopair_and_ptc/figure3_panelC_event_prevalence.py` + `data_export.R`. New Rmd embeds via `sec1-fig3abc`. Scope: pop_BC (n=1,548). | **Verifiable locally.** |
| 4.15 | "Terminal events (Alt TSS / TES) and skipped exons more frequently led to sequence GAIN in NMD; intron retention led to sequence gain in Controls" *(now in supplement per Pete's plan)* | `isopair:R/event-detection.R` (GAIN/LOSS semantics — clarified in commit `d242798`). Moved out of main Figure 3 in the 2026-06-15 refactor. | Verifiable locally; supplementary. |
| 4.16 | ~~Figure 3D (gain/loss direction)~~ — **retired from main Figure 3.** Current Figure 3 Panel D is the stop-to-EJC-distance panel (was old Figure 4A). See claim 4.19. | n/a | superseded |

### Paragraph 3 — PTC enrichment + mechanism breakdown (Figure 3 D, E, F at n=130)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 4.17 | "24-fold enrichment of PTCs in NMD-susceptible isoforms (36.9% vs 1.5% PTC rate, p = 1.58×10⁻¹⁴)" at n=130 (Subset 1, gencode_all3) | `isopair:R/ptc-attribution.R` + `ptc.R` (50-nt rule). New Rmd `sec2a-ptc` chunk computes per-pair PTC determination at the strict GENCODE-restricted scope. | **Verifiable locally.** Pinned by `verify_pass7_new_rmd.R` (36.9% / 1.5% / OR=37.06 / p=1.58e-14). *(Legacy "44% vs 5%, 15-fold" was pre-bias-correction.)* |
| 4.18 | "Dose-response relationship between upstream distance of the stop codon and the magnitude of the NMD response" | `isopair:R/spatial.R` — stop-codon-to-last-EJC distance computation. New Rmd `sec2a-ptc`. | Verifiable locally. |
| 4.19 | **Figure 3D** (stop-codon-to-last-EJC distance, NMD vs Control, n=130) | `figures/multipanel/figure3_isopair_and_ptc/figure3_panelD_stop_codon_distance.py` + `panel_e_compute.R`. New Rmd embeds via `sec2a-fig3-panelD`. | **Verifiable locally.** |
| 4.20 | "NMD response peaked at 4–5 downstream EJCs from PTC" | New Rmd §2a — joins mashr posterior logFC with downstream EJC count. | Verifiable locally. |
| 4.21 | "Frameshift events were the leading mechanism (55%)" — 38/69 PTC+ pairs | New Rmd `sec2b-attribution` + `sec2b-mech-table`. Frameshift n=38. | **Verifiable locally.** Pinned by verifier (55.1%). |
| 4.22 | "In-frame stop events (33%)" — 23/69 | Same. | Pinned by verifier (33.3%). |
| 4.23 | "3'UTR splicing (12%)" — 8/69 | Same. | Pinned by verifier (11.6%). |
| 4.24 | "Skipped exon was the most common event type (47% of PTC-causing events vs 15% in Controls, Fisher p < 10⁻⁶)" | New Rmd §2b — per-event Fisher table at n=130 PTC+. Panel E TSV: SE 47.9% / 14.4% / p=6.88e-7. | **Verifiable locally.** Pinned by verifier. |
| 4.25 | "A5SS also significantly enriched (13% vs 4.5%, p < 10⁻²)" | Same — Panel E TSV: A5SS 13% / 4.5% / p=0.00888. | **Verifiable locally.** Pinned by verifier. |
| 4.26 | **Figure 3E** (PTC-causing event attribution, per-event Fisher) + **Figure 3F** (mechanism breakdown × event type) | `figures/multipanel/figure3_isopair_and_ptc/figure3_panelE_ptc_event_attribution.py` + `figure3_panelF_mechanism_breakdown.py` + `panel_e_compute.R`. New Rmd embeds via `sec2b-fig3-panelE-F`. | **Verifiable locally.** |

### Paragraph 4 — 5'UTR / uORF analysis + occult-PTC + TD2 bias (Figure 4 + TD2BiasEvidence supplement)

> Note: The old "Figure 5" structure (TD2 vs ref-AUG occult-PTC panels) has been retired. Occult-PTC + TD2 bias evidence are now in the **TD2BiasEvidence supplement**; 5'UTR / longest uORF analysis at n=130 and n=819 is **Figure 4**.

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 4.27 | "1,116 NMD and Control pairs" *(typo in manuscript prose — should read 1,166; Pete corrected 2026-06-16)*; expanded ref-AUG-traceable pair set | New Rmd `sec3-scope` — `mechanism_class_4` filter yields n=819 NMD+Control matched pairs (was 1,166 pre-floor; manuscript prose to find/replace in Phase 7). | **Verifiable locally.** Pinned by verifier (819). |
| 4.28 | "50% of these cases, the TD2-called CDS used the same AUG as the reference isoform CDS" | New Rmd §4 (`sec4-compute`, TD2 vs ref-AUG ORF length). At Subset 2 (n=819): TD2 == ref AUG in 429/819 = 52.4%. | **Verifiable locally.** Pinned by verifier. |
| 4.29 | "492 novel NMD+ isoforms where the ORF defined by the reference AUG included a PTC" — occult-PTC scope | New Rmd §4 — within Subset 2 NMD+, count of pairs where ref-AUG defines a PTC but TD2 misses it = 348 (manuscript prose 492 to find/replace). Mechanism via `05r_ref_atg_analysis.R::traceReferenceAtg()` (`effectively_ptc` category). | **Verifiable locally.** Pinned by verifier. |
| 4.30 | "99.1% (345/348) of the TD2-called CDS were downstream from the reference AUG" — Figure S TD2BiasEvidence panel | New Rmd `sec4c` (TD2 ATG position) + `figures/SupplementalFigures/SF34-SF35_TD2BiasEvidence/data_export.R`. | **Verifiable locally.** Pinned by verifier (99.1%). |
| 4.31 | "Reference AUG was stronger in 82% of cases (286/348)" — paired Kozak PWM, occult-PTC scope | New Rmd `sec4b` (paired Kozak) + TD2BiasEvidence supplement. Wilcoxon paired test p = 9.34e-36. | **Verifiable locally.** Pinned by verifier (82.2%). |
| 4.32 | **Figure 4 Panel A** (5'UTR length, NMD+/PTC+ vs NMD+/PTC− vs Control at n=130) | `figures/multipanel/figure4_ptcneg_and_model/figure4_panelA_5utr_length_all3enst.py` + `data_export.R`. New Rmd embeds via `sec2c-fig4AB`. NMD+/PTC− 5'UTR 426.5 nt vs PTC+ 87.5 vs Control 145. | **Verifiable locally.** Pinned by verifier. |
| 4.33 | **Figure 4 Panel B** (longest 5'UTR ORF, n=130) | `figures/multipanel/figure4_ptcneg_and_model/figure4_panelB_longest_5utr_orf_all3enst.py`. New Rmd embeds via `sec2c-fig4AB`. | **Verifiable locally.** |
| 4.34 | **Figure 4 Panel C** (5'UTR length, ref-AUG-projected, n=819) | `figures/multipanel/figure4_ptcneg_and_model/figure4_panelC_5utr_length_refaug.py`. Uses `05k_b_utr5_refaug.R` (NEW ref-AUG-projected scan, TD2-bias-free). New Rmd embeds via `sec3a-fig4CD`. | **Verifiable locally.** |
| 4.35 | **Figure 4 Panel D** (longest 5'UTR ORF, ref-AUG-projected, n=819) | `figures/multipanel/figure4_ptcneg_and_model/figure4_panelD_longest_5utr_orf_refaug.py`. New Rmd embeds via `sec3a-fig4CD`. | **Verifiable locally.** |
| 4.36 | TD2BiasEvidence supplement (2×3 layout: ORF length, paired Kozak, TD2 ATG position; rows = n=819 broad and n=348 occult-PTC) | `figures/SupplementalFigures/SF34-SF35_TD2BiasEvidence/figure_s_td2_bias.py` + `data_export.R`. New Rmd embeds via `sec4-supp`. | **Verifiable locally.** |
| 4.37 | "PTCs in 92% of the NMD+ isoforms (756/819)" — cumulative PTC accounting at expanded scope | New Rmd `sec5-accounting` + `sec5-acct-table`. Composite of original TD2-detected PTC + occult-PTC under ref-AUG tracing. | **Verifiable locally.** Pinned by verifier (92.1% = 756/819). *(Legacy "85%" was at the prior 2,289-pair scope.)* |

### Paragraph 5 — 3'UTR length comparisons (n=130 and n=819)

> Note: the previous "498 isoforms remaining / three subgroups" framing was at the old 2,289-pair scope before the 2026-06-15 refactor reframed to Subset 1 / Subset 2 / Occult-PTC. The 3'UTR + 5'UTR / uORF analyses now sit in Figure 4 (claims 4.32–4.35 above) + the CDSand3UTR_GENCODEonly supplement; the per-subgroup mechanism breakdown sits in Figure 3 Panels D/E/F (claims 4.19–4.26).

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 4.38 | "NMD+/PTC− isoforms had longer 5'UTRs and uORFs than both other groups" at n=130 | New Rmd §2c (`sec2c-secA-groups`, `sec2c-panelA-tables`, `sec2c-panelB-tables`). 5'UTR: PTC+ 87.5 / PTC− 426.5 / Ctrl 145 nt. | **Verifiable locally.** Pinned by verifier. |
| 4.39 | "No significant 3'UTR length difference at n=130 (GENCODE-restricted set)" — bias-corrected 3'UTR | New Rmd §2d (`sec2d-utr3np`). Bias-corrected 3'UTR (non-PTC-stop): PTC+ 550 / PTC− 353.5 / Ctrl 581.5, p=0.15. | **Verifiable locally.** Pinned by verifier. |
| 4.40 | "Significantly longer 3'UTR in NMD+/PTC+ at expanded set (median 1311 / 834.5 / 948 nt)" | New Rmd §3a (`sec3a-secC-groups`) + §2d (`sec2d-utr3np` at n=819 row). Bias-corrected 3'UTR at Section C scope. | **Verifiable locally.** Pinned by verifier (1311 / 834.5 / 948 exact match). |
| 4.41 | "Mechanism classification by subgroup at expanded scope" (effectively_ptc, truncated_no_ejc, no_downstream_ejc, ref_atg_lost) | `05r_ref_atg_analysis.R` (via `Isopair::traceReferenceAtg()`); subgroup class table built inline in the canonical Rmd `sec3-class-table` chunk (no longer via the removed `05l_unified_model.R`). | **Verifiable locally.** |
| 4.42 | "5'UTR longer for NMD+ at Section C (ref-AUG-projected, n=819)" — Figure 4 Panel C | New Rmd §3a (`sec3a-panelC-tables`) + Figure 4 Panel C. Uses ref-AUG-projected scan from `05k_b_utr5_refaug.R` (NEW, TD2-bias-free) — replaces the legacy TD2-anchored scan from `05k_utr5_all_isoforms.R`. | **Verifiable locally.** |
| 4.43 | "Longer 5'UTR uORFs in NMD+ at Section C" — Figure 4 Panel D | New Rmd §3a (`sec3a-panelD-tables`) + Figure 4 Panel D. | **Verifiable locally.** |
| 4.44 | Helper / package functions backing 4.38–4.43 | `isopair:R/uorf.R` (uORF detection); `nmd:figures/lib/mechanism_class.R` (derived helper, replaces legacy cached column); `nmd:results/.../isopair_wrapper/05k_b_utr5_refaug.R` (ref-AUG-projected 5'UTR features). | **Verifiable locally.** |

### Paragraph 6 — CDS / 3'UTR sanity check on GENCODE-CDS-restricted scope (CDSand3UTR_GENCODEonly supplement)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 4.45 | CDS length comparison at n=130 (strict GENCODE-restricted scope) | New Rmd `sec2d-cds` chunk + `figures/SupplementalFigures/SF33-SF36_CDSand3UTR_GENCODEonly/data_export.R`. 2×3 supplement: row 1 = n=130 (panels A CDS / B 3'UTR translation-based / C 3'UTR bias-corrected); row 2 = n=819 ref-AUG-projected (panels D/E/F). | **Verifiable locally.** Replaces the old "n=1,904 / 22,335" claim (which was at the pre-refactor scope). |
| 4.46 | 3'UTR length comparison — translation-based vs bias-corrected (non-PTC-stop), demonstrating that the bias-correction pattern replicates across both n=130 strict and n=819 broad scopes | Same. The 2×3 supplement shows the same bias-correction collapse pattern at both scopes, validating the methodological choice. | **Verifiable locally.** |

### Section 4 figure render-script summary

In contrast to §1–3, **all §4 figure render scripts are local.** Figures 3 and 4 + the three supplements live under `figures/multipanel/` and `figures/SupplementalFigures/`; the new canonical Rmd `05_final_report_gencode_scope_2026-07-11.Rmd` embeds each panel via dedicated chunks. The Isopair package provides:

- `isopair:R/visualization.R::plotIsoformPair()` — used for Figure 2E (Yul-side §3) and any isoform-structure renders.

Figure 3 panels A–F are in `figures/multipanel/figure3_isopair_and_ptc/`; Figure 4 panels A–D are in `figures/multipanel/figure4_ptcneg_and_model/`. The three supplements are in `figures/SupplementalFigures/{CDSand3UTR_GENCODEonly, TD2BiasEvidence, PairAnalysisFlowchart}/`. Figure 5 (DL model) is in `figures/multipanel/figure5_dl_model/` — panels A–E exist; composite + RATIONALE.md pending (Task #43).

### Section 4 verifiable-locally summary

**All of §4 is locally verified** (status 2026-06-15). Both verifiers PASS:

- `reproducibility/verify_pass7_new_rmd.R` (37 checks) — pins every figure number to its canonical source TSV
- `reproducibility/verify_cross_check_new_rmd_vs_figures.R` (57 checks) — binds Rmd HTML claims to figure-side TSVs

Pipeline + Rmd are in `~/claude_projects/Isopair/` and `~/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/`. Headline verification anchors:

- **4.3** (1,548 pop_BC genes) — pinned by verifier
- **4.17** (36.9% vs 1.5% / 24-fold / p<10⁻¹³ at n=130) — pinned
- **4.21–4.23** (mechanism breakdown 52% / 38% / 10%) — pinned
- **4.24–4.25** (Panel E Fisher tests: SE 48%/14% p<10⁻⁶, A5SS 13%/4% p<0.05) — pinned
- **4.28** (52% TD2 == ref AUG at n=819) — pinned
- **4.30–4.31** (348 occult-PTC, 99% TD2 downstream, 82% ref-AUG Kozak stronger) — pinned
- **4.37** (90% = 756/819 cumulative PTC) — pinned
- **4.40** (3'UTR medians 1311 / 834.5 / 948 at n=819) — pinned

Re-run verifiers from repo root: `Rscript reproducibility/verify_pass7_new_rmd.R` and `Rscript reproducibility/verify_cross_check_new_rmd_vs_figures.R`. Both should report PASS.

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
| Architecture schematic (supplementary / Methods, NOT a Figure 5 panel) | `model:make_architecture_figure.R` | Architecture schematic |
| Figure 5 SHAP panels (legacy upstream) | `model:make_shap_interpretation_figure.R` | Legacy DeepSHAP-based interpretation panels — now rebuilt as `figures/multipanel/figure5_dl_model/figure5_panel{A,B,C,D,E}_*.py` |
| Analysis report | `model:orf_model_report_v5.Rmd` | End-to-end analysis report consuming the TSV exports above |

### Paragraph 1 — model architecture + input encoding + train/test split

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 5.1 | "Five candidate ORFs per transcript, priority: (1) reference CDS ORF when dominant non-NMD isoform's start is present; (2) TD2 CDS when different from (1); (3) remaining slots filled by top Kozak-scoring ORFs from ORFik scan" | `model:data_prep.py` + `model:METHODS.md`. Implementation in data_prep.py — ORF ranking logic. | **Verifiable locally.** |
| 5.2 | "9-channel encoding per position: A/C/G/T (4) + EEJ position + 50-bp rolling GC + 3 reading-frame channels relative to AUG" | `model:data_prep.py` + `model:METHODS.md`. | Verifiable locally. |
| 5.3 | "5 per-ORF structural features: frac_start, frac_stop, is_ref_cds, is_sqanti_cds, n_downstream_ejc" | `model:data_prep.py` + `model:METHODS.md`. | Verifiable locally. |
| 5.4 | "Two shared-weight CNN branches (AUG window, stop window) + linear branch (structural) → fused per-ORF embedding × 5 ORFs → attention aggregation → classification head" | `model:model.py`. | Verifiable locally. |
| 5.5 | "Chr 2, 4 for val; Chr 1, 3, 5, 7 for test; paralogs excluded from test set" | `model:data_prep.py` + paralog list in the `NMD_orf_model_v5_4ct` repo. | **Verifiable locally.** The paralog exclusion is the model's frozen H5 `test_paralog` split (linked repo); the isopair-wrapper generator `05u_paralog_annotation.R` was removed 2026-07-11 (its output fed only the now-frozen legacy report). |
| 5.6 | "BCEWithLogitsLoss with pos_weight = n_neg/n_pos; Adam (lr=1e-3, differential weight decay 1e-3 CNN / 1e-4 elsewhere); batch=256; fp16 AMP; ReduceLROnPlateau factor=0.5 patience=5; early-stop patience=10 on val AUC; max 100 epochs; seed=42" | `model:03_train.py` + `model:config.yaml`. | Verifiable locally — open `config.yaml` and confirm values match Methods text. |
| 5.7 | "Sweep over AUG ∈ {100, 500, 1000} × stop ∈ {100, 500, 1000, 2000}; selected AUG=500 / stop=500 by held-out AUC" | Sweep driver script (TBD — `03_train.py` may run individual configs while a sweep orchestrator like SLURM array job picks combinations). The selected config is in `model:config.yaml`. | Cite-only for the sweep results table (config-by-config AUC); the chosen config is verifiable locally. |
| 5.8 | Architecture schematic (likely supplementary / Methods illustration; NOT a panel in Figure 5 per the 2026-06-15 5-panel reframing) | `model:make_architecture_figure.R`. Final placement TBD. | **Verifiable locally.** |

### Paragraph 2 — performance + block-level importance (Figure 5 Panels A, B, C)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 5.9 | "AUC=0.93, AUPRC=0.833 on held-out test set" | `model:evaluate.py` — produces test-set AUC / AUPRC. The canonical AUPRC is 0.833 (4ct model), NOT 0.781 (predecessor v5 mirror); see [[feedback_nmd_orf_model_4ct_canonical]]. | **Verifiable locally** (assuming trained-model weights are in `~/claude_projects/NMD_orf_model_v5_4ct/results_4ct/` or cluster equivalent). |
| 5.10 | **Figure 5 Panel A** — ROC / PR curve at AUC=0.93 / AUPRC=0.833 | `figures/multipanel/figure5_dl_model/figure5_panelA_roc_curve.py`. (Per the 2026-06-15 figure numbering, the DL-model figure is Figure 5, not Figure 6.) | **Verifiable locally.** |
| 5.11 | "Roughly ⅔ of predictive information from ORF structural data, ⅓ from START+STOP sequence" — Shapley block-level decomposition | `model:11_kernel_shap_branches.py` — KernelSHAP at embedding level, exact additive Shapley across the 3 sub-encoders (AUG / STOP / structural). | **Verifiable locally.** |
| 5.12 | "STOP sequence ~3× as important as START sequence" | Same `11_kernel_shap_branches.py` — ratio of STOP to START sub-encoder attributions. | Verifiable locally. |
| 5.13 | **Figure 5 Panel B** — block-level Shapley decomposition (3-branch) | `figures/multipanel/figure5_dl_model/figure5_panelB_branch_importance.py`. | **Verifiable locally.** |
| 5.14 | "Of individual ORF structural features, EJC count was by far the most important" | `model:05_interpret_structural.py` + `model:06_export_deepshap_tsv.py` (per-feature DeepSHAP). | Verifiable locally. |
| 5.15 | **Figure 5 Panel C** — ranked per-feature structural importance | `figures/multipanel/figure5_dl_model/figure5_panelC_structural_features.py`. | **Verifiable locally.** |

### Paragraph 3 — nucleotide-level attributions (Figure 5 Panels D, E + SFx)

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 5.16 | "Importance concentrated within 50 nt preceding the start or stop site" (SFx — Shapley value profile across windows) | `model:06_export_deepshap_tsv.py` — per-position attribution. | Verifiable locally. |
| 5.17 | "Both windows: more importance for NMD than Control isoforms" | Same — stratify per-position attributions by NMD label. | Verifiable locally. |
| 5.18 | "START window: model learned Kozak sequence; no weight on AUG itself because invariant across training ORFs" | `model:07_motif_analysis.py` + `scripts/export_joint_motif_logos.py`. | Verifiable locally. |
| 5.19 | **Figure 5 Panel D** — nucleotide-level attribution around AUG (Kozak motif logo) | `figures/multipanel/figure5_dl_model/figure5_panelD_atg_logo.py`. | **Verifiable locally.** |
| 5.20 | "STOP window: UGA shifts predictions toward NMD; U at +4 has largest importance for that position (matches readthrough biology)" | `model:07_motif_analysis.py` + driver. **This is the headline biological finding for §5.** | **Verifiable locally — high priority verification target given the lit-review §6 dependence on this claim.** |
| 5.21 | **Figure 5 Panel E** — nucleotide-level attribution around stop codon (UGA + U+4 highlighted) | `figures/multipanel/figure5_dl_model/figure5_panelE_stop_logo.py`. | **Verifiable locally.** |

### Paragraph 4 — attention distribution

> **Open question (§5 figure-naming).** Paragraph 4 in the manuscript references "Fig X panel A/B" for the attention-distribution comparison (NMD vs Control). Per the 2026-06-15 reframing, these attention panels are NOT in the 5-panel Figure 5 (A=ROC, B=branch importance, C=structural features, D=ATG logo, E=stop logo). They likely belong to a supplementary attention figure. Pinning still pending.

| # | Claim | Source code | Notes / status |
|---|---|---|---|
| 5.22 | "ORFs numbered by CDS likelihood; ORF0 = most likely CDS; most attention on ORF0" | `model:04_interpret_attention.py` — exports per-transcript attention vectors over 5 ORFs. | Verifiable locally. |
| 5.23 | "Attention more broadly distributed in NMD than in Control isoforms" | Same — entropy or top-1-share of attention vector, stratified by NMD label. | Verifiable locally. |
| 5.24 | **"Fig X panel A and B"** (attention distribution NMD vs Control) — supplementary, NOT in Figure 5 | Render likely in `model:make_shap_interpretation_figure.R` or `model:orf_model_report_v5.Rmd`. Final panel identity TBD. | Verifiable locally once panel identity is pinned down. |

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

Local. The current Figure 5 (5 panels A–E) is rendered from `figures/multipanel/figure5_dl_model/figure5_panel{A,B,C,D,E}_*.py` against TSV exports from `model:` (the linked DL repo). Architecture schematic (`model:make_architecture_figure.R`) is supplementary or Methods illustration, NOT a Figure 5 panel. The attention panels in Paragraph 4 are a separate (supplementary) figure — final panel identity still TBD. The composite render + RATIONALE.md for Figure 5 is the remaining open task (#43).

### §4/§5 Supplemental Figures — 2026-07-08 build

This block records the §4 and §5 supplemental-figure work started 2026-07-08 (post Yul's paper draft), adding new SFs and porting existing SFs to the ggplot-mimic matplotlib style. Every SF entry names the render script, its input TSVs, and the upstream code/commit that produced those TSVs. The traceability audit for §4/§5 SFs runs from these entries.

**Manuscript SF numbering was corrected 2026-07-08** to insert a new SF34 slot (TD2 CDS length + downstream AUG); everything from old-SF34 onward shifted +1. Figure-side dirs and filenames now match the corrected manuscript numbering.

| SF number in manuscript | Render script | Input data | Upstream (data provenance) |
|---|---|---|---|
| SF25 (Twelve splice event categories detected by Isopair) | `nmd:figures/SupplementalFigures/SF25_SpliceEventCategories/figure_sf25_splice_event_categories.py` | (diagrammatic — coordinates embedded in the script) | Canonical event definitions and illustrative coordinates from `nmd:results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_mashr.Rmd`, § Event Type Glossary. Isopair implementation lives in the linked `peter4244/Isopair` R package. |
| SF27 (Isoforms per gene, n=1,548 pair-set cohort) | `nmd:figures/SupplementalFigures/SF27_IsoformsPerGene/figure_sf27_isoforms_per_gene.py` | `data/isoforms_per_gene.tsv`, `data/descriptives_summary.tsv` | Copied from `SupplementalFigures/PairSetDescriptives/data/` (canonical Rmd source: `05_final_report_gencode_scope_2026-07-11.Rmd`) |
| SF28 (Reference share of gene expression, n=1,548 pair-set cohort) | `nmd:figures/SupplementalFigures/SF28_ReferenceShare/figure_sf28_reference_share.py` | `data/ref_expression_fraction.tsv`, `data/descriptives_summary.tsv` | Copied from `SupplementalFigures/PairSetDescriptives/data/` (canonical Rmd source: `05_final_report_gencode_scope_2026-07-11.Rmd`). Reference isoform = highest-expressed non-NMD isoform per gene under DMSO, 4-CT mean. Median 64.3%; 67.5% at ≥50%. |
| SF30 (GAIN direction by splice event type, n=1,548 pair-set triplets) | `nmd:figures/SupplementalFigures/SF30_GainDirectionByEvent/figure_sf30_gain_direction_by_event.py` | `data/table1c_gain_loss_events.csv` | Isopair pipeline canonical table generated by `nmd:results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_mashr.Rmd`, § "Gain/Loss Direction by Event Type" (chunk `goal1-table1c-gain-loss`). Fisher's exact per event type on the 2×2 (GAIN/LOSS × NMD/Control) counts. |
| SF31 (NMD response magnitude vs stop→last-EJC distance, n=2,790 gene-matched C2 comparators) | `nmd:figures/SupplementalFigures/SF31_PTCDistanceDoseResponse/figure_sf31_ptc_distance_dose_response.py` | `data/sf31_ptc_distance_logfc.tsv` | Produced by `data_export.R` (same dir) from `profiles_c2_allsamples.rds` + `ptc.rds` + `mashr_isoform_model_2026.3.10.rds`. Matches the Rmd chunk `goal1-fig2-logfc-dist` in `05_final_report_mashr.Rmd` § "NMD Response vs Stop Codon Distance". Spearman ρ = 0.182, p = 3.9 × 10⁻²². |
| SF32 (NMD effect size by downstream EJC count, PTC+ subset n=1,211) | `nmd:figures/SupplementalFigures/SF32_NMDEffectByEJCCount/figure_sf32_nmd_effect_by_ejc_count.py` | `data/sf32_ejc_count_logfc.tsv` | Produced by `nmd:figures/SupplementalFigures/SF31_PTCDistanceDoseResponse/data_export.R` (same run as SF31). Matches the Rmd chunk `goal1-fig3-ejc-boxplot` in `05_final_report_mashr.Rmd` § "NMD Strength by Downstream EJC Count". PTC+ subset (`has_ptc == TRUE`) of SF31 population. |
| SF29 (Transcript length by pair role, n=1,548 triplets) | `nmd:figures/SupplementalFigures/SF29_TranscriptLengthByRole/figure_sf29_transcript_length_by_role.py` | `data/tx_length_by_role_long.tsv`, `data/descriptives_summary.tsv` | Copied from `SupplementalFigures/PairSetDescriptives/data/` (canonical Rmd source: `05_final_report_gencode_scope_2026-07-11.Rmd`). Medians: NMD comparator 3,018 nt, Reference 2,896 nt, Control comparator 2,742 nt. |
| SF37 (Total \|SHAP\| across AUG + stop windows) | `nmd:figures/SupplementalFigures/SF37_ShapAcrossWindows/figure_sf37_shap_across_windows.py` | `data/shap_profile_{atg,stop}_joint_atg500_stop500.tsv` | `model:09b_export_subgroup_profiles.py` (sbatch `model:slurm_export_subgroup_profiles_09b.sh`), tag `atg500_stop500`. NPZ ancestry: `slurm_deepshap_seq_500bg.sh` 5-run array. Cluster export commit `model:5c19591`. |
| SF37 supplemental (end-of-window channel decomposition) | `nmd:figures/SupplementalFigures/SF37_ShapAcrossWindows/edge_decomposition_notes.md` (analysis note) | `data/sf37_edge_channel_decomposition.tsv` (90 rows: window × zone × channel) | Same as SF37 input; decomposition tabulated locally. Right-edge stop-window elevation dominated by T channel (consistent with poly-A signal machinery: AAUAAA + downstream U-rich element). Left-edge AUG elevation dominated by G/C nucleotide channels with a zero-padding contribution in short-5′UTR transcripts; frame channels are not the drivers at either edge (corrects an earlier assumption). |
| SF42 (Rolling GC content across AUG + stop windows) | `nmd:figures/SupplementalFigures/SF42_GCcontentStopWindow/figure_sf42_gc_content_stop_window.py` | `data/gc_content_across_{atg,stop}_window_refaug_only_atg500_stop500.tsv` | `model:09d_export_gc_content_refaug_only.py` (sbatch `model:slurm_export_gc_refaug_09d.sh`), tag `atg500_stop500` run 1. Restricts NMD to `is_ref_cds=1` on ORF0 (n = 1,743 / 2,268 = 76.9%); Controls unfiltered (n = 7,863). Cluster export commit `model:4476c7f`. |

**Pending §4/§5 SFs (tracked here at planning time; entries will be filled in as each ships):**

| Manuscript SF | Content | Status |
|---|---|---|
| SF26 (Isopair pair-analysis cohort flowchart) | Existing dir `SupplementalFigures/SF26_PairAnalysisFlowchart/` (Graphviz render, not matplotlib) | Existing render; legend format check pending (Task #148); no matplotlib restyle applicable |
| SF33 (CDS/3'UTR at n=130 GENCODE-restricted; cited by manuscript prose "no significant difference in CDS or 3'UTR ... in NMD+/PTC+ vs NMD+/PTC-") | Existing dir `SupplementalFigures/SF33-SF36_CDSand3UTR_GENCODEonly/` (matplotlib); row 1 of that figure at n=130 | Existing render restyled to ggplot-mimic 2026-07-08 (Task #148 in-place). Same dir also serves SF36 (row 2 at n=1166). |
| SF34 (TD2 CDS length + downstream AUG position vs ref-AUG-anchored ORF; cited by "the TD2 CDS was usually longer than the reference AUG-tranced CDS and used an AUG downstream from the reference AUG") | Existing dir `SupplementalFigures/SF34-SF35_TD2BiasEvidence/` (matplotlib); length + position panels A + C (broad n=819) and D + F (occult n=348 subset) | Existing render restyled to ggplot-mimic 2026-07-08 (Task #148 in-place). Same dir also serves SF35 (Kozak panels B, E). |
| SF35 (Kozak PWM at ref-AUG vs TD2-called AUG; cited by "reference AUG was stronger in 82% of cases (286/348)") | Existing dir `SupplementalFigures/SF34-SF35_TD2BiasEvidence/` (matplotlib); Kozak panels B (broad n=819) and E (occult n=348) | See SF34 row. |
| SF36 (CDS/3'UTR at n=819 ref-AUG-traceable; cited by manuscript prose about 3'UTR length in GENCODE-restricted vs expanded sets, medians 1324/722/922 nt for NMD+/PTC+ / NMD+/PTC- / Control) | Existing dir `SupplementalFigures/SF33-SF36_CDSand3UTR_GENCODEonly/` (matplotlib); row 2 of that figure at n=819 | See SF33 row. |
| SF38 (UGA usage NMD vs Control) | Existing dir `SupplementalFigures/SF38_StopCodonUsage/` (matplotlib) | Existing render; restyle pending (Task #153) |
| SF39 A/B (Attention distribution NMD vs Control) | Existing dir `SupplementalFigures/SF39_AttentionDistribution/` (matplotlib) | Existing render; restyle pending (Task #153) |
| SF40 (Branch SHAP by PTC subclass) | Existing dir `SupplementalFigures/SF40_PTCSubclassBranchSHAP/` (matplotlib) | Existing render; restyle pending (Task #153) |
| SF41 (Discrimination performance by PTC subclass) | Existing dir `SupplementalFigures/SF41_PTCSubclassPerformance/` (matplotlib) | Existing render; restyle pending (Task #153) |
| SF43 (NMDetective-B / NMDEP comparison) | Render: `SupplementalFigures/SF43_ModelComparison/` (matplotlib). **Upstream analysis: `nmd:code/nmd_predictor_comparison/`** — `01_extract_our_isoforms.R` → `02_score_nmdetective_b.R` / `03_score_nmdep_rule_baseline.R` → `04_compute_metrics.R`, report `nmd_predictor_comparison.Rmd`, methods `METHODS.md`, cluster inference under `explorer_run/`. Outputs `metrics_summary_2026.7.11.tsv`, `per_isoform_scores_2026.7.11.tsv`. | Existing render; restyle pending (Task #153). **Producer mapped 2026-07-21** — this directory previously had **zero** map citations despite producing SF43 and being the source for the manuscript's NMDetective-B / NMDEP comparison (cited 4×). It is the third demonstration that the map is not a keep-list: uncited ≠ unused. |

### Section 5 verifiable-locally summary

**Essentially all of §5 is locally verifiable**, with these caveats:

- **Training itself requires GPU + HDF5 input.** Pete already has the trained-model weights cached in `model:results_4ct/`; verification of AUC/AUPRC is from running `evaluate.py` against those weights, not retraining.
- **KernelSHAP and DeepSHAP runs already executed.** The TSV exports under `model:results_4ct/` (per the `06_export_deepshap_tsv.py` and `11_kernel_shap_branches.py` outputs) are what feed the figure-render scripts. Verification is reproducing the chain from TSV → figure.
- **Cross-repo subgroup dependency:** the PTC+ / PTC− ref-AUG retained / lost labels come from the Isopair wrapper (`05r_ref_atg_analysis.R`). Verification of subgroup attributions requires that those labels match between the two repos.

Headline biological-claim verification targets:
- **5.9** AUC=0.93 / AUPRC=0.833 (4ct model — see [[feedback_nmd_orf_model_4ct_canonical]])
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

### Unmapped Rmd / script gaps (Yul-side)

5. **DMSO-only one-vs-rest cell-type-marker Rmd** (Methods M5; Section 1 numbers 1.9–1.14) — code path TBD.
6. **Tan et al. mashr reanalysis Rmd** (Methods M9; Section 2 claim 2.18) — code path TBD; likely in Yul's `final/`.
7. **Sequencing-summary aggregation script** (Section 1 claim 1.1: "26 samples / 13M FLNC reads / 341,638,920 total reads") — code path TBD.
8. **Section 5 sweep results table** (claim 5.7: AUG × stop window-size sweep results table) — orchestrator script and recorded results TBD.

### Figure-render-script gaps (Yul-side)

9. **Figure 1 panels** (A, C, D, F) — render scripts TBD; nominally Yul-side.
10. **Figure 2 panels** — ⚠️ **re-scoped 2026-07-24.** Figure 2 is now A–E, not A–F. **2A** producer identified: `yul:transcriptional_output.Rmd` chunk `barplot-global`. **2B** (productive logFC distribution) TBD in `yul:productive_response.Rmd`. **2C/2D/2E** (*SHMT2*, *SRSF2*, *PCNA* structure + logFC) rendered via `isopair:R/visualization.R::plotIsoformPair()` — locally available. The former **GPR180** 2E/2F panels are retired with the PCI analysis.

### Methods-text / claim verification flags (Yul-side)

11. ~~**SR↔LR coverage threshold** (old §3 claim 3.10, "~40% lacked sufficient LR isoform coverage")~~ — **RETIRED from §3 2026-07-24** (the SR↔LR concordance paragraph is no longer in §3). ⚠️ Before closing, check whether this claim migrated to §1/§2; if it did, the undefined "sufficient coverage" threshold still needs pinning from `yul:comparison_analysis.Rmd` (now local).
12. ~~**PCI 63–90% denominator** (old §3 claim 3.7)~~ — **RETIRED 2026-07-24**; the PCI analysis is not in the manuscript. Superseded by §3 finding **F-3** (the "% output lost" metric is a total variation distance) as the open §3 methods question for Yul.

### Figure-numbering inconsistencies

13. **§4 paragraph 5 panel-letter orphan** — manuscript still says "Fig X panel A/B/C/D" for the 5'UTR / 3'UTR / uORF subgroup panels. Under the 2026-06-15 reframing these map to Figure 4 A/B/C/D + the CDSand3UTR_GENCODEonly supplement. Manuscript text needs find/replace to make panel references explicit; current §4 audit (Pete applied 2026-06-16) handled the headline numbers but not the panel letters.
14. **§5 attention panels overlap** (claim 5.24) — paragraph 4 cites "Fig X panel A and B" for attention distribution. NOT in the 5-panel Figure 5 (A=ROC / B=branch importance / C=structural features / D=ATG logo / E=stop logo). Likely a supplementary attention figure. Final panel identity TBD.

### Resolved since 2026-06-15

- ✓ Pete's `isocall_limma_dge_fullmodel_2026.3.1.Rmd` is QC-only — confirmed per the explicit Pete vs Yul ownership model in `isocall_dge/limma/README.md` (Yul owns limma + mashr DGE; Pete owns Isopair + DL model). Pete's parallel `pete/` snapshot was deleted on 2026-06-16.
- ✓ "85% explained by PTCs" (old claim 4.37) — superseded by the current 90% (756/819) at the n=819 expanded scope; pinned by the cross-check verifier.
- ✓ Figure 3 / Figure 4 panel-drop / renumbering — final 6-panel Figure 3 (A–F) + 4-panel Figure 4 (A–D) per the 2026-06-15 plan; methodology files, RATIONALE.md, and find/replace pairs all in place.

---

## Next steps (after this map)

In rough priority order:

1. **Yul-repo access setup** — Yul to add `repjc` as a member of her changit project (or `scp` `final/` from Channing). Unblocks §1, most of §2, and all of §3.
2. **Build Figure 5 composite + RATIONALE.md** (Task #43) — the 5 panels exist; composite render + pre-registration doc pending.
3. **Manuscript-wide TD2 bias audit** (Task #49) — `TD2_BIAS_AUDIT.md` at repo root is opened; needs to be carried through to identify any §1/§2/§3 analyses still consuming TD2-derived CDS boundaries.
4. **Pin §1b descriptives** — manuscript §4 ¶1 still cites "7 isoforms/gene" + "70% of parent gene expression" + "75% accounting for >50%" from the legacy Rmd; these aren't in the current verifier suite. Pete deferred recomputation 2026-06-16; worth adding to the verifier if they get revisited.
5. **Decide on legacy Rmd removal** — `05_final_report_mashr.Rmd` (317KB) still on disk with deprecation banners. Could be deleted once we're confident nothing else depends on it.
6. **Yul-side verification** — once Yul access is in place, run the 5-step scientific-report protocol on §1, §2, §3 against Yul's Rmds.

---

*Revised 2026-06-16. §4 now pinned to the new canonical Rmd `05_final_report_gencode_scope_2026-07-11.Rmd` + both verifiers (pass-7 + cross-check). 144 numbered claims; §4 fully verified locally.*
