# ONBOARDING — NMD long-read paper repo (Castaldi + Randell labs)

**Audience:** Claude Code (or any AI coding assistant) joining this project on behalf of a co-author. If you are a human, you can read this too — but the level of detail is calibrated for an AI agent that has just cloned the repo and needs to be productive on its first task.

Repo: `peter4244/nmd_isopair_analysis` (public). Canonical local path conventions assume `~/claude_projects/nmd/`.

---

## 1. Project mission and current phase

This is the analysis repo for the NMD long-read paper:

> Leshem, Kasai, Thakur, Paul, Ziniti, Boueiz, Saferali, DeMarzio, Laederach, Randell, Castaldi. **Long-read RNA sequencing in primary lung cell types reveals principles of nonsense-mediated decay.** (In preparation, May 2026.)

**Current phase:** manuscript finalization. The biology is settled. Active work is split between (a) finalizing the manuscript draft itself in a shared Google Doc, and (b) finalizing figures, tables, and methods text in the repo. Pete and Yul are the two principal users of the repo. The two cluster homes are Channing (BWH) for Pete and Northeastern Discovery (`/home/p.castaldi/cc/`) for Yul.

**Headline scientific claims** (so you can sanity-check your contributions against them):

- NMD in primary human lung cells is widespread, isoform-level, cell-type-invariant quantitative regulation, not a binary disposal pathway. ~90–92% of significant isoforms are NMD-responsive across the 4 cell types in scope.
- The deep-learning model (CNN + ORF features hybrid) reaches AUC ≈ 0.93 / AUPRC ≈ 0.83 on held-out chromosomes.
- The headline methodological finding: **NMD prediction is bottlenecked by ORF identification, not by feature engineering.** Length-priority and splicing-aware ORF callers give indistinguishable AUC, but only the splicing-aware caller attributes calls to the correct biological mechanism (PTC vs spurious "long uORF"). Quantitatively, splicing-aware ref-ATG tracing lifts the PTC-attributable fraction of gene-matched NMD pairs from 43% to 77%.

---

## 2. Cell-type scope (CRITICAL)

The lab generated long-read + short-read data on six primary lung cell types under SMG1 inhibitor (Smg1i) vs DMSO control:

| Code | Cell type |
|---|---|
| **AT** (sometimes AT2) | Alveolar type 2 cells |
| **DD** | Day-differentiated large airway epithelial cells (LAE), submerged 2D primary culture |
| **DD_ALI** | LAE at air-liquid interface |
| **DO** | Donor cells / primary SAEC, submerged |
| **DO_ALI** | Donor cells at ALI = **SAEC ALI** |
| **FB** | Lung fibroblasts |
| **MV** (or MVE) | Microvascular endothelial cells |

**The manuscript is scoped to four cell types: AT, DD, FB, MV.** DD_ALI was excluded due to near-zero SR-vs-LR logFC correlation (r = 0.002). DO_ALI was excluded for similar reliability reasons (low effect-size correlation; small donor n). Note that DO_ALI is still important for grant-related work (it's the SAEC ALI preparation the NMD 2026 grant's Aim 1 uses; ATF4 is strongly NMD-targeted there at adj.P = 1.7×10⁻¹⁷), so per-cell-type gene-level mashr CSVs for all 6 cell types still exist in the repo.

**Cell-type label rename (April 2026):** `AT → AT2` and `DO → DO_ALI` in the canonical labeling. Files pre-dating commit `1cd5328` may still use the older `AT` / `DO` labels — be careful when joining across vintages.

---

## 3. Repository layout (and what is and isn't tracked)

```
~/claude_projects/nmd/
├── paper/                              # Manuscript draft (markdown export of Google Doc)
│   └── NMD manuscript 2026.2.5.md     # Current. THE GOOGLE DOC IS THE SOURCE OF TRUTH.
├── code/                               # Analysis scripts (124 tracked R/Rmd/py files)
├── pheno/                              # Sample / donor metadata
├── shortread_dge/mashr/                # Gene-level mashr DGE per cell type
│   ├── nmd_mashr_dge_{at,dd,ddali,doali,fb,mv}_2026.3.10.csv   # canonical per-CT
│   └── mashr_{lfsr,posterior_means,...}_2026.3.10.csv          # shared 4-CT
├── isocall_dge/mashr/                  # Isoform-level mashr DIE (4 CT only: at, dd, fb, mv)
├── isocall_dge/old/mashr/              # Archived 6-CT isoform mashr (DO_ALI still in here)
├── results/                            # gitignored EXCEPT for one subtree (see below)
│   └── isoform_transitions/Version_6.0/isopair_wrapper/   # force-tracked
│       ├── 01_..05_*.R                    # The Isopair analysis pipeline
│       ├── 05_final_report_mashr.Rmd      # Main analysis report
│       └── paper_outputs/                 # Older manuscript drafts (now in paper/)
├── tmp/                                # gitignored. Regenerable analysis outputs land here.
├── data/                               # gitignored. Raw counts.
├── resources/, sqanti/, reference_files/  # gitignored. Large reference data.
├── CLAUDE.md, README.md                # Project-level docs (CLAUDE.md is partly stale)
└── ONBOARDING.md                       # This file.
```

Key gitignore rules from `.gitignore`:

- **`results/` is gitignored** — only `results/isoform_transitions/Version_6.0/` is force-tracked. **Anything you save under any other `results/` subdirectory is invisible to version control.** Don't put work there expecting it to persist.
- **`tmp/` is gitignored.** Conventional place for analysis outputs (TSVs, intermediate RDS files, GSEA tables). Regenerable from code.
- **`*.csv` under `longread_dge/`, `shortread_dge/`, `isocall_dge/mashr/`** are gitignored individually as well — DGE outputs are large and regenerable.
- `.rds`/`.RDS`/`.RData`, R cache dirs, `archive/` trees, lit-review PDFs, generated HTML reports are all gitignored.

---

## 4. THE BIG PUBLICATION WARNING ⚠️

This repo has a **dual-push** setup. The local `origin` remote has two push URLs:

1. `changit.bwh.harvard.edu:repjc/nmd_lungcells_2026.git` (Channing GitLab, private)
2. `https://github.com/peter4244/nmd_isopair_analysis.git` (**PUBLIC GitHub**)

**Every `git push` goes to both at the same time.** A commit on `main` becomes publicly visible on GitHub instantly. There is no staging gate. Before committing anything that contains:

- patient identifiers,
- unpublished phenotype tables,
- credentials, API keys, or cluster paths with embedded usernames anyone could grep for,

…verify it's safe for public release. The `.gitignore` is already broad to prevent accidental leaks of phenotype CSVs and large data files, but it does not check the contents of new files you create.

If you are working on something you would prefer to review locally first, work on a feature branch and don't push it until ready.

---

## 5. Linked repos — the NMD work spans multiple GitHub repos

Don't conflate them.

| Repo | Purpose | Public? |
|---|---|---|
| `peter4244/nmd_isopair_analysis` | THIS REPO. Analysis pipeline, manuscript markdown, code. | ✓ |
| `peter4244/Isopair` | The Isopair R/Bioconductor package itself (gene-matched isoform-pair analysis, splicing-event detection, PTC attribution). | ✓ |
| `peter4244/NMD_orf_model_v5_4ct` | **Canonical source for the manuscript's deep-learning model** (training code, config, METHODS.md, DeepSHAP scripts, SLURM wrappers). The model was trained on Northeastern Discovery; this repo is the version-controlled source. | ✓ |
| `peter4244/NMD_orf_model_v5` | Predecessor of v5_4ct, pre-4-cell-type-scope migration. | ✓ |
| `peter4244/Isoscope` | Per-gene isoform annotation/visualization (long-read + short-read). | ✓ |
| `peter4244/Isovar` | Variant → isoform splicing functions for COPD GWAS sQTL work. | ✓ |
| `peter4244/copd-nmd-sqtl-airway-epithelial` | Shiny app for the GWAS-sQTL × NMD integration. | private |
| `peter4244/nmd-2026-grant` | NMD R01 grant materials. | private |

The manuscript's "Isopair" methods text refers to the Isopair package (v0.99.2 at time of writing); the "deep learning model" methods text refers to `NMD_orf_model_v5_4ct`'s `03_train.py` and `METHODS.md`. If you're writing or auditing methods prose, **pull from those repos' own METHODS.md / vignettes as the source of truth**, not from anything in this repo.

---

## 6. Methodological conventions (the non-obvious rules)

These are things that have bitten people before. Internalize before writing analysis code.

1. **Baseline expression is from DMSO samples only.** Never use limma's `AveExpr` or any pooled `cpm_*` column as a baseline — Smg1i inflates NMD targets and overestimates baseline. Filter your DGEList to `treatment == "DMSO"` before computing baseline statistics. Label baseline columns explicitly ("DMSO baseline CPM") so the convention is unambiguous downstream.

2. **SQANTI/TD2 CDS calls are biased against PTC-containing ORFs** — they preferentially pick non-PTC ORFs, which is the *opposite* of what an NMD experiment needs. For main-ORF PTC NMD substrates, use the reference-projected ORFs from `NMD_orf_model_v5_4ct` (e.g., `nmd_orf_model_v5_4ct/tmp/ref_orf_ptc_cache_with_nmd_2026.5.12.rds`). For uORF-mediated NMD (ATF4 is the canonical case), the SQANTI main-CDS is usually fine because the NMD-triggering termination is in the uORF, not the main ORF.

3. **NMD-responsive definition (manuscript convention):** `mashr lfsr < 0.05` AND `posterior_mean > 0`. Don't substitute `adj.P.Val < 0.05` from the limma step underlying mashr without checking — they're not the same denominator.

4. **Non-NMD definition (for Isopair pair construction):** `adj.P.Val > 0.30` in the limma step underlying mashr (this threshold was lowered from 0.50 to 0.30 in the 4-CT mashr refit because the new mashr distribution is tighter). Specific to gene-matched-pair analyses.

5. **Filename date-stamp convention:** `yyyy.m.d` (e.g., `2026.3.10`, `2026.5.18`). The date is when the file was *generated*, not when it was last edited. When you write a new analysis output, use today's date.

6. **Speculation vs observation discipline.** Pete enforces this on every scientific deliverable. In prose, *observation* language ("X is +1.19 in DO_ALI with adj.P = 1.7×10⁻¹⁷") must be clearly separated from *interpretation* language ("consistent with NMD stabilizing ATF4-axis transcripts," "we hypothesize…"). Mechanism / causal phrasing ("is a compensatory response to") needs the strongest hedging. When in doubt, hedge.

7. **The 5-step scientific-report verification protocol.** Before declaring any publication-quality report "done," run these as five separate focused passes — do **not** combine. The canonical anti-example is publishing OR = 0 for all domains because verification was combined with development.

   1. **Factual accuracy.** Recompute every stated number independently against the raw data; verify inline R values reference the correct objects.
   2. **Result correctness.** For each statistical test, check the 2×2 table construction, plausibility of ORs/p-values, denominators. Actively try to disprove each result.
   3. **Documentation accuracy.** Every prose statement matches its corresponding table/figure exactly.
   4. **Reproducibility & completeness.** A reader can trace every result to its source code + data. Every chunk cites its script, function, and METHODS section.
   5. **METHODS.md verification.** METHODS.md text matches the actual code implementations and current parameters.

---

## 7. Figure conventions (the immediate collaboration focus)

For publication / grant-quality figures the lab follows a specific set of conventions. The full principles document lives at `~/.claude/skills/figures/references/figures_principles.md` on Pete's machine (which is not in this repo — see follow-up note below). The headline rules:

- **Validator-clean output as the bar.** A figure isn't "done" until a layout-validator pass returns zero warnings.
- **Two font sizes only** in a figure; body text never bold.
- **Structural symmetry.** Panels of the same type share dimensions, margins, and named coordinate constants.
- **Reference-image-first** for biology panels — start from the published figure you're emulating.
- **Snapshot triplets** at close-but-not-final milestones (overview + detail + comparison) so Pete can review evolution.
- **Data-integrity rules:** every number in a figure is recomputable from a script in the repo; no manual values.
- **Empty space → compress, don't pad.** When Pete flags whitespace, the move is to shrink panels / trim margins / remove titles, not add content.

**Brainstorming / mockup stage = relaxed rigor.** At early stages skip validator-clean output and snapshot triplets; the figure is for team alignment. Bring full rigor back when the figure is heading into the submitted draft.

**Follow-up needed for shared figure styles:** the figure helpers (`figure_primitives.py`, `figure_geometry.py`, `validate_figure_layout.py`) currently live in `~/.claude/utils/` on Pete's machine — they are *not* in this repo. To share styles with Yul, those need to be either (a) copied into `figures/lib/` here and tracked, or (b) packaged as a separate shared `nmd_figures` repo. This is the next infrastructure task.

---

## 8. Manuscript workflow

- **Source of truth:** the shared Google Doc (Pete + Yul are co-editing).
- **Repo copy:** `paper/NMD manuscript 2026.2.5.md` is a markdown export periodically synced from the Google Doc. The repo copy gives Claude / scripts something to read and diff against. **Do not** assume the repo markdown is canonical for prose changes — when you want to suggest an edit to the manuscript, *return a list of manual find/replace pairs the user can apply in Google Docs*, not commits to the markdown.
- **Manuscript methods reconciliation:** the deep-learning model methods text in the manuscript should track `NMD_orf_model_v5_4ct/METHODS.md`. The Isopair methods text should track `Isopair/vignettes/NMD-attribution.Rmd`.

---

## 9. Common analysis tasks — pointers

| Task | Where to start |
|---|---|
| Per-cell-type gene-level mashr DGE | `shortread_dge/mashr/nmd_mashr_dge_{ct}_2026.3.10.csv` |
| Per-cell-type isoform-level mashr DIE (4-CT) | `isocall_dge/mashr/nmd_mashr_die_{at,dd,fb,mv}_2026.3.10.csv` |
| Per-cell-type isoform-level mashr DIE (6-CT archived) | `isocall_dge/old/mashr/nmd_mashr_die_{at,dd,ddali,doali,fb,mv}_2026.3.10.csv` |
| GSEA on March 2026 mashr | `code/gsea_mashr_2026.3.10.R` → `tmp/gsea_mashr_gene_2026.3.10_run2026-05-18.tsv` |
| Isopair gene-matched pair analysis | `results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_mashr.Rmd` |
| Deep-learning model training | `peter4244/NMD_orf_model_v5_4ct/03_train.py` (clone that repo separately) |
| Deep-learning model methods text | `peter4244/NMD_orf_model_v5_4ct/METHODS.md` |
| Isopair package vignette | `peter4244/Isopair/vignettes/NMD-attribution.Rmd` |
| Manuscript markdown | `paper/NMD manuscript 2026.2.5.md` |

---

## 10. If you only read three things

1. Section 4 — the public-by-default warning.
2. Section 6 — methodological conventions.
3. Section 8 — the manuscript-edits-as-find/replace-pairs rule.
