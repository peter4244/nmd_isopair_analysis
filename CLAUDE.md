# CLAUDE.md

This file provides guidance to Claude Code when working in this repository.

**First-time reader: start with [`ONBOARDING.md`](ONBOARDING.md).** It covers the project mission, public-by-default warning, cell-type scope, linked-repo separation, methodological conventions (DMSO-only baseline, SQANTI PTC bias, mashr conventions, speculation-vs-observation discipline, 5-step verification), figure conventions, and the manuscript-edits-as-find/replace-pairs rule. CLAUDE.md (this file) covers repo-internal navigation only.

## Project at a glance

NMD long-read paper repo. Manuscript phase (May 2026). Lead authors: Leshem, Kasai, ..., Randell, Castaldi. The paper characterizes nonsense-mediated decay across primary human lung cell types using paired short-read + PacBio Iso-Seq under SMG1 inhibition, and trains a CNN + ORF-features hybrid to predict NMD from sequence. Headline finding: NMD prediction is bottlenecked by ORF identification, not feature engineering.

**Manuscript scope:** 4 cell types — AT (alveolar type 2), DD (large airway epithelial, submerged), FB (fibroblast), MV (microvascular endothelial). DD_ALI and DO_ALI are in the data and per-cell-type mashr outputs but excluded from primary manuscript analyses due to low SR-vs-LR effect-size correlation.

## Repository structure (current)

```
nmd/
├── ONBOARDING.md                       # Start here for context
├── README.md                           # Public-facing repo description
├── CLAUDE.md                           # This file — navigation
├── paper/
│   └── NMD manuscript 2026.2.5.md     # Manuscript markdown export (Google Doc is source of truth)
├── code/                               # Analysis scripts (124 tracked R/Rmd/py files)
├── figures/
│   ├── README.md
│   └── lib/                            # Shared figure primitives, geometry, validators (see figures/README.md)
├── pheno/                              # Donor / sample metadata
├── shortread_dge/mashr/                # Gene-level mashr DGE per cell type (6 CTs available)
├── isocall_dge/mashr/                  # Isoform-level mashr DIE (4 CTs: at, dd, fb, mv)
├── isocall_dge/old/mashr/              # Archived 6-CT isoform mashr (DO_ALI, DD_ALI still here)
├── results/                            # gitignored
│   └── isoform_transitions/Version_6.0/   # FORCE-TRACKED exception; the Isopair analysis pipeline lives here
├── tmp/                                # gitignored; regenerable analysis outputs
└── data/, resources/, sqanti/, reference_files/   # gitignored; large reference / raw data
```

## Cell-type abbreviations

| Code | Cell type | In manuscript scope? |
|---|---|---|
| AT (or AT2) | Alveolar type 2 | ✓ |
| DD | LAE, submerged | ✓ |
| DD_ALI | LAE, air-liquid interface | ✗ |
| DO | SAEC, submerged (rare) | ✗ |
| DO_ALI | SAEC ALI (= grant Aim 1 prep) | ✗ |
| FB | Lung fibroblasts | ✓ |
| MV (or MVE) | Microvascular endothelial | ✓ |

**Label rename (April 2026):** `AT → AT2`, `DO → DO_ALI` in canonical labeling. Files pre-dating commit `1cd5328` use the older labels.

## Data file conventions

**Filename date-stamp:** `yyyy.m.d` (e.g., `2026.3.10`). When writing a new analysis output, use today's date.

**Per-cell-type mashr CSVs** (`shortread_dge/mashr/`, `isocall_dge/mashr/`, `isocall_dge/old/mashr/`):
- `ensembl_gene_id_version` (or `txid`), `logFC`, `adj.P.Val`, `nmd_responsive`, `hgnc_symbol` (and `gene_id`)
- `logFC` is the **mashr posterior mean**.
- `adj.P.Val` is the **limma-voom FDR** underlying mashr (NOT the mashr `lfsr`).
- `nmd_responsive == TRUE` means `mashr lfsr < 0.05 AND posterior mean > 0`.

**Shared mashr objects** (`mashr_{lfsr, posterior_means, ...}_2026.3.10.csv`):
- Wide format with one column per cell type (`Smg1i_in_DD`, `Smg1i_in_AT`, etc.)
- 4-CT scope (`DD, AT, FB, MV`). DO_ALI and DD_ALI are NOT in the shared 4-CT objects.

**Phenotype** (`pheno/`):
- `nmd_pheno_longreadbamids_YYYY.M.D.csv`
- Columns: `lib.size`, `norm.factors`, `id` (donor), `sample`, `treatment` (Smg1i / DMSO), `ct` (cell type), `sample_num`, `bam`

## Loading data (quick reference)

```r
library(data.table)
# Gene-level mashr (DO_ALI included)
gene <- fread("shortread_dge/mashr/nmd_mashr_dge_doali_2026.3.10.csv")
sig  <- gene[adj.P.Val < 0.05]
nmd  <- gene[nmd_responsive == TRUE]

# Isoform-level mashr (4-CT only)
iso <- fread("isocall_dge/mashr/nmd_mashr_die_dd_2026.3.10.csv")

# Shared mashr (wide; 4-CT)
lfsr <- fread("shortread_dge/mashr/mashr_lfsr_2026.3.10.csv")
```

```python
import pandas as pd
gene = pd.read_csv("shortread_dge/mashr/nmd_mashr_dge_doali_2026.3.10.csv")
sig  = gene[gene["adj.P.Val"] < 0.05]
nmd  = gene[gene["nmd_responsive"] == True]
```

## Reminders that bite

- **Baseline expression = DMSO samples only.** Never use limma's `AveExpr` or any pooled CPM as a baseline.
- **SQANTI/TD2 CDS calls are biased against PTC-containing ORFs.** For main-ORF PTC substrates, use reference-projected ORFs from the `NMD_orf_model_v5_4ct` repo. Exception: uORF-mediated NMD (e.g., ATF4) is fine with SQANTI CDS.
- **Every `git push` goes to BOTH Channing GitLab and PUBLIC GitHub.** No staging gate. Review before committing anything sensitive. (See ONBOARDING.md §4.)
- **`results/` is gitignored except for `results/isoform_transitions/Version_6.0/`.** Anything you save elsewhere under `results/` is invisible to git.
- **Manuscript prose changes return as Google Docs find/replace pairs, not commits to the markdown.** The Google Doc is the source of truth.

## Linked repos (the model and the package live elsewhere)

| Repo | What lives there |
|---|---|
| `peter4244/Isopair` | The Isopair R/Bioconductor package (gene-matched isoform pair analysis, splicing event detection, PTC attribution) |
| `peter4244/NMD_orf_model_v5_4ct` | Canonical source for the manuscript's deep-learning model — `03_train.py`, `config.yaml`, `model.py`, `METHODS.md`, DeepSHAP scripts |
| `peter4244/NMD_orf_model_v5` | Predecessor of v5_4ct; logic-identical training code, pre-4-CT-scope |
| `peter4244/Isoscope`, `peter4244/Isovar` | Per-gene isoform annotation/visualization; variant→isoform sQTL tooling |

When writing or auditing methods text, pull from the linked repos' METHODS.md / vignettes as the source of truth — not from anything in this repo's `results/`.

## Working in this repo

- **Branch from `main`** for feature work; `main` is the dual-push target.
- **Use the figure tooling in `figures/lib/`** (not `~/.claude/`) so figures stay reproducible for both Pete and Yul.
- **Save analysis outputs to `tmp/`** (regenerable, gitignored) or to `code/` (tracked) — never to `results/` outside the force-tracked Version_6.0 subtree.
- **For deep-learning model work, clone `NMD_orf_model_v5_4ct` separately.** This repo doesn't contain training code; a stale local rsync mirror exists at `results/isoform_transitions/NMD_orf_model_v5/` but is not version-controlled and should not be edited.

## Where to go next

- Project context, conventions, and collaboration rules → [`ONBOARDING.md`](ONBOARDING.md)
- Figure tooling and style guide → [`figures/README.md`](figures/README.md)
- Manuscript markdown → [`paper/NMD manuscript 2026.2.5.md`](paper/NMD%20manuscript%202026.2.5.md)
- Isopair analysis pipeline → [`results/isoform_transitions/Version_6.0/isopair_wrapper/`](results/isoform_transitions/Version_6.0/isopair_wrapper/)
