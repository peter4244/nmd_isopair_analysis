# Reproducibility plan — code repos + a defined starting-file deposit

**Drafted 2026-07-24.** Proposal, not executed. Nothing moves without sign-off.

**Goal.** A reader with (a) the GitHub/Zenodo code repos and (b) one deposited starting-file
set can reproduce every number and figure in §1–§5 without access to Channing, without editing
paths, and without retraining the model.

**Design.** Code repos stay code-only (unchanged deposit philosophy). A **separate Zenodo data
record** holds the starting files. The chain from raw reads to those starting files is
*documented but not required* — that is where the accepted M2 gap lives.

```
GEO GSE329233 (raw reads)
      │   [documented, NOT required for reproduction — M2 gap lives here]
      ▼
★ STARTING-FILE DEPOSIT (new Zenodo data record)  ←── reproduction starts HERE
      │
      ▼
code repos (nmd_isopair_analysis · Isopair · Isocall_v1 · NMD_orf_model_v5_4ct)
      │
      ▼
every §1–§5 number, table and figure
```

---

## 1. The starting-file set

Pete's proposed four are in **bold**; everything else is an addition this audit identified.

### 1.1 Required

| # | File(s) | Size | Why required | Status today |
|---|---|---|---|---|
| **A** | **Complete mashr results** — per-CT `nmd_mashr_dge_{at,dd,fb,mv}` + `nmd_mashr_die_{…}`, shared `mashr_{,isoform_}{lfsr,posterior_means,sharing_2fold,sharing_sign,null_correlation}`, plus the fitted `mashr_{,isoform_}model_*.rds` | **166 MB** | Backs nearly every §1–§3 claim | ⚠️ **gitignored** (`.gitignore:65,69`) |
| **B** | **Raw isoform count matrix** | see C | §1 landscape, §3 output-lost | ⚠️ not tracked |
| **C** | **The four DGELists** — `dge_isoform_longread{,_filtered}_2026.3.3`, `dge_shortread_gene{,_filtered}_2026.3.2`, `dge_gene_unfiltered_2026.1.2` | 48 MB | **B alone is insufficient**: the DGEList also carries `$samples` (donor, treatment, ct, lib.size, **norm.factors**) and `$genes` (txid→gene_id→symbol). Shipping only a count matrix forces the reader to recompute norm factors and rebuild the isoform→gene map, which will not reproduce exactly | ⚠️ not tracked |
| **D** | **SQANTI GFF** — `nmd_lungcells_corrected.cds.gff3`, `_corrected.fasta`, `_filtered.gtf`, `nmd_isocall.isoforms.gtf` | ~GB | §4 structures, ORF projection | ⚠️ not tracked; "no identified home" |
| **E** | **SQANTI `classification.txt`** | **368 MB** (slim ≈5 MB) | **Distinct from the GFF and not covered by it.** Supplies `structural_category` (claim 1.4), `coding` (protein-coding gene definition in `productive_response`), `associated_gene`. Read by at least 3 pipelines | ⚠️ not tracked |
| **F** | **Short-read gene counts** — `salmon.merged.gene_counts_length_scaled.rds` | 17 MB | Fig 1A SR↔LR correlation (1.6–1.8), §2 gene-level, and the SR induction signal in `productive_response` | ⚠️ not tracked |
| **G** | Sample/phenotype metadata — `pheno/*.csv/.tsv/.fofn` | <1 MB | Maps count columns → donor/treatment/cell type | ✅ **already tracked** |
| **H** | Gene/transcript annotation maps — `gmap_ENSGv115`, `gmap_txlevel_ENSGv115` | 6.8 MB | Symbol mapping; pinning beats "regenerate from Ensembl v115" | ⚠️ not tracked |
| **I** | RBP rosters — `gerstberger_2014_rbp_census.csv`, `encode_rbp_roster_vannostrand2020.csv` | 135 KB | §2 RBP enrichment (M8a) | ⚠️ **on disk but gitignored** by a bare `data/` rule (`.gitignore:31`) — see §4 |
| **J** | **NMD model** — `best_model_atg500_stop500.pt` (453 KB) **+ the interpretation export TSVs** (predictions, attention, DeepSHAP, motif, GC, subgroup, kernel-SHAP) | ~50–100 MB | Weights alone do **not** reproduce Fig 5 / SF37–43 — those panels read the export TSVs, and **we stripped them from the model repo in the v2.0.0 simplification.** Under a no-retraining rule they must live in the data deposit | ⚠️ **stripped from the model repo** |

**Subtotal ≈ 700 MB** with the full classification file, ≈340 MB with a slimmed one.

### 1.2 Decisions needed

| # | Item | Size | The question |
|---|---|---|---|
| **K** | Isopair `data_mashr/analysis_cache` (83 objects) | **149 MB** | **Is §4 reproduced from the deposit, or regenerated?** The Version_6.0 *code* is tracked (54 files) and in principle rebuilds the cache from D+E+C+A via `scripts/core/02–05` → `isopair_wrapper/01–06`. Cheaper to deposit the cache; more honest to prove regeneration. Recommend: **test regeneration; deposit the cache as a fallback.** |
| **L** | `nmd_orf_data.h5` | **4 GB** | Needed only to *re-run* inference/SHAP. Under "no retraining", J's TSVs suffice. Recommend **exclude**, and document how to rebuild it via `data_prep.py`. |
| **M** | Tan et al. supplementary tables | small | Third-party; redistribution is questionable and they were deliberately gitignored. Recommend **document the download**, do not redistribute — 2.18 then reproduces from a documented fetch. |
| **N** | `tan_tx_mashr_model.rds` | small | Our fitted object over Tan's data — depositing it makes 2.18 reproducible even if their download moves. Recommend **include**. |
| **O** | Example-gene annotation (`srsf.*`, `sub.isoforms.gtf`, …) | small | Only if figure panels read them. **Audit, then include if used.** |

### 1.3 Deliberately excluded (documented, not deposited)

- **Reference genome / GENCODE / Ensembl v115** — large, public, versioned. Pin versions in Methods.
- **Raw reads** — GEO GSE329233.
- **Intermediate figures / rendered reports** — regenerable.
- **The nf-core/rnaseq launcher (M2)** — your standing accepted gap. State it plainly: the deposit starts *downstream* of M2, so M2 is not on the reproduction path.

---

## 2. Where each piece goes

| Destination | Contents | Rationale |
|---|---|---|
| **GitHub (code repos)** | all code, `pheno/`, RBP rosters (I), the RUNBOOK, config template | Small, diffable, already the deposit philosophy |
| **New Zenodo *data* record** | A–F, H, J, and whichever of K/N/O are approved | **Required**: the 368 MB classification file exceeds GitHub's 100 MB per-file hard limit; the full set (~700 MB) is unreasonable in git |

Cite the data DOI in Data Availability alongside GEO, and reference it from each code repo's README.

---

## 3. De-hardcoding the code (the portability work)

**Measured:** 17 of 59 tracked files in `code/upstream/` hardcode `/udd/reyle/...`. My own
`verify_*.R` scripts hardcode `~/claude_projects/nmd` and are equally guilty.

**Approach — one config, no path literals:**

1. Add `config/paths.yml` (or `.Renviron` + `here::here()`) with one root per input class:
   `MASHR_DIR`, `DGE_DIR`, `SQANTI_DIR`, `MODEL_DIR`, `OUT_DIR`.
2. Add `R/load_config.R` (and a Python equivalent) resolving those roots, defaulting to
   `./data_deposit/` so an unzipped Zenodo record "just works".
3. Rewrite the 17 files (+ my verify scripts) to read from config. Mechanical.
4. **Add a linter to CI** that fails on `/udd/`, `/proj/`, `/home/`, `~/claude_projects` in
   tracked code — this is the "rules in code, not docs" pattern; without it the paths creep back.

**Also fix the `.gitignore` hazards.** The bare `data/` rule (`:31`) silently excludes the RBP
rosters, and bare `shortread_dge/mashr/` (`:69`) and `nmd_fig_data/` (`:103`) hide required
inputs. This is the fourth instance of the bare-directory-rule hazard in this project. Anchor
them and add explicit `!` exemptions.

---

## 4. RUNBOOK (new, tracked)

A single `REPRODUCTION.md` giving, in order: environment (R 4.5.2 / Bioc 3.22 / msigdbr 26.1.0
/ **seed 42**; Python 3.14.4), how to fetch and unpack the data record, then a numbered table of
*claim/figure → script → expected output*. It must state plainly what is **not** reproducible
from the deposit: M2, model retraining, and the Tan download.

---

## 5. The verification that actually proves it

Everything above is unproven until this runs. **Clean-room test:**

1. Fresh clone of all four repos into an empty directory on a machine with no Channing access.
2. Download the Zenodo data record; unpack to `./data_deposit/`.
3. Run the RUNBOOK top to bottom.
4. Compare against the manuscript: **numbers exact; figures by PNG byte-compare** (your
   established N-2 finding — PDFs embed `/CreationDate` and never match).

**Success = a fresh clone reproduces §1–§5 with zero path edits.** Until that passes, the
correct claim is "the code is archived", not "the paper is reproducible".

---

## 6. Suggested phasing

| Phase | Work | Effort | Value |
|---|---|---|---|
| **1** | Confirm the file list (§1.1) and settle K/L/M/N/O | — | Unblocks everything |
| **2** | Fix `.gitignore`; track I; assemble the deposit; mint the data DOI | ~half day | Makes inputs obtainable at all |
| **3** | Config system + rewrite the 17 files + CI path linter | ~1 day | Removes the biggest barrier |
| **4** | Write `REPRODUCTION.md` | ~half day | Makes the chain followable |
| **5** | Clean-room test; fix what breaks | ~1–2 days | **The only step that proves the claim** |

Phases 2–4 are worth doing before submission regardless. Phase 5 is what converts
"comprehensive repos" into "verified reproducible".

---

## 7. Open questions for Pete

1. **§4 scope** — deposit the Isopair cache (K, 149 MB), or require regeneration and test it?
2. **SQANTI classification** — ship the full 368 MB, a slimmed column subset (~5 MB), or both?
3. **Model H5** (L, 4 GB) — exclude and document, as recommended?
4. **Tan tables** (M) — document-the-download, or check whether the licence permits redistribution?
5. **Does the data record get its own DOI**, or become a new version of `nmd_isopair_analysis`?
   (Separate record is cleaner: data and code version independently.)
6. **Is Yul's sign-off needed** on depositing her derived intermediates (the DGELists are her outputs)?
