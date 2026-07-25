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

**Revised after Pete's correction (2026-07-24):** *"the DGELists would be produced by our code."*
Correct — and verified, with a consequence that goes further:

- `Isoform_Level_Quantification.Rmd` reads the SQANTI count matrix + GTF + gmap and writes
  `dge_isoform_longread_2026.3.3.rds` (`:275`), the filtered/TMM version (`:302`), **and**
  `mashr_isoform_model_*.rds` (`:565`).
- `NMD_shortread_dge_fullmodel_2026.5.5.Rmd` reads the salmon gene counts + gmap and writes
  `dge_shortread_gene_2026.3.2.rds` (`:212`), the filtered version (`:227`), **and**
  `mashr_model_*.rds` (`:628`).

So **both the DGELists and the mashr results are intermediates with tracked producers.**
Depositing them as "starting files" would violate ship-producers-not-products.

**Determinism checked**, since regeneration only helps if it reproduces: mashr is seeded
(`set.seed(42)` before `cov_pca`; `cov_pca` + `cov_canonical` only, no stochastic extreme
deconvolution), and `calcNormFactors(method="TMM")` is deterministic. Regeneration should
reproduce exactly, modulo package versions.

This yields a **two-tier deposit**.

### 1.1 Tier 1 — irreducible starting files (no in-repo producer)

| # | File(s) | Size | Why irreducible |
|---|---|---|---|
| **A** | **Isoform count matrix** — `nmd_lungcells_filtered.count_matrix.txt` | 65 MB | SQANTI/isocall output; regenerating needs raw reads + heavy compute |
| **B** | **SQANTI `classification.txt`** | 368 MB (slim ≈5 MB) | Same pipeline. **Distinct from the GFF** — supplies `structural_category` (1.4), `coding` (protein-coding definition in `productive_response`), `associated_gene` |
| **C** | **SQANTI structures** — `_corrected.cds.gff3` (1.1 GB), `_filtered.gtf` (1.1 GB), `_corrected.fasta` (1.6 GB) | **3.8 GB** | §4 isoform structures + ORF/PTC projection. **See scoping question §7.2** — the FASTA may be droppable if the model is not re-run |
| **D** | Short-read gene counts — `salmon.merged.gene_counts_length_scaled` | 17 MB | **nf-core/rnaseq (M2) output, and M2 is the accepted gap** — so this is genuinely a starting point, not an intermediate |
| **E** | Gene/transcript annotation maps — `gmap_ENSGv115`, `gmap_txlevel_ENSGv115` | 6.8 MB | Built outside this project (`/udd/repjc/RESEARCH/REUSABLE_CODE/`), so no producer here |
| **F** | Sample/phenotype metadata — `pheno/*` | <1 MB | ✅ already tracked; keep in GitHub |
| **G** | RBP rosters | 135 KB | Third-party. ⚠️ on disk but silently gitignored by a bare `data/` rule (`.gitignore:31`) |
| **H** | **NMD model** — `best_model_atg500_stop500.pt` + **the interpretation export TSVs** | ~50–100 MB | Retraining is out of scope by your rule. Weights alone do **not** reproduce Fig 5 / SF37–43 — those panels read the export TSVs, and **we stripped them from the model repo in the v2.0.0 simplification**, so they currently exist in neither repo |

**Tier 1 ≈ 4.3 GB** (≈0.5 GB if the three large SQANTI structure files are scoped out).

### 1.2 Tier 2 — regenerable checkpoints (deposit as baselines, not requirements)

Everything here has a tracked producer, so it is **not required**. Depositing it anyway is
cheap and buys two things: a reader can start midway instead of re-running the whole chain, and
— more valuably — it becomes the **regression baseline** for the clean-room test (§5), letting
us prove regeneration matches rather than assuming it.

| Item | Size | Producer |
|---|---|---|
| The 4 DGELists (+ `dge_gene_unfiltered`) | 48 MB | `Isoform_Level_Quantification.Rmd`, `NMD_shortread_dge_fullmodel_2026.5.5.Rmd` |
| Complete mashr results — per-CT CSVs, shared lfsr/posterior-means/sharing objects, fitted `mashr_*_model.rds` | 166 MB | same two scripts + `ct_de.Rmd` |
| Isopair `data_mashr/analysis_cache` (83 objects) | 149 MB | `scripts/core/02–05` → `isopair_wrapper/01–06` |

**Tier 2 ≈ 360 MB.** Label it plainly in the deposit as *derived checkpoints, regenerable from
Tier 1*, so nobody mistakes it for primary data.

### 1.3 Remaining scoping decisions

| # | Item | Size | The question |
|---|---|---|---|
| **K** | Isopair `analysis_cache` | 149 MB | Now Tier 2 (regenerable). Question is only whether we **test** that regeneration works, or ship it and skip the test. Recommend: test it. |
| **L** | `nmd_orf_data.h5` | **4 GB** | Needed only to *re-run* inference/SHAP. Under "no retraining", J's TSVs suffice. Recommend **exclude**, and document how to rebuild it via `data_prep.py`. |
| **M** | Tan et al. supplementary tables | small | Third-party; redistribution is questionable and they were deliberately gitignored. Recommend **document the download**, do not redistribute — 2.18 then reproduces from a documented fetch. |
| **N** | `tan_tx_mashr_model.rds` | small | Our fitted object over Tan's data — depositing it makes 2.18 reproducible even if their download moves. Recommend **include**. |
| **O** | Example-gene annotation (`srsf.*`, `sub.isoforms.gtf`, …) | small | Only if figure panels read them. **Audit, then include if used.** |

### 1.4 Deliberately excluded (documented, not deposited)

- **Reference genome / GENCODE / Ensembl v115** — large, public, versioned. Pin versions in Methods.
- **Raw reads** — GEO GSE329233.
- **Intermediate figures / rendered reports** — regenerable.
- **The nf-core/rnaseq launcher (M2)** — your standing accepted gap. State it plainly: the deposit starts *downstream* of M2, so M2 is not on the reproduction path.

---

## 2. Where each piece goes

| Destination | Contents | Rationale |
|---|---|---|
| **GitHub (code repos)** | all code, `pheno/` (F), RBP rosters (G), the RUNBOOK, config template | Small, diffable, already the deposit philosophy |
| **New Zenodo *data* record** | **Tier 1** (A–E, H) and, clearly labelled as regenerable, **Tier 2** | **Required**: the 368 MB classification file alone exceeds GitHub's 100 MB per-file hard limit, and Tier 1 runs to ~4.3 GB. Zenodo allows 50 GB/record |

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

1. **Tier 2** — deposit the regenerable checkpoints as baselines (recommended, +360 MB), or omit them and rely purely on regeneration?
2. **SQANTI structures (C, 3.8 GB)** — is the 1.6 GB `corrected.fasta` needed if the model is not re-run, or do §4/ORF analyses require sequence? This one decision moves Tier 1 between ~0.5 GB and ~4.3 GB.
3. **SQANTI classification** — full 368 MB, a slimmed column subset (~5 MB), or both?
4. **Model H5** (4 GB) — exclude and document, as recommended?
5. **Tan tables** — document-the-download, or check whether the licence permits redistribution?
6. **Does the data record get its own DOI**, or become a new version of `nmd_isopair_analysis`?
   (Separate record is cleaner: data and code version independently.)
7. **Yul's sign-off** — less pressing now that the DGELists are Tier 2 rather than required inputs, but still worth asking before depositing her derived objects.
