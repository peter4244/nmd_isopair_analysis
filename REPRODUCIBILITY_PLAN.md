# Reproducibility plan — NMD long-read lung study

**Current as of 2026-07-25.** Supersedes the 2026-07-24 draft; settled decisions are folded in
rather than left as open questions.

**Goal.** A reader with the Zenodo source-data record and the code repositories can reproduce
every number and figure in §1–§5 — on a machine with no Channing access, with no path edits,
and without retraining the model.

```
GEO GSE329233 (raw reads)
      │   documented, NOT on the reproduction path (the M2 gap lives here)
      ▼
★ ZENODO SOURCE-DATA RECORD  ←── reproduction starts here
      │
      ▼
nmd_lung_longread_2026  (§1–§5)  ·  Isopair  ·  Isocall_v1
      │
      ▼
every §1–§5 number, table and figure
```

---

## 1. Status

| Phase | State |
|---|---|
| 1. Scope + decisions | ✅ complete (§2) |
| 2. Source-data deposit | ✅ built, trimmed, validated (§3) — Pete uploading to Zenodo |
| 3. New repository | ✅ created, 395 files (§4) |
| 4. Regression baseline | ✅ captured and committed **before** any edit (§6.1) |
| 5. **The rewrite** | ⬜ **next — the subject of this plan** (§5) |
| 6. Regression check | ⬜ after the rewrite (§6.1) |
| 7. Clean-room test | ⬜ the gate before minting the code DOI (§6.2) |

---

## 2. Settled decisions

| # | Decision |
|---|---|
| 1 | **Tier 1 only.** Deposit the irreducible starting files; anything a shipped script regenerates stays out. |
| 2 | SQANTI `classification.txt` deposited in full (no slimmed copy). |
| 3 | Isopair `analysis_cache` **not** deposited — regeneration is tested instead. |
| 4 | Example-gene annotation included **if** a figure reads it (see §7 B). |
| 5 | **Tan et al. tables not redistributed** — the README cites the article DOI and names tables/sheets. |
| 6 | ORF-scan objects are 4-CT-valid (the scan is sequence-derived; only labels are cell-type-specific). |
| 7 | Source data gets **its own Zenodo DOI**, separate from the code record. |
| 8 | No sign-off required to deposit derived objects. |
| 9 | **One consolidated code repo** (`nmd_lung_longread_2026`) for §1–§5; Isopair and Isocall stay separate. |
| 10 | **Option C — the code is rewritten to be 4-cell-type native.** No shipped script loads other conditions and drops them. |
| 11 | **The repo holds the code that produced the paper's results, and not more.** History and internal process docs are excluded; correctness rests on verification, not provenance. |

---

## 3. Source-data deposit — final contents

`~/claude_projects/nmd_deposit_2026/` — **10 files, 4.0 GB**

| Path | Contents |
|---|---|
| `nmd_lungcells_counts_4ct.csv` | Long-read isoform counts, 614,992 × 26 |
| `salmon_gene_counts_4ct.csv` | Short-read gene counts, 46,571 × 26 |
| `sqanti/` | classification · corrected CDS GFF3 · filtered GTF · corrected FASTA |
| `annotation/` | Ensembl v115 gene + transcript maps |
| `model/` | trained model weights |
| `MANIFEST.sha256`, `README.md`, `.zenodo.json` | checksums, contents, metadata |

**Validated:** the count matrix, classification and FASTA carry an **identical isoform ID set**
(614,992) — not merely equal counts.

**Restricted to observed features.** Features with zero counts in every sample were removed
(isoforms 645,272 → 614,992; genes 78,899 → 46,571) after verifying it is result-neutral: the
`filterByExpr` sets are identical and percent-output-lost is unchanged to four decimals.

**Excluded because shipped code regenerates them:** DGELists, mashr results, the Isopair cache,
ORF-scan objects, the model HDF5, and all interpretability exports.

---

## 4. The repository

`nmd_lung_longread_2026` — 395 files.

```
analysis/{upstream,isopair,predictor_comparison}   §1–§3 · §4 · §5 benchmark
model/                                             §5 training + interpretability
figures/{lib,multipanel,supplemental}              Figures 1–5, SF1–43
verification/                                      verifier suite + regression baseline
metadata/, docs/, config/
```

Two defects fixed by construction: the §4 pipeline no longer sits under `results/` (whose bare
gitignore rule caused four near-misses), and the model is in-tree, removing the absolute
cross-repo `source()` that made Figure 5 panel A unrunnable off one laptop.

**Producer recovery.** `05s_orfik_scan.R`, `05t_ref_cds_features.R` and
`05u_paralog_annotation.R` were restored from history after a 2026-07-11 cleanup deleted them as
"legacy report feeders". They are the sole producers of the model's ORF feature inputs; without
them the model chain had no producer at all. A sweep of all 11 scripts named in `PIPELINE.md`
found no further gaps — `05l`/`05v` are correctly absent (deprecated-report feeders only) and
`04` was a stale reference, superseded by the present `03b_rebuild_cache.R`.

---

## 5. The rewrite — execution plan

Two changes applied in **one pass**, because they touch the same files and separating them would
force two full re-verifications.

### 5.1 Path portability

**Measured:** 17 of 59 files under `analysis/upstream/` carry `/udd/reyle/...`; `05s` carries an
absolute FASTA path; the `verify_section*.R` scripts carry `~/claude_projects/nmd`.

1. `config/paths.yml` — one root per input class (`DEPOSIT`, `CACHE`, `OUT`), defaulting to
   `./data_deposit/` so an unpacked Zenodo record works with no edits.
2. `R/load_config.R` plus a Python equivalent to resolve them.
3. Rewrite every hardcoded path to resolve through config.
4. **CI linter** failing on `/udd/`, `/proj/`, `/home/`, `~/claude_projects` in tracked code —
   without it, absolute paths creep back.

### 5.2 4-cell-type native

1. Delete the sample-dropping blocks (`Isoform_Level_Quantification.Rmd`,
   `NMD_shortread_dge_fullmodel_2026.5.5.Rmd`). The deposited matrices contain only the four
   cell types, so these are already no-ops.
2. Keep `CELL_TYPES` as a **scope declaration** (it orders factors and documents scope); remove
   the *drop* framing and references to other conditions across ~32 files.
3. Normalise cell-type naming to `AT2`/`LAE` throughout. The deposit already uses canonical
   names, removing the alias-patching that `resolve_mash_col()`, `resolve_sr()` and
   `CELLTYPE_MAP` currently need.
4. Update READMEs and in-repo docs to state four-cell-type scope.

### 5.3 Pre-declared expected differences

Recorded in advance so they are predictions, not post-hoc excuses:

- `verify_section3_p1.txt` — **exactly two lines** change (`DGEList dimensions:` and
  `CPM matrix:`), because the deposit is trimmed to observed isoforms; values become 614992.
- **All nine numeric-result lines must be identical.** Anything else differing is a defect.

---

## 6. Verification gates

### 6.1 Regression check — gate on the rewrite

`verification/capture_verification_baseline.sh --check`. Baseline captured 2026-07-24 before any
edit: 9 scripts, 358 lines, covering ~47 verified claims across §1–§3. Requirement: every number
identical apart from §5.3's two declared lines.

⚠️ The script's baseline path must be updated for the new layout (`reproducibility/` →
`verification/`) before `--check` will run there.

### 6.2 Clean-room test — gate on the DOI

Fresh clone into an empty directory on a machine with no Channing access; download the Zenodo
record; unpack to `./data_deposit/`; run `REPRODUCTION.md` end to end. Compare **numbers exactly**
and **figures by PNG byte-compare** (PDFs embed a creation timestamp and never match).

**Success = a fresh clone reproduces §1–§5 with zero path edits.** Until it passes, the honest
claim is "the code is archived", not "the paper is reproducible". Under decision 11 there is no
history to fall back on, so this gate is not optional.

---

## 7. Open items

| # | Item |
|---|---|
| A | **`05_final_report_mashr.Rmd` — recommend excluding.** Deprecated, cannot run from a clean checkout, its feeders (`05l`/`05v`) are intentionally absent, and no shipped figure reads its outputs. Pete's call. |
| B | **Example-gene annotation** (`srsf.*`, `sub.isoforms.gtf`) is not on this laptop. Audit whether any figure reads it; if none does, close decision 4 as "not needed". |
| C | **`REPRODUCTION.md` not yet written** — the runbook mapping claim/figure → script → expected output, stating plainly what is *not* reproducible (M2, model retraining, the Tan download). |
| D | **Zenodo metadata before publishing:** expand author given-names, add ORCIDs, add related-identifier links to the code DOI once minted. Immutable at mint. |
| E | **ORFik version risk.** 1.30.2 is pinned, but that is what is installed now, not verifiably what produced the original scan. If a regenerated scan drifts, §5 figures will not match — the clean-room test is the detector. |
| F | **GEO GSE329233** must be reduced to the four cell types' reads, after verification (Pete, 2026-07-25). |
