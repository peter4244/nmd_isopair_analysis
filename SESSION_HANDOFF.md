# Session handoff — NMD repo consolidation

Written 2026-07-22 at the end of a long session, for a fresh-context successor.
**Read this first, then the four docs listed under "Orientation".**

---

## ⏭ NEXT (updated 2026-07-24 — read this before the older sections below)

**Deposit workstream is DONE.** Model repo simplified → `v2.0.0` released → all four Zenodo
DOIs minted (see DOI table) → Code Availability section drafted at `paper/code_availability.md`.

**Active next task: §3 manuscript verification — 26 unverified claims (parked).** Run them
through the 5-step process one at a time. Also open: §1 (8/15) and §2 (15/26 + 3 robust)
stragglers; Claim 1.11 (AT2/LAE Spearman range) pending Yul's rerun.

**Small closeout items (Yul-dependent):** add full author metadata + ORCIDs to the 4 Zenodo
records; Yul's sign-off on the deposits; paste Code Availability into the Google Doc (+ make
"Figure 5" a live cross-ref).

---

## Where things stand

The consolidation is ~80% done. **Phases 0, 1, 3, 4, 5.2, 6.1, 6.5 complete.**
Remaining: Phase 7 (citable snapshot) and Phase 8 (submission compliance).

| | |
|---|---|
| Active tracked code | **301 → 157** |
| Archived (reversible, `superseded/<original/path>`) | 158 |
| Rollback branch | `archive/pre-consolidation` (pushed to GitHub) |
| Working tree | clean; GitHub `main` current |
| Channing (`changit`) | **~70 commits behind — needs a VPN push from Pete** |

### Orientation — read in this order
1. `REPO_CONSOLIDATION_PLAN.md` — the plan (v3.1), all rulings folded in
2. `KEEP_LIST.md` §E–I — keep-set, trace findings, the review's catch, producer audit
3. `paper/results_to_code_map.md` — per-claim verification status
4. Memory: `project_nmd_repo_merge_citable`, `feedback_citable_repo_code_and_starting_data_only`

---

## ✅ DECIDED (Pete, 2026-07-22): external repo packaging = **Option B**

**Each external repo gets its own Zenodo DOI**, cited from Code Availability. Submodules
are ruled out; vendoring is used only for the single cross-repo file (see "Recommended"
below). **This decision is made — do not re-litigate it.** What remains is execution:

1. Add a **LICENSE** to `NMD_orf_model_v5_4ct` and `Isocall_v1` (both currently unlicensed
   = all-rights-reserved; Zenodo/journals will object). `Isopair` already has MIT.
2. **Tag a release** in each of the three (none has any tag today, so nothing is pinnable).
3. Mint a Zenodo DOI per repo from that tag.
4. **Vendor** `make_architecture_figure.R` into this repo (or document a clone-as-siblings
   step) to kill the cross-repo `source()` in `figure5_panelA_architecture.R:33`.
5. Cite all three DOIs in the manuscript's Code Availability section — **which does not yet
   exist and must be written.**

Supporting detail and the options considered are kept below for the record.

### ⚠️ The repo split is by KIND, not by manuscript section

A natural but **wrong** assumption is "§1–4 in this repo, §5 in the model repo." Measured
2026-07-22, this repo holds **25 §5 files**:

- `figures/multipanel/figure5_dl_model/` (10) — all Figure 5 panels
- `figures/SupplementalFigures/SF37–SF43` (9) — SHAP-across-windows, stop-codon usage,
  attention distribution, PTC-subclass branch SHAP + performance, GC content, model comparison
- `code/nmd_predictor_comparison/` (6) — the NMDetective-B / NMDEP benchmark behind SF43

`NMD_orf_model_v5_4ct` holds the model itself: `03_train.py`, `model.py`, `config.yaml`,
the DeepSHAP/KernelSHAP scripts, its own `METHODS.md`.

| Layer | Where |
|---|---|
| §1–3 analysis | this repo (`code/upstream/`, Yul) |
| §4 Isopair pipeline + figures | this repo |
| §5 **model training + interpretability computation** | external model repo |
| §5 **figure rendering + predictor benchmark** | **this repo** |

**The model repo computes; this repo visualizes and benchmarks.** Figure 5 panels read TSV
exports from the model repo, and `figure5_panelA` additionally `source()`s a script from it.

**Consequence for Code Availability:** do NOT write "sections 1–4 here, section 5 there" —
reproducing **Figure 5 requires both repos**. Phrase by role, e.g.:

> Analysis and figure-generation code: [this repo DOI]. Deep-learning model training and
> interpretability: [`NMD_orf_model_v5_4ct` DOI]. Long-read processing pipeline:
> [`Isocall_v1` DOI]. Isoform-pair analysis package: [`Isopair` DOI].

---

## Background — external repo packaging (options as considered)

Three external repos are needed to reproduce the paper. **Facts (checked 2026-07-22):**

| Repo | Public | License | Tags | Size | Role |
|---|---|---|---|---|---|
| `peter4244/Isopair` | ✅ | MIT | **0** | 967 KB | §4 splice-event/PTC package |
| `peter4244/NMD_orf_model_v5_4ct` | ✅ | **none** | **0** | 1.2 MB | §5 DL model |
| `peter4244/Isocall_v1` | ✅ | **none** | **0** | 23 KB | M3 long-read pipeline |

**Two blockers:** two repos are unlicensed (legally all-rights-reserved — Zenodo/journals
will object), and **none has a tagged release**, so nothing is currently pinnable.

### Options
- **A. Cite URL only.** Zero work; `main` moves, so a future reader gets different code.
  Weakest reproducibility.
- **B. Own Zenodo DOI each.** ⭐ *Recommended.* Standard for research software with
  independent identity. Immutable + versioned + citable. Needs: licenses on two repos,
  one tagged release each.
- **C. Git submodule.** ❌ **Rule out.** GitHub release tarballs and Zenodo archives do
  **not** include submodule contents — the deposit would ship three empty directories.
- **D. Vendor into the snapshot.** Self-contained and cheap (~2.2 MB), but duplicates a
  package with its own life and muddies attribution.

### Recommended: B + one targeted vendor
DOI all three. **Exception:** `figures/multipanel/figure5_dl_model/figure5_panelA_architecture.R:33`
does `source(UPSTREAM)` where `UPSTREAM = ../NMD_orf_model_v5_4ct/make_architecture_figure.R`
— a relative path *outside* the repo that breaks under A, B **and** C. Vendor that single
file (or document a "clone as siblings" setup step).

---

## DOIs — collected as Pete mints them (Code Availability)

| Repo | DOI | status |
|---|---|---|
| `Isocall_v1` | **10.5281/zenodo.21536486** | ✅ minted 2026-07-24 |
| `Isopair` | **10.5281/zenodo.21536495** | ✅ minted 2026-07-24 |
| `NMD_orf_model_v5_4ct` v1.0.0 | **10.5281/zenodo.21536502** | history (unsimplified) |
| `NMD_orf_model_v5_4ct` **v2.0.0** | **10.5281/zenodo.21539601** | ✅ minted 2026-07-24 (code-only simplified; **cite this**) |
| `NMD_orf_model_v5_4ct` concept | **10.5281/zenodo.21536501** | always-latest |
| `nmd_isopair_analysis` v1.0.0 | **10.5281/zenodo.21539735** (concept **21539734**) | ✅ minted 2026-07-24 |

**Code Availability section drafted** → `paper/code_availability.md` (concept DOIs, all four repos + GEO GSE329233, verified via Zenodo API). Paste into the Google Doc; make "Figure 5" a live cross-ref. All four Zenodo records still need full author metadata added (with Yul).

⚠️ v2 Zenodo record author metadata is currently minimal (repo owner only — no `.zenodo.json`/
`CITATION.cff` at mint). Edit the Zenodo record to add the full author list + ORCIDs with Yul;
the DOI string is unaffected.

⚠️ **The model DOI (…502) was minted BEFORE the model-repo simplification**, so `v1.0.0`
archives the current *unsimplified* state — including the **broken `nmd_orf_data.h5`
symlink**. Two ways forward, Pete's call:
- **(a) Accept …502 as the citable version.** The code is all present; simplification
  becomes optional. But the archived tarball contains a dangling symlink and the one-off
  scripts. Also fix the broken h5 before relying on it for reproduction.
- **(b) Simplify, then mint a `v2.0.0`** and cite that version DOI instead (…502 stays as
  v1 history). Recommended if we want the deposited model artifact clean + runnable.
- Either way, the **concept DOI** (Zenodo's "Cite all versions") always resolves to the
  latest — citing the concept DOI in the paper future-proofs a later v2.

*Confirm version-vs-concept DOI for each; cite the **version** DOI (or the concept DOI if
you want it to track a future v2). Draft Code Availability text goes here once all four exist.*

## ▶ ACTIVE NEXT TASK (decided 2026-07-24): simplify the model repo → v2

Pete chose **option (b)**: simplify `NMD_orf_model_v5_4ct`, then publish a **`v2.0.0`**
release and cite that. `v1.0.0` = DOI `…502` stays as history.

### ✅ EXECUTED 2026-07-24 on branch `simplify-v2` (local, UNPUSHED — awaiting Pete's review then tag)

Pete signed off (4-question gate) on: commit the uORF-attention analysis into v2; strict
code-only strip of ALL `results_4ct/` data; archive old-v5 slurm variants; FF + branch.

Reconciled dirty state: fast-forwarded `master` 75599c2→755244f (= `v1.0.0`, the "Create
LICENSE" commit pulled via `git fetch --tags`); branched `simplify-v2` from there. Five commits:
1. `5303f37` add uORF-attention analysis (4 code files) + METHODS section.
2. `4f22eac` SF37 report edit (drop chi-sq subtitle).
3. `299b64e` **code-only strip** — deleted 8 dangling symlinks, untracked (kept on disk) 8
   real product TSVs, removed dead `nmd_orf_data.h5`/`selected_orfs.tsv` symlinks; `.gitignore`
   now blanket-ignores `results_4ct/`.
4. `df836b9` archived 14 files to `superseded/` (5-step QC logs + 7 verified-superseded slurm
   wrappers) + `superseded/README.md`.
5. `35924ad` README: corrected input provenance (`export_rds.R` produces ref_cds/td2 from the
   Isopair cache — the "upstream 05t_*.R" claim was stale) + documented code-only regeneration.

Tracked 88→77 (62 active + 15 in `superseded/`); **zero data files, zero dangling symlinks**.
Adversarial pre-prune review (subagent stalled at 600s watchdog → completed inline): split is
SAFE — every report/exporter input has an in-repo KEEP producer; **catch:** `BUGFIX_STOP_CODON`
must stay KEEP (cited by METHODS + report + patch script). **The whole symlink layer (not just
h5) was broken** — 10 dangling symlinks, all products of `export_rds.R`/`data_prep.py`; step 4
(h5) + step 5 (code-only) collapsed into the one strip commit.

**Render-verify (2026-07-24, Explorer):** rendered `orf_model_report_v5.Rmd` on the live
Explorer copy (master@75599c2, all inputs present) — **171/171 chunks, all 42 figures, 0
missing, no errors** ⇒ KEEP set complete, nothing stripped was load-bearing. Byte-compare
vs the Apr-30 baseline showed 42/42 diffs but that is **graphics-env drift, not content**
(inputs unchanged since April ⇒ content is deducibly identical); pixel-compare abandoned as
confounded. Pete satisfied; no `export_rds.R` smoke test run.

**✅ PUSHED (Mac, 2026-07-24):** `origin/master` = **35924ad** (the 5 simplification commits);
**`v2.0.0` tag pushed** (points at 35924ad). A draft `CITATION.cff` was briefly committed then
**dropped** at Pete's request (author metadata needs Yul; Zenodo captures its own metadata at
mint and it's editable post-publish). So v2.0.0 is the clean code-only simplification, no CFF.

**REMAINING for Pete:**
1. **Mint the v2 DOI:** create a **GitHub Release from tag `v2.0.0`** (repo → Releases → Draft
   new release → tag v2.0.0 → Publish). v1.0.0 already has a GitHub Release + Zenodo …502, so
   the webhook auto-mints a new *version* DOI under the same concept. ⚠️ No `.zenodo.json`/
   `CITATION.cff` ⇒ Zenodo captures **minimal metadata (repo owner only)** — grab the DOI now,
   then **edit the Zenodo record's authors/ORCIDs later with Yul** (DOI stays fixed; same as v1).
2. Cite the v2 **version** DOI (or the **concept** DOI for always-latest) in **Code Availability**
   — section still does not exist in the manuscript; write it (Phase 7 blocker).
3. Yul's sign-off on the deposit (Phase 7 prereq; can follow the mint since metadata is editable).

Explorer working copy (separate, master@75599c2 + 2 uncommitted slurm edits) will need a
`git fetch` to sync — not the deposit source, low priority.

**⚠️ Original note (still applies as history):** start fresh-window work here if resuming.

**Before touching anything in the model repo, reconcile its dirty state:**
- It has **4 uncommitted changes** — inspect and commit/stash first.
- **No local tag** for v1.0.0 (the release was made via GitHub web; `git fetch --tags`
  to pull it, so the prune branches off the right point).
- Remote: `github.com/peter4244/NMD_orf_model_v5_4ct` (single remote, not dual-push).

**Method — identical discipline to this repo's prune (which worked; see §§E–I and the
pre-prune review that caught a real break):**
1. Seed from manuscript deliverables: `03_train.py`, `model.py`, `utils.py`, `config.yaml`,
   the interpretability exporters (`04_interpret_attention.py`, `06`/`08_export_deepshap`,
   `11_kernel_shap_branches.py`, the `09_*` exporters that feed Fig 5 + SF37–43), METHODS.md.
2. Build a keep-set; get Pete's ruling on ambiguous files; **run an adversarial pre-prune
   review before moving anything** (it earned its keep here).
3. Archive the rest to `superseded/` (NOT `archive/` — check that repo's `.gitignore` too).
4. **Resolve the broken `nmd_orf_data.h5` symlink** — the deposited repo must be runnable,
   or document how the `.h5` is produced/fetched. This is the `cds_exons` lesson; do not
   deposit around it.
5. **Code-only rule:** `results_4ct/` derived TSVs are intermediates — ship producers, not
   products.
6. Render/verify a Fig-5 family output if feasible; then Pete cuts `v2.0.0`.

## NEW WORKSTREAM — background/scope: simplify the model repo (`NMD_orf_model_v5_4ct`)

Pete (2026-07-24): the model-development repo needs the **same** code-simplification
pass this repo got. Scope: **87 tracked files, 56 code**; prefixes show the same
canonical-pipeline-mixed-with-one-offs pattern (`05`, `09`×3, `patch_`, `relabel_`,
`make_shap_interpretation_figure`, `run_infer_all`). Same method: seed from manuscript
deliverables (Fig 5 + SF37–43 model exports, `03_train.py`, `model.py`, `11_kernel_shap_branches.py`,
`04_interpret_attention.py`, `06/08_export_deepshap`, the `09_*` exporters), build a keep-set,
archive the rest.

**Two ways it is NOT a copy-paste of this repo's prune:**
1. **It gets a DOI directly** — its keep-set *is* the deposited artifact, higher stakes
   than this repo's internal archive.
2. **The `nmd_orf_data.h5` feature file is a broken symlink** (`results_4ct/nmd_orf_data.h5`
   → a deleted `NMD_orf_model_v5/` path; found during the SF40 investigation). The `.h5`
   defines the scoreable cohort and is required to run inference / SHAP / attention. A
   deposited repo that can't produce or fetch it cannot be run — the `cds_exons` lesson
   again. **Resolve h5 provenance (regenerate from a tracked producer, or document a
   fetch) as part of the model-repo simplification, before its DOI.**

Also relevant: apply the **code-only** rule there too — `results_4ct/` holds many derived
TSVs (predictions, SHAP, attention exports) that are intermediates; the deposit ships the
producing scripts + a path to the `.h5`, not the derived TSVs.

## Phase 7 prerequisites (beyond the decision)

**Pete's, needed regardless:**
- ❗ **The manuscript has NO Code Availability section.** Phase 7 assumes one exists to cite
  the DOI into. Blocks the deposit.
- **Yul's sign-off** on the deposit (work is Leshem et al. + Castaldi); `CITATION.cff` must
  carry the full author list — Zenodo DOIs are immutable, so metadata must be right at mint.
- Add licenses to the two unlicensed repos; tag a release in each of the three.

**Technical, optional but valuable:**
- 4-CT-correct `isopair_wrapper/archive/02_build_profiles.R` (config at `:47-48` is still
  6-CT) and validate it as the canonical infrastructure builder. That would make Ruling 1
  true in practice and let `scripts/core/02–05` archive after all.

---

## Hard-won constraints — do NOT relearn these

1. **`.gitignore` is hazardous.** Three near-misses this session. Fixed: `results/` is now
   anchored to `/results/`, and `**/archive/` has exemptions. **54 tracked files are still
   force-tracked under `/results/`** — a fresh `git init` + `git add .` drops them silently.
   **Phase 7.1 must use `git add --pathspec-from-file` (fails loudly), never `git add .`,
   never a shell loop, never `-f`.** A directory named `archive/` is unusable as a
   destination (that is why the archive lives at `superseded/`).
2. **Static scanning does not work on this codebase.** Only 13–32% of reads use literal
   paths, and writes are just as dynamic (`sprintf`-constructed). A producer scan reported
   ~12 panel TSVs as producerless when their producer was directly readable. **Trust
   hand-inspection and actually running things.**
3. **The map is not a keep-list.** Uncited ≠ unused — proven three times
   (`ggplot_style.py`: 0 citations/21 importers; `scripts/core`; `nmd_predictor_comparison`
   produces SF43 with 0 citations).
4. **Snapshot ships CODE only.** Starting data → GEO **GSE329233**. No intermediates, no
   figures, no rendered reports. Consequence: every unshipped intermediate needs an in-repo
   producer.
5. **PNG byte-compare is the valid render check; PDF is not** (PDFs embed wall-clock
   `/CreationDate`). Verified three renders byte-identical post-archive.
6. **Donor codes are fine.** They match the public GEO/SRA submission
   (`~/claude_projects/ncbi_submissions/nmd_lung_cells/`). Not a blocker — do not re-open.

---

## Known-acceptable gaps (documented, not bugs)

- **The deprecated `05_final_report_mashr.Rmd` cannot run from a clean checkout** — its
  6-CT-era cache inputs (`orf_landscape.rds`, `orfik_scan.rds`, `unified_model.rds`) have no
  producer and exist only in the 6-CT backup. Acceptable because it is retained as SF25/30/31/32
  provenance, **not** as a runnable step. The snapshot README must say so.
- `utr5_features_all.rds`, `tx_summary.tsv` — same class, lower confidence. **Verify by
  running before the snapshot.**
- **M2 (nf-core/rnaseq launcher) is an accepted gap** (Pete's call) — the pipeline is
  standard and versioned in the Methods; do not re-open as a task.

---

## Separately open on the paper (independent of repo work)

- **§3 is parked** — 26 unverified claims.
- **Claim 1.11** — Yul is rerunning the AT2/LAE Spearman range and will update the text.
- Verification standing: §1 8/15 · §2 15/26 + 3 robust · **Tables 1–4 all exact** ·
  §4 complete · §5 near-complete. Two manuscript errors found and fixed this session
  (`9,161 → 9,154`; the 1.11 range).
