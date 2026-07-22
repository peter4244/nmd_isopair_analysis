# Session handoff — NMD repo consolidation

Written 2026-07-22 at the end of a long session, for a fresh-context successor.
**Read this first, then the four docs listed under "Orientation".**

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
