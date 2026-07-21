# Plan v3 — Consolidate Pete + Yul code into one curated repo, then a citable snapshot

Supersedes v2 (2026-07-20, same day). v2 was critiqued by two independent agents;
both found fatal flaws, verified by hand before this rewrite. v3 is a structural
rewrite, not a patch.

---

## Goal

One curated repo in Pete's GitHub containing only code needed to reproduce the paper
+ Supplement, with deprecated and manuscript-unrelated code removed; then a clean
citable snapshot (Zenodo DOI) for Code Availability. **Reproducibility target = FULL.**

## What killed v2 (all verified against the repo, not taken on faith)

- **F-1 — the snapshot silently deleted the Isopair pipeline.** v2's Phase 7.3 said
  "`.gitignore` must carry the force-track exception for `Version_6.0/`." **No such
  exception exists.** `.gitignore:27` is a bare `results/`; there is no negation
  pattern in the file. Those 99 files are tracked only because they were `git add -f`'d
  — force-tracking lives in the **index**, not `.gitignore`. Measured: a fresh
  `git init` + `git add .` silently drops **118 tracked files, 99 under
  `Version_6.0/`, including `05_final_report_mashr.Rmd`** — the file v2's own Phase 4.3
  went out of its way to protect. `git add .` is silent on ignored files and `git
  status` comes back clean, so **no v2 gate could detect it**: Phases 1–5 all ran
  against the old repo where the files were tracked.
- **F-2 — the "dependency closure" was not a closure.** Measured on the real code:
  only **13%** of `readRDS()`, **32%** of `fread()`, and ~**50%** of `source()` calls
  use literal paths; the rest sit behind variables. And `figures/lib/ggplot_style.py`
  — v2's flagship example — is imported as a **bare module name** after
  `sys.path.insert(0, ...)` with `Path(__file__).parents[N]` arithmetic, so v2's
  described parser fails on its own headline case. The hand-reviewed residual is the
  majority of the graph, not a tail.
- **F-3 — v2 was code-first; completeness requires manuscript-first.** SF1–24 have
  **zero** map coverage (SF3/SF12/SF20/SF23 = 0 mentions each); only **19** SF dirs
  exist on disk for 43 SFs; SF20–23 have no render script in either repo; SF1–6 are
  produced in the **Isopair** repo; Tables 1–4 and ST1–ST4 have no producer and no
  plan entry (ST4 is cited **7×** in the manuscript). v2's Phase 5.1 gate would have
  reported "no orphans" only because the enumeration never happened.

## Corrections carried into v3
Yul's branch has **34** files (not 31) — the 3 READMEs carry her panel→script→input
contract and must be imported deliberately. **Three** linked repos, not two
(`Isopair`, `NMD_orf_model_v5_4ct`, `NMD_orf_model_v5`), plus an unversioned rsync
mirror at `results/isoform_transitions/NMD_orf_model_v5/` that must not ship.
`code/nmd_atlas/` has its own `.gitignore` that would apply inside a fresh init.

---

## D-5 (Pete, 2026-07-20) — REPRODUCIBILITY STANDARD, supersedes "FULL"

> "We need to make available the starting files (raw reads, which are provided in GEO)
> and the code that processes our data and produces the results in the paper, but we do
> not need to provide all interim data files."

**Standard = raw reads (GEO GSE329233) + complete code chain. Interim files NOT deposited.**
This is the conventional genomics standard and replaces the earlier "FULL reproducibility"
framing that drove v1–v3.1.

**What this REMOVES (large simplification):**
- The `nmd_fig_data/` bundle no longer needs a deposit → **group A of
  `code/upstream/DATA_INPUTS_NEEDED.md` is no longer a deliverable.**
- The SQANTI/isocall products no longer need a deposit → **group B's "no identified home"
  problem is MOOT.** (Do not resurrect it.)
- Most of workstream **D2** (deposit strategy) drops away.
- Phase 6.5's "large `.rds`" audit is moot — no intermediates ship.

**What this ELEVATES to critical — the catch:**
Without deposited intermediates, the code chain must be **unbroken from raw reads
forward**. A missing pipeline step used to be patchable by depositing its output; now it
is a hard break at step one.

**Ownership corrected 2026-07-20 (Pete): neither upstream gap is Yul's.**
- **M3 — CLOSED.** The long-read pipeline is **Pete's own `peter4244/Isocall_v1`**
  (verified reachable; ships `main.nf`, `Dockerfile`, `environment.yml`). ⚠ The working
  copy at `isocall_pipeline/` is **gitignored, not a submodule** — invisible to this repo
  and to any snapshot. Must be **cited by URL/DOI or made a real submodule**.
  **External repos the deposit must cover = THREE: `Isopair`, `NMD_orf_model_v5_4ct`,
  `Isocall_v1`.** (Corrected 2026-07-20 — I had briefly said four by including
  `NMD_orf_model_v5`. **v5 is DEPRECATED and must NOT be included:** it is never cited as
  a result source anywhere in the map, `model:` resolves exclusively to `v5_4ct`, and
  v5's AUPRC 0.781 does not match the manuscript's 0.833. Same reasoning excludes
  `Isoscope` / `Isovar`.)
  **Root cause of that error, worth not repeating:** I built the list from `CLAUDE.md`'s
  *Linked repos* table, which is a **navigation** aid listing repos Pete has — not a
  citation list of repos the paper needs. Same failure class as "the map is a keep-list":
  using a document for a purpose it was not built for. Also still true: the stale
  unversioned mirror at `results/isoform_transitions/NMD_orf_model_v5/` must not ship.
- **M2 — ACCEPTED GAP, closed by decision (Pete, 2026-07-20). NOT a task.** The
  nf-core/rnaseq run was John Ziniti's. His launcher/config will **not** be obtained.
  *(Retracting my earlier framing of this as "the single highest-priority gap" — that was
  written before Pete's call and is superseded.)*
  **Defensible because** nf-core/rnaseq is a standardized, versioned community pipeline
  with **v3.14.0 + Nextflow 24.04.4** pinned in the Methods; the launcher is
  parameterization, not novel logic. **Consequence to state plainly:** the chain is not
  literally one-command runnable from GEO → salmon counts; that step is reconstructed
  from the Methods, which therefore must stay accurate. **⇒ There is now NO blocking
  upstream code gap** (M3 closed via `Isocall_v1`, M2 accepted).

**What is UNCHANGED — do not conflate deposit with verification:**
We still cannot *run* Yul's §1–3 scripts locally without her intermediates, so:
- Phase **1.1b (empirical trace)** remains blocked for §1–3 → those stay hand-reviewed.
- Phase **5.2 (re-render)** can only verify §4/§5 here. **Someone must still verify the
  §1–3 chain actually runs** — most practically Yul, on Channing. Shipping code we have
  never confirmed executes end-to-end is the residual risk of this standard, and it is
  a real one.

## Decisions settled

- **D-1: prune in place on a branch, THEN fresh-init a snapshot.** `.git` is 1.1 GB
  vs a 74 MB tracked tree; in-place deletion frees none of it, a snapshot alone leaves
  the working repo uncleaned.
- **D-2: `code/nmd_atlas/` → exclude from snapshot, do NOT delete** (8 files, 0 map
  citations, 0 manuscript mentions; a real artifact Pete may publish separately).
- **D-3: `code/nmd_predictor_comparison/` → KEEP** (produces SF43; NMDetective-B /
  NMDEP cited 4× in the manuscript). Note it has **0 map citations** — a third
  independent proof that the map is not a keep-list.
- **D-4 (Pete, 2026-07-20): the manuscript/supplement `.md` files are NOT canonical
  and do not belong in the repo.** The Google Doc is the source of truth. **Consequence:
  the deliverables ledger (Phase 0.5) becomes the only in-repo record of what the paper
  requires, so it must be a tracked, durable artifact — not a scratch file.**
- Home = Pete's repo. Import Yul's files, not her git history (unrelated histories).
- Keep necessary infrastructure, not just result-generating scripts.

---

## Plan

### Phase 0 — Unblock — ✅ FULLY RESOLVED 2026-07-20

All six questions answered by Pete + Yul, each verified against `yul/main`
(now at commit `78f1a8c`, 11 commits). **Phase 0 is closed; Phase 3 (import) unblocked.**

- **O-1 — Yul's GitHub IS canonical.** The `/udd/reyle/.../code/final/` blocker is DEAD.
- **O-2 — download, and now COMMITTED.** Both rosters live in `yul:data/`:
  `encode_rbp_roster_vannostrand2020.csv` (357 rows),
  `gerstberger_2014_rbp_census.csv` (1,543 rows). Verified.
- **O-3 — `NMD_shortread_dge_fullmodel_2026.5.5.Rmd` IS the script.** Yul confirmed the
  later date is fine. Provenance-vintage worry (R-E) CLOSED. No duplicate to reconcile.
- **O-4 — FB excluded DELIBERATELY** (outlier in one analysis). `exclude_fb = TRUE` is
  intentional. Documentation fix only: manuscript must state the 86.2%/13.5% split is
  AT2/LAE/MV while the same section reports FB in the logFC medians. (§3 parked.)
- **O-5 — SF20–23 producers PUSHED.** `yul:Figures/render_output_lost_gene.R`,
  `render_output_lost_per_isoform.R`, `render_output_lost_threshold.R`,
  `render_sr_lr_correlation.R`. Verified. **R-F (largest completeness hole) CLOSED.**
- **O-6 — SF1–6 now render IN YUL'S REPO, not Isopair-only.** 6 new `render_*.R`
  (`render_isoform_length.R`, `render_pairwise_expression.R`, `render_sqanti_categories.R`,
  `render_mashr_sharing.R`, `render_proportion_vs_expression.R`, `render_sr_lr_correlation.R`),
  all wired into `make_supplemental_figures.Rmd`. **DECISION FLIP (good direction):
  Isopair can now be CITED BY DOI, not vendored** — Q4 simplifies. Verified.

**Yul's 2026-07-20 push also added TWO NEW WRINKLES (see W-1, W-2 below).**

Original list (all now answered — kept for provenance):
- **O-1** Is her GitHub now canonical, or does `/udd/reyle/.../code/final/` hold more?
  *(Phases 3–5 are unsafe if not canonical.)*
- **O-2** Is `encode_rbp_roster_vannostrand2020.csv` curated or a straight download?
  *(deposit vs cite)*
- **O-3** Which script emitted the manuscript's `2026.3.10` CSVs?
  `NMD_shortread_dge_fullmodel_2026.5.5.Rmd` post-dates them. **Do not resolve from
  the map alone** — picking wrong silently changes provenance.
- **O-4** Why is FB excluded as a "provisional CT" in `productive_response.Rmd`?
  *(affects §3 text; parked by Pete but must not be lost)*
- **O-5** Do producers exist anywhere for **SF20, SF21, SF22, SF23**? Neither repo has one.
- **O-6** Are SF1–SF6 in fact rendered in the Isopair repo, per her `Figures/README.md`?

### Phase 0.5 — Deliverables ledger (NEW; the plan's spine)
**This inverts v2's direction of travel.** v2 derived the keep-set from code and
checked it against the manuscript at the very end — backwards for a completeness
guarantee.

- **0.5.1** Mechanically extract every deliverable from the Google Doc export +
  `Supplemental Figures NMD.md`: Fig 1–5 panels, SF1–SF43, Table 1–4, ST1–ST4,
  13 Supplemental Methods subsections, 10 Supplemental Results subsections (~90 rows).
- **0.5.2** Columns: deliverable · producing script (or `NONE`) · repo (pete / yul /
  Isopair / v5_4ct / external) · inputs · status.
- **0.5.3** Commit as `paper/DELIVERABLES_LEDGER.md`. Per **D-4** this is the repo's
  only record of the paper's requirements — it must be tracked and kept current.
- **0.5.4 (NEW, v3.1) PROVENANCE STAMP — required, or the ledger is unverifiable.**
  D-4 means the ledger is derived from inputs that will never exist in the repo, so
  nothing in-repo can regenerate it, diff it against the paper, or detect drift — yet
  5.1 makes it the authoritative gate. Mitigate by committing, alongside it: the
  **Google Doc revision ID**, the **export date**, and the **extraction script**. Then
  "is the ledger current?" is answerable by re-exporting and diffing rather than by
  memory. Without this the ledger is an assertion, not a record — and since the Doc
  keeps moving while the repo freezes for a DOI mint, **drift is the default outcome,
  not a risk**.
- **0.5.5** Keep the map as an independent cross-check. v3 said "the ledger wins on any
  conflict," which collapses two imperfect records into one unchallenged authority
  (**N-5**). A single Phase-0.5 extraction error would otherwise propagate through 5.1
  and 5.5 into the DOI. Conflicts should be *surfaced*, and resolved case by case.
- **Expected day-one output:** SF1–24, SF20–23, ST1–ST4, Tables 1–4 land as `NONE`.
  That is the point — surface the holes *before* the prune, not after.
- **This is a prerequisite for Phases 1, 2, and 5, not a parallel task.**

### Phase 1 — Keep-set: seed → hints → reviewed residual (no writes)
- **1.1 Seed** = every producer named in the ledger + figure/SF render scripts +
  `paper/build_supplemental_figures_docx.js` + reproducibility verifiers +
  `code/nmd_predictor_comparison/` + Yul's 34 files. **Seeded from the ledger, not
  the map** — the map resolves only 27 of its 117 cited paths exactly.
- **1.1b (NEW, v3.1) EMPIRICAL TRACE — run this BEFORE any hand review.** v3 asked
  Pete to hand-classify 500+ files; that is multi-day work whose failure mode is
  exactly the silent one (a file `source()`d from a single Rmd chunk). A runtime trace
  dominates it: re-render the figure families under an `Rprofile` / `sitecustomize.py`
  shim (or `fs_usage -w -f filesys`) that logs every `file()` / `readRDS()` / `open()`
  argument. This yields the **true** dependency set empirically and subsumes every
  channel 1.2 parses — constant-variable paths, globs, CWD conventions, `parents[N]` —
  because it does not care how the path was computed.
  **Limit:** the trace only covers stages that actually execute. Per 6.1, §1–3
  "reproduce from nothing today," so those still need hand review — but for a small
  identified subset, not for 500 files.
  *(Corollary: v3 over-stated F-2. 28 of 30 `sys.path` manipulations are two literal
  forms — `HERE.parents[2]/"lib"` and `HERE.parents[1]/"lib"` — followed by a
  stereotyped `sys.path.insert` + `from ggplot_style import`. That is two regexes, not
  an intractable dynamic-import problem. Static resolution is weak, not useless.)*
- **1.2 Dependency hint generator** (demoted from "crawler"/"closure" per F-2). Parse
  `source()`, `import`/`from`, `library()`, literal paths — AND, because they carry
  most of the graph: constant-variable path assignments (`X <- "..."` then
  `readRDS(X)`), `sys.path.insert` + `parents[N]` module resolution, glob reads, and
  per-file CWD convention (repo-root-relative vs file-relative). Emit
  `keep_hints.tsv`: file → reached-from → channel → confidence.
- **1.3 Reviewed residual is the PRIMARY artifact.** Every file not in seed or hints
  gets a hand class with evidence: **(i)** missed dependency, **(ii)** deprecated
  (per the A–D methodology in `DEPRECATED_CODE_INVENTORY.md`), **(iii)** alive but
  manuscript-unrelated, **(iv)** Pete rules.
- **1.4 Classify NON-code too.** v2 covered only the 257 code files; **222 tracked
  non-code files** (479 total) had no class — including map data inputs, `figures/lib`
  style SSOTs, and ~13 MB of TSV/`.md` under `Version_6.0/`.
- **1.5 Deliverable:** `KEEP_LIST.md`, every one of the 479 + 34 files in exactly one
  class. **GATE: Pete reviews before any deletion.**

### Phase 2 — Close map holes (no code writes)
- **2.1** Add `nmd_predictor_comparison` → SF43.
- **2.2** Add the **tables → code** section (absent today, so v2's table gate was vacuous).
- **2.3** Reconcile map ↔ ledger; the ledger wins on any conflict.
- **2.4** Re-run the "every cited file exists" check after the Yul import.
- **Parked (Pete, 2026-07-20):** §3 / M10 / the wrong `productive_compensation.Rmd`
  filename. **Must be revisited before the snapshot** — do not ship a known-wrong entry.

### Phase 3 — Import Yul's 34 files (writes, additive only)
- **3.1** Copy into `code/upstream/`, preserving her `Figures/` and `tan_reanalysis/`
  substructure and all 3 READMEs. Files only, no history.
- **3.2** Rewrite `yul:` citations to in-repo paths.
- **3.3** Reconcile the duplicate SR DGE analyses per **O-3**.
- Additive only — cannot break anything.

### Phase 4 — Prune (writes, reversible)
- **4.0** Pin the remote by name in every push command. **`origin` has two push URLs**
  (public GitHub + Channing GitLab) and the **`yul` remote has a push URL to the
  collaborator's repo** — remove it before any push step.
- **4.1** Branch `repo-consolidation`; push `archive/pre-consolidation` as the rollback
  point. Resolve the 4 pending working-tree deletions first.
- **4.2** Delete only classes (ii) and (iii), after Pete's gate.
- **4.3** Do NOT delete `05_final_report_mashr.Rmd` (SF25/30/31/32 provenance).
- **4.4** Update `README`, `CLAUDE.md`, `ONBOARDING.md`, map, ledger to the new layout.

### Phase 5 — Verify (GATE) — runs TWICE
- **5.1 Trace against the ledger** (not the map): every one of the ~90 deliverables →
  exactly one surviving producer, or an explicit accepted `NONE`.
- **5.2 Re-render check — split by format (v3.1 correction).** v3 claimed byte-compare
  "fails 100% of the time." That is true for **PDF** (each embeds a wall-clock
  `/CreationDate`) but **false for PNG**, the format the docx builder actually consumes:
  tested three tracked PNGs, **0 timestamp markers** in each. So:
  - **PNG → plain byte-compare.** Cheap, exact, no infrastructure. Valid once 6.5 pins
    the matplotlib version (the only non-determinism is a version string).
  - **PDF → `SOURCE_DATE_EPOCH` / strip `pdf.infoDict`, or compare rasterized content.**
  Re-render one figure per family (Yul panel, multipanel, SF).
- **5.3** `figures/lib/lint_sf_legends.py` + pass-7 / cross-check verifiers.
- **5.4** 5-step verification discipline applied to the ledger and map.
- **5.5 RE-RUN 5.1–5.3 AFTER Phase 7.1**, against the snapshot repo. v2's gate
  validated a tree that Phase 7 then mutated — that structural gap is what made F-1
  invisible.

### Phase 6 — Data reproducibility
- **6.1 D1 inventory:** every kept script → input, size, source. Known: Yul's
  `nmd_fig_data/` bundle is in neither repo; her Rmds read `/udd/reyle/...`; the
  in-repo `shortread_dge/` + `isocall_dge/` CSVs are her *outputs*, not inputs — §1–3
  reproduce from nothing today.
- **6.2 D2 deposit strategy.** Note **GEO GSE329233 is an expression record and will
  not hold** DGEList `.rds`, mashr posterior CSVs, or `srsf.gtf`/`.cds.gff3`/`.fasta`.
  A companion Zenodo **data** record is likely required. Name an owner.
- **6.3 D3 path normalization — scope corrected 20× (v3.1).** v3 named 5 files.
  **Measured: 99 tracked files contain `/Users/petecastaldi` or `/udd/reyle`, 35 of
  them in keep-likely dirs** (`figures/`, `code/nmd_predictor_comparison/`,
  `reproducibility/`, `paper/`). That includes
  `code/nmd_predictor_comparison/01_extract_our_isoforms.R` (a D-3 KEEP) and the
  SF25–39 render/export scripts — i.e. the very files 5.2 re-renders. It also includes
  **`paper/build_supplemental_figures_docx.js:29`**
  (`const SFROOT = '/Users/petecastaldi/claude_projects/nmd/...'`), which is in the 1.1
  seed and was absent from v3's list. This is a real work item, not a half-day.
  **Must land BEFORE 5.2**, or the deposited code is not the code that was verified.
- **6.4** Per stage: "run from raw" vs "run from deposited intermediate" (nf-core,
  isocall, SQANTI3, DL training are compute-heavy).
- **6.5 (corrected, v3.1) The size risk is the opposite of what v3 stated.** The large
  root `.rds` (`all_proportions_2026.1.15.rds`, `all_proportions_2026.1.20.rds`) are
  **untracked and ignored** — verified — and 7.1 transfers by explicit `git ls-files`,
  so they can never reach the snapshot. There is no 100 MB breach: the largest *tracked*
  file is ~10 MB. The real issue is ~30 MB of **tracked, regenerable** HTML/docx/PDF
  that probably should not ship. Audit for shippability, not for size limits.

### Phase 6.5 — Environment capture (NEW; gates the snapshot)
There is **no `renv.lock`, `requirements.txt`, conda env, `DESCRIPTION`, or captured
`sessionInfo()`** anywhere in the tree. Without this, "FULL reproducibility" means
"runs on Pete's laptop in July 2026."
- `renv::init()` + committed `renv.lock` (R); pinned `requirements.txt` (Python:
  `figures/lib/` + Yul's three matplotlib stitchers).
- Pin versions/DOIs for nf-core/rnaseq v3.14.0, Nextflow 24.04.4, PacBio Isocall,
  SQANTI3, TD2.
- **Cache the `pathview` KEGG `hsa04141` template into the repo** — currently fetched
  from the network at runtime and will drift.
- Capture `sessionInfo()` from the run that produced the final figures.

### Phase 7 — Citable snapshot

> **v3.1 correction (tested, 2026-07-20).** v3's original 7.0 pattern
> (`results/` + `!results/isoform_transitions/Version_6.0/**`) **DOES NOT WORK** —
> measured in a scratch repo it staged **0** files. Git cannot re-include a file whose
> **parent directory** is excluded; a bare `results/` prunes the directory and no `!`
> beneath it is ever consulted. The **primary F-1 control is 7.1, not 7.0.**

- **7.1 (PRIMARY CONTROL) Transfer by explicit file list — NEVER `git add .`.**
  Measured: `git add` of an ignored path **exits 1 and stages nothing**, both for a
  literal path and for `--pathspec-from-file`. `git add .` is silent. So an
  explicit-list transfer makes F-1 *structurally impossible*, whatever `.gitignore` says.
  - Use `git add --pathspec-from-file=<list>` built from `git ls-files` on the curated
    branch. **Do NOT use a shell loop** — each iteration exits 1 but the loop continues
    and the operator watches hints scroll past.
  - **`-f` / `--force` is FORBIDDEN in this phase.** It is the obvious reflex when
    `git add` complains, and it restores F-1's full blast radius while looking like progress.
  - Then assert the file count matches exactly and diff the two lists; fail loudly.
- **7.0 (now belt-and-braces) `.gitignore` rewrite.** Correct form unexcludes each
  ancestor one level at a time, using `dir/*` not `dir/`:
  ```gitignore
  results/*
  !results/isoform_transitions/
  results/isoform_transitions/*
  !results/isoform_transitions/Version_6.0/
  ```
  Verified: stages the target file. **Second-order trap:** `.gitignore:70-71`
  (`**/archive/`, `**/archive_*/`) is a *later* rule that still wins — **29 tracked
  files sit under an `archive/` dir**, including
  `Version_6.0/scripts/archive/event_detection_functions.R` and `ptc_analysis.R`,
  names that are not obviously disposable. These are currently masked because
  directory-pruning short-circuits before `**/archive/` is consulted, so **fixing the
  `results/` chain unmasks a previously invisible ignore rule** — the exact class of
  bug that produced F-1. Also re-check the global `*.rds` / `*_proportions_*.tsv`
  patterns, which only bind once the directory is re-included.
  Drop the `code/nmd_atlas/` `.gitignore` line — it ignores only caches, and per D-2
  `nmd_atlas` is excluded from the snapshot anyway.
- **7.2 Phase 7.1 is IRREVERSIBLE** — it discards the `.git` that holds
  `archive/pre-consolidation`. The rollback survives only if pushed to a reachable
  remote first. Verify the push before proceeding.
- **7.3** README (repro steps + deliverable↔script table from the ledger), **LICENSE**,
  data-availability pointer.
- **7.4** Decide inclusion for the three linked repos; ensure the unversioned
  `results/isoform_transitions/NMD_orf_model_v5/` mirror does **not** ship.
- **7.5** Re-run Phase 5 against the snapshot (= 5.5).

### Phase 8 — Submission compliance (NEW; blocks Zenodo)
- **8.1 IRB / data-steward decision** on donor codes (`DD027U`, `DD001V`, …) in
  `pheno/nmd_pheno_longreadbamids_2026.1.18.csv` and `pheno/fastq_sample_map3.tsv`,
  joined to treatment and cell type. Whether these may appear in a permanent, public,
  DOI-minted archive is an ethics call, not a prune decision. Sweep all tracked
  `.csv`/`.tsv`/`.rds`/`.txt` for the same pattern.
- **8.2 `CITATION.cff`** with the full author list + **Yul's sign-off** — Zenodo DOIs
  are immutable, so the metadata must be right at mint time.
- **8.3 Manuscript edits** (Pete, in the Doc): add a **Code Availability** section —
  there isn't one, and Phase 7.5 assumes it exists; add `NMD_orf_model_v5_4ct` to Data
  Availability; confirm ethics statement + Author Contributions.
- **8.4** Untrack the Excel lock file `PJC_Tables/~$NMD_Tables.xlsx`.

---

## NEW WRINKLES from Yul's 2026-07-20 push (verified)

- **W-1 — Git LFS.** `yul:.gitattributes` routes `tan_reanalysis/data/*.xlsx` through
  Git LFS (`filter=lfs`). The four Tan tables are LFS objects (S2 alone = 16 MB). **LFS
  objects do NOT travel through a plain `git init` snapshot** — a fresh init sees only
  the 3-line pointer stubs unless LFS is installed and the objects re-pushed. Phase 7.1's
  explicit-list transfer does not currently handle LFS.
  → **Phase 7 action:** decide per file — either de-LFS (commit the real bytes, fine if
  <100 MB/file; S2 at 16 MB is OK) or carry LFS config into the snapshot and re-upload
  objects. Simplest: de-LFS, since these are external inputs we may not ship at all (W-2).
- **W-2 — Tan tables are now redistributed.** Yul committed the four Tan et al. (2025)
  `Supplementary_Table_*.xlsx` into her repo — but her own earlier README said they are
  *"not redistributed here — the authors' published data."* Committing another paper's
  supplementary tables into a **public, DOI-minted** archive is a redistribution-rights
  question. Fine in her working repo; **must be resolved before the snapshot.**
  → **Phase 8 action:** default to NOT shipping them — replace with a download manifest
  (the file→sheet table already in `tan_reanalysis/README.md`). Removes W-1 entirely.
- **W-3 — path collision.** Yul's data now lives at `yul:data/`, which will collide with
  our layout on import. → **Phase 3 placement decision** (`code/upstream/` per Q1 keeps
  it namespaced).

## EXECUTION ORDER (v3.1 — the phase numbers are NOT the run order)

v3's numbering had four forward dependencies, two of which it admitted in the step text
and left unreordered. Run in **this** order:

```
0     Ask Yul (O-1…O-6)
0.5   Deliverables ledger  (+ provenance stamp)
3     Import Yul's 34 files      ← MOVED UP. Additive-only, "cannot break anything",
                                   and 1.2/1.5 both reason over her files
1     Keep-set: 1.1 seed → 1.1b EMPIRICAL TRACE → 1.2 hints → 1.3-1.5 reviewed residual
6.3   Path normalization (99 files)   ← must precede 5.2
6.5   Environment capture             ← must precede 5.2 (PNG byte-compare needs a
                                        pinned matplotlib; else 5.2 tests this laptop)
8.1   IRB / donor-code decision        ← gates the keep-list: whether pheno/ ships is a
                                        keep/drop input, and v3 had it 4 phases late
2     Close map holes (incl. 2.4, which explicitly needs Phase 3 done)
4     Prune  →  5  Verify  →  6  Data  →  7  Snapshot  →  5.5  RE-VERIFY  →  8  Rest
```

**Moving Phase 3 ahead of Phase 1 is the single highest-leverage change** — it is
explicitly additive-only, and three later phases silently assume it has already run.

**8.1 de-scoped (N-4):** v3 priced the donor-code sweep as an unbounded audit. Measured,
the `DD###X` pattern appears in exactly **2 tracked files**, both already named. It is a
one-line check plus one ethics decision — do not let it become a blocker by inertia.

## Risks
- **R-A: the hint generator misses a live file.** Measured coverage is low (F-2).
  Mitigations: reviewed residual is primary, not a tail; Phase 5.2 deterministic
  re-render is the empirical backstop; archive branch + explicit-list transfer.
- **R-B: "manuscript-unrelated" is a judgment call** (`nmd_atlas` is the type case).
  Pete rules; never automate.
- **R-C: Phase 7 file loss** — the F-1 mechanism. Mitigated by 7.0 + 7.1's count
  assertion + 5.5's post-snapshot re-verification. Three independent controls because
  this failure is silent.
- **R-D: §3 parked with a known-wrong map entry.** Must not ship.
- **R-E: duplicate SR DGE** — gated on O-3.
- **R-F: SF1–24 have no in-repo producers and no map coverage.** The single largest
  completeness hole. Gated on O-1/O-5/O-6.
- **R-G: `nmd_fig_data/` may have no valid deposit target** (6.2). If GEO can't hold
  it and no Zenodo data record is created, §1–3 never reproduce and FULL is unmet.

## Open questions for Pete
- Q1 Location for imported Yul files (`code/upstream/`?).
- Q2 Archive branch vs `archive/` directory for pruned code.
- Q3 Which derived data ships with the citable repo vs external deposit.
- Q4 Isopair / `v5_4ct` / `v5`: submodule, vendored, or cited by DOI? (Q4 interacts
  with R-F: if SF1–6 render in Isopair, "cited by DOI" loses six supplements.)
- Q5 Does `nmd_atlas` get its own repo?
- Q6 Who owns creating the `nmd_fig_data` deposit — Pete or Yul?
