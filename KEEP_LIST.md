# KEEP_LIST — Phase 1 deliverable (repo consolidation)

Draft 2026-07-21. **Gate document: Pete reviews before ANY deletion.** Produced per
`REPO_CONSOLIDATION_PLAN.md` Phase 1. Nothing here has been deleted or moved.

## Method

**Seed** = every area that produces a manuscript deliverable or is proven load-bearing
infrastructure: `figures/` (Fig 3–5 multipanel, SF25–43, `figures/lib`), `code/upstream/`
(Yul: Fig 1–2, SF1–24, §1–3, Tables 1–4), `code/nmd_predictor_comparison/` (SF43),
`reproducibility/` (verifiers), `results/isoform_transitions/Version_6.0/` (Isopair §4),
`paper/` (docx builder, map).

**Candidates** = everything else. For each, test whether its basename is referenced by
any seed file *or* by `results_to_code_map.md`.

| | Count |
|---|---|
| Tracked code files (`.R/.Rmd/.py/.js/.sh`) | 301 |
| — Seed (producers + infrastructure) | **193** |
| — Candidates | 108 |
| —— referenced by seed or map | **2** |
| —— no inbound reference | **106** |

⚠️ **Static analysis only.** "No inbound reference" is evidence, not proof — the plan (F-2)
measured that only 13% of `readRDS()` and 32% of `fread()` calls use literal paths. The
empirical runtime trace (1.1b) has **not** been run. Treat the 106 as *candidates for
archive*, never as confirmed-dead. **Archive, don't hard-delete** (plan Phase 4.2).

---

## A. KEEP — seed (193 files)
All of `figures/`, `code/upstream/`, `code/nmd_predictor_comparison/`, `reproducibility/`,
`Version_6.0/`, `paper/`. Includes items with **zero** map citations that are nonetheless
required — the three proven cases: `figures/lib/ggplot_style.py` (0 citations, 21
importers), `Version_6.0/scripts/core/` (0 citations, sourced by the wrapper),
`code/nmd_predictor_comparison/` (0 citations, produces SF43). **This is why the map is
not the keep-list.**

## B. KEEP but reclassify — referenced candidates (2)
| File | Status |
|---|---|
| `code/gsea_mashr_2026.3.10.R` | Cited in map M8, but **verified 2026-07-21 as NOT the producer of Tables 3–4** (Yul's Rmds are). Superseded; keep for provenance, mark non-canonical. |
| `code/isocall_limma_dge_fullmodel_2026.3.1.Rmd` | Map explicitly records it as Pete's **parallel QC implementation**, not canonical (Yul's is). Keep as a documented cross-check. |

## C. ARCHIVE CANDIDATES — 106 files, no inbound reference

| Category | n | Notes |
|---|---|---|
| Pre-Isopair event-detection lineage | ~49 | `detect_events_*`, `fix_event_vectors_v3.0*`, `v3.0_distance_analysis*`, `v4.0_topology_analysis`, union-exon builders, `extract_exon_structures`, `recalculate_*`. The development ancestry of splice-event detection, **superseded by the Isopair package + `Version_6.0/`**. Coherent group — archive together. |
| Exploratory / one-off | 30 | `analyze_*`, `check_*`, `compare_*`, `investigate_*`, `debug_*`. |
| Explicitly archived already | 11 | `code/archive/` — already marked. |
| Dated analysis scripts | 5 | superseded by later dated versions. |
| Utility / shell | 4 | `cleanup_*`, `add_headers.sh`, monitors. |
| `code/nmd_atlas/` | 4 | **Pete rules (D-2)** — web app, 0 map citations, 0 manuscript mentions. Decision: exclude from snapshot, do NOT delete. |
| `code/ptcneg_go_handoff/` | 1 | check before archiving. |
| knitr `_files/` dir | 2 | regenerable render artifacts. |

---

## D. TRACE RESULTS (1.1b, run 2026-07-21) — three findings

**D-1 ✅ The original 106 are confirmed safe on the `source()` channel.** Audited all 37
`source()` calls in the keep-set: 21 literal, 16 dynamic. Every dynamic one resolves —
`VIS` → `isopair_wrapper/visualization_functions.R`, `mc_path` → `figures/lib/mechanism_class.R`,
`UPSTREAM` → the external `NMD_orf_model_v5_4ct` repo, the rest to `figures/lib`,
`isopair_wrapper`, `scripts/core`. **None points at any of the 106.**

**D-2 ⚠️ The seed was TOO COARSE — this is the trace's main catch.** Treating all of
`Version_6.0/` as seed swept in exploratory clusters. Map-citation by subdirectory:

| Subdir | files | cited | verdict |
|---|---|---|---|
| `isopair_wrapper/` | 13 | **9** | **KEEP** — the §4 pipeline |
| `scripts/core/` | 9 | 0 | **KEEP** — load-bearing (`event_detection_functions.R` sourced by `scripts/tests`) |
| `scripts/tests/` | 3 | 0 | KEEP — test infrastructure |
| `results/ptc/scripts/` | 10 | 0 | **→ ARCHIVE.** rmats appears **0 times** in the manuscript, supplement, *and* map |
| `scripts/dev/` | 7 | 0 | **→ ARCHIVE** (dev) |
| `scripts/archive/` | 10 | 0 | **→ ARCHIVE** (already so named) |
| `scripts/nmd/` | 11 | 0 | **→ REVIEW** — check before archiving |

⇒ **+27 archive candidates** beyond the original 106 (ptc 10 + dev 7 + archive 10), with
`scripts/nmd/` (11) needing a look. Found via the data-flow channel: these read the
**deprecated 2026.1.18 oarfish DGE outputs**, which the map records as superseded by the
2026.3.10 isocall pipeline — a self-consistent dead cluster.

**D-3 ⚠️ Divergent duplicate.** `visualization_functions.R` exists in **both**
`scripts/core/` (530 lines) and `isopair_wrapper/` (538 lines) and they **differ**. Every
cited reference points at the **wrapper** copy (confirmed via SF25's `VIS` absolute path),
so the `scripts/core/` copy is stale. Do not let the prune keep the wrong one.
Related: `scripts/core/` is the Isopair *development* source (event detection) and may be
superseded by the `peter4244/Isopair` package — same pattern as the ~49 `code/` files.
**Pete decision.**

**Also recorded:** `figures/multipanel/figure5_dl_model/figure5_panelA_architecture.R`
sources `../NMD_orf_model_v5_4ct/make_architecture_figure.R` — a **cross-repo dependency**
that must be resolved for the citable snapshot (plan Q4).

---

## E. PETE'S RULINGS (2026-07-21) — keep-set FINAL

**Ruling 1: the `peter4244/Isopair` package supersedes the in-repo event-detection code.**
**Ruling 2: `scripts/nmd/` reviewed → archive.**

Safety checks run before accepting both:

- **`isopair_wrapper` is self-contained.** Its `source("visualization_functions.R")` and
  `source("analysis_functions.R")` are **bare relative paths**, resolving to the wrapper's
  own directory — *not* `scripts/core/`. Verified in `05_final_report_mashr.Rmd:42` and
  `06_orf_analysis_mashr.Rmd:37`. Archiving `scripts/core/` therefore cannot break §4.
  ⇒ **Resolves D-3**: the wrapper's 538-line copy is live; `scripts/core`'s 530-line copy
  is stale and archives with it.
- **`scripts/tests/`** sources `scripts/core/event_detection_functions.R`, so it travels
  with core (test harness for superseded code).
- **`scripts/nmd/`** — zero inbound references from any keep-set file, and reads the
  deprecated `dge_isoform_nofilter_2026.2.7.rds`. It is the **predecessor numbered
  pipeline** (01→15), superseded by `isopair_wrapper`'s own 01→06. ⇒ archive.

### Final tally

| Class | files |
|---|---|
| **KEEP** | **~143** |
| **ARCHIVE** | **~158** |

Archive set = 106 original + `results/ptc/scripts` 10 + `scripts/dev` 7 +
`scripts/archive` 10 + **`scripts/core` 9** + **`scripts/tests` 3** + **`scripts/nmd` 11**.

Keep set = `figures/` (multipanel, SF25–43, `lib`), `code/upstream/` (Yul),
`code/nmd_predictor_comparison/`, `reproducibility/`, `Version_6.0/isopair_wrapper/`,
`paper/`, plus the 2 reclassified provenance files (§B).

**Coherent story:** everything archived is one of — (a) the pre-Isopair event-detection
lineage now owned by the `Isopair` package, (b) the predecessor numbered pipeline
superseded by `isopair_wrapper`, (c) exploratory/rmats work absent from the manuscript, or
(d) already-archived. Nothing archived is cited by the map or reachable from a keep-set
`source()`.

---

## F. DONOR CODES — RESOLVED, NOT A BLOCKER (2026-07-21)

Plan item **8.1 is closed.** The `DD###X` / `001V`-style codes in the repo are the **same
de-identified study codes already submitted for public release** in the GEO/SRA
submission (`~/claude_projects/ncbi_submissions/nmd_lung_cells/`, GSE329233, public
2026-09-01):

- `sample_metadata_longread.csv` → `sample_alias` = `AT2_001V_DMSO_LR`,
  `sample_title` = "AT2 donor 001V DMSO (long-read)", plus an explicit `donor_id` column.
- `sra/sra_biosample_attributes.tsv` → the **required, public** BioSample `*isolate`
  field is the donor code itself; `*age` = "not collected".

⇒ The repo does **not** exceed the disclosure already made to a public archive. Keeping
the codes is in fact **desirable** — a reader needs them to map repo files to GEO samples.
**No scrub, no IRB escalation, no snapshot blocker.**

**Corrections to my earlier flag, recorded so the mistakes aren't repeated:**
- I said "2 files"; it is **8** — I missed `pheno/fastq_sample_map2.fofn`, three rendered
  `code/*.html` reports, `isopair_wrapper/05_final_report_gencode_scope_2026-07-11.html`
  (in the **keep** set), and `REPO_CONSOLIDATION_PLAN.md` (where I wrote them myself).
- I framed it as "will they ship?" when `peter4244/nmd_isopair_analysis` is **public**
  (`"private": false`) — they already had. The right question was never *whether to
  publish* but *whether publication was already authorized*. It was.
- **Lesson:** check the existing disclosure record before treating something as a new
  disclosure decision. The answer was on disk the whole time.

---

## G. 🛑 PRE-PRUNE REVIEW FOUND A REAL BREAK — keep-set CORRECTED (2026-07-21)

**Verdict: the split as specified in §E was NOT safe.** An adversarial review found one
concrete way the prune breaks Figures 3, 4, 5 and SF26–34 + SF43. **Verified independently
before accepting.**

### The break
`isopair_wrapper/02_build_profiles_mashr.R:89–94` (**KEEP**, the §4 pipeline) loads seven
isoform-intrinsic objects and **builds none of them**:
`structures.rds`, `union_exons.rds`, `isoform_union_mapping.rds`, `cds.rds`,
`annotated_ue.rds`, `ptc.rds`.

A repo-wide `saveRDS` scan over all tracked code returns **four** writers — **every one in
the ARCHIVE set**:

| Producer | writes |
|---|---|
| `scripts/core/02_extract_isoform_structures.R` | `isoform_structures.rds` |
| `scripts/core/03_build_union_exons.R` | `union_exons.rds`, `isoform_union_mapping.rds` |
| `scripts/core/04_extract_cds_annotations.R` | CDS metadata |
| `scripts/core/05_annotate_region_types.R` | `isoform_union_exons_annotated.rds` |
| `scripts/nmd/06_filter_to_analysis_subset.R` | the 7 `*_filtered.rds` variants *(missed by the review too — found on independent check)* |

⇒ Under D-5 a reader gets through `01_prepare_data_mashr.R`, then
`02_build_profiles_mashr.R` dies on a missing `structures.rds`. **29 KEEP files** read the
downstream `profiles_c2/c4_*` or infra objects, including all three multipanel
`data_export.R` scripts and `nmd_predictor_comparison/01_extract_our_isoforms.R` (SF43).

### Why my Ruling 1 was wrong
"The Isopair package supersedes the in-repo event-detection code" is true of the
**functions** — `Isopair` exports `parseIsoformStructures`, `buildUnionExons`,
`extractCdsAnnotations`, `annotateRegionTypes`, `computePtcStatus` — but **not of the
driver scripts** that call them to produce these specific artifacts. No KEEP file calls
those functions against `data_mashr`. **Having the library ≠ having the pipeline step.**
The static scan missed it because the dependency travels through **data files**, not
script basenames — the same channel that caught the `scripts/nmd` cluster in §D-2.

### Corrections
**ARCHIVE → KEEP:** `scripts/core/02`, `03`, `04`, `05`; `scripts/nmd/01_prepare_expression_data.R`
(its `--source isocall` branch feeds `core/02`+`04`; §E's rationale cited the wrong
`--source oarfish` branch); `scripts/nmd/06_filter_to_analysis_subset.R`.
Remainder of `scripts/core` (`08`, `09`, `event_detection_functions.R`,
`reconstruction_functions.R`), `scripts/tests`, `scripts/nmd/09–15` stay ARCHIVE —
confirmed no KEEP file calls any of their 30 defined functions.

**Revised tally: ~151 KEEP / ~150 ARCHIVE** (also reconcile §E's "~143/~158" — a
`git ls-files` count gives 145/156; fix before using as a checksum).

### Confirmed-safe by the review (no change)
The ~49 `code/` lineage (spot-checked 4: all read the deprecated
`v4.0_reference_based/` tree, zero tracked files there); `scripts/dev`; `scripts/archive`;
`results/ptc/scripts`; `scripts/tests`; D-3's divergent duplicate; the `source()` and glob
channels; zero knitr `child=` in the keep set.

### Mechanism risks to honor at execution
1. **`.gitignore:27` is a bare `results/`.** Any archive move to a path *under* `results/`
   drops files silently. Use `git mv` to a path **outside** `results/`, and assert
   `git ls-files results/isoform_transitions/Version_6.0 | wc -l` before/after.
2. **Preserve full relative paths** — flat archive collides: `event_detection_functions.R`
   in two dirs, `visualization_functions.R` in core *and* kept wrapper, three distinct
   `01_prepare_*.R`.
3. Move the 5 orphaned rendered `.html` and `code/isoform_transition_splicing_analysis_files/`
   with their sources.

### Pre-existing hole (not caused by the prune, do not certify around it)
`cds_exons.rds` is read by KEEP `03b_rebuild_cache.R:54` and `05_final_report_mashr.Rmd:383`
but **nothing in the repo writes it**, tracked or untracked. `Isopair::extractCdsExons` is
exported and never called. Separately, `isopair_wrapper/archive/02_build_profiles.R` — the
*actual* working producer, calling the Isopair functions — is **untracked** (`.gitignore:27`)
and would not survive a fresh-init snapshot either. Resolve both before Phase 7.

### Documentation debt
`Version_6.0/METHODS.md` (§4 methods source of truth, incl. the provenance table at
:1427–1439) and `Version_6.0/README.md` document much of the archive set. Update them in
the same commit as the prune, or the §4 methods description goes stale.

## Open before Phase 4 (prune)
1. **Pete's review of category C** — especially the ~49 event-detection lineage: confirm
   the Isopair package fully supersedes them.
2. **Empirical trace (1.1b)** — not yet run; would upgrade C from "no static reference"
   to "not opened at runtime." The data bundle is now local, so §1–3 scripts can finally
   be traced.
3. ~~`pheno/` IRB decision (8.1)~~ — **RESOLVED, see §F.** Codes match the public GEO/SRA submission; not a blocker.
4. Non-code tracked files (526 total − 301 code = **225**) are **not yet classified** —
   plan 1.4 requires this before the snapshot.

---

## H. ✅ PHASE 4 EXECUTED + RENDER TEST PASSED (2026-07-21)

**Archive executed:** 158 files → `superseded/<original/path>`. Active tracked code
**301 → 157**. Nothing deleted; rollback branch `archive/pre-consolidation` (05713a4)
pushed to GitHub.

⚠️ **`archive/` was unusable** — `.gitignore:70`'s `**/archive/` would have silently
untracked all 158. Caught with `git check-ignore` pre-move; used `superseded/`.
**Third** near-miss from these globs (after the modelled fresh-init dropping 118 files,
and the bare `results/` rule hiding the real infra builder). **Fix `.gitignore` before
Phase 7.**

### Render test (plan 5.2) — the empirical backstop

Re-rendered one figure per family and byte-compared PNGs against the committed versions:

| Render | Family | Exercises | Result |
|---|---|---|---|
| `SF33_CdsAnd3UTR_GENCODE/figure_sf33.py` | SF, Python | `figures/lib`, six-validator gate | **byte-identical** |
| `figure3_isopair_and_ptc/figure3_panelC_event_prevalence.py` | multipanel, Python | multipanel data export | **byte-identical** |
| `SF25_SpliceEventCategories/render_sf25_canonical.R` | SF, **R** | sources `isopair_wrapper/visualization_functions.R` | **byte-identical** |

The third is the decisive one: it is the **D-3 divergent-duplicate** case, and it renders
correctly against the *kept* wrapper copy — empirically confirming the right one of the two
was archived.

Only the `.pdf` outputs differed (wall-clock `/CreationDate`), confirming the plan's N-2
finding that **PNG byte-compare is valid and PDF is not**. Tree restored clean.

⇒ **The archive is verified in practice, not just statically.**

---

## I. INTERMEDIATE-PRODUCER AUDIT (2026-07-22)

Driven by Pete's rule: **the citable repo ships CODE only** — starting data lives in GEO
(GSE329233); no intermediate data, figures, or rendered reports. That rule makes code
completeness the binding constraint: any intermediate that isn't shipped MUST have a
producer in the repo, or the chain breaks there.

### Method note — automated scanning FAILED, twice over
Regex-based read/write matching produced an unusable false-negative rate: writes are
constructed dynamically (`fwrite(file.path(OUT, sprintf("panel%s_long.tsv", slug)))` at
`figure4/data_export.R:240`), so ~12 panel TSVs were flagged as producerless when their
producer was directly readable. Consistent with the plan's F-2 finding (only 13–32% of
*reads* use literal paths) — writes are no better. **Do not trust a static producer scan
on this codebase.** The reliable checks are hand-inspection and actually running the chain.

### Hand-checked result — no `cds_exons`-class hole remains

| Input | Verdict |
|---|---|
| `ref_expression_fraction.tsv` | ✅ produced — `SF28_ReferenceShare/data_export.R:55` (scan false positive) |
| panel*/`profiles_c*`/`isoforms_per_gene` etc. | ✅ produced via constructed paths (scan false positives) |
| `orf_landscape.rds`, `orfik_scan.rds`, `unified_model.rds` | ⚠️ **no producer anywhere** — but **absent from the live 4-CT `data_mashr/analysis_cache/` (83 objects)** and present *only* in `data_mashr.bak_6ct_2026-07-11/`. Read solely by `05_final_report_mashr.Rmd` — **explicitly DEPRECATED** (map:258, "retained for reference", carries deprecation banners) — and `06_orf_analysis_mashr.Rmd`. ⇒ **vestigial reads of pre-rescope 6-CT artifacts, not live gaps.** |
| `utr5_features_all.rds`, `tx_summary.tsv` | ⚠️ same class — no writer found; lower confidence than the three above, **verify by running before the snapshot**. |

**Externally sourced by design** (correctly absent from this repo): `nmd_mashr_*`, `dge_*`,
`salmon.*` (Yul upstream) · `nmd_lungcells_*`, `sqanti3_corrected.cds.gff3` (`Isocall_v1`)
· `kernel_shap_*`, `predictions_*`, `model_comparison.rds` (`NMD_orf_model_v5_4ct`) ·
`gmap_*`, RBP rosters (third-party reference) · `pheno/*.csv` (starting metadata).

### Consequence for Phase 7
The live 4-CT chain has no known producer gap. The deprecated `05_final_report_mashr.Rmd`
**cannot run** from a clean checkout (its 6-CT cache inputs have no producer) — that is
acceptable *only* because it is retained as provenance for SF25/30/31/32, not as a runnable
step. **Say so in the snapshot README rather than implying the whole tree executes.**
