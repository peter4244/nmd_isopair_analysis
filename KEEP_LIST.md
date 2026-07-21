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

## Open before Phase 4 (prune)
1. **Pete's review of category C** — especially the ~49 event-detection lineage: confirm
   the Isopair package fully supersedes them.
2. **Empirical trace (1.1b)** — not yet run; would upgrade C from "no static reference"
   to "not opened at runtime." The data bundle is now local, so §1–3 scripts can finally
   be traced.
3. **`pheno/` IRB decision (8.1)** — 2 files with donor codes; gates whether they ship.
4. Non-code tracked files (526 total − 301 code = **225**) are **not yet classified** —
   plan 1.4 requires this before the snapshot.
