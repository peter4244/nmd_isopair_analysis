# Deprecated / dead / superseded code & data inventory — NMD repo

## ✅ CLEANUP EXECUTED 2026-07-11 (branch `cleanup-deprecated-code`, commits c738a0a / a534c3c / 8ff209e)
Done AFTER the re-scope landed on `main`. Every deletion was re-verified against the post-re-scope tree (zero live readers) before removal.
- **Category B (deleted):** `05d_protein_domain_analysis.R` (git rm) + its cache; 2 April `data_mashr_p030_*` backups (397 M); 6-CT prefloor backup (144 M); 2 `.html.bak`, `Rplots.pdf`, `.wrangler`, `.DS_Store`. **~568 M reclaimed.**
- **Category D — Pete ruled 2026-07-11:** (1) **Freeze legacy reports + delete the feeder cluster** → removed `04`,`04b`,`05l`,`05s`,`05m`,`05t`,`05u`,`05v` (+ caches), `paper_figures/`, `figures_mashr/` (legacy `.html` reports kept as frozen records). (2) **Delete protein track** → `05f`,`05g` (+ caches). (3) **Delete `deep_nmd_model/`** (superseded by linked repo; SF37/38 read the linked repo). **Another ~137 M.**  Total ≈ **705 M**, isopair_wrapper 1.3 G → 612 M.
- **Map hygiene:** fixed 4 stale refs in `paper/results_to_code_map.md` (claims 4.41/5.5 verified independent of the deleted scripts).
- **Corrections to this inventory (found on re-verify):** `05t`/`05v` feed the *retained* legacy report (they were D, not A auto-drop — the re-scope correctly left them, then Pete's freeze-decision removed them); `paper_figures/` was TRACKED (not untracked); SF42's `model_comparison` matches were the figure name, not an `.rds` read (05v confirmed dead).
- **STILL FLAGGED (not acted on):**
  - 5 per-figure dev docs carry stale provenance to deleted files (SF37/SF38 `README.md`, `figure3_panelA_methodology.md`, `figure4 RATIONALE.md`, `TD2_BIAS_AUDIT.md`) — several pre-date the deletion (SF37 already re-rendered natively per `d1d97f7`); need a figure-context pass, not a mechanical fix.
  - `isocall_pipeline.channing-backup/` (7.5 M) — a **separate git repo** at repo root; not touched (not mine to delete).
  - `data_mashr.bak_6ct_2026-07-11` — kept as the re-scope rollback backup; delete once fully confident.
  - `CLEANUP_CHECKLIST.md`-governed archive trees (`isopair_wrapper/archive/`, `Version_6.0/archive/`) — end-of-project only.

---


**Purpose:** Evidence-based inventory of dead, superseded, or orphaned code and data artifacts in
`~/claude_projects/nmd`, so cleanup can be done deliberately and safely. Every "dead" claim below is
backed by concrete evidence (grep output, banner text, dates, git tracking status). Being conservative
and correct is the goal — anything that cannot be *proven* dead is left in category D.

**Date context:** 2026-07-11. The 4-CT reference-floor **re-scope is planned but NOT yet executed**
(see `results/isoform_transitions/Version_6.0/isopair_wrapper/REFERENCE_FLOOR_PLAN.md`, Amendment +
phases R1–R8, all unchecked). The re-scope itself regenerates/renames/drops a specific set of items.

**⚠ GOLDEN RULE:** Cleanup happens **AFTER the re-scope lands**. Act only on category **B** (independently
dead) and **D** (after Pete rules). Category **A** is owned by the re-scope — do not touch it now;
**re-verify A** once R1–R8 complete. Category **C** is deliberately retained — never propose deleting.

---

## How "live" was defined

The live manuscript pipeline = everything reachable from:
- **Canonical report** `isopair_wrapper/05_final_report_gencode_scope_2026-07-10.Rmd` (to be renamed
  `_2026-07-11.Rmd` by the re-scope). Its `readRDS`/`source` (lines 41–52, 798–799) load:
  `analysis_functions.R`; `data_mashr/{profiles_c2,profiles_c4,cds,structures,expression_data,nmd_classification,dmso_samples}.rds`;
  `data_mashr/analysis_cache/{ref_atg_analysis,utr5_features_refaug,utr5_features_all,fc_c2_allsamples,fw_c2_allsamples}.rds`.
- **Live pipeline scripts** that produce those: `01_prepare_data_mashr.R`, `02_build_profiles_mashr.R`,
  `03b_rebuild_cache.R`, `05r_ref_atg_analysis.R`, `05k_utr5_all_isoforms.R`, `05k_b_utr5_refaug.R`,
  `analysis_functions.R`.
- **Live figures**: repo-root `figures/multipanel/figure{3,4,5}_*` + `figures/SupplementalFigures/SF24…SF42`.
- **Predictor comparison**: `code/nmd_predictor_comparison/01→04` → SF42.
- **Trace map**: `paper/results_to_code_map.md`.

A script/data artifact is a **dead candidate** if no live artifact reads it. References found only in
documentation (`results_to_code_map.md`, per-figure `RATIONALE.md`/`*_methodology.md`/`README.md`,
`TD2_BIAS_AUDIT.md`) are **provenance, not a live dependency**.

All isopair_wrapper paths below are relative to
`/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/`.

---

## Category A — owned by the 4-CT re-scope (leave to it; re-verify after R1–R8)

| path | kind | evidence | recommendation |
|---|---|---|---|
| `05t_ref_cds_features.R` | script (tracked) | PLAN Amendment + ROUND-3: "05t/05v DROPPED — dead/legacy". Output `data_mashr/analysis_cache/ref_cds_features_all.rds` read only by `05v` (`05v_model_comparison.R:47`). | Leave to re-scope (it drops 05t). Re-verify after R8. |
| `05v_model_comparison.R` | script (tracked) | PLAN ROUND-3: "05v feeds only `05_final_report_mashr.Rmd` (LEGACY) → dead for the paper". Reads 05l/05t/05u caches; writes `model_comparison.rds`. | Leave to re-scope (it drops 05v). |
| `data_mashr/analysis_cache/ref_cds_features_all.rds` | data | 05t output; consumer 05v being dropped. | Leave to re-scope. |
| `data_mashr/analysis_cache/model_comparison.rds` | data | PLAN: "non-paper cache". 05v output → read only by legacy Rmd (`05_final_report_mashr.Rmd:388`). | Leave to re-scope. |
| `05_final_report_gencode_scope_2026-07-10.Rmd` (+ `.html`, `_cache/`) | report | PLAN R3: `git mv` → `..._2026-07-11.Rmd`, re-render 4-CT. This is the *current* canonical report, only transiently "deprecated" by rename. | Leave to re-scope (rename, not delete). |
| `paper/section4_findreplace_2026-07-10_referencefloor.md` | doc | PLAN R6: `git rm` (never applied, 6-CT); replaced by fresh `section4_findreplace_2026-07-11_4ct.md`. | Leave to re-scope (it `git rm`s this). |
| `data_mashr/*` caches (profiles_c2/c4, fc/fw, ref_atg, utr5_*) | data | PLAN R2: regenerated 4-CT in place. | NOT deletion candidates — regenerated, not removed. |
| `data_mashr_prefloor_backup_2026-07-10/` (144 M, gitignored) | dir | 6-CT floor-pass before-snapshot; R2 will take a fresh backup, superseding this. | Leave to re-scope; safe local-delete after R8. |

---

## Category B — independently dead (safe to clean later; unrelated to re-scope)

Each row cites the grep showing zero *live* references.

| path | kind | evidence (zero live refs) | recommendation |
|---|---|---|---|
| `05d_protein_domain_analysis.R` (tracked) + `data_mashr/analysis_cache/protein_domain_full_analysis.rds` | script + data | Output basename `protein_domain_full_analysis` → `grep -rIln` across repo (excl. archive) = **zero hits anywhere** (not even in the orphan protein report 05g, which reads `protein_analysis_v2` not this). Protein/domain track has **zero** refs in `paper/`, `figures/`, `code/`. No SF for protein/domain. | `git rm` after re-scope. Highest confidence. |
| `05_final_report_mashr_p030_postfix.html.bak` (14 MB) | render backup | `git check-ignore` = ignored (untracked); Apr-29 `.bak` of a superseded render. | local `rm`. |
| `05_final_report_mashr_p050threshold_4CT.html.bak` (13 MB) | render backup | untracked `.bak`, Apr-29, p050 threshold (superseded by p030). | local `rm`. |
| `Rplots.pdf` | stray | untracked (gitignored); default R device auto-output, Apr-29. | local `rm`. |
| `data_mashr_p030_2026-04-29/` (133 MB) | dir (data) | `git check-ignore` = ignored (untracked); April pre-06-15 cache backup, superseded by live `data_mashr/`. | local `rm` after re-scope. |
| `data_mashr_p030_2026-04-29_postrender/` (264 MB) | dir (data) | untracked April cache backup; superseded. | local `rm` after re-scope. |
| `isopair_wrapper/figures_mashr/` (336 K) | dir (figures) | legacy-report figure outputs (`05g`/`06` `ggsave` target `fig_dir`); referenced only in a `*_methodology.md`/`README.md` (provenance). Current manuscript figures live at repo-root `figures/`. | local `rm` after confirming not needed for legacy .html. |
| `isopair_wrapper/paper_figures/` (15 MB) | dir (figures) | legacy figure outputs; no live reference (only doc). | local `rm` after re-scope. |
| `.DS_Store` files (wrapper root, Version_6.0) | junk | macOS metadata. | local `rm`. |

---

## Category C — deliberately retained (do NOT propose deleting)

| path | kind | evidence |
|---|---|---|
| `05_final_report_mashr.Rmd` (+ `.html`) | report | Banner L2–3: **"[LEGACY 2026-06-13] … Retained for sensitivity / reproducibility of prior submission."** ALSO still live-read: `figures/multipanel/figure3_isopair_and_ptc/verify_pass5_methods.R:163` does `readLines(...05_final_report_mashr.Rmd)`. Do not delete. |
| `06_orf_analysis_mashr.Rmd` (+ `.html`) | report | Banner L2–3: **"[LEGACY 2026-06-13] … Retained for … reproducibility."** |
| `05_final_report_gencode_scope_2026-06-15.Rmd` (+ `.html`, `_cache/`) | report | Banner L2–3: **"[LEGACY 2026-06-15 — pre-reference-share-floor] … retained for reproducibility of the prior state only. DO NOT cite its numbers as canonical."** (task-declared retained). |
| `isocall_dge/old/mashr/` | dir (data) | 6-CT isoform-mashr quarantine kept on purpose (ONBOARDING §3; CLAUDE.md). |
| `paper/archive/` | dir | Just-archived manuscript audits/plans (2026-06-15/18). |
| `isopair_wrapper/archive/` and `Version_6.0/archive/` (115 MB / 716 MB) | dir | Deliberate dev quarantine. `CLEANUP_CHECKLIST.md` governs their deletion **at end-of-project only**. |
| `CLEANUP_CHECKLIST.md` | doc | The end-of-project deletion checklist itself. |

---

## Category D — uncertain, needs Pete's judgment (do NOT auto-delete)

These are coupled to **retained legacy reports** (so deleting them breaks re-rendering the retained
`.html`s) or are entry-point reports (which legitimately have no inbound refs). Conservative call: leave
for Pete.

| path | kind | evidence | why uncertain |
|---|---|---|---|
| `04_productive_frameshift_precompute.R` + `data_mashr/analysis_cache/productive_frameshift_allsamples.rds` | script + data | Output read **only** by retained legacy `05_final_report_mashr.Rmd:596`. No live/canonical reader (PLAN Phase-0: "OUT of scope … read only by the deprecated legacy Rmd"). | Feeds a *retained* report. Droppable only if that report never re-renders. |
| `04b_all_c4_protein_comparison.R` + `all_c4_protein_comparison_allsamples.rds` | script + data | Output read only by legacy `05_final_report_mashr.Rmd:585`. | same coupling. |
| `05l_unified_model.R` + `unified_model.rds` | script + data | Read by legacy `05_final_report_mashr.Rmd:386` **and** by (dropped) `05v:46`. | feeds retained legacy report directly. |
| `05s_orfik_scan.R` + `orfik_scan.rds` | script + data | Read only by legacy `05_final_report_mashr.Rmd:4330,6477`. | same coupling. |
| `05m_orf_landscape.R` + `orf_landscape.rds` (+ `orf_landscape_per_orf.rds`) | script + data | Read by **retained** legacy `06_orf_analysis_mashr.Rmd:88`. | feeds a retained report. |
| `05u_paralog_annotation.R` + `paralog_genes.rds` | script + data | In-repo consumer is **only** `05v:48` (which the re-scope drops). ROUND-4 confirms it is not read by `code/nmd_predictor_comparison/`; the model's "paralog" split is a frozen H5 label in the `NMD_orf_model_v5_4ct` repo. | Becomes fully dead once 05v is dropped, **but** 05u may have been the manual upstream generator of the model's paralog split. Confirm with Pete. |
| `05g_protein_report_mashr.Rmd` (+ `.html`) | report | **Untracked** (`git ls-files` = none). Reads `protein_analysis_v2.rds`. No inbound refs; not in `results_to_code_map.md`; not bannered; protein track absent from manuscript. | Entry-point report — likely abandoned, but intent unconfirmed. |
| `05f_protein_analysis_v2.R` + `protein_analysis_v2.rds` | script + data | Output read **only** by orphan `05g` (`05g:101`). | dies with 05g if protein track is abandoned. |
| `deep_nmd_model/` (43 MB) | dir | Old in-repo CNN dev reports (`orf_model_report_v3/v4/v5.html`, bugfix notes). Superseded by linked repo `NMD_orf_model_v5_4ct`. Referenced only in `SF37/SF38 README.md` (provenance). | Whole-dir; provenance-referenced. Confirm superseded before removing. |

---

## Safe-deletion batch (category B, highest confidence)

Act on these **after the re-scope lands**. All are either untracked local clutter or a code/data pair
with proven zero live references.

1. `05d_protein_domain_analysis.R` + `data_mashr/analysis_cache/protein_domain_full_analysis.rds` — `git rm` (output read by nothing, incl. the protein report).
2. `05_final_report_mashr_p030_postfix.html.bak` — untracked; `rm`.
3. `05_final_report_mashr_p050threshold_4CT.html.bak` — untracked; `rm`.
4. `Rplots.pdf` — untracked stray; `rm`.
5. `data_mashr_p030_2026-04-29/` (133 MB) — untracked April cache backup; `rm`.
6. `data_mashr_p030_2026-04-29_postrender/` (264 MB) — untracked April cache backup; `rm`.
7. `isopair_wrapper/paper_figures/` (15 MB) + `figures_mashr/` (336 K) — legacy figure outputs; `rm` after confirming legacy `.html`s don't need re-render.
8. `.DS_Store` files — `rm`.

## Needs-Pete list (category D — one decision each)

1. **Legacy-report feeder scripts** `04`, `04b`, `05l`, `05s`, `05m` (+ their `.rds`): keep only if the
   retained legacy reports (`05_final_report_mashr.Rmd`, `06_orf_analysis_mashr.Rmd`) must stay
   *re-renderable*. If the retained `.html`s are frozen artifacts (no re-render), these can be dropped
   like 05t/05v.
2. **`05u_paralog_annotation.R` + `paralog_genes.rds`**: was this the manual generator of the frozen
   model's paralog split, or fully superseded? Its only in-repo reader (05v) is being dropped.
3. **Protein track `05g` (untracked report) + `05f` (+ `protein_analysis_v2.rds`)**: is the
   protein-domain / mass-spec analysis abandoned for the manuscript? (Zero manuscript refs.)
4. **`deep_nmd_model/` (43 MB)**: confirm fully superseded by the linked `NMD_orf_model_v5_4ct` repo
   before removing.

---

## Notes on directories & dated-file families

- **Dated-file family (canonical report):** three vintages coexist — `..._2026-06-15.Rmd` (C, retained
  legacy), `..._2026-07-10.Rmd` (A, current canonical → to be renamed `_2026-07-11`), and the pre-gencode
  `05_final_report_mashr.Rmd` (C, retained legacy). Only the newest is canonical; the older two are
  bannered and kept.
- **`data_mashr` family:** live `data_mashr/` (242 MB) is canonical; three sibling backups
  (`_p030_2026-04-29`, `_p030_2026-04-29_postrender`, `_prefloor_backup_2026-07-10`) are untracked and
  superseded (B / A).
- **Archive trees are OUT of scope for deletion now:** `isopair_wrapper/archive/`, `Version_6.0/archive/`,
  and `paper/archive/` are deliberate quarantines governed by `CLEANUP_CHECKLIST.md` at end-of-project.
- **`.html.bak` files** are the clearest independently-dead artifacts (untracked, superseded renders).
