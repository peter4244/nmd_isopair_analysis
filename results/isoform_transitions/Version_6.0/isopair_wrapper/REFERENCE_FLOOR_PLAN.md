# Reference-share floor — implementation plan & check-off

**Branch:** `reference-floor-25pct` · **Started:** 2026-07-10 · **Owner:** Pete (Castaldi)

## Decision (locked)
- **Family A** — reference *selection* unchanged: rank-1 by DMSO mean CPM within the **strict non-NMD** pool (`nmd_class[[ct]]$non_nmd`). Do NOT re-anchor to the gene dominant.
- **25% floor** — additionally drop any gene whose selected reference accounts for **< 25% of total gene expression**.
- **all-isoform denominator** — the 25% is `mean_dmso(reference) / Σ_{all isoforms of gene} mean_dmso`, DMSO basis = `dmso_samp[[ct]]` per profile. NB "all isoforms" = all isoforms **passing the 5% condition-stratified filter** in `01_prepare_data_mashr.R` (expr_mat is already filtered); Methods must not overstate this.
- **Rationale:** analyze only genes whose reference is *truly non-NMD* AND *truly dominant*. Displacement cases (true dominant is a "neither"-gap isoform) are dropped, not re-anchored. See `REFERENCE_FLOOR_RATIONALE.md`.

---

## ⚠ AMENDMENT 2026-07-11 — FULL 4-CT re-scope (decision **b**; supersedes the 6-CT floor pass entirely)

**Trigger:** the SF26 legend audit exposed that the `all_samples` profile used a 6-CT sample basis (18 DMSO / 18 Smg1i libs across AT, DD, **DD_ALI, DO**, FB, MV) — DD_ALI + DO are out of manuscript scope. An independent adversarial review (2026-07-11) then found the leak is **deeper than the two DMSO/Smg1i pointers**: the isoform *universe itself* is built 6-CT. Full leak map:

| # | Leak | Site | affects |
|---|---|---|---|
| 1 | **Isoform filter / universe** (`expr_mat`) | `01` L130–131 (`dmso_cols`/`smg1i_cols`), L146–147 (5% max-prop filter), L152–155 (`filterByExpr` on 6-CT design) | which isoforms exist at all; SF28 discloses 36→26 sample drop |
| 2 | Reference selection | `02` L154–160 ← `dmso_samples[["all_samples"]]` (`01` L204) | which reference isoform |
| 3 | 25% floor denominator | `02` L173 ← same | which genes retained |
| 4 | C2 NMD-partner | `02` L196–202 ← `smg1i_samples[["all_samples"]]` (`01` L205) | which NMD isoform analyzed |
| 5 | SF42 predictor-comparison cohort | `code/nmd_predictor_comparison/01_extract_our_isoforms.R:45-46` reads `profiles_c2/c4_allsamples.rds` (floor+rescope-affected) | SF42 (NMDetective-B/NMDEP). Re-run sub-pipeline (round 3); **05t/05v dropped — dead/legacy** |

**Clean (verified by review):** `all_samples` NMD *classification* (`01` L198–203, union/intersect over AT/DD/FB/MV) is already 4-CT; `buildProfiles` is structural (no per-sample leak); per-CT profiles never leaked; `05r`/`05k`/`05k_b` do **not** reconstruct a sample basis (structural / classification-keyed → inherit 4-CT free); pop_BC producers (`figure5_dl_model/data_export*.R`) and SF39/40 are basis-free and inherit the re-run profiles.

**Decision (Pete, 2026-07-11): option (b) — EVERYTHING 4-CT, including the filter. TRUE zero-trace.** The whole `all_samples` analysis — isoform universe, normalization, reference, floor, partner — is AT/DD/FB/MV only. No DD_ALI/DO anywhere; no 6-CT numbers, no "all sequenced libraries" wording, no sensitivity note in any prose/legend/Methods/figure. 6-CT provenance lives ONLY in this doc + git (squashed at merge). SF28 shows 4-CT sample counts with **no** 36→26 drop.

**Implementation (source edit, `01_prepare_data_mashr.R`):** drop DD_ALI + DO at the **source** — restrict `sample_metadata`/`count_mat` to `ct %in% c("AT","DD","FB","MV")`. **Placement is critical (M-1): insert STRICTLY AFTER L103** (after both `count_mat` and `sample_metadata` are subset by the `exclude_samples` block) and before the `dge`(L119)/`cpm_mat`(L123) build — NOT into/before the L94 block, or the L96–98 `sum(excluded)!=length(exclude_samples)` guard aborts once DO-029T is already gone. The DO-029T exclusion then becomes redundant but harmless. Everything downstream (`calcNormFactors`/`cpm` on 26 libs, `dmso_cols`/`smg1i_cols`, the 5% + `filterByExpr` filters, per-CT lists, and `all_samples` L204/205) becomes 4-CT **by construction** — verified `dge`/`cpm_mat` build after the exclude block. **05t/05v: DROPPED (round 3)** — both are dead for the paper (05v feeds only the legacy Rmd); do NOT re-run them. **SF42 un-freeze (Pete 2026-07-11) goes via its REAL upstream `code/nmd_predictor_comparison/` (01→04)**, re-run after R2 — see the round-3 correction block below. This OVERRIDES the earlier "SF42 frozen" exemption; SF42 becomes a regenerated 4-CT figure.
**This supersedes** the Decision block's "DMSO basis = `dmso_samp[[ct]]` per profile" clause (per-CT unchanged; `all_samples` is now 4-CT throughout, and the isoform universe is 4-CT-normalized/filtered).

### New acceptance gates (4-CT; exact numbers from the data-layer run)
`expr_mat` **will change** (fewer samples → different normalization + filter → different isoform set), so ALL isopair numbers move — more than the earlier L204-only "~1,519" guess, which is now **void** (it ignored reference-reselection, C2-partner churn, AND the universe change). Pre-commit *direction + band*, not a point:
| Set | pre-floor (6-CT) | **4-CT floor (target)** | **4-CT MEASURED (2026-07-11)** |
|---|---|---|---|
| pop_BC (all_samples C2) | 3,009 | **~1,400–1,550** | **1,548** ✓ (top of band; 6-CT floor was 1,585) |
| n=190 GENCODE-restricted | 190 | **TBD (< 136)** | **130** ✓ |
| n=1,166 ref-AUG | 1,166 | **TBD (< 888)** | **819** ✓ |
| occult-PTC | 492 | **TBD (< 380)** | **348** ✓ |
| per-CT C2 (AT/DD/FB/MV) | 2,583–2,907 | all ≫ MIN_PAIRS=50 | **AT 1,432 · DD 1,476 · FB 1,340 · MV 1,466** — all ≫ 50 ✓ |

**Isoform universe (M-2 decomposition):** 6-CT filtered 123,790 → 4-CT filtered 95,623 (net −28,167).
Removed 28,278 = (a) 1,790 DD_ALI/DO-exclusive (0 counts in all 26 4-CT samples) + (b) 26,488
re-normalization/threshold churn (5% max-prop + `filterByExpr` flips from dropping 10 libs); added 111
(flipped in). Recomputed filter reproduces 95,623 exactly → decomposition faithful. Gate: **PASS** (no per-CT < 50; pop_BC in band). AWAITS Pete's confirmation before propagating to figures/report.
**Hard-halt** only if any per-CT C2 < MIN_PAIRS=50. **Advisory** (investigate, don't auto-halt): pop_BC expected ~1,400–1,550; flag if outside ~1,300–1,585. **Do NOT gate on "isoform drop explained by DD_ALI/DO-only" (M-2)** — dropping 8–10 libs also shifts TMM norm factors + `filterByExpr` + the 5% max-prop, flipping borderline 4-CT isoforms pass↔fail, so the retained-set delta is genuine filter/normalization churn, not decomposable as DD_ALI/DO-only. At R2 instead **decompose** the isoform-set delta into (a) DD_ALI/DO-exclusive isoforms removed, (b) borderline 4-CT isoforms flipped on the 5%/`filterByExpr` thresholds by re-normalization, (c) net pop_BC delta — and sanity-check each is plausible.

### Safeguard — 6-CT dependency scan (run at EVERY stage before trusting output)
Grep each stage's code for: `DD_ALI` · `DO_ALI` · `\bDO\b` · `Sample11|Sample15|Sample16|Sample21|Sample23` · `dmso_cols` · `smg1i_cols` · `all_samples` · `union(` · `Reduce(intersect` · `colnames(expr` · `rowMeans(.*dmso` · `treatment == "DMSO"` · `treatment == "Smg1i"` · hardcoded `\b18\b`/`\b36\b`/`\b26\b`/6-CT counts · **stale floor literals — BOTH post-floor** `1585` `136` `888` `380` `818` `70` `508` `345` `54` `82` `72` `69` **AND pre-floor (these ALSO move under b — C-1)** `3009`/`3,009` `8323`/`8,323` `190` `1166`/`1,166` `492` `2583`–`2907`. Scan `.py` **docstrings/titles + `data_export.R`/`.Rmd` cat-strings/comments + per-figure `README.md`** (m-3), not just `.R` assert lines. NB `figure3.../data_export.R` already prints a stale `pop_BC = 3,009` label in the shipped 6-CT state — direct proof the label scan (not just asserts) is required. Any hit → review. **Positive assertions** (`stopifnot(length(s)==13, !any(grepl("DD_ALI|_DO_", s)))`) into the `02` floor block, `code/nmd_predictor_comparison/01`, every figure `data_export.R` touching a basis (SF26 first), and `SF30/data_export.R` (`stopifnot(length(smg1i_cols)==4)` — m1). Stages: (1) `01`/`02`, (2) caches `03b`/`05r`/`05k`/`05k_b`, (3) `code/nmd_predictor_comparison/` (01→04) + Rmd, (4) every figure `data_export.R`, (5) verifiers.

### Deprecation of the 6-CT floor pass (Option A + squash)
Overwrite the 6-CT artifacts in place; `git mv` `05_final_report_gencode_scope_2026-07-10.Rmd` → `..._2026-07-11.Rmd` (+ .html). Regenerate figures/verifiers/map with 4-CT numbers. **Find/replace — ELIMINATE, don't keep (Pete 2026-07-11):** `git rm paper/section4_findreplace_2026-07-10_referencefloor.md` (never applied, 6-CT) and generate a fresh `paper/section4_findreplace_2026-07-11_4ct.md` with 4-CT numbers + 4-CT basis wording (Methods M-1 "all sequenced libraries" → "the four cell types"). **Filename bump surface (M3/m3 — widen beyond .html):** every `2026-07-10` ref in `figure5_dl_model/data_export.R`, `data_export_n1166.R`, `SF26/data_export.R`, `SF26/*.py`, `SF39/*.py`, `figures/lib/ggplot_style.py`, `paper/results_to_code_map.md`, and the three tracking docs (the 6-CT find/replace is **deleted**, not renamed). Squash-merge at Phase 8 → no 6-CT report file, no 6-CT commit in `main`. 2026-06-15 pre-floor report stays LEGACY (shipped state); 6-CT 07-10 report does NOT survive.

### Review findings folded in — ROUND 1 (2026-07-11 adversarial review)
- **C1** (isoform universe 6-CT) → resolved by decision (b): filter/normalize on 4-CT at source.
- **M1** — Phase-0 inventory below **mis-lists 05t as a Rmd dependency — it is NOT** (Rmd L47–49 load `ref_atg_analysis`/`utr5_features_refaug`/`utr5_features_all` = 05r/05k_b/05k; 05t feeds `05v` predictor-comparison only). Correct the inventory. (05t *basis* fix superseded by round-2 M-3 — no edit, just re-run.)
- **M2** — halt gate has a band + "churn expected" (further reframed by round-2 M-2 below).
- **M3** — scan tokens include floor literals + docstring/cat-string scanning (above).
- **m1** — `SF30/data_export.R` `stopifnot(length(smg1i_cols)==4)` (round 2 confirms `smg1i_cols` = the 4 mashr posterior-mean columns, an external basis-free object — the guard is correct, not a sample-count leak).
- **m2** — reconcile stray pre-floor `Control n=1166` at `figure4…/verify_pass1_factual.R:83` before R4.
- **m3** — filename-bump surface widened (above).

### ⚠ CORRECTION — ROUND 3 (2026-07-11): SF42 data source was MISIDENTIFIED
The "un-freeze SF42 via `05t→05v→SF42`" wording (leak table #5, Implementation ¶, R2/R3/R4/R6) is **factually wrong** and must be re-targeted:
- **SF42 reads `code/nmd_predictor_comparison/per_isoform_scores_2026.6.20.tsv` + `metrics_summary_2026.6.20.tsv`** (`SF42_ModelComparison/figure_s_model_comparison.py:49-50`), NOT `05v`'s `model_comparison.rds`. SF42 = NMDetective-B / NMDEP predictor comparison (`results_to_code_map.md:647`), NOT the 05v TD2/Ref/Combined comparison.
- **05v feeds only `05_final_report_mashr.Rmd` (LEGACY/superseded)** → 05t/05v are **dead for the paper**; DROP them from the re-scope (note `model_comparison.rds` as a non-paper cache; round-2 M-3 "re-run 05t" is moot).
- SF42's real upstream is `code/nmd_predictor_comparison/` (01→02→03→04), which reads `profiles_c2/c4_allsamples.rds` (`01_extract_our_isoforms.R:45-46`, floor+rescope-affected) + `ref_atg_analysis.rds` (05r) + a frozen model-global `predictions_all_atg500_stop500.tsv`. Its TSVs are dated `2026.6.20` — **pre-floor**, so SF42 is stale against BOTH the floor and the rescope.
- **To genuinely 4-CT SF42:** move `code/nmd_predictor_comparison/` from "Separate deliverables / not regenerated here" (L99) INTO scope; re-run it after R2; re-render SF42; update its hardcoded cohort numbers (n=1,166/561/293, 255/30/276) in the SF42 legend, `code/nmd_predictor_comparison/METHODS.md` + `README.md` + sub-Rmd; add those files to the scan surface.
- **CONFIRMED (Pete 2026-07-11):** bring `code/nmd_predictor_comparison/` into scope, regenerate SF42 4-CT; drop 05t/05v. R2/R3/R4/R6 + leak table + Implementation ¶ + L108 updated to match.
- Also (m-4, round 3): sequence **R5 (map) before R6** so the fresh 4-CT find/replace is built against the updated map; R6 should reference predictor-comparison cohort/Spearman numbers, not "05v model-comparison."

### Review findings folded in — ROUND 4 (2026-07-11, targeted pass on `code/nmd_predictor_comparison/`)
**Verdict: the DATA path is genuinely 4-CT with zero surviving 6-CT trace** — only `01_extract_our_isoforms.R` reads external data (02/03/04 read only 01's TSVs); its inputs are all class (a) R2-refreshed (`profiles_c2/c4_allsamples`, `ref_atg_analysis`) or (b) legitimately basis-independent: the gold-standard mashr is an **Apr-29 4-condition refit** (6-CT quarantined in `isocall_dge/old/mashr/`; the `2026.3.10` is a stale label, not the gen date), the model is the **frozen v5_4ct** (not retrained), and `cds.rds`/`structures.rds` are transcript-structure annotations invariant to CT scope. **No category-(c) 6-CT dependency.** Two premise corrections: **R2 does NOT refresh `cds.rds`/`structures.rds`** (Mar-dated core outputs, untouched — fine, they're basis-independent); **`paralog_genes.rds` is NOT an input** (only an H5 split label; SF42 filters to `test`, excluding `test_paralog`). Gaps are DOC-ONLY (folded into R2/R3 above): SF42 legend counts (CRITICAL), DATESTAMP+path sync (MAJOR), METHODS/README/Rmd literals (MAJOR), var-name/comment counts (MINOR), + assert post-merge model-prob coverage (MINOR). Do those and SF42 is genuinely 4-CT.

### Review findings folded in — ROUND 2 (2026-07-11, re-review of decision b)
Round 2 verified (b) is mechanically sound (`dge`/`cpm_mat` build after the exclude block; classification already 4-CT; model not retrained; 03b/05r/05k/05k_b inherit 4-CT free). Fixes folded in:
- **C-1** — scan tokens now include **pre-floor** literals (3009/8323/190/1166/492/2583–2907), which also move under (b); scan prose/cat-strings/comments/README, not just asserts. (`figure3/data_export.R` already prints a stale `pop_BC = 3,009`.)
- **M-1** — 01 source edit pinned **strictly after L103** (else the L96–98 `sum(excluded)` guard aborts). R1 updated.
- **M-2** — halt reframed: hard-halt only on per-CT < 50; pop_BC band advisory; **decompose** the isoform-set delta rather than gate on "DD_ALI/DO-only explains it" (re-normalization flips borderline 4-CT isoforms too). Gate section + R2 updated.
- **M-3** — SUPERSEDED by the round-3 correction: 05t/05v are dead/legacy → **DROPPED** (not re-run). SF42 un-freeze goes via `code/nmd_predictor_comparison/`, not 05t/05v.
- **M-4** — SF28 is a **full manual DOT rebuild** (~20 hardcoded floor literals) + fix the `N_C2` node mislabeled "before matching" + delete dead 36/10 vars. R3 updated. (The "no 36→26" is already satisfied — those vars aren't referenced in the DOT.)
- **m-1/m-2/m-3** — stale task-list item corrected; halt band = advisory vs hard split stated; per-figure `README.md` added to scan surface.

---

## ~~Acceptance gates (SUPERSEDED — see Amendment above)~~
| Set | pre-floor | ~~post-floor (25% all-iso, 6-CT)~~ |
|---|---|---|
| pop_BC (all_samples C2) | 3,009 | ~~**~1,585**~~ |
| n=190 GENCODE-restricted | 190 | ~~**~134**~~ |
| n=1,166 ref-AUG | 1,166 | ~~**~861**~~ |
| per-CT C2 (AT/DD/FB/MV) | 2,583–2,907 | 1,358–1,487 (all ≫ MIN_PAIRS=50) |

## Impacted-artifact inventory (Phase 0)
**Injection point:** `02_build_profiles_mashr.R` (post-`generatePairsExpression` filter; C2 inherits via gene_id join).

**Caches to rebuild:**
- `02` → `data_mashr/pairs_c{2,4}_*.rds`, `profiles_c{2,4}_*.rds` (pop_BC derives from these).
- `03b_rebuild_cache.R --force` → `analysis_cache/{div,ri,cooc,er,ptc,fw,fc,cps}_*.rds` (**--force required**; `cached_compute` returns stale otherwise). Rmd reads `fc_c2`,`fw_c2` (lines 798–799).
- Feature-cache producers the Rmd loads (lines 47–49): `05r_ref_atg_analysis.R` → `ref_atg_analysis.rds`; `05k_utr5_all_isoforms.R` → `utr5_features_all.rds`; `05k_b_utr5_refaug.R` → `utr5_features_refaug.rds`; `05t_ref_cds_features.R`. **Rerun OR prove each is inner-joined to floored pop_BC (superset-safe).**

**OUT of scope:** `04_productive_frameshift_precompute.R` / `04b` — outputs read only by the *deprecated* legacy `05_final_report_mashr.Rmd` + archive; canonical Rmd & figures never read them. Confirmed by grep.

**Report:** copy `05_final_report_gencode_scope_2026-06-15.Rmd` → `..._2026-07-10.Rmd`, edit, render.

**Figures (regenerate one-at-a-time; ordering matters):**
- Producers first: `figures/multipanel/figure5_dl_model/data_export.R` + `data_export_n1166.R` emit `gencode_all3_n190_isoforms.tsv` / `subset2_n1166_isoforms.tsv`.
- Then consumers of those TSVs: `SF30_PTCDistanceDoseResponse`, `SF39_PTCSubclassBranchSHAP`.
- Multipanel: `figure3_isopair_and_ptc` (A/B/C pop_BC; D/E/F n=190), `figure4_ptcneg_and_model` (A/B n=190; C/D n=1,166), `figure5_dl_model`.
- Supplements: SF24, SF25, SF26 (**export logic deleted with PairSetDescriptives — recover via `git show f6d96bc^:…`; rebuild on all-iso basis**), SF27, SF28 flowchart, SF29, SF30, SF31, SF32/SF35 (companion — edit together), SF33/SF34, SF39/SF40, SF41.

**Verifiers (repo-root `reproducibility/` + per-figure):**
- Central: `reproducibility/verify_pass7_new_rmd.R` (37 `expected=`), `reproducibility/verify_cross_check_new_rmd_vs_figures.R` (~12 `exp=` — reconcile vs map's "57"; line 22 hardcodes `..._2026-06-15.html` → update for new date).
- Per-figure: `figure3_isopair_and_ptc/verify_pass{1_factual,2_correctness,5_methods}.R`; `figure4_ptcneg_and_model/verify_pass4_reproducibility.R`. **New expected values independently recomputed (5-step step 2), NOT pasted from the new render.**

**Separate deliverables:** `code/nmd_atlas/export_atlas_data.R` — still OUT of scope (not regenerated). ~~`code/nmd_predictor_comparison/`~~ → **NOW IN SCOPE (round 3, Pete 2026-07-11):** re-run 01→04 after R2 to regenerate SF42 on 4-CT.

## Phase check-off
- [x] **Phase 0** — branch, tracking docs, inventory, before-state snapshot
- [x] **Phase 1** — floor in `02_build_profiles_mashr.R` (`REF_SHARE_FLOOR<-0.25`, all-iso denom)
- [x] **Phase 2** — caches rebuilt (02 → 03b --force → 05r → 05k_b; 05k/05t floor-independent, left stale); pop_BC=1,585 confirmed, no CT < MIN_PAIRS
- [x] **Phase 3** — Rmd date-bumped to 2026-07-10 + rendered (all floored numbers verified); figures done
      - [x] data layer: all `data_export.R` regenerated (guards 190→136, 1166→888, occult 492→380); see `REFERENCE_FLOOR_NUMBERS_DELTA.md`
      - [x] **Main figures 3, 4, 5** — rendered + visually inspected + CORRECT (data-driven; n's updated; model panels unchanged). Python = `/opt/homebrew/bin/python3` (system python3 lacks pandas).
      - [x] **Supplements — ALL 13 floor-affected DONE** (SF25–35, 39, 40). Full audit in `REFERENCE_FLOOR_SF_INVENTORY.md`: SF24/36/37/38/41 N/A (model-global/static), ~~SF42 left as-is~~ (SF42 UN-FROZEN at R3 → regenerate on 4-CT, per 2026-07-11 decision). Fixes covered: path-splits (SF32/33/34/35), hardcoded stats (SF33/34), layout-clips (SF30/31/40), rebuilt producers (SF25/26/27/28/29).
      - [x] Rmd 2026-07-10 rendered clean; ref-share now 67% (matches SF26 66.7%); repointed 2 removed composites to split SFs
- [x] **Phase 4** — all 6 verifiers PASS (252 checks); expecteds independently re-derived; caught+fixed Fisher one-sided/two-sided + pass2 N_BC scope. 5-step: steps 1-4 covered by suite+independent derivations; step 5 (document floor in METHODS) rolls into Phase 6/7
- [~] **Phase 5** — 2026-06-15 Rmd bannered LEGACY/superseded (done); SF26 supersession handled in rebuild
- [x] **Phase 6** — results_to_code_map.md: §4 claims 4.1-4.46 + scope table + verifiable-summary all floored; M11 documents the floor; 4.5 drift RESOLVED (70/75 was dominant-share mislabel -> 67/71); all 2026-06-15->2026-07-10 filename bumps
- [x] **Phase 7** — manuscript find/replace pairs drafted (paper/section4_findreplace_2026-07-10_referencefloor.md: 26 §4/legend + 5 Methods incl. floor documentation); SF legends updated; docx rebuilt. AWAITS Pete applying pairs to the Google Doc + verifying vs live Doc.
- [ ] **Phase 8** — one coherent commit; dual-push

---

## Phase check-off — 4-CT RE-SCOPE (2026-07-11); supersedes the 6-CT pass above
The phases above completed under the 6-CT `all_samples` basis and are now **superseded** (never shipped). Each phase re-opens against the 4-CT basis. **Run the 6-CT dependency scan (Amendment) at every phase.**

- [x] **R1 — Source edit** — `01`: restrict `sample_metadata`/`count_mat` to `ct %in% c("AT","DD","FB","MV")` **strictly AFTER L103** (M-1: NOT in/before the L94 exclude block — the L96–98 guard aborts). Confirm `dge`(L119)/`cpm_mat`(L123) build after → `dmso_cols`/`smg1i_cols`=13, `all_samples` 4-CT by construction (no L204/L205 edit). Scan `01`/`02` + `code/nmd_predictor_comparison/` for residual 6-CT deps (+ pre-floor literals). (05t/05v dropped — dead/legacy.)
- [ ] **R2 — Data layer + gate** — backup `data_mashr`; re-run `01`. **Gate (b):** assert `sample_metadata`/`expr_mat` have only 4 CTs (26 cols, no DD_ALI/DO). **Decompose the isoform-set delta (M-2)** into (a) DD_ALI/DO-exclusive isoforms removed, (b) borderline 4-CT isoforms flipped by re-normalization on the 5%/`filterByExpr` thresholds, (c) net — sanity-check each; do NOT gate on "(a) explains all." Re-run `02` (+ `stopifnot(length==13,…)` floor-block guard) → `03b --force` → `05r`/`05k`/`05k_b` (**05t/05v DROPPED — dead/legacy**). Then re-run **`code/nmd_predictor_comparison/` 01→04** on the refreshed profiles (bump `DATESTAMP` `2026.6.20`→run date; **assert post-merge `our_model_prob` coverage** — joins are `all.x=TRUE`, so new-cohort isoforms absent from the frozen v5_4ct H5 silently drop; don't trust the old 95%/114-dropped) → new 4-CT TSVs for SF42. Record pop_BC, n=190→?, n=1166→?, occult→?, per-CT. **Hard-halt** only if any per-CT C2 < 50; pop_BC ~1,400–1,550 expected, investigate if outside ~1,300–1,585 (advisory).
- [ ] **R3 — Report + figures** — `git mv` Rmd → `..._2026-07-11.Rmd`; re-render; strip ALL 6-CT basis wording → "four cell types". Regenerate figures ONE AT A TIME (Figs 3/4/5, SF25–35, SF39, SF40, **SF42 — UN-FROZEN (round-4 review: data path is genuinely 4-CT, gaps are DOC-ONLY). After re-running the pipeline: (i) bump `DATESTAMP` in `01-04` + sync the 4 hardcoded paths (`figure_s_model_comparison.py:49-50`, `nmd_predictor_comparison.Rmd:23-24`); (ii) re-derive SF42 **legend** counts (n=561/255/30/276) from the fresh `metrics_summary` `head-to-head:test:*` rows — CRITICAL, legend prose is frozen while the PNG auto-derives n; (iii) re-derive cohort literals in `METHODS.md`/`README.md`/`.Rmd`; (iv) rename `c2_n1166`/`c4_n1166` (`01:78-79`) + drop count literals from comments**) + 6-CT scan (pre+post-floor literals in docstrings/cat-strings/README) each `data_export.R`/`.py` + basis assertion. **SF28 = full manual DOT rebuild (M-4):** ~20 hand-typed floor literals in `build_flowchart.R` (L135/144/157-159/176/190-192/209/212 + filter counts) re-derived by hand; fix the `N_C2` node mislabeled "before matching" (it reads the *floored* profiles_c2); delete dead `N_AFTER_OUT`/`N_FILTER_SMP`/`N_DROPPED_SMP` (36/10) vars — render already shows 26 samples / no 36→26 node. Fix SF26 legend "all sequenced libraries" → "four cell types".
- [ ] **R4 — Verifiers** — expecteds INDEPENDENTLY re-derived (not pasted from render); resolve m2 (`verify_pass1_factual.R:83` stray 1166); update SF42 predictor-comparison expecteds/cohort (un-frozen); date-bump `2026-07-10` → 07-11 across ALL refs (m3 surface, not just .html); all pass.
- [ ] **R5 — Map** — results_to_code_map.md → 4-CT numbers + 07-11 filenames; M11 basis wording → 4-CT; correct the 05t-as-Rmd-dep misattribution (M1).
- [ ] **R6 — Manuscript** — `git rm paper/section4_findreplace_2026-07-10_referencefloor.md` (never applied, 6-CT — **eliminate, don't supersede**) and generate a fresh `paper/section4_findreplace_2026-07-11_4ct.md` with 4-CT numbers + basis wording (Methods M-1 "all sequenced libraries" → "the four cell types"). **Run R5 (map) BEFORE this (m-4)** so find/replace builds against the updated map. Include any SF42 predictor-comparison cohort/Spearman number changes if cited in §4/§5 prose.
- [ ] **R7 — Legend tasks** (after numbers settle): SF28 (center "Construct paired comparisons" text + genes-dropped count/reason from new floor delta), SF29 + all star-using figs significance-star key, SF30 cross-CT-averaging note + `stopifnot(len==4)` (m1) + NMD-susceptible-in-≥1-CT verification, "exon junction complex"→EJC + Abbreviations everywhere.
- [ ] **R8 — Squash-merge** to one 4-CT commit; NO 6-CT report/commit survives; dual-push on Pete's go.
