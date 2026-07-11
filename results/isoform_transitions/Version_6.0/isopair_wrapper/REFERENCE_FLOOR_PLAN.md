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
| 5 | 05t reference CDS features | `05t_ref_cds_features.R:78` rebuilds 6-CT DMSO basis from `sample_metadata` | predictor-comparison reference (via `05v`) |

**Clean (verified by review):** `all_samples` NMD *classification* (`01` L198–203, union/intersect over AT/DD/FB/MV) is already 4-CT; `buildProfiles` is structural (no per-sample leak); per-CT profiles never leaked; `05r`/`05k`/`05k_b` do **not** reconstruct a sample basis (structural / classification-keyed → inherit 4-CT free); pop_BC producers (`figure5_dl_model/data_export*.R`) and SF39/40 are basis-free and inherit the re-run profiles.

**Decision (Pete, 2026-07-11): option (b) — EVERYTHING 4-CT, including the filter. TRUE zero-trace.** The whole `all_samples` analysis — isoform universe, normalization, reference, floor, partner — is AT/DD/FB/MV only. No DD_ALI/DO anywhere; no 6-CT numbers, no "all sequenced libraries" wording, no sensitivity note in any prose/legend/Methods/figure. 6-CT provenance lives ONLY in this doc + git (squashed at merge). SF28 shows 4-CT sample counts with **no** 36→26 drop.

**Implementation (source edit, `01_prepare_data_mashr.R`):** drop DD_ALI + DO at the **source** — restrict `sample_metadata`/`count_mat` to `ct %in% c("AT","DD","FB","MV")` right after parsing (~L100, folding into/next to the existing `exclude_samples` block L94). Everything downstream (cpm/`calcNormFactors`, `dmso_cols`/`smg1i_cols`, the 5% + `filterByExpr` filters, per-CT lists, and `all_samples` at L204/205) then becomes 4-CT **by construction** — L204/L205 need no separate edit once `dmso_cols`/`smg1i_cols` are 13 libs. Also fix **`05t:78`** → `& ct %in% main_cts` (or confirm 05v/predictor-comparison is out-of-scope and freeze with a documented caveat — Pete to confirm at R6).
**This supersedes** the Decision block's "DMSO basis = `dmso_samp[[ct]]` per profile" clause (per-CT unchanged; `all_samples` is now 4-CT throughout, and the isoform universe is 4-CT-normalized/filtered).

### New acceptance gates (4-CT; exact numbers from the data-layer run)
`expr_mat` **will change** (fewer samples → different normalization + filter → different isoform set), so ALL isopair numbers move — more than the earlier L204-only "~1,519" guess, which is now **void** (it ignored reference-reselection, C2-partner churn, AND the universe change). Pre-commit *direction + band*, not a point:
| Set | pre-floor (6-CT) | **4-CT floor (target)** |
|---|---|---|
| pop_BC (all_samples C2) | 3,009 | **modest decrease from 1,585; expect ~1,400–1,550 — TBD, churn both ways is EXPECTED not a bug** |
| n=190 GENCODE-restricted | 190 | **TBD (< 136)** |
| n=1,166 ref-AUG | 1,166 | **TBD (< 888)** |
| occult-PTC | 492 | **TBD (< 380)** |
| per-CT C2 (AT/DD/FB/MV) | 2,583–2,907 | all ≫ MIN_PAIRS=50 (**halt if any < 50**) |
**Halt** if: any per-CT < 50; pop_BC outside ~1,300–1,585; or the isoform-count drop is not explained by DD_ALI/DO-only isoforms.

### Safeguard — 6-CT dependency scan (run at EVERY stage before trusting output)
Grep each stage's code for: `DD_ALI` · `DO_ALI` · `\bDO\b` · `Sample11|Sample15|Sample16|Sample21|Sample23` · `dmso_cols` · `smg1i_cols` · `all_samples` · `union(` · `Reduce(intersect` · `colnames(expr` · `rowMeans(.*dmso` · `treatment == "DMSO"` · `treatment == "Smg1i"` · hardcoded `\b18\b`/`\b36\b`/`\b26\b`/6-CT counts · **and the stale 6-CT floor literals** `1585` `136` `888` `380` `818` `70` `508` `345` `54` `82` `72` `69`. Scan `.py` **docstrings/titles + `data_export.R` cat-strings**, not just `.R` assert lines (M3). Any hit → review. **Positive assertions** (`stopifnot(length(s)==13, !any(grepl("DD_ALI|_DO_", s)))`) into the `02` floor block, `05t`, every figure `data_export.R` touching a basis (SF26 first), and `SF30/data_export.R` (`stopifnot(length(smg1i_cols)==4)` — m1). Stages: (1) `01`/`02`, (2) caches `03b`/`05r`/`05k`/`05k_b`/**`05t`**, (3) Rmd, (4) every figure `data_export.R`, (5) verifiers.

### Deprecation of the 6-CT floor pass (Option A + squash)
Overwrite the 6-CT artifacts in place; `git mv` `05_final_report_gencode_scope_2026-07-10.Rmd` → `..._2026-07-11.Rmd` (+ .html). Regenerate figures/verifiers/map/find-replace with 4-CT numbers. **Filename bump surface (M3/m3 — widen beyond .html):** every `2026-07-10` ref in `figure5_dl_model/data_export.R`, `data_export_n1166.R`, `SF26/data_export.R`, `SF26/*.py`, `SF39/*.py`, `figures/lib/ggplot_style.py`, `paper/results_to_code_map.md`, `paper/section4_findreplace_2026-07-10_referencefloor.md`, and the three tracking docs. Regenerate `section4_findreplace_*` (never applied) with 4-CT numbers + 4-CT basis wording (M-1 "all sequenced libraries" → "the four cell types"). Squash-merge at Phase 8 → no 6-CT report file, no 6-CT commit in `main`. 2026-06-15 pre-floor report stays LEGACY (shipped state); 6-CT 07-10 report does NOT survive.

### Review findings folded in (2026-07-11 adversarial review)
- **C1** (isoform universe 6-CT) → resolved by decision (b): filter/normalize on 4-CT at source.
- **M1** — `05t:78` reconstructs 6-CT DMSO; and the Phase-0 inventory below **mis-lists 05t as a Rmd dependency — it is NOT** (Rmd L47–49 load `ref_atg_analysis`/`utr5_features_refaug`/`utr5_features_all` = 05r/05k_b/05k; 05t feeds `05v` predictor-comparison only). Fix 05t basis + correct the inventory.
- **M2** — halt gate now has a band + explicit "churn expected" note (above).
- **M3** — scan tokens now include floor literals + docstring/cat-string scanning (above).
- **m1** — `SF30/data_export.R` assert exactly 4 Smg1i cols (folded into R7 note).
- **m2** — reconcile stray pre-floor `Control n=1166` at `figure4…/verify_pass1_factual.R:83` (intentional intermediate vs stale) before R4.
- **m3** — filename-bump surface widened (above).

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

**Separate deliverables (flag to Pete; not regenerated here):** `code/nmd_atlas/export_atlas_data.R`, `code/nmd_predictor_comparison/01_extract_our_isoforms.R`.

## Phase check-off
- [x] **Phase 0** — branch, tracking docs, inventory, before-state snapshot
- [x] **Phase 1** — floor in `02_build_profiles_mashr.R` (`REF_SHARE_FLOOR<-0.25`, all-iso denom)
- [x] **Phase 2** — caches rebuilt (02 → 03b --force → 05r → 05k_b; 05k/05t floor-independent, left stale); pop_BC=1,585 confirmed, no CT < MIN_PAIRS
- [x] **Phase 3** — Rmd date-bumped to 2026-07-10 + rendered (all floored numbers verified); figures done
      - [x] data layer: all `data_export.R` regenerated (guards 190→136, 1166→888, occult 492→380); see `REFERENCE_FLOOR_NUMBERS_DELTA.md`
      - [x] **Main figures 3, 4, 5** — rendered + visually inspected + CORRECT (data-driven; n's updated; model panels unchanged). Python = `/opt/homebrew/bin/python3` (system python3 lacks pandas).
      - [x] **Supplements — ALL 13 floor-affected DONE** (SF25–35, 39, 40). Full audit in `REFERENCE_FLOOR_SF_INVENTORY.md`: SF24/36/37/38/41 N/A (model-global/static), SF42 left as-is (frozen June comparison). Fixes covered: path-splits (SF32/33/34/35), hardcoded stats (SF33/34), layout-clips (SF30/31/40), rebuilt producers (SF25/26/27/28/29).
      - [x] Rmd 2026-07-10 rendered clean; ref-share now 67% (matches SF26 66.7%); repointed 2 removed composites to split SFs
- [x] **Phase 4** — all 6 verifiers PASS (252 checks); expecteds independently re-derived; caught+fixed Fisher one-sided/two-sided + pass2 N_BC scope. 5-step: steps 1-4 covered by suite+independent derivations; step 5 (document floor in METHODS) rolls into Phase 6/7
- [~] **Phase 5** — 2026-06-15 Rmd bannered LEGACY/superseded (done); SF26 supersession handled in rebuild
- [x] **Phase 6** — results_to_code_map.md: §4 claims 4.1-4.46 + scope table + verifiable-summary all floored; M11 documents the floor; 4.5 drift RESOLVED (70/75 was dominant-share mislabel -> 67/71); all 2026-06-15->2026-07-10 filename bumps
- [x] **Phase 7** — manuscript find/replace pairs drafted (paper/section4_findreplace_2026-07-10_referencefloor.md: 26 §4/legend + 5 Methods incl. floor documentation); SF legends updated; docx rebuilt. AWAITS Pete applying pairs to the Google Doc + verifying vs live Doc.
- [ ] **Phase 8** — one coherent commit; dual-push

---

## Phase check-off — 4-CT RE-SCOPE (2026-07-11); supersedes the 6-CT pass above
The phases above completed under the 6-CT `all_samples` basis and are now **superseded** (never shipped). Each phase re-opens against the 4-CT basis. **Run the 6-CT dependency scan (Amendment) at every phase.**

- [ ] **R1 — Source edit** — `01`: restrict `sample_metadata`/`count_mat` to `ct %in% c("AT","DD","FB","MV")` at source (~L100, with `exclude_samples`). Confirm this makes `dmso_cols`/`smg1i_cols` = 13 each and L204/L205 4-CT by construction (no separate edit). Also `05t:78` → 4-CT (or freeze per R6). Scan `01`/`02`/`05t` for residual 6-CT deps.
- [ ] **R2 — Data layer + HALT gate** — backup current `data_mashr`; re-run `01`. **Gate (redefined for b — `expr_mat` WILL change):** assert `sample_metadata`/`expr_mat` have only 4 CTs (26 cols, no DD_ALI/DO); the isoform-count drop is fully explained by DD_ALI/DO-only isoforms; nothing non-`all_samples`/non-`expr_mat` shifts unexpectedly. Re-run `02` (+ `stopifnot(length==13,…)` guard in floor block) → `03b --force` → `05r`/`05k`/`05k_b`/`05t`. Record exact pop_BC, n=190→?, n=1166→?, occult→?, per-CT. **Halt if outside the Amendment's band** (pop_BC ~1,300–1,585, per-CT ≥ 50, drop explained). Churn in/out is expected — not a bug.
- [ ] **R3 — Report + figures** — `git mv` Rmd → `..._2026-07-11.Rmd`; re-render; strip ALL 6-CT basis wording → "four cell types". Regenerate figures ONE AT A TIME (Figs 3/4/5, SF25–35, **SF28 flowchart: 4-CT counts, NO 36→26**, SF39, SF40) + 6-CT scan (incl. floor-literal docstrings) each `data_export.R`/`.py` + add basis assertion. Fix SF26 legend "all sequenced libraries" → "four cell types".
- [ ] **R4 — Verifiers** — expecteds INDEPENDENTLY re-derived (not pasted from render); resolve m2 (`verify_pass1_factual.R:83` stray 1166); date-bump `2026-07-10` → 07-11 across ALL refs (m3 surface, not just .html); all pass.
- [ ] **R5 — Map** — results_to_code_map.md → 4-CT numbers + 07-11 filenames; M11 basis wording → 4-CT; correct the 05t-as-Rmd-dep misattribution (M1).
- [ ] **R6 — Manuscript** — regenerate `section4_findreplace_*` with 4-CT numbers + basis wording (M-1 "all sequenced libraries" → "the four cell types"); supersede the never-applied 6-CT draft. **Confirm with Pete: is 05v/predictor-comparison in manuscript scope?** (drives whether 05t must be 4-CT or frozen-with-caveat).
- [ ] **R7 — Legend tasks** (after numbers settle): SF28 (center "Construct paired comparisons" text + genes-dropped count/reason from new floor delta), SF29 + all star-using figs significance-star key, SF30 cross-CT-averaging note + `stopifnot(len==4)` (m1) + NMD-susceptible-in-≥1-CT verification, "exon junction complex"→EJC + Abbreviations everywhere.
- [ ] **R8 — Squash-merge** to one 4-CT commit; NO 6-CT report/commit survives; dual-push on Pete's go.
