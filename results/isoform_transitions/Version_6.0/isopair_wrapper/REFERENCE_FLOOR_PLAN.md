# Reference-share floor — implementation plan & check-off

**Branch:** `reference-floor-25pct` · **Started:** 2026-07-10 · **Owner:** Pete (Castaldi)

## Decision (locked)
- **Family A** — reference *selection* unchanged: rank-1 by DMSO mean CPM within the **strict non-NMD** pool (`nmd_class[[ct]]$non_nmd`). Do NOT re-anchor to the gene dominant.
- **25% floor** — additionally drop any gene whose selected reference accounts for **< 25% of total gene expression**.
- **all-isoform denominator** — the 25% is `mean_dmso(reference) / Σ_{all isoforms of gene} mean_dmso`, DMSO basis = `dmso_samp[[ct]]` per profile. NB "all isoforms" = all isoforms **passing the 5% condition-stratified filter** in `01_prepare_data_mashr.R` (expr_mat is already filtered); Methods must not overstate this.
- **Rationale:** analyze only genes whose reference is *truly non-NMD* AND *truly dominant*. Displacement cases (true dominant is a "neither"-gap isoform) are dropped, not re-anchored. See `REFERENCE_FLOOR_RATIONALE.md`.

## Acceptance gates (predicted; halt if off)
| Set | pre-floor | post-floor (25% all-iso) |
|---|---|---|
| pop_BC (all_samples C2) | 3,009 | **~1,585** |
| n=190 GENCODE-restricted | 190 | **~134** |
| n=1,166 ref-AUG | 1,166 | **~861** |
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
- [ ] **Phase 3** — Rmd date-bump + render; figures one-at-a-time
      - [x] data layer: all `data_export.R` regenerated (guards 190→136, 1166→888, occult 492→380); see `REFERENCE_FLOOR_NUMBERS_DELTA.md`
      - [x] **Main figures 3, 4, 5** — rendered + visually inspected + CORRECT (data-driven; n's updated; model panels unchanged). Python = `/opt/homebrew/bin/python3` (system python3 lacks pandas).
      - [ ] **Supplements — need per-figure work (NOT batch):**
            - PATH SPLIT: `figure_sfNN.py` reads `SFNN_*/data/` (stale, 06-15) but `data_export.R` writes to shared dir (e.g. `TD2BiasEvidence/data/`, `CDSand3UTR_GENCODEonly/data/`, fresh 07-10). Confirmed SF33; check SF32/34/35/39. → fix figure DATA path (or export target) per figure, then re-render.
            - HARDCODED stat annotations: e.g. SF33 `figure_sf33.py:201` p-value `"1.2×10^-38"` → `2.9×10^-39`. Grep each SF for hardcoded n/p/median.
            - LAYOUT-CLIP errors (validator caught, new data ranges): SF30, SF31, SF40 — adjust top margin / ylim.
            - SF26: `data_export` deleted with PairSetDescriptives — recover via `git show f6d96bc^:…`, rebuild on all-iso basis.
            - ⚠ SF33/34/35/39 PNGs currently rendered from STALE data — do NOT use until path fixed + re-rendered + re-inspected.
      - [ ] date-bump Rmd → 2026-07-10, update hardcoded prose/DOT/captions, render
- [ ] **Phase 4** — update + independently re-derive all verifier expecteds; PASS; full 5-step
- [ ] **Phase 5** — deprecate 2026-06-15 Rmd (banner); SF26 supersession
- [ ] **Phase 6** — results_to_code_map.md 4.1–4.46 + M11 + filename/verifier-path bumps
- [ ] **Phase 7** — manuscript find/replace (grep Abstract + legends + supplement, not only §4); Methods (Isopair vignette + manuscript); docx rebuild
- [ ] **Phase 8** — one coherent commit; dual-push
