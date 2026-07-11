# Plan v4 — branch + linear workflow (executable)

**Status:** v4 — ready to execute after Plan-agent pass #3 found and the user approved 5 critical fixes.
**Supersedes:** v1 (in-place rewrite, rejected), v2 (figure-TSV consumption, rejected for workflow tangling), v3 (linear workflow but contained provenance errors and helper-classification errors).

## Architecture (settled across v2/v3/v4)

**Branch-don't-rewrite + linear workflow:**

```
primary RDS caches  →  new Rmd (loads caches, computes scopes,
  (isopair_wrapper       runs Isopair + analysis_functions.R, statistical tests,
   analysis output)       renders all tables + bound inline-R values)
                       →  inline-R values feed paper prose
                       →  Figure scripts independently compute the same scopes
                           from the same caches (downstream renderers)
```

The new Rmd is the canonical analysis report. The figure scripts continue as independent downstream renderers from the same upstream RDS caches. Both must produce identical numbers; drift is caught by Pass-1 verification (Phase 4).

## Design decisions (locked)

| ID | Decision |
|---|---|
| **RD1** | New Rmd filename: `05_final_report_gencode_scope_2026-06-15.Rmd` |
| **RD2** | Single file covering all of §3–§4 at current scopes |
| **RD3** | Rmd is canonical source of truth. Loads primary RDS caches; computes everything inline; no reading of figure-folder TSVs. |
| **RD4** | Legacy markers on old Rmds: (a) visible banner chunk + (c) YAML title flag |
| **RD5** | Iteration cadence: Plan-agent reviews continued through pass #3; v4 ready to execute |
| **RD6** | Shared libraries: **Isopair** (canonical workflow library) + **`isopair_wrapper/analysis_functions.R`** (workflow-internal helper, sourced by the legacy Rmd + figure3 panel_e_compute.R + the new Rmd). Inline-replicate anything else. |

## Critical fixes from Plan-agent pass #3 (applied in v4)

| # | Fix | Where |
|---|---|---|
| **C1** | `attribute_ptc_events()` and `attribute_3utr_splice()` are in `isopair_wrapper/analysis_functions.R` — NOT Isopair. The new Rmd `source()`s that file. | Setup chunk + all chunks that call these functions |
| **C2** | `figures/lib/mechanism_class.R` is figure-side. Inline-replicate the 50-line `mechanism_class()` + `mechanism_class_4()` helpers in the new Rmd, NOT source. | New Rmd helper chunk |
| **C3** | Replace `here::here()` with an explicit `NMD_ROOT` constant defined at top of setup chunk. The `nmd/` repo has no `.git` / `.Rproj` / `.here` sentinel and `here::here()` would resolve to Isopair's `DESCRIPTION` instead. | Setup chunk |
| **C4** | §1 numbers (3,067 / 3,107 / 2,906 nt; 70%; 75%) are NOT upstream Isopair constants — they're inline computations. The new Rmd recomputes them from `structures.rds` + `profiles_list`. | §1 |
| **C5** | Setup chunk asserts `packageVersion("Isopair") >= "X.Y.Z"` + checks that `utr5_features_refaug.rds` has the columns the Rmd binds against. | Setup chunk |

## Should-fix items from review #3

| # | Item | Where |
|---|---|---|
| **S1** | Phase 1 smoke is hybrid: `run_pandoc=FALSE` on the legacy Rmds + a separate 3-chunk minimal banner-test Rmd that knits through pandoc to confirm the HTML `<div>` renders. | Phase 1 |
| **S2** | Phase 4 verification: concrete `verify_pass7_new_rmd.R` skeleton with an expected-value manifest (`expected_values <- list(nmd_ptc_pos_rate_secA = 37.9, ...)`) with per-value tolerance (counts exact; percentages ±0.1; p-values same order of magnitude). | Phase 4 |
| **S3** | Move legacy verifier `verify_pass6_rmd_source.R` out of `figures/multipanel/figure3_isopair_and_ptc/` to a top-level `reproducibility/` folder and rename to `verify_legacy_rmd_reproducibility.R`. Workflow-coherence: a verifier of the Rmd doesn't belong in a figure folder. | Phase 4.5 |
| **S4** | Rollback policy in Phase 4: if new Rmd ≠ figure number, do not modify either until root cause identified. Do not paper over off-by-one discrepancies. | Phase 4 |
| **S5** | Phase 6 (follow-up task, not blocking the v4 execution): promote `build_section_A_scope()`, `build_section_C_scope()`, `compute_own_gencode_ptc()` to Isopair so the next scope refinement is a one-place change. | Phase 6 (deferred) |

## Nice-to-have items from review #3

| # | Item |
|---|---|
| **N1** | §4 (TD2 bias) has a §3-only scope dependency. It could be authored immediately after §3 scope is built, before §2 subsections, if a reviewer wants the supplement validated first. The default ordering keeps §2 before §4 for narrative flow, but the dependency is explicitly noted. |
| **N2** | §5 cumulative accounting is flagged as a **verification chokepoint** (not just a write-down): it's the natural place to spot internal inconsistency between §2 and §3 results. Inline-R values in §5 should be derived from the §2/§3 result objects, not re-computed. |
| **N3** | Realistic re-knit time: **8–15 min cold-start**. First debugging iterations slower. |
| **N4** | §6 cross-references the legacy Rmd content (banner-only is fine; no need to also list every legacy chunk by name). |

## What the new Rmd does

```
05_final_report_gencode_scope_2026-06-15.Rmd

  YAML header:
    title: "NMD analysis: current GENCODE-only scopes (2026-06-15)"
    output: html_document { toc, toc_depth: 3, toc_float, code_folding: hide }

  Setup chunk (include=FALSE):
    NMD_ROOT <- "/Users/petecastaldi/claude_projects/nmd"  # explicit, NOT here::here()
    ISOPAIR_WRAPPER <- file.path(NMD_ROOT, "results/isoform_transitions/Version_6.0/isopair_wrapper")
    DM <- file.path(ISOPAIR_WRAPPER, "data_mashr")

    library(Isopair)
    stopifnot(packageVersion("Isopair") >= "X.Y.Z")  # set to current
    source(file.path(ISOPAIR_WRAPPER, "analysis_functions.R"))
      # source: attribute_ptc_events(), attribute_3utr_splice(), helpers

    # Load primary RDS caches
    profiles_c2 <- readRDS(file.path(DM, "profiles_c2_allsamples.rds"))
    profiles_c4 <- readRDS(file.path(DM, "profiles_c4_allsamples.rds"))
    cds         <- as.data.table(readRDS(file.path(DM, "cds.rds")))
    structures  <- as.data.table(readRDS(file.path(DM, "structures.rds")))
    ref_atg     <- readRDS(file.path(DM, "analysis_cache/ref_atg_analysis.rds"))
    utr5_refaug <- readRDS(file.path(DM, "analysis_cache/utr5_features_refaug.rds"))
    expr_mat    <- readRDS(file.path(DM, "expression_data.rds"))
    nmd_class   <- readRDS(file.path(DM, "nmd_classification.rds"))

    # Schema validation (cheap insurance per C5)
    stopifnot(all(c("utr5_length", "longest_orf_nt", "isoform_features")
                  %in% c(names(utr5_refaug), names(utr5_refaug$isoform_features))))

  Helper chunk (include=FALSE):
    # mechanism_class() + mechanism_class_4() — inline-replicated from
    # figures/lib/mechanism_class.R (50 lines, no source). Kept inline per C2
    # so the Rmd is fully self-contained. If this helper ever drifts, the
    # two copies are simple enough to compare in one diff.
    mechanism_class <- function(...) { ... }
    mechanism_class_4 <- function(...) { ... }

  §1 Isopair construction and pair sets (inline computations per C4)
    - pop_BC = 3,009/3,009 (re-implement the Stage-2 gene-match filter inline)
    - Pair-set descriptive stats (median isoforms/gene, expression percentiles,
      transcript length comparisons) recomputed from structures.rds + expr_mat +
      profiles_list. Numbers (3,067 / 3,107 / 2,906 / 70% / 75%) are inline outputs.
    - Reference Figure 3 Panels A/B/C (PNG embedded via include_graphics)

  §2 PTC analysis at the GENCODE-only scope (n=190)
    Scope construction (inline, ~25 lines):
      pop_BC ∩ ENST-reference ∩ ENST-NMD-comp ∩ ENST-Control-comp
              ∩ coding-everywhere ∩ re-intersect on (gene, ref)
              → gencode_all3_c2, gencode_all3_c4 (190 / 190)
    §2a PTC determination
      - Own GENCODE stop via Isopair::genomicToTranscript + 50-nt rule
      - 72 / 4 PTC+; 18× enrichment, p<10⁻²⁰
      - Figure 3 Panel D (PNG embedded)
    §2b PTC-causing event attribution
      - source()ed attribute_ptc_events() + attribute_3utr_splice() on 72 PTC+
      - 66 direct + 8 same-stop = 74 attributions → 69 unique → 6 unresolved
      - Per-event-type Fisher table (Figure 3 Panel E, recomputed inline)
      - Mechanism breakdown 38/23/8 (Figure 3 Panel F)
    §2c 5'UTR length + longest 5'UTR ORF at Section A
      - Look up by comparator_isoform_id in utr5_refaug$isoform_features
      - Pairwise Wilcoxon + HL + Cliff's δ (functions defined inline or in
        analysis_functions.R; verify before relying on it)
      - Figure 4 Panels A/B (PNG embedded)
    §2d 3'UTR / CDS sanity at n=190
      - CDS length from cds$cds_stop - cds$cds_start
      - 3'UTR translation-based + non-PTC-stop based from ref_atg fields
      - CDSand3UTR_GENCODEonly supplement (PNG embedded)

  §3 Larger-scope replication (Section C, n=1,166)
    Scope construction (inline):
      pop_BC ∩ ENST-reference ∩ ref-AUG-traceable categories (each side)
              ∩ re-intersect on (gene, ref) AFTER category filter
              → refaug_matched_c2, refaug_matched_c4 (1,166 / 1,166)
    Filter cascade documented inline:
      ENST-ref pop_BC: 2,098 / 2,098
      → ref-AUG-traceable: 1,659 / 1,286
      → re-intersect: 1,166 / 1,166
    Classification via inline mechanism_class_4 → 1,050 PTC+ / 113 PTC− / 3 Ref-AUG-absent
    §3a 5'UTR length + longest 5'UTR ORF at Section C
      - utr5_refaug lookups; same pairwise tests as §2c
      - Figure 4 Panels C/D (PNG embedded)
    §3b PTC rate (1,050/1,166 = 90%)

  §4 TD2 bias evidence (TD2BiasEvidence supplement)
    [N1: this section has a §3-only scope dependency. Can be authored
     immediately after §3 scope build if convenient.]
    Scope:
      refaug_matched_c2 ∩ TD2 CDS available ∩ ref-AUG ORF computable ∩ ref AUG exonic
        → n=1,166 broad
      ∩ effectively_ptc ∩ original_ptc == FALSE
        → n=492 occult-PTC
    §4a TD2 vs ref-AUG ORF length
    §4b Paired Kozak PWM (Isopair::scoreKozakPWM)
    §4c TD2 ATG position vs ref AUG (Isopair::genomicToTranscript)
    TD2BiasEvidence supplement composite (PNG embedded)

  §5 Cumulative accounting [N2: verification chokepoint]
    - Section A: 72 PTC+ + 118 PTC− = 190 (recompute from §2 result objects)
    - Section C: 1,050 PTC+ + 113 PTC− + 3 Ref-AUG-absent = 1,166 (from §3)
    - Why Section A 38% vs Section C 90% — classifier methodology

  §6 What is intentionally NOT in this Rmd (legacy banner pointers)
    - Reclassification two-pass: 05_final_report_mashr.Rmd
    - Per-CT PTC stratification at old scope: 05_final_report_mashr.Rmd
    - "85% combined" framework: 05_final_report_mashr.Rmd
    - "n=1,904 / 22,335" GENCODE-only 3'UTR outside-set: 05_final_report_mashr.Rmd
    - ORF / model analyses: 06_orf_analysis_mashr.Rmd (also bannered)
```

## Open content questions — resolved by branch architecture

| Original Q | Resolution under v4 |
|---|---|
| Q1 — per-CT PTC table | Not in new Rmd; legacy keeps it |
| Q2 — TD2 length chunk | Compute inline (no TSV reading) |
| Q3 — PTC-attribution-on-reclassified | Not in new Rmd; legacy keeps it |
| Q4 — "remaining ~15%" subgroup | Not in new Rmd |
| Q5 — `06_orf_analysis_mashr.Rmd` §1 | Legacy banner only; no new ORF Rmd |
| Q6 — naming | Distinct: `gencode_all3_*`, `refaug_matched_*` for new |
| Q7 — cross-CT | Not in new Rmd; legacy keeps it |

## Phases

### Phase 1 — Legacy banners on old Rmds (hybrid smoke test per S1)

1. Edit `05_final_report_mashr.Rmd`:
   - YAML title: `"[LEGACY 2026-06-13] NMD as Post-Transcriptional Reading Frame Selection"`
   - YAML subtitle: `"Superseded by 05_final_report_gencode_scope_2026-06-15.Rmd as of 2026-06-15. Retained for sensitivity / reproducibility of prior submission."`
   - Insert banner chunk after setup chunk (visible HTML `<div>` with warning style).
2. Same for `06_orf_analysis_mashr.Rmd`.
3. **Smoke test (hybrid):**
   - `rmarkdown::render("05_final_report_mashr.Rmd", run_pandoc = FALSE)` — confirms R chunks still parse + run.
   - Build a 3-chunk minimal test Rmd containing JUST the YAML + banner chunk + a placeholder body chunk; knit it through pandoc to confirm the HTML `<div>` renders correctly.
   - **Do not attempt a full pandoc knit of the legacy Rmds** — they may fail for unrelated upstream-schema-drift reasons.

### Phase 2 — Skeleton new Rmd

1. Write `05_final_report_gencode_scope_2026-06-15.Rmd` with:
   - YAML header.
   - Setup chunk loading primary RDS caches + sourcing `analysis_functions.R` + asserting Isopair version + schema-checking `utr5_features_refaug.rds`.
   - Helper chunk with inline `mechanism_class()` / `mechanism_class_4()` (per C2).
   - Section headers + subsection placeholders matching the outline.
   - One inline-R value bound per section to confirm setup works.
2. **Pass-knit gate:** Rmd must knit without errors. All caches load, all libraries resolve, all schema checks pass.

### Phase 3 — Fill section by section

For each subsection:
- Add inline scope-construction code with comments documenting each filter.
- Add the statistical computations (Fisher, Wilcoxon, attribution).
- Bind inline-R values.
- Embed the relevant figure PNG via `include_graphics(file.path(NMD_ROOT, "figures/..."))`.
- Re-knit; spot-check against the corresponding figure number.

Execution order (build subset objects first, then dependent analyses):
1. §1 (pop_BC) — foundation
2. §2 scope (gencode_all3_c2/c4) — Section A subset
3. §3 scope (refaug_matched_c2/c4) — Section C subset (also §4's foundation per N1)
4. §2 subsections 2a → 2b → 2c → 2d
5. §3 subsections 3a → 3b
6. §4 subsections 4a → 4b → 4c
7. §5 cumulative accounting [verification chokepoint per N2]
8. §6 cross-reference list

### Phase 4 — Verification (concrete per S2 + rollback policy per S4)

1. Re-knit `05_final_report_gencode_scope_2026-06-15.Rmd` end-to-end.
2. Write `verify_pass7_new_rmd.R`:

```r
# verify_pass7_new_rmd.R
# Runs after the new Rmd knits; validates inline-R values against canonical figures.

expected_values <- list(
  # Section A
  nmd_ptc_pos_n_secA      = list(value = 72L,    tol = 0L),
  ctrl_ptc_pos_n_secA     = list(value = 4L,     tol = 0L),
  nmd_ptc_pos_rate_secA   = list(value = 37.9,   tol = 0.1),
  ctrl_ptc_pos_rate_secA  = list(value = 2.1,    tol = 0.1),
  fold_enrichment_secA    = list(value = 18,     tol = 1),
  attributed_events_n     = list(value = 69L,    tol = 0L),
  mechanism_fs_n          = list(value = 38L,    tol = 0L),
  mechanism_inframe_n     = list(value = 23L,    tol = 0L),
  mechanism_utr3_n        = list(value = 8L,     tol = 0L),

  # Section C
  ptc_pos_n_secC          = list(value = 1050L,  tol = 0L),
  ptc_neg_n_secC          = list(value = 113L,   tol = 0L),
  ctrl_n_secC             = list(value = 1166L,  tol = 0L),
  ptc_pos_rate_secC       = list(value = 90.1,   tol = 0.1),

  # TD2 bias supplement (broad scope)
  td2_orf_ratio_broad     = list(value = 4.69,   tol = 0.05),
  td2_kozak_p_broad       = list(value = 1.2e-38,tol_order = 1),
  ref_kozak_higher_broad  = list(value = 446L,   tol = 0L),

  # TD2 bias supplement (occult-PTC)
  td2_orf_ratio_occult    = list(value = 7.73,   tol = 0.05),
  td2_kozak_p_occult      = list(value = 2.6e-38,tol_order = 1),
  td2_downstream_pct_occult = list(value = 99.0, tol = 0.1)
)

# Source: figure3_composite_legend.md, figure4_composite_legend.md,
#         RATIONALE.md, TD2BiasEvidence/README.md.
```

3. The verifier extracts inline-R-bound values from the knitted .R (`knitr::purl()`) and validates against `expected_values`. Counts must match exactly; percentages within tolerance; p-values within order of magnitude.
4. **Rollback policy (S4):** if the new Rmd produces a value outside tolerance, DO NOT modify either the Rmd or the figure script until the root cause is identified. Off-by-one discrepancies always have a cause — either a filter-cascade bug, a deduplication bug, or upstream-cache drift. Investigate before papering over.

### Phase 4.5 — Verify-script housekeeping (S3)

- Rename `figures/multipanel/figure3_isopair_and_ptc/verify_pass6_rmd_source.R` to `reproducibility/verify_legacy_rmd_reproducibility.R` and move it out of the figure folder.
- Add deprecation header explaining it asserts numbers from the prior submission scope; reader should run the new Rmd + `verify_pass7_new_rmd.R` instead.
- Create `verify_pass7_new_rmd.R` in the same `reproducibility/` folder.

### Phase 5 — Cross-references + workflow doc

- Add cross-reference from `figure3/RATIONALE.md` and `figure4/RATIONALE.md` pointing at the new Rmd.
- Add a short section to `figures/README.md` titled "Workflow path":
  - The Rmd is the canonical analysis report.
  - Figure scripts are downstream renderers from the same upstream RDS caches.
  - Drift is verified by `verify_pass7_new_rmd.R`.

### Phase 6 — Promote scope construction to Isopair (DEFERRED, S5)

- Not blocking v4 execution.
- After v4 lands and the new Rmd is stable, propose adding `build_section_A_scope()`, `build_section_C_scope()`, `compute_own_gencode_ptc()`, `mechanism_class()`, `mechanism_class_4()` to Isopair.
- Once in Isopair, both the Rmd and figure scripts can `library(Isopair)` and call these directly — eliminating the inline-cascade duplication permanently.

## Risk assessment (v4)

| Risk | v4 status |
|---|---|
| Cascade breakage from `coding_pairs_c2` changes | NONE (old object untouched) |
| Workflow tangling between Rmd and figure scripts | NONE (Rmd doesn't read figure TSVs) |
| Wrong helper sourced (figure-side leaking into Rmd) | Mitigated by C2 (inline-replicate mechanism_class.R) |
| Wrong path resolution | Mitigated by C3 (explicit NMD_ROOT) |
| `attribute_ptc_events` provenance error | Fixed by C1 (source analysis_functions.R) |
| Silent number drift between Rmd and figures | LOW (both compute from same upstream RDS via same functions; verified by Phase 4) |
| Replication of scope-construction code | LOW-MEDIUM (~100-150 lines per Rmd; managed for now; promoted to Isopair in Phase 6) |
| Re-knit time | ~8–15 min cold-start (N3) |
| Schema change in upstream RDS | Mitigated by C5 (setup-chunk schema check) |
| Rollback policy violation | Mitigated by S4 (no off-by-one paper-over) |

## Time estimate

- Phase 1: 2 hours (banner addition + hybrid smoke test).
- Phase 2: 1 day (skeleton, all caches loading, all schema checks passing).
- Phase 3: 1.5 days (section-by-section fill + spot-check at each subsection).
- Phase 4 + 4.5: 0.5 day (verifier + housekeeping).
- Phase 5: 0.5 day (cross-references + workflow doc).
- **Total: 3–4 days of focused work.**

## Provenance

- v4 written 2026-06-15 by Claude after Plan-agent review pass #3.
- Critical fixes C1–C5 + Should-fix S1–S5 + Nice-to-have N1–N4 applied.
- Plan-agent's verdict: needs 5 specific revisions, then ready to execute. Revisions are in v4.
- Pete confirmed: write v4 with all critical fixes.
- Ready to execute. Phase 1 can start immediately upon Pete's go-ahead.
