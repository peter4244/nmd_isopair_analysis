# Reference-floor — Phase 4 verification record (2026-07-10)

## Automated verifier suite — 252 checks, all PASS
| Verifier | Result |
|---|---|
| `reproducibility/verify_pass7_new_rmd.R` | 37 / 37 |
| `reproducibility/verify_cross_check_new_rmd_vs_figures.R` | 57 / 57 |
| `figure3_isopair_and_ptc/verify_pass1_factual.R` | 57 / 57 |
| `figure3_isopair_and_ptc/verify_pass2_correctness.R` | 38 / 38 |
| `figure3_isopair_and_ptc/verify_pass5_methods.R` | 23 / 23 |
| `figure4_ptcneg_and_model/verify_pass4_reproducibility.R` | 40 / 0 |

Expected values were **independently re-derived** (not pasted from the render). The process caught
**two real errors**, both fixed:
1. Fisher p was one-sided (2.58e-18); the report + manifest use **two-sided → 5.15e-18** (confirmed:
   pre-floor 2×2 two-sided = 1.88e-20 matches the old manifest exactly). Fixed the DOT flowchart p.
2. `verify_pass2` hardcoded `N_BC = 3009` (pre-floor) as the 2×2 denominator → wrong OR; the panel
   TSVs were correct. Updated `N_BC→1585`, `N_DEF→136`, `CTRL_TOTAL_FULL→297`.

## Pete's 5-step scientific verification (independent adversarial pass)
An independent verifier recomputed every headline number from the raw `.rds` caches with its own code:

- **Step 1 — Factual:** pop_BC 1,585; all-3-ENST+CDS 136; ref-AUG 888; occult 380 (effectively_ptc ∩
  original_ptc==FALSE); PTC 54/1 of 136 (full independent strand-aware stop→EJC recompute); tx-length
  2,958/2,991/2,808 (inclusive end−start+1). **All match.**
- **Step 2 — Correctness:** 2×2 [[54,82],[1,135]] → OR 88.11, two-sided p 5.155e-18, fold 54×. Direction
  biologically sensible (NMD arm PTC-enriched; Control ~0 by construction). **No floor-induced selection
  artifact** — the floor drops on a gene-level property (reference share) and removes both arms together,
  so it cannot inflate the paired NMD-vs-Control asymmetry; reference share is mechanistically unrelated
  to comparator PTC status.
- **Step 3 — Documentation:** report internally consistent; the one-sided 2.58e-18 appears nowhere;
  ref-share 67% (complement-of-NMD, 4-CT-mean) ≈ SF26 66.7% (all-iso, all_samples). No number two ways.
- **Step 4 — Reproducibility:** every §4 figure/number traces to a committed `data_export.R` / script.
- **Step 5 — METHODS:** ⚠ **the floor is NOT yet documented** in the Isopair vignette, manuscript Methods,
  or the canonical Rmd. **Action for Phase 6/7.** Precise definition to document: a gene is retained iff
  `mean_DMSO(reference) / Σ_{i ∈ structures∩expressed isoforms of gene} mean_DMSO(i) ≥ 0.25`, DMSO basis =
  all_samples (18 DMSO libraries). NB the denominator is the structures-mapped, 5%-condition-filtered
  isoform set — not literally every annotated isoform (a naive all-annotated basis puts ~5 genes at ~24%).

## Outcome
All load-bearing numbers verified correct. Remaining: document the floor in METHODS (Phase 6/7).
