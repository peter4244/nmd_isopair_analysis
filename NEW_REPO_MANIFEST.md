# `nmd_lung_longread_2026` — file manifest

**Drafted 2026-07-24. Decisions settled by Pete the same day — see §3.** Repo name:
**`nmd_lung_longread_2026`**. Nothing created or moved yet.

**Principle (Pete):** the repo contains the code that produced the paper's results — **and not
more**. Project history, exploratory lineage, and internal process documents are out of scope,
on the condition that the reproducibility and correctness work is thorough.

**Sources:** `nmd_isopair_analysis` (561 tracked files) + `NMD_orf_model_v5_4ct` v2.0.0 (62
active files). Curation follows `KEEP_LIST.md` §§A–I.

---

## 1. Proposed layout (~385 files)

```
<new-repo>/
├── README.md                      what this is · scope · how to cite            [new]
├── REPRODUCTION.md                the runbook: deposit → every number/figure    [new]
├── ENVIRONMENT.md                 R 4.5.2 / Bioc 3.22 / Python 3.14.4 · seed 42 [from root]
├── config/paths.yml               single source of input roots                  [new]
│
├── analysis/
│   ├── upstream/          50   §1–§3 — Yul's DGE/DIE/mashr/landscape/GSEA/Tan
│   ├── isopair/           32   §4 — isopair_wrapper (26) + scripts/core (6)
│   └── predictor_comparison/ 16 §5 — NMDetective-B / NMDEP benchmark (SF43)
│
├── model/                 62   §5 — training, evaluation, interpretability
│                               (NMD_orf_model_v5_4ct v2.0.0, minus superseded/)
│
├── figures/
│   ├── lib/               15   ggplot_style (SSOT), validators, geometry
│   ├── multipanel/        90   Figures 1–5
│   └── supplemental/      93   SF1–SF43
│
├── verification/          22   reproducibility verifiers (13) + verify_section*.R (9)
│   └── baseline/           9   regression baseline for the rewrite
│
├── metadata/
│   ├── pheno/              3   donor / sample / treatment map
│   └── reference/          2   RBP rosters (Gerstberger, ENCODE eCLIP)
│
└── docs/
    ├── results_to_code_map.md  every claim → producing script
    └── METHODS/                §4 + model METHODS.md
```

Two structural fixes are folded in:

- **`results/isoform_transitions/Version_6.0/` → `analysis/isopair/`.** The old path is an
  accident of history and the direct cause of three `.gitignore` near-misses (a bare `results/`
  rule silently untracking files). A normal path removes that hazard permanently.
- **The model merges in as `model/`**, killing the absolute cross-repo `source()` at
  `figure5_panelA_architecture.R:29` that hard-fails for anyone who is not Pete.

---

## 2. Excluded, with reasons

| Excluded | n | Why |
|---|---|---|
| `superseded/` | 158 | Already-archived dev lineage. By definition not paper code. |
| `paper/` manuscripts + `paper/archive/` | 24 | Manuscript drafts are not code; the Doc is canonical. |
| `code/nmd_atlas/` | 8 | Web app. 0 map citations, 0 manuscript mentions. Pete ruled exclude (KEEP_LIST D-2). |
| `Version_6.0/results/ptc/` | 19 | rmats/PTC **data products** + docs; rmats appears 0 times in manuscript, supplement and map. Violates code-only. |
| `code/` root provenance files | 5 | `gsea_mashr_2026.3.10.R` (verified **not** the Tables 3–4 producer) and `isocall_limma_dge_fullmodel` (Pete's parallel QC, non-canonical). Kept previously *for provenance* — which "not more than that" excludes. |
| `code/ptcneg_go_handoff/` | 5 | A GO handoff bundle (CSVs + METHODS), not a producer of any manuscript number. **Confirm.** |
| Internal process docs | 8 | `KEEP_LIST`, `REPO_CONSOLIDATION_PLAN`, `SESSION_HANDOFF`, `DEPRECATED_CODE_INVENTORY`, `TD2_BIAS_AUDIT`, `REPRODUCIBILITY_PLAN`, `NEW_REPO_MANIFEST`, `ONBOARDING.md`. |
| `CLAUDE.md` | 1 | Internal working guidance; also carries the dual-push warning and lab-internal conventions. Not for a public deposit. |
| `paper/VERIFICATION_ISSUES.md`, `paper/MANUSCRIPT_REVIEW_FOR_YUL.md` | 2 | ⚠️ **Must not ship** — these enumerate unresolved manuscript problems by gene and cell type. Internal only. |
| `PJC_Tables/`, `lit_review/`, `isocall_dge/limma/README.md` | 4 | Manuscript tables (`.xlsx`, plus an `~$` lock-file artifact), a literature review, and a pointer doc. Not code. |

**Excluded total ≈ 234 files.** (ptc 19 and ptcneg 5 confirmed excluded.)

---

## 3. Decisions — SETTLED (Pete, 2026-07-24)

1. ✅ **`verification/` — INCLUDE.** The verifier suite + `verify_section*.R`. Makes the deposit
   self-checking, which is what Pete's "provided the correctness work is thorough" condition
   requires.
2. ✅ **`docs/results_to_code_map.md` — INCLUDE**, after stripping the internal
   verification-status annotations (which reference open manuscript issues).
3. ✅ **`Version_6.0/results/ptc/` (19) — EXCLUDE.**
4. ✅ **`code/ptcneg_go_handoff/` (5) — EXCLUDE.**
5. ✅ **Repo name: `nmd_lung_longread_2026`.**
6. ✅ **Isopair split confirmed** — `analysis/isopair/` holds the wrapper + drivers; the
   `Isopair` **package** remains its own repo with its own DOI.

---

## 4. What must be true before the DOI is minted

Under "not more than that", completeness has no safety net — there is no history to fall back
on. So the gate is:

1. `capture_verification_baseline.sh --check` passes: every number identical after the rewrite.
2. The clean-room test passes: fresh clone + Zenodo deposit, zero path edits, §1–§5 reproduce.
3. The old repos are **archived, not deleted**, so the development record survives if a gap
   appears later. Their existing DOIs remain valid historical records.
