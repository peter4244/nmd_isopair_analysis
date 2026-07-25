# Reproducibility plan — NMD long-read lung study

**v3, 2026-07-25.** Restructured after Pete's observation that "rewrite" was too dramatic a
framing: the work splits into a mechanical path change and a semantic cell-type change, and
bundling them is what made v1/v2 hard to verify. Two independent reviews informed this
(findings referenced as R1-x / R2-x).

---

## 0. Scope — Task 1 of two

| | Task 1 — **fidelity** (this plan) | Task 2 — **correctness** (deferred) |
|---|---|---|
| Question | Same numbers as before? | Are the numbers right? |
| Success | Every number identical | A wrong number is found and changed |

Success criteria contradict, so they must not interleave. Task 1 runs against frozen analysis
semantics. 📌 `CORRECTNESS_PLAN.md` to be drafted once Task 1 completes.

---

## 1. Why phased

v1 and v2 both treated this as one "rewrite" and then tried to build an instrument capable of
verifying it. That failed twice — the instrument measured cached outputs of the files being
edited (R1-B1), and it sat inside its own blast radius (R1-B2).

Splitting the work gives each phase an invariant simple enough to check directly:

| Phase | Change | Invariant |
|---|---|---|
| **0** | Repoint the verification suite | It runs at all |
| **1** | Paths only — where files are read from | **Byte-identical outputs, same inputs** |
| **2a** | Cell-type handling, known sites | Still byte-identical on *current* data |
| **2b** | Point the config at the deposit | Only pre-declared differences |
| **3** | Cleanup — naming, dead code | Byte-identical to 2b |

Phase 1 needs no deposit, no pending decisions, and no elaborate harness: for a path-only
change the check is running each edited script and diffing its own outputs. That is a complete
test of inertness, and it dissolves R1-B1 — I am no longer asking a downstream verifier to
detect a break in the file that produced its input.

It also dissolves R2's circularity finding. The verifiers currently patch nine things in
`productive_response.Rmd`; **seven are pure path assignments**. Once paths resolve through
config, the verifier sets the config instead of rewriting source lines, so the fragile regex
anchors evaporate as a side effect rather than needing a special refactor first.

---

## 2. Phase 0 — make the verification suite runnable (hours)

Mechanical, independently checkable, prerequisite to everything (R2-7).

- `capture_verification_baseline.sh:20` — `BASE` points at `reproducibility/baseline_verification`;
  the baseline is at `verification/baseline_verification`. `verification/baseline/` is empty.
- All 15 verifier scripts hardcode `~/claude_projects/nmd` (the *excluded* old repo), plus three
  `setwd()` calls to it.
- Harness hygiene (R1-S2/S3): stop stripping warnings — the later phases change factor handling
  and join keys, exactly what warnings surface. Fix the two dead filters
  (`select\(\) returned` never matches the real `'select()' returned …`) and the unanchored
  terms that can delete legitimate result lines.

---

## 3. Phase 1 — paths only

**Rule: change where a file is read from; change nothing else.** No renaming, no dropping, no
restructuring. If a diff does anything but move a path, it belongs to a later phase.

1. `config/paths.yml` (currently `config/` is empty — R2-14) with one root per input class.
   At least six root families exist in shipped code (R2-12): `/udd/reyle/`, `/udd/repjc/`,
   `/proj/regeps/`, `/projects/talisman/`, `/home/p.castaldi/`, `~/claude_projects/`.
2. `R/load_config.R` + Python equivalent.
3. Convert path literals script by script, **verifying each by output diff before moving on**.
4. Linter: `^[^#]*["'](/|~/)` plus an explicit `setwd(` ban. Scope measured at ~99 files with an
   absolute-path literal, 54 containing `/Users/petecastaldi`, 15 `setwd()` calls (R2-11).

**Gate:** every edited script reproduces its own outputs byte-for-byte from unchanged inputs.

---

## 4. Phase 2a — cell-type handling, known sites

A focused pass over the sites we can already name, *before* pointing at the deposit — so the
empirical run in 2b is a safety net rather than a serial debugging loop.

**Still verified against current data**, so this stays byte-identical and any surprise is
attributable to the edit.

### 4.1 Will error loudly (easy)

| Site | Issue |
|---|---|
| `01_prepare_data_mashr.R:97` | `stop()` when two named `DO` samples are absent |
| `01_prepare_data_mashr.R:118` | `stopifnot(setequal(ct, c("AT","DD","FB","MV")))` |

### 4.2 Will fail *silently* — the dangerous class

| Site | Issue |
|---|---|
| `Isoform_Level_Quantification.Rmd:219-227` | `Sample##_CT_ID_Treatment` parser; canonical input scrambles treatment/donor/CT and **passes its own duplicate check** (R1-B4) |
| `01_prepare_data_mashr.R` `parse_sample()` | same |
| `render_sr_lr_correlation.R:61-71`, `correlation_analysis.Rmd:126-139` | DD_ALI disambiguation via `dge_gene_unfiltered`. Used **only** as a lookup, not as expression data, and guarded — so on canonical input it is a no-op. **Delete the dependency with the block** (resolves R2-1 without depositing DD_ALI data). |
| `NMD_shortread_dge_fullmodel:108-135` | de-canonicalises AT2→AT, LAE→DD to match count columns |
| `Isoform-Level_DIE_Summary_p1.Rmd:129-135` | drops unmapped cell-type rows **with a message** — a half-finished rename yields tables over two cell types |
| `CELLTYPE_MAP` sites | `Gene-Level_DGE_Summary_mashR.Rmd:64`, `die_mashr_enrichment_part2:56`, `Isoform-Level_DIE_Summary_p1:80` expect `at2`/`lae`; `model/relabel_tx_summary_4ct.R:14` expects `at`/`dd`; tracked CSVs use `AT`/`DD` |

### 4.3 Convert silent to loud — do this first

Because 2b's empirical run only finds what *breaks*:

- `stopifnot(all(codes %in% names(CELLTYPE_MAP)))` before any rename
- assert exactly four cell types at each load site
- assert expected sample count after every subset

Note two existing assertions already guard this path and should be kept:
`transcriptional_output.Rmd:355` and `comparison_analysis.Rmd:125-126`.

### 4.4 Do not touch

`relevel(factor(pheno$ct), ref = "LAE")` is **load-bearing** (R1-B8). Design order is
`LAE, AT2, FB, MV`; display constants are `AT2, LAE, FB, MV`; `makeContrasts` hardcodes the
former and `NMD_shortread_dge_fullmodel:296-297` pairs names **positionally**. Unifying them
silently swaps the AT2 and LAE result tables. Keep `CT_DESIGN_REF` and `CT_DISPLAY_ORDER`
separate and never substitute one for the other.

---

## 5. Phase 2b — point the config at the deposit

Now inputs genuinely change (4-CT, trimmed, canonical names). Run the pipeline and fix what
surfaces; 4.3's assertions make silent misbehaviour loud.

**Wiring, not blockers** (R2-5 corrected my triage):

| Item | Reality |
|---|---|
| `nmd_isocall.isoforms.gtf.gz` | repoint — the deposit GTF carries the same 614,992 transcripts and the needed `transcript_id`→`gene_id` |
| `gencode.v49…`, `merge-collapsed.gff` | in the **oarfish** branch of `02_extract_isoform_structures.R:85-91`; the manuscript path is isocall. Not required. |
| `tx_summary_6ct.tsv` | `export_rds.R`'s own output under an old name. Fix: write `tx_summary_unlabeled.tsv`, have the relabel step read that. |
| `sub.isoforms.gtf` | a **3-gene subset** (SHMT2/SRSF2/PCNA) for Fig 1 C/D/E — write a ~5-line producer over the deposit GTF |
| pheno CSVs | on disk (Dropbox); add to the deposit |

**Expected differences: compute mechanically.** Do not hand-predict — v1 and v2 both did and
both were wrong (R2-9). Run the verifiers against deposit-derived inputs and paste the literal
diff.

---

## 6. Phase 3 — cleanup

Delete the now-dead drop blocks, normalise naming to `AT2`/`LAE`, update docs. Verified
byte-identical against phase 2b.

---

## 7. Gates

| Gate | Requirement |
|---|---|
| Phase 1 | Each edited script's outputs byte-identical from unchanged inputs |
| Phase 2a | Still byte-identical on current data |
| Phase 2b | Only the mechanically-computed expected differences |
| Stage-wise | Regenerate each stage from the deposit, diff against retained references (5 DGELists, 20 mashr CSVs, 83 Isopair objects, 3 ORF-scan objects). **Typed comparison, not bit-compare** — references were produced on Channing against unknown package versions (R2-3): exact for IDs/factors/logicals/integers, tolerance for doubles, `setequal` on the mashr `strong` index and `filterByExpr` keep-sets. |
| Clean-room | Fresh clone + deposit, zero path edits. **Gate on the DOI.** |

⚠️ §1–§3 numbers were produced in Yul's Channing environment, so "reproduce exactly" holds
against local re-derivations, not against the manuscript. State this in `REPRODUCTION.md`.

⚠️ No reference PNGs exist for Figures 1–2 or SF1–24 (R2-13) — the code is present but the
images are not, so the clean-room image check has no §1–§2 targets. Either commit references or
scope the image check.

---

## 8. Open items

| # | Item |
|---|---|
| A | `05_final_report_mashr.Rmd` — recommend excluding (deprecated, non-runnable, no shipped figure reads it). |
| B | `REPRODUCTION.md` and repo `README.md` do not exist; both must precede the clean-room test. |
| C | Zenodo metadata before publishing: author given-names, ORCIDs, related-identifier links, and **clear the literal `TODO BEFORE PUBLISHING` string from `.zenodo.json`'s notes field** (R2). |
| D | Pin `useEnsembl(version = 115)` in `05u_paralog_annotation.R:90` — an unpinned live Ensembl query feeds `paralog_genes` → the model's `test_paralog` split → §5 AUPRC (R2-4). Same for the biotype fallback in `Isoform_Level_Quantification.Rmd:161-170`, and the `pathview` KEGG template fetched at runtime. |
| E | `ENVIRONMENT.md` lists msigdbr **25.1.1** and **26.1.0** in different tables (R2-16). |
| F | `DATA_INPUTS_NEEDED.md` and `IMPORT_PROVENANCE.md` shipped but carry internal plan references, stale claims and a PHI question — remove or rewrite before the DOI (R2-6). |
| G | `model/relabel_tx_summary_4ct.R:6` docstring says `> 0.50`; code sets `0.30`. Code is right (documented at `01_prepare_data_mashr.R:52-58`); fix the comment (R2-15). |
| H | Ask Yul: does a producer for `sub.isoforms.gtf` / `gene_ex/` exist? Only if the 5-line subset producer proves insufficient. |
| I | GEO GSE329233 reduced to four cell types, after verification. |
