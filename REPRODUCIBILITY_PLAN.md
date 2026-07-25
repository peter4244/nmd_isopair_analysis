# Reproducibility plan — NMD long-read lung study

**v2, 2026-07-25.** Rewritten after an independent adversarial review found v1 unsafe to
execute. Findings are referenced by their review IDs (B1–B8, S1–S8, M1–M5).

---

## 0. Scope — this plan is Task 1 of two

| | Task 1 — **fidelity** (this plan) | Task 2 — **correctness** (`CORRECTNESS_PLAN.md`, deferred) |
|---|---|---|
| Question | Does the rewritten code, run from the deposit, produce **the same** numbers? | Are those numbers **right**? |
| Success | Every number identical | A wrong number is found and **changed** |
| Status | Instrument being rebuilt; rewrite blocked | §1–§3 done (14 findings); §4 46/46; §5 ~30/31 |

Their success criteria contradict, so they must not interleave: a regression would be
indistinguishable from a correction. Task 1 runs against **frozen analysis semantics** — the
rewrite changes paths, sample-name parsing and cell-type scope, never a computation.

---

## 1. Status

| Phase | State |
|---|---|
| Source-data deposit | ⚠️ built and uploaded as a **draft**, DOI `10.5281/zenodo.21544337` reserved — **but its column format is incompatible with the code (B4)** |
| New repository | ✅ `nmd_lung_longread_2026`, 395 files |
| Regression baseline v1 | ❌ **invalid** — measures cached outputs of the files the rewrite edits (B1) |
| **Instrument rebuild** | ⬜ **next — §5** |
| The rewrite | ⬜ blocked on §5 and §4 |
| Verification gates | ⬜ §6 |

**v1's central error.** The baseline read `nmd_fig_data/dge_*.rds` and the mashr CSVs — the
*outputs* of `Isoform_Level_Quantification.Rmd` and `NMD_shortread_dge_fullmodel_2026.5.5.Rmd`,
which are exactly the files §5.2 edits. A rewrite that broke the design matrix, reference level
or mashr fit would have produced a green check. The instrument measured something other than
the thing being changed.

---

## 2. Settled decisions

1. **Tier 1 only** in the deposit. 2. Full SQANTI classification. 3. Isopair cache regenerated, not
deposited. 4. Example-gene annotation — **now a blocker, see §4**. 5. Tan tables cited, not
redistributed. 6. ORF-scan objects are 4-CT-valid. 7. Source data gets its own DOI
(`…21544337`, reserved). 8. No sign-off needed to deposit derived objects. 9. One consolidated
code repo. 10. **Option C** — code is 4-cell-type native. 11. **The repo holds the code that
produced the paper's results, and not more**; correctness rests on verification, not provenance.

**New (v2):**

12. **The deposit is an interface, not just a file set.** Its column format is part of the
    contract with the code, and changing one requires changing the other (B4).
13. **`relevel(..., ref = "LAE")` is load-bearing** and must never be unified with display
    order (B8).
14. **No rewrite begins until the instrument can detect its failure** (B1, B2).

---

## 3. Deposit — status and the interface defect

`10.5281/zenodo.21544337` (reserved, draft). 9 files, 4.0 GB, at
`~/claude_projects/nmd_deposit_2026/`.

Validated: counts, classification, FASTA, GTF and GFF3 all carry an **identical 614,992 isoform
ID set**; genes 46,571. Trimming to observed features is result-neutral — `filterByExpr` sets
identical, percent-output-lost identical to 4 dp, and **TMM `norm.factors` bit-identical**
(verified 2026-07-25, closing review item M5).

Fixed after review: `.DS_Store` removed from the record (M2); `MANIFEST.sha256` entries
re-rooted so the documented `shasum -c` works (M1). **Both need re-uploading to the draft.**

### 3.1 ⛔ The column-format defect (B4)

The deposit uses `001V_DMSO_AT2`. Three parsers expect `Sample##_CT_ID_Treatment`:

| parser | effect when fed the deposit |
|---|---|
| `Isoform_Level_Quantification.Rmd:219-227` | parses treatment=`AT2`, donor=`DMSO`, ct=`001V`; duplicate check passes → **silent corruption** until `relevel(ref="LAE")` errors |
| `NMD_shortread_dge_fullmodel_2026.5.5.Rmd:108-135` | de-canonicalises AT2→AT, LAE→DD to match count columns; every sample "not found", dropped |
| `01_prepare_data_mashr.R` `parse_sample()` | same break |

**Resolution: change the parsers, not the deposit.** The canonical format is the point of
decision 10; reverting the deposit to `Sample##_` would reintroduce internal codes and defeat
it. This is **declared work in §5.3**, not cleanup.

---

## 4. ⛔ Missing inputs — triaged 2026-07-25

These are **missing producers or missing files**, not merely undeposited data. Investigated
before escalating (the `05s`/`05t`/`05u` recovery showed history often holds them).

### 4.1 Found locally — deposit gap only, no code missing

| Input | Located at | Action |
|---|---|---|
| `nmd_isocall.isoforms.gtf.gz` | `nmd/isocall/nmd_lungcells/results/call/` | add to the deposit |
| `FB_MV_pheno_2025.8.2.csv` | Dropbox `…/NMD/full_dataset/pheno/` | add to the deposit (or `metadata/pheno/`) |
| `DD_DO_AT_pheno_2025.8.15.csv` | Dropbox `…/NMD/full_dataset/pheno/` | same |
| `gencode.v49…annotation.sorted.gtf.gz` | `nmd/reference_files/` | external reference — **document the download**, do not redistribute |

### 4.2 Genuinely absent — must be sourced

| Input | Needed by | Where to look |
|---|---|---|
| `sub.isoforms.gtf` (+ `gene_ex/`) | `Figures/make_gene_examples.R:26` → **Figure 1 C/D/E** | **Ask Yul.** No producer anywhere: nothing in the repo, `superseded/`, or git history writes a GTF. Her repo was imported as *files, not history* (59 files), so a producer may exist there and never have come across. |
| ~~`tx_summary_6ct.tsv`~~ | ~~`model/relabel_tx_summary_4ct.R:45`~~ | ✅ **RESOLVED 2026-07-25 — not missing.** It is `export_rds.R`'s own `tx_summary.tsv` under an old name. The `_6ct` suffix records *when* it was generated, not 6-CT-specific content: the ORF summary is per-transcript ORF structure, and the README confirms the scan is "shared with original v5", i.e. cell-type-independent. The only gap was an undocumented manual rename, which existed so the relabel step would not read and overwrite its own output filename. **Fix in the rewrite:** have `export_rds.R` write `tx_summary_unlabeled.tsv` and `relabel_tx_summary_4ct.R` read that, writing `tx_summary.tsv`. Removes the manual step and the misleading name. |
| `merge-collapsed.gff` | `scripts/core/02_extract_isoform_structures.R:90` | PacBio isoseq collapse output. Likely Channing/Randell-side; confirm whether §4 still needs it or whether the SQANTI GTF supersedes it. |

⚠️ **Channing is unreachable from this machine** (`changit.bwh.harvard.edu` does not resolve
without VPN), so Yul's repo could not be inspected directly.

### 4.3 The ask for Yul

1. Does a script producing **`sub.isoforms.gtf`** / the `gene_ex/` example-gene annotation exist
   in her repo? It is the only input behind three main-figure panels with no producer at all.
2. Was anything else under `/udd/reyle/nmd_lungcells_2026/` **not** included in the 59-file
   import — particularly figure-support or data-prep scripts?

## 5. Rebuilding the instrument — before any rewrite

### 5.1 Complete the input audit

Enumerate every input of every shipped script; classify as deposited / produced-by-shipped-code
/ **missing**. Every missing entry is resolved, or named in `REPRODUCTION.md` as out of scope.

### 5.2 Make the verifiers survive the rewrite (B2)

`verify_section3_p2.R` and `verify_section2_claims_210_211.R` purl their target Rmds and patch
them by regex (`^proj <-`, `^resolve_mash_col <- function`, the `saveRDS(m, file.path(OUT_DIR`
truncation anchor). §5.3 destroys all of those, so the measuring instrument breaks — and once
edited to work again, the baseline it produced is no longer a baseline.

**Fix:** give the Rmds a stable, documented entry point (a `params:` block or sourced config
object) that the verifiers target by name and the rewrite may not rename. Refactor the
verifiers onto it, **then** re-baseline.

### 5.3 Rebuild the baseline so it measures the rewrite (B1)

New **gate 6.0**: regenerate the DGELists and mashr CSVs *from the deposit with the rewritten
code*, and bit-compare against the pre-rewrite copies retained locally. Claim-level verifiers
run only after that passes, and against deposit-derived inputs rather than `nmd_fig_data/`.

Also fix the harness (`capture_verification_baseline.sh`):

- **stop stripping warnings** (S2) — the rewrite changes factor handling and join keys, exactly
  what warnings would surface; sort/dedupe them instead of discarding
- fix the two dead filters (S3) — `select\(\) returned` never matches the real
  `'select()' returned …`, and unanchored terms (`tidyverse`, `geom_smooth`) can delete
  legitimate result lines
- update `SCRIPTS` for the new layout (`code/upstream/` → `verification/`)

### 5.4 Add §4/§5 fixtures (S4)

The nine verifiers cover §1–§3 only; §5's headline AUPRC 0.833 is unguarded. Pin, at minimum:
Isopair `data_mashr/*.rds` dimensions and column checksums; `tx_summary.tsv` label counts;
and `predictor_comparison/metrics_summary_2026.7.11.tsv` as an expected-output fixture.

### 5.5 Compute the expected diff mechanically (B3)

v1 predicted "exactly two lines". At least a third changes —
`isoform x ct x donor rows: 8388536` → `7994896` (645,272 × 13 → 614,992 × 13) — and
`gene-level CPM matrix: 30355` and the `unfiltered DGEList:` line are also at risk.
**Generate the literal expected diff by running the current verifiers against deposit-derived
inputs, and paste it into the plan.** A wrong prediction is worse than none: several
unexplained diffs arriving together invites absorbing them all as "expected".

---

## 6. The rewrite — declared work

### 6.1 Sample-name interface (B4) — newly declared

Rewrite the three parsers to consume `<donor>_<treatment>_<CT>` canonical columns directly,
deleting the `Sample##` reformat and the de-canonicalisation in
`NMD_shortread_dge_fullmodel:108-135`.

### 6.2 Hard-stop guards (B5) — newly declared

v1 asserted the sample-dropping blocks were "already no-ops". They are not:
`01_prepare_data_mashr.R:93-98` hard-`stop()`s when it cannot find two named `DO` samples, and
`:106-113` asserts `setequal(unique(ct), c("AT","DD","FB","MV"))` — both fatal against the
deposit. **Audit every `stop()`/`stopifnot()` referencing sample counts or cell-type codes.**

### 6.3 Factor-order contract (B8)

`relevel(factor(pheno$ct), ref="LAE")` gives design order `LAE, AT2, FB, MV`; display
constants use `AT2, LAE, FB, MV`; `makeContrasts` hardcodes the former and
`NMD_shortread_dge_fullmodel:296-297` pairs names **positionally**. Unifying them under one
constant would silently swap the AT2 and LAE result tables. Keep two named constants —
`CT_DESIGN_REF <- "LAE"` and `CT_DISPLAY_ORDER <- c("AT2","LAE","FB","MV")` — and never
substitute one for the other.

### 6.4 Cell-type code mapping must fail loudly (S1)

Three regimes coexist: consumers expect `at2`/`lae`, `model/relabel_tx_summary_4ct.R:14`
expects `at`/`dd`, and the tracked CSVs use `AT`/`DD`.
`Isoform-Level_DIE_Summary_p1.Rmd:129-135` **drops unmapped rows with a message**, so a
half-finished rename yields tables silently computed over two cell types. Add
`stopifnot(all(codes %in% names(CELLTYPE_MAP)))` and assert four cell types at each load site
**before** renaming anything.

### 6.5 Path portability — rescoped (S5, S6)

v1 said "17 of 59 files". Measured: **119 files carry an absolute path**; 54 carry
`/Users/petecastaldi/...`, of which 51 match none of v1's patterns; **15 `setwd()` calls** point
at a directory that no longer exists; and `model/` uses a fifth root
(`/projects/talisman/...`, `/home/p.castaldi/cc/...`).

- Linter regex `^[^#]*["'](/|~/)` plus an explicit `setwd(` ban, with a small allowlist.
- Config needs a **fourth root** beyond `DEPOSIT`/`CACHE`/`OUT`: the §1–§3 → §4/§5 handoff
  intermediates (`nmd_fig_data/`, `isocall_dge/mashr/`, `shortread_dge/mashr/`), with a
  documented build order.

### 6.6 4-cell-type native

Delete the drop blocks (after 6.2's audit), keep `CELL_TYPES` as a scope declaration only
(never as factor order — 6.3), normalise naming to `AT2`/`LAE`, update in-repo docs.

---

## 7. Verification gates

| Gate | Requirement |
|---|---|
| **6.0 Intermediate regeneration** *(new)* | DGELists and mashr CSVs regenerated from the deposit **bit-compare** to the retained pre-rewrite copies |
| **6.1 Regression check** | Every number identical apart from §5.5's mechanically-computed expected diff |
| **6.1b Stage-wise comparison** | Regenerate each stage from the deposit and diff against retained references — 5 DGELists, 20 mashr CSVs, 83 Isopair cache objects, 3 ORF-scan objects, model inputs. Localises a defect to a stage instead of a final number. These references are ours, not a reader's. |
| **6.2 Clean-room test** — *gate on the DOI* | Fresh clone, no Channing access, unpack the Zenodo record, run `REPRODUCTION.md`. Numbers exact; figures by **PNG** byte-compare (PDFs embed timestamps). |

⚠️ **`REPRODUCTION.md` must be written before 6.2, not after** (M4) — it is the artifact the
clean-room test executes. v1 sequenced it wrongly.

⚠️ **§1–§3 numbers were produced in Yul's Channing environment** (`/udd/reyle/Rlibs`), so
"reproduce exactly" is achievable against local re-derivations, not against the manuscript
itself. State this plainly in `REPRODUCTION.md` (S8).

---

## 8. Open items

| # | Item |
|---|---|
| A | `05_final_report_mashr.Rmd` — recommend excluding (deprecated, non-runnable, no shipped figure reads it). |
| C | `REPRODUCTION.md` unwritten — resequenced ahead of 6.2. |
| D | Zenodo metadata before publishing: author given-names, ORCIDs, related-identifier links. Immutable at mint. Re-upload the fixed `MANIFEST.sha256`. |
| E | ORFik version — **pre-check now**: regenerate `05s` under 1.30.2 and compare row counts/hashes to the cached `orfik_scan.rds`, rather than waiting for 6.2 (M4). |
| F | GEO GSE329233 reduced to four cell types, after verification. |
| G | 📌 Draft `CORRECTNESS_PLAN.md` (Task 2) once the rewrite is complete. |
| H | `ENVIRONMENT.md` lists msigdbr **25.1.1** in one table and **26.1.0** in another; these cannot both have produced Tables 3–4 (S8). Also records a `pathview` KEGG template fetched from the network at runtime — a live dependency inside the reproduction path. |
| I | `model/relabel_tx_summary_4ct.R` docstring says non-NMD is `adj.P.Val > 0.50`; the code sets `0.30` (M3). This defines the model's negative class. Reconcile before touching the file. |
