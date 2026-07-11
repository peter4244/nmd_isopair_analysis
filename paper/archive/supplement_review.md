# Supplementary Material — review notes

**Source:** Google Doc fetched 2026-06-13. 6 tables + 31 SFx headers (28 with image markers, 3 without).

**Status:** First-pass review. Issues numbered S1–S20 below; priorities for Pete to triage at the bottom.

---

## Critical issues (fix before submission)

### S1. DIE limma table — column 5 is mislabeled

The DIE limma table (manuscript line 36–45) has a header that doesn't match its data:

| Cell Type | Tested | Sig | NMD-Responsive | "% of Sig that are NMD-Responsive" | Mean logFC |
|---|---|---|---|---|---|
| AT | 210,484 | 14,438 | 13,519 | **6.86** | 2.90 |

The column reads "Percent of Significant Transcripts that are NMD-Responsive" but the values shown are actually **% of TESTED transcripts that are SIGNIFICANT**:
- AT: 14,438 / 210,484 = **6.86%** ✓ (matches column 5)
- DD: 38,084 / 210,484 = **18.09%** ✓
- DD_ALI: 45,838 / 210,484 = **21.77%** ✓
- FB: 7,222 / 210,484 = **3.43%** ✓
- MV: 23,872 / 210,484 = **11.34%** ✓

The other three tables (SR Limma+Voom, LR Gene mashr, LR limma) use the correct %-NMD-of-significant computation (e.g., AT SR: 1,019 / 1,160 = 87.84%). **Either fix the header in the DIE limma table to "Percent of Tested Transcripts that are Significant" or recompute the column to actually mean %-NMD-of-significant** (which for AT DIE limma would be 13,519 / 14,438 = 93.6%). Recommend the latter for consistency with the other tables.

### S2. Three missing images

Three SFx have captions but no image markers in the doc:

- Line 83: **SFx — Expression Levels of Isoform Pair Sets** (cited in main text §4 paragraph 1)
- Line 89: **SFx — NMD Effect Size by Number of EJCs** (cited in main text §4 paragraph 3)
- Line 91: **SFx — NMD Pairs Flowchart** (cited in main text §4 paragraph 1)

All three are referenced in §4 text and need images before submission.

### S3. Cell-type labels are internal (pre-rename) — not manuscript-facing

The supplement tables use `AT`, `DD`, `DD_ALI`, `DO_ALI`, `FB`, `MV` (the internal labels). The manuscript main text uses **AT2** for `AT` and **LAE** for `DD` (per ONBOARDING §2 — April 2026 rename: `AT → AT2`, `DD = LAE` externally). Need to translate before submission. This also affects the Top Genes / Top Isoforms tables (lines 49, 64) which use `DD` instead of `LAE`.

---

## Substantive inconsistencies

### S4. Six-CT supplement vs four-CT main text scope

Supplement SR Limma+Voom, LR Gene mashr, LR limma, and DIE limma tables show all 6 CTs (`AT, DD, DD_ALI, DO_ALI, FB, MV`). Manuscript main text is scoped to 4 CTs (`AT2, LAE, FB, MV` — DD_ALI excluded for near-zero SR-vs-LR logFC correlation per ONBOARDING §2, DO_ALI excluded for reliability/donor n). Two options:

- **Keep 6-CT in supplement** — supplement carries the full data, main text reports the 4-CT scope. Document this scope choice in a single sentence at the top of the supplement.
- **Drop DD_ALI / DO_ALI from supplement** — match main-text scope.

Recommend the former; per ONBOARDING the 6-CT per-CT CSVs still exist in the repo precisely for this reason. But it should be stated.

### S5. Top Genes / Top Isoforms tables are 4-CT — inconsistent with other tables

While the SR/LR limma/mashr tables show 6 CTs, the Top Genes (line 49) and Top Isoforms (line 64) tables only show 4 CTs (AT, DD, FB, MV). Either include DD_ALI and DO_ALI in the Top tables, or apply the same 4-CT scope to all tables.

### S6. Supplement table counts don't match main-text-quoted counts

| Quantity | Supplement (6-CT) | Main text §2 (4-CT) | Likely cause |
|---|---|---|---|
| SR genes tested | 28,822 | 25,955 | Different filterByExpr cutoff because of different cell-type set |
| LR genes tested | 21,055 | 19,056 | Same |
| LR isoforms tested | 210,484 | 162,800 | Same |
| SR significant range across CTs | 14 – 4,075 | 3,122 – 6,753 | Limma adj.P.Val (supplement) vs **mashr posterior** (main text) |
| LR isoform significant range | 8 – 45,838 | 24,803 – 35,336 | Same: limma vs mashr |

The 4-CT vs 6-CT scope explains the testing-pool size differences. The significant-count differences are because **the supplement reports limma adj.P.Val while the manuscript main text reports mashr lfsr** (mashr borrows strength across CTs and shifts more features to significance for some CTs while pulling others toward shared signal). This is fine but the column header in the supplement says "FDR<0.05" which is ambiguous — mashr's column is lfsr, not FDR. **Clarify "FDR" vs "lfsr" in the supplement table headers.**

### S7. No DIE mashr table in the supplement

The supplement has DIE limma but **no DIE mashr** table. Manuscript main text reports DIE *mashr* numbers (24,803–35,336 sig per CT, 90.9–92.0% NMD susceptible). If the policy is "supplement carries the limma-only counts as a companion to the main-text mashr counts," fine — but a brief note at the top of the supplement explaining this would help readers.

### S8. SR / LR limma "FDR" column heading needs clarification

Same as S6 — the "FDR<0.05" column heading is used in tables that are limma (where FDR ≡ adj.P.Val) and tables that are mashr (where the equivalent is lfsr). The headers should distinguish. Specifically:

- "Short-Read Limma+Voom" — keep "FDR<0.05" (limma; correct)
- "Long-Read Gene mashr" — change to "lfsr<0.05" (mashr; current "FDR<0.05" is wrong)
- "Long-Read limma" — keep "FDR<0.05" (limma; correct)
- "DIE limma" — keep "FDR<0.05" (limma; correct)

### S9. NMD-susceptible definition wording is ambiguous in the mashr table

For the LR Gene mashr table, the column "NMD-Responsive Genes (FDR<0.05 & logFC>0)" should read "NMD-Responsive Genes (lfsr<0.05 & posterior_mean>0)" — the convention defined in main text Methods and ONBOARDING §6. Same fix as S8.

---

## Main-text → supplement reference gaps

These SFx are cited in the manuscript main text but **don't appear in the supplement**:

### S10. Section 1 — SFx that may or may not be in supplement

Manuscript §1 cites:
- "SF: correlation by ct" → matches **SFx Correlation by Cell Type** (supplement line 144) ✓
- "SFx — Pairwise Expression Similarity" → ambiguous; could match **SFx Pairwise Isoform Expression** (line 136) or **SFx Pairwise Similarity Across Cell Types for Isoform** (line 148). The naming should converge.

### S11. Section 3 — SFx compensation pathways MISSING

Manuscript §3 paragraph 7 cites: "(SF X — compensation pathways)" for the PCI pathway enrichments (oxphos, complex I, ATP synthesis). **No matching SFx in the supplement.** Either add it or remove the SFx reference and inline the pathway p-values in a main-text table.

### S12. Section 5 — GC content SFx MISSING

Manuscript §5 paragraph 5 cites GC content channel results: "GC content is a strong signal differentiating NMD from Control isoforms in the STOP but not the START window (SFx)." No matching SFx in supplement. The model repo has `model:09_export_gc_content.py` which produces the underlying data — the figure rendering may exist but the supplement entry is missing.

### S13. Section 2 — SFx volcano plot reference incomplete

Manuscript main text doesn't explicitly cite a volcano-plot SFx, but the supplement has **SFx Volcano Plot SR** (line 156). Either (a) add a main-text reference where appropriate, or (b) reclassify as supplement-only / data-availability.

---

## Supplement → main-text reference gaps (supplement-only figures)

These SFx are in the supplement but not directly referenced in main text §1–5:

### S14. SFx PCA Plot SR (line 120)
### S15. SFx PCA LR (line 126)
### S16. SFx Isoforms per Gene (line 123) — partially referenced; main-text says "median gene expressed five isoforms... 13,941 (76.3%) genes expressing two or more isoforms, 5,963 (32.6%) expressing ten or more" so this could fit the citation if explicit
### S17. SFx Change in Dominant Isoform Between Conditions (line 129)
### S18. SFx Volcano Plot SR (line 156) — see S13

For PCA plots, suggestion: add a brief mention in §1 or §2 ("Sample-level PCA showed clean separation by cell type and treatment, SFx PCA Plot SR / SFx PCA LR") so they're not orphaned in the supplement.

---

## Structural / editorial

### S19. Figure ordering

The 31 supplementary figures appear in an order that doesn't strictly follow main-text citation order. For final submission, recommend renumbering as `SF1, SF2, ...` and ordering them to match main-text citation order. This is also necessary to replace the placeholder `SFx` references in the main text with concrete numbers.

### S20. The supplement has no table titles using a numbering scheme

The 6 tables are headed "Short-Read Limma+Voom (SF)" etc. — "(SF)" suggests they're supplementary but they aren't numbered. Recommend `ST1, ST2, ...` (Supplementary Tables) numbering.

---

## Priorities for Pete to triage

**Block-or-fix-immediately** (would catch reviewer attention):
- **S1** — DIE limma column header / data mismatch
- **S2** — three missing images
- **S3** — cell-type labels need rename to manuscript-facing (`AT → AT2`, `DD → LAE`)
- **S11** — missing SFx for PCI compensation pathways

**Substantive, fix before final submission:**
- **S4** — document the 6-CT-supplement vs 4-CT-main-text scope choice
- **S6 / S8 / S9** — clarify limma "FDR<0.05" vs mashr "lfsr<0.05" in table headers
- **S7** — decide on DIE mashr table in supplement
- **S12** — GC content SFx needs to be added (Section 5)

**Editorial, do before submission:**
- **S5** — Top tables 4-CT vs other tables 6-CT (resolve once S4 is decided)
- **S10 / S13** — name reconciliation
- **S14–S18** — decide which orphan supplement figures need main-text references
- **S19 / S20** — concrete SF / ST numbering when freezing

---

*Review v0.1 — 2026-06-13. Pete to review and prioritize before applying to the Google Doc.*
