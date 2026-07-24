# Verification issues — running list

Open items surfaced by the 5-step verification passes, with routing. **Living document**
— append as issues are found, strike through when resolved. Started 2026-07-24 during the
§3 pass.

**Routing key:** 🔵 **PETE** (decision/judgment call) · 🟣 **YUL** (owns the analysis) ·
🟢 **CLAUDE** (mechanical, no decision needed) · ⚪ **INFO** (recorded, no action)

**Severity:** 🔴 blocks submission · 🟡 should fix before submission · 🟢 minor/polish

---

## Open

### F-3 🔴 🟣 YUL (+ 🔵 Pete) — "% transcriptional output lost" is a total variation distance

**Where it lands:** Abstract ("NMD degraded 11–18% of transcriptional output"), §3 ¶1,
Fig 2A, SF20, SF21, Discussion (comparison to Fair et al. 15%).

**Status of the numbers:** ✅ reproduce **exactly** — this is not a reproducibility problem.
It is an interpretation problem.

**The finding.** Every CPM column sums to exactly 1e6 (closed composition), so `sum(Δ) = 0`
and therefore `Σmax(Δ,0) ≡ Σmax(−Δ,0) ≡ ½Σ|Δ|`. The statistic **is** the total variation
distance between the DMSO and Smg1i isoform-composition vectors. Measured consequences:

| Diagnostic | Result |
|---|---|
| Reverse direction (DMSO − Smg1i) | **Identical** in all 4 CTs (14.40 / 14.56 / 18.23 / 15.53) — the number does not encode direction |
| Share from NMD-susceptible isoforms (lfsr<0.05 & PM>0) | **20.7–54.5%** (AT2 4.58 of 14.40; FB 3.01 of 14.56; LAE 9.93 of 18.23; MV 6.03 of 15.53 pts) |
| Same statistic, two **DMSO** samples, different donors | 13.4–33.3% — **exceeds** the treatment value in 3 of 4 CTs |

**Caveats (stated so this isn't overread).** The equal forward/reverse totals are a
*normalization* consequence; the isoform *sets* moving up vs down differ, so this is not
evidence the underlying biology is symmetric. The cross-donor DMSO baseline carries genuine
donor biology, so it is an upper bound on a no-treatment baseline, not a matched technical
null.

**Options to resolve:** (a) restrict the numerator to NMD-susceptible isoforms → 3.0–9.9%;
(b) relabel as a compositional-shift / TVD measure rather than "output lost"; (c) add a
permutation baseline demonstrating the treatment effect exceeds donor-level variability;
(d) argue the current framing is defensible and document why.

**Evidence:** `code/upstream/verify_section3_p1_metric_diagnostics.R` · commit `c74bf2c`

---

### F-8 🟡 🟣 YUL — §3 ¶3's pooled 13.5% composition share is almost entirely an LAE phenomenon

**Claim:** "…the level term was dominant in 86.2% of genes and the composition term in only
13.5%…" and "Among significant productive-down NMD genes, 32.1% were composition-dominated."

**All three percentages reproduce exactly.** The issue is again pooling (same class as
[F-6](#f-6)). Per cell type, composition-dominance:

| CT | composition-dominated | n genes |
|---|---|---|
| AT2 | **1.1%** | 871 |
| MV | **5.9%** | 853 |
| LAE | **21.3%** | 2,225 |

LAE supplies **56% of the pooled population**, so the 13.5% headline is driven by it.
"Composition matters in only 13.5% of genes" reads as a transcriptome-wide constant, when it
ranges from ~1% (AT2) to ~21% (LAE). Since the paragraph's conclusion is mechanistic — that
isoform-pool shifts are "a meaningful driver of productive loss" — the cell-type dependence
seems worth stating rather than averaging away.

**Also:** ¶3's population **excludes FB** (provisional), but the prose does not restate this;
a reader would reasonably assume all four cell types. Including FB gives 87.1% / 12.5%.

**Evidence:** `code/upstream/verify_section3_p3.R`

---

### F-6 🟡 🟣 YUL — §3 ¶2's "roughly 83%" is a median across cell types, not a pooled proportion

**Claim:** "…with no significant change in roughly 83% of genes **(Fig 2B; ST4)**."

Recomputed from Yul's own pipeline (NMD genes, mashr lfsr on the adjusted fit):

| CT | down | ns | up | total | % no change |
|---|---|---|---|---|---|
| AT2 | 88 | 4233 | 783 | 5104 | **82.9** |
| FB | 11 | 4784 | 306 | 5101 | **93.8** |
| LAE | 1337 | 3310 | 888 | 5535 | **59.8** |
| MV | 218 | 4448 | 635 | 5301 | **83.9** |

Candidate summaries: pooled all-4 **79.7%** · pooled excl. FB **75.2%** · mean of per-CT
**80.1%** · **median of per-CT 83.4%** ← the only one matching "roughly 83%".

**Why it matters:** the single figure hides a 59.8–93.8% spread. In **LAE** — the cell type
with the strongest NMD response, and the one supplying the top of most §3 ranges —
**40% of NMD genes change significantly**, not 17%. Suggest reporting the range
("60–94% depending on cell type") or the pooled value, rather than one number.

**Evidence:** `code/upstream/verify_section3_p2.R`

---

### F-5 🟡 🟣 YUL — §3 ¶2 DE-count ranges silently exclude FB

**Claim:** "We found 88 (AT2) to 1337 (LAE) downregulated genes and 635 (MV) to 888 (LAE)
upregulated genes **(ST4)**."

**All four quoted numbers are exactly correct.** The problem is the min-to-max framing:

| | claimed range | actual min | actual max |
|---|---|---|---|
| downregulated | 88 (AT2) – 1337 (LAE) | **FB 11** | LAE 1337 ✅ |
| upregulated | 635 (MV) – 888 (LAE) | **FB 306** | LAE 888 ✅ |

FB falls **below both stated minima**, so the ranges are not min-to-max across cell types.
FB is declared provisional, so excluding it is defensible — but **the same sentence's median
list includes FB** ("AT2 +0.14, LAE −0.04, **FB +0.10**, MV +0.05"), making the treatment
inconsistent within a single paragraph. Fix: say "excluding FB (provisional)" on the counts,
or report FB's values.

---

### F-7 🟡 🟣 YUL — the stated reason FB is provisional does not match the analysis code

**Manuscript:** "…and considered FB provisional due to **donor confounding** (see methods)."

**`productive_response.Rmd:61`:** "FB is treated as **provisional** (smaller effect sizes;
fails part of the reproducibility battery)."

These are different rationales. `confound` appears only 3× in the Rmd, none about FB's status
(all in Step 14, about gene-level predictors). Additionally the "(see methods)" pointer looks
unfulfilled — no discussion of donor confounding appears in `NMD manuscript 2026.7.17.md`.

**Caveat:** that file is a Google Doc export with an abbreviated Methods section, so the
explanation may exist in the full Methods. Needs confirmation of which reason is correct, and
that Methods actually covers it.

---

### F-1 🟡 🟣 YUL — SF21 quartile range is wrong for MV

"IQR 30–100%" does not hold for MV: observed quartiles **Q1 25.7 / Q3 96.2**
(AT2 31.4/100, FB 30.2/100, LAE 31.6/100). Either widen to ~26–100% or scope the statement
to AT2/FB/LAE. Separately, the manuscript uses "IQR" to mean the Q1–Q3 *interval*, not the
IQR width (68.4–70.5) — worth making explicit.

Medians are **exact** as published (60.2–67.1 vs claimed 60–67).

**Evidence:** `code/upstream/verify_section3_p1.R`

---

### F-4 🟢 🟣 YUL — SRSF2 example does not state its cell type in prose

§3 ¶5 gives *SHMT2* (MV) and *PCNA* (LAE) with cell types in the prose, but *SRSF2* has none
— only the Fig 2D legend supplies "in MV". Add for parallelism.

---

### O-1 🔵 PETE — Zenodo author metadata is minimal on all four records

All four deposits were minted without a `.zenodo.json` / `CITATION.cff`, so each record
lists only the repo owner. Metadata is editable post-publication and **the DOI does not
change**. Needs the full author list + ORCIDs, ideally with Yul.

Records: `nmd_isopair_analysis` 10.5281/zenodo.21539735 · `Isopair` …21536495 ·
`Isocall_v1` …21536486 · `NMD_orf_model_v5_4ct` v2.0.0 …21539601

---

### O-2 🔵 PETE — Code Availability section not yet in the manuscript

Drafted and DOI-verified at `paper/code_availability.md`; still needs pasting into the
Google Doc, with "Figure 5" as a live cross-reference field.

---

### O-3 🟢 🟢 CLAUDE — check whether the retired SR↔LR coverage claim migrated

Old §3 claim 3.10 ("~40% of SR NMD-susceptible genes lacked sufficient LR isoform coverage")
left §3 in the rewrite. If it now appears in §1/§2, its undefined "sufficient coverage"
threshold still needs pinning from `yul:comparison_analysis.Rmd` (now local). If it is gone
from the manuscript entirely, close this.

---

## Resolved

### F-2 ✅ — map described the output-lost metric as NMD-susceptible-only

It sums over **all** 645,272 isoforms with no susceptibility filter. Map corrected
2026-07-24 (commit `c74bf2c`).

### §3 map staleness ✅ — PCI / GPR180 framework retired

`results_to_code_map.md` §3 mapped an analysis absent from the manuscript (old claims
3.6–3.26). Rewritten against manuscript 2026.7.17; consistency fixes applied to the summary
table, M10, and open items 10–12. Commit `c74bf2c`.
