# Manuscript verification — notes for Yul (Results §1–§3)

**From:** Pete (verification pass, 2026-07-21 and 2026-07-24)
**Scope:** Results §1, §2 and §3 of manuscript version 2026.7.17. §4 and §5 were verified
separately and are not covered here.

---

## Headline

**Every computation reproduces.** Across all three sections I could not find an arithmetic,
pipeline, or code error. Where numbers are quoted, they are the numbers the code produces —
usually to the decimal (median adjusted log2FC 0.141 / −0.040 / 0.098 / 0.051; Kitagawa
86.2% / 13.5%; splicing OR 2.52; SR/hnRNP OR 6.76; PCNA 251.5 → 129.9; 25,955 genes;
162,800 isoforms; 34,387 unique NMD-susceptible isoforms).

The items below are therefore **not** about whether the analysis is right. They fall into three
kinds:

1. **Two places where the prose says more than the statistic shows** (§1.2, §1.3). A third
   (§1.1) I raised and then withdrew — recorded there so the reasoning is on file.
2. **A recurring pooling pattern** — cell-type-specific spread averaged into one number, with
   LAE the outlier every time (§2).
3. **Definitions a reader cannot reconstruct from Methods** (§3).

**How this was checked.** Rather than reimplement anything, your Rmds were purled and executed
against local copies of the same inputs, and the claims read off your own objects (`base`,
`prod_dir_mashr`, `kit_nmd_table()`, `theme_enrichment()`, `all_proportions`). One disclosed
deviation throughout: the mashr CSVs in Pete's repo predate the AT→AT2 / DD→LAE rename, so
`resolve_mash_col()`, `resolve_sr()` and `CELLTYPE_MAP` were widened to accept both spellings.
Column lookup only — no data, thresholds, or computations altered.

Scripts live in `code/upstream/`: `verify_section3_p{1,2,3,4,5}.R`,
`verify_section3_p1_metric_diagnostics.R`, `verify_section3_p2_claim311.R`,
`verify_section2_claims_210_211.R`, `verify_section1_claim113.R`, and the shared harness
`_run_productive_response.R`.

**One thing worth flagging in your favour:** `interpret_isoform_patterns_mashr…Rmd`'s
`q3-identify` chunk prints "(lfsr < 0.05)" while filtering on `adj.P.Val`. That looked like a
mismatch, but it is **correct** — in these files `adj.P.Val` *is* the mashr lfsr and `logFC`
*is* the mashr posterior mean (verified identical, max |difference| = 0, across all 8 per-CT
`2026.3.10` files). Pete's `CLAUDE.md` had this documented backwards; it has been fixed.

---

## 1. Needs a decision before submission

### 1.1 "% transcriptional output lost" — **no action needed** (an earlier concern of mine, withdrawn)

An earlier draft of these notes raised this as a serious issue. Having read the Methods in full
and examined SF22, **I withdraw it**, and record why in case it comes up in review.

The numbers reproduce exactly (gene-level 11.06–12.50; isoform-level 14.40–18.23).

**What I raised, and why it does not hold.** Because CPM columns each sum to 1e6,
`sum(Δ) = 0` for every donor pair, so `Σ max(Δ,0) ≡ Σ max(−Δ,0) ≡ ½ Σ|Δ|` — the measure equals
the **total variation distance** between the DMSO and SMG1i composition vectors. I read that as
a defect. It is not: it is a property of compositional data. The quantity being measured — what
fraction of the SMG1i library consists of material absent under DMSO — is meaningful; the
isoform *sets* contributing on each side differ; and the experimental design supplies the
direction.

I also argued the `max(Δ,0)` operator could inflate the estimate with noise. **SF22 already
answers this, and more directly than the check I ran** — its floor is applied to the per-feature
CPM *difference* (`SMG1i − DMSO ≥ t`), which is exactly where noise lives. Even at a ≥10 CPM
difference floor, isoform-level remains ~7–9% and gene-level ~9–11%. My independent sweep on
the expression level agrees: the gene-level measure is nearly invariant (11.06 → 10.89 from a
0 to 5 CPM floor), and the isoform-level measure holds at 10–13% when only 2.2% of isoforms
survive. **The signal is not a noise artifact.**

Finally, computing over *all* isoforms — rather than only lfsr<0.05 ones — is the right choice:
thresholding would make the measure power-dependent (n=3–4 donors) and discard genuine
sub-threshold targets.

**Optional, at most:** one Methods sentence noting the measure is compositional (relative to
library composition) and equals the total variation distance between conditions. That framing
is arguably a strength, not a caveat.

### 1.2 The *SHMT2* example (§3 ¶5, Fig 2C) — wrong cell type, and "significantly" is unsupported

> "…*SHMT2* was NMD susceptible yet **significantly** increased its productive output upon NMD
> inhibition (**MV**, 86.5 to 111.1 CPM; adjusted logFC +0.42) … (short-read logFC +0.56)."

**All four quoted values are AT2's, not MV's:**

| quantity | quoted | **AT2** | MV |
|---|---|---|---|
| productive DMSO CPM | 86.5 | **86.5** | 44.1 |
| productive SMG1i CPM | 111.1 | **111.1** (raw) | 54.8 |
| adjusted logFC | +0.42 | **0.421** | 0.391 |
| short-read logFC | +0.56 | **0.561** | 1.161 |

Fig 2C's legend also says "in MV", so if the panel plots MV, the figure and the prose show
different cell types.

**Separately, SHMT2 is `ns` in every cell type** by the paper's own rule: AT2 lfsr = **0.0582**,
MV 0.114, FB 0.162, LAE 0.376. AT2 is marginally above threshold. "Significantly" should come
out, or the sentence be reframed as an illustrative trend.

**Minor third point:** SRSF2 and PCNA quote the *adjusted* SMG1i CPM (50.4→"50"; 129.9→"130");
SHMT2's "111.1" is the *raw* value (adjusted = 116.2) while the "+0.42" beside it comes from the
adjusted figure.

The other two examples are exact, so this looks like a one-gene slip.

### 1.3 "Transcriptional-repression genes were **instead** enriched for DNA replication and repair (OR 2.1)"

Recomputed with your `theme_enrichment()` (transcriptional n=1,631; composition n=430):

| theme | trans_OR | trans_p | comp_OR | comp_p |
|---|---|---|---|---|
| RNA splicing | 1.19 | 0.10 (ns) | **2.52** | 9.7e-07 |
| SR/hnRNP factors | 0.76 | 0.76 (ns) | **6.76** | 2.6e-04 |
| DNA replication | **2.12** | 9.3e-07 | **2.03** | 4.5e-03 |
| DNA repair | **1.57** | 3.4e-05 | **1.88** | 4.2e-04 |

The quoted **OR 2.1 matches DNA replication only** (2.12) — DNA repair is 1.57. And **"instead"
implies a specificity the data do not show**: composition genes have nearly identical
DNA-replication enrichment (2.03 vs 2.12) and *higher* DNA-repair enrichment (1.88 vs 1.57).

The **composition half of the same sentence is excellent** — splicing 1.19 (ns) → 2.52 and
SR/hnRNP 0.76 (ns) → 6.76, exactly the double-dissociation the paragraph wants. Only the
transcriptional half overreaches. Simplest fix: quote DNA replication with its own OR and
soften "instead" for the DNA themes.

---

## 2. Pooled numbers that hide cell-type spread

The same issue three times, so it may warrant one consistent decision rather than three edits.
LAE is the outlier in every case.

### 2.1 "no significant change in roughly 83% of genes" (§3 ¶2)

Not a pooled proportion. Per cell type (NMD genes):

| AT2 | FB | LAE | MV |
|---|---|---|---|
| 82.9% | 93.8% | **59.8%** | 83.9% |

Pooled across all four is **79.7%** (75.2% excluding FB); mean of the per-CT values is 80.1%.
The only summary giving 83% is the **median across cell types (83.4%)**. The figure that
matters: in **LAE, 40% of NMD genes change significantly** — not 17%.

### 2.2 The DE-count ranges exclude FB, but the median list does not (§3 ¶2)

> "We found 88 (AT2) to 1337 (LAE) downregulated genes and 635 (MV) to 888 (LAE) upregulated."

All four numbers are exactly right, but FB is **down 11 / up 306** — below both stated minima,
so these are not min-to-max ranges. Excluding FB as provisional is reasonable; the inconsistency
is that the *same sentence* includes FB in the median list ("AT2 +0.14, LAE −0.04, **FB +0.10**,
MV +0.05").

### 2.3 The 13.5% composition share is essentially an LAE phenomenon (§3 ¶3)

| AT2 | MV | LAE |
|---|---|---|
| **1.1%** | **5.9%** | **21.3%** |

LAE contributes 56% of the pooled population (2,225 of 3,949), so it drives the 13.5%.
"Composition matters in only 13.5% of genes" reads as a transcriptome-wide constant when the
range is ~1% to ~21%. Since the paragraph's claim is mechanistic, the cell-type dependence may
be the more interesting result. (¶3's FB exclusion is also not restated in the prose; including
FB gives 87.1% / 12.5%.)

---

## 3. Definitions a reader cannot reconstruct (§1)

Two §1 claims quote percentages that do not reproduce, while their **counts are exact**.

### 3.1 "0.36–2.24% of each cell type's expressed isoforms" (1.12)

Counts are exact — **LAE 8,087 · AT2 2,267 · MV 2,015 · FB 1,131**. The percentages recompute
as **0.88–5.48%**, so the "expressed isoforms" denominator differs from the one I can build
from Methods.

### 3.2 "Only 27.8% (LAE) to 51.7% (MV) … passed the additional expression filter" (1.13)

Pete has already resolved the conceptual half: cell-type-restricted expression is defined over
a **larger universe than the 162,800 DIE isoforms**, since a restricted isoform may by
definition fail the expression filter. That is confirmed — defining restricted isoforms on the
**unfiltered** DGEList and intersecting with the DIE universe reproduces 1.12's counts exactly
(LAE 72,508 → 8,087; AT2 30,240 → 2,267; MV 39,456 → 2,015; FB 25,361 → 1,131).

What still does not reproduce is the percentages themselves:

| denominator | AT2 | LAE | FB | MV | ordering vs claim |
|---|---|---|---|---|---|
| unfiltered restricted set | 7.5 | 11.2 | 4.5 | 5.1 | ✗ (FB lowest, LAE highest) |
| filtered + per-CT `filterByExpr` | 31.7 | **25.6** | 37.2 | **45.0** | ✅ matches (LAE lowest, MV highest) |
| claimed | — | **27.8** | — | **51.7** | |

The per-CT variant has the right structure but runs 2–7 points low. Since
`Isoform_Landscape.Rmd` reads the **filtered** DGEList with `get_expressed = rowSums > 0`, the
landscape numbers are post-filter — which suggests the "additional expression filter" is a
*further* filter on top of the 162,800, and that "all reached significance" refers to that step
rather than the mashr DIE.

**Ask:** what exactly is the "expressed isoforms" denominator (1.12), and which filter and
which testing step does 1.13 refer to? Neither is currently specified in Methods, so a reader
cannot reproduce either percentage. The substantive counts are exact regardless.

---

## 4. Smaller items

| # | Item |
|---|---|
| **4.1** | **Why FB is provisional.** The manuscript says "due to **donor confounding** (see methods)". `productive_response.Rmd:61` says "smaller effect sizes; fails part of the reproducibility battery". `confound` never appears in the Rmd about FB, and no confounding discussion appears in the manuscript export (its Methods are abbreviated — it may be in the full version). Which is correct? |
| **4.2** | **SF21 quartiles.** "IQR 30–100%" does not hold for MV (Q1 **25.7**, Q3 **96.2**); AT2/FB/LAE are 30.2–31.6 / 100. Medians are exact (60.2–67.1). Also "IQR" is used for the Q1–Q3 *interval*, not the IQR width. |
| **4.3** | **SRSF2's cell type** is not given in the prose (SHMT2 and PCNA both name theirs). The values confirm **MV**, so the Fig 2D legend is right — just add it to the text. |
| **4.4** | **"the 12 canonical SR proteins"** — only **11** are in the data (SRSF1–SRSF11; SRSF12 absent). Counts unaffected. The Discussion already says "11 of the 12", so the two sections disagree. |
| **4.5** | **The UPR enrichment claim (§3 ¶2).** The NF-κB/TNFα half verifies cleanly (UP in NMD genes, all 4 CTs, p 4.3e-05…0.008). For UPR: all three named genes appear, but **no single cell type has all three** (AT2 HSPA5+EIF2AK3; LAE/MV HSPA5+XBP1; FB HSPA5 only). And `leading-edge-check` lists membership without testing enrichment — a direct Fisher test is significant in LAE (p=0.004) and MV (p=0.034) but not FB (0.053) or **AT2 (0.44)**. Since "leading-edge" implies a GSEA I could not locate (the claim cites ST4), **which analysis backs the "enriched" wording?** My Fisher test is a different test and should not be read as refuting it. |
| **4.6** | **§2's cell-type-specific pathways.** Two claims tracked internally — AT2-specific translation-regulatory / xenobiotic metabolism, and AT2+MV regulation of RNA splicing — could not be checked because the specifics are in **neither** manuscript version; the text now defers to "(supplemental results)". **Is there a supplemental results document to verify against?** |
| **4.8** | **SF23 vs the text.** §3 ¶2 claims "Spearman ≥ 0.96" for non-NMD genes, but SF23 itself prints LAE = **0.95** (AT2/FB/MV 0.97). My recomputation sides with the text (LAE **0.9643**, which rounds to 0.96), so SF23 was most likely rendered from an earlier run and not refreshed. Worth re-rendering, since a reviewer comparing the sentence to the figure it cites would see the stated floor violated. |
| **4.7** | The Kitagawa percentages sum to 99.7 / 99.5 / 99.9 because of **NA rows** (12 / 11 / 1) — genes in `prod_dir_mashr` dropped from `base` by the expression floor. Not ties (there are **zero**) and not an interaction term. Cosmetic; noted so it isn't mistaken for a bug. |

---

## 5. Reproduced exactly — no action needed

Recorded so the scope of what held up is clear.

**§1** — 162,800 isoforms and 18,270 genes; SQANTI categories (FSM 55,770 · ISM 31,667 · NIC …);
isoforms per gene (median 5 · mean 8.9 · IQR 2–13 · max 97); SR↔LR correlations (Pearson
0.83–0.91, Spearman 0.849–0.901); 105,938 filter-passing isoforms; the mashr marker counts in
all four cell types; CT-restricted counts LAE 8,087 · AT2 2,267 · MV 2,015 · FB 1,131; pairwise
Jaccard 0.863–0.923.

**§2** — 25,955 genes and 162,800 isoforms tested; significant genes 3,122–6,753 and isoforms
24,803–35,336; 34,387 unique NMD-susceptible isoforms; 19,803 (57.6%) core-shared and LAE
83.4% of CT-specific; SR-vs-LR effect-size ratio 2.0–3.2×; per-CT specificity LAE 34.2% /
23.7% vs 1.9–5.3% and 0.8–3.4%; pairwise sharing at both levels; the Tan et al. reanalysis
end-to-end (**3,069 shared · 49.8% · 40.6%**). Isoform proportions: median DMSO proportion
**0.1%**, 99.0–99.2% below 50%, expression spanning **4.21–4.94** log10; DMSO→SMG1i shift up in
**96.8–97.9%** (proportion) and **99.7–100%** (expression), paired Wilcoxon p ≈ 0.

**§3** — gene-level 11.06 / 11.83 / 12.03 / 12.50 and isoform-level 14.40 / 14.56 / 15.53 /
18.23, including which cell type is min and max; per-isoform medians 60.2–67.1; Spearman ≥0.99
(0.9980–0.9997) and ≥0.96 (0.9643–0.9705); median adjusted log2FC 0.141 / −0.040 / 0.098 /
0.051; DE counts 88 / 1337 / 635 / 888; Kitagawa 86.2% / 13.5%, up-arm 99.2% / 0.3%, down-arm
32.1% / 67.9% (identity holds to 9.1e−13); splicing OR 2.52 and SR/hnRNP OR 6.76; three-way
split transcriptional 70.1%; **SR proteins — 13 significant changes across AT2/LAE/MV, every one
a decrease, zero increases**; SRSF2 (MV) 84.1 → 50.4, −0.728, 91.0% → 26.4%; PCNA (LAE)
251.5 → 129.9, −0.947, 99.5% → 92.5%, SR −0.811.

**Producers traced:** Fig 2B → `productive_response.Rmd:505` (`dist-others`); Fig 2C/2D/2E →
`Isopair::plotIsoformPair()` (`Isopair/R/visualization.R:47`).

---

## 6. Questions

1. **§1.2** — is the *SHMT2* panel AT2 or MV, and can "significantly" be dropped?
2. **§1.3** — happy to quote DNA replication alone and soften "instead"?
3. **§2** — report per-cell-type ranges rather than pooled/median summaries, given LAE is
   consistently the outlier?
4. **§3** — what are the "expressed isoforms" denominator (1.12) and the "additional expression
   filter" plus its testing step (1.13)?
5. **§4.1** — which is the correct reason FB is provisional?
6. **§4.5** — which analysis produced the UPR "enriched" claim (ST4)?
7. **§4.6** — is there a supplemental results document holding the cell-type-specific pathway
   findings?
8. **§4.8** — can SF23 be re-rendered so its LAE Spearman matches the text?
