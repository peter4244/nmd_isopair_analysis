# Section 3 verification — notes for Yul

**From:** Pete (verification pass run 2026-07-24)
**Scope:** Results §3, "Transcriptional Output Lost and Productive Isoform Response",
manuscript version 2026.7.17 — all 24 claims.

---

## Headline

**Every computation in §3 reproduces.** I re-ran the analysis end to end and could not find a
single arithmetic, pipeline, or code error. Where numbers are quoted, they are the numbers
the code produces — usually to the decimal (median adjusted log2FC 0.141 / −0.040 / 0.098 /
0.051; Kitagawa 86.2% / 13.5%; splicing OR 2.52; SR/hnRNP OR 6.76; PCNA 251.5 → 129.9).

The issues below are therefore **not** about whether the analysis is right. They are about
places where the **prose says something slightly different from what the statistic shows** —
mostly pooled summaries that hide cell-type spread, one attribution error, and one metric
whose name promises more than the math delivers.

**How this was checked.** Rather than reimplement anything, I purled
`productive_response.Rmd` and `transcriptional_output.Rmd` and executed **your code** against
local copies of the same inputs, then read the claims off your own objects (`base`,
`prod_dir_mashr`, `kit_nmd_table()`, `theme_enrichment()`). One disclosed deviation: the
mashr CSVs in Pete's repo predate the AT→AT2 / DD→LAE rename, so `resolve_mash_col()` and
`resolve_sr()` were widened to accept both spellings. Column lookup only — no data,
thresholds, or computations were altered.

Scripts: `code/upstream/verify_section3_p{1,2,3,4,5}.R`, plus
`verify_section3_p1_metric_diagnostics.R` and the shared harness
`_run_productive_response.R`.

---

## 1. Needs a decision before submission

### 1.1 "% transcriptional output lost" is mathematically a total variation distance

This one affects an **Abstract** claim ("NMD degraded 11–18% of transcriptional output"),
Fig 2A, SF20/21, and the Discussion's comparison to Fair et al.

The numbers reproduce exactly (gene-level 11.06–12.50; isoform-level 14.40–18.23). The
concern is what the statistic measures. Because CPM columns each sum to 1e6, `sum(Δ) = 0`
for every donor pair, and therefore

```
Σ max(Δ,0)  ≡  Σ max(−Δ,0)  ≡  ½ Σ|Δ|
```

so the quantity is the **total variation distance** between the DMSO and SMG1i composition
vectors. Three consequences I measured:

| diagnostic | result |
|---|---|
| Recompute in reverse (DMSO − SMG1i) | **identical** in all 4 CTs (14.40 / 14.56 / 18.23 / 15.53) |
| Share contributed by NMD-susceptible isoforms (lfsr<0.05 & PM>0) | **20.7–54.5%** (AT2 4.58 of 14.40 pts; FB 3.01 of 14.56) |
| Same statistic between two **DMSO** samples, different donors | 13.4–33.3%, exceeding the treatment value in 3 of 4 CTs |

**Two caveats stated up front**, because neither of these makes the biology wrong: the equal
forward/reverse totals are a *normalization* consequence — the isoform *sets* going up and
down are different, so this is not evidence the biology is symmetric. And the cross-donor
DMSO baseline includes real donor biology, so it is an upper bound, not a matched technical
null.

**Options:** (a) restrict the numerator to NMD-susceptible isoforms — gives 3.0–9.9%;
(b) describe it explicitly as a compositional-shift / TVD measure; (c) add a permutation
baseline showing the treatment effect exceeds donor-level variability; (d) keep the framing
and document why it is defensible. Any of these works — what would be uncomfortable is a
reviewer deriving the symmetry themselves.

*Evidence: `verify_section3_p1_metric_diagnostics.R`*

### 1.2 The *SHMT2* example (Fig 2C) — wrong cell type, and "significantly" is not supported

> "…*SHMT2* was NMD susceptible yet **significantly** increased its productive output upon NMD
> inhibition (**MV**, 86.5 to 111.1 CPM; adjusted logFC +0.42) … (short-read logFC +0.56)."

**All four quoted values are AT2's, not MV's:**

| quantity | quoted | **AT2** | MV |
|---|---|---|---|
| productive DMSO CPM | 86.5 | **86.5** | 44.1 |
| productive SMG1i CPM | 111.1 | **111.1** (raw) | 54.8 |
| adjusted logFC | +0.42 | **0.421** | 0.391 |
| short-read logFC | +0.56 | **0.561** | 1.161 |

Fig 2C's legend also says "in MV", so if the panel plots MV, the figure and the prose are
showing different cell types.

**Separately, SHMT2 is `ns` in every cell type** by the paper's own rule
(`dir = case_when(lfsr<0.05 & pm>0 ~ "up", …)`): AT2 lfsr = **0.0582**, MV 0.114, FB 0.162,
LAE 0.376. AT2 is marginally above threshold. The word "significantly" needs to come out, or
the sentence reframed as an illustrative trend.

**Minor third point:** SRSF2 and PCNA both quote the *adjusted* SMG1i CPM (50.4→"50";
129.9→"130"), but SHMT2's "111.1" is the *raw* value (adjusted = 116.2) while the "+0.42"
beside it comes from the adjusted figure.

The other two examples are exact — see §4 — so this looks like a one-gene slip.

*Evidence: `verify_section3_p5.R`*

### 1.3 "Transcriptional-repression genes were instead enriched for DNA replication and repair (OR 2.1)"

Recomputed with your `theme_enrichment()` (transcriptional n=1,631; composition n=430):

| theme | trans_OR | trans_p | comp_OR | comp_p |
|---|---|---|---|---|
| RNA splicing | 1.19 | 0.10 (ns) | **2.52** | 9.7e-07 |
| SR/hnRNP factors | 0.76 | 0.76 (ns) | **6.76** | 2.6e-04 |
| DNA replication | **2.12** | 9.3e-07 | **2.03** | 4.5e-03 |
| DNA repair | **1.57** | 3.4e-05 | **1.88** | 4.2e-04 |

Two issues. The quoted **OR 2.1 matches DNA replication only** (2.12) — DNA repair is 1.57.
And **"instead" implies a specificity the data do not show**: composition genes have nearly
identical DNA-replication enrichment (2.03 vs 2.12) and *higher* DNA-repair enrichment
(1.88 vs 1.57).

Worth emphasizing that the **composition half of the same sentence is excellent** — splicing
goes 1.19 (ns) → 2.52 and SR/hnRNP 0.76 (ns) → 6.76, which is exactly the clean
double-dissociation the paragraph wants. Only the transcriptional half overreaches. Simplest
fix: quote DNA replication with its own OR and soften "instead" for the DNA themes.

*Evidence: `verify_section3_p4.R`*

---

## 2. Pooled numbers that hide cell-type spread

These are the same issue three times, so they may be worth one consistent decision rather
than three separate edits.

### 2.1 "no significant change in roughly 83% of genes"

Not a pooled proportion. Per cell type (NMD genes):

| AT2 | FB | LAE | MV |
|---|---|---|---|
| 82.9% | 93.8% | **59.8%** | 83.9% |

Pooled across all four is **79.7%** (75.2% excluding FB); the mean of the per-CT values is
80.1%. The only summary that gives 83% is the **median across cell types (83.4%)**. The
figure that matters: in **LAE**, 40% of NMD genes change significantly — not 17%.

### 2.2 The DE-count ranges exclude FB, but the median list does not

> "We found 88 (AT2) to 1337 (LAE) downregulated genes and 635 (MV) to 888 (LAE) upregulated."

All four numbers are exactly right, but FB is **down 11 / up 306** — below both stated minima,
so these are not min-to-max ranges. Excluding FB as provisional is reasonable; the
inconsistency is that the *same sentence* includes FB in the median list ("AT2 +0.14, LAE
−0.04, **FB +0.10**, MV +0.05").

### 2.3 The 13.5% composition share is essentially an LAE phenomenon

| AT2 | MV | LAE |
|---|---|---|
| **1.1%** | **5.9%** | **21.3%** |

LAE contributes 56% of the pooled population (2,225 of 3,949), so it drives the 13.5%.
"Composition matters in only 13.5% of genes" reads as a transcriptome-wide constant when the
range is ~1% to ~21%. Since the paragraph's claim is mechanistic, the cell-type dependence
may be the more interesting result. (¶3's FB exclusion is also not restated in the prose;
including FB gives 87.1% / 12.5%.)

---

## 3. Smaller items

| # | Item |
|---|---|
| **3.1** | **Why FB is provisional.** The manuscript says "due to **donor confounding** (see methods)". `productive_response.Rmd:61` says "smaller effect sizes; fails part of the reproducibility battery". `confound` never appears in the Rmd about FB, and I found no confounding discussion in the manuscript export (though its Methods are abbreviated — it may be in the full version). Which is correct? |
| **3.2** | **SF21 quartiles.** "IQR 30–100%" does not hold for MV (Q1 **25.7**, Q3 **96.2**); AT2/FB/LAE are 30.2–31.6 / 100. Medians are exact (60.2–67.1). Also "IQR" is being used for the Q1–Q3 *interval*, not the IQR width. |
| **3.3** | **SRSF2's cell type** is not given in the prose (SHMT2 and PCNA both name theirs). The values confirm **MV**, so the Fig 2D legend is right — just add it to the text. |
| **3.4** | **"the 12 canonical SR proteins"** — only **11** are present in the data (SRSF1–SRSF11; SRSF12 absent). Counts are unaffected. The Discussion already says "11 of the 12", so the two sections disagree. |
| **3.5** | **The UPR enrichment claim.** The NF-κB/TNFα half verifies cleanly (UP in NMD genes, all 4 CTs, p 4.3e-05…0.008). For UPR: all three named genes appear, but **no single cell type has all three** (AT2 HSPA5+EIF2AK3; LAE/MV HSPA5+XBP1; FB HSPA5 only). And `leading-edge-check` lists membership without testing enrichment — a direct Fisher test is significant in LAE (p=0.004) and MV (p=0.034) but not FB (0.053) or **AT2 (0.44)**. Since "leading-edge" implies a GSEA I could not locate (the claim cites ST4), **which analysis backs the "enriched" wording?** My Fisher test is a different test and should not be read as refuting it. |
| **3.6** | The Kitagawa percentages sum to 99.7 / 99.5 / 99.9 because of **NA rows** (12 / 11 / 1) — genes in `prod_dir_mashr` dropped from `base` by the expression floor. Not ties (there are **zero** exact ties) and not an interaction term. Cosmetic; noted only so it isn't mistaken for a bug. |

---

## 4. Reproduced exactly — no action needed

Recorded so the scope of what held up is clear:

- **¶1** gene-level 11.06 / 11.83 / 12.03 / 12.50; isoform-level 14.40 / 14.56 / 15.53 / 18.23; per-isoform medians 60.2–67.1 — all matching, including which cell type is min and max.
- **¶2** Spearman ≥0.99 (observed 0.9980–0.9997) and ≥0.96 (0.9643–0.9705); median adjusted log2FC 0.141 / −0.040 / 0.098 / 0.051; DE counts 88 / 1337 / 635 / 888.
- **¶3** Kitagawa 86.2% / 13.5%; up-arm 99.2% / 0.3%; down-arm 32.1% / 67.9%. Decomposition identity holds to 9.1e−13.
- **¶4** splicing OR 2.52; SR/hnRNP OR 6.76; three-way split transcriptional 70.1% (consistent with Kitagawa's 67.9%). **SR proteins: 13 significant changes across AT2/LAE/MV, every one a decrease, zero increases** — that result is as clean as it reads.
- **¶5** SRSF2 (MV) 84.1 → 50.4, −0.728, 91.0% → 26.4%; PCNA (LAE) 251.5 → 129.9, −0.947, 99.5% → 92.5%, SR −0.811.
- **Producers traced:** Fig 2B → `productive_response.Rmd:505` (`dist-others`); Fig 2C/2D/2E → `Isopair::plotIsoformPair()` (`Isopair/R/visualization.R:47`).

---

## 5. Questions

1. **§1.1** — which framing do you want for "% transcriptional output lost"?
2. **§1.2** — is the *SHMT2* panel AT2 or MV, and can "significantly" be dropped?
3. **§1.3** — happy to quote DNA replication alone and soften "instead"?
4. **§2** — report per-cell-type ranges instead of pooled/median summaries, given LAE is consistently the outlier?
5. **§3.1** — which is the correct reason FB is provisional?
6. **§3.5** — which analysis produced the UPR "enriched" claim (ST4)?
