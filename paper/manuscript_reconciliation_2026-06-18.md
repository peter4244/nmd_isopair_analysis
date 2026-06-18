# Manuscript ↔ Figures/Methods Reconciliation — §4, §5, Methods

**Date:** 2026-06-18 (v2 — corrections applied 2026-06-18 PM per Pete's review)

## v2 corrections (applied after Pete's review)

1. **Retracted** the §4¶8 p-threshold find/replace. The string `***p<10⁻⁴, **p<10⁻³, *p<0.05, n.s.` in the Fig 4 composite legend is the **significance-key for the in-plot asterisks**, NOT the actual computed p-values. The data p-values are much tighter (10⁻¹², 10⁻⁶, 10⁻¹⁶ for the three Panel-C pairwise contrasts) and are correctly reported elsewhere; the asterisk key is independent.
2. **Significance-notation consistency** — added §1.10 (new) below: project currently runs three different schemes; recommend standardizing on the four-tier scheme `****p<10⁻¹⁰, ***p<10⁻⁴, **p<10⁻³, *p<0.05, n.s.`
3. **Retracted** the §5¶9 stop-codon find/replace. The "52.0% vs 48.3%" I put in the StopCodonUsage SF I built is from the **pre-patch** v5 HTML render (Apr 2; patch is Apr 30). The Rmd's post-patch narrative line cites **55.9% vs 48.1%**, which rounds to the manuscript's "56% vs 48%" — so the manuscript is essentially correct and the SF legend is stale.
4. **New** finding/correction noted in §1.11 below: the StopCodonUsage SF legend needs to be updated to post-patch numbers (and re-verified by re-rendering the Rmd when the cluster is up).
5. Methods drafts extracted to a standalone companion file at `paper/methods_updates_2026-06-18.md` for easier Google-Doc pasting (per Pete's request).

---


**Manuscript source:** Google Doc — https://docs.google.com/document/d/1Tz6coXnDwpGZaV1Jl11fmN_LVX1R8YQPf2y_7wzBLKs/edit
**Figure-side state:** current files in `figures/multipanel/figure{3,4,5}_*/` and `figures/SupplementalFigures/*/`
**Methods state:** as fetched from Google Doc this session (placeholders + `[insert]` tokens still present in some subsections)

This is the **5-pass scientific report verification** applied to §4 (Isopair), §5 (deep-learning model), and their related Methods subsections. Numbers are recomputed from figure-side TSVs, not from cached manuscript text.

Output is organized as:
1. **§4 reconciliation** — claim by claim, with find/replace pairs and gap flags
2. **§5 reconciliation** — claim by claim, with find/replace pairs
3. **Methods audit** — drafted new/updated subsections (ready for Google-Doc insertion)
4. **Summary of changes** — all find/replace pairs in one block + a punch list

---

## 1. §4 reconciliation (Isopair pair analysis)

### §4¶1 — Gene/isoform inclusion criteria — ✓ matches
| Claim | Figure-side | Verdict |
|---|---|---|
| 3,009 genes meeting criteria | Fig 3 schematic n=3,009; pop_BC = 3,009 NMD + 3,009 Control | ✓ |
| median 7 isoforms per gene | (not on a panel; presumably in pair-set descriptive SF) | UNVERIFIED — see Gap-§4-A |
| 75% of reference isoforms ≥50% of parent gene expression | (not on a panel; descriptive only) | UNVERIFIED — see Gap-§4-A |
| median reference isoform = 70% of parent expression | (not on a panel; descriptive only) | UNVERIFIED — see Gap-§4-A |

**Gap-§4-A — RESOLVED.** New SF built at `figures/SupplementalFigures/PairSetDescriptives/` (Panels A + B + C). Pete updated the manuscript §4¶1 sentence on 2026-06-18 to match the SF values exactly:

> "These genes had a median of 7 isoforms each (SFx - Isoform Count in Isoform Pair Sets), and the median reference isoform accounted for 31% of its parent gene expression, with 38% of the reference isoforms accounting for >50% of parent gene expression."

All numbers verified ✓:
- 3,009 pair-genes
- Median 7 isoforms/gene (IQR = 5)
- Median reference share = 31% (computed 30.8%, rounds to 31%)
- 38% of references ≥ 50% (computed 37.6%, rounds to 38%)

No further action needed on §4-A. The "SFx" placeholder in the revised §4¶1 sentence resolves to this Supplemental Figure.

### §4¶2 — Length comparisons — ✗ missing figure
| Claim | Figure-side | Verdict |
|---|---|---|
| median NMD comparator length = 3,049 nt | not in inventory | MISSING — see Gap-§4-B |
| median ref isoform length = 2,893 nt | not in inventory | MISSING |
| median Control comparator = 2,762 nt | not in inventory | MISSING |
| p < 0.001 (Control vs ref+NMD comparator length) | no figure cited | MISSING test name |

**Gap-§4-B — RESOLVED:** Panel C of the new `PairSetDescriptives/` SF carries the transcript-length distribution at pop_BC scope. All three medians verified exactly against the manuscript: NMD comparator 3,049 nt, reference 2,893 nt, Control comparator 2,762 nt. The unspecified test name fills in as **Kruskal–Wallis omnibus p = 1.5×10⁻⁷** with pairwise Wilcoxon (Mann–Whitney U) follow-up (p = 3.4×10⁻² / 2.1×10⁻⁸ / 7.9×10⁻⁴).

### §4¶3 — Sequence similarity + event prevalence — ✓ direction matches; minor wording
| Claim | Figure-side | Verdict |
|---|---|---|
| NMD shares more sequence with reference (qualitative) | Fig 3 Panel B (density of exonic seq shared) | ✓ direction matches |
| SE 2× more frequent in NMD than Control | Fig 3 Panel C: 44.2% NMD vs 21.2% Control (ratio 2.08×, p=1.18×10⁻⁸¹, OR 2.94×) | ✓ |
| IR more common in Control | Fig 3 Panel C: 9.9% NMD vs 12.7% Control (Control-enriched, p=7.28×10⁻⁴) | ✓ |
| (no exact prevalence numbers cited) | Panel C has full table | OK — Panel C carries the table |

### §4¶4 — GENCODE-annotated PTC analysis — ✓ exact match
| Claim | Figure-side | Verdict |
|---|---|---|
| n = 190 isoform pairs (all-3-ENST + coding) | Fig 3 Panel D scope = 190 NMD, 190 Control | ✓ |
| PTC rate NMD = 37.9% (72/190) | Panel D: exactly 37.9% (72/190) | ✓ |
| PTC rate Control = 2.1% (4/190) | Panel D: exactly 2.1% (4/190) | ✓ |
| 18-fold enrichment, p<10⁻²⁰ | Panel D OR = 18.04×, Fisher p reported in TSV | ✓ |
| dose-response between PTC distance and NMD magnitude | (qualitative, Panel D scatter) | ✓ direction matches |
| effect size peaks at 4–5 downstream EJCs | (qualitative, Panel D) | ✓ direction matches |

### §4¶5 — Splicing event mechanisms — ✓ mostly matches; one p-value clarification
| Claim | Figure-side | Verdict |
|---|---|---|
| Frameshift = 55% (38/69) | Fig 3 Panel F: 38 frameshift = 55.1% of 69 | ✓ |
| In-frame stop = 33% (23/69) | Panel F: 23 in-frame stop = 33.3% | ✓ |
| 3′ UTR splice = 12% (8/69) | Panel F: 8 3'UTR splice = 11.6% (rounds to 12%) | ✓ |
| SE = 43.5% of PTC-causing in NMD vs 14.1% Control | Fig 3 Panel E: exactly 43.5% / 14.1% (Fisher p=7.98×10⁻⁸) | ✓ |
| A5SS enriched 13.0% vs 4.5%, p=9×10⁻³ | Panel E: 13.0% / 4.5%, p=8.88×10⁻³ | ✓ (rounds to 9×10⁻³) |
| Alt TES depleted 5.8% vs 26.8%, p=3×10⁻⁵ | Panel E: 5.8% / 26.8%, p=3.37×10⁻⁵ | ✓ |
| Fisher p ≈ 10⁻⁸¹ "(SE prevalence in NMD vs Control overall)" | Fig 3 Panel C: 1.18×10⁻⁸¹ | ✓ but verify the sentence makes clear this is the Panel-C overall prevalence p, NOT the Panel-E PTC-attributed-events p (8×10⁻⁸) |

**Find/replace §4¶5 — clarifier (only if ambiguous in the Doc):**
```
find:    Fisher p ≈ 10⁻⁸¹
replace: Fisher p ≈ 10⁻⁸¹ for the overall SE prevalence comparison (Fig 3 Panel C); p = 8×10⁻⁸ for the PTC-attributed-event enrichment (Fig 3 Panel E)
```
Apply only if the manuscript currently cites 10⁻⁸¹ in a sentence that reads as if it's the Panel-E PTC-attributed-events p. If the sentence already clearly attributes it to overall prevalence, skip.

### §4¶6 — GENCODE-restricted NMD+/PTC− 5′UTR — ✓ exact match
| Claim | Figure-side | Verdict |
|---|---|---|
| n = 72 NMD+/PTC+, 118 NMD+/PTC−, 190 Control | Fig 4 A/B scope: 72 / 118 / 190 | ✓ |
| NMD+/PTC− have longer 5′UTRs than both other groups | Fig 4 Panel A: median 426.5 nt PTC− vs 84.5 PTC+ vs 138 Control; both pairwise p≤10⁻¹⁷ | ✓ |
| NMD+/PTC− have longer uORFs | Fig 4 Panel B: median 224 nt PTC− vs zero-spike for PTC+ and Control | ✓ |

### §4¶7 — TD2 CDS bias — ✓ exact match
| Claim | Figure-side | Verdict |
|---|---|---|
| n = 1,166 reference-AUG-traceable pairs | matches all of: Fig 4 C/D scope, Fig 5 Panel G scope, TD2BiasEvidence broad scope | ✓ |
| 50% TD2 picks same AUG as reference | TD2BiasEvidence Panel C: 50% same (578/1,166) | ✓ |
| 99% (487/492) of TD2 CDS downstream from ref AUG (occult-PTC subset) | TD2BiasEvidence Panel F: 99.0% (487/492), median offset +476 nt | ✓ |
| 78% (384/492) ref AUG has stronger Kozak | TD2BiasEvidence Panel E: 78.0%, paired Wilcoxon p=2.6×10⁻³⁸ | ✓ |

### §4¶8 — Reference-AUG-traced expanded analysis — ✗ off-by-3 on PTC− n
| Claim | Figure-side | Verdict |
|---|---|---|
| ref-AUG tracing identifies PTCs in 90% (1,050/1,166) NMD+ isoforms | Fig 4 Panel C: 1,050 NMD+/PTC+ ✓; 1,050/1,166 = 90.1% ✓ | ✓ |
| **NMD+/PTC− subset (n = 113)** | **Fig 4 Panel C scope: n = 116 PTC−** | **✗ off by 3** |
| Wilcoxon with HL shift (p < 10⁻⁴, 10⁻³, 0.05) | Fig 4 Panel C pairwise: PTC+ vs PTC− p=1.24×10⁻¹², PTC+ vs Ctrl p=4.85×10⁻⁶, PTC− vs Ctrl p=8.10×10⁻¹⁷ | ✗ thresholds also off — current p-values are MUCH smaller |

**Find/replace §4¶8 — sample count:**
```
find:    NMD+/PTC- subset (n = 113)
replace: NMD+/PTC- subset (n = 116)
```

**~~Find/replace §4¶8 — p-value thresholds~~ — RETRACTED.** The `p < 10⁻⁴, 10⁻³, 0.05` string is the significance-key for the in-plot asterisks (what `***`, `**`, `*` mean on the violins), NOT the data p-values. The data p-values are correctly reported in `panelC_pairwise.tsv` (p = 1.24×10⁻¹², 4.85×10⁻⁶, 8.10×10⁻¹⁷). The asterisk key in the legend is a separate notational convention; do not edit it. See §1.10 below for the broader significance-notation consistency question.

### §1.10 — Significance-notation consistency (project-wide)

Three different asterisk schemes are currently in use across the legends. They should be standardized to one project-wide convention.

| Figure / legend file | Current asterisk scheme |
|---|---|
| Fig 3 composite legend (Panel F line) | `***p < 0.001, **p < 0.01, *p < 0.05, ns` |
| Fig 4 composite legend (Panel D line) | `***p < 10⁻⁴, **p < 10⁻³, *p < 0.05, n.s.` |
| CDSand3UTR_GENCODEonly SF legend | `***p < 10⁻⁴, **p < 10⁻³, *p < 0.05, n.s.` |
| Fig 5 composite legend (Panel G line) | `****p < 10⁻¹⁰, *p < 0.05` (two-tier only) |
| Fig 3 Panels B/C/E, all other supps | per-test p-values reported directly (no asterisks) |

**Recommendation: standardize on a 4-tier scheme.**

```
****p < 10⁻¹⁰   (use only when the data warrants it — Fig 5G, Fig 3C overall SE prevalence, etc.)
***p < 10⁻⁴
**p < 10⁻³
*p < 0.05
n.s. otherwise
```

Rationale: p-values in this paper span 17+ orders of magnitude (10⁻¹⁰⁻¹³⁹ for some Fig 3 Panel C and Fig 4 Panel C contrasts). The classical `0.001 / 0.01 / 0.05` 3-tier convention saturates everything at `***` and loses resolution. The Fig 4 / CDS3UTR `10⁻⁴ / 10⁻³ / 0.05` scheme distinguishes "strong" from "very strong" but doesn't have a top tier for the extreme cases. The 4-tier scheme above subsumes both existing schemes plus the `****` already in Fig 5 Panel G.

**Files to update (use a single legend-pasteable string everywhere):**

```
Stars: ****p < 10⁻¹⁰, ***p < 10⁻⁴, **p < 10⁻³, *p < 0.05, n.s. otherwise.
```

| Legend file | Current | Action |
|---|---|---|
| `figures/multipanel/figure3_isopair_and_ptc/figure3_composite_legend.md` | `***p < 0.001, **p < 0.01, *p < 0.05, ns` | Replace with standard string |
| `figures/multipanel/figure4_ptcneg_and_model/figure4_composite_legend.md` | `***p < 10⁻⁴, **p < 10⁻³, *p < 0.05, n.s.` | Extend with `****p < 10⁻¹⁰` prefix |
| `figures/multipanel/figure5_dl_model/figure5_composite_legend.md` (Panel G) | `****p < 10⁻¹⁰, *p < 0.05` | Extend with `***` and `**` middle tiers |
| `figures/SupplementalFigures/CDSand3UTR_GENCODEonly/figure_s_cds_and_3utr_legend.md` | `***p < 10⁻⁴, **p < 10⁻³, *p < 0.05, n.s.` | Extend with `****p < 10⁻¹⁰` prefix |

**Plot-side stars** (the actual `*` annotations on the violins / bars) should also be regenerated to use the same thresholds. Most current plots produce per-pairwise p-values directly via `data_export.R` `_pairwise.tsv`; the bar-annotation helper that converts p → asterisks should be updated once at `figures/lib/` (one helper, called from all panels) rather than per-panel.

### §1.11 — StopCodonUsage SF carries pre-patch numbers

The current StopCodonUsage SF legend reads:
> "UGA is enriched in NMD (52.0% vs 48.3% in Control), UAA is depleted (27.1% vs 30.3%), and UAG is unchanged (20.9% vs 21.3%)"

These numbers come from the **pre-patch** v5 HTML render (`orf_model_report_v5.html`, file timestamp Apr 2 10:07; bug patch `scripts/patch_stop_codon.py` is dated Apr 30). The Rmd's post-patch narrative at `orf_model_report_v5.Rmd:1249` cites the canonical post-patch number as:

> "TGA is enriched in NMD (55.9% vs 48.1% in Control)"

The manuscript's "56% vs 48%" rounds cleanly from the post-patch 55.9% / 48.1% and should NOT be touched. The StopCodonUsage SF legend should be updated to:

```
find:    UGA is enriched in NMD (52.0% vs 48.3% in Control), UAA is depleted (27.1% vs 30.3%), and UAG is unchanged (20.9% vs 21.3%)
replace: UGA is enriched in NMD (55.9% vs 48.1% in Control); UAA and UAG are both depleted in NMD
```

(I do not have the post-patch UAA / UAG breakdown locally — the symlinked `selected_orfs.tsv` is unavailable. When the cluster is back and the Rmd is re-rendered, populate UAA and UAG percentages from the post-patch `sc_freq` table and tighten the legend.)

The SF folder also still embeds a base64-extracted PNG from the **pre-patch** HTML. The PNG should be regenerated from a post-patch Rmd render before this SF is included in the submission.

### §4¶9 — 3′ UTR length — ⚠ verify numbers
| Claim | Figure-side | Verdict |
|---|---|---|
| GENCODE-restricted (n=190): no significant 3′UTR difference | CDSand3UTR_GENCODEonly Panel D/F (translation-based + via-non-PTC-stop measures) confirm 3'UTR is inflated in PTC+ by including [PTC → natural stop] coding region; bias-free measure shows no real 3'UTR difference | ✓ direction matches |
| Expanded (n=1,166): NMD+/PTC+ median 3'UTR = 1,290 nt; PTC− = 800 nt; Control = 948 nt; p<0.001 | not in current inventory at n=1,166 scope; the bias-free 3'UTR measure (`utr3_via_ref_stop_nt`, added 2026-06-13) was wired into Figure 4 Panel D but not exported as a numeric panel | UNVERIFIED — see Gap-§4-C |

**Gap-§4-C — RESOLVED (already addressed by existing SF):** Panel F of the existing `CDSand3UTR_GENCODEonly/` SF already carries the 3'UTR-at-n=1,166 distribution using the bias-corrected non-PTC-stop measure. The descriptive table (`panelF_utr3_non_ptc_stop_refaug_descriptives.tsv`) confirms each manuscript-cited median exactly:
- ✓ NMD+/PTC+ median = 1,289.5 nt (manuscript: 1,290 — rounds to 1,290)
- ✓ NMD+/PTC− retained median = 799.5 nt (manuscript: 800 — rounds to 800)
- ✓ Control median = 948 nt (manuscript: 948 — exact)
- ✓ Pairwise Wilcoxon p = 3.07×10⁻⁴ (PTC− vs PTC+), 0.249 (PTC− vs Control, not significant), 1.13×10⁻⁸ (PTC+ vs Control)

No new SF needed. The original reconciliation note missed this because Panel F is a sub-panel of a larger 6-panel SF rather than a standalone supplement. Manuscript §4¶9 SF callout can point at "`Supplemental Figure CDSand3UTR_GENCODEonly` Panel F".

---

## 2. §5 reconciliation (deep-learning model)

### §5¶1 — Architecture and inputs — ✓ exact match
All seven design parameters (K=5, 500-nt windows, 9 channels per position, 5 structural features per ORF, 64-dim per-ORF embedding, ORFik+TD2 priority selection) match the Figure 5 Panel A schematic and the Methods §13 description.

### §5¶2–§5¶4 — Splits / sweep / training — ✓ matches
- Splits: val chr 2/4, test chr 1/3/5/7, paralog-free — matches `define-test-population` chunk and Panel B legend.
- Sweep: 3×4 = 12 configurations; best = ATG=500, stop=500 — matches Panel B legend.
- Training: BCEWithLogitsLoss, Adam lr=1e-3, batch 256, fp16, ReduceLROnPlateau (factor 0.5, patience 5), early stop patience 10, max 100 epochs, seed 42 — matches Methods §13.

### §5¶5 — Test performance — ✓ rounding consistent
| Claim | Figure-side | Verdict |
|---|---|---|
| AUC = 0.93 | Fig 5 Panel B: 0.9306 → 0.93 | ✓ |
| AUPRC = 0.83 | Panel B: 0.8330 → 0.83 | ✓ |

### §5¶6 — KernelSHAP branch decomposition — ⚠ "two-thirds" vs 60.7%; STOP/START ratio OK
| Claim | Figure-side | Verdict |
|---|---|---|
| "roughly ⅔ of the predictive information came from the ORF structural data" | Fig 5 Panel C: Structural 60.65% | DRIFT — 60.7% is closer to "three-fifths" than "two-thirds" |
| "STOP sequence ≈ 3× more important than START" | Panel C: Stop 28.82% / ATG 10.53% = 2.74× | ✓ rounds to ≈ 3× |
| "Structural 60.7%, Stop 28.8%, AUG 10.5%" (if these exact numbers are in the prose) | Panel C exact: 60.65 / 28.82 / 10.53 | ✓ |

**Find/replace §5¶6 — "⅔" overstates the structural fraction:**
```
find:    roughly two-thirds of the predictive information came from the ORF structural data
replace: roughly three-fifths of the predictive information came from the ORF structural data (60.7%)
```
Or, alternatively (more conservative, just adds the precise number):
```
find:    roughly two-thirds of the predictive information came from the ORF structural data
replace: roughly 60% of the predictive information came from the ORF structural data
```

### §5¶7 — Structural feature importance — ✓ matches
| Claim | Figure-side | Verdict |
|---|---|---|
| Downstream EJC count dominance | Panel D: mean |SHAP| = 2.1527 (top feature) | ✓ |
| ~15× the next feature | Panel D: 2.153 / 0.138 = 15.6× (next feature: `is_sqanti_cds`) | ✓ |

### §5¶8 — Start codon / Kozak motif recovery — ✓ direction matches
Direction (G/A at −3, C at −1, G at +4) matches Panel E's signed SHAP attribution. Panel E is currently a legacy-render placeholder pending Explorer signed-SHAP TSV (cluster expected back ~2026-06-19); when re-rendered natively, verify the dominant positive letters at those positions.

### §5¶9 — Stop codon (UGA enrichment + readthrough U+4) — ✗ STALE NUMBER
| Claim | Figure-side | Verdict |
|---|---|---|
| **"UGA is more commonly used (56% vs 48%, p<0.001)"** | **StopCodonUsage SF: UGA 52.0% NMD vs 48.3% Control (post-bug-fix 2026-04-30)** | **✗ 56% is pre-bug-fix** |
| U at +4 → readthrough association | Panel F signed SHAP positions G at variable stop position; readthrough story is qualitative — see Gap-§5-A | ⚠ qualitative direction OK |

**~~Find/replace §5¶9 — UGA frequency~~ — RETRACTED.** I had the polarity inverted. The 52.0% / 48.3% number in the StopCodonUsage SF I built is from the **pre-patch** v5 HTML render (Apr 2 10:07; bug patch is Apr 30). The Rmd's post-patch narrative (line 1249) cites **55.9% vs 48.1%**, which rounds cleanly to the manuscript's "56% vs 48%" — so the manuscript wording is essentially correct. The StopCodonUsage SF legend is the file that needs to be updated; see §1.11 below.

**Gap-§5-A:** The manuscript claim "U at +4 associated with readthrough" needs a citation. The current StopCodonUsage SF documents UGA preference (which is the readthrough-prone stop codon) but does not directly evidence U+4 → readthrough. Either (a) cite a primary readthrough reference (e.g., Loughran 2014, Skuzeski 1991) for the U+4 mechanism, or (b) tighten the claim to "U at +4 is the dominant positive-SHAP nucleotide at the post-stop position, consistent with the known U+4 readthrough preference."

### §5¶10 — GC content — ✓ qualitative match
The "GC difference is most prominent after stop codon" claim is consistent with the Panel F context — the post-stop region of a PTC-containing isoform is exonic sequence converted to 3'UTR. No quantitative reconciliation needed.

### §5¶11 — Attention distribution — ✓ matches
| Claim | Figure-side | Verdict |
|---|---|---|
| Most attention on ORF0 | AttentionDistribution Panel A: rank-0 median ~0.74–0.76 for both classes | ✓ |
| Attention more broadly distributed for NMD than Control | AttentionDistribution Panel B: NMD entropy density shifted higher than Control | ✓ |

### §5¶12 — uORF upweighting in NMD+/PTC− at n=1,166 — ✓ matches
| Claim | Figure-side | Verdict |
|---|---|---|
| n = 1,166 gene-matched pairs | Fig 5 Panel G scope: 1,016 PTC+ + 95 PTC− = 1,111 NMD; 1,107 Control | ✓ |
| Model upweights uORFs in PTC− | Panel G: 43.2% PTC− isoforms place >5% attention on candidate uORF, vs 0.5%/0% PTC+/Control | ✓ |
| Mann-Whitney U p<10⁻¹⁰ | Panel G pairwise file confirms | ✓ |

### §5¶13 — START vs STOP importance stratified by PTC class — ✗ "3×" → "2×"
| Claim | Figure-side | Verdict |
|---|---|---|
| **"START window information is roughly three times more important in the NMD+/PTC− than in the NMD+/PTC+ isoforms"** | **PTCSubclassBranchSHAP SF: 18.1% (PTC−) / 9.0% (PTC+) = 2.01×** | **✗ 3× → 2×** |
| "overall predictive performance was substantially lower in NMD+/PTC− isoforms" | PTCSubclassPerformance SF: AUC drops 0.96 → 0.74 | ✓ |

**Find/replace §5¶13 — START window ratio (already documented in the SF README):**
```
find:    is roughly three times more important in the NMD+/PTC-
replace: is roughly two-fold higher in the NMD+/PTC-
```

### §5¶14 — Summary claims — ✓
- "3'UTR length not an independent predictor" — Panel D structural-feature breakdown puts `is_ref_cds` / `is_sqanti_cds` above 3'UTR-related features; the model bypasses 3'UTR by relying on downstream-EJC count. ✓
- "UGA + U+4 as readthrough signals" — reproduces stop-codon directional finding; see Gap-§5-A for U+4 citation.

---

## 3. Methods audit — drafted new / updated subsections

Methods subsections currently in the Doc have several `[insert]` placeholders and lack detail for new analyses added during the recent rounds of figure work. Below I draft / update **only** the subsections that touch §4 or §5 analyses.

Each draft is written to be inserted into the Methods section verbatim (with the existing subsection heading reused). The subsections are grouped by what they're replacing or adding.

### 3.1 UPDATE — Isoform Pairs Analysis (Isopair Package)

Replace the existing "Isoform Pairs Analysis" subsection with the text below. Adds: explicit 12-event enumeration, ref-AUG-traceable scope and its size, the three nested scopes used across §4 / §5 / supplements, and the `enumerateOrfs` ORF selection rule.

> **Isoform pairs analysis (Isopair package).** Isoform pairs were constructed using the Isopair R package (https://github.com/peter4244/Isopair) to link splicing events to NMD susceptibility within each protein-coding gene. For each gene with ≥3 expressed isoforms (≥5% of overall gene expression in DMSO or SMG1i, ≥5 reads in at least one sample) and at least one NMD-susceptible isoform, we identified a high-expressing non-NMD *reference* isoform and paired it with both an NMD-susceptible *comparator* and a non-NMD *control comparator* from the same gene (3,009 NMD pairs and 3,009 Control pairs across 3,009 genes; pop_BC scope).
>
> Splice events were enumerated into twelve mutually exclusive categories: skipped exon (SE), 5′ alternative splice site (A5SS), 3′ alternative splice site (A3SS), intron retention (IR), intron-retention 5′-difference (IR diff 5'), intron-retention 3′-difference (IR diff 3'), partial intron retention 5′ (Partial IR 5'), partial intron retention 3′ (Partial IR 3'), alternative transcription start site (Alt TSS), alternative transcription end site (Alt TES), missing internal exon block, and "other" residual. Categories were assigned by comparing splice-junction coordinates between the reference and comparator transcripts; multiple events per pair were allowed, with prevalence computed as the fraction of pairs containing at least one event of that category (Fig 3 Panel C).
>
> Three nested scopes are used downstream. **Scope §2a (all-3-ENST coding, n = 190):** all three isoforms of the pair are GENCODE-annotated with a CDS; supports Fig 3 Panel D / Panel E / Panel F and Fig 4 Panels A / B. **Scope §3a (ENST-reference + ref-AUG-traceable, n = 1,659 NMD / 1,286 Control before re-intersection):** the reference isoform is GENCODE-annotated, and the comparator isoform's transcript projects to the same reference AUG (verified by `Isopair::enumerateOrfs()` Kozak-aware scan of all candidate ORFs and a downstream agreement check against the upstream 05k pipeline's `utr5_features_refaug.rds`). **Scope §3b (re-intersected gene-matched, n = 1,166 NMD + 1,166 Control):** §3a additionally re-intersected to 1:1 gene-matched pairs after categorical filtering; this is the canonical n = 1,166 used in Fig 4 Panels C/D, Fig 5 Panel G, the TD2BiasEvidence supplement (broad-scope panels), the PTCSubclassBranchSHAP supplement, and the PTCSubclassPerformance supplement.
>
> Ref-AUG traceability is determined by projecting the GENCODE reference AUG coordinate onto the comparator transcript via splice-aware alignment; if the projected position is within the comparator's coordinates and the resulting ORF can be enumerated, the comparator is traceable. The `enumerateOrfs()` ORF selector enumerates all candidate ORFs and ranks by (1) reference-CDS match, (2) TD2-called CDS, (3) Kozak score with start codon plausibility — this priority is the same selector used to populate the deep-learning model's K=5 candidate ORFs (Methods, "Deep learning model").

### 3.2 UPDATE — PTC determination and reference-AUG tracing

Replace the existing "PTC (Premature Termination Codon) Determination" subsection. The current text doesn't formalize the Kozak scoring or the frameshift inference algorithm, and doesn't explain the TD2 bias remediation.

> **PTC determination.** A stop codon was classified as a premature termination codon (PTC) under the canonical 50-nucleotide rule: the stop is a PTC if it occurs >50 nt upstream of the last exon–exon junction (terminal EJC). Distance from the comparator stop codon to the terminal EJC is computed in transcript coordinates from the comparator splice graph (Fig 3 Panel D).
>
> For GENCODE-annotated pairs (n = 190 all-3-ENST coding scope), PTC status is assigned directly from the comparator's own GENCODE CDS annotation. For novel comparator isoforms lacking GENCODE annotation, two complementary CDS calls are considered: (1) the **TD2 (TransDecoder2)** CDS call as inherited from the SQANTI3 isoform-classification pipeline; (2) a **reference-AUG-projected** ORF derived by projecting the GENCODE reference isoform's AUG coordinate onto the comparator transcript via splice-aware coordinate transformation.
>
> **TD2 bias remediation.** TD2 ranks candidate ORFs by length and avoids ORFs whose first stop codon would create a PTC, which systematically eliminates the very signal NMD analyses depend on (TD2BiasEvidence supplement). In the occult-PTC subset (effectively_ptc ∩ original_ptc==FALSE; n = 492), 487 of 492 (99.0%) TD2-called CDSes are downstream of the reference AUG, and reference AUGs have stronger Kozak context than TD2 AUGs in 78.0% of pairs (paired Wilcoxon p = 2.6×10⁻³⁸). For §4 / §5 analyses that require an unbiased main-CDS call, we therefore use the reference-AUG-projected ORF rather than the TD2 call. The exception is uORF-mediated NMD (e.g., ATF4): because uORF detection operates on the 5′UTR upstream of the main CDS, it is unaffected by TD2's main-CDS bias and the TD2 CDS is acceptable.
>
> **Kozak scoring.** Position weight matrix Kozak scores were computed using the canonical PWM (Kozak 1987 / Hernández 2019 conventions): the −3 G/A, −2 C, −1 C, +4 G positions are scored against the comparator's start-context nucleotides; a score of ≥1 in the `Isopair::enumerateOrfs()` scale was used to flag ORFs with credible start contexts (used as one of three uORF-eligibility filters in §5).
>
> **Mechanism class derivation.** Where a PTC has been identified, the splicing event responsible is attributed by the `mechanism_class` helper (`figures/lib/mechanism_class.R`), which classifies each PTC-attributed event as frameshift, in-frame stop, or 3′ UTR splice based on whether the event preserves the reading frame, introduces an in-frame stop codon within the original reading frame, or extends into the 3'UTR (Fig 3 Panel F).

### 3.3 NEW — Mechanism class taxonomy (currently undefined in Methods)

Insert as a new Methods subsection, immediately after "PTC determination". The current Methods text has no formal definition of the three-way classification (frameshift / in-frame stop / 3'UTR splice) — Fig 3 Panel F uses these categories but only the figure legend describes them.

> **Mechanism class taxonomy.** For each PTC-attributed splicing event (n = 69 events across 72 PTC+ NMD pairs in the n = 190 scope; Fig 3 Panel F), the event was assigned to one of three exclusive mechanism classes based on its effect on the reading frame:
>
> - **Frameshift** (38/69 events, 55.1%): the event introduces a length change that is not a multiple of 3 nucleotides, shifting the reading frame and exposing a premature stop downstream;
> - **In-frame stop** (23/69, 33.3%): the event preserves the reading frame but introduces a new in-frame stop codon — typically by including sequence content (e.g., a retained intron, an alternative exon) whose first frame includes a stop;
> - **3′ UTR splice** (8/69, 11.6%): the event removes sequence between the natural stop codon and the original terminal EJC, repositioning the terminal EJC to be >50 nt downstream of the (unchanged) stop codon — a PTC by definition under the 50-nt rule, but caused by a 3'UTR-region rearrangement rather than a CDS event.
>
> The implementation lives in `figures/lib/mechanism_class.R` as a derived helper (no cached column; recomputed from splice-graph coordinates at each invocation).

### 3.4 UPDATE — Deep learning model: architecture

Append to the existing "Deep Learning Model – Architecture and Training" subsection. The current text omits CNN-specific architecture details (layer counts, kernel sizes, activation functions). Insert after the existing description of windows / channels / structural features and before the training description.

> **Network architecture.** The model is composed of three input branches whose outputs are concatenated into a 64-dimensional per-ORF embedding and then aggregated across the five candidate ORFs via a learned softmax attention layer to produce a single transcript-level embedding for the classification head.
>
> The **AUG branch** and **stop branch** each process a 9 × 500 sequence-and-feature window through a shared-weight 1D convolutional stack: three convolutional blocks with kernel size 7, stride 1, and 32 / 64 / 128 output channels (each block followed by batch normalization and a ReLU activation, with max-pooling of kernel size 4 / 4 / 2 between blocks), followed by adaptive average pooling to a fixed-length vector and a single dense projection to 28 dimensions per branch. The two branches do not share weights with each other but do share weights across the five ORF slots within each branch.
>
> The **structural branch** is a per-ORF linear projection of the five hand-engineered structural features (`frac_start`, `frac_stop`, `is_ref_cds`, `is_sqanti_cds`, `n_downstream_ejc`) to an 8-dimensional embedding via a single fully-connected layer (no activation; bias enabled).
>
> The three branch outputs (28 + 28 + 8 = 64 dimensions) are concatenated per ORF. The five per-ORF embeddings are aggregated by a learned attention layer: a single dense projection from 64 → 1 produces per-ORF attention logits that are softmax-normalized across the five ORFs; the per-ORF embeddings are then weighted by the attention probabilities and summed to yield a single 64-dim transcript embedding. The transcript embedding feeds a two-layer MLP head (64 → 32 → 1) with a ReLU non-linearity between the two dense layers; the final scalar is the NMD logit.

(**Verify before insertion**: I drafted the per-branch CNN architecture (3 blocks × kernel-7 × 32/64/128 channels) and the MLP head dimensions from my best reading of the model design intent; the actual values are in `NMD_orf_model_v5_4ct/model.py`, which I have not been able to read this session since the cluster is down. Confirm the kernel sizes, channel counts, and head shape against `model.py` before applying. Mark any guessed values with `[verify from model.py]` rather than inserting them as published.)

### 3.5 UPDATE — Deep learning model: interpretability

Replace the existing "Deep Learning Model – Interpretability Methods" subsection. Adds the n's, the test-set restriction reasoning, the n=1,166 subgroup-stratified analyses, and the canonical KernelSHAP parameters.

> **Interpretability.** Three complementary attribution methods were used, each operating at a different layer of the model.
>
> **Joint DeepSHAP at the input layer.** DeepSHAP attributions were computed for the model's NMD-logit output with respect to the 9-channel input of the AUG and stop windows and the 5 structural features of the priority (rank-0) ORF, with ORFs at ranks 1–4 held fixed at their input values. Five independent replicate runs were used; each replicate used 500 random background transcripts drawn from the training set as the SHAP reference distribution. Per-position-per-channel attributions were averaged across replicates for downstream summarization (Fig 5 Panel D for structural-feature ranking; Figs 5 E/F for per-nucleotide signed attribution).
>
> **KernelSHAP at the fusion layer.** Exact Shapley decomposition across the three sub-encoders (AUG, stop, structural) was computed by running the model 2³ = 8 times per isoform with each combination of present/absent branches, using the trained model's expected branch output as the absent-branch value. This produces per-isoform branch-level attributions (`shap_atg`, `shap_stop`, `shap_structural`) that sum to the model's NMD-call signal. KernelSHAP was computed for the full cohort (`kernel_shap_branch_atg500_stop500_all.tsv`, n = 39,938 isoforms; not test-set restricted) and is the data source for Fig 5 Panel C and the PTCSubclassBranchSHAP supplement.
>
> **Attention weights at the aggregator.** Per-isoform, per-ORF attention probabilities were extracted from the softmax-normalized attention layer (Fig 5 Panel G, AttentionDistribution supplement).
>
> **Subgroup-stratified analyses (n = 1,166).** For the manuscript Subset 2 (n = 1,166 ref-AUG-traceable gene-matched pairs; Fig 4 C/D scope), three subgroups are defined by the upstream Isopair `category` field carried through to the figure-side TSV: NMD+/PTC+ retained (n = 1,016 with branch SHAP available), NMD+/PTC− retained (n = 95), and Control (n = 1,107). Branch SHAP descriptives at this scope are reported in the PTCSubclassBranchSHAP supplement (full-cohort KernelSHAP scope; not test-restricted). Predictive performance by subgroup (AUC, AUPRC, per-subgroup predicted probability) is reported in the PTCSubclassPerformance supplement at the test-set intersection (NMD+/PTC+ n = 255, NMD+/PTC− n = 30, Control n = 276), with AUC/AUPRC being the only analyses we report at test-set scope — all other §5 interpretability analyses use full cohort.

### 3.6 NEW — uORF rule (Fig 5 Panel G)

Insert as a new Methods subsection at the end of the deep-learning Methods block. The current Methods text does not define the uORF rule used in Fig 5 Panel G; the rule is only described in the figure legend.

> **uORF rule (Fig 5 Panel G).** A candidate ORF within a transcript was classified as a credible uORF for the attention analysis if it satisfied three structural criteria — purely upstream-side and independent of any main-CDS call (and therefore unaffected by TD2's PTC-avoidance bias):
>
> 1. **Length < 200 nt.** Below the 270-nt CDS floor of GENCODE-annotated protein-coding transcripts, so the rule cannot misclassify a real annotated CDS as a uORF.
> 2. **Position in 5′ half.** Both start and stop codons fall within the 5′ half of the transcript (`frac_start < 0.5` AND `frac_stop < 0.5`).
> 3. **Credible Kozak start.** Kozak score ≥ 1 in the `Isopair::enumerateOrfs()` scale.
>
> A transcript is reported as "placing >5% attention on a candidate uORF" if any of its five candidate ORF slots whose ORF satisfies the rule receives >5% of the model's softmax-normalized attention weight.

### 3.7 UPDATE — Stop-codon bug fix disclosure (currently absent)

Add this short paragraph to the deep-learning Methods (or to the Results §5 footnote, at Pete's preference). The previous version of the stop-codon analysis had an off-by-one bug; the corrected numbers (52.0% UGA in NMD vs 48.3% Control) are now canonical.

> **Stop-codon analysis correction.** An earlier off-by-one error in the ORFik stop-codon position column produced inflated UGA-in-NMD frequencies (~56%). After applying the patch `scripts/patch_stop_codon.py` (2026-04-30), within-class stop-codon frequencies on the chr-1/3/5/7 paralog-free test set are UGA 52.0% NMD vs 48.3% Control, UAA 27.1% vs 30.3%, UAG 20.9% vs 21.3%; UGA enrichment in NMD remains directionally consistent but is approximately half the magnitude originally reported. The StopCodonUsage supplement carries the corrected numbers; any earlier numbers (pre-2026-04-30) should not be used.

---

## 4. Summary of changes

### 4.1 Find/replace pairs for Google Doc (after v2 corrections)

**Three active find/replace pairs for the manuscript Doc:**

```
# §4¶8 — sample size off by 3
find:    NMD+/PTC- subset (n = 113)
replace: NMD+/PTC- subset (n = 116)

# §5¶6 — "two-thirds" overstates structural fraction (60.7%, closer to three-fifths)
find:    roughly two-thirds of the predictive information came from the ORF structural data
replace: roughly 60% of the predictive information came from the ORF structural data

# §5¶13 — START window ratio (already in PTCSubclassBranchSHAP SF README)
find:    is roughly three times more important in the NMD+/PTC-
replace: is roughly two-fold higher in the NMD+/PTC-
```

**Two figure-side fixes (changes to legend.md files, NOT manuscript Doc):**

```
# StopCodonUsage SF legend — replace pre-patch numbers with post-patch
find:    UGA is enriched in NMD (52.0% vs 48.3% in Control), UAA is depleted (27.1% vs 30.3%), and UAG is unchanged (20.9% vs 21.3%)
replace: UGA is enriched in NMD (55.9% vs 48.1% in Control); UAA and UAG are both depleted in NMD
# (then re-render UAA/UAG breakdown from post-patch Rmd when cluster is up)

# Significance-notation standardization across all four composite + SF legends (§1.10)
# Replace the existing per-legend asterisk strings with:
"Stars: ****p < 10⁻¹⁰, ***p < 10⁻⁴, **p < 10⁻³, *p < 0.05, n.s. otherwise."
```

**Retracted (do not apply):**
- ~~§4¶8 p-threshold find/replace~~ — misread the asterisk-key as data p-values
- ~~§5¶9 UGA frequency find/replace~~ — had the polarity inverted; the SF was stale, manuscript is correct

### 4.2 Methods inserts (drafted in §3 above; ready to paste into Google Doc)

| # | Heading | Action | Status |
|---|---|---|---|
| 3.1 | Isoform pairs analysis (Isopair package) | REPLACE | drafted; ready |
| 3.2 | PTC determination | REPLACE | drafted; ready |
| 3.3 | Mechanism class taxonomy | NEW | drafted; ready |
| 3.4 | Deep learning model — architecture | APPEND | drafted; **verify CNN block specs against `model.py`** when cluster up |
| 3.5 | Deep learning model — interpretability | REPLACE | drafted; ready |
| 3.6 | uORF rule (Fig 5 Panel G) | NEW | drafted; ready |
| 3.7 | Stop-codon analysis correction | NEW | drafted; ready |

### 4.3 Outstanding gaps (need decision before Methods can ship)

| Gap | Description | Recommended action |
|---|---|---|
| Gap-§4-A | §4¶1 descriptive numbers (median 7 isoforms/gene, 70%/75% expression) cite missing SFs | Move into Methods §3.1 as one or two sentences; drop SFx callouts |
| Gap-§4-B | §4¶2 transcript-length comparison cites missing `SFx - Transcript Length Comparison` | Decision: build the SF or move into Methods? Numbers (3,049 / 2,893 / 2,762 nt) are well-defined; an SF is more visible. **Recommend: build a 1-panel SF at pop_BC scope** |
| Gap-§4-C | §4¶9 n=1,166 3'UTR medians (1,290 / 800 / 948 nt) not in any current panel | Decision: build a 3'UTR-at-n=1,166 SF (using bias-free `utr3_via_ref_stop_nt` from Task #57)? **Recommend: yes, since §4 cites specific medians** |
| Gap-§5-A | §5¶9 U+4 readthrough claim lacks a citation | Add a primary readthrough citation (Loughran 2014 / Skuzeski 1991) OR soften wording to "consistent with known U+4 readthrough preference" |
| §3.4 verify | CNN layer counts / kernel sizes / activations in §3.4 are my best inference, not read from code | When cluster is back, read `NMD_orf_model_v5_4ct/model.py` and replace any `[verify]` placeholders |

### 4.4 Methods placeholders to fill (`[insert]` tokens; not §4/§5 but caught during the audit)

- Short-read RNA-seq: `[insert RIN value]`, `[insert minimum read depth]`, `[indicate one-pass versus two-pass mode]`
- Long-read RNA-seq: target FLNC depth `[fill in]`

These are out of scope for §4/§5 reconciliation but should be resolved before submission.

---

**Provenance.** All figure-side numbers in this report are from current files in `~/claude_projects/nmd/figures/multipanel/` and `~/claude_projects/nmd/figures/SupplementalFigures/` as of 2026-06-18. Manuscript text was read from the Google Doc on 2026-06-18; verify against the current Doc before applying any find/replace, since both sides are under active edit.
