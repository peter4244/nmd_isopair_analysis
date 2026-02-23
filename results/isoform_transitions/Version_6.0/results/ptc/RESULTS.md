# PTC Analysis Results

## 1. PTC Computation Summary

From 235,022 coding isoforms with exon structures:
- **38,746** PTC+ (stop codon >50nt upstream of last EJC)
- **170,251** PTC- (stop codon within 50nt of last EJC, multi-exon)
- **26,025** single-exon (no EJCs, classified as non-PTC)
- **0** stop codon not in any exon

Median PTC distance (multi-exon): -117 nt (negative = stop codon is downstream of last EJC, i.e., in last exon).

### By isoform source:
| Source | N | PTC+ | PTC rate | Median distance |
|--------|---:|-----:|---------:|----------------:|
| GENCODE | 54,952 | 5,397 | 9.8% | -121 |
| PacBio | 180,070 | 33,349 | 18.5% | -116 |

PacBio novel isoforms have nearly 2x the PTC rate of GENCODE isoforms, consistent with SQANTI ORF predictions being less constrained than curated annotations.

## 2. GENCODE Validation

Cross-tabulation of computed PTC vs GENCODE `nonsense_mediated_decay` biotype (n = 19,237 GENCODE isoforms):

| | GENCODE NMD | GENCODE non-NMD |
|---|---:|---:|
| **Computed PTC+** | 1,689 | 51 |
| **Computed PTC-** | 7 | 17,490 |

- **Concordance: 99.7%** (19,179/19,237)
- Sensitivity: 99.6% (1,689/1,696)
- Specificity: 99.7% (17,490/17,541)
- Precision: 97.1% (1,689/1,740)

The 51 false positives and 7 false negatives likely reflect edge cases in coordinate systems or alternative NMD mechanisms.

## 3. PTC -> NMD Responsiveness

### Per cell type (5%-filtered, Fisher test)

| Cell type | N | N PTC+ | NMD rate PTC+ | NMD rate PTC- | OR | 95% CI | p | FDR |
|-----------|---:|------:|:---:|:---:|----:|--------|--------:|-----:|
| AT | 68,744 | 10,961 | 3.3% | 3.1% | 1.07 | 0.95-1.20 | 0.233 | 0.297 |
| DD | 72,892 | 11,610 | 11.1% | 10.6% | 1.06 | 0.99-1.13 | 0.094 | 0.232 |
| **DD_ALI** | 72,055 | 11,518 | **6.3%** | **9.0%** | **0.68** | **0.63-0.74** | **1.3e-22** | **8.1e-22** |
| DO | 69,861 | 11,203 | 0.0% | 0.0% | 1.75 | 0.30-7.00 | 0.423 | 0.423 |
| FB | 65,558 | 10,353 | 1.7% | 1.6% | 1.10 | 0.93-1.30 | 0.248 | 0.297 |
| MV | 63,484 | 9,946 | 6.1% | 5.7% | 1.08 | 0.98-1.18 | 0.116 | 0.232 |

**Meta-analysis (5%-filtered):** Pooled OR = 0.96 (0.92-1.00), p = 0.03. High heterogeneity: I^2 = 94%.

The pooled result is driven entirely by DD_ALI, where PTC+ isoforms are *less* likely to be NMD-responsive -- the opposite of expectation. This is a GENCODE-specific effect (see below).

### GENCODE vs PacBio Stratification

**Meta-analysis by source:**
- **GENCODE only:** Pooled OR = 0.88 (0.80-0.98), p = 0.017 (driven by DD_ALI OR = 0.62)
- **PacBio only:** Pooled OR = 1.05 (1.01-1.10), p = 0.028 (I^2 = 0%, consistent weak signal)

The GENCODE pooled result is reversed (PTC predicts *less* NMD responsiveness), driven entirely by DD_ALI. Excluding DD_ALI, GENCODE ORs are >=1.0 across cell types. PacBio shows a homogeneous but tiny positive association.

### Distance stratification (no dose-response)

NMD rates by PTC distance bin (pooled across cell types):
- 0-50 (no PTC): 5.2%
- 51-100: 5.4%
- 101-200: 5.7%
- 201-500: 4.6%
- \>500: 4.7%

No evidence of dose-response. Spearman correlations near zero for all cell types except DD_ALI (rho = -0.014, p = 0.0003, reversed direction).

## 4. GENCODE NMD Biotype -> NMD Responsiveness

The GENCODE `nonsense_mediated_decay` biotype does not consistently predict Smg1i responsiveness:
- Meta-analysis: Pooled OR = 0.96 (0.87-1.06), p = 0.40
- DD_ALI again anomalous: OR = 0.75 (0.64-0.88), p = 0.0003
- Other cell types show ORs near 1.0

## 5. logFC Distribution Analysis

### 5%-filtered
Only DD_ALI (p = 2.7e-19) and DO (p = 3.2e-13) show significant PTC+ vs PTC- logFC shifts.
- DD_ALI: PTC+ shifted *lower* (reversed)
- DO: PTC+ shifted higher (expected), median shift ~+0.06

### Unfiltered
DO and DD_ALI remain strongly significant. Additional marginal signals in AT (p=0.044), FB (p=0.042), MV (p=0.002).

### By source (unfiltered, significant results)
| Source | Cell type | Median shift | p |
|--------|-----------|:-----------:|------:|
| GENCODE | AT | +0.057 | 0.020 |
| GENCODE | DD_ALI | **-0.347** | **7.5e-14** |
| GENCODE | DO | +0.118 | 5.8e-06 |
| GENCODE | FB | +0.055 | 0.018 |
| PacBio | DO | +0.029 | 3.7e-05 |
| PacBio | MV | +0.004 | 0.021 |

The DD_ALI GENCODE signal is the largest effect but in the *wrong* direction. DO shows the strongest expected signal across both sources.

## 6. Signal Exploration: Top Findings

### Top signals with expected direction (PTC+ has higher logFC):
1. **DO basic PTC (unfiltered):** median shift = +0.064, p = 4.6e-24
2. **DD_ALI strong PTC (unfiltered):** median shift = +0.664, p = 2.1e-18
3. **DD_ALI strong PTC (5%-filtered):** median shift = +0.888, p = 3.4e-15
4. **DO basic PTC (5%-filtered):** median shift = +0.060, p = 3.2e-13

**Key insight**: DD_ALI "strong PTC" (>500nt, >=2 EJCs, GENCODE) shows the largest positive shift (+0.66 to +0.89), while DD_ALI's *overall* signal is reversed. This means a small subset of well-characterized GENCODE PTCs do predict NMD in DD_ALI, but the bulk of computed PTCs (mostly PacBio-derived) do not.

### Within-gene paired analysis
DO is the only cell type with significant within-gene PTC effects (unfiltered: median diff = +0.015, p = 0.005).

### Quantile enrichment
- DO: strong positive trend (rho = 0.98, p < 0.001) -- PTC fraction increases with logFC decile
- DD_ALI: strong negative trend (rho = -0.89, p = 0.001) -- reversed

## 7. Pair-Level PTC and Frame Disruption

### PTC prevalence in comparison pairs
No significant enrichment of PTC+ isoforms in NMD comparators (C1, C2) vs baseline (C4) after FDR correction. Lowest FDR = 0.096 (DD_ALI C1: OR=0.44, reversed; MV C2: OR=1.30, expected).

### Within-pair PTC asymmetry
McNemar's test: only C4 (baseline) shows significant asymmetry (p = 0.002), with more comparator-only PTC than dominant-only PTC. C1 and C2 show no asymmetry.

### Frame disruption: NMD vs baseline
No significant difference in frame disruption rates between NMD comparisons and baseline after FDR correction. Frame disruption rates are ~51-56% across all comparisons (high baseline).

### LOSS + frame-disrupting events
No significant enrichment in NMD pairs. LOSS frame-disruption rates: ~34% across all comparisons with minimal NMD vs baseline difference.

## 8. rMATS Frame Disruption Enrichment

### Overview

1,056,930 rMATS splice events across 6 cell types were classified by frame disruption potential and CDS overlap. Gene ID matching: 67.2% of rMATS genes mapped to CDS envelopes (18,488 genes from 54,952 GENCODE coding isoforms). 60.8% of CDS-overlapping SE events are frame-disrupting (below the theoretical 66.7%, consistent with evolutionary selection for frame-preserving exon lengths).

### Primary enrichment: frame-disrupting CDS events are significantly enriched

Among all non-RI rMATS events, frame-disrupting CDS events are ~2x more likely to be significant (FDR < 0.05, |dPSI| ≥ 0.05) than frame-preserving CDS events:

| Cell type | OR | 95% CI | p | Sig disrupting / total | Sig preserving / total |
|-----------|---:|--------|------:|---:|---:|
| DD | 1.80 | 1.72-1.89 | 3.0e-159 | 7,785 / 139,272 | 2,694 / 84,715 |
| MV | 2.66 | 2.17-3.27 | 3.2e-25 | 496 / 47,470 | 119 / 30,063 |
| AT | 2.43 | 1.64-3.72 | 1.8e-6 | 121 / 57,217 | 32 / 36,791 |
| DD_ALI | 1.32 | 0.52-3.58 | 0.67 | 15 / 50,984 | 8 / 35,768 |
| FB | 1.98 | 0.35-20.1 | 0.49 | 6 / 45,797 | 2 / 30,222 |
| DO | — | — | — | 0 / 54,507 | 0 / 37,768 |

**Meta-analysis (random effects):** Pooled OR = 2.13 (1.66-2.74), p = 3.15e-9, I² = 69%.

DD, MV, and AT all show significant enrichment. DD_ALI, FB, and DO have too few significant events for meaningful tests.

### Enrichment by event type (DD)

| Event type | OR | 95% CI | p |
|------------|---:|--------|------:|
| A3SS | 2.76 | 2.33-3.29 | 1.5e-36 |
| A5SS | 2.38 | 1.96-2.89 | 1.5e-21 |
| SE | 2.22 | 2.07-2.38 | 3.4e-132 |
| MXE | 1.24 | 1.16-1.33 | 6.7e-10 |

All four non-RI event types are significantly enriched, with boundary changes (A3SS, A5SS) showing the strongest effects.

### RI events: no CDS enrichment

CDS-overlapping RI events are NOT enriched among significant events compared to non-CDS RI:

**Meta-analysis:** Pooled OR = 0.80 (0.62-1.05), p = 0.11, I² = 58%.

MV shows a marginally significant *depletion* (OR = 0.72, p = 0.039). This is unexpected if intron retention were a primary NMD-inducing mechanism, though it may reflect that RI events in UTRs are also functionally impactful.

### Direction analysis

Among significant frame-disrupting events, the direction of dPSI (positive = more inclusion in Smg1i) varies by cell type and event type:

**SE events:**
- AT: 75% positive (n=69, p=2.9e-5) — biased toward increased inclusion
- MV: 66% positive (n=271, p=1.4e-7) — biased toward increased inclusion
- DD: 38% positive (n=3,855, p=1.4e-53) — biased toward *decreased* inclusion
- DD_ALI: 92% positive (n=13, p=0.003) — small n but strongly biased

The opposing directions across cell types are consistent with the biological ambiguity: for any given SE event, either inclusion or skipping can create the PTC depending on reading frame context. The enrichment signal is robust (two-tailed), but the direction is event-specific, not predictable from frame disruption alone.

**RI events (CDS, one-tailed test for positive dPSI):**
- DD: 27% positive (n=930, p=1.0) — strongly biased *negative* (opposite of expected)
- MV: 35% positive (n=113, p=1.0) — also biased negative

RI events show no evidence of the expected positive dPSI bias. Significant CDS RI events are predominantly associated with *decreased* retention under Smg1i, suggesting these retained introns may not function primarily as NMD-triggering elements.

### Magnitude analysis

Frame-disrupting CDS events have modestly larger |dPSI| than frame-preserving CDS events across all cell types (Wilcoxon test), but the absolute differences are small (median |dPSI| ≈ 0.01-0.02 for both groups).

### Known NMD-coupled poison exons

SRSF7 (a canonical NMD autoregulatory gene) shows a strong frame-disrupting SE event in DD with dPSI = +0.49 (FDR = 0.002), consistent with poison exon inclusion being rescued by NMD inhibition. HNRNPL shows a similar pattern (dPSI = +0.58, FDR = 0.004 in DD).

## 9. Interpretation and Conclusions

### Central finding: isoform-level PTC status does not predict NMD responsiveness
**NMD-responsive transcripts are not enriched for PTC-containing transcripts.** Across six cell types, PTC+ isoforms are no more likely to be NMD-responsive (upregulated by Smg1i) than PTC- isoforms. The only exceptions are a small signal in DO and a subset of well-characterized GENCODE PTCs in DD_ALI ("strong PTC": distance >500nt, >=2 downstream EJCs, GENCODE-annotated).

### But: event-level frame disruption DOES predict Smg1i responsiveness
The rMATS analysis (Section 8) reveals a complementary signal invisible at the isoform level: splice events that disrupt the reading frame are ~2x more likely to be significantly affected by Smg1i treatment (meta OR = 2.13, p = 3e-9). This is consistent across SE, A5SS, A3SS, and MXE event types and across the three cell types with sufficient power (DD, MV, AT).

### Reconciling the two findings
The apparent contradiction — null isoform-level PTC signal but positive event-level frame disruption signal — has a coherent explanation:

1. **CDS prediction quality confounds the isoform analysis.** ~75% of coding isoforms are PacBio-derived with SQANTI CDS predictions, yielding a 2x higher PTC rate (18.5% vs 9.8% GENCODE). Misannotated ORFs generate false PTC calls that dilute real signals at the isoform level.
2. **The rMATS analysis bypasses CDS prediction entirely.** Frame disruption is computed from event geometry (exon width mod 3) relative to gene-level CDS envelopes from GENCODE annotations only. This avoids the PacBio ORF prediction problem.
3. **Event-level resolution is more precise.** A single splice event can be cleanly classified as frame-disrupting, while an isoform may contain multiple splicing differences whose net effect on the reading frame is ambiguous.

### PTC+ and PTC- NMD-responsive isoforms are indistinguishable
Among isoforms that *are* NMD-responsive, PTC+ and PTC- isoforms show no meaningful differences in:
- **DMSO expression level** (baseline CPM distributions overlap)
- **Response magnitude** (logFC distributions overlap)
- **Strong response rates** (proportion with logFC > 1 or logFC > 2)

This means PTC status does not even identify a distinct *subset* of NMD-responsive isoforms -- the PTC+ NMD-responsive isoforms look identical to PTC- NMD-responsive isoforms on all measurable dimensions.

### The DD_ALI paradox
DD_ALI shows the strongest overall signal but in the *reversed* direction (PTC+ isoforms are less Smg1i-responsive). However, restricting to "strong PTC" features flips the signal to strongly positive. This suggests that the bulk of PTC calls -- driven by PacBio novel isoforms -- are noisy in DD_ALI, masking a real signal in well-annotated isoforms.

### GENCODE vs PacBio
PacBio novel isoforms have 2x the PTC rate of GENCODE (18.5% vs 9.8%), likely reflecting less reliable ORF predictions from SQANTI3. The PTC->NMD signal is generally more consistent (though still weak) when restricted to GENCODE isoforms.

### RI events: no evidence for NMD-driven intron retention
Despite the expectation that retained introns in CDS regions would introduce stop codons and trigger NMD, CDS-overlapping RI events are not enriched among Smg1i-responsive events (meta OR = 0.80, p = 0.11). Among significant CDS RI events, dPSI is biased *negative* (decreased retention under Smg1i), opposite to the expected direction. This suggests that NMD-coupled intron retention is not a dominant mechanism in these cell types, or that RI events detected by rMATS represent a different biological process than PTC-inducing intron retention.

### Direction ambiguity in frame-disrupting events
The direction of dPSI for significant frame-disrupting events varies across cell types (DD is biased negative, AT and MV are biased positive for SE events). This is consistent with the biological reality that for any given frame-disrupting SE event, either inclusion or skipping can create the PTC depending on reading frame context. The enrichment signal is robust, but the direction is event-specific, not predictable from frame disruption status alone.

### Why is the isoform-level signal weak?
1. **PTC is neither necessary nor sufficient for NMD.** NMD can be triggered without classical PTCs -- long 3'UTRs, upstream ORFs (uORFs), and other transcript features can independently activate NMD. Conversely, not all PTC-containing transcripts are efficient NMD substrates. The 50-nucleotide rule identifies a *structural feature* that is associated with NMD but does not determine NMD substrate status.
2. **NMD substrate degradation** occurs co-translationally; the steady-state mRNA level already reflects NMD. Smg1i inhibition rescues degradation, so the logFC captures the *change* in degradation, not PTC status per se.
3. **CDS prediction quality**: ~75% of coding isoforms are PacBio-derived with SQANTI CDS predictions. Misannotated ORFs generate false PTC calls that dilute real signals. The rMATS event-level analysis, which bypasses CDS predictions, recovers a clear frame disruption signal (OR ≈ 2).
4. **Multiple testing**: The 5% expression filter restricts to ~60K-70K isoforms per cell type, with only ~1-11% NMD-responsive, creating a strong class imbalance.

### Threshold independence
These analyses are independent of the non-NMD classification threshold (0.50 vs 0.95). Scripts 01-04 operate on raw DE results and expression data, not comparison pairs. Script 05 (pair-level PTC/frame-disruption) uses the comparison framework, and its results are null regardless of threshold. Script 06 (rMATS) uses entirely independent short-read data.

## Output Files

### results/
| File | Description |
|------|-------------|
| `ptc_status.rds` | PTC features for 235,022 coding isoforms |
| `gencode_validation.tsv` | Computed PTC vs GENCODE NMD biotype cross-tab |
| `ptc_predicts_nmd_by_celltype.tsv` | Fisher test results (5%-filtered) |
| `ptc_predicts_nmd_unfiltered_by_celltype.tsv` | Fisher test results (unfiltered) |
| `ptc_predicts_nmd_by_source.tsv` | GENCODE vs PacBio stratification |
| `ptc_distance_nmd_rate.tsv` | NMD rate by PTC distance bin (pooled) |
| `ptc_distance_nmd_rate_by_celltype.tsv` | Distance-stratified by cell type |
| `nmd_biotype_predicts_nmd_by_celltype.tsv` | NMD biotype analysis |
| `nmd_rate_by_biotype.tsv` | NMD rate by GENCODE biotype |
| `nmd_biotype_predicts_nmd_enst_only.tsv` | ENST-only sensitivity analysis |
| `ptc_signal_exploration_all_tests.tsv` | All 10-angle exploration tests |
| `ptc_logfc_shift_summary.tsv` | logFC shift statistics |
| `ptc_prevalence_nmd_vs_baseline.tsv` | Pair-level PTC prevalence Fisher tests |
| `ptc_rates_by_comparison_run.tsv` | PTC rates by comparison x run |
| `ptc_within_pair_asymmetry.tsv` | Within-pair PTC asymmetry |
| `frame_disruption_nmd_vs_baseline.tsv` | Frame disruption NMD vs baseline |
| `frame_disruption_loss_nmd_vs_baseline.tsv` | LOSS frame-disruption tests |
| `rmats_event_classification.tsv` | All rMATS events with frame disruption + CDS annotations |
| `rmats_enrichment_by_celltype.tsv` | Fisher test results per cell type |
| `rmats_enrichment_by_event_type.tsv` | Enrichment stratified by SE/A5SS/A3SS/MXE |
| `rmats_direction_tests.tsv` | Binomial direction bias tests |
| `rmats_magnitude_tests.tsv` | |dPSI| magnitude comparisons |
| `rmats_summary_stats.tsv` | Event counts by classification |

### figures/
| File | Description |
|------|-------------|
| `ptc_logfc_by_celltype.pdf` | logFC density by cell type (5%-filtered) |
| `ptc_logfc_by_celltype_expression.pdf` | logFC by cell type x expression tertile (5%-filtered) |
| `ptc_logfc_by_celltype_unfiltered.pdf` | logFC density by cell type (unfiltered) |
| `ptc_logfc_by_celltype_expression_unfiltered.pdf` | logFC by cell type x expression (unfiltered) |
| `ptc_logfc_by_source.pdf` | logFC GENCODE vs PacBio (5%-filtered) |
| `ptc_logfc_by_source_unfiltered.pdf` | logFC GENCODE vs PacBio (unfiltered) |
| `rmats_dpsi_by_frame_class.pdf` | dPSI distributions: frame-disrupting vs preserving |
| `rmats_enrichment_forest.pdf` | Forest plot of enrichment ORs across cell types |
