# PTC Analysis Results

## Overview

Two complementary approaches test whether reading-frame disruption predicts NMD responsiveness (Smg1i upregulation):

1. **Event-level (rMATS, short-read):** Frame-disrupting splice events are ~2x enriched among Smg1i-responsive events (meta OR = 2.13, p = 3e-9). This signal is robust across event types and cell types with adequate power.

2. **Isoform-level (PTC status, long-read):** PTC+ isoforms are NOT enriched among NMD-responsive transcripts, even when restricted to GENCODE-annotated isoforms with validated CDS. This null holds across all six cell types.

The discrepancy reflects a fundamental difference in what these analyses measure, not a technical artifact of CDS annotation quality. Canonical case studies (SRSF1, SRSF7) confirm that rMATS detects poison exon splicing changes that are invisible at the isoform level due to PacBio's limited ability to quantify lowly-expressed NMD substrates.

---

## 1. rMATS Frame Disruption Enrichment

### 1.1 Event classification

1,056,930 rMATS splice events (SE, A5SS, A3SS, MXE, RI) across 6 cell types were classified by frame disruption potential (exon width mod 3) and CDS overlap (gene-level GENCODE CDS envelope). 67.2% of rMATS genes matched CDS envelopes. 60.8% of CDS-overlapping SE events are frame-disrupting (slightly below the theoretical 66.7%, consistent with evolutionary selection for frame-preserving exon lengths).

### 1.2 Primary result: frame-disrupting enrichment

Among non-RI events, frame-disrupting CDS events are ~2x more likely to be significant (FDR < 0.05, |dPSI| >= 0.05) than frame-preserving CDS events:

| Cell type | OR | 95% CI | p | Sig disrupting / total | Sig preserving / total |
|-----------|---:|--------|------:|---:|---:|
| DD | 1.80 | 1.72-1.89 | 3.0e-159 | 7,785 / 139,272 | 2,694 / 84,715 |
| MV | 2.66 | 2.17-3.27 | 3.2e-25 | 496 / 47,470 | 119 / 30,063 |
| AT | 2.43 | 1.64-3.72 | 1.8e-6 | 121 / 57,217 | 32 / 36,791 |
| DD_ALI | 1.32 | 0.52-3.58 | 0.67 | 15 / 50,984 | 8 / 35,768 |
| FB | 1.98 | 0.35-20.1 | 0.49 | 6 / 45,797 | 2 / 30,222 |
| DO | -- | -- | -- | 0 / 54,507 | 0 / 37,768 |

**Meta-analysis (random effects):** Pooled OR = 2.13 (1.66-2.74), p = 3.15e-9, I^2 = 69%.

DD, MV, and AT all show significant enrichment. DD_ALI, FB, and DO have too few significant events for meaningful tests.

### 1.3 Enrichment by event type (DD)

| Event type | OR | 95% CI | p |
|------------|---:|--------|------:|
| A3SS | 2.76 | 2.33-3.29 | 1.5e-36 |
| A5SS | 2.38 | 1.96-2.89 | 1.5e-21 |
| SE | 2.22 | 2.07-2.38 | 3.4e-132 |
| MXE | 1.24 | 1.16-1.33 | 6.7e-10 |

All four non-RI event types are significantly enriched. Boundary changes (A3SS, A5SS) show the strongest effects.

### 1.4 RI events: no CDS enrichment

CDS-overlapping RI events are NOT enriched among significant events vs non-CDS RI. Meta-analysis: OR = 0.80 (0.62-1.05), p = 0.11, I^2 = 58%. Among significant CDS RI events, dPSI is biased negative (DD: 27% positive, MV: 35% positive), opposite to the expected positive direction for NMD-coupled intron retention.

### 1.5 Direction analysis

Among significant frame-disrupting SE events, the direction of dPSI varies by cell type:
- AT: 75% positive (n=69, p=2.9e-5)
- MV: 66% positive (n=271, p=1.4e-7)
- DD: 38% positive (n=3,855, p=1.4e-53)
- DD_ALI: 92% positive (n=13, p=0.003)

The opposing directions reflect biological ambiguity: for any given frame-disrupting SE event, either inclusion or skipping can create the PTC depending on reading frame context. The enrichment is robust regardless of direction.

### 1.6 Magnitude analysis

Frame-disrupting CDS events have modestly larger |dPSI| than frame-preserving CDS events (Wilcoxon, all cell types significant), but absolute differences are small (median |dPSI| ~ 0.01-0.02 for both groups).

### 1.7 Known NMD-coupled poison exons

Multiple canonical NMD autoregulatory genes show strong rMATS signals:
- **SRSF1**: 3 significant frame-disrupting SE events in DD (dPSI ~ -0.31 to -0.33) at the classic poison exon locus (chr17:58005398-58005600)
- **SRSF6**: SE events in DD (dPSI = +0.65) and MV (dPSI = +0.63) -- strong inclusion increase
- **SRSF7**: SE event in DD (dPSI = +0.49, FDR = 0.002) -- poison exon inclusion rescued
- **HNRNPL**: 11 significant frame-disrupting events including SE (dPSI = +0.58) and MXE (dPSI = +0.54) in DD
- **PTBP1**: 12 significant frame-disrupting events across DD and MV
- **SRSF5**: 15 significant frame-disrupting events across A3SS, A5SS, MXE, SE in DD

---

## 2. Isoform-Level PTC Analysis (Summary)

### 2.1 PTC computation

From 235,022 coding isoforms: 38,746 PTC+ (16.5%), 170,251 PTC- (multi-exon), 26,025 single-exon.

| Source | N | PTC+ | PTC rate |
|--------|---:|-----:|---------:|
| GENCODE | 54,952 | 5,397 | 9.8% |
| PacBio | 180,070 | 33,349 | 18.5% |

PTC computation validated against GENCODE NMD biotype: 99.7% concordance (19,179/19,237).

### 2.2 PTC does not predict NMD responsiveness

| Cell type | NMD rate PTC+ | NMD rate PTC- | OR | p |
|-----------|:---:|:---:|----:|------:|
| AT | 3.3% | 3.1% | 1.07 | 0.233 |
| DD | 11.1% | 10.6% | 1.06 | 0.094 |
| **DD_ALI** | **6.3%** | **9.0%** | **0.68** | **1.3e-22** |
| DO | 0.0% | 0.0% | 1.75 | 0.423 |
| FB | 1.7% | 1.6% | 1.10 | 0.248 |
| MV | 6.1% | 5.7% | 1.08 | 0.116 |

Meta-analysis: OR = 0.96 (0.92-1.00), p = 0.03, I^2 = 94%. The only significant cell type is DD_ALI, which goes in the *wrong* direction (PTC+ less NMD-responsive).

### 2.3 GENCODE-only analysis: still null

Restricting to GENCODE isoforms (validated CDS) does not rescue the signal:
- **GENCODE meta-analysis:** OR = 0.88 (0.80-0.98), p = 0.017, driven by DD_ALI (OR = 0.62)
- Excluding DD_ALI, GENCODE ORs are all non-significant and near 1.0

This is a critical result: if the isoform-level null were simply due to noisy PacBio CDS predictions, GENCODE-only analysis should show strong PTC enrichment. It does not.

### 2.4 GENCODE NMD biotype: also null

GENCODE's own `nonsense_mediated_decay` biotype annotation does not predict Smg1i responsiveness (meta OR = 0.96, p = 0.40).

### 2.5 No dose-response

NMD rates are flat across PTC distance bins (0-50: 5.2%, 51-100: 5.4%, 101-200: 5.7%, 201-500: 4.6%, >500: 4.7%). No dose-response relationship.

### 2.6 logFC distributions

Only DD_ALI and DO show significant PTC+ vs PTC- logFC shifts. DD_ALI is in the wrong direction; DO shows a small expected signal (median shift +0.06). "Strong PTC" features (>500nt, >=2 EJCs, GENCODE) show a positive shift in DD_ALI (+0.66 to +0.89), representing a small well-characterized subset.

### 2.7 Pair-level analysis

No significant PTC enrichment in NMD comparison pairs (C1, C2) vs baseline (C4). No significant frame disruption differences between NMD and baseline comparisons.

---

## 3. Bridging the Two Analyses: Case Studies

### 3.1 SRSF1 -- the classic poison exon

SRSF1 auto-regulates via a well-characterized poison exon whose inclusion introduces a PTC, targeting the transcript for NMD.

**rMATS (short-read):** 3 significant frame-disrupting SE events in DD at chr17:58005398-58005600 (the known poison exon locus). dPSI ~ -0.31 to -0.33 (FDR ~ 0.02), meaning the poison exon is less included under Smg1i. This is consistent with NMD rescue: when NMD is inhibited, the poison exon-containing transcript accumulates, providing negative feedback that reduces poison exon inclusion.

**Isoform-level (PacBio):** Only a single SRSF1 isoform passes expression filters across all 6 cell types: ENST00000583741.1, annotated as `nonsense_mediated_decay` biotype with PTC+ (distance = 389 nt). This isoform is **not significant in any cell type** (adj.P.Val >> 0.05). The logFC is actually negative in 5 of 6 cell types.

**Interpretation:** PacBio long-read sequencing lacks the depth to quantify the poison exon-containing transcript, which is being actively degraded by NMD. The splicing change is clearly visible at the event level (rMATS) but invisible at the isoform level (PacBio DGE).

### 3.2 SRSF7

**rMATS:** Significant SE event in DD (dPSI = +0.49, FDR = 0.002) and MV (dPSI = +0.32) -- poison exon inclusion increases under Smg1i, consistent with NMD rescue.

**Isoform-level:** 6 isoforms per cell type, only 1 significant hit in DD_ALI (ENST00000477635.5, logFC = +2.89), but PTC status is NA for this isoform.

### 3.3 Linking rMATS genes to isoform-level PTC

As a direct bridge test: among DD GENCODE isoforms in genes with significant rMATS frame-disrupting events (3,237 genes), PTC+ NMD rate = 8.7% vs PTC- NMD rate = 7.6% (Fisher OR = 1.16, p = 0.28). Even restricting to the genes where rMATS independently confirms frame-disrupting splicing, isoform-level PTC status still does not predict NMD responsiveness.

---

## 4. Interpretation and Conclusions

### The two analyses test fundamentally different things

The rMATS and isoform-level analyses are not the same question measured with different precision. They test genuinely different biological quantities:

- **rMATS** asks: does the *usage* of frame-disrupting exons shift under NMD inhibition? This measures splicing changes at individual events, using short-read junction evidence with replicate-aware statistics.
- **Isoform-level PTC** asks: are full-length transcripts with annotated PTCs upregulated under NMD inhibition? This measures transcript-level abundance changes from PacBio long-read DGE.

### Why rMATS detects a signal that isoform-level analysis does not

1. **Sensitivity to lowly-expressed NMD substrates.** NMD substrates are being actively degraded and are therefore lowly expressed. PacBio long-read sequencing has limited depth; NMD target isoforms may not accumulate enough reads for quantification or statistical testing. The SRSF1 case demonstrates this directly: the known NMD target isoform is detectable but not statistically significant in any cell type. Short-read rMATS, with much higher sequencing depth, detects the splicing shift without needing to quantify full-length isoforms.

2. **Event-level resolution avoids isoform model complexity.** rMATS classifies individual splice events by a simple geometric criterion (exon width mod 3). The isoform-level analysis requires correct CDS annotation for the entire transcript, correct exon-exon junction identification, and sufficient expression for statistical power -- each adding potential failure points.

3. **The rMATS signal may partly reflect non-NMD biology.** SMG1 phosphorylates UPF1, which participates in pathways beyond classical NMD (Staufen-mediated decay, histone mRNA turnover, other RNA surveillance). Frame-disrupting exons may be enriched among Smg1i-responsive events because they engage UPF1-dependent pathways more broadly, not exclusively through PTC-triggered NMD.

### The CDS quality explanation is insufficient

The original hypothesis -- that noisy PacBio CDS predictions mask a real PTC signal -- is directly contradicted by the GENCODE-only analysis:

- GENCODE isoforms have validated CDS annotations (99.7% concordance with GENCODE NMD biotype)
- Yet GENCODE-only meta-analysis shows OR = 0.88 (PTC predicts *less* NMD responsiveness)
- Even GENCODE's own `nonsense_mediated_decay` biotype annotation shows OR = 0.96 (p = 0.40)

PacBio CDS quality may add noise, but it is not the primary explanation for the null isoform-level result. The null holds even with the best available CDS annotations.

### PTC is not a useful predictor of Smg1i responsiveness

Across all analysis angles, isoform-level PTC status fails to predict NMD responsiveness:
- **Overall:** OR near 1.0 across 5 of 6 cell types
- **GENCODE-only:** OR = 0.88 (wrong direction)
- **GENCODE NMD biotype:** OR = 0.96 (null)
- **No dose-response** with PTC distance
- **No enrichment** in NMD comparison pairs vs baseline
- **PTC+ and PTC- NMD-responsive isoforms are indistinguishable** in expression level, response magnitude, and strong response rates

This does not mean PTCs are biologically irrelevant to NMD -- it means that static PTC annotation does not capture the dynamic, context-dependent process of NMD substrate recognition as measured by Smg1i treatment.

### The frame disruption enrichment is real but its mechanism is uncertain

The rMATS OR ~ 2 enrichment for frame-disrupting events is statistically robust, but several features complicate a pure NMD interpretation:

- **Direction ambiguity:** In DD (highest power), significant frame-disrupting SE events are biased toward *decreased* inclusion (62% negative dPSI), while AT and MV show the opposite. A pure NMD model does not predict this heterogeneity.
- **RI events show no CDS enrichment** and have negative dPSI bias, contradicting the expectation that intron retention introduces PTCs rescued by Smg1i.
- **The enrichment may partly reflect non-NMD UPF1-dependent pathways** or evolutionary selection on exon lengths that correlates with regulatory complexity.

### DD dominance

DD (4 replicates, differentiated cells) has ~10x more significant events than other cell types (3 replicates each). All per-event-type rMATS results are essentially DD-specific. The meta-analysis provides formal pooling, but DD drives the estimate.

### DD_ALI paradox

DD_ALI consistently shows the strongest isoform-level signal but in the reversed direction (PTC+ less NMD-responsive). This holds for GENCODE-only (OR = 0.62) and GENCODE NMD biotype (OR = 0.75). The reversal is unexplained and may reflect cell-type-specific biology. A small subset of "strong PTC" features (>500nt, >=2 EJCs, GENCODE) does show a positive signal in DD_ALI, but this represents ~3% of GENCODE isoforms and emerges from exploratory analysis.

---

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
| `rmats_case_studies.pdf` | SRSF1/SRSF7 rMATS vs isoform-level comparison |
