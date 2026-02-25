# NMD Isoform Transitions — Results Summary

**Date:** 2026-02-21
**Pipeline:** Version 6.0
**Primary threshold:** non-NMD adj.P.Val > 0.50

---

## 1. Dataset Overview

### 1.1 Input Data

The pipeline analyzed isoform-level expression from PacBio long-read RNA sequencing across 6 cell types (AT2, DD, DD_ALI, DO, FB, MV) under two conditions: DMSO (vehicle control) and Smg1i (SMG1 inhibitor that blocks NMD). Isoforms were classified as NMD-sensitive (adj.P.Val < 0.05, logFC > 0 under Smg1i) or non-NMD (adj.P.Val > 0.50) based on differential expression.

### 1.2 Comparison Framework

Four comparison types were defined to distinguish NMD-triggering structural differences from baseline splicing variation:

| Comparison | Dominant | Comparator | Pairs/gene | Purpose |
|:---:|---|---|:---:|---|
| **C1** | Dominant non-NMD (DMSO) | Dominant NMD-sensitive (Smg1i) | 1 | True NMD transitions (expression-based dominance) |
| **C2** | Top non-NMD by CPM | Top NMD-sensitive by CPM | 1 | NMD transitions (abundance-based selection) |
| **C3** | Dominant non-NMD | All other non-NMD | Multiple | Baseline variation (full range) |
| **C4** | Dominant non-NMD | Next-best non-NMD by CPM | 1 | Baseline variation (matched 1:1 design) |

**Scope:** C3 was excluded from cross-comparison analysis (Scripts 14-15) because C4 is a strict subset of C3, creating pseudo-replication. DO cell type was excluded from per-cell-type analyses due to insufficient NMD pairs (C1: 1 pair, C2: 6 pairs at 0.50 threshold).

### 1.3 Splicing Event Detection

Pairs were deduplicated across all comparison x cell-type sets, yielding **39,647 unique pairs** covering **10,688 genes**. Hierarchical event detection identified **12 event types** across these pairs.

---

## 2. Pooled Splicing Profile Characterization (Scripts 09-13)

These analyses describe the overall architecture of splicing differences between isoform pairs, pooled across all comparisons.

### 2.1 Event Complexity (Script 09)

Isoform pairs differ by a median of **3 events** (interquartile range: 2-4). Comparator isoform exon count (a proxy for transcript complexity) shows a weak but significant association with the number of detected events (R^2 ~ 0.06). Profiles were binned into quartiles for downstream stratification:

| Bin | n_exons range | Profiles | % Combined pattern |
|-----|:---:|:---:|:---:|
| Q1 (low) | Fewest exons | 12,471 | 65.4% |
| Q2 (medium-low) | — | 7,353 | 81.8% |
| Q3 (medium-high) | — | 10,884 | 86.9% |
| Q4 (high) | Most exons | 8,939 | 90.1% |

Higher-complexity transcripts are more likely to exhibit combined (multi-category) splicing patterns (Chi-square = 3,410, Cramer's V = 0.169, p < 0.001).

### 2.2 Event Co-occurrence (Script 10)

Analysis of 28 pairwise event combinations across 39,647 profiles revealed systematic co-occurrence and mutual exclusion patterns:

**Co-occurring event pairs (OR > 1.3):**

| Event Pair | Odds Ratio | 95% CI | Interpretation |
|---|:---:|---|---|
| A5SS + A3SS | 1.45 | 1.37–1.54 | Boundary events co-occur |
| A5SS + SE | 1.45 | 1.36–1.54 | Boundary + inclusion co-occur |
| A3SS + SE | 1.41 | 1.34–1.48 | Boundary + inclusion co-occur |
| Alt_TSS + Alt_TES | 1.39 | 1.32–1.47 | Terminal events co-occur |
| Alt_TES + Partial_IR | 1.65 | 1.57–1.73 | Terminal + partial retention |

**Mutually exclusive event pairs (OR < 0.5):**

| Event Pair | Odds Ratio | 95% CI | Interpretation |
|---|:---:|---|---|
| Partial_IR + IR | 0.40 | 0.37–0.43 | Full vs. partial retention mutually exclusive |
| IR + Missing_Internal | 0.37 | 0.34–0.40 | Retention vs. missing exons exclusive |
| Alt_TES + SE | 0.45 | 0.43–0.47 | Terminal vs. inclusion exclusive |
| Alt_TES + IR | 0.46 | 0.43–0.49 | Terminal vs. retention exclusive |
| Partial_IR + Missing_Internal | 0.43 | 0.41–0.46 | Retention vs. missing exons exclusive |
| SE + Missing_Internal | 0.69 | 0.65–0.73 | Two inclusion types weakly exclusive |

### 2.3 Spatial Organization (Script 11)

**Proximity analysis:**
Events cluster significantly closer together than expected by chance. Observed mean inter-event distance = **13,540 bp** vs. permutation expectation = **15,004 bp** (z = -10.1, p < 0.001, 10,000 permutations). This 1,464 bp reduction suggests that splicing differences between isoforms tend to be spatially concentrated.

**Event topology:**
Among profiles with exactly 2 internal events, face-to-face (F2F) topology was significantly enriched:

| Topology | Observed % | Permutation p |
|---|:---:|:---:|
| Interleaved | 69.5% | < 0.001 |
| Back-to-back (B2B) | 16.5% | < 0.001 |
| Face-to-face (F2F) | 14.0% | < 0.001 |

**Positional bias:**
All event types showed a central/uniform distribution along the transcript (KS test vs. uniform, all p < 0.001 except IR_diff_5_3 with n = 83). No event type showed a strong 5' or 3' positional bias.

### 2.4 Functional Context (Script 12)

**Regional enrichment/depletion of events:**

| Event Type | 5'UTR | CDS | 3'UTR | ORF Start | ORF Stop |
|---|:---:|:---:|:---:|:---:|:---:|
| **Alt_TSS** | **Enriched** (2.2x) | Depleted (0.29x) | **Enriched** (2.0x) | Depleted (0.56x) | Depleted (0.65x) |
| **Alt_TES** | Depleted (0.71x) | Depleted (0.27x) | Depleted (0.62x) | Depleted (0.79x) | Depleted (0.89x) |
| **A5SS** | Depleted (0.62x) | Depleted (0.45x) | Depleted (0.59x) | Depleted (0.13x) | Depleted (0.14x) |
| **A3SS** | Depleted (0.27x) | Depleted (0.59x) | Depleted (0.28x) | Depleted (0.18x) | Depleted (0.17x) |
| **SE** | Depleted (0.78x) | Depleted (0.65x) | Depleted (0.75x) | Depleted (0.05x) | Depleted (0.08x) |
| **Missing_Internal** | **Enriched** (1.21x) | Depleted (0.43x) | **Enriched** (1.12x) | Depleted (0.11x) | Depleted (0.14x) |
| **IR** | NS | Depleted (0.78x) | Depleted (0.81x) | Depleted (0.37x) | Depleted (0.44x) |

*Note: Most events are intronic (41-77% depending on type). Enrichment ratios are relative to the expected proportion of each region among non-intronic exons. Alt_TSS is the only event type enriched in both 5'UTR and 3'UTR.*

**ORF boundary susceptibility:**
A5SS and A3SS events are depleted at ORF boundary exons (OR = 0.69, p < 10^-20), suggesting splice site variation avoids positions critical for reading frame maintenance.

**ORF impact of terminal events:**

| Event | % Affecting ORF Boundary | n Events |
|---|:---:|:---:|
| Alt_TSS affecting ORF start | 12.5% | 4,011 / 32,204 |
| Alt_TES affecting ORF stop | 15.3% | 4,173 / 27,332 |

### 2.5 Profile Pattern Classification (Script 13)

Splicing choice profiles were classified into 5 primary categories based on the event types present:

1. **Terminal:** Alt_TSS, Alt_TES
2. **Boundary:** A5SS, A3SS
3. **Inclusion:** SE, Missing_Internal
4. **Partial Retention:** Partial_IR_5, Partial_IR_3
5. **Full Retention:** IR, IR_diff_5/3/5_3

Profiles containing events from multiple categories were classified as **Combined**.

**Top 10 profile signatures:**

| Rank | Signature | n | % | Cumulative % |
|:---:|---|:---:|:---:|:---:|
| 1 | Terminal + Inclusion | 6,400 | 20.2% | 20.2% |
| 2 | Terminal + Partial Retention | 5,313 | 16.8% | 37.0% |
| 3 | Terminal + Boundary + Inclusion | 4,020 | 12.7% | 49.7% |
| 4 | Terminal + Boundary | 3,779 | 11.9% | 61.6% |
| 5 | Terminal + Boundary + Partial Retention | 2,839 | 9.0% | 70.6% |
| 6 | Terminal + Inclusion + Partial Retention | 2,606 | 8.2% | 78.8% |
| 7 | Terminal + Full Retention | 1,941 | 6.1% | 84.9% |
| 8 | Terminal + Boundary + Inclusion + Partial Retention | 875 | 2.8% | 87.7% |
| 9 | Terminal + Inclusion + Full Retention | 706 | 2.2% | 89.9% |
| 10 | Terminal + Boundary + Full Retention | 703 | 2.2% | 92.1% |

Terminal events are nearly ubiquitous — **all top 10 signatures include Terminal events**. The most common combined pattern (Terminal + Inclusion, 20.2%) reflects isoform pairs that differ at both their TSS/TES boundaries and through exon skipping or missing internal exons.

---

## 3. Per-Comparison Analysis (Script 14)

Script 14 filtered the 39,647 deduplicated profiles to each comparison x cell-type combination and ran Scripts 09-13 on each subset.

### 3.1 Profile Counts

| Run | C1 | C2 | C4 |
|---|:---:|:---:|:---:|
| all_samples | 29* | 676 | 1,156 |
| AT2 | 98 | 640 | 5,533 |
| DD | 250 | 1,425 | 4,611 |
| DD_ALI | 254 | 880 | 2,348 |
| FB | 51 | 303 | 5,694 |
| MV | 178 | 939 | 5,097 |

*Skipped (n < 50 minimum). **17 of 18 runs completed successfully.**

C1 has substantially fewer pairs than C2 or C4 because it requires a gene to have both a dominant non-NMD isoform AND a dominant NMD-sensitive isoform in the same cell type — a restrictive criterion. C4 has the most pairs because it only requires two non-NMD isoforms per gene.

---

## 4. Cross-Comparison Statistical Framework (Script 15)

This analysis asks: **Do NMD-triggering isoform transitions (C1, C2) differ structurally from baseline splicing variation (C4)?**

### 4.1 Phase 1: NMD vs. Baseline — Per-Cell-Type Tests

#### 4.1.1 Event Complexity (Unpaired)

| Comparison | Run | NMD n | Baseline n | NMD Median | Baseline Median | Cliff's delta | p |
|:---:|---|:---:|:---:|:---:|:---:|:---:|:---:|
| C1 | AT2 | 98 | 5,533 | 3 | 3 | 0.016 | 0.782 |
| C1 | DD | 250 | 4,611 | 3 | 3 | -0.059 | 0.106 |
| C1 | DD_ALI | 254 | 2,348 | 3 | 3 | 0.035 | 0.347 |
| C1 | FB | 51 | 5,694 | 3 | 3 | 0.100 | 0.203 |
| C1 | MV | 178 | 5,097 | 3 | 3 | 0.005 | 0.908 |
| C2 | all_samples | 676 | 1,156 | 3 | 3 | 0.059 | **0.029** |
| C2 | AT2 | 640 | 5,533 | 3 | 3 | 0.007 | 0.756 |
| C2 | DD | 1,425 | 4,611 | 3 | 3 | 0.019 | 0.259 |
| C2 | DD_ALI | 880 | 2,348 | 3 | 3 | 0.078 | **< 0.001** |
| C2 | FB | 303 | 5,694 | 3 | 3 | 0.102 | **0.002** |
| C2 | MV | 939 | 5,097 | 3 | 3 | 0.067 | **< 0.001** |

Both NMD and baseline pairs show a median of 3 events. C1 shows no significant differences in any cell type. C2 shows small but significant positive Cliff's deltas in 4 of 6 runs, indicating slightly more events in NMD pairs.

#### 4.1.2 Paired Analysis (Same Genes in Both Comparisons)

For genes appearing in both NMD and baseline comparisons, paired tests were performed:

| Comparison | Run | n Paired | Median Diff | Wilcoxon p | % NMD More | % Profile Type Change |
|:---:|---|:---:|:---:|:---:|:---:|:---:|
| C1 | DD | 107 | 0 | 0.235 | 54.4% | 36.4% |
| C1 | DD_ALI | 94 | 0 | 0.127 | 59.7% | 26.6% |
| C1 | MV | 85 | 0 | 0.155 | 42.9% | 34.1% |
| C2 | all_samples | 444 | 0 | 0.252 | 53.2% | 24.1% |
| C2 | AT2 | 497 | 0 | 0.687 | 51.2% | 29.2% |
| C2 | DD | 1,043 | 0 | 0.613 | 51.1% | 26.7% |
| C2 | DD_ALI | 520 | 0 | **0.045** | 56.5% | 29.2% |
| C2 | FB | 251 | 0 | **0.043** | 56.8% | 29.9% |
| C2 | MV | 684 | 0 | 0.089 | 54.0% | 27.0% |

Paired analyses confirm the same pattern: median event count difference is 0 in all cases. C2 DD_ALI and FB show marginal significance (p ~ 0.04), with ~57% of NMD pairs having more events than their baseline counterparts for the same gene. Approximately 25-36% of gene pairs change their profile type classification between NMD and baseline comparisons.

### 4.2 Phase 2: Meta-Analysis Across Cell Types

Random-effects meta-analysis (REML) pooled effect sizes across 5 cell types (AT2, DD, DD_ALI, FB, MV):

#### 4.2.1 Event Complexity and Profile Types

| Metric | C1 Pooled | C1 95% CI | C1 p | C2 Pooled | C2 95% CI | C2 p |
|---|:---:|---|:---:|:---:|---|:---:|
| Cliff's delta (events) | 0.002 | -0.048, 0.051 | 0.942 | 0.050 | 0.017, 0.082 | **0.003** |
| Cramer's V (profile types) | 0.026 | 0.008, 0.044 | **0.004** | 0.031 | 0.014, 0.048 | **< 0.001** |

C1 shows no event complexity difference (delta ~ 0). C2 shows a small but significant excess of events in NMD pairs (pooled delta = 0.05). Both show statistically significant but very small profile type shifts (Cramer's V ~ 0.03, negligible effect size).

#### 4.2.2 Event Type Prevalence (Odds Ratios, NMD vs. Baseline)

| Event Type | C1 OR | C1 95% CI | C1 p | C1 I^2 | C2 OR | C2 95% CI | C2 p | C2 I^2 |
|---|:---:|---|:---:|:---:|:---:|---|:---:|:---:|
| **Alt_TSS** | 1.03 | 0.86–1.22 | 0.775 | 0% | 1.08 | 0.99–1.18 | 0.074 | 0% |
| **Alt_TES** | **1.24** | **1.06–1.44** | **0.008** | 0% | 0.98 | 0.87–1.10 | 0.679 | 61% |
| **A5SS** | 1.03 | 0.85–1.24 | 0.786 | 0% | 1.09 | 0.98–1.21 | 0.113 | 22% |
| **A3SS** | 0.97 | 0.81–1.16 | 0.754 | 9% | 1.07 | 0.98–1.16 | 0.124 | 12% |
| **Partial_IR** | 0.92 | 0.79–1.07 | 0.282 | 0% | 1.05 | 0.97–1.13 | 0.203 | 1% |
| **IR** | 1.06 | 0.86–1.29 | 0.590 | 0% | 0.98 | 0.89–1.08 | 0.703 | 0% |
| **SE** | 0.96 | 0.72–1.27 | 0.748 | 62% | **1.13** | **1.03–1.24** | **0.013** | 29% |
| **Missing_Internal** | 1.13 | 0.97–1.33 | 0.124 | 0% | 1.07 | 0.98–1.17 | 0.124 | 18% |

**Key findings:**
- **C1 — Alt_TES enrichment:** The only significant event-level signal for C1. NMD-sensitive dominant isoforms are 24% more likely to differ from non-NMD dominants at the transcription end site (OR = 1.24, p = 0.008). Low heterogeneity (I^2 = 0%) indicates this is consistent across cell types.
- **C2 — SE enrichment:** NMD transitions identified by top CPM are 13% more likely to involve skipped exons (OR = 1.13, p = 0.013). Moderate heterogeneity (I^2 = 29%).
- **All other event types** show odds ratios close to 1.0 with no significant differences.

### 4.3 Phase 3: Sensitivity Analyses

#### 4.3.1 C1 vs. C2 Concordance

Correlation of event-level Cliff's deltas between C1 and C2 across 5 cell types:

| Cell Type | Pearson r | p | Direction Concordance |
|---|:---:|:---:|:---:|
| AT2 | 0.22 | 0.606 | Yes |
| DD | 0.51 | 0.194 | No |
| DD_ALI | **0.78** | **0.024** | Yes |
| FB | 0.06 | 0.892 | Yes |
| MV | -0.45 | 0.265 | Yes |

C1 and C2 agree on the direction of overall event complexity effects in 4 of 5 cell types. However, event-level correlations are variable — only DD_ALI shows significant concordance (r = 0.78). This suggests that while C1 and C2 capture related but not identical aspects of NMD biology, their specific event-type enrichments may differ.

#### 4.3.2 Complexity Confound

A potential confound: do NMD comparator isoforms differ in baseline transcript complexity from baseline comparators?

| Comparison | Pattern | Significant Runs |
|---|---|---|
| **C1** | NMD comparators have slightly fewer exons (median 4 vs. 5) | DD only (p = 0.005) |
| **C2** | NMD comparators have **more** exons (median 5-6 vs. 5) | AT2, DD_ALI, FB, MV (p < 0.01) |

The C2 complexity confound (NMD comparators are more complex) may partially explain C2's event count excess. C1 shows the opposite pattern (NMD comparators slightly simpler) but this is mostly non-significant.

#### 4.3.3 Event Direction (GAIN vs. LOSS)

The proportion of LOSS events (comparator lost sequence relative to dominant) is ~54% in both NMD and baseline comparisons across all cell types. No systematic directional asymmetry distinguishes NMD transitions from baseline variation.

| Comparison | LOSS Proportion (NMD) | LOSS Proportion (Baseline) | Largest |diff|| Significant Runs |
|---|:---:|:---:|:---:|---|
| C1 | 0.47–0.56 | 0.54–0.55 | 0.079 (FB) | FB only (p = 0.034) |
| C2 | 0.54–0.57 | 0.54–0.54 | 0.032 (DD_ALI) | DD_ALI only (p = 0.003) |

---

## 5. PTC and Frame Disruption Analysis

A parallel analysis tested whether premature termination codons (PTCs) and reading-frame disruption predict NMD responsiveness. This connects the isoform comparison framework (Sections 2-4) to the broader question of what makes an isoform an NMD substrate.

### 5.1 rMATS frame-disruption enrichment

Short-read rMATS differential splicing analysis found that frame-disrupting splice events are ~2x enriched among Smg1i-responsive events (meta OR = 2.13, p = 3e-9) across SE, A5SS, A3SS, and MXE event types. DD, MV, and AT showed significant enrichment.

### 5.2 Isoform-level PTC null

Isoform-level PTC status (computed via the 50nt EJC rule) does NOT predict NMD responsiveness (meta OR = 0.96, p = 0.03, with DD_ALI driving the signal in the wrong direction). This null holds for GENCODE-only analysis (OR = 0.88) and GENCODE NMD biotype (OR = 0.96). No dose-response with PTC distance.

### 5.3 Cross-platform concordance resolves the discrepancy

Systematic junction matching between rMATS events and PacBio isoforms (Script 11, results/ptc/) explains why the rMATS signal does not propagate to isoform-level PTC enrichment:

| Bucket | Description | % of events |
|--------|-------------|---:|
| 1 | No PacBio gene coverage | 2.8% |
| 2 | No junction match (forms missing from catalog) | 69.1% |
| 3 | Both forms present, no PTC difference | 17.0% |
| 4 | PTC differs, but not NMD-sensitive | 8.5% |
| 5 | Fully concordant | 2.7% |

**The primary explanation is isoform catalog incompleteness** (71.8% of events). PacBio does not have isoforms representing most of the splice variants rMATS detects. Among events with both forms present, frame disruption does not reliably predict PTC status (42.5% PTC-difference rate, identical to frame-preserving controls).

For full details, see `results/ptc/RESULTS.md` and `results/ptc/METHODS.md`.

---

## 6. Summary of Key Findings

### 6.1 Splicing Architecture is Remarkably Conserved

The most striking finding is the **similarity** between NMD-triggering transitions and baseline splicing variation:
- Median event count is **3** in both NMD and baseline pairs
- Profile type distributions are nearly identical (Cramer's V ~ 0.03)
- Event direction balance (~54% LOSS) is the same
- No dramatic enrichment or depletion of any event type

This suggests that the structural mechanisms distinguishing NMD-sensitive isoforms from non-NMD isoforms are drawn from the same repertoire of splicing variation that distinguishes any two isoforms of the same gene.

### 6.2 Subtle but Consistent Signals

Two event-level signals emerge from meta-analysis:

1. **Alt_TES enrichment in C1 (OR = 1.24, p = 0.008):** Isoform pairs defined by expression-based NMD dominance are more likely to differ at the transcription end site. This is consistent across all 5 cell types (I^2 = 0%). Alternative polyadenylation or 3' end variation may be mechanistically relevant to NMD triggering in these pairs.

2. **SE enrichment in C2 (OR = 1.13, p = 0.013):** Isoform pairs defined by top CPM NMD selection are more likely to involve skipped exons. However, this signal should be interpreted cautiously given the complexity confound — C2 NMD comparators tend to be more complex transcripts, which may increase the baseline probability of observing SE events.

### 6.3 C1 and C2 Capture Different Facets

C1 (dominant NMD-sensitive) and C2 (top NMD-sensitive by CPM) do not produce identical results:
- C1 detects Alt_TES enrichment; C2 does not
- C2 detects overall event complexity increase and SE enrichment; C1 does not
- Their event-level concordance is variable across cell types

This likely reflects their different selection criteria: C1 requires the NMD isoform to be the single most expressed isoform under Smg1i treatment, while C2 selects the highest-CPM NMD isoform regardless of relative expression rank.

### 6.4 Cell-Type Consistency

Meta-analysis heterogeneity (I^2) is generally low (0-30%) for most metrics, indicating that NMD-vs-baseline differences are consistent across the 5 analyzed cell types. The main exception is Alt_TES in C2 (I^2 = 61%) and SE in C1 (I^2 = 62%), where cell-type-specific effects are present.

### 6.5 Cross-Platform Concordance Reveals Isoform Catalog Limitations

The rMATS frame-disruption enrichment (OR = 2.13) does not propagate to isoform-level PTC analysis because 69% of significant frame-disrupting events have no matching junction pair in the PacBio isoform catalog. This identifies a fundamental limitation of long-read isoform-level analysis: many biologically relevant splice variants detected by short-read event-level methods were not assembled as distinct isoforms by PacBio IsoSeq or annotated in GENCODE. The 2.7% fully-concordant rate quantifies the narrow overlap between these two analytical frameworks.

---

## 7. Generated Figures

### 7.1 Pooled Downstream Analysis Figures

| Figure | Location | Description |
|---|---|---|
| Complexity vs. Events | `deduplicated/results/figures/complexity_vs_events.pdf` | Scatter/density plot of transcript complexity vs. event count |
| Co-occurrence Heatmap | `deduplicated/results/figures/cooccurrence_heatmap.pdf` | Heatmap of pairwise event co-occurrence odds ratios |
| Pattern Classification | `deduplicated/results/figures/pattern_classification.pdf` | Distribution of splicing profile pattern categories |
| Regional Enrichment | `deduplicated/results/figures/regional_enrichment.pdf` | Event enrichment/depletion across genomic regions |
| Spatial Patterns | `deduplicated/results/figures/spatial_patterns.pdf` | Event proximity and topology analyses |

### 7.2 Cross-Comparison Analysis Figures

| Figure | Location | Description |
|---|---|---|
| Profile Types by Comparison | `cross_comparison/figures/profile_types_by_comparison.pdf` | Side-by-side profile type distributions (C1, C2, C4) |
| Event Prevalence Comparison | `cross_comparison/figures/event_prevalence_comparison.pdf` | Event type proportions across comparisons |
| Paired Event Differences | `cross_comparison/figures/paired_event_differences.pdf` | Distribution of within-gene event count differences |
| Forest Plot — C1 Event ORs | `cross_comparison/figures/forest_event_ORs_C1.pdf` | Meta-analysis forest plots for C1 event type odds ratios |
| Forest Plot — C2 Event ORs | `cross_comparison/figures/forest_event_ORs_C2.pdf` | Meta-analysis forest plots for C2 event type odds ratios |
| Forest Plot — C1 Events | `cross_comparison/figures/forest_events_C1.pdf` | Cell-type-specific event ORs for C1 |
| Forest Plot — C2 Events | `cross_comparison/figures/forest_events_C2.pdf` | Cell-type-specific event ORs for C2 |

All figure paths are relative to `comparisons/nonNMD_0.50/`.

---

## 8. Output File Inventory

All output files are located under `comparisons/nonNMD_0.50/`.

### 8.1 Pooled Results (`deduplicated/results/`)

| File | Description |
|---|---|
| `splicing_choice_profiles_with_bins.rds` | 39,647 profiles with complexity quartile bins |
| `complexity_relationship_results.rds` | Linear model results |
| `complexity_bins_definition.rds` | Quartile boundary definitions |
| `cooccurrence_crude_results.tsv` | 28 pairwise event co-occurrence tests |
| `cooccurrence_stratified_results.tsv` | Stratified co-occurrence with ORs and CIs |
| `positional_bias_results.tsv` | KS test results for 12 event types |
| `proximity_analysis_results.tsv` | Permutation test for event clustering |
| `topology_enrichment_results.tsv` | F2F/B2B/Interleaved enrichment |
| `regional_distribution_results.tsv` | Event enrichment across genomic regions |
| `orf_boundary_susceptibility.tsv` | A5SS/A3SS depletion at ORF boundaries |
| `orf_impact_summary.tsv` | Alt_TSS/Alt_TES impact on ORF start/stop |
| `pattern_frequencies_by_complexity.tsv` | Pattern types stratified by complexity quartile |
| `pattern_complexity_association.tsv` | Chi-square test of pattern x complexity |
| `top_patterns.tsv` | Top 10 profile pattern signatures |

### 8.2 Per-Comparison Results

`per_comparison_run_summary.tsv` — Summary of 18 comparison x run combinations (17 completed, 1 skipped).

Per-run results stored in `{C1,C2,C4}/{run}/results/` with the same file structure as pooled results.

### 8.3 Cross-Comparison Results (`cross_comparison/`)

| File | Description |
|---|---|
| `phase1_unpaired_results.tsv` | 11 unpaired NMD vs. baseline tests |
| `phase1_paired_results.tsv` | 9 paired within-gene tests |
| `phase1_event_prevalence.tsv` | 160 event-type prevalence tests |
| `phase1_regional_bootstrap.tsv` | ~700 regional enrichment bootstrap CIs |
| `phase2_meta_analysis.tsv` | 20 random-effects meta-analysis models |
| `phase3_c1c2_concordance.tsv` | C1-C2 agreement across 5 cell types |
| `phase3_complexity_confound.tsv` | Complexity confound analysis (33 tests) |
| `phase3_direction.tsv` | GAIN/LOSS balance comparison (11 tests) |

### 8.4 PTC and Concordance Results (`results/ptc/results/`)

| File | Description |
|---|---|
| `rmats_junction_concordance.tsv` | Per-event junction concordance classification |
| `rmats_ptc_verification.tsv` | PTC status comparison for events with both PacBio forms |
| `rmats_psi_concordance.tsv` | rMATS vs PacBio dPSI comparison |
| `rmats_discrepancy_waterfall.tsv` | Sequential bucket assignment explaining signal loss |
| `rmats_concordance_summary.tsv` | Summary statistics from concordance analysis |

See `results/ptc/RESULTS.md` for complete output file listing including all PTC and rMATS analysis outputs.

---

## 9. Validation Status

- **Event detection test suite:** 126/126 PASS (100%) — 44 synthetic + 80 real-data + 2 deduplication tests
- **GENCODE validation:** 4,274/4,274 PASS (100%) — random dominant assignment, seed = 42
- **Reconstruction verification:** All detected events can be used to perfectly reconstruct the dominant isoform from the comparator + events, confirming detection accuracy
