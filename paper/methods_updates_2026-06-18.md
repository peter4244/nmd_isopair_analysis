# Methods updates for §4 + §5 analyses

**Date:** 2026-06-18
**Companion to:** `paper/manuscript_reconciliation_2026-06-18.md`
**Purpose:** Drafted Methods subsections ready to paste into the manuscript Google Doc. Seven inserts; three REPLACE existing subsections, four are NEW.

Reading order:

| # | Heading | Action | Notes |
|---|---|---|---|
| 1 | Isoform pairs analysis (Isopair package) | REPLACE | adds 12-event taxonomy + three nested scopes + `enumerateOrfs` priority |
| 2 | PTC determination | REPLACE | adds Kozak scoring, TD2 bias remediation, ref-AUG-projection algorithm |
| 3 | Mechanism class taxonomy | NEW | formal definitions of frameshift / in-frame stop / 3′UTR splice |
| 4 | Deep learning model — network architecture | APPEND | CNN layers, kernel sizes, attention aggregator. **Verify against `model.py` when cluster is up** — placeholders marked `[verify]` |
| 5 | Deep learning model — interpretability methods | REPLACE | DeepSHAP / KernelSHAP / attention; n=1,166 subgroup-stratified analyses |
| 6 | uORF rule (Fig 5 Panel G) | NEW | three structural criteria used in the attention analysis |
| 7 | Stop-codon analysis correction | NEW | discloses the 2026-04-30 patch and notes that post-patch UGA = ~55.9% NMD vs ~48.1% Control |

---

## 1. REPLACE — Isoform pairs analysis (Isopair package)

> **Isoform pairs analysis (Isopair package).** Isoform pairs were constructed using the Isopair R package (https://github.com/peter4244/Isopair) to link splicing events to NMD susceptibility within each protein-coding gene. For each gene with ≥3 expressed isoforms (≥5% of overall gene expression in DMSO or SMG1i, ≥5 reads in at least one sample) and at least one NMD-susceptible isoform, we identified a high-expressing non-NMD *reference* isoform and paired it with both an NMD-susceptible *comparator* and a non-NMD *control comparator* from the same gene. This yielded 3,009 NMD pairs and 3,009 Control pairs across 3,009 genes (this universe is referred to as pop_BC). Across the 3,009 pair-genes, the median number of expressed non-NMD isoforms per gene was 7 (IQR = 5), and the median reference isoform captured 31% of the parent gene's non-NMD DMSO 4-CT (AT, DD, FB, MV) expression (38% of references ≥ 50%) — see Supplemental Figure *PairSetDescriptives* Panels A / B. The median spliced transcript lengths were 3,049 nt for the NMD comparator, 2,893 nt for the reference, and 2,762 nt for the Control comparator (Kruskal–Wallis p = 1.5×10⁻⁷; pairwise Mann–Whitney U p ≤ 0.034 for all three contrasts; Supplemental Figure *PairSetDescriptives* Panel C).
>
> Splice events were enumerated into twelve mutually exclusive categories: skipped exon (SE), 5′ alternative splice site (A5SS), 3′ alternative splice site (A3SS), intron retention (IR), intron-retention 5′-difference (IR diff 5′), intron-retention 3′-difference (IR diff 3′), partial intron retention 5′ (Partial IR 5′), partial intron retention 3′ (Partial IR 3′), alternative transcription start site (Alt TSS), alternative transcription end site (Alt TES), missing internal exon block, and a residual "other" category. Categories were assigned by comparing splice-junction coordinates between the reference and comparator transcripts; multiple events per pair were allowed, with prevalence (Fig 3 Panel C) computed as the fraction of pairs containing at least one event of that category.
>
> Three nested scopes are used downstream. **Scope §2a (all-3-ENST coding, n = 190):** all three isoforms of the pair are GENCODE-annotated with a CDS; supports Fig 3 Panels D / E / F and Fig 4 Panels A / B. **Scope §3a (ENST-reference + ref-AUG-traceable):** the reference isoform is GENCODE-annotated, and the comparator's transcript projects to the same reference AUG (verified by `Isopair::enumerateOrfs()` Kozak-aware scan of all candidate ORFs and a downstream agreement check against the 5′UTR-features pipeline's `utr5_features_refaug.rds`). **Scope §3b (re-intersected gene-matched, n = 1,166 NMD + 1,166 Control):** §3a additionally re-intersected to 1:1 gene-matched pairs after categorical filtering — this is the canonical n = 1,166 used in Fig 4 Panels C / D, Fig 5 Panel G, the TD2BiasEvidence supplement (broad-scope panels), the PTCSubclassBranchSHAP supplement, and the PTCSubclassPerformance supplement.
>
> Ref-AUG traceability is determined by projecting the GENCODE reference AUG coordinate onto the comparator transcript via splice-aware coordinate transformation; if the projected position is within the comparator's coordinates and the resulting ORF can be enumerated, the comparator is traceable. The `enumerateOrfs()` ORF selector enumerates all candidate ORFs and ranks them by (1) reference-CDS match, (2) TD2-called CDS, (3) Kozak score with start-codon plausibility — this priority is the same selector used to populate the deep-learning model's K = 5 candidate ORFs (Methods, "Deep learning model").


---

## 2. REPLACE — PTC determination

> **PTC determination.** A stop codon was classified as a premature termination codon (PTC) under the canonical 50-nucleotide rule: the stop is a PTC if it occurs >50 nt upstream of the last exon–exon junction (terminal EJC). Distance from the comparator stop codon to the terminal EJC is computed in transcript coordinates from the comparator splice graph (Fig 3 Panel D).
>
> For GENCODE-annotated pairs (n = 190 all-3-ENST coding scope), PTC status is assigned directly from the comparator's own GENCODE CDS annotation. For novel comparator isoforms lacking GENCODE annotation, two complementary CDS calls are considered: (1) the **TD2 (TransDecoder2)** CDS call as inherited from the SQANTI3 isoform-classification pipeline; (2) a **reference-AUG-projected** ORF derived by projecting the GENCODE reference isoform's AUG coordinate onto the comparator transcript via splice-aware coordinate transformation.
>
> **TD2 bias remediation.** TD2 ranks candidate ORFs by length and avoids ORFs whose first stop codon would create a PTC, systematically eliminating the signal that NMD analyses depend on (TD2BiasEvidence supplement). In the occult-PTC subset (effectively_ptc ∩ original_ptc==FALSE; n = 492), 487 of 492 (99.0%) TD2-called CDSes are downstream of the reference AUG, and reference AUGs have stronger Kozak context than TD2 AUGs in 78.0% of pairs (paired Wilcoxon p = 2.6×10⁻³⁸). For §4 / §5 analyses that require an unbiased main-CDS call, we therefore use the reference-AUG-projected ORF rather than the TD2 call. The exception is uORF-mediated NMD (e.g., ATF4): because uORF detection operates on the 5′UTR upstream of the main CDS, it is unaffected by TD2's main-CDS bias and the TD2 CDS is acceptable.
>
> **Kozak scoring.** Position weight matrix Kozak scores were computed using the canonical PWM (Kozak 1987 / Hernández 2019 conventions): the −3 G/A, −2 C, −1 C, +4 G positions are scored against the comparator's start-context nucleotides. A score ≥ 1 in the `Isopair::enumerateOrfs()` scale was used to flag ORFs with credible start contexts (one of three uORF-eligibility filters in §5).
>
> **Mechanism class derivation.** Where a PTC has been identified, the splicing event responsible is attributed by the `mechanism_class` helper (`figures/lib/mechanism_class.R`), which classifies each PTC-attributed event as frameshift, in-frame stop, or 3′ UTR splice — see "Mechanism class taxonomy" below for the exclusive definitions.

---

## 3. NEW — Mechanism class taxonomy

*(Insert immediately after "PTC determination". Current Methods text has no formal definition of the three-way classification used in Fig 3 Panel F — only the figure legend describes them.)*

> **Mechanism class taxonomy.** For each PTC-attributed splicing event (n = 69 events across 72 PTC+ NMD pairs in the n = 190 scope; Fig 3 Panel F), the event is assigned to one of three exclusive mechanism classes based on its effect on the reading frame:
>
> - **Frameshift** (38 / 69 events, 55.1%): the event introduces a length change that is not a multiple of 3 nucleotides, shifting the reading frame and exposing a premature stop downstream;
> - **In-frame stop** (23 / 69, 33.3%): the event preserves the reading frame but introduces a new in-frame stop codon — typically by including sequence content (e.g., a retained intron, an alternative exon) whose first frame includes a stop;
> - **3′ UTR splice** (8 / 69, 11.6%): the event removes sequence between the natural stop codon and the original terminal EJC, repositioning the terminal EJC to be >50 nt downstream of the (unchanged) stop codon — a PTC by definition under the 50-nt rule, but caused by a 3′UTR-region rearrangement rather than a CDS event.
>
> The implementation lives in `figures/lib/mechanism_class.R` as a derived helper (no cached column; recomputed from splice-graph coordinates at each invocation).

---

## 4. APPEND — Deep learning model: network architecture

*(Append to the existing "Deep Learning Model – Architecture and Training" subsection, after the description of windows / channels / structural features and BEFORE the training description. The current Methods text omits CNN-specific architecture details.)*

> **Network architecture.** The model is composed of three input branches whose outputs are concatenated into a 64-dimensional per-ORF embedding and then aggregated across the five candidate ORFs via a learned softmax attention layer to produce a single transcript-level embedding for the classification head.
>
> The **AUG branch** and **stop branch** each process a 9 × 500 sequence-and-feature window through a 1D convolutional stack: two convolutional blocks with kernel sizes 15 (first) and 7 (second), 32 output channels in each block, padding chosen to preserve length (kernel-size/2), each block followed by batch normalization and a ReLU activation. A max-pooling layer of kernel size 4 between the two convolutions reduces sequence length. After the second convolution an adaptive max-pooling layer collapses the variable-length output to a fixed-length vector of 32 dimensions, which is then linearly projected to a 32-dimensional per-branch embedding. The two branches do not share weights with each other but DO share weights across the five ORF slots within each branch (the same encoder is applied independently to each of the five candidate ORFs).
>
> The **structural branch** is a per-ORF linear projection of five hand-engineered structural features (`frac_start`, `frac_stop`, `is_ref_cds`, `is_sqanti_cds`, `n_downstream_ejc`) to a 32-dimensional embedding via a single fully-connected layer followed by ReLU activation.
>
> The three per-branch embeddings (32 + 32 + 32 = 96 dimensions) are concatenated per ORF and fused into a 64-dimensional per-ORF embedding through a single fully-connected layer with ReLU activation and dropout (rate 0.2). The five per-ORF embeddings are aggregated by a learned attention layer: a single linear projection (64 → 1) produces per-ORF attention scores that are masked-softmax normalized across the five ORFs (masking pads ORFs absent from a given transcript and gracefully handles transcripts with zero ORFs); the per-ORF embeddings are then weighted by the attention probabilities and summed to yield a single 64-dimensional transcript-level embedding. The transcript embedding feeds a two-layer multilayer perceptron classification head (64 → 32, ReLU, dropout 0.3, → 1) whose final scalar is the NMD logit.

---

## 5. REPLACE — Deep learning model: interpretability methods

> **Interpretability.** Three complementary attribution methods were used, each operating at a different layer of the model.
>
> **Joint DeepSHAP at the input layer.** DeepSHAP attributions were computed for the model's NMD-logit output with respect to the 9-channel input of the AUG and stop windows and the 5 structural features of the priority (rank-0) ORF, with ORFs at ranks 1–4 held fixed at their input values. Five independent replicate runs were used; each replicate drew 500 random background transcripts from the training set as the SHAP reference distribution. Per-position-per-channel attributions were averaged across replicates for downstream summarization (Fig 5 Panel D for structural-feature ranking; Figs 5 E / F for per-nucleotide signed attribution).
>
> **KernelSHAP at the fusion layer.** Exact Shapley decomposition across the three sub-encoders (AUG, stop, structural) was computed by running the model 2³ = 8 times per isoform with each combination of present/absent branches, using the trained model's expected branch output as the absent-branch value. This produces per-isoform branch-level attributions (`shap_atg`, `shap_stop`, `shap_structural`) that sum to the model's NMD-call signal. KernelSHAP was computed for the full cohort (`kernel_shap_branch_atg500_stop500_all.tsv`, n = 39,938 isoforms; not test-set restricted) and is the data source for Fig 5 Panel C and the PTCSubclassBranchSHAP supplement.
>
> **Attention weights at the aggregator.** Per-isoform, per-ORF attention probabilities were extracted from the softmax-normalized attention layer (Fig 5 Panel G; AttentionDistribution supplement).
>
> **Subgroup-stratified analyses (n = 1,166).** For the manuscript Subset 2 (n = 1,166 ref-AUG-traceable gene-matched pairs; Fig 4 C/D scope), three subgroups are defined by the upstream Isopair `category` field carried through to the figure-side TSV: NMD+/PTC+ retained (n = 1,016 with branch SHAP available), NMD+/PTC− retained (n = 95), and Control (n = 1,107). Branch SHAP descriptives at this scope are reported in the PTCSubclassBranchSHAP supplement at full-cohort KernelSHAP scope (not test-restricted). Predictive performance by subgroup (AUC, AUPRC, per-subgroup predicted probability) is reported in the PTCSubclassPerformance supplement at the test-set intersection (NMD+/PTC+ n = 255, NMD+/PTC− n = 30, Control n = 276); AUC and AUPRC are the only analyses we report at test-set scope — all other §5 interpretability analyses use the full cohort.

---

## 6. NEW — uORF rule (Fig 5 Panel G)

*(Insert at the end of the deep-learning Methods block. Currently the uORF rule is defined only in the figure legend.)*

> **uORF rule (Fig 5 Panel G).** A candidate ORF within a transcript was classified as a credible uORF for the attention analysis if it satisfied three structural criteria — purely upstream-side and independent of any main-CDS call (and therefore unaffected by TD2's PTC-avoidance bias):
>
> 1. **Length < 200 nt.** Below the 270-nt CDS floor of GENCODE-annotated protein-coding transcripts, so the rule cannot misclassify a real annotated CDS as a uORF.
> 2. **Position in 5′ half.** Both start and stop codons fall within the 5′ half of the transcript (`frac_start < 0.5` AND `frac_stop < 0.5`).
> 3. **Credible Kozak start.** Kozak score ≥ 1 in the `Isopair::enumerateOrfs()` scale.
>
> A transcript is reported as "placing >5% attention on a candidate uORF" if any of its five candidate ORF slots whose ORF satisfies the rule receives >5% of the model's softmax-normalized attention weight.

---

## 7. NEW — Stop-codon analysis correction

*(Insert as a short methods-side disclosure, OR as a footnote to the §5 paragraph that cites the UGA frequency. Pete's preference.)*

> **Stop-codon analysis correction.** An earlier off-by-one error in the ORFik stop-codon position column produced inflated UGA-in-NMD frequencies in renders generated before April 30, 2026. The column was patched on 2026-04-30 (`scripts/patch_stop_codon.py`) by regenerating the stop-codon column from the upstream HDF5 one-hot encoding. Post-patch, within-class stop-codon frequencies on the chr-1/3/5/7 paralog-free test set are UGA ~55.9% NMD vs ~48.1% Control [exact UAA / UAG breakdown to be populated from post-patch Rmd render]. The StopCodonUsage supplement carries the corrected figure; any numbers from pre-April-2026 renders (e.g., the prior 52.0% / 48.3% figures embedded in the pre-patch `orf_model_report_v5.html`) should be considered superseded.

---

## How to apply

1. **Open the Google Doc** at https://docs.google.com/document/d/1Tz6coXnDwpGZaV1Jl11fmN_LVX1R8YQPf2y_7wzBLKs/edit
2. For each subsection above with action **REPLACE**, locate the existing heading and replace the body with the drafted text.
3. For each **NEW** subsection, insert at the location indicated in the parenthetical guidance above.
4. For **APPEND**, paste at the indicated insertion point within the existing subsection.
5. **Before pasting subsection 4 (network architecture)**: read `NMD_orf_model_v5_4ct/model.py` when cluster is up, replace `[verify]` placeholders with actual values, or rewrite if the architecture diverges from this skeleton.
6. **After all pastes**: run a final pass to ensure cross-references (Fig 3, Fig 4, Fig 5 Panel X, supplement names) match the manuscript's current figure numbering.

Outstanding placeholders elsewhere in Methods (not §4 / §5 but caught during the audit) — these need filling before submission:

- Short-read RNA-seq: `[insert RIN value]`, `[insert minimum read depth]`, `[indicate one-pass versus two-pass mode]`
- Long-read RNA-seq: target FLNC depth `[fill in]`

These are out of scope for this Methods update but should be resolved before submission.
