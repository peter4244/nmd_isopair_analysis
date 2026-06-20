# NMD predictor comparison — our model vs Lindeboom 2019 (NMDetective) vs Saadat 2025 (NMDEP)

**Status:** Pre-registration, signed off 2026-06-20. Decisions locked in §2.

**Folder convention.** Following the existing analysis convention in `code/ptcneg_go_handoff/` — the analysis lives in `code/`, the rendered paper figure lives separately in `figures/SupplementalFigures/`, and the figure script consumes from the analysis outputs:

```
code/nmd_predictor_comparison/                       ← analysis pipeline (this folder)
├── METHODS.md                                       ← pre-registration; canonical
├── README.md                                        ← how to reproduce, dependencies
├── 01_extract_our_isoforms.R                        ← gather n=1,166 scope + features + gold standard
├── 02_score_nmdetective_b.R                         ← apply 4-rule decision tree to our isoforms
├── 03_score_nmdep_rule_baseline.R                   ← apply NMDEP-re-thresholded rules
├── 04_compute_metrics.R                             ← Spearman / Pearson / R² / MAE
├── nmd_predictor_comparison.Rmd                     ← reproducible report
├── nmd_predictor_comparison_2026.6.20.html          ← rendered report
├── isoforms_2026.6.20.tsv                           ← committed: per-isoform input table
├── per_isoform_scores_2026.6.20.tsv                 ← committed: each model's prediction
├── metrics_summary_2026.6.20.tsv                    ← committed: Spearman/Pearson/R²/MAE
└── external/                                        ← gitignored; only if NMDetective-A is run

figures/SupplementalFigures/ModelComparison/         ← rendered SF only after analysis is run
├── figure_s_model_comparison.py
├── figure_s_model_comparison.{pdf,png}
└── figure_s_model_comparison_legend.md
```

---

## 1. Question

How do existing variant-level NMD prediction models — NMDetective (Lindeboom et al., 2019) and NMDEP (Saadat & Fellay, 2025) — perform on splicing-derived NMD substrates measured in primary cells under controlled SMG1 inhibition, compared to our deep-learning isoform-level model?

This is fundamentally a **cross-dataset generalization test**. The prior models were trained on per-PTC NMD efficacy estimated from somatic-mutation allele-specific expression in TCGA tumor data. We test them on a different gold standard: per-isoform mashr-shrunk log fold-change under SMG1i in lung primary cells.

---

## 2. Locked-in design decisions

The following were agreed in conversation on 2026-06-20.

### 2.1 Cohort

Use the full **n = 1,166 ref-AUG-traceable cohort** (the same one supporting Figure 4 Panels C/D and Figure 5 Panel G):
- 1,050 NMD+/PTC+ comparator isoforms
- 116 NMD+/PTC− comparator isoforms
- 1,166 matched non-NMD Control comparator isoforms

Total: **2,332 isoforms.** No train/test split filtering; the comparison is "given these isoforms, what do the various models predict."

### 2.2 Gold standard

**Per-isoform mashr posterior-mean log fold-change of isoform abundance under SMG1i versus DMSO, 4-CT scope (AT, DD, FB, MV).** This is the same quantity that underlies our binary `nmd_responsive == TRUE` label but kept continuous to allow regression metrics that match what NMDetective and NMDEP report.

Higher mashr log-FC = stronger NMD substrate. Controls cluster near zero; PTC+ NMD substrates have high log-FC; PTC− NMD substrates also have high log-FC by construction (they're in the mashr-NMD set).

### 2.3 Scoring approach

**Apply each model's published feature definitions to features computed from OUR isoform structures**, not to canonical Ensembl/UCSC transcripts at lookup positions. This eliminates the "canonical-transcript-only" limitation of variant-level models and makes the comparison apples-to-apples at the isoform level.

For each of the 2,332 isoforms we compute four NMDetective-B features from its long-read-observed transcript structure:
- `InLastExon` — is the comparator stop codon in the last exon?
- `DistanceToStart` — nt from the comparator's start codon to its stop codon
- `ExonLength` — length of the exon containing the comparator stop codon
- `50ntToLastEJ` — is the stop ≤ 50 nt upstream of the last exon-exon junction?

For Controls (no premature stop) and NMD+/PTC− isoforms (no canonical PTC by the 50-nt rule), the features describe the natural stop codon. The rule-based models will therefore predict "no NMD" for these by construction — that's the expected behaviour and the comparison reveals it as a real limitation.

### 2.4 Models scored

**Included in v1:**
- **NMDetective-B** (Lindeboom 2019) — 4-rule decision tree, published thresholds 50 / 150 / 407 nt
- **NMDEP rule baseline** (Saadat 2025, Table 3) — same 4 rules re-thresholded to 49 / 120 / 355 nt
- **Our model** — predicted NMD probability per isoform (already in `predictions_atg500_stop500.tsv`)

**Deferred unless trivially runnable:**
- **NMDetective-A** (full random forest, ~71% R² in their hands) — requires their trained model object. Easy to score from the figshare resource IF that resource ships the model weights, not just the genome-wide predictions. Check first; if it's only the canonical-transcript predictions, defer.
- **Full NMDEP** (60-feature MLP with sequence embeddings + imputed half-life / ribosome loading / localization auxiliary models) — requires their GitHub code (must verify exists) plus running Orthrus sequence embeddings on our isoform sequences plus running their 3 auxiliary models. High infrastructure cost. Defer unless their code repo provides an inference API for arbitrary mRNA sequences.

### 2.5 Metrics

Match the conventions reported in NMDEP Table 4 (they report against the same NMD-efficacy gold standard from TCGA, so this is the natural comparison frame):
- Spearman correlation
- Pearson correlation
- R²
- Mean Absolute Error
- Root Mean Squared Error

Stratified by isoform subclass (PTC+, PTC−, Control) and pooled across all three.

### 2.6 Figure deliverable

**Two-panel supplemental figure:**
- **Panel A** — scatter plot per model, one subplot each: predicted NMD score (x) vs observed mashr posterior-mean log-FC (y), points coloured by subclass. Spearman correlation annotated.
- **Panel B** — bar chart of Spearman correlation per model, stratified by subclass (3 colours per model: PTC+ / PTC− / Control) plus pooled.

---

## 3. Pipeline

### 3.1 Inputs

| Source | What | Where |
|---|---|---|
| `ref_atg_analysis.rds` | Per-isoform `category`, `comp_stop_genomic`, `n_downstream_ejc` for the 1,166 scope | `results/.../data_mashr/analysis_cache/` |
| `structures.rds` | Per-isoform exon coordinates, junction positions | `results/.../data_mashr/` |
| `cds.rds` | Per-isoform start codon position | `results/.../data_mashr/` |
| Mashr per-cell-type DIE | Per-isoform mashr posterior mean + lfsr | `isocall_dge/mashr/nmd_mashr_die_*_2026.3.10.csv` |
| Our model predictions | Per-isoform predicted NMD probability | `NMD_orf_model_v5_4ct/results_4ct/predictions_atg500_stop500.tsv` |

### 3.2 Steps

1. **Cohort + gold-standard table** (`01_extract_our_isoforms.R`). For each of the 2,332 isoforms: gather `gene_id`, `comparator_isoform_id`, `subclass` (PTC+ / PTC− / Control), `mashr_posterior_mean_logfc` (the gold standard), and the four feature columns needed for NMDetective-B (`in_last_exon`, `distance_to_start_nt`, `exon_length_at_stop_nt`, `within_50nt_of_last_ej`). Write `isoforms_2026.6.20.tsv`.

2. **NMDetective-B scoring** (`02_score_nmdetective_b.R`). Apply the published 4-rule decision tree with thresholds 50 / 150 / 407 nt to each isoform. Output binary trigger/escape call AND a continuous NMD-efficacy estimate via the leaf-level NMD score from their Figure 1c (last-exon: 0.00; long-exon: 0.41; start-proximal: 0.12; 50nt-rule: 0.20; trigger: 0.65). Higher = stronger NMD trigger.

3. **NMDEP rule-baseline scoring** (`03_score_nmdep_rule_baseline.R`). Same 4-rule structure, re-thresholded to 49 / 120 / 355 nt (Table 3 of Saadat 2025). Same leaf-level NMD-efficacy assignment as NMDetective-B unless their paper specifies different leaf scores (re-check; if not, use the same leaf values).

4. **Our model scores**. Pull `prob` from `predictions_atg500_stop500.tsv` for each of our 2,332 isoforms.

5. **Metrics** (`04_compute_metrics.R`). For each model × stratum (PTC+ / PTC− / Control / pooled): Spearman, Pearson, R², MAE, RMSE against `mashr_posterior_mean_logfc`. Write `metrics_summary_2026.6.20.tsv`.

6. **Reproducible report** (`nmd_predictor_comparison.Rmd`). Inline R values pulled from the TSVs above; rendered HTML. Output: `nmd_predictor_comparison_2026.6.20.html`.

### 3.3 Outputs

Committed (datestamped per CLAUDE.md `yyyy.m.d` convention):
- `isoforms_2026.6.20.tsv` — per-isoform feature + gold-standard input
- `per_isoform_scores_2026.6.20.tsv` — predictions from each model side-by-side
- `metrics_summary_2026.6.20.tsv` — pooled and stratified Spearman / Pearson / R² / MAE / RMSE
- `nmd_predictor_comparison_2026.6.20.html` — rendered analysis report

Not committed:
- `external/` — only if NMDetective-A or full NMDEP are run

---

## 4. Open risks

1. **Cross-dataset domain shift.** TCGA-trained models tested on SMG1i-derived primary-cell data is a real domain shift — performance differences may reflect dataset rather than method. The paired Lindeboom-vs-Saadat-vs-ours comparison helps separate these (the two prior models are on the same training domain; we're the cross-domain target).

2. **NMD+/PTC− subclass behaviour.** Rule-based models will systematically miss the PTC− subgroup by construction (no canonical PTC → no rule fires → predicted no NMD). This is expected; the stratified-by-subclass metric (§2.5) will quantify it. The interpretation is "rule-based models cannot recover non-canonical NMD substrates" rather than "models perform poorly."

3. **NMDetective-A model availability.** If the trained random-forest weights aren't bundled in the figshare resource, NMDetective-A is deferred. Final v1 will state explicitly which variant was scored and which were deferred.

4. **NMDEP full-model deferral.** If their GitHub code isn't released or doesn't support inference on arbitrary mRNA sequences, the rule-only baseline is the only NMDEP variant we'll include in v1. State this explicitly in the supplement.

5. **Gold-standard noise.** mashr posterior means are shrunk estimates; the SE is non-trivial for low-abundance isoforms. Consider whether to weight points by isoform expression level, or to filter to a minimum-expression cohort. Defer this decision to a sensitivity analysis after the primary results land.

---

## 5. Sign-off

All decisions locked in 2026-06-20:
- §2.1 Cohort: all 2,332 isoforms in the n = 1,166 ref-AUG-traceable scope ✓
- §2.2 Gold standard: per-isoform mashr posterior mean log-FC under SMG1i (4-CT) ✓
- §2.3 Scoring: features from our isoforms, NOT canonical-transcript lookup ✓
- §2.4 Models: NMDetective-B and NMDEP rule baseline in v1; defer full models unless easy ✓
- §2.5 Metrics: Spearman / Pearson / R² / MAE / RMSE, pooled + stratified by subclass ✓
- §2.6 Figure: 2-panel (scatter + stratified-Spearman bars) ✓

Ready to write the README and `01_extract_our_isoforms.R`.
