# NMD predictor comparison — README

How to reproduce the head-to-head comparison of our deep-learning NMD model
against the variant-level NMD predictors NMDetective-B (Lindeboom et al.,
*Nat Genet* 2019) and NMDEP-rule-baseline (Saadat & Fellay, *arXiv* 2025).

See `METHODS.md` for the pre-registered analysis design and decision log.

## Dependencies

R packages: `data.table`, `dplyr`, `pROC` (for AUC sanity checks).

No external downloads required for v1. NMDetective-A and full NMDEP are
deferred unless their model weights / inference code prove easy to run
(see `METHODS.md` §2.4).

## Run order

```bash
cd code/nmd_predictor_comparison

# 1. Build the per-isoform input table — features and gold standard
Rscript 01_extract_our_isoforms.R

# 2. Score each isoform with the published 4-rule NMDetective-B decision tree
Rscript 02_score_nmdetective_b.R

# 3. Score each isoform with the NMDEP-re-thresholded 4-rule baseline
Rscript 03_score_nmdep_rule_baseline.R

# 4. Compute pooled and stratified metrics
Rscript 04_compute_metrics.R
```

Each step writes one datestamped TSV that the next step reads. The full
pipeline runs in under a minute on a laptop — there's no external data
download, no model inference, and no heavy compute.

## Outputs

Committed (datestamped per CLAUDE.md `yyyy.m.d` convention):

| File | Produced by | Contents |
|---|---|---|
| `isoforms_2026.6.20.tsv` | 01 | Per-isoform input row: `gene_id`, `comparator_isoform_id`, `subclass` (PTC+ / PTC− / Control), `mashr_posterior_mean_logfc`, the four NMDetective rule features, model probability |
| `per_isoform_scores_2026.6.20.tsv` | 03 | Per-isoform scores from each model side-by-side (`nmdetective_b_score`, `nmdep_rule_score`, `our_model_prob`) and the gold standard column |
| `metrics_summary_2026.6.20.tsv` | 04 | Per-model × per-stratum Spearman, Pearson, R², MAE, RMSE — matches NMDEP Table 4 columns |

## What "gold standard" means here

The gold standard is the **mashr posterior-mean log fold-change under SMG1i
vs DMSO**, 4-CT scope (alveolar type 2, large-airway epithelial, fibroblast,
microvascular endothelial). Higher values indicate stronger NMD substrates:
the isoform's abundance recovers more under SMG1 inhibition. This is the
continuous quantity that underlies our binary `nmd_responsive == TRUE` label
but kept continuous so we can use the regression metrics NMDetective and
NMDEP report.

This is a different gold standard from the one the prior models were trained
on — they used per-PTC NMD efficacy estimated from somatic-mutation
allele-specific expression in TCGA tumors. Applying their models to ours is
therefore a cross-dataset generalization test, not a head-to-head on the
training domain.

## Cross-references

- Analysis pre-registration and decision log: `METHODS.md`
- Final paper figure: `figures/SupplementalFigures/ModelComparison/`
  (created after the analysis is run and we've agreed on the panel layout)
- Companion mainline figure: `figures/multipanel/figure5_dl_model/` (our
  model's training, architecture, and other interpretability analyses)
