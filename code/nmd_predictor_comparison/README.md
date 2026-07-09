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

## Prerequisite — full-cohort model predictions

Step 01 requires `predictions_all_atg500_stop500.tsv` (one row per isoform in the deep-learning model's H5 universe, split="all"). The trained model + H5 only live on Explorer (Northeastern GPU cluster), so this file is produced there and copied back. See `explorer_run/README.md`. Briefly:

```bash
# On Explorer
cd ~/cc/nmd_orf_model_v5_4ct
git pull
sbatch ... run_infer_all.py --config config.yaml --atg-window 500 --stop-window 500
# Output: results_4ct/predictions_all_atg500_stop500.tsv (~40k rows)

# From laptop, scp back into this folder
scp p.castaldi@login.explorer.northeastern.edu:~/cc/nmd_orf_model_v5_4ct/results_4ct/predictions_all_atg500_stop500.tsv .
```

## Run order

```bash
cd code/nmd_predictor_comparison

# 1. Build the per-isoform input table — features + gold standard + model prob
Rscript 01_extract_our_isoforms.R

# 2. Score each isoform with the published 4-rule NMDetective-B decision tree
Rscript 02_score_nmdetective_b.R

# 3. Score each isoform with the NMDEP-re-thresholded 4-rule baseline
Rscript 03_score_nmdep_rule_baseline.R

# 4. Compute pooled, stratified, and test-only-sensitivity metrics
Rscript 04_compute_metrics.R

# 5. Knit the reproducible report (HTML, self-contained, figure embedded)
RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64 \
  Rscript -e "rmarkdown::render('nmd_predictor_comparison.Rmd', \
                                 output_file='nmd_predictor_comparison_2026.6.20.html')"
```

Each step writes one datestamped TSV that the next step reads. The four R scripts together run in under a minute on a laptop. Heavy compute (model inference) is the prerequisite step on Explorer.

## Outputs

Committed (datestamped per CLAUDE.md `yyyy.m.d` convention):

| File | Produced by | Contents |
|---|---|---|
| `predictions_all_atg500_stop500.tsv` | `explorer_run/run_infer_all.py` (on Explorer) | Per-isoform deep-learning model score across the full H5 universe (~40k isoforms; `isoform_id`, `chr`, `h5_split`, `label`, `logit`, `prob`) |
| `isoforms_2026.6.20.tsv` | 01 | Per-isoform input row: `gene_id`, `comparator_isoform_id`, `subclass` (PTC+ / PTC− / Control), `stop_tx_pos_used`, the four NMDetective rule features, `mashr_posterior_mean_logfc`, `chr`, `h5_split`, `our_model_prob` |
| `nmdetective_b_scores_2026.6.20.tsv` | 02 | Per-isoform NMDetective-B leaf assignment (`nmdetective_b_score`, `_leaf`, `_call`) |
| `nmdep_rule_baseline_scores_2026.6.20.tsv` | 03 | Per-isoform NMDEP rule baseline leaf assignment (same columns, NMDEP thresholds) |
| `per_isoform_scores_2026.6.20.tsv` | 04 | Joined per-isoform scores across all three models + gold standard (one row per cohort isoform) |
| `metrics_summary_2026.6.20.tsv` | 04 | Per-model × per-stratum Spearman / Pearson / R² / MAE / RMSE — pooled, by-subclass, head-to-head intersection, and test-only sensitivity slice (`head-to-head:test:*`) |
| `nmd_predictor_comparison.Rmd` + `nmd_predictor_comparison_YYYY.M.D.html` | 5 (manual knit) | Reproducible HTML report with inline R values, methods, metrics table, and embedded figure. Primary analytical frame is test-only (chr1/3/5/7 held-out split) for fair comparison against the rule-based models which never saw any of our data. |

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
- Final paper figure: `figures/SupplementalFigures/SF42_ModelComparison/`
  (created after the analysis is run and we've agreed on the panel layout)
- Companion mainline figure: `figures/multipanel/figure5_dl_model/` (our
  model's training, architecture, and other interpretability analyses)
