#!/usr/bin/env Rscript
# ==============================================================================
# 05v_model_comparison.R
#
# Combined NMD prediction model: TD2 CDS vs Reference CDS vs Combined.
# Uses XGBoost gradient-boosted trees with native NA handling.
#
# Loads pre-computed features from:
#   - unified_model.rds (TD2 features, ~60K isoforms)
#   - ref_cds_features_all.rds (ref-CDS features, ~60K isoforms)
#   - paralog_genes.rds (expressed paralogs for test set filtering)
#
# Constructs matched train/test sets and fits:
#   Model Set 1 — TD2 progressive (steps 1-4) — all isoforms
#   Model Set 2 — Ref-CDS progressive (steps 1-4) — ref_atg_available only
#   Model Set 3 — Combined (all TD2 + all ref-CDS features) — all isoforms
#
# XGBoost handles NA natively, so the combined model includes all isoforms
# (ref-CDS columns are NA for ref_atg_lost isoforms).
#
# Output:
#   - data_mashr/analysis_cache/model_comparison.rds
#
# Usage:
#   Rscript 05v_model_comparison.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(xgboost)
  library(pROC)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

OUTPUT_PATH <- "data_mashr/analysis_cache/model_comparison.rds"

cat("=== Combined Prediction Model: TD2 vs Reference CDS (XGBoost) ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Load pre-computed data
# ==============================================================================
cat("Loading data...\n")

unified <- readRDS("data_mashr/analysis_cache/unified_model.rds")
ref_cds <- readRDS("data_mashr/analysis_cache/ref_cds_features_all.rds")
paralogs <- readRDS("data_mashr/analysis_cache/paralog_genes.rds")

td2_df <- unified$df  # Full TD2 feature matrix
ref_features <- ref_cds$features

cat("  TD2 feature matrix:", nrow(td2_df), "isoforms\n")
cat("  Ref-CDS features:", nrow(ref_features), "isoforms\n")
cat("  Paralog leakage genes:", length(paralogs$leakage_genes), "\n")

# ==============================================================================
# 2. Construct populations
# ==============================================================================
cat("\nConstructing populations...\n")

# TD2 features needed (from unified_model.rds$df)
td2_cols <- c("isoform_id", "is_nmd", "chr", "gene_id",
              "downstream_ejc",
              "atg_count", "atg_density", "atg_strong_kozak", "atg_orphan_count",
              "uorf_count_overlapping", "uorf_count_inframe", "uorf_count_outframe",
              "uorf_longest_nt", "utr5_orf_coverage", "stop_density",
              "log_utr3_length")

# Ref-CDS features needed
ref_cols <- c("isoform_id", "ref_downstream_ejc", "ref_atg_available",
              "ref_atg_count", "ref_atg_density", "ref_atg_strong_kozak",
              "ref_atg_orphan_count",
              "ref_uorf_count_overlapping", "ref_uorf_count_inframe",
              "ref_uorf_count_outframe", "ref_uorf_longest_nt",
              "ref_utr5_orf_coverage", "ref_stop_density",
              "ref_log_utr3_length", "ref_utr5_excluded")

# Merge TD2 and ref-CDS features — LEFT JOIN to keep all TD2 isoforms
td2_sub <- td2_df[, td2_cols]
ref_sub <- ref_features[, intersect(ref_cols, names(ref_features))]

merged <- left_join(td2_sub, ref_sub, by = "isoform_id")

cat("  After left join:", nrow(merged), "isoforms\n")

# Full population: all isoforms with TD2 features (ref-CDS may be NA)
merged_all <- merged
cat("  Full population (TD2 complete):", nrow(merged_all), "\n")

# Ref-complete population: for ref-CDS-only models
merged_complete <- merged %>%
  filter(
    ref_atg_available == TRUE,
    !is.na(ref_downstream_ejc),
    !is.na(ref_utr5_excluded) & !ref_utr5_excluded,
    !is.na(ref_log_utr3_length)
  )

cat("  With complete ref-CDS features:", nrow(merged_complete), "\n")
cat("  Full population NMD:", sum(merged_all$is_nmd), "\n")
cat("  Full population Non-NMD:", sum(!merged_all$is_nmd), "\n")

# ==============================================================================
# 3. Train/test split with paralog removal from test
# ==============================================================================
cat("\nSplitting train/test...\n")

holdout_chrs <- c("chr1", "chr3", "chr5", "chr7")
paralog_gene_set <- paralogs$leakage_genes

# Full population splits (for TD2 and combined models)
train_all <- merged_all %>% filter(!chr %in% holdout_chrs)
test_all_full <- merged_all %>% filter(chr %in% holdout_chrs)
test_all_paralogs <- test_all_full %>% filter(gene_id %in% paralog_gene_set)
test_all_clean <- test_all_full %>% filter(!gene_id %in% paralog_gene_set)

# Ref-complete splits (for ref-CDS-only models)
train_ref <- merged_complete %>% filter(!chr %in% holdout_chrs)
test_ref_all <- merged_complete %>% filter(chr %in% holdout_chrs)
test_ref_clean <- test_ref_all %>% filter(!gene_id %in% paralog_gene_set)

n_paralog_removed_genes <- length(unique(test_all_paralogs$gene_id))
n_paralog_removed_isoforms <- nrow(test_all_paralogs)

cat("  Train (full):", nrow(train_all), "(NMD:", sum(train_all$is_nmd), ")\n")
cat("  Test (full, before paralog):", nrow(test_all_full),
    "(NMD:", sum(test_all_full$is_nmd), ")\n")
cat("  Paralog genes removed from test:", n_paralog_removed_genes,
    "(", n_paralog_removed_isoforms, "isoforms)\n")
cat("  Test (full, final):", nrow(test_all_clean),
    "(NMD:", sum(test_all_clean$is_nmd), ")\n")
cat("  Train (ref-complete):", nrow(train_ref),
    "(NMD:", sum(train_ref$is_nmd), ")\n")
cat("  Test (ref-complete, final):", nrow(test_ref_clean),
    "(NMD:", sum(test_ref_clean$is_nmd), ")\n")

# ==============================================================================
# 4. Feature definitions
# ==============================================================================

# TD2 feature sets (progressive)
td2_step1 <- "downstream_ejc"
td2_step2 <- c(td2_step1, "atg_density", "atg_count", "atg_strong_kozak")
td2_step3 <- c(td2_step2,
               "uorf_count_overlapping", "uorf_longest_nt",
               "uorf_count_inframe", "uorf_count_outframe",
               "utr5_orf_coverage", "stop_density", "atg_orphan_count")
td2_step4 <- c(td2_step3, "log_utr3_length")

# Ref-CDS feature sets (progressive, same structure)
ref_step1 <- "ref_downstream_ejc"
ref_step2 <- c(ref_step1, "ref_atg_density", "ref_atg_count", "ref_atg_strong_kozak")
ref_step3 <- c(ref_step2,
               "ref_uorf_count_overlapping", "ref_uorf_longest_nt",
               "ref_uorf_count_inframe", "ref_uorf_count_outframe",
               "ref_utr5_orf_coverage", "ref_stop_density", "ref_atg_orphan_count")
ref_step4 <- c(ref_step3, "ref_log_utr3_length")

# Combined: all TD2 + all ref-CDS features
combined_features <- c(td2_step4, ref_step4)

# Feature display labels
feature_labels <- c(
  # TD2
  downstream_ejc = "TD2: Downstream EJC count",
  atg_density = "TD2: 5\u2032UTR ATG density",
  atg_count = "TD2: 5\u2032UTR ATG count",
  atg_strong_kozak = "TD2: Strong Kozak ATG count",
  atg_orphan_count = "TD2: Orphan ATG count",
  uorf_count_overlapping = "TD2: Overlapping uORF count",
  uorf_longest_nt = "TD2: Longest uORF (nt)",
  uorf_count_inframe = "TD2: In-frame uORF count",
  uorf_count_outframe = "TD2: Out-of-frame uORF count",
  utr5_orf_coverage = "TD2: 5\u2032UTR ORF coverage",
  stop_density = "TD2: 5\u2032UTR stop density",
  log_utr3_length = "TD2: 3\u2032UTR length (log bp)",
  # Ref-CDS
  ref_downstream_ejc = "Ref: Downstream EJC count",
  ref_atg_density = "Ref: 5\u2032UTR ATG density",
  ref_atg_count = "Ref: 5\u2032UTR ATG count",
  ref_atg_strong_kozak = "Ref: Strong Kozak ATG count",
  ref_atg_orphan_count = "Ref: Orphan ATG count",
  ref_uorf_count_overlapping = "Ref: Overlapping uORF count",
  ref_uorf_longest_nt = "Ref: Longest uORF (nt)",
  ref_uorf_count_inframe = "Ref: In-frame uORF count",
  ref_uorf_count_outframe = "Ref: Out-of-frame uORF count",
  ref_utr5_orf_coverage = "Ref: 5\u2032UTR ORF coverage",
  ref_stop_density = "Ref: 5\u2032UTR stop density",
  ref_log_utr3_length = "Ref: 3\u2032UTR length (log bp)"
)

# ==============================================================================
# 5. XGBoost model fitting helper
# ==============================================================================

fit_xgboost <- function(train_df, test_df, features, step_name) {
  cat("\n--- ", step_name, " ---\n")

  # Build matrices (XGBoost handles NA natively via missing param)
  x_train <- as.matrix(train_df[, features, drop = FALSE])
  y_train <- as.integer(train_df$is_nmd)
  x_test <- as.matrix(test_df[, features, drop = FALSE])
  y_test <- as.integer(test_df$is_nmd)

  dtrain <- xgb.DMatrix(data = x_train, label = y_train, missing = NA)
  dtest <- xgb.DMatrix(data = x_test, label = y_test, missing = NA)

  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 4,
    eta = 0.1,
    min_child_weight = 10,
    subsample = 0.8,
    colsample_bytree = 0.8,
    seed = 42
  )

  # Find optimal nrounds via CV with early stopping
  set.seed(42)
  cv_result <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 500,
    nfold = 10,
    early_stopping_rounds = 20,
    verbose = 0
  )

  best_nrounds <- cv_result$early_stop$best_iteration
  cat("  Best nrounds (CV):", best_nrounds, "\n")

  # Train final model
  set.seed(42)
  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_nrounds,
    evals = list(train = dtrain, test = dtest),
    verbose = 0
  )

  # Predictions
  pred_train <- predict(xgb_model, dtrain)
  pred_test <- predict(xgb_model, dtest)

  # ROC
  roc_train <- roc(y_train, pred_train, quiet = TRUE)
  roc_test <- roc(y_test, pred_test, quiet = TRUE)

  cat("  Train AUC:", round(auc(roc_train), 4), "\n")
  cat("  Test AUC:", round(auc(roc_test), 4), "\n")

  # SHAP values (predcontrib = TRUE returns n_features + 1 columns, last is BIAS)
  shap_raw <- predict(xgb_model, dtest, predcontrib = TRUE)
  # Drop BIAS column (last column)
  shap_test <- shap_raw[, seq_len(ncol(shap_raw) - 1), drop = FALSE]
  colnames(shap_test) <- features

  # Store raw feature values for beeswarm color scale
  shap_feature_values <- x_test

  list(
    xgb_model = xgb_model,
    roc_train = roc_train, roc_test = roc_test,
    auc_train = as.numeric(auc(roc_train)),
    auc_test = as.numeric(auc(roc_test)),
    pred_train = pred_train, pred_test = pred_test,
    features = features,
    shap_test = shap_test,
    shap_feature_values = shap_feature_values,
    best_nrounds = best_nrounds,
    n_train = nrow(train_df), n_test = nrow(test_df)
  )
}

# ==============================================================================
# 6. Model Set 1 — TD2 CDS progressive (all isoforms)
# ==============================================================================
cat("\n========================================")
cat("\n  MODEL SET 1: TD2 CDS Annotations")
cat("\n========================================\n")

td2_m1 <- fit_xgboost(train_all, test_all_clean, td2_step1,
                       "TD2 Step 1: Downstream EJC")
td2_m2 <- fit_xgboost(train_all, test_all_clean, td2_step2,
                       "TD2 Step 2: + Start codon context")
td2_m3 <- fit_xgboost(train_all, test_all_clean, td2_step3,
                       "TD2 Step 3: + ORF features")
td2_m4 <- fit_xgboost(train_all, test_all_clean, td2_step4,
                       "TD2 Step 4: + 3'UTR length")

# ==============================================================================
# 7. Model Set 2 — Reference CDS progressive (ref-complete only)
# ==============================================================================
cat("\n========================================")
cat("\n  MODEL SET 2: Reference CDS")
cat("\n========================================\n")

ref_m1 <- fit_xgboost(train_ref, test_ref_clean, ref_step1,
                       "Ref Step 1: Downstream EJC")
ref_m2 <- fit_xgboost(train_ref, test_ref_clean, ref_step2,
                       "Ref Step 2: + Start codon context")
ref_m3 <- fit_xgboost(train_ref, test_ref_clean, ref_step3,
                       "Ref Step 3: + ORF features")
ref_m4 <- fit_xgboost(train_ref, test_ref_clean, ref_step4,
                       "Ref Step 4: + 3'UTR length")

# ==============================================================================
# 8. Model Set 3 — Combined (TD2 + Ref-CDS, all isoforms)
# ==============================================================================
cat("\n========================================")
cat("\n  MODEL SET 3: Combined (TD2 + Ref-CDS)")
cat("\n========================================\n")

combined_m <- fit_xgboost(train_all, test_all_clean, combined_features,
                           "Combined: All TD2 + all ref-CDS features")

# Build SHAP importance table (replaces combined_coefficients)
combined_shap_importance <- data.frame(
  feature = combined_features,
  mean_abs_shap = colMeans(abs(combined_m$shap_test)),
  source = ifelse(grepl("^ref_", combined_features), "Ref-CDS", "TD2"),
  stringsAsFactors = FALSE
) %>% arrange(desc(mean_abs_shap))

cat("\n  Combined model feature importance (mean |SHAP|):\n")
cat("  TD2 features with non-zero importance:",
    sum(combined_shap_importance$mean_abs_shap[combined_shap_importance$source == "TD2"] > 0.001), "\n")
cat("  Ref-CDS features with non-zero importance:",
    sum(combined_shap_importance$mean_abs_shap[combined_shap_importance$source == "Ref-CDS"] > 0.001), "\n")

for (j in seq_len(nrow(combined_shap_importance))) {
  cat(sprintf("    %s: %.4f [%s]\n",
              combined_shap_importance$feature[j],
              combined_shap_importance$mean_abs_shap[j],
              combined_shap_importance$source[j]))
}

# ==============================================================================
# 9. AUC summary table
# ==============================================================================
cat("\n\n=== AUC SUMMARY ===\n\n")

auc_summary <- data.frame(
  model_set = c(rep("TD2", 4), rep("Ref-CDS", 4), "Combined"),
  step = c("Step 1: PTC only", "Step 2: + Start codon", "Step 3: + ORF features",
           "Step 4: + 3'UTR",
           "Step 1: PTC only", "Step 2: + Start codon", "Step 3: + ORF features",
           "Step 4: + 3'UTR",
           "All features"),
  n_features = c(1, length(td2_step2), length(td2_step3), length(td2_step4),
                  1, length(ref_step2), length(ref_step3), length(ref_step4),
                  length(combined_features)),
  n_train = c(rep(nrow(train_all), 4), rep(nrow(train_ref), 4), nrow(train_all)),
  n_test = c(rep(nrow(test_all_clean), 4), rep(nrow(test_ref_clean), 4),
             nrow(test_all_clean)),
  auc_train = round(c(td2_m1$auc_train, td2_m2$auc_train, td2_m3$auc_train, td2_m4$auc_train,
                       ref_m1$auc_train, ref_m2$auc_train, ref_m3$auc_train, ref_m4$auc_train,
                       combined_m$auc_train), 4),
  auc_test = round(c(td2_m1$auc_test, td2_m2$auc_test, td2_m3$auc_test, td2_m4$auc_test,
                      ref_m1$auc_test, ref_m2$auc_test, ref_m3$auc_test, ref_m4$auc_test,
                      combined_m$auc_test), 4),
  stringsAsFactors = FALSE
)

print(auc_summary, right = FALSE)

# ==============================================================================
# 10. Save all results
# ==============================================================================
cat("\nSaving results...\n")

results <- list(
  # Model objects
  td2_models = list(step1 = td2_m1, step2 = td2_m2, step3 = td2_m3, step4 = td2_m4),
  ref_models = list(step1 = ref_m1, step2 = ref_m2, step3 = ref_m3, step4 = ref_m4),
  combined_model = combined_m,

  # Summary tables
  auc_summary = auc_summary,
  combined_shap_importance = combined_shap_importance,

  # Data — use full population as train/test (report accesses mc$train, mc$test)
  train = train_all,
  test = test_all_clean,
  merged_complete = merged_all,  # full population (preserved for compatibility)

  # Feature metadata
  feature_labels = feature_labels,
  feature_sets = list(
    td2 = list(step1 = td2_step1, step2 = td2_step2, step3 = td2_step3, step4 = td2_step4),
    ref = list(step1 = ref_step1, step2 = ref_step2, step3 = ref_step3, step4 = ref_step4),
    combined = combined_features
  ),

  # Filtering metadata
  metadata = list(
    run_timestamp = Sys.time(),
    holdout_chrs = holdout_chrs,
    n_merged = nrow(merged),
    n_complete = nrow(merged_all),
    n_ref_complete = nrow(merged_complete),
    n_train = nrow(train_all),
    n_train_ref = nrow(train_ref),
    n_test_before_paralog = nrow(test_all_full),
    n_paralog_removed_genes = n_paralog_removed_genes,
    n_paralog_removed_isoforms = n_paralog_removed_isoforms,
    n_test_final = nrow(test_all_clean),
    n_test_ref_final = nrow(test_ref_clean),
    n_nmd_train = sum(train_all$is_nmd),
    n_nmd_test = sum(test_all_clean$is_nmd),
    model_type = "xgboost",
    note = "TD2 and combined models use full population; ref-CDS models restricted to ref_atg_available. XGBoost handles NA natively."
  )
)

# Sanity checks
stopifnot(nrow(results$combined_model$shap_test) == nrow(test_all_clean))
stopifnot(ncol(results$combined_model$shap_test) == length(combined_features))
stopifnot(all(colnames(results$combined_model$shap_test) == combined_features))
stopifnot(nrow(results$combined_model$shap_feature_values) == nrow(test_all_clean))
stopifnot(length(results$combined_model$pred_test) == nrow(test_all_clean))

saveRDS(results, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
