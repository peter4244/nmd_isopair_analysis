#!/usr/bin/env Rscript
# ==============================================================================
# 05v_model_comparison.R
#
# Combined NMD prediction model: TD2 CDS vs Reference CDS vs Combined.
#
# Loads pre-computed features from:
#   - unified_model.rds (TD2 features, ~60K isoforms)
#   - ref_cds_features_all.rds (ref-CDS features, ~60K isoforms)
#   - paralog_genes.rds (expressed paralogs for test set filtering)
#
# Constructs matched train/test sets and fits:
#   Model Set 1 — TD2 progressive (steps 1-4)
#   Model Set 2 — Ref-CDS progressive (steps 1-4)
#   Model Set 3 — Combined (all TD2 + all ref-CDS features)
#
# All models trained and evaluated on the SAME isoform population: the
# intersection of isoforms with both TD2 and ref-CDS features.
# Test set: holdout chromosomes 1, 3, 5, 7 with paralog genes removed.
#
# Output:
#   - data_mashr/analysis_cache/model_comparison.rds
#
# Usage:
#   Rscript 05v_model_comparison.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(pROC)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

OUTPUT_PATH <- "data_mashr/analysis_cache/model_comparison.rds"

cat("=== Combined Prediction Model: TD2 vs Reference CDS ===\n")
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
# 2. Construct matched population
# ==============================================================================
cat("\nConstructing matched population...\n")

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

# Merge TD2 and ref-CDS features
td2_sub <- td2_df[, td2_cols]
ref_sub <- ref_features[, intersect(ref_cols, names(ref_features))]

merged <- inner_join(td2_sub, ref_sub, by = "isoform_id")

cat("  After inner join:", nrow(merged), "isoforms\n")

# Filter to isoforms with ref-CDS features available:
# ref ATG must be available AND 5'UTR features not excluded
merged_complete <- merged %>%
  filter(
    ref_atg_available == TRUE,
    !is.na(ref_downstream_ejc),
    !is.na(ref_utr5_excluded) & !ref_utr5_excluded,
    !is.na(ref_log_utr3_length)
  )

cat("  With complete ref-CDS features:", nrow(merged_complete), "\n")
cat("    NMD:", sum(merged_complete$is_nmd), "\n")
cat("    Non-NMD:", sum(!merged_complete$is_nmd), "\n")

# ==============================================================================
# 3. Train/test split with paralog removal from test
# ==============================================================================
cat("\nSplitting train/test...\n")

holdout_chrs <- c("chr1", "chr3", "chr5", "chr7")
paralog_gene_set <- paralogs$leakage_genes

train <- merged_complete %>% filter(!chr %in% holdout_chrs)
test_all <- merged_complete %>% filter(chr %in% holdout_chrs)

# Remove paralog genes from test set only
test_with_paralogs <- test_all %>% filter(gene_id %in% paralog_gene_set)
test <- test_all %>% filter(!gene_id %in% paralog_gene_set)

n_paralog_removed_genes <- length(unique(test_with_paralogs$gene_id))
n_paralog_removed_isoforms <- nrow(test_with_paralogs)

cat("  Train:", nrow(train), "(NMD:", sum(train$is_nmd), ")\n")
cat("  Test (before paralog removal):", nrow(test_all),
    "(NMD:", sum(test_all$is_nmd), ")\n")
cat("  Paralog genes removed from test:", n_paralog_removed_genes,
    "(", n_paralog_removed_isoforms, "isoforms)\n")
cat("  Test (final):", nrow(test), "(NMD:", sum(test$is_nmd), ")\n")

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
# 5. Model fitting helpers
# ==============================================================================

make_model_matrix <- function(data, features) {
  f <- as.formula(paste("~ 0 +", paste(features, collapse = " + ")))
  model.matrix(f, data = data)
}

fit_glm <- function(train_df, test_df, features, step_name) {
  cat("\n--- ", step_name, " (logistic regression) ---\n")

  f <- as.formula(paste("is_nmd ~", paste(features, collapse = " + ")))
  fit <- glm(f, data = train_df, family = binomial())

  pred_train <- predict(fit, train_df, type = "response")
  pred_test <- predict(fit, test_df, type = "response")

  roc_train <- roc(train_df$is_nmd, pred_train, quiet = TRUE)
  roc_test <- roc(test_df$is_nmd, pred_test, quiet = TRUE)

  cat("  Train AUC:", round(auc(roc_train), 4), "\n")
  cat("  Test AUC:", round(auc(roc_test), 4), "\n")

  list(
    cv_fit = fit,
    roc_train = roc_train, roc_test = roc_test,
    auc_train = as.numeric(auc(roc_train)),
    auc_test = as.numeric(auc(roc_test)),
    pred_train = pred_train, pred_test = pred_test,
    features = features,
    scaling = NULL,
    n_train = nrow(train_df), n_test = nrow(test_df)
  )
}

fit_elastic_net <- function(train_df, test_df, features, step_name, alpha = 0.5) {
  cat("\n--- ", step_name, " ---\n")

  x_train <- make_model_matrix(train_df, features)
  y_train <- train_df$is_nmd
  x_test <- make_model_matrix(test_df, features)
  y_test <- test_df$is_nmd

  means <- colMeans(x_train, na.rm = TRUE)
  sds <- apply(x_train, 2, sd, na.rm = TRUE)
  sds[sds == 0] <- 1

  x_train_s <- scale(x_train, center = means, scale = sds)
  x_test_s <- scale(x_test, center = means, scale = sds)
  x_train_s[is.na(x_train_s)] <- 0
  x_test_s[is.na(x_test_s)] <- 0

  set.seed(42)
  cv_fit <- cv.glmnet(x_train_s, y_train, family = "binomial",
                       alpha = alpha, nfolds = 10, type.measure = "auc")

  pred_train <- predict(cv_fit, x_train_s, s = "lambda.min", type = "response")[, 1]
  pred_test <- predict(cv_fit, x_test_s, s = "lambda.min", type = "response")[, 1]

  roc_train <- roc(y_train, pred_train, quiet = TRUE)
  roc_test <- roc(y_test, pred_test, quiet = TRUE)

  cat("  Train AUC:", round(auc(roc_train), 4), "\n")
  cat("  Test AUC:", round(auc(roc_test), 4), "\n")

  # Extract coefficients
  coefs <- as.matrix(coef(cv_fit, s = "lambda.min"))
  nonzero <- coefs[coefs[, 1] != 0 & rownames(coefs) != "(Intercept)", , drop = FALSE]
  cat("  Non-zero features:", nrow(nonzero), "/", length(features), "\n")

  list(
    cv_fit = cv_fit,
    roc_train = roc_train, roc_test = roc_test,
    auc_train = as.numeric(auc(roc_train)),
    auc_test = as.numeric(auc(roc_test)),
    pred_train = pred_train, pred_test = pred_test,
    features = features,
    scaling = list(means = means, sds = sds),
    n_train = nrow(train_df), n_test = nrow(test_df)
  )
}

# ==============================================================================
# 6. Model Set 1 — TD2 CDS progressive
# ==============================================================================
cat("\n========================================")
cat("\n  MODEL SET 1: TD2 CDS Annotations")
cat("\n========================================\n")

td2_m1 <- fit_glm(train, test, td2_step1, "TD2 Step 1: Downstream EJC")
td2_m2 <- fit_elastic_net(train, test, td2_step2, "TD2 Step 2: + Start codon context")
td2_m3 <- fit_elastic_net(train, test, td2_step3, "TD2 Step 3: + ORF features")
td2_m4 <- fit_elastic_net(train, test, td2_step4, "TD2 Step 4: + 3'UTR length")

# ==============================================================================
# 7. Model Set 2 — Reference CDS progressive
# ==============================================================================
cat("\n========================================")
cat("\n  MODEL SET 2: Reference CDS")
cat("\n========================================\n")

ref_m1 <- fit_glm(train, test, ref_step1, "Ref Step 1: Downstream EJC")
ref_m2 <- fit_elastic_net(train, test, ref_step2, "Ref Step 2: + Start codon context")
ref_m3 <- fit_elastic_net(train, test, ref_step3, "Ref Step 3: + ORF features")
ref_m4 <- fit_elastic_net(train, test, ref_step4, "Ref Step 4: + 3'UTR length")

# ==============================================================================
# 8. Model Set 3 — Combined (TD2 + Ref-CDS)
# ==============================================================================
cat("\n========================================")
cat("\n  MODEL SET 3: Combined (TD2 + Ref-CDS)")
cat("\n========================================\n")

combined_m <- fit_elastic_net(train, test, combined_features,
                               "Combined: All TD2 + all ref-CDS features")

# Extract combined model coefficients for diagnostic
combined_coefs <- as.matrix(coef(combined_m$cv_fit, s = "lambda.min"))
combined_nonzero <- combined_coefs[
  combined_coefs[, 1] != 0 & rownames(combined_coefs) != "(Intercept)", ,
  drop = FALSE]

cat("\n  Combined model surviving features:\n")
coef_df <- data.frame(
  feature = rownames(combined_nonzero),
  coefficient = combined_nonzero[, 1],
  stringsAsFactors = FALSE
) %>% arrange(-abs(coefficient))

# Tag CDS source
coef_df$source <- ifelse(grepl("^ref_", coef_df$feature), "Ref-CDS", "TD2")

cat("  TD2 features surviving:", sum(coef_df$source == "TD2"), "\n")
cat("  Ref-CDS features surviving:", sum(coef_df$source == "Ref-CDS"), "\n")

for (j in seq_len(nrow(coef_df))) {
  cat(sprintf("    %s: %.4f [%s]\n",
              coef_df$feature[j], coef_df$coefficient[j], coef_df$source[j]))
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
  n_train = c(rep(nrow(train), 9)),
  n_test = c(rep(nrow(test), 9)),
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
  combined_coefficients = coef_df,

  # Data
  train = train,
  test = test,
  merged_complete = merged_complete,

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
    n_complete = nrow(merged_complete),
    n_train = nrow(train),
    n_test_before_paralog = nrow(test_all),
    n_paralog_removed_genes = n_paralog_removed_genes,
    n_paralog_removed_isoforms = n_paralog_removed_isoforms,
    n_test_final = nrow(test),
    n_nmd_train = sum(train$is_nmd),
    n_nmd_test = sum(test$is_nmd),
    note = "All models trained/tested on same population (intersection of TD2 and ref-CDS features)"
  )
)

saveRDS(results, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
