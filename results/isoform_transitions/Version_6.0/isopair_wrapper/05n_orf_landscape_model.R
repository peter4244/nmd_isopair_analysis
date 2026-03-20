#!/usr/bin/env Rscript
# ==============================================================================
# 05n_orf_landscape_model.R
#
# ORF landscape NMD prediction models.
# Compares three approaches to predicting NMD from the ORF landscape:
#   1. Simple threshold: has any NMD-competent ORF?
#   2. Elastic net with ORF landscape features
#   3. Combined: ORF landscape + existing 5'UTR features
#
# All evaluated on holdout chr 1, 3, 5, 7 (same as unified model).
#
# Inputs:
#   - data_mashr/analysis_cache/orf_landscape.rds (from 05m)
#   - data_mashr/analysis_cache/utr5_features_all.rds (existing 5'UTR features)
#   - data_mashr/analysis_cache/unified_model.rds (for comparison)
#
# Output:
#   - data_mashr/analysis_cache/orf_landscape_model.rds
#
# Usage:
#   Rscript 05n_orf_landscape_model.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(pROC)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

OUTPUT_PATH <- "data_mashr/analysis_cache/orf_landscape_model.rds"

cat("=== ORF Landscape NMD Prediction Models ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Load data
# ==============================================================================
cat("Loading data...\n")

orf_df <- readRDS("data_mashr/analysis_cache/orf_landscape.rds")
unified <- readRDS("data_mashr/analysis_cache/unified_model.rds")

cat("  ORF landscape isoforms:", nrow(orf_df), "\n")
cat("  NMD:", sum(orf_df$is_nmd), "  Non-NMD:", sum(!orf_df$is_nmd), "\n")

# ==============================================================================
# 2. Train/test split (same holdout as unified model)
# ==============================================================================
holdout_chrs <- c("chr1", "chr3", "chr5", "chr7")

train <- orf_df[!orf_df$chr %in% holdout_chrs, ]
test <- orf_df[orf_df$chr %in% holdout_chrs, ]

cat("\n  Train:", nrow(train), "(NMD:", sum(train$is_nmd), ")\n")
cat("  Test:", nrow(test), "(NMD:", sum(test$is_nmd), ")\n")

# ==============================================================================
# 3. Model 1: Simple threshold — has any NMD-competent ORF?
# ==============================================================================
cat("\n=== Model 1: Simple Threshold ===\n")

pred_threshold_test <- as.integer(test$has_any_nmd_competent)
roc_threshold <- roc(test$is_nmd, pred_threshold_test, quiet = TRUE)
cat("  Test AUC:", round(auc(roc_threshold), 4), "\n")
thresh_coords <- pROC::coords(roc_threshold, x = "best", best.method = "youden")
cat("  Sensitivity:", round(thresh_coords$sensitivity, 4), "\n")
cat("  Specificity:", round(thresh_coords$specificity, 4), "\n")

# ==============================================================================
# 4. Model 2: Elastic net with ORF landscape features
# ==============================================================================
cat("\n=== Model 2: Elastic Net (ORF Landscape) ===\n")

orf_features <- c(
  "n_nmd_competent_capped",    # NMD-competent ORF count (capped at 20)
  "max_downstream_ejc_capped", # max downstream EJCs across all ORF stops (capped at 5)
  "n_nmd_competent_strong_kozak", # NMD-competent ORFs with strong Kozak
  "best_kozak_nmd",            # best Kozak score among NMD-competent ATGs
  "frac_nmd_competent",        # fraction of ORFs that are NMD-competent
  "n_orfs",                    # total ORFs (proxy for transcript complexity)
  "n_atgs",                    # total ATGs
  "n_junctions",               # total exon-exon junctions
  "tx_length"                  # transcript length
)

# Build model matrices
make_matrix <- function(data, features) {
  f <- as.formula(paste("~ 0 +", paste(features, collapse = " + ")))
  model.matrix(f, data = data)
}

x_train <- make_matrix(train, orf_features)
x_test <- make_matrix(test, orf_features)
y_train <- train$is_nmd
y_test <- test$is_nmd

# Scale
means <- colMeans(x_train, na.rm = TRUE)
sds <- apply(x_train, 2, sd, na.rm = TRUE)
sds[sds == 0] <- 1

x_train_s <- scale(x_train, center = means, scale = sds)
x_test_s <- scale(x_test, center = means, scale = sds)
x_train_s[is.na(x_train_s)] <- 0
x_test_s[is.na(x_test_s)] <- 0

set.seed(42)
cv_orf <- cv.glmnet(x_train_s, y_train, family = "binomial",
                     alpha = 0.5, nfolds = 10, type.measure = "auc")

pred_orf_train <- predict(cv_orf, x_train_s, s = "lambda.min", type = "response")[, 1]
pred_orf_test <- predict(cv_orf, x_test_s, s = "lambda.min", type = "response")[, 1]

roc_orf_train <- roc(y_train, pred_orf_train, quiet = TRUE)
roc_orf_test <- roc(y_test, pred_orf_test, quiet = TRUE)

cat("  Train AUC:", round(auc(roc_orf_train), 4), "\n")
cat("  Test AUC:", round(auc(roc_orf_test), 4), "\n")

# Feature importance
coefs <- as.matrix(coef(cv_orf, s = "lambda.min"))
nonzero <- coefs[coefs[, 1] != 0 & rownames(coefs) != "(Intercept)", , drop = FALSE]
cat("\n  Non-zero coefficients:\n")
for (j in order(-abs(nonzero[, 1]))) {
  cat(sprintf("    %s: %.4f\n", rownames(nonzero)[j], nonzero[j, 1]))
}

# ==============================================================================
# 5. Model 3: ORF landscape + downstream EJC from annotated CDS
# ==============================================================================
cat("\n=== Model 3: ORF Landscape + Annotated CDS EJC ===\n")

# Add the annotated-CDS downstream EJC count from the unified model's feature matrix
ptc <- readRDS("data_mashr/ptc.rds")
ptc_dejc <- ptc[, c("isoform_id", "n_downstream_ejcs")]
ptc_dejc$annotated_dejc <- pmin(ptc_dejc$n_downstream_ejcs, 5L)

orf_df_combined <- merge(orf_df, ptc_dejc[, c("isoform_id", "annotated_dejc")],
                          by = "isoform_id", all.x = TRUE)
orf_df_combined$annotated_dejc[is.na(orf_df_combined$annotated_dejc)] <- 0L

train_c <- orf_df_combined[!orf_df_combined$chr %in% holdout_chrs, ]
test_c <- orf_df_combined[orf_df_combined$chr %in% holdout_chrs, ]

combined_features <- c(orf_features, "annotated_dejc")

x_train_c <- make_matrix(train_c, combined_features)
x_test_c <- make_matrix(test_c, combined_features)

means_c <- colMeans(x_train_c, na.rm = TRUE)
sds_c <- apply(x_train_c, 2, sd, na.rm = TRUE)
sds_c[sds_c == 0] <- 1

x_train_cs <- scale(x_train_c, center = means_c, scale = sds_c)
x_test_cs <- scale(x_test_c, center = means_c, scale = sds_c)
x_train_cs[is.na(x_train_cs)] <- 0
x_test_cs[is.na(x_test_cs)] <- 0

set.seed(42)
cv_combined <- cv.glmnet(x_train_cs, train_c$is_nmd, family = "binomial",
                          alpha = 0.5, nfolds = 10, type.measure = "auc")

pred_comb_train <- predict(cv_combined, x_train_cs, s = "lambda.min", type = "response")[, 1]
pred_comb_test <- predict(cv_combined, x_test_cs, s = "lambda.min", type = "response")[, 1]

roc_comb_train <- roc(train_c$is_nmd, pred_comb_train, quiet = TRUE)
roc_comb_test <- roc(test_c$is_nmd, pred_comb_test, quiet = TRUE)

cat("  Train AUC:", round(auc(roc_comb_train), 4), "\n")
cat("  Test AUC:", round(auc(roc_comb_test), 4), "\n")

coefs_c <- as.matrix(coef(cv_combined, s = "lambda.min"))
nonzero_c <- coefs_c[coefs_c[, 1] != 0 & rownames(coefs_c) != "(Intercept)", , drop = FALSE]
cat("\n  Non-zero coefficients:\n")
for (j in order(-abs(nonzero_c[, 1]))) {
  cat(sprintf("    %s: %.4f\n", rownames(nonzero_c)[j], nonzero_c[j, 1]))
}

# ==============================================================================
# 6. Comparison summary
# ==============================================================================
cat("\n=== MODEL COMPARISON ===\n\n")

comparison <- data.frame(
  model = c("Threshold (any NMD-competent ORF)",
            "Elastic net (ORF landscape)",
            "Elastic net (ORF + annotated CDS EJC)",
            "Unified model (Section 2a)"),
  auc_test = round(c(
    as.numeric(auc(roc_threshold)),
    as.numeric(auc(roc_orf_test)),
    as.numeric(auc(roc_comb_test)),
    unified$elastic_net$step4$auc_test  # step4 = full model with 3'UTR
  ), 4),
  stringsAsFactors = FALSE
)

print(comparison)

# ==============================================================================
# 7. Save results
# ==============================================================================
results <- list(
  threshold = list(
    roc_test = roc_threshold,
    auc_test = as.numeric(auc(roc_threshold))
  ),
  orf_landscape = list(
    cv_fit = cv_orf,
    roc_train = roc_orf_train, roc_test = roc_orf_test,
    auc_train = as.numeric(auc(roc_orf_train)),
    auc_test = as.numeric(auc(roc_orf_test)),
    features = orf_features,
    scaling = list(means = means, sds = sds),
    n_train = nrow(train), n_test = nrow(test)
  ),
  combined = list(
    cv_fit = cv_combined,
    roc_train = roc_comb_train, roc_test = roc_comb_test,
    auc_train = as.numeric(auc(roc_comb_train)),
    auc_test = as.numeric(auc(roc_comb_test)),
    features = combined_features,
    scaling = list(means = means_c, sds = sds_c),
    n_train = nrow(train_c), n_test = nrow(test_c)
  ),
  comparison = comparison,
  orf_landscape_data = orf_df,
  metadata = list(
    holdout_chrs = holdout_chrs,
    run_timestamp = Sys.time(),
    population = nrow(orf_df)
  )
)

saveRDS(results, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
