#!/usr/bin/env Rscript
# ==============================================================================
# 05l_unified_model.R
#
# Unified prediction model: can isoform structural features predict NMD status?
# Progressive elastic net with chromosome-based holdout validation.
#
# Model steps (progressive — each adds a feature category):
#   Step 1: Stop codon position — downstream_ejc_category (0 / 1 / 2–3 / ≥4)
#   Step 2: + Start codon context — atg_density, atg_count, atg_strong_kozak
#   Step 3: + ORF features — uorf_count_overlapping, uorf_longest_nt, etc.
#
# Population: ALL mashr-classified NMD and non-NMD coding isoforms with
# complete structural features (~60K isoforms).
#
# Train/test split: holdout chromosomes 1, 3, 5, 7 (~26% of data).
# Feature selection and model tuning use ONLY training data.
#
# Outputs elastic net models with regularized coefficients and AUC progression.
#
# Feature name mapping (R variable → paper display name):
#   downstream_ejc_category  → Downstream EJC count (0 / 1 / 2–3 / ≥4)
#   uorf_count_overlapping   → Overlapping uORF count
#   uorf_longest_nt          → Longest uORF (nt)
#   uorf_count_inframe       → In-frame uORF count
#   uorf_count_outframe      → Out-of-frame uORF count
#   atg_density              → 5'UTR ATG density (per 100 bp)
#   atg_count                → 5'UTR ATG count
#   atg_strong_kozak         → Strong Kozak ATG count
#   atg_orphan_count         → Orphan ATG count
#   utr5_orf_coverage        → 5'UTR ORF coverage (%)
#   stop_density             → 5'UTR stop density (per 100 bp)
#
# Input files:
#   - data_mashr/nmd_classification.rds
#   - data_mashr/cds.rds
#   - data_mashr/structures.rds
#   - data_mashr/ptc.rds
#   - data_mashr/analysis_cache/utr5_features_all.rds
#
# Output:
#   - data_mashr/analysis_cache/unified_model.rds
#
# Usage:
#   Rscript 05l_unified_model.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(pROC)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

OUTPUT_PATH <- "data_mashr/analysis_cache/unified_model.rds"

cat("=== Unified NMD Prediction Model ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# Feature display names (for figures and tables)
# ==============================================================================
feature_labels <- c(
  downstream_ejc_category = "Downstream EJC count",
  # Factor levels get their own labels in model output:
  "downstream_ejc_category1"   = "Downstream EJC: 1",
  "downstream_ejc_category2-3" = "Downstream EJC: 2\u20133",
  "downstream_ejc_category4+"  = "Downstream EJC: \u22654",
  uorf_count_overlapping = "Overlapping uORF count",
  uorf_longest_nt        = "Longest uORF (nt)",
  uorf_count_inframe     = "In-frame uORF count",
  uorf_count_outframe    = "Out-of-frame uORF count",
  atg_density            = "5\u2032UTR ATG density (per 100 bp)",
  atg_count              = "5\u2032UTR ATG count",
  atg_strong_kozak       = "Strong Kozak ATG count",
  atg_orphan_count       = "Orphan ATG count",
  utr5_orf_coverage      = "5\u2032UTR ORF coverage (%)",
  stop_density           = "5\u2032UTR stop density (per 100 bp)"
)

# ==============================================================================
# 1. Load data and build feature matrix
# ==============================================================================
cat("Loading data...\n")

nmd_class <- readRDS("data_mashr/nmd_classification.rds")
cds <- readRDS("data_mashr/cds.rds")
structures <- readRDS("data_mashr/structures.rds")
ptc <- readRDS("data_mashr/ptc.rds")
utr5_all <- readRDS("data_mashr/analysis_cache/utr5_features_all.rds")

# --- Define NMD outcome ---
nmd_ids <- nmd_class$all_samples$nmd
non_nmd_ids <- nmd_class$all_samples$non_nmd
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]

nmd_coding <- intersect(nmd_ids, coding_ids)
non_nmd_coding <- intersect(non_nmd_ids, coding_ids)

cat("  NMD coding isoforms:", length(nmd_coding), "\n")
cat("  Non-NMD coding isoforms:", length(non_nmd_coding), "\n")

outcome <- data.frame(
  isoform_id = c(nmd_coding, non_nmd_coding),
  is_nmd = c(rep(1L, length(nmd_coding)), rep(0L, length(non_nmd_coding))),
  stringsAsFactors = FALSE
)

# --- Stop codon features: downstream EJC category ---
ptc_features <- ptc %>%
  filter(isoform_id %in% outcome$isoform_id) %>%
  select(isoform_id, n_downstream_ejcs) %>%
  mutate(
    downstream_ejc_category = factor(
      case_when(
        n_downstream_ejcs == 0 ~ "0",
        n_downstream_ejcs == 1 ~ "1",
        n_downstream_ejcs %in% 2:3 ~ "2-3",
        n_downstream_ejcs >= 4 ~ "4+"
      ),
      levels = c("0", "1", "2-3", "4+")
    )
  )

cat("  PTC features:", nrow(ptc_features), "isoforms\n")
cat("  EJC category distribution:\n")
print(table(ptc_features$downstream_ejc_category))

# --- 5'UTR features (full population from 05k) ---
utr5_feat <- utr5_all$isoform_features %>%
  filter(!excluded, isoform_id %in% outcome$isoform_id) %>%
  transmute(
    isoform_id,
    atg_count            = n_atg,
    atg_density          = atg_density,
    atg_strong_kozak     = n_strong_kozak_atg,
    atg_orphan_count     = n_atg_without_orf,
    uorf_count_overlapping = n_orfs_overlapping,
    uorf_count_inframe   = n_orfs_inframe,
    uorf_count_outframe  = n_orfs_outframe,
    uorf_longest_nt      = longest_orf_nt,
    utr5_orf_coverage    = pct_utr5_in_orfs,
    stop_density         = stop_density
  )

cat("  5'UTR features:", nrow(utr5_feat), "isoforms\n")

# --- Chromosome for train/test split ---
chr_map <- structures %>%
  filter(isoform_id %in% outcome$isoform_id) %>%
  select(isoform_id, chr, gene_id)

# ==============================================================================
# 2. Assemble feature matrix
# ==============================================================================
cat("\nAssembling feature matrix...\n")

df <- outcome %>%
  inner_join(ptc_features, by = "isoform_id") %>%
  inner_join(utr5_feat, by = "isoform_id") %>%
  inner_join(chr_map, by = "isoform_id")

cat("  Complete feature matrix:", nrow(df), "isoforms\n")
cat("    NMD:", sum(df$is_nmd), "  Non-NMD:", sum(!df$is_nmd), "\n")

# Population coverage
cat("\n  Population coverage:\n")
cat("    Classified coding:", nrow(outcome), "\n")
cat("    With PTC features:", nrow(ptc_features), "\n")
cat("    With 5'UTR features:", nrow(utr5_feat), "\n")
cat("    Complete (intersection):", nrow(df),
    sprintf("(%.1f%%)\n", 100 * nrow(df) / nrow(outcome)))
cat("    Excluded (no 5'UTR or ATG validation):",
    nrow(outcome) - nrow(df), "\n")

# ==============================================================================
# 3. Train/test split
# ==============================================================================
holdout_chrs <- c("chr1", "chr3", "chr5", "chr7")

train <- df %>% filter(!chr %in% holdout_chrs)
test <- df %>% filter(chr %in% holdout_chrs)

cat("\n  Train:", nrow(train), "(NMD:", sum(train$is_nmd), ")\n")
cat("  Test:", nrow(test), "(NMD:", sum(test$is_nmd), ")\n")

# ==============================================================================
# 4. Progressive models
# ==============================================================================
cat("\n=== Progressive Model Building ===\n")

# Feature sets for each step
step1_features <- "downstream_ejc_category"  # factor
step2_features <- c(step1_features, "atg_density", "atg_count", "atg_strong_kozak")
step3_features <- c(step2_features,
                    "uorf_count_overlapping", "uorf_longest_nt",
                    "uorf_count_inframe", "uorf_count_outframe",
                    "utr5_orf_coverage", "stop_density", "atg_orphan_count")

# 5'UTR-only (no stop codon features) for complementarity
utr5_only_features <- setdiff(step3_features, step1_features)

# --- Helper: build model matrix handling factors ---
make_model_matrix <- function(data, features) {
  f <- as.formula(paste("~ 0 +", paste(features, collapse = " + ")))
  model.matrix(f, data = data)
}

# --- Helper: fit elastic net on model matrix ---
fit_elastic_net <- function(train_df, test_df, features, step_name, alpha = 0.5) {
  cat("\n--- ", step_name, " ---\n")

  x_train <- make_model_matrix(train_df, features)
  y_train <- train_df$is_nmd
  x_test <- make_model_matrix(test_df, features)
  y_test <- test_df$is_nmd

  # Scale numeric columns (not factor dummies)
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

# Progressive elastic net
en_step1 <- fit_elastic_net(train, test, step1_features, "Step 1: Stop codon position")
en_step2 <- fit_elastic_net(train, test, step2_features, "Step 2: + Start codon context")
en_step3 <- fit_elastic_net(train, test, step3_features, "Step 3: + ORF features (unified)")

# Complementarity
en_utr5_only <- fit_elastic_net(train, test, utr5_only_features, "5'UTR only (no stop codon)")

# ==============================================================================
# 5. Score all isoforms with unified elastic net model
# ==============================================================================
cat("\n=== Scoring all isoforms ===\n")

# Score all isoforms with the unified elastic net (step 3)
x_all <- make_model_matrix(df, step3_features)
x_all_s <- scale(x_all,
                 center = en_step3$scaling$means,
                 scale = en_step3$scaling$sds)
x_all_s[is.na(x_all_s)] <- 0

pred_all <- predict(en_step3$cv_fit, x_all_s,
                    s = "lambda.min", type = "response")[, 1]

scored_isoforms <- data.frame(
  isoform_id = df$isoform_id,
  gene_id = df$gene_id,
  is_nmd = df$is_nmd,
  pred_prob = pred_all,
  chr = df$chr,
  in_holdout = df$chr %in% holdout_chrs,
  stringsAsFactors = FALSE
)

cat("  Scored", nrow(scored_isoforms), "isoforms\n")
cat("  Mean pred_prob (NMD):",
    round(mean(scored_isoforms$pred_prob[scored_isoforms$is_nmd == 1]), 3), "\n")
cat("  Mean pred_prob (non-NMD):",
    round(mean(scored_isoforms$pred_prob[scored_isoforms$is_nmd == 0]), 3), "\n")

# ==============================================================================
# 7. Summary
# ==============================================================================
cat("\n=== SUMMARY ===\n\n")

auc_summary <- data.frame(
  step = c("1_stop_codon", "2_start_codon", "3_unified", "utr5_only"),
  description = c("Downstream EJC count",
                   "+ ATG density/count/Kozak",
                   "+ uORF/ORF features",
                   "5'UTR features only (no stop codon)"),
  n_train = c(en_step1$n_train, en_step2$n_train,
              en_step3$n_train, en_utr5_only$n_train),
  n_test = c(en_step1$n_test, en_step2$n_test,
             en_step3$n_test, en_utr5_only$n_test),
  auc_train_en = c(en_step1$auc_train, en_step2$auc_train,
                   en_step3$auc_train, en_utr5_only$auc_train),
  auc_test_en = c(en_step1$auc_test, en_step2$auc_test,
                  en_step3$auc_test, en_utr5_only$auc_test),
  stringsAsFactors = FALSE
)

print(auc_summary, digits = 4, right = FALSE)

cat("\nPopulation:", nrow(df), "isoforms with complete features\n")
cat("Holdout chromosomes: 1, 3, 5, 7\n")

# ==============================================================================
# 8. Save results
# ==============================================================================
results <- list(
  # Elastic net models
  elastic_net = list(
    step1 = en_step1,
    step2 = en_step2,
    step3 = en_step3,
    utr5_only = en_utr5_only
  ),

  # AUC summary table
  auc_summary = auc_summary,

  # Scored isoforms
  scored_isoforms = scored_isoforms,

  # Feature matrix
  df = df,

  # Feature metadata
  feature_labels = feature_labels,

  # Metadata
  metadata = list(
    holdout_chrs = holdout_chrs,
    run_timestamp = Sys.time(),
    step_features = list(
      step1 = step1_features,
      step2 = step2_features,
      step3 = step3_features,
      utr5_only = utr5_only_features
    ),
    population = nrow(df),
    note = "Full population: all classified coding isoforms with complete features"
  )
)

saveRDS(results, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
