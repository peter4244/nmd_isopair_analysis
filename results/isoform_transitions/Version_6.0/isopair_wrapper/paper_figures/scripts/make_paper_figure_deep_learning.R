#!/usr/bin/env Rscript
# Paper Figure: Deep Learning Model Performance
# Panels: A) Progressive XGBoost AUC across feature steps (TD2 vs Ref-CDS)
#         B) ROC overlay (Combined, TD2 step3, Ref-CDS step3, Ref EJC-only)
#         C) Combined model SHAP top features
#         D) Matched-population AUC comparison

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
source("visualization_functions.R")
source("analysis_functions.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(pROC)
})

paper_theme <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        plot.margin = margin(5, 8, 5, 5))

pal <- c("NMD" = "#ef8a62", "Control" = "#67a9cf")

# Family palette: TD2 vs Ref-CDS
family_pal <- c("TD2" = "#91bfdb", "Ref-CDS" = "#b2182b")

# ---- Load data ----
mc <- readRDS("data_mashr/analysis_cache/model_comparison.rds")
# Note: unified_model.rds has elastic-net AUCs but only for TD2 features.
# model_comparison.rds$auc_summary holds the matched progressive AUCs for
# TD2 and Ref-CDS feature families (XGBoost, paralog-free holdout) — the
# only place both families are progressively scored on identical splits.

# ---- Panel A: Progressive AUC ----
auc_summary <- mc$auc_summary
# Drop the "Combined" row from the progressive plot
auc_progressive <- auc_summary[auc_summary$model_set %in% c("TD2", "Ref-CDS"), ]

step_labels <- c(
  "Step 1: PTC only"      = "Step 1\nStop codon",
  "Step 2: + Start codon" = "Step 2\n+ Start codon",
  "Step 3: + ORF features"= "Step 3\n+ ORF features",
  "Step 4: + 3'UTR"       = "Step 4\n+ 3' UTR length"
)
auc_progressive$step_short <- step_labels[auc_progressive$step]
auc_progressive$step_short <- factor(auc_progressive$step_short, levels = unname(step_labels))
auc_progressive$model_set <- factor(auc_progressive$model_set, levels = c("TD2", "Ref-CDS"))

# Long format with split-eval column
panelA_df <- rbind(
  data.frame(model_set = auc_progressive$model_set,
             step_short = auc_progressive$step_short,
             auc = auc_progressive$auc_test,
             split = "Test",
             stringsAsFactors = FALSE),
  data.frame(model_set = auc_progressive$model_set,
             step_short = auc_progressive$step_short,
             auc = auc_progressive$auc_train,
             split = "Train",
             stringsAsFactors = FALSE)
)
panelA_df$split <- factor(panelA_df$split, levels = c("Train", "Test"))

n_td2_full <- mc$td2_models$step3$n_test
n_ref_match <- mc$ref_models$step3$n_test

pA <- ggplot(panelA_df,
             aes(x = step_short, y = auc,
                 color = model_set, group = interaction(model_set, split),
                 shape = split, linetype = split)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 3.2) +
  scale_color_manual(values = family_pal) +
  scale_shape_manual(values = c("Train" = 1, "Test" = 16)) +
  scale_linetype_manual(values = c("Train" = "dashed", "Test" = "solid")) +
  scale_y_continuous(limits = c(0.5, 1.0),
                      breaks = seq(0.5, 1.0, 0.1)) +
  labs(title = "Progressive XGBoost AUC by Feature Step",
       subtitle = sprintf("Test sets: TD2 n = %s, Ref-CDS n = %s (paralog-free holdout)",
                          format(n_td2_full, big.mark = ","),
                          format(n_ref_match, big.mark = ",")),
       x = NULL, y = "AUC",
       color = "Feature family", shape = "Split", linetype = "Split") +
  paper_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 9))

# ---- Panel B: ROC overlay ----
roc_specs <- list(
  list(roc = mc$combined_model$roc_test,
       name = "Combined (24 feat.)",
       color = "#1a1a1a"),
  list(roc = mc$td2_models$step3$roc_test,
       name = "TD2 step3 (11 feat.)",
       color = "#91bfdb"),
  list(roc = mc$ref_models$step3$roc_test,
       name = "Ref-CDS step3 (11 feat.)",
       color = "#b2182b"),
  list(roc = mc$ref_models$step1$roc_test,
       name = "Ref-CDS EJC-only (1 feat.)",
       color = "#4393c3")
)

roc_df <- do.call(rbind, lapply(roc_specs, function(spec) {
  auc_val <- as.numeric(auc(spec$roc))
  data.frame(
    fpr = 1 - spec$roc$specificities,
    tpr = spec$roc$sensitivities,
    model = sprintf("%s (AUC = %.3f)", spec$name, auc_val),
    stringsAsFactors = FALSE
  )
}))
roc_df$model <- factor(roc_df$model, levels = unique(roc_df$model))
roc_colors <- setNames(sapply(roc_specs, `[[`, "color"), levels(roc_df$model))

pB <- ggplot(roc_df, aes(x = fpr, y = tpr, color = model)) +
  geom_abline(linetype = "dashed", color = "grey70") +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = roc_colors) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.02, 0)) +
  labs(title = "Holdout ROC: Model Comparison",
       subtitle = sprintf("TD2 / Combined evaluated on n = %s; Ref-CDS on n = %s",
                          format(mc$td2_models$step3$n_test, big.mark = ","),
                          format(mc$ref_models$step3$n_test, big.mark = ",")),
       x = "False Positive Rate", y = "True Positive Rate", color = NULL) +
  paper_theme +
  theme(legend.position = c(0.62, 0.22),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.key.width = unit(1.0, "cm"),
        legend.text = element_text(size = 8.5))

# ---- Panel C: Combined model SHAP top features ----
shap_mat <- mc$combined_model$shap_test
mean_abs <- colMeans(abs(shap_mat))
shap_df <- data.frame(
  feature = names(mean_abs),
  mean_abs_shap = unname(mean_abs),
  stringsAsFactors = FALSE
)
shap_df$source <- ifelse(grepl("^ref_", shap_df$feature), "Ref-CDS", "TD2")
shap_df$label <- mc$feature_labels[shap_df$feature]
shap_df <- shap_df[order(-shap_df$mean_abs_shap), ]
top_n <- 10
shap_top <- head(shap_df, top_n)
shap_top$label <- factor(shap_top$label, levels = rev(shap_top$label))
shap_top$source <- factor(shap_top$source, levels = c("TD2", "Ref-CDS"))

n_shap <- nrow(shap_mat)

pC <- ggplot(shap_top, aes(x = mean_abs_shap, y = label, fill = source)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", mean_abs_shap)),
            hjust = -0.15, size = 3, color = "grey20") +
  scale_fill_manual(values = family_pal) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = sprintf("Combined Model SHAP: Top %d Features", top_n),
       subtitle = sprintf("Mean |SHAP| over n = %s holdout isoforms",
                          format(n_shap, big.mark = ",")),
       x = "Mean |SHAP value|", y = NULL, fill = "Feature family") +
  paper_theme +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 9))

# ---- Panel D: Matched-population AUC comparison ----
# TD2 step3 evaluated on:
#   (i) full test set (all classified coding isoforms)
#   (ii) matched subset where Ref-CDS features are available
# Ref-CDS step3 evaluated on its native subset (same as matched).
matched_idx <- which(!is.na(mc$test$ref_utr5_excluded) &
                       mc$test$ref_utr5_excluded == FALSE)
pred_matched <- mc$td2_models$step3$pred_test[matched_idx]
truth_matched <- mc$test$is_nmd[matched_idx]
roc_td2_matched <- roc(truth_matched, pred_matched, quiet = TRUE)
auc_td2_matched <- as.numeric(auc(roc_td2_matched))

panelD_df <- data.frame(
  bar_label = factor(
    c("TD2 (full pop.)", "TD2 (matched pop.)", "Ref-CDS (matched pop.)"),
    levels = c("TD2 (full pop.)", "TD2 (matched pop.)", "Ref-CDS (matched pop.)")),
  family = factor(c("TD2", "TD2", "Ref-CDS"), levels = c("TD2", "Ref-CDS")),
  auc = c(mc$td2_models$step3$auc_test,
          auc_td2_matched,
          mc$ref_models$step3$auc_test),
  n = c(mc$td2_models$step3$n_test,
        length(matched_idx),
        mc$ref_models$step3$n_test)
)
panelD_df$bar_text <- sprintf("AUC = %.3f\nn = %s",
                                panelD_df$auc,
                                format(panelD_df$n, big.mark = ","))

pD <- ggplot(panelD_df, aes(x = bar_label, y = auc, fill = family)) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_text(aes(label = bar_text), vjust = -0.25, size = 3.4, lineheight = 0.95) +
  scale_fill_manual(values = family_pal) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.2),
                      expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Matched-Population Comparison (Step-3 Models)",
       subtitle = "Most of the apparent Ref-CDS advantage is population selection;\non matched populations the gap is small",
       x = NULL, y = "Test AUC", fill = "Feature family") +
  paper_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 9))

# ---- Compose ----
combined <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A",
    theme = theme(plot.margin = margin(10, 10, 10, 10))) &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# ---- Validate composite layout ----
fig_w <- 14; fig_h <- 10
source("~/.claude/utils/validate_composite_layout.R")
validate_composite_layout(combined, fig_width = fig_w, fig_height = fig_h)

# ---- Save ----
ggsave("paper_figures/FigX-DeepLearning_generated.png", combined,
       width = fig_w, height = fig_h, dpi = 300, bg = "white")
ggsave("paper_figures/FigX-DeepLearning_generated.pdf", combined,
       width = fig_w, height = fig_h, bg = "white", device = cairo_pdf)

cat("Saved: paper_figures/FigX-DeepLearning_generated.png and .pdf\n")
