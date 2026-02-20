#!/usr/bin/env Rscript
################################################################################
# Analyze Failures from 08_validate_reconstruction.R
################################################################################
#
# Investigates the 2,570 FAILs and 53 ERRORs from the NMD dataset validation.
# Categorizes failures by:
#   - Failure type (exon count mismatch direction, coordinate mismatch)
#   - Event type composition of failing pairs
#   - Novel vs GENCODE isoform involvement
#   - Exon count difference distribution
#   - Specific event types associated with failures
#
################################################################################

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   Analyze Reconstruction Failures                            ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"
setwd(base_dir)

# ==============================================================================
# 1. Load Data
# ==============================================================================

cat("Loading data...\n")
verification <- readRDS("data/reconstruction_verification.rds")
profiles <- readRDS("data/splicing_choice_profiles.rds")

cat(sprintf("  Verification results: %d\n", nrow(verification)))
cat(sprintf("  Profiles: %d\n", nrow(profiles)))

# ==============================================================================
# 2. Overall Summary
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SECTION 1: Overall Summary\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

status_summary <- verification %>%
  count(status) %>%
  mutate(pct = 100 * n / sum(n))

for (i in 1:nrow(status_summary)) {
  cat(sprintf("  %s: %d (%.1f%%)\n", status_summary$status[i],
              status_summary$n[i], status_summary$pct[i]))
}

# ==============================================================================
# 3. Failure Classification
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SECTION 2: Failure Classification\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

failures <- verification %>% filter(status == "FAIL")

# Parse failure reason
failures <- failures %>%
  mutate(
    failure_type = case_when(
      grepl("Exon count mismatch", reason) ~ "exon_count",
      grepl("mismatch", reason) ~ "coordinate",
      TRUE ~ "other"
    ),
    exon_diff = n_exons_reconstructed - n_exons_original,
    exon_diff_direction = case_when(
      exon_diff > 0 ~ "too_many_exons",
      exon_diff < 0 ~ "too_few_exons",
      exon_diff == 0 ~ "same_count_wrong_coords",
      TRUE ~ "unknown"
    )
  )

cat("Failure types:\n")
failures %>%
  count(failure_type, sort = TRUE) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  {for (i in 1:nrow(.)) cat(sprintf("  %s: %d (%.1f%%)\n", .$failure_type[i], .$n[i], .$pct[i]))}

cat("\nExon count direction:\n")
failures %>%
  count(exon_diff_direction, sort = TRUE) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  {for (i in 1:nrow(.)) cat(sprintf("  %s: %d (%.1f%%)\n", .$exon_diff_direction[i], .$n[i], .$pct[i]))}

cat("\nExon count difference distribution:\n")
diff_summary <- failures %>%
  filter(failure_type == "exon_count") %>%
  pull(exon_diff)

cat(sprintf("  Mean: %.1f\n", mean(diff_summary)))
cat(sprintf("  Median: %.0f\n", median(diff_summary)))
cat(sprintf("  Range: %d to %d\n", min(diff_summary), max(diff_summary)))
cat(sprintf("  IQR: %d to %d\n", quantile(diff_summary, 0.25), quantile(diff_summary, 0.75)))

cat("\nExon difference breakdown:\n")
failures %>%
  filter(failure_type == "exon_count") %>%
  count(exon_diff, sort = FALSE) %>%
  arrange(exon_diff) %>%
  {for (i in 1:min(nrow(.), 30)) cat(sprintf("  diff = %+d: %d pairs\n", .$exon_diff[i], .$n[i]))}

# ==============================================================================
# 4. Novel vs GENCODE Isoform Analysis
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SECTION 3: Novel (PB.) vs GENCODE (ENST) Isoforms\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

# Join verification with profile info
verification_enriched <- verification %>%
  left_join(
    profiles %>%
      select(gene_id, dominant_isoform_id, non_dominant_isoform_id,
             n_events, n_a5ss, n_a3ss, n_ir, n_ir_diff, n_partial_ir,
             n_se, n_missing_internal, n_alt_tss, n_alt_tes,
             tss_changed, tes_changed),
    by = c("gene_id",
           "dominant_isoform_id" = "dominant_isoform_id",
           "comparator_isoform_id" = "non_dominant_isoform_id")
  )

# Classify dominant and comparator isoform types
verification_enriched <- verification_enriched %>%
  mutate(
    dom_type = if_else(grepl("^PB\\.", dominant_isoform_id), "novel_PB", "GENCODE"),
    comp_type = if_else(grepl("^PB\\.", comparator_isoform_id), "novel_PB", "GENCODE"),
    pair_type = paste0(dom_type, " vs ", comp_type)
  )

cat("Pair type distribution (all):\n")
verification_enriched %>%
  count(pair_type) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  {for (i in 1:nrow(.)) cat(sprintf("  %s: %d (%.1f%%)\n", .$pair_type[i], .$n[i], .$pct[i]))}

cat("\nPass rate by pair type:\n")
verification_enriched %>%
  group_by(pair_type) %>%
  summarise(
    total = n(),
    pass = sum(status == "PASS"),
    fail = sum(status == "FAIL"),
    error = sum(status == "ERROR"),
    pass_rate = 100 * pass / total,
    .groups = "drop"
  ) %>%
  arrange(desc(pass_rate)) %>%
  {for (i in 1:nrow(.)) cat(sprintf("  %s: %d/%d PASS (%.1f%%), %d FAIL, %d ERROR\n",
    .$pair_type[i], .$pass[i], .$total[i], .$pass_rate[i], .$fail[i], .$error[i]))}

# GENCODE-only pair comparison (to compare with 99.6% benchmark)
gencode_only <- verification_enriched %>%
  filter(dom_type == "GENCODE" & comp_type == "GENCODE")

cat(sprintf("\nGENCODE-only pairs: %d total, %d PASS (%.1f%%), %d FAIL, %d ERROR\n",
            nrow(gencode_only),
            sum(gencode_only$status == "PASS"),
            100 * sum(gencode_only$status == "PASS") / nrow(gencode_only),
            sum(gencode_only$status == "FAIL"),
            sum(gencode_only$status == "ERROR")))

# ==============================================================================
# 5. Event Type Composition of Failures
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SECTION 4: Event Types in Failing vs Passing Pairs\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

event_cols <- c("n_a5ss", "n_a3ss", "n_ir", "n_ir_diff", "n_partial_ir",
                "n_se", "n_missing_internal", "n_alt_tss", "n_alt_tes")

# Compare event type prevalence
for (col in event_cols) {
  pass_mean <- mean(verification_enriched[[col]][verification_enriched$status == "PASS"], na.rm = TRUE)
  fail_mean <- mean(verification_enriched[[col]][verification_enriched$status == "FAIL"], na.rm = TRUE)
  pass_pct <- 100 * mean(verification_enriched[[col]][verification_enriched$status == "PASS"] > 0, na.rm = TRUE)
  fail_pct <- 100 * mean(verification_enriched[[col]][verification_enriched$status == "FAIL"] > 0, na.rm = TRUE)

  cat(sprintf("  %-20s  PASS: %.2f mean (%.0f%% have)  FAIL: %.2f mean (%.0f%% have)\n",
              col, pass_mean, pass_pct, fail_mean, fail_pct))
}

# Mean total events
cat(sprintf("\n  Mean n_events:  PASS = %.1f,  FAIL = %.1f\n",
            mean(verification_enriched$n_events[verification_enriched$status == "PASS"], na.rm = TRUE),
            mean(verification_enriched$n_events[verification_enriched$status == "FAIL"], na.rm = TRUE)))

# ==============================================================================
# 6. Event types that predict failure (logistic regression)
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SECTION 5: Event Types Predictive of Failure\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

# Binary flags for event presence
model_data <- verification_enriched %>%
  filter(status %in% c("PASS", "FAIL")) %>%
  mutate(
    is_fail = as.integer(status == "FAIL"),
    has_a5ss = as.integer(n_a5ss > 0),
    has_a3ss = as.integer(n_a3ss > 0),
    has_ir = as.integer(n_ir > 0),
    has_ir_diff = as.integer(n_ir_diff > 0),
    has_partial_ir = as.integer(n_partial_ir > 0),
    has_se = as.integer(n_se > 0),
    has_missing_internal = as.integer(n_missing_internal > 0),
    has_alt_tss = as.integer(n_alt_tss > 0),
    has_alt_tes = as.integer(n_alt_tes > 0),
    is_novel_dom = as.integer(dom_type == "novel_PB"),
    is_novel_comp = as.integer(comp_type == "novel_PB")
  ) %>%
  filter(!is.na(n_events))  # drop rows that didn't join

cat("Univariate logistic regression (event type → failure):\n\n")

predictors <- c("has_a5ss", "has_a3ss", "has_ir", "has_ir_diff", "has_partial_ir",
                "has_se", "has_missing_internal", "has_alt_tss", "has_alt_tes",
                "is_novel_dom", "is_novel_comp", "n_events")

for (pred in predictors) {
  fit <- glm(as.formula(paste("is_fail ~", pred)), data = model_data, family = "binomial")
  coef_summary <- summary(fit)$coefficients
  if (nrow(coef_summary) >= 2) {
    or <- exp(coef_summary[2, 1])
    p_val <- coef_summary[2, 4]
    sig <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else ""
    cat(sprintf("  %-25s OR = %6.2f  p = %.2e  %s\n", pred, or, p_val, sig))
  }
}

# ==============================================================================
# 7. Detailed failure examples (first 10)
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SECTION 6: Detailed Failure Examples (first 10)\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

fail_examples <- verification_enriched %>%
  filter(status == "FAIL") %>%
  head(10)

for (i in seq_len(nrow(fail_examples))) {
  f <- fail_examples[i, ]
  cat(sprintf("─── %d. %s ───\n", i, f$gene_id))
  cat(sprintf("  DOM: %s (%s)  COMP: %s (%s)\n",
              f$dominant_isoform_id, f$dom_type,
              f$comparator_isoform_id, f$comp_type))
  cat(sprintf("  Exons: %d original → %d reconstructed (diff = %+d)\n",
              f$n_exons_original, f$n_exons_reconstructed,
              f$n_exons_reconstructed - f$n_exons_original))
  cat(sprintf("  Events: %d total | A5SS:%d A3SS:%d IR:%d IR_diff:%d PIR:%d SE:%d MI:%d TSS:%d TES:%d\n",
              f$n_events, f$n_a5ss, f$n_a3ss, f$n_ir, f$n_ir_diff,
              f$n_partial_ir, f$n_se, f$n_missing_internal, f$n_alt_tss, f$n_alt_tes))
  cat(sprintf("  Reason: %s\n\n", f$reason))
}

# ==============================================================================
# 8. Zero-event failures
# ==============================================================================

cat("═══════════════════════════════════════════════════════════════════\n")
cat("SECTION 7: Zero-Event Pairs\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

zero_event_results <- verification_enriched %>%
  filter(n_events == 0 | is.na(n_events))

cat(sprintf("Pairs with zero events: %d\n", nrow(zero_event_results)))
if (nrow(zero_event_results) > 0) {
  cat("Status breakdown:\n")
  zero_event_results %>%
    count(status) %>%
    {for (i in 1:nrow(.)) cat(sprintf("  %s: %d\n", .$status[i], .$n[i]))}
}

# ==============================================================================
# 9. Multi-event complexity analysis
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SECTION 8: Failure Rate by Event Complexity\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

complexity_analysis <- verification_enriched %>%
  filter(status %in% c("PASS", "FAIL"), !is.na(n_events)) %>%
  mutate(
    event_bin = case_when(
      n_events == 0 ~ "0 events",
      n_events == 1 ~ "1 event",
      n_events == 2 ~ "2 events",
      n_events == 3 ~ "3 events",
      n_events <= 5 ~ "4-5 events",
      n_events <= 10 ~ "6-10 events",
      TRUE ~ "11+ events"
    ),
    event_bin = factor(event_bin, levels = c("0 events", "1 event", "2 events",
                                             "3 events", "4-5 events", "6-10 events", "11+ events"))
  ) %>%
  group_by(event_bin) %>%
  summarise(
    total = n(),
    fail = sum(status == "FAIL"),
    fail_rate = 100 * fail / total,
    .groups = "drop"
  )

for (i in 1:nrow(complexity_analysis)) {
  cat(sprintf("  %-14s  %5d pairs  %4d fail  (%.1f%%)\n",
              complexity_analysis$event_bin[i],
              complexity_analysis$total[i],
              complexity_analysis$fail[i],
              complexity_analysis$fail_rate[i]))
}

# ==============================================================================
# 10. ERROR analysis
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SECTION 9: ERROR Analysis\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

errors <- verification_enriched %>% filter(status == "ERROR")

cat(sprintf("Total errors: %d\n\n", nrow(errors)))

cat("Error reasons:\n")
errors %>%
  count(reason, sort = TRUE) %>%
  {for (i in 1:nrow(.)) cat(sprintf("  %s: %d\n", .$reason[i], .$n[i]))}

cat("\nError pair types:\n")
errors %>%
  count(pair_type, sort = TRUE) %>%
  {for (i in 1:nrow(.)) cat(sprintf("  %s: %d\n", .$pair_type[i], .$n[i]))}

# ==============================================================================
# Summary
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("KEY FINDINGS\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

# Top 3 most enriched event types in failures
cat("Top predictors of failure (by odds ratio):\n")
or_results <- list()
for (pred in c("has_a5ss", "has_a3ss", "has_ir", "has_ir_diff", "has_partial_ir",
               "has_se", "has_missing_internal", "has_alt_tss", "has_alt_tes",
               "is_novel_dom", "is_novel_comp")) {
  fit <- glm(as.formula(paste("is_fail ~", pred)), data = model_data, family = "binomial")
  coef_summary <- summary(fit)$coefficients
  if (nrow(coef_summary) >= 2) {
    or_results[[pred]] <- tibble(
      predictor = pred,
      or = exp(coef_summary[2, 1]),
      p_val = coef_summary[2, 4]
    )
  }
}
or_df <- bind_rows(or_results) %>% arrange(desc(or))
for (i in 1:min(5, nrow(or_df))) {
  cat(sprintf("  %d. %s (OR = %.2f, p = %.2e)\n", i, or_df$predictor[i], or_df$or[i], or_df$p_val[i]))
}

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("Analysis complete.\n")
