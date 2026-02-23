#!/usr/bin/env Rscript
# 06_rmats_frame_disruption.R
#
# rMATS frame disruption and PTC-inducing splice event analysis.
#
# Are PTC-inducing splice events (those that disrupt the reading frame
# within CDS) preferentially affected by Smg1i treatment?
#
# Uses short-read rMATS differential splicing results (Smg1i vs DMSO,
# 6 cell types) as a complement to the long-read PTC analysis (Scripts 01-05).
#
# Input:
#   results/rds/allrmats_jc_2025.8.20.rds   (junction counts, primary)
#   data/isoform_cds_metadata.rds            (CDS boundaries)
#   data/isoform_structures.rds              (isoform → gene_id, strand)
#
# Output (results/ptc/results/):
#   rmats_event_classification.tsv           (all events with annotations)
#   rmats_enrichment_by_celltype.tsv         (Fisher tests per cell type)
#   rmats_enrichment_by_event_type.tsv       (stratified by event type)
#   rmats_direction_tests.tsv                (direction bias tests)
#   rmats_magnitude_tests.tsv                (|dPSI| comparisons)
#   rmats_summary_stats.tsv                  (event counts by class)
#
# Output (results/ptc/figures/):
#   rmats_dpsi_by_frame_class.pdf
#   rmats_enrichment_forest.pdf
#
# Run from: results/isoform_transitions/Version_6.0/

suppressPackageStartupMessages({
  library(tidyverse)
  library(maser)
  library(GenomicRanges)
  library(metafor)
})

cat("\n")
cat("==================================================================\n")
cat("   06: rMATS Frame Disruption & PTC-Inducing Splice Events\n")
cat("==================================================================\n\n")

# Configurable thresholds
FDR_THRESH <- 0.05
DPSI_THRESH <- 0.05

out_dir <- "results/ptc/results"
fig_dir <- "results/ptc/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Cell type name mapping: rMATS → pipeline
CT_MAP <- c(AT2 = "AT", DD = "DD", DD_ALI = "DD_ALI",
            DO_ALI = "DO", FB = "FB", MV = "MV")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Load and Preprocess rMATS Data
# ═══════════════════════════════════════════════════════════════════════════════

cat("── Section 1: Loading rMATS data ──────────────────────────────────\n\n")

rmats_jc <- readRDS("../../../results/rds/allrmats_jc_2025.8.20.rds")
cat(sprintf("  Cell types: %s\n", paste(names(rmats_jc), collapse = ", ")))

# Verify dPSI direction convention from first object
obj0 <- rmats_jc[[1]]
cond_labels <- slot(obj0, "conditions")
cat(sprintf("  Conditions: cond1 = %s, cond2 = %s\n", cond_labels[1], cond_labels[2]))
cat("  IncLevelDifference = PSI(cond1) - PSI(cond2)\n")

# Extract all events into tidy tibbles
EVENT_TYPES <- c("SE", "A5SS", "A3SS", "MXE", "RI")

extract_events <- function(maser_obj, cell_type) {
  results <- list()
  for (etype in EVENT_TYPES) {
    events <- slot(maser_obj, paste0(etype, "_events")) %>% as_tibble()
    stats  <- slot(maser_obj, paste0(etype, "_stats"))  %>% as_tibble()
    gr     <- slot(maser_obj, paste0(etype, "_gr"))

    # Merge events + stats on ID
    merged <- inner_join(events, stats, by = "ID")

    # Extract coordinates from GRanges (already 1-based from Maser)
    if (etype == "SE") {
      target <- gr[["exon_target"]]
      merged$event_chr   <- as.character(seqnames(target))
      merged$event_start <- start(target)
      merged$event_end   <- end(target)
      merged$event_width <- width(target)
      # Flanking exons for context
      merged$upstream_end   <- end(gr[["exon_upstream"]])
      merged$downstream_start <- start(gr[["exon_downstream"]])
    } else if (etype %in% c("A5SS", "A3SS")) {
      long_gr  <- gr[["exon_long"]]
      short_gr <- gr[["exon_short"]]
      merged$event_chr   <- as.character(seqnames(long_gr))
      merged$long_start  <- start(long_gr)
      merged$long_end    <- end(long_gr)
      merged$long_width  <- width(long_gr)
      merged$short_start <- start(short_gr)
      merged$short_end   <- end(short_gr)
      merged$short_width <- width(short_gr)
      merged$event_start <- pmin(start(long_gr), start(short_gr))
      merged$event_end   <- pmax(end(long_gr), end(short_gr))
      merged$width_diff  <- abs(merged$long_width - merged$short_width)
    } else if (etype == "MXE") {
      e1 <- gr[["exon_1"]]
      e2 <- gr[["exon_2"]]
      merged$event_chr   <- as.character(seqnames(e1))
      merged$exon1_start <- start(e1)
      merged$exon1_end   <- end(e1)
      merged$exon1_width <- width(e1)
      merged$exon2_start <- start(e2)
      merged$exon2_end   <- end(e2)
      merged$exon2_width <- width(e2)
      merged$event_start <- pmin(start(e1), start(e2))
      merged$event_end   <- pmax(end(e1), end(e2))
      merged$width_diff  <- abs(merged$exon1_width - merged$exon2_width)
    } else if (etype == "RI") {
      up_gr   <- gr[["exon_upstream"]]
      down_gr <- gr[["exon_downstream"]]
      merged$event_chr   <- as.character(seqnames(up_gr))
      # Intron = gap between upstream exon end and downstream exon start
      merged$intron_start <- end(up_gr) + 1L
      merged$intron_end   <- start(down_gr) - 1L
      merged$intron_width <- merged$intron_end - merged$intron_start + 1L
      merged$event_start  <- merged$intron_start
      merged$event_end    <- merged$intron_end
    }

    merged$event_type <- etype
    merged$cell_type  <- CT_MAP[cell_type]
    # Strip gene ID version for matching
    merged$gene_id_unversioned <- sub("[.].*", "", merged$GeneID)

    results[[etype]] <- merged
  }
  bind_rows(results)
}

all_events <- map2_dfr(rmats_jc, names(rmats_jc), extract_events)
cat(sprintf("  Total events extracted: %s across %d cell types\n",
            format(nrow(all_events), big.mark = ","), length(rmats_jc)))

# Event counts per cell type
ct_counts <- all_events %>% count(cell_type, event_type) %>%
  pivot_wider(names_from = event_type, values_from = n)
cat("\n  Events per cell type:\n")
print(as.data.frame(ct_counts), row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Build Gene-Level CDS Reference
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Section 2: Building gene-level CDS reference ────────────────────\n\n")

cds_meta <- readRDS("data/isoform_cds_metadata.rds")
iso_str  <- readRDS("data/isoform_structures.rds")

# Join CDS + structure, GENCODE coding only
gencode_coding <- cds_meta %>%
  filter(coding_status == "coding", grepl("^ENST", isoform_id)) %>%
  inner_join(
    iso_str %>% select(isoform_id, gene_id, seqnames, strand),
    by = "isoform_id"
  )
cat(sprintf("  GENCODE coding isoforms: %s\n",
            format(nrow(gencode_coding), big.mark = ",")))

# Gene-level CDS envelope: min(cds_start) to max(cds_stop)
gene_cds <- gencode_coding %>%
  group_by(gene_id, seqnames) %>%
  summarise(
    cds_min = min(cds_start),
    cds_max = max(cds_stop),
    n_isoforms = n(),
    .groups = "drop"
  )
cat(sprintf("  Gene-level CDS envelopes: %s genes\n",
            format(nrow(gene_cds), big.mark = ",")))

# Match rate with rMATS genes
rmats_genes <- unique(all_events$gene_id_unversioned)
matched <- sum(rmats_genes %in% gene_cds$gene_id)
cat(sprintf("  rMATS gene match rate: %d / %d (%.1f%%)\n",
            matched, length(rmats_genes),
            100 * matched / length(rmats_genes)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Frame Disruption Classification
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Section 3: Classifying frame disruption ─────────────────────────\n\n")

# Join events with CDS envelope
all_events <- all_events %>%
  left_join(gene_cds, by = c("gene_id_unversioned" = "gene_id",
                              "event_chr" = "seqnames"))

# Determine CDS overlap: event coordinates within gene CDS envelope
all_events <- all_events %>%
  mutate(
    overlaps_cds = !is.na(cds_min) &
                   event_start <= cds_max &
                   event_end >= cds_min
  )

# Frame disruption classification per event type
# - SE: exon_target width mod 3 != 0 → frame-disrupting
# - A5SS/A3SS: |long_width - short_width| mod 3 != 0
# - MXE: |exon1_width - exon2_width| mod 3 != 0
# - RI: all CDS-overlapping RI treated as PTC-inducing (introns contain stops)
all_events <- all_events %>%
  mutate(
    frame_disrupting = case_when(
      event_type == "SE"  ~ (event_width %% 3) != 0,
      event_type %in% c("A5SS", "A3SS") ~ (width_diff %% 3) != 0,
      event_type == "MXE" ~ (width_diff %% 3) != 0,
      event_type == "RI"  ~ overlaps_cds,  # all CDS RI are PTC-inducing
      TRUE ~ NA
    )
  )

# Three-way classification
all_events <- all_events %>%
  mutate(
    frame_class = case_when(
      !overlaps_cds                      ~ "non_cds",
      event_type == "RI" & overlaps_cds  ~ "cds_ptc_inducing",
      overlaps_cds & frame_disrupting    ~ "cds_frame_disrupting",
      overlaps_cds & !frame_disrupting   ~ "cds_frame_preserving",
      TRUE                               ~ "unclassified"
    )
  )

# Flag significance
all_events <- all_events %>%
  mutate(
    is_significant = !is.na(FDR) & FDR < FDR_THRESH &
                     !is.na(IncLevelDifference) & abs(IncLevelDifference) >= DPSI_THRESH
  )

# Summary
class_summary <- all_events %>%
  count(event_type, frame_class) %>%
  pivot_wider(names_from = frame_class, values_from = n, values_fill = 0)
cat("  Classification summary:\n")
print(as.data.frame(class_summary), row.names = FALSE)

sig_summary <- all_events %>%
  filter(overlaps_cds | frame_class == "non_cds") %>%
  group_by(event_type, frame_class) %>%
  summarise(total = n(), significant = sum(is_significant),
            pct_sig = round(100 * significant / total, 2), .groups = "drop")
cat("\n  Significance rates by class:\n")
print(as.data.frame(sig_summary), row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Enrichment Tests
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Section 4: Enrichment analysis ──────────────────────────────────\n\n")

# Helper: run Fisher test on 2x2 table
run_fisher <- function(df, group_col, positive_val, negative_val) {
  tbl <- df %>%
    filter(!!sym(group_col) %in% c(positive_val, negative_val)) %>%
    mutate(is_target = !!sym(group_col) == positive_val) %>%
    count(is_target, is_significant) %>%
    pivot_wider(names_from = is_significant, values_from = n, values_fill = 0)

  if (nrow(tbl) < 2 || !all(c(TRUE, FALSE) %in% tbl$is_target)) {
    return(tibble(OR = NA_real_, OR_lower = NA_real_, OR_upper = NA_real_,
                  p_value = NA_real_, n_target = NA_integer_,
                  n_control = NA_integer_,
                  sig_target = NA_integer_, sig_control = NA_integer_))
  }

  mat <- matrix(0, nrow = 2, ncol = 2)
  target_row  <- tbl %>% filter(is_target)
  control_row <- tbl %>% filter(!is_target)

  mat[1, 1] <- if ("TRUE" %in% names(target_row))  target_row[["TRUE"]]  else 0
  mat[1, 2] <- if ("FALSE" %in% names(target_row))  target_row[["FALSE"]] else 0
  mat[2, 1] <- if ("TRUE" %in% names(control_row))  control_row[["TRUE"]]  else 0
  mat[2, 2] <- if ("FALSE" %in% names(control_row))  control_row[["FALSE"]] else 0

  ft <- fisher.test(mat)
  tibble(
    OR = ft$estimate,
    OR_lower = ft$conf.int[1],
    OR_upper = ft$conf.int[2],
    p_value = ft$p.value,
    n_target = mat[1, 1] + mat[1, 2],
    n_control = mat[2, 1] + mat[2, 2],
    sig_target = mat[1, 1],
    sig_control = mat[2, 1]
  )
}

# --- 4a: Frame-disrupting vs frame-preserving CDS events, per cell type ---
# (Excludes RI — uses separate classification; and non-CDS events)
enrichment_by_ct <- all_events %>%
  filter(event_type != "RI",
         frame_class %in% c("cds_frame_disrupting", "cds_frame_preserving")) %>%
  group_by(cell_type) %>%
  group_modify(~run_fisher(.x, "frame_class",
                            "cds_frame_disrupting", "cds_frame_preserving")) %>%
  ungroup() %>%
  mutate(comparison = "frame_disrupting_vs_preserving")

# RI-specific: CDS RI vs non-CDS RI
enrichment_ri <- all_events %>%
  filter(event_type == "RI") %>%
  mutate(ri_class = if_else(overlaps_cds, "cds_ri", "non_cds_ri")) %>%
  group_by(cell_type) %>%
  group_modify(~run_fisher(.x, "ri_class", "cds_ri", "non_cds_ri")) %>%
  ungroup() %>%
  mutate(comparison = "ri_cds_vs_non_cds")

enrichment_combined <- bind_rows(enrichment_by_ct, enrichment_ri)

cat("  Enrichment: frame-disrupting CDS vs frame-preserving CDS (excl. RI)\n")
enrichment_by_ct %>%
  mutate(across(c(OR, OR_lower, OR_upper), ~round(., 3)),
         p_value = signif(p_value, 3)) %>%
  print(n = Inf)

cat("\n  Enrichment: CDS RI vs non-CDS RI\n")
enrichment_ri %>%
  mutate(across(c(OR, OR_lower, OR_upper), ~round(., 3)),
         p_value = signif(p_value, 3)) %>%
  print(n = Inf)

# --- 4b: Stratified by event type ---
enrichment_by_etype <- all_events %>%
  filter(event_type != "RI",
         frame_class %in% c("cds_frame_disrupting", "cds_frame_preserving")) %>%
  group_by(cell_type, event_type) %>%
  group_modify(~run_fisher(.x, "frame_class",
                            "cds_frame_disrupting", "cds_frame_preserving")) %>%
  ungroup()

cat("\n  Enrichment by event type (showing DD only for brevity):\n")
enrichment_by_etype %>%
  filter(cell_type == "DD") %>%
  mutate(across(c(OR, OR_lower, OR_upper), ~round(., 3)),
         p_value = signif(p_value, 3)) %>%
  print(n = Inf)

# --- 4c: Meta-analysis across cell types (random effects) ---
run_meta <- function(enrich_df, label) {
  df <- enrich_df %>% filter(!is.na(OR), OR > 0)
  if (nrow(df) < 2) {
    cat(sprintf("  %s: <2 cell types with valid OR, skipping meta-analysis\n", label))
    return(NULL)
  }
  # Log OR and SE for meta-analysis
  df <- df %>%
    mutate(
      log_or = log(OR),
      se_log_or = (log(OR_upper) - log(OR_lower)) / (2 * 1.96)
    )
  fit <- tryCatch(
    rma(yi = df$log_or, sei = df$se_log_or, method = "REML"),
    error = function(e) {
      cat(sprintf("  %s meta-analysis error: %s\n", label, e$message))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)

  cat(sprintf("\n  Meta-analysis: %s\n", label))
  cat(sprintf("    Pooled OR: %.3f (%.3f-%.3f), p=%.2e\n",
              exp(fit$beta), exp(fit$ci.lb), exp(fit$ci.ub), fit$pval))
  cat(sprintf("    I² = %.1f%%, tau² = %.4f, Q p=%.3f\n",
              fit$I2, fit$tau2, fit$QEp))
  fit
}

meta_frame <- run_meta(enrichment_by_ct, "Frame-disrupting vs preserving")
meta_ri    <- run_meta(enrichment_ri, "CDS RI vs non-CDS RI")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5: Direction Analysis
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Section 5: Direction analysis ───────────────────────────────────\n\n")

# For SE, A5SS, A3SS, MXE: we cannot predict direction (inclusion vs skipping
# could create the PTC), so test for enrichment at BOTH tails (|dPSI|).
# For RI: intron retention almost always introduces stop codons, so we expect
# the retention form to be rescued by Smg1i → positive dPSI.

direction_tests <- list()

# --- 5a: SE direction among significant frame-disrupting events ---
# Two-tailed: test if dPSI is biased in either direction
for (ct in unique(all_events$cell_type)) {
  for (etype in c("SE", "A5SS", "A3SS", "MXE")) {
    sig_fd <- all_events %>%
      filter(cell_type == ct, event_type == etype,
             frame_class == "cds_frame_disrupting", is_significant)
    n_events <- nrow(sig_fd)
    if (n_events < 5) next

    n_pos <- sum(sig_fd$IncLevelDifference > 0)
    n_neg <- sum(sig_fd$IncLevelDifference < 0)
    bt <- binom.test(n_pos, n_pos + n_neg, p = 0.5)

    direction_tests[[paste(ct, etype, "frame_disrupting")]] <- tibble(
      cell_type = ct, event_type = etype,
      frame_class = "cds_frame_disrupting",
      n_events = n_events, n_positive = n_pos, n_negative = n_neg,
      frac_positive = round(n_pos / (n_pos + n_neg), 3),
      binom_p = bt$p.value,
      test_type = "two_tailed"
    )
  }

  # Frame-preserving control
  for (etype in c("SE", "A5SS", "A3SS", "MXE")) {
    sig_fp <- all_events %>%
      filter(cell_type == ct, event_type == etype,
             frame_class == "cds_frame_preserving", is_significant)
    n_events <- nrow(sig_fp)
    if (n_events < 5) next

    n_pos <- sum(sig_fp$IncLevelDifference > 0)
    n_neg <- sum(sig_fp$IncLevelDifference < 0)
    bt <- binom.test(n_pos, n_pos + n_neg, p = 0.5)

    direction_tests[[paste(ct, etype, "frame_preserving")]] <- tibble(
      cell_type = ct, event_type = etype,
      frame_class = "cds_frame_preserving",
      n_events = n_events, n_positive = n_pos, n_negative = n_neg,
      frac_positive = round(n_pos / (n_pos + n_neg), 3),
      binom_p = bt$p.value,
      test_type = "two_tailed"
    )
  }

  # --- 5b: RI direction (one-tailed: expected positive dPSI) ---
  sig_ri <- all_events %>%
    filter(cell_type == ct, event_type == "RI",
           frame_class == "cds_ptc_inducing", is_significant)
  n_ri <- nrow(sig_ri)
  if (n_ri >= 5) {
    n_pos_ri <- sum(sig_ri$IncLevelDifference > 0)
    n_neg_ri <- sum(sig_ri$IncLevelDifference < 0)
    bt_ri <- binom.test(n_pos_ri, n_pos_ri + n_neg_ri,
                         p = 0.5, alternative = "greater")
    direction_tests[[paste(ct, "RI_cds")]] <- tibble(
      cell_type = ct, event_type = "RI",
      frame_class = "cds_ptc_inducing",
      n_events = n_ri, n_positive = n_pos_ri, n_negative = n_neg_ri,
      frac_positive = round(n_pos_ri / (n_pos_ri + n_neg_ri), 3),
      binom_p = bt_ri$p.value,
      test_type = "one_tailed_positive"
    )
  }
}

direction_df <- bind_rows(direction_tests)

if (nrow(direction_df) > 0) {
  cat("  Direction tests (showing cell types with results):\n")
  direction_df %>%
    mutate(binom_p = signif(binom_p, 3)) %>%
    arrange(event_type, frame_class, cell_type) %>%
    print(n = 40)
} else {
  cat("  No direction tests had sufficient events (n >= 5)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6: Magnitude Analysis
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Section 6: Magnitude analysis ───────────────────────────────────\n\n")

# Compare |dPSI| between frame-disrupting and frame-preserving CDS events
# (not limited to significant — tests all events with CDS overlap)
magnitude_tests <- list()

for (ct in unique(all_events$cell_type)) {
  # All non-RI event types together
  fd <- all_events %>%
    filter(cell_type == ct, event_type != "RI",
           frame_class == "cds_frame_disrupting",
           !is.na(IncLevelDifference))
  fp <- all_events %>%
    filter(cell_type == ct, event_type != "RI",
           frame_class == "cds_frame_preserving",
           !is.na(IncLevelDifference))

  if (nrow(fd) >= 10 && nrow(fp) >= 10) {
    wt <- wilcox.test(abs(fd$IncLevelDifference),
                       abs(fp$IncLevelDifference))
    magnitude_tests[[paste(ct, "all_non_RI")]] <- tibble(
      cell_type = ct, event_type = "all_non_RI",
      n_disrupting = nrow(fd), n_preserving = nrow(fp),
      median_abs_dpsi_disrupting = round(median(abs(fd$IncLevelDifference)), 4),
      median_abs_dpsi_preserving = round(median(abs(fp$IncLevelDifference)), 4),
      wilcox_p = wt$p.value
    )
  }

  # Per event type
  for (etype in c("SE", "A5SS", "A3SS", "MXE")) {
    fd_e <- fd %>% filter(event_type == etype)
    fp_e <- fp %>% filter(event_type == etype)
    if (nrow(fd_e) >= 10 && nrow(fp_e) >= 10) {
      wt_e <- wilcox.test(abs(fd_e$IncLevelDifference),
                           abs(fp_e$IncLevelDifference))
      magnitude_tests[[paste(ct, etype)]] <- tibble(
        cell_type = ct, event_type = etype,
        n_disrupting = nrow(fd_e), n_preserving = nrow(fp_e),
        median_abs_dpsi_disrupting = round(median(abs(fd_e$IncLevelDifference)), 4),
        median_abs_dpsi_preserving = round(median(abs(fp_e$IncLevelDifference)), 4),
        wilcox_p = wt_e$p.value
      )
    }
  }

  # RI: CDS vs non-CDS
  ri_cds <- all_events %>%
    filter(cell_type == ct, event_type == "RI", overlaps_cds,
           !is.na(IncLevelDifference))
  ri_non <- all_events %>%
    filter(cell_type == ct, event_type == "RI", !overlaps_cds,
           !is.na(IncLevelDifference))
  if (nrow(ri_cds) >= 10 && nrow(ri_non) >= 10) {
    wt_ri <- wilcox.test(abs(ri_cds$IncLevelDifference),
                          abs(ri_non$IncLevelDifference))
    magnitude_tests[[paste(ct, "RI")]] <- tibble(
      cell_type = ct, event_type = "RI",
      n_disrupting = nrow(ri_cds), n_preserving = nrow(ri_non),
      median_abs_dpsi_disrupting = round(median(abs(ri_cds$IncLevelDifference)), 4),
      median_abs_dpsi_preserving = round(median(abs(ri_non$IncLevelDifference)), 4),
      wilcox_p = wt_ri$p.value
    )
  }
}

magnitude_df <- bind_rows(magnitude_tests)
cat("  Magnitude tests (|dPSI| frame-disrupting vs frame-preserving):\n")
magnitude_df %>%
  mutate(wilcox_p = signif(wilcox_p, 3)) %>%
  print(n = 40)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7: Summary and Output
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Section 7: Saving outputs ───────────────────────────────────────\n\n")

# --- Summary stats table ---
summary_stats <- all_events %>%
  group_by(cell_type, event_type, frame_class) %>%
  summarise(
    total = n(),
    significant = sum(is_significant),
    pct_sig = round(100 * significant / total, 2),
    median_abs_dpsi = round(median(abs(IncLevelDifference), na.rm = TRUE), 4),
    .groups = "drop"
  )

# --- Save all outputs ---
# Event classification (large file — write key columns only)
event_out <- all_events %>%
  select(cell_type, event_type, ID, GeneID, gene_id_unversioned, geneSymbol,
         event_chr, event_start, event_end,
         FDR, IncLevelDifference, is_significant,
         overlaps_cds, frame_disrupting, frame_class)

write_tsv(event_out, file.path(out_dir, "rmats_event_classification.tsv"))
cat(sprintf("  Wrote rmats_event_classification.tsv (%s events)\n",
            format(nrow(event_out), big.mark = ",")))

write_tsv(enrichment_combined, file.path(out_dir, "rmats_enrichment_by_celltype.tsv"))
cat("  Wrote rmats_enrichment_by_celltype.tsv\n")

write_tsv(enrichment_by_etype, file.path(out_dir, "rmats_enrichment_by_event_type.tsv"))
cat("  Wrote rmats_enrichment_by_event_type.tsv\n")

if (nrow(direction_df) > 0) {
  write_tsv(direction_df, file.path(out_dir, "rmats_direction_tests.tsv"))
  cat("  Wrote rmats_direction_tests.tsv\n")
}

write_tsv(magnitude_df, file.path(out_dir, "rmats_magnitude_tests.tsv"))
cat("  Wrote rmats_magnitude_tests.tsv\n")

write_tsv(summary_stats, file.path(out_dir, "rmats_summary_stats.tsv"))
cat("  Wrote rmats_summary_stats.tsv\n")

# --- Figures ---
cat("\n  Generating figures...\n")

# Figure 1: dPSI distributions by frame class
fig1_data <- all_events %>%
  filter(event_type != "RI",
         frame_class %in% c("cds_frame_disrupting", "cds_frame_preserving"),
         !is.na(IncLevelDifference))

p1 <- ggplot(fig1_data, aes(x = IncLevelDifference,
                              fill = frame_class)) +
  geom_density(alpha = 0.5) +
  facet_grid(event_type ~ cell_type, scales = "free_y") +
  scale_fill_manual(values = c("cds_frame_disrupting" = "#E74C3C",
                                "cds_frame_preserving" = "#3498DB"),
                    labels = c("Frame-disrupting", "Frame-preserving")) +
  geom_vline(xintercept = 0, lty = 2, color = "grey40") +
  labs(x = "IncLevelDifference (Smg1i - DMSO)",
       y = "Density",
       fill = "Frame class",
       title = "dPSI distributions: frame-disrupting vs frame-preserving CDS events") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8))

ggsave(file.path(fig_dir, "rmats_dpsi_by_frame_class.pdf"),
       p1, width = 14, height = 10)
cat("  Wrote rmats_dpsi_by_frame_class.pdf\n")

# Figure 2: Forest plot of enrichment ORs
# Filter out cell types with OR=0 or Inf (no significant events)
forest_data <- enrichment_by_ct %>%
  filter(!is.na(OR), OR > 0, is.finite(OR_upper))

if (nrow(forest_data) > 0) {
  # Add meta-analysis result if available
  if (!is.null(meta_frame)) {
    meta_row <- tibble(
      cell_type = "Meta-analysis",
      OR = exp(as.numeric(meta_frame$beta)),
      OR_lower = exp(meta_frame$ci.lb),
      OR_upper = exp(meta_frame$ci.ub),
      is_meta = TRUE
    )
    forest_data2 <- bind_rows(
      forest_data %>% select(cell_type, OR, OR_lower, OR_upper) %>%
        mutate(is_meta = FALSE),
      meta_row
    ) %>%
      mutate(cell_type = factor(cell_type,
                                levels = c("Meta-analysis",
                                           rev(sort(setdiff(unique(cell_type),
                                                            "Meta-analysis"))))))
  } else {
    forest_data2 <- forest_data %>%
      select(cell_type, OR, OR_lower, OR_upper) %>%
      mutate(is_meta = FALSE,
             cell_type = factor(cell_type,
                                levels = rev(sort(unique(cell_type)))))
  }

  p2 <- ggplot(forest_data2, aes(x = OR, y = cell_type)) +
    geom_point(aes(shape = is_meta), size = 3) +
    geom_errorbar(aes(xmin = OR_lower, xmax = OR_upper),
                  width = 0.2, orientation = "y") +
    geom_vline(xintercept = 1, lty = 2, color = "grey50") +
    scale_x_log10() +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    labs(x = "Odds Ratio (log scale)",
         y = "",
         title = "Enrichment: frame-disrupting CDS events among significant rMATS events",
         subtitle = sprintf("FDR < %.2f, |dPSI| > %.2f | Random-effects meta (excl. RI, cell types w/o sig events)",
                            FDR_THRESH, DPSI_THRESH)) +
    theme_bw(base_size = 12)

  ggsave(file.path(fig_dir, "rmats_enrichment_forest.pdf"),
         p2, width = 9, height = 5)
  cat("  Wrote rmats_enrichment_forest.pdf\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# CONSOLE SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("==================================================================\n")
cat("   SUMMARY\n")
cat("==================================================================\n\n")

total_events <- nrow(all_events)
cds_events <- sum(all_events$overlaps_cds, na.rm = TRUE)
cat(sprintf("  Total rMATS events: %s\n", format(total_events, big.mark = ",")))
cat(sprintf("  CDS-overlapping: %s (%.1f%%)\n",
            format(cds_events, big.mark = ","),
            100 * cds_events / total_events))

# Non-RI frame disruption
non_ri_cds <- all_events %>%
  filter(event_type != "RI",
         frame_class %in% c("cds_frame_disrupting", "cds_frame_preserving"))
cat(sprintf("  CDS non-RI events: %s (disrupting: %s, preserving: %s)\n",
            format(nrow(non_ri_cds), big.mark = ","),
            format(sum(non_ri_cds$frame_class == "cds_frame_disrupting"), big.mark = ","),
            format(sum(non_ri_cds$frame_class == "cds_frame_preserving"), big.mark = ",")))

# Key enrichment results
cat("\n  Key enrichment results (OR > 1 = frame-disrupting enriched):\n")
enrichment_by_ct %>%
  mutate(summary = sprintf("    %s: OR=%.3f (%.3f-%.3f), p=%.2e [%d/%d vs %d/%d]",
                           cell_type, OR, OR_lower, OR_upper, p_value,
                           sig_target, n_target, sig_control, n_control)) %>%
  pull(summary) %>%
  cat(sep = "\n")

# Meta-analysis summary
if (!is.null(meta_frame)) {
  cat(sprintf("\n  Meta-analysis (random effects): OR=%.3f (%.3f-%.3f), p=%.2e, I²=%.1f%%\n",
              exp(as.numeric(meta_frame$beta)), exp(meta_frame$ci.lb),
              exp(meta_frame$ci.ub), meta_frame$pval, meta_frame$I2))
}

# RI summary
cat("\n  RI enrichment (CDS vs non-CDS):\n")
enrichment_ri %>%
  mutate(summary = sprintf("    %s: OR=%.3f (%.3f-%.3f), p=%.2e",
                           cell_type, OR, OR_lower, OR_upper, p_value)) %>%
  pull(summary) %>%
  cat(sep = "\n")

# Direction summary
if (nrow(direction_df) > 0) {
  cat("\n  Direction bias (significant frame-disrupting events):\n")
  direction_df %>%
    filter(frame_class %in% c("cds_frame_disrupting", "cds_ptc_inducing")) %>%
    mutate(summary = sprintf("    %s %s: %.0f%% positive (n=%d), p=%.3f",
                             cell_type, event_type,
                             100 * frac_positive, n_events, binom_p)) %>%
    pull(summary) %>%
    cat(sep = "\n")
}

cat("\n\n  Done.\n")
