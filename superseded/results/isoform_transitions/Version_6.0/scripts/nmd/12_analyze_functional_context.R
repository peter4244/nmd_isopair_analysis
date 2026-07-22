#!/usr/bin/env Rscript
#
# Script: 12_analyze_functional_context.R
# Version: 6.0
# Purpose: Analyze regional distribution (UTR vs CDS) and ORF impact of splicing events
#
# Input:
#   - splicing_choice_profiles.rds (from Script 08)
#   - data/isoform_union_exons_annotated_filtered.rds (region types per isoform x union exon)
#   - data/union_exons_filtered.rds (union exon coordinates for bp calculations)
#
# Output:
#   - {output_dir}/event_regions.rds (per-event region mappings for cross-comparison bootstrap)
#   - {output_dir}/regional_distribution_results.tsv
#   - {output_dir}/orf_impact_summary.tsv
#   - {output_dir}/figures/regional_enrichment.pdf
#
# Flags:
#   --input path.rds       Input profiles (default: data/splicing_choice_profiles.rds)
#   --output-dir dir/      Output directory (default: results)
#

library(tidyverse)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 12: Analyze Functional Context                         ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Paths
base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"
setwd(base_dir)

# ═══════════════════════════════════════════════════════════════════
# Parse CLI arguments
# ═══════════════════════════════════════════════════════════════════

args <- commandArgs(trailingOnly = TRUE)

input_path <- "data/splicing_choice_profiles.rds"
if ("--input" %in% args) {
  idx <- which(args == "--input")
  input_path <- args[idx + 1]
}

output_dir <- "results"
if ("--output-dir" %in% args) {
  idx <- which(args == "--output-dir")
  output_dir <- args[idx + 1]
}

# Parse --source (default: oarfish)
source_type <- "oarfish"
if ("--source" %in% args) {
  src_idx <- which(args == "--source")
  if (src_idx < length(args)) {
    source_type <- args[src_idx + 1]
    if (!source_type %in% c("oarfish", "isocall")) {
      stop("--source must be 'oarfish' or 'isocall'")
    }
  }
}

# Set data directory based on source
data_dir <- if (source_type == "isocall") "data/isocall" else "data"

cat(sprintf("  Input:      %s\n", input_path))
cat(sprintf("  Output dir: %s\n", output_dir))
cat(sprintf("  Source:     %s\n\n", source_type))

# ═══════════════════════════════════════════════════════════════════
# SECTION 1: Load Data
# ═══════════════════════════════════════════════════════════════════

cat("Loading data...\n")

if (!file.exists(input_path)) {
  stop(sprintf("ERROR: Input file '%s' not found.", input_path))
}

profiles <- readRDS(input_path)
cat(sprintf("  Loaded %d profiles\n", nrow(profiles)))

# Load annotated union exons (region type is per isoform x union_exon pair)
annotated_ue <- readRDS(file.path(data_dir, "isoform_union_exons_annotated_filtered.rds"))
cat(sprintf("  Loaded %d annotated union exon mappings\n", nrow(annotated_ue)))

# Load union exons (for bp calculations)
union_exons <- readRDS(file.path(data_dir, "union_exons_filtered.rds"))
cat(sprintf("  Loaded %d union exons\n", nrow(union_exons)))

# ═══════════════════════════════════════════════════════════════════
# SECTION 2: Extract Events and Map to Regions
# ═══════════════════════════════════════════════════════════════════

cat("\nExtracting events and mapping to region types...\n")

# Extract all events
# Note: detailed_events contains gene_id internally, so don't select it from outer
all_events <- profiles %>%
  select(dominant_isoform_id, non_dominant_isoform_id, detailed_events) %>%
  filter(map_int(detailed_events, nrow) > 0) %>%
  unnest(detailed_events) %>%
  mutate(
    genomic_start = pmin(five_prime, three_prime),
    genomic_end = pmax(five_prime, three_prime),
    event_row_id = row_number()
  )

cat(sprintf("  Total events: %d\n", nrow(all_events)))

# For each event, find overlapping union exons from the DOMINANT isoform's annotations
# Filter annotated_ue to only dominant isoforms in our profiles
dom_isoform_ids <- unique(all_events$dominant_isoform_id)

dom_annotations <- annotated_ue %>%
  filter(isoform_id %in% dom_isoform_ids) %>%
  select(gene_id, isoform_id, union_exon_id, union_exon_start, union_exon_end, region_type)

cat(sprintf("  Dominant isoform annotations: %d rows\n", nrow(dom_annotations)))

# Join events to overlapping union exons by coordinate overlap within same gene
# An event overlaps a union exon if: event_genomic_start <= ue_end AND event_genomic_end >= ue_start
# Report ALL region types touched (per plan decision)
event_regions <- all_events %>%
  select(event_row_id, gene_id, dominant_isoform_id, non_dominant_isoform_id,
         event_type, direction, genomic_start, genomic_end) %>%
  inner_join(
    dom_annotations %>% rename(dominant_isoform_id = isoform_id),
    by = c("gene_id", "dominant_isoform_id"),
    relationship = "many-to-many"
  ) %>%
  filter(genomic_start <= union_exon_end, genomic_end >= union_exon_start) %>%
  distinct(event_row_id, gene_id, dominant_isoform_id, non_dominant_isoform_id,
           event_type, direction, region_type)

cat(sprintf("  Event-region mappings: %d (events can map to multiple regions)\n", nrow(event_regions)))

# Events that didn't overlap any union exon (intronic regions)
events_no_region <- all_events %>%
  filter(!event_row_id %in% event_regions$event_row_id)

cat(sprintf("  Events with no overlapping union exon (intronic): %d\n", nrow(events_no_region)))

# Add intronic events as "intronic" region
if (nrow(events_no_region) > 0) {
  intronic_regions <- events_no_region %>%
    select(event_row_id, gene_id, dominant_isoform_id, non_dominant_isoform_id,
           event_type, direction) %>%
    mutate(region_type = "intronic")
  event_regions <- bind_rows(event_regions, intronic_regions)
}

# ═══════════════════════════════════════════════════════════════════
# SECTION 3: Regional Distribution
# ═══════════════════════════════════════════════════════════════════

cat("\nAnalyzing regional distribution of events...\n")

# Compute expected proportions from dominant isoform bp per region
dom_bp_by_region <- annotated_ue %>%
  filter(isoform_id %in% dom_isoform_ids) %>%
  # Use one row per (isoform, union_exon) — avoid double counting
  distinct(isoform_id, union_exon_id, .keep_all = TRUE) %>%
  mutate(bp = union_exon_end - union_exon_start + 1) %>%
  group_by(region_type) %>%
  summarise(total_bp = sum(bp), .groups = "drop") %>%
  mutate(expected_prop = total_bp / sum(total_bp))

cat("\nExpected proportions (dominant isoform bp):\n")
print(dom_bp_by_region)

# Count events by region type for each event type
regional_counts <- event_regions %>%
  count(event_type, region_type, name = "observed_count")

# Total events per event type (for proportions)
event_totals <- event_regions %>%
  count(event_type, name = "total_events")

regional_counts <- regional_counts %>%
  left_join(event_totals, by = "event_type") %>%
  mutate(observed_prop = observed_count / total_events) %>%
  left_join(dom_bp_by_region %>% select(region_type, expected_prop), by = "region_type") %>%
  mutate(
    expected_prop = replace_na(expected_prop, 0),
    enrichment_ratio = ifelse(expected_prop > 0, observed_prop / expected_prop, NA_real_)
  )

# Binomial test for each event_type x region_type
regional_distribution <- regional_counts %>%
  rowwise() %>%
  mutate(
    binom_test = list(
      tryCatch(
        binom.test(observed_count, total_events, p = max(expected_prop, 1e-6)),
        error = function(e) NULL
      )
    ),
    binom_p = ifelse(!is.null(binom_test), binom_test$p.value, NA_real_),
    binom_ci_lower = ifelse(!is.null(binom_test), binom_test$conf.int[1], NA_real_),
    binom_ci_upper = ifelse(!is.null(binom_test), binom_test$conf.int[2], NA_real_)
  ) %>%
  ungroup() %>%
  select(-binom_test) %>%
  mutate(
    binom_adj_p = p.adjust(binom_p, method = "fdr"),
    significant = binom_adj_p < 0.05,
    direction_enrich = case_when(
      !significant ~ "NS",
      enrichment_ratio > 1 ~ "Enriched",
      enrichment_ratio < 1 ~ "Depleted",
      TRUE ~ "NS"
    )
  )

cat("\nRegional distribution results (significant only):\n")
sig_results <- regional_distribution %>%
  filter(significant) %>%
  select(event_type, region_type, observed_count, enrichment_ratio, binom_adj_p, direction_enrich) %>%
  arrange(event_type, desc(enrichment_ratio))
print(sig_results, n = 30)

# ═══════════════════════════════════════════════════════════════════
# SECTION 4: ORF Boundary Susceptibility
# ═══════════════════════════════════════════════════════════════════

cat("\nTesting ORF boundary exon susceptibility to A5SS/A3SS...\n")

# Identify union exons with ORF boundary region types (for dominant isoforms)
orf_boundary_ues <- annotated_ue %>%
  filter(isoform_id %in% dom_isoform_ids) %>%
  distinct(isoform_id, union_exon_id, .keep_all = TRUE) %>%
  mutate(
    is_orf_boundary = region_type %in% c("contains_orf_start", "contains_orf_stop",
                                          "contains_orf_start_stop"),
    is_regular_cds = region_type == "CDS"
  )

n_orf_boundary_ues <- sum(orf_boundary_ues$is_orf_boundary)
n_regular_cds_ues <- sum(orf_boundary_ues$is_regular_cds)

cat(sprintf("  ORF boundary union exons: %d\n", n_orf_boundary_ues))
cat(sprintf("  Regular CDS union exons: %d\n", n_regular_cds_ues))

# Count union exons with at least one A5SS/A3SS event overlapping them
# First get all A5SS/A3SS events with their genomic coordinates
a5a3_events <- all_events %>%
  filter(event_type %in% c("A5SS", "A3SS"))

# Find union exons hit by A5SS/A3SS events (via coordinate overlap)
if (nrow(a5a3_events) > 0) {
  a5a3_event_ue_hits <- a5a3_events %>%
    select(gene_id, dominant_isoform_id, genomic_start, genomic_end) %>%
    inner_join(
      dom_annotations %>% rename(dominant_isoform_id = isoform_id),
      by = c("gene_id", "dominant_isoform_id"),
      relationship = "many-to-many"
    ) %>%
    filter(genomic_start <= union_exon_end, genomic_end >= union_exon_start) %>%
    distinct(dominant_isoform_id, union_exon_id, region_type)

  # Count UEs with A5SS/A3SS hits by region class
  ues_hit_orf_boundary <- a5a3_event_ue_hits %>%
    filter(region_type %in% c("contains_orf_start", "contains_orf_stop",
                               "contains_orf_start_stop")) %>%
    nrow()
  ues_hit_regular_cds <- a5a3_event_ue_hits %>%
    filter(region_type == "CDS") %>%
    nrow()
} else {
  ues_hit_orf_boundary <- 0L
  ues_hit_regular_cds <- 0L
}

if (n_orf_boundary_ues > 0 && n_regular_cds_ues > 0) {
  # Ensure counts don't exceed totals (cap at total)
  ues_hit_orf_boundary <- min(ues_hit_orf_boundary, n_orf_boundary_ues)
  ues_hit_regular_cds <- min(ues_hit_regular_cds, n_regular_cds_ues)

  # Fisher's exact test: are ORF boundary exons more likely to have A5SS/A3SS events?
  contingency <- matrix(
    c(ues_hit_orf_boundary, n_orf_boundary_ues - ues_hit_orf_boundary,
      ues_hit_regular_cds, n_regular_cds_ues - ues_hit_regular_cds),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("ORF_boundary", "Regular_CDS"), c("Has_A5SS_A3SS", "No_event"))
  )

  fisher_result <- fisher.test(contingency)

  orf_susceptibility <- tibble(
    exon_class = c("ORF_boundary", "Regular_CDS"),
    n_exons = c(n_orf_boundary_ues, n_regular_cds_ues),
    n_with_a5ss_a3ss = c(ues_hit_orf_boundary, ues_hit_regular_cds),
    rate = n_with_a5ss_a3ss / n_exons,
    odds_ratio = fisher_result$estimate,
    or_ci_lower = fisher_result$conf.int[1],
    or_ci_upper = fisher_result$conf.int[2],
    fisher_p = fisher_result$p.value
  )

  cat(sprintf("  ORF boundary exons with A5SS/A3SS: %d/%d (%.4f)\n",
              ues_hit_orf_boundary, n_orf_boundary_ues,
              ues_hit_orf_boundary / n_orf_boundary_ues))
  cat(sprintf("  Regular CDS exons with A5SS/A3SS: %d/%d (%.4f)\n",
              ues_hit_regular_cds, n_regular_cds_ues,
              ues_hit_regular_cds / n_regular_cds_ues))
  cat(sprintf("  Fisher's OR: %.2f (%.2f - %.2f), p = %.2e\n",
              fisher_result$estimate, fisher_result$conf.int[1],
              fisher_result$conf.int[2], fisher_result$p.value))
} else {
  cat("  Insufficient data for ORF boundary susceptibility test\n")
  orf_susceptibility <- tibble(
    exon_class = character(0),
    n_exons = integer(0),
    n_with_a5ss_a3ss = integer(0),
    rate = numeric(0),
    odds_ratio = numeric(0),
    or_ci_lower = numeric(0),
    or_ci_upper = numeric(0),
    fisher_p = numeric(0)
  )
}

# ═══════════════════════════════════════════════════════════════════
# SECTION 5: ORF Impact Summary
# ═══════════════════════════════════════════════════════════════════

cat("\nSummarizing ORF impact...\n")

# For Alt_TSS events: what fraction affect exons containing ORF start?
alt_tss_events <- event_regions %>%
  filter(event_type == "Alt_TSS")

alt_tss_total <- n_distinct(alt_tss_events$event_row_id)
alt_tss_orf_start <- alt_tss_events %>%
  filter(region_type %in% c("contains_orf_start", "contains_orf_start_stop")) %>%
  pull(event_row_id) %>% n_distinct()

# For Alt_TES events: what fraction affect exons containing ORF stop?
alt_tes_events <- event_regions %>%
  filter(event_type == "Alt_TES")

alt_tes_total <- n_distinct(alt_tes_events$event_row_id)
alt_tes_orf_stop <- alt_tes_events %>%
  filter(region_type %in% c("contains_orf_stop", "contains_orf_start_stop")) %>%
  pull(event_row_id) %>% n_distinct()

orf_impact <- tibble(
  event_type = c("Alt_TSS", "Alt_TES"),
  orf_region = c("contains_orf_start", "contains_orf_stop"),
  n_events = c(alt_tss_total, alt_tes_total),
  n_affecting_orf = c(alt_tss_orf_start, alt_tes_orf_stop),
  pct_affecting_orf = ifelse(c(alt_tss_total, alt_tes_total) > 0,
                              100 * c(alt_tss_orf_start, alt_tes_orf_stop) /
                                c(alt_tss_total, alt_tes_total),
                              NA_real_)
)

cat("\nORF Impact Summary:\n")
print(orf_impact)

# ═══════════════════════════════════════════════════════════════════
# SECTION 6: Visualization
# ═══════════════════════════════════════════════════════════════════

cat("\nCreating regional enrichment plots...\n")

fig_dir <- file.path(output_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(fig_dir, "regional_enrichment.pdf"), width = 14, height = 10)

# Plot 1: Enrichment ratios by event type and region (grouped barplot)
plot_data <- regional_distribution %>%
  filter(!is.na(enrichment_ratio), total_events >= 10) %>%
  mutate(
    log2_enrich = log2(enrichment_ratio),
    sig_label = ifelse(significant, "*", "")
  )

if (nrow(plot_data) > 0) {
  p1 <- ggplot(plot_data, aes(x = region_type, y = log2_enrich, fill = event_type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_text(aes(label = sig_label, y = log2_enrich + sign(log2_enrich) * 0.1),
              position = position_dodge(width = 0.8), size = 5) +
    labs(
      title = "Regional Enrichment of Splicing Events",
      subtitle = "log2(observed/expected) by region type; * = FDR < 0.05",
      x = "Region type",
      y = "log2(Enrichment ratio)",
      fill = "Event type"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))

  print(p1)

  # Plot 1b: Faceted by event type for readability
  p1b <- ggplot(plot_data, aes(x = region_type, y = log2_enrich, fill = direction_enrich)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    facet_wrap(~event_type, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("Enriched" = "firebrick", "Depleted" = "steelblue", "NS" = "grey60")) +
    labs(
      title = "Regional Enrichment by Event Type",
      subtitle = "Red = enriched, Blue = depleted, Grey = NS",
      x = "Region type",
      y = "log2(Enrichment ratio)",
      fill = "Direction"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          plot.title = element_text(face = "bold"))

  print(p1b)
}

# Plot 2: ORF boundary susceptibility (forest plot style)
if (nrow(orf_susceptibility) > 0 && !is.na(orf_susceptibility$odds_ratio[1])) {
  p2 <- ggplot(orf_susceptibility, aes(x = exon_class, y = odds_ratio)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = or_ci_lower, ymax = or_ci_upper), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    labs(
      title = "ORF Boundary Susceptibility to A5SS/A3SS",
      subtitle = sprintf("Fisher's exact OR = %.2f, p = %.2e",
                         fisher_result$estimate, fisher_result$p.value),
      x = "Exon class",
      y = "Odds ratio (with 95% CI)"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"))

  print(p2)
}

# Plot 3: ORF impact summary
if (all(orf_impact$n_events > 0)) {
  p3 <- ggplot(orf_impact, aes(x = event_type, y = pct_affecting_orf, fill = event_type)) +
    geom_col(width = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct_affecting_orf,
                                   n_affecting_orf, n_events)),
              vjust = -0.3, size = 3.5) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "ORF Impact of Terminal Events",
      subtitle = "Fraction of terminal events affecting ORF boundary exons",
      x = "Event type",
      y = "% affecting ORF boundary"
    ) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))

  print(p3)
}

dev.off()

cat(sprintf("  ✓ %s\n", file.path(fig_dir, "regional_enrichment.pdf")))

# ═══════════════════════════════════════════════════════════════════
# SECTION 7: Save Outputs
# ═══════════════════════════════════════════════════════════════════

cat("\nSaving results...\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

out_path <- file.path(output_dir, "event_regions.rds")
saveRDS(event_regions, out_path)
cat(sprintf("  ✓ %s\n", out_path))

out_path <- file.path(output_dir, "regional_distribution_results.tsv")
write_tsv(regional_distribution, out_path)
cat(sprintf("  ✓ %s\n", out_path))

out_path <- file.path(output_dir, "orf_boundary_susceptibility.tsv")
write_tsv(orf_susceptibility, out_path)
cat(sprintf("  ✓ %s\n", out_path))

out_path <- file.path(output_dir, "orf_impact_summary.tsv")
write_tsv(orf_impact, out_path)
cat(sprintf("  ✓ %s\n", out_path))

cat("\n✓ Step 12 complete\n\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Events analyzed: %d\n", nrow(all_events)))
cat(sprintf("Event-region mappings: %d\n", nrow(event_regions)))
cat(sprintf("Significant enrichments: %d\n", sum(regional_distribution$significant, na.rm = TRUE)))
cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")
