#!/usr/bin/env Rscript
#
# Script: 11_analyze_spatial_patterns.R
# Version: 6.0
# Purpose: Analyze spatial organization of splicing events (positional bias, topology, proximity)
#
# Input:
#   - splicing_choice_profiles.rds (from Script 08)
#   - data/isoform_structures_filtered.rds (for gene span and exon boundaries)
#
# Output:
#   - {output_dir}/positional_bias_results.tsv
#   - {output_dir}/topology_enrichment_results.tsv
#   - {output_dir}/proximity_analysis_results.tsv
#   - {output_dir}/figures/spatial_patterns.pdf
#
# Flags:
#   --input path.rds       Input profiles (default: data/splicing_choice_profiles.rds)
#   --output-dir dir/      Output directory (default: results)
#

library(tidyverse)

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 11: Analyze Spatial Organization                       ║\n")
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

# Load isoform structures for gene span and exon boundary mapping
structures <- readRDS(file.path(data_dir, "isoform_structures_filtered.rds"))
cat(sprintf("  Loaded %d isoform structures\n", nrow(structures)))

# ═══════════════════════════════════════════════════════════════════
# SECTION 2: Extract All Events with Genomic Coordinates
# ═══════════════════════════════════════════════════════════════════

cat("\nExtracting events with genomic coordinates...\n")

# Extract all detailed_events into a single data frame
# Note: detailed_events contains gene_id, dominant_transcript_id, comparator_transcript_id
# which would collide with outer columns. Keep outer IDs, drop gene_id from select
# since it's already inside detailed_events.
all_events <- profiles %>%
  select(dominant_isoform_id, non_dominant_isoform_id, detailed_events) %>%
  filter(map_int(detailed_events, nrow) > 0) %>%
  unnest(detailed_events)

cat(sprintf("  Total events: %d\n", nrow(all_events)))

# Convert strand-aware five_prime/three_prime to genomic coordinates
all_events <- all_events %>%
  mutate(
    genomic_start = pmin(five_prime, three_prime),
    genomic_end = pmax(five_prime, three_prime),
    event_midpoint = (genomic_start + genomic_end) / 2
  )

# Compute gene span from dominant isoform structure
# Keep dominant_isoform_id to avoid row multiplication when a gene has multiple dominants
gene_spans <- structures %>%
  select(isoform_id, gene_id, tx_start, tx_end) %>%
  inner_join(
    profiles %>% distinct(gene_id, dominant_isoform_id),
    by = c("isoform_id" = "dominant_isoform_id", "gene_id")
  ) %>%
  select(dominant_isoform_id = isoform_id, gene_id,
         gene_start = tx_start, gene_end = tx_end)

# Join gene spans to events by both gene_id and dominant_isoform_id
all_events <- all_events %>%
  left_join(gene_spans, by = c("gene_id", "dominant_isoform_id" = "dominant_isoform_id")) %>%
  filter(!is.na(gene_start), !is.na(gene_end))

# Compute relative position (0 = gene start, 1 = gene end, in genomic coords)
all_events <- all_events %>%
  mutate(
    gene_length = gene_end - gene_start,
    relative_pos = ifelse(gene_length > 0,
                          (event_midpoint - gene_start) / gene_length,
                          0.5)
  ) %>%
  # Clamp to [0,1] in case events extend slightly beyond dominant structure
  mutate(relative_pos = pmax(0, pmin(1, relative_pos)))

cat(sprintf("  Events with valid positions: %d\n", sum(!is.na(all_events$relative_pos))))

# ═══════════════════════════════════════════════════════════════════
# SECTION 3: Positional Bias Analysis
# ═══════════════════════════════════════════════════════════════════

cat("\nTesting positional bias by event type...\n")

event_type_counts <- all_events %>%
  count(event_type) %>%
  arrange(desc(n))

cat("  Event type counts:\n")
print(event_type_counts, n = Inf)

# KS test vs uniform(0,1) for each event type with n >= 30
positional_results <- event_type_counts %>%
  filter(n >= 30) %>%
  pull(event_type) %>%
  map_dfr(function(et) {
    positions <- all_events %>%
      filter(event_type == et) %>%
      pull(relative_pos)

    ks <- ks.test(positions, "punif", 0, 1)

    tibble(
      event_type = et,
      n_events = length(positions),
      median_position = median(positions),
      mean_position = mean(positions),
      sd_position = sd(positions),
      ks_D = ks$statistic,
      ks_p_value = ks$p.value
    )
  }) %>%
  mutate(
    ks_adj_p = p.adjust(ks_p_value, method = "fdr"),
    significant = ks_adj_p < 0.05,
    bias = case_when(
      !significant ~ "None",
      median_position < 0.3 ~ "5' biased",
      median_position > 0.7 ~ "3' biased",
      TRUE ~ "Central/other"
    )
  )

cat("\nPositional bias results:\n")
print(positional_results %>% select(event_type, n_events, median_position, ks_D, ks_adj_p, bias), n = Inf)

# ═══════════════════════════════════════════════════════════════════
# SECTION 4: Topology Enrichment (A5SS + A3SS)
# ═══════════════════════════════════════════════════════════════════

cat("\nTesting topology enrichment (A5SS + A3SS configurations)...\n")

# Find profiles with both A5SS and A3SS events
profiles_both <- profiles %>%
  filter(n_a5ss > 0, n_a3ss > 0)

cat(sprintf("  Profiles with both A5SS and A3SS: %d\n", nrow(profiles_both)))

# Build lookup of exon boundaries per isoform
# For each dominant isoform, we need its exon starts/ends
dom_structures <- structures %>%
  filter(isoform_id %in% profiles_both$dominant_isoform_id)

# Function to classify an A5SS-A3SS pair as F2F, B2B, or Interleaved
# F2F: A5SS and A3SS flank the same intron (A5SS on exon i 3' boundary, A3SS on exon i+1 5' boundary)
# B2B: A5SS and A3SS on the same exon (A3SS on 5' boundary, A5SS on 3' boundary)
classify_topology <- function(a5ss_event, a3ss_event, exon_starts, exon_ends, gene_strand) {

  # Get genomic coordinates of event midpoints
  a5ss_mid <- (min(a5ss_event$five_prime, a5ss_event$three_prime) +
                max(a5ss_event$five_prime, a5ss_event$three_prime)) / 2
  a3ss_mid <- (min(a3ss_event$five_prime, a3ss_event$three_prime) +
                max(a3ss_event$five_prime, a3ss_event$three_prime)) / 2

  # Find which exon boundary each event is nearest to
  # A5SS should be near a 3' exon boundary (exon_ends on + strand, exon_starts on - strand)
  # A3SS should be near a 5' exon boundary (exon_starts on + strand, exon_ends on - strand)
  if (gene_strand == "+") {
    a5ss_exon_idx <- which.min(abs(exon_ends - a5ss_mid))
    a3ss_exon_idx <- which.min(abs(exon_starts - a3ss_mid))
  } else {
    # On minus strand, 5' end is at higher genomic coordinate
    a5ss_exon_idx <- which.min(abs(exon_starts - a5ss_mid))
    a3ss_exon_idx <- which.min(abs(exon_ends - a3ss_mid))
  }

  if (a5ss_exon_idx == a3ss_exon_idx) {
    return("B2B")  # Same exon
  } else if (gene_strand == "+" && a5ss_exon_idx == a3ss_exon_idx - 1) {
    return("F2F")  # A5SS on exon i, A3SS on exon i+1 (flanking same intron)
  } else if (gene_strand == "-" && a5ss_exon_idx == a3ss_exon_idx + 1) {
    return("F2F")  # Minus strand: reversed genomic order
  } else {
    return("Interleaved")
  }
}

# Classify each A5SS-A3SS pair
topology_results_list <- list()

for (i in seq_len(nrow(profiles_both))) {
  prof <- profiles_both[i, ]
  events <- prof$detailed_events[[1]]

  a5ss_events <- events %>% filter(event_type == "A5SS")
  a3ss_events <- events %>% filter(event_type == "A3SS")

  # Get dominant isoform structure
  dom_struct <- dom_structures %>% filter(isoform_id == prof$dominant_isoform_id)
  if (nrow(dom_struct) == 0) next

  exon_starts <- dom_struct$exon_starts[[1]]
  exon_ends <- dom_struct$exon_ends[[1]]
  gene_strand <- dom_struct$strand

  # Classify each A5SS-A3SS pair
  for (a5 in seq_len(nrow(a5ss_events))) {
    for (a3 in seq_len(nrow(a3ss_events))) {
      topo <- classify_topology(
        a5ss_events[a5, ], a3ss_events[a3, ],
        exon_starts, exon_ends, gene_strand
      )
      topology_results_list[[length(topology_results_list) + 1]] <- tibble(
        gene_id = prof$gene_id,
        dominant_isoform_id = prof$dominant_isoform_id,
        topology = topo
      )
    }
  }
}

topology_classified <- bind_rows(topology_results_list)

if (nrow(topology_classified) > 0) {
  topo_counts <- topology_classified %>%
    count(topology) %>%
    mutate(prop = n / sum(n))

  cat("\nTopology classification:\n")
  print(topo_counts)

  # Permutation test: shuffle A5SS/A3SS event positions within each gene's exon boundaries
  # and re-classify to build null distribution
  n_perm <- 1000
  cat(sprintf("\nRunning %d permutations for topology null distribution...\n", n_perm))

  observed_f2f_prop <- topo_counts %>% filter(topology == "F2F") %>% pull(prop)
  if (length(observed_f2f_prop) == 0) observed_f2f_prop <- 0

  perm_f2f_props <- numeric(n_perm)

  for (perm in seq_len(n_perm)) {
    perm_topos <- character(0)

    for (i in seq_len(nrow(profiles_both))) {
      prof <- profiles_both[i, ]
      events <- prof$detailed_events[[1]]
      n_a5 <- sum(events$event_type == "A5SS")
      n_a3 <- sum(events$event_type == "A3SS")

      dom_struct <- dom_structures %>% filter(isoform_id == prof$dominant_isoform_id)
      if (nrow(dom_struct) == 0) next

      exon_starts <- dom_struct$exon_starts[[1]]
      exon_ends <- dom_struct$exon_ends[[1]]
      n_exons <- length(exon_starts)
      gene_strand <- dom_struct$strand

      if (n_exons < 2) next

      # Randomly assign A5SS to exon 3' boundaries and A3SS to exon 5' boundaries
      if (gene_strand == "+") {
        a5ss_indices <- sample(seq_len(n_exons), n_a5, replace = TRUE)
        a3ss_indices <- sample(seq_len(n_exons), n_a3, replace = TRUE)
      } else {
        a5ss_indices <- sample(seq_len(n_exons), n_a5, replace = TRUE)
        a3ss_indices <- sample(seq_len(n_exons), n_a3, replace = TRUE)
      }

      for (a5_idx in a5ss_indices) {
        for (a3_idx in a3ss_indices) {
          if (a5_idx == a3_idx) {
            perm_topos <- c(perm_topos, "B2B")
          } else if (gene_strand == "+" && a5_idx == a3_idx - 1) {
            perm_topos <- c(perm_topos, "F2F")
          } else if (gene_strand == "-" && a5_idx == a3_idx + 1) {
            perm_topos <- c(perm_topos, "F2F")
          } else {
            perm_topos <- c(perm_topos, "Interleaved")
          }
        }
      }
    }

    perm_f2f_props[perm] <- sum(perm_topos == "F2F") / max(length(perm_topos), 1)
  }

  # Compute permutation p-value (is observed F2F higher than expected?)
  topo_perm_p <- mean(perm_f2f_props >= observed_f2f_prop)

  topology_enrichment <- tibble(
    topology = topo_counts$topology,
    observed_count = topo_counts$n,
    observed_prop = topo_counts$prop,
    f2f_permutation_p = topo_perm_p,
    n_permutations = n_perm
  )

  cat(sprintf("  F2F observed proportion: %.3f\n", observed_f2f_prop))
  cat(sprintf("  F2F permutation p-value: %.4f\n", topo_perm_p))
} else {
  cat("  No profiles with both A5SS and A3SS — skipping topology analysis\n")
  topology_enrichment <- tibble(
    topology = character(0),
    observed_count = integer(0),
    observed_prop = numeric(0),
    f2f_permutation_p = numeric(0),
    n_permutations = integer(0)
  )
}

# ═══════════════════════════════════════════════════════════════════
# SECTION 5: Proximity Analysis
# ═══════════════════════════════════════════════════════════════════

cat("\nAnalyzing event proximity...\n")

# For profiles with >= 2 internal events (excluding Alt_TSS/Alt_TES),
# test whether events are closer together than expected
internal_events <- all_events %>%
  filter(!event_type %in% c("Alt_TSS", "Alt_TES"))

profiles_multi_internal <- internal_events %>%
  group_by(gene_id, dominant_isoform_id, non_dominant_isoform_id) %>%
  summarise(n_internal = n(), .groups = "drop") %>%
  filter(n_internal >= 2)

cat(sprintf("  Profiles with >= 2 internal events: %d\n", nrow(profiles_multi_internal)))

if (nrow(profiles_multi_internal) > 0) {

  # Compute observed mean pairwise distance per profile
  compute_mean_pairwise_dist <- function(midpoints) {
    if (length(midpoints) < 2) return(NA_real_)
    dists <- as.numeric(dist(midpoints))
    mean(dists)
  }

  observed_distances <- internal_events %>%
    inner_join(profiles_multi_internal,
               by = c("gene_id", "dominant_isoform_id", "non_dominant_isoform_id")) %>%
    group_by(gene_id, dominant_isoform_id, non_dominant_isoform_id) %>%
    summarise(
      n_events = n(),
      mean_dist = compute_mean_pairwise_dist(event_midpoint),
      gene_start = first(gene_start),
      gene_end = first(gene_end),
      .groups = "drop"
    ) %>%
    filter(!is.na(mean_dist))

  overall_observed_mean <- mean(observed_distances$mean_dist, na.rm = TRUE)

  # Permutation test: for each profile, permute event positions within gene span
  # Vectorized approach: pre-compute all random positions, then compute distances
  n_perm <- 10000
  cat(sprintf("  Running %d permutations for proximity null...\n", n_perm))

  set.seed(42)

  # Pre-extract vectors for speed
  n_events_vec <- observed_distances$n_events
  gene_start_vec <- observed_distances$gene_start
  gene_end_vec <- observed_distances$gene_end
  n_profiles <- length(n_events_vec)

  # For profiles with exactly 2 events, mean pairwise distance = |x1 - x2|
  # For profiles with more events, use dist()
  # Pre-compute per-profile permuted means in a matrix (n_perm x n_profiles)
  perm_matrix <- matrix(NA_real_, nrow = n_perm, ncol = n_profiles)

  for (j in seq_len(n_profiles)) {
    ne <- n_events_vec[j]
    gs <- gene_start_vec[j]
    ge <- gene_end_vec[j]
    gene_span <- ge - gs

    if (ne == 2) {
      # Special case: distance = |x1 - x2|, mean of uniform order stats
      # Generate all random pairs at once
      r1 <- runif(n_perm, 0, gene_span)
      r2 <- runif(n_perm, 0, gene_span)
      perm_matrix[, j] <- abs(r1 - r2)
    } else {
      # General case: compute mean pairwise distance for ne random points
      for (perm in seq_len(n_perm)) {
        pts <- runif(ne, min = gs, max = ge)
        perm_matrix[perm, j] <- mean(as.numeric(dist(pts)))
      }
    }
  }

  perm_means <- rowMeans(perm_matrix, na.rm = TRUE)

  # One-sided test: are observed distances smaller than expected?
  prox_p_value <- mean(perm_means <= overall_observed_mean)

  proximity_results <- tibble(
    n_profiles = nrow(observed_distances),
    n_events_total = sum(observed_distances$n_events),
    observed_mean_distance = overall_observed_mean,
    permutation_mean_distance = mean(perm_means),
    permutation_sd = sd(perm_means),
    z_score = (overall_observed_mean - mean(perm_means)) / sd(perm_means),
    p_value = prox_p_value,
    n_permutations = n_perm,
    interpretation = ifelse(prox_p_value < 0.05,
                           "Events are closer than expected by chance",
                           "No significant clustering")
  )

  cat(sprintf("  Observed mean distance: %.0f bp\n", overall_observed_mean))
  cat(sprintf("  Permutation mean: %.0f bp (sd = %.0f)\n",
              mean(perm_means), sd(perm_means)))
  cat(sprintf("  z-score: %.2f, p-value: %.4f\n", proximity_results$z_score, prox_p_value))

} else {
  cat("  Insufficient profiles with multiple internal events — skipping proximity analysis\n")
  proximity_results <- tibble(
    n_profiles = 0L,
    n_events_total = 0L,
    observed_mean_distance = NA_real_,
    permutation_mean_distance = NA_real_,
    permutation_sd = NA_real_,
    z_score = NA_real_,
    p_value = NA_real_,
    n_permutations = 0L,
    interpretation = "Insufficient data"
  )
}

# ═══════════════════════════════════════════════════════════════════
# SECTION 6: Visualization
# ═══════════════════════════════════════════════════════════════════

cat("\nCreating spatial plots...\n")

fig_dir <- file.path(output_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(fig_dir, "spatial_patterns.pdf"), width = 14, height = 10)

# Plot 1: Density of relative positions by event type
if (nrow(all_events) > 0) {
  event_order <- positional_results %>%
    arrange(median_position) %>%
    pull(event_type)

  plot_events <- all_events %>%
    filter(event_type %in% event_order) %>%
    mutate(event_type = factor(event_type, levels = event_order))

  p1 <- ggplot(plot_events, aes(x = relative_pos, fill = event_type)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
    facet_wrap(~event_type, scales = "free_y", ncol = 3) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("5'", "25%", "50%", "75%", "3'")) +
    labs(
      title = "Positional Bias: Event Distribution Along Gene Body",
      subtitle = "Dashed line = gene midpoint; x-axis = relative position (0=5' start, 1=3' end)",
      x = "Relative position within gene",
      y = "Density"
    ) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))

  print(p1)
}

# Plot 2: Topology classification barplot
if (nrow(topology_classified) > 0) {
  p2 <- ggplot(topo_counts, aes(x = reorder(topology, -n), y = n, fill = topology)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f%%", prop * 100)), vjust = -0.3, size = 4) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "A5SS-A3SS Topology Classification",
      subtitle = sprintf("F2F permutation p = %.4f (n_perm = %d)",
                         topo_perm_p, n_perm),
      x = "Configuration",
      y = "Count"
    ) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))

  print(p2)

  # Permutation null distribution
  perm_df <- tibble(f2f_prop = perm_f2f_props)
  p2b <- ggplot(perm_df, aes(x = f2f_prop)) +
    geom_histogram(bins = 50, fill = "grey70", color = "white") +
    geom_vline(xintercept = observed_f2f_prop, color = "red", linewidth = 1.2) +
    annotate("text", x = observed_f2f_prop, y = Inf, label = "Observed",
             color = "red", vjust = 2, hjust = -0.1) +
    labs(
      title = "Topology Permutation Null Distribution",
      subtitle = sprintf("Observed F2F proportion = %.3f", observed_f2f_prop),
      x = "F2F proportion (permuted)",
      y = "Count"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"))

  print(p2b)
}

# Plot 3: Proximity analysis — observed vs permuted
if (nrow(profiles_multi_internal) > 0 && exists("perm_means")) {
  perm_dist_df <- tibble(mean_dist = perm_means)
  p3 <- ggplot(perm_dist_df, aes(x = mean_dist)) +
    geom_histogram(bins = 80, fill = "grey70", color = "white") +
    geom_vline(xintercept = overall_observed_mean, color = "red", linewidth = 1.2) +
    annotate("text", x = overall_observed_mean, y = Inf,
             label = sprintf("Observed = %.0f bp", overall_observed_mean),
             color = "red", vjust = 2, hjust = -0.1) +
    labs(
      title = "Event Proximity: Observed vs Permutation Null",
      subtitle = sprintf("p = %.4f (n_perm = %d, n_profiles = %d)",
                         prox_p_value, n_perm, nrow(observed_distances)),
      x = "Mean pairwise distance (bp)",
      y = "Permutation count"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"))

  print(p3)
}

dev.off()

cat(sprintf("  ✓ %s\n", file.path(fig_dir, "spatial_patterns.pdf")))

# ═══════════════════════════════════════════════════════════════════
# SECTION 7: Save Outputs
# ═══════════════════════════════════════════════════════════════════

cat("\nSaving results...\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

out_path <- file.path(output_dir, "positional_bias_results.tsv")
write_tsv(positional_results, out_path)
cat(sprintf("  ✓ %s\n", out_path))

out_path <- file.path(output_dir, "topology_enrichment_results.tsv")
write_tsv(topology_enrichment, out_path)
cat(sprintf("  ✓ %s\n", out_path))

out_path <- file.path(output_dir, "proximity_analysis_results.tsv")
write_tsv(proximity_results, out_path)
cat(sprintf("  ✓ %s\n", out_path))

cat("\n✓ Step 11 complete\n\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Events analyzed: %d\n", nrow(all_events)))
cat(sprintf("Positional bias tests: %d event types (n >= 30)\n", nrow(positional_results)))
cat(sprintf("Topology pairs classified: %d\n", nrow(topology_classified)))
cat(sprintf("Proximity profiles: %d\n", nrow(profiles_multi_internal)))
cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")
