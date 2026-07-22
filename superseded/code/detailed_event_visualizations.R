#!/usr/bin/env Rscript
# Detailed Event Visualizations
# Creates publication-quality figures for detailed event analysis

library(tidyverse)
library(patchwork)
library(scales)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║           DETAILED EVENT VISUALIZATIONS                        ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

output_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions"

# Load data
cat("Loading data...\n")
transitions <- readRDS(file.path(output_dir, "detailed_event_vectors.rds"))
incidence <- read_tsv(file.path(output_dir, "q1_event_incidence_by_celltype.tsv"),
                      show_col_types = FALSE)

# ============================================================================
# Figure 1: Event Incidence by Cell Type (Detailed)
# ============================================================================

cat("Creating Figure 1: Event incidence by cell type...\n")

# Reshape for plotting
incidence_long <- incidence %>%
  select(cell_type, n_alt_tss, n_alt_tes, n_se, n_a5ss, n_a3ss, n_ir) %>%
  pivot_longer(-cell_type, names_to = "event_type", values_to = "count") %>%
  mutate(
    event_type = str_remove(event_type, "n_"),
    event_type = factor(event_type,
                        levels = c("alt_tss", "alt_tes", "se", "a5ss", "a3ss", "ir"),
                        labels = c("Alt TSS", "Alt TES", "SE", "A5SS", "A3SS", "IR"))
  )

fig1 <- ggplot(incidence_long, aes(x = cell_type, y = count, fill = event_type)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Detailed Event Incidence by Cell Type",
    subtitle = "Count of each specific event type across isoform transitions",
    x = "Cell Type",
    y = "Event Count",
    fill = "Event Type"
  ) +
  scale_y_continuous(labels = comma) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "fig1_detailed_event_incidence.pdf"),
       fig1, width = 10, height = 6)
cat("  Saved: fig1_detailed_event_incidence.pdf\n")

# ============================================================================
# Figure 2: Event Proportion (% of transitions with each event)
# ============================================================================

cat("Creating Figure 2: Event proportions...\n")

incidence_pct <- incidence %>%
  select(cell_type, pct_with_alt_tss, pct_with_alt_tes, pct_with_se,
         pct_with_a5ss, pct_with_a3ss, pct_with_ir) %>%
  pivot_longer(-cell_type, names_to = "event_type", values_to = "pct") %>%
  mutate(
    event_type = str_remove(event_type, "pct_with_"),
    event_type = factor(event_type,
                        levels = c("alt_tss", "alt_tes", "se", "a5ss", "a3ss", "ir"),
                        labels = c("Alt TSS", "Alt TES", "SE", "A5SS", "A3SS", "IR"))
  )

fig2 <- ggplot(incidence_pct, aes(x = cell_type, y = pct, fill = event_type)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Proportion of Transitions with Each Event Type",
    subtitle = "% of isoform transitions containing at least one event of each type",
    x = "Cell Type",
    y = "% of Transitions",
    fill = "Event Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "fig2_event_proportions.pdf"),
       fig2, width = 10, height = 6)
cat("  Saved: fig2_event_proportions.pdf\n")

# ============================================================================
# Figure 3: Event Complexity Distribution
# ============================================================================

cat("Creating Figure 3: Event complexity distribution...\n")

complexity_data <- transitions %>%
  mutate(
    n_base_events_total = n_alt_tss + n_alt_tes + n_se + n_a5ss + n_a3ss + n_ir,
    complexity_class = case_when(
      n_base_events_total == 0 ~ "0 events",
      n_base_events_total == 1 ~ "1 event",
      n_base_events_total == 2 ~ "2 events",
      n_base_events_total == 3 ~ "3 events",
      n_base_events_total >= 4 ~ "4+ events"
    ),
    complexity_class = factor(complexity_class,
                              levels = c("0 events", "1 event", "2 events", "3 events", "4+ events"))
  ) %>%
  group_by(cell_type, complexity_class) %>%
  summarize(n_transitions = n(), .groups = "drop") %>%
  group_by(cell_type) %>%
  mutate(pct = 100 * n_transitions / sum(n_transitions))

fig3 <- ggplot(complexity_data, aes(x = cell_type, y = pct, fill = complexity_class)) +
  geom_col(position = "stack") +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(
    title = "Transition Complexity by Cell Type",
    subtitle = "Distribution of number of base splicing events per transition",
    x = "Cell Type",
    y = "% of Transitions",
    fill = "Complexity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "fig3_complexity_distribution.pdf"),
       fig3, width = 10, height = 6)
cat("  Saved: fig3_complexity_distribution.pdf\n")

# ============================================================================
# Figure 4: Co-occurrence Heatmap
# ============================================================================

cat("Creating Figure 4: Event co-occurrence heatmap...\n")

if (file.exists(file.path(output_dir, "q2_event_cooccurrence_matrix.tsv"))) {
  cooccurrence <- read_tsv(file.path(output_dir, "q2_event_cooccurrence_matrix.tsv"),
                           show_col_types = FALSE)
  
  # Create symmetric matrix for heatmap
  event_types <- c("alt_tss", "alt_tes", "se", "a5ss", "a3ss", "ir")
  cooccur_matrix <- matrix(0, nrow = length(event_types), ncol = length(event_types))
  rownames(cooccur_matrix) <- event_types
  colnames(cooccur_matrix) <- event_types
  
  for (i in 1:nrow(cooccurrence)) {
    row <- cooccurrence[i, ]
    r_idx <- which(event_types == row$event_A)
    c_idx <- which(event_types == row$event_B)
    cooccur_matrix[r_idx, c_idx] <- log10(row$odds_ratio)
    cooccur_matrix[c_idx, r_idx] <- log10(row$odds_ratio)
  }
  diag(cooccur_matrix) <- 0
  
  # Convert to long format
  cooccur_long <- as_tibble(cooccur_matrix, rownames = "event_A") %>%
    pivot_longer(-event_A, names_to = "event_B", values_to = "log_odds_ratio") %>%
    mutate(
      event_A = factor(event_A, levels = rev(event_types),
                       labels = rev(c("Alt TSS", "Alt TES", "SE", "A5SS", "A3SS", "IR"))),
      event_B = factor(event_B, levels = event_types,
                       labels = c("Alt TSS", "Alt TES", "SE", "A5SS", "A3SS", "IR"))
    )
  
  fig4 <- ggplot(cooccur_long, aes(x = event_B, y = event_A, fill = log_odds_ratio)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      name = "Log10\nOdds Ratio"
    ) +
    labs(
      title = "Event Co-occurrence Pattern",
      subtitle = "Fisher's exact test odds ratios for event pair co-occurrence",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    ) +
    coord_fixed()
  
  ggsave(file.path(output_dir, "fig4_cooccurrence_heatmap.pdf"),
         fig4, width = 8, height = 7)
  cat("  Saved: fig4_cooccurrence_heatmap.pdf\n")
} else {
  cat("  Skipping Figure 4: Q2 results not yet available\n")
}

# ============================================================================
# Figure 5: Event Distribution Comparison
# ============================================================================

cat("Creating Figure 5: Event type comparison across cell types...\n")

# Calculate mean events per transition
mean_events <- transitions %>%
  group_by(cell_type) %>%
  summarize(
    mean_alt_tss = mean(n_alt_tss),
    mean_alt_tes = mean(n_alt_tes),
    mean_se = mean(n_se),
    mean_a5ss = mean(n_a5ss),
    mean_a3ss = mean(n_a3ss),
    mean_ir = mean(n_ir),
    .groups = "drop"
  ) %>%
  pivot_longer(-cell_type, names_to = "event_type", values_to = "mean_count") %>%
  mutate(
    event_type = str_remove(event_type, "mean_"),
    event_type = factor(event_type,
                        levels = c("alt_tss", "alt_tes", "se", "a5ss", "a3ss", "ir"),
                        labels = c("Alt TSS", "Alt TES", "SE", "A5SS", "A3SS", "IR"))
  )

fig5 <- ggplot(mean_events, aes(x = event_type, y = mean_count, fill = cell_type)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Paired") +
  labs(
    title = "Mean Event Count per Transition by Cell Type",
    subtitle = "Average number of each event type per isoform transition",
    x = "Event Type",
    y = "Mean Events per Transition",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "fig5_mean_events_by_celltype.pdf"),
       fig5, width = 10, height = 6)
cat("  Saved: fig5_mean_events_by_celltype.pdf\n")

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║           VISUALIZATION GENERATION COMPLETE                    ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("Generated 5 figures:\n")
cat("  1. fig1_detailed_event_incidence.pdf\n")
cat("  2. fig2_event_proportions.pdf\n")
cat("  3. fig3_complexity_distribution.pdf\n")
cat("  4. fig4_cooccurrence_heatmap.pdf\n")
cat("  5. fig5_mean_events_by_celltype.pdf\n")
cat("\n")
