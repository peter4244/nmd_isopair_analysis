#!/usr/bin/env Rscript
# Generate a single-panel supplemental figure showing all 12 Isopair splice event types
# Output: figures/supplemental_event_glossary.png / .pdf

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
source("visualization_functions.R")
library(ggplot2)
library(patchwork)

theme_set(theme_bw(base_size = 16))

# Helper: create a compact event diagram
event_panel <- function(ref_s, ref_e, comp_s, comp_e, title, subtitle) {
  plot_ucsc_pair(ref_s, ref_e, comp_s, comp_e, strand = "+",
                 ref_label = "Ref", comp_label = "Comp",
                 title = title, subtitle = subtitle) +
    theme(plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 13),
          plot.margin = margin(2, 4, 2, 4),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

# Use wider intron gaps (>=300 bp) so directional arrows render consistently.
# Exon width ~150 bp, intron gaps ~300 bp.

# ---- Terminal Events ----
p_tss <- event_panel(
  c(100, 550, 1000), c(250, 700, 1150),
  c(350, 550, 1000), c(500, 700, 1150),
  "Alt_TSS", "Different transcription start site")

p_tes <- event_panel(
  c(100, 550, 1000), c(250, 700, 1250),
  c(100, 550, 1000), c(250, 700, 1100),
  "Alt_TES", "Different transcription end site")

# ---- Boundary Events ----
p_a5ss <- event_panel(
  c(100, 550, 1000), c(250, 700, 1150),
  c(100, 550, 1000), c(250, 780, 1150),
  "A5SS", "Different 5' splice site (donor)")

p_a3ss <- event_panel(
  c(100, 550, 1000), c(250, 700, 1150),
  c(100, 470, 1000), c(250, 700, 1150),
  "A3SS", "Different 3' splice site (acceptor)")

p_pir5 <- event_panel(
  c(100, 550, 1000), c(250, 700, 1150),
  c(100, 550, 1000), c(250, 800, 1150),
  "Partial_IR_5", "Partial intron retention (5' side)")

p_pir3 <- event_panel(
  c(100, 550, 1000), c(250, 700, 1150),
  c(100, 550, 850),  c(250, 700, 1150),
  "Partial_IR_3", "Partial intron retention (3' side)")

# ---- Exon-Level Events ----
p_se <- event_panel(
  c(100, 550, 900, 1250), c(250, 700, 1050, 1400),
  c(100, 550, 1250),      c(250, 700, 1400),
  "SE", "Skipped exon (comparator LOSS)")

p_mi <- event_panel(
  c(100, 450, 750, 1050, 1400), c(250, 600, 900, 1200, 1550),
  c(100, 1400),                  c(250, 1550),
  "Missing_Internal", "Multiple internal exons skipped")

# ---- Intron Retention Events ----
p_ir <- event_panel(
  c(100, 550, 1000), c(250, 700, 1150),
  c(100, 550),       c(250, 1150),
  "IR", "Full intron retention (exact boundaries)")

p_ird5 <- event_panel(
  c(100, 550, 1000), c(250, 700, 1150),
  c(100, 600),       c(250, 1150),
  "IR_diff_5", "Intron retention, 5' boundary differs")

p_ird3 <- event_panel(
  c(100, 550, 1000), c(250, 700, 1150),
  c(100, 550),       c(250, 1100),
  "IR_diff_3", "Intron retention, 3' boundary differs")

p_ird53 <- event_panel(
  c(100, 550, 1000), c(250, 700, 1150),
  c(100, 600),       c(250, 1100),
  "IR_diff_5_3", "Intron retention, both boundaries differ")

# ---- Compose into a 4x3 grid ----
combined <- (p_tss | p_tes | p_a5ss) /
            (p_a3ss | p_pir5 | p_pir3) /
            (p_se | p_mi | p_ir) /
            (p_ird5 | p_ird3 | p_ird53) +
  plot_annotation(
    title = "Isopair Splice Event Types",
    subtitle = "Reference (top) vs Comparator (bottom). Events described from comparator perspective: GAIN = comparator has additional sequence; LOSS = comparator missing sequence.",
    theme = theme(plot.title = element_text(size = 20, face = "bold"),
                  plot.subtitle = element_text(size = 14)))

ggsave("figures/supplemental_event_glossary.png", combined,
       width = 16, height = 16, dpi = 300)
ggsave("figures/supplemental_event_glossary.pdf", combined,
       width = 16, height = 16)

cat("Saved: figures/supplemental_event_glossary.png and .pdf\n")
