#!/usr/bin/env Rscript
# SF24 — Twelve splice event categories detected by Isopair.
#
# Renders the 12-panel schematic using the canonical `plot_ucsc_pair` function
# from the Isopair pipeline's visualization module. Matches Yul's transcript-
# model display standard (UCSC-browser chevron introns, exon rectangles with
# CDS/UTR coloring, Reference-on-top / Comparator-on-bottom pair convention).
#
# Coordinates are the illustrative examples from the Rmd's Event Type Glossary
# section (`05_final_report_mashr.Rmd`, chunks `glossary-*`), reproduced verbatim.
#
# Output: figure_sf24_splice_event_categories.{pdf,png}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# Source the canonical visualization functions
VIS <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/visualization_functions.R"
source(VIS)

# Convenience helper. Per Yul-style convention, no titles or subtitles on the
# figure itself — the composite caption identifies the panels.
glossary_plot <- function(ref_s, ref_e, comp_s, comp_e, title, subtitle) {
  plot_ucsc_pair(ref_s, ref_e, comp_s, comp_e, strand = "+",
                 ref_label = "Reference", comp_label = "Comparator",
                 title = NULL, subtitle = NULL)
}

# 12 event definitions — coordinates verbatim from the Rmd Event Type Glossary.
p_alt_tss <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(250, 400, 700), c(350, 550, 850),
  "Alt_TSS", "Alternative transcription start site")

p_alt_tes <- glossary_plot(
  c(100, 400, 700), c(200, 550, 900),
  c(100, 400, 700), c(200, 550, 800),
  "Alt_TES", "Alternative transcription end site")

p_a5ss <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 400, 700), c(200, 600, 850),
  "A5SS", "Alternative 5' splice site")

p_a3ss <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 350, 700), c(200, 550, 850),
  "A3SS", "Alternative 3' splice site")

p_se <- glossary_plot(
  c(100, 400, 600, 800), c(200, 500, 700, 900),
  c(100, 400, 800),      c(200, 500, 900),
  "SE", "Skipped exon")

p_mi <- glossary_plot(
  c(100, 350, 550, 750, 950), c(200, 450, 650, 850, 1050),
  c(100, 950),                c(200, 1050),
  "Missing_Internal", "Missing internal exon(s)")

p_pir5 <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 400, 700), c(200, 620, 850),
  "Partial_IR_5", "Partial intron retention (5' side)")

p_pir3 <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 400, 640), c(200, 550, 850),
  "Partial_IR_3", "Partial intron retention (3' side)")

p_ir <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 400),      c(200, 850),
  "IR", "Full intron retention")

p_ir_d5 <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 430),      c(200, 850),
  "IR_diff_5", "IR with 5' boundary difference")

p_ir_d3 <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 400),      c(200, 820),
  "IR_diff_3", "IR with 3' boundary difference")

p_ir_d53 <- glossary_plot(
  c(100, 400, 700), c(200, 550, 850),
  c(100, 430),      c(200, 820),
  "IR_diff_5_3", "IR with 5' and 3' boundary differences")

# 3 rows × 4 cols composite via patchwork.
composite <- (p_alt_tss | p_alt_tes | p_a5ss | p_a3ss) /
             (p_se      | p_mi      | p_pir5 | p_pir3) /
             (p_ir      | p_ir_d5   | p_ir_d3| p_ir_d53) +
             patchwork::plot_annotation(
               theme = theme(plot.margin = margin(6, 6, 6, 6))
             )

HERE <- "/Users/petecastaldi/claude_projects/nmd/figures/SupplementalFigures/SF24_SpliceEventCategories"
ggsave(file.path(HERE, "figure_sf24_splice_event_categories.pdf"),
       composite, width = 15, height = 8.4, device = cairo_pdf)
ggsave(file.path(HERE, "figure_sf24_splice_event_categories.png"),
       composite, width = 15, height = 8.4, dpi = 200)
message("wrote figure_sf24_splice_event_categories.{pdf,png}")
