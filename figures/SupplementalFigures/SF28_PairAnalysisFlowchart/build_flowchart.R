#!/usr/bin/env Rscript
# Pair-analysis cohort flowchart — standalone supplement render.
#
# Plain-language flow chart of how the pair analysis arrives at the 3,009
# matched NMD-sensitive / Control pairs (upstream chain) and how those pairs
# then narrow into the two analysis subsets used in the manuscript
# (downstream cascade), with kept / dropped counts and the reason for each
# dropout at every filter step.
#
# Outputs (in this folder):
#   figure_s_pair_analysis_flowchart.dot   — raw Graphviz DOT spec
#   figure_s_pair_analysis_flowchart.html  — standalone htmlwidget; opens
#                                            in any browser
#
# For PNG / PDF: open the HTML in Chrome and File → Print → "Save as PDF",
# or paste the .dot file into https://dreampuf.github.io/GraphvizOnline/.
#
# Numbers match the inline §1 chunk of
#   results/isoform_transitions/Version_6.0/isopair_wrapper/
#     05_final_report_gencode_scope_2026-06-15.Rmd
# and the figure-side TSVs; cross-pipeline drift is guarded by the verifiers
# at reproducibility/verify_pass7_new_rmd.R and
# verify_cross_check_new_rmd_vs_figures.R.

suppressPackageStartupMessages({
  library(DiagrammeR)
  library(htmlwidgets)
})

if (Sys.getenv("RSTUDIO_PANDOC") == "" && file.exists(
      "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/pandoc")) {
  Sys.setenv(RSTUDIO_PANDOC =
             "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64")
}

NMD_ROOT <- Sys.getenv("NMD_ROOT",
                        "/Users/petecastaldi/claude_projects/nmd")
DM <- file.path(NMD_ROOT,
                 "results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr")
HERE <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- args[grep("--file=", args)]
  if (length(m) > 0) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  normalizePath(getwd())
})()

# ─────────────────────────────────────────────────────────────────────
# Live counts pulled from upstream RDS caches (kept transparent).
# ─────────────────────────────────────────────────────────────────────
expr_mat    <- readRDS(file.path(DM, "expression_data.rds"))
nmd_class   <- readRDS(file.path(DM, "nmd_classification.rds"))
profiles_c2 <- readRDS(file.path(DM, "profiles_c2_allsamples.rds"))
profiles_c4 <- readRDS(file.path(DM, "profiles_c4_allsamples.rds"))

# Upstream constants — fixed (raw isocall matrix dimensions are constants of
# the input data and are not stored in any RDS; recorded in the legacy Rmd's
# flowchart-data chunk).
N_RAW_ISO     <- 645273L
N_RAW_SMP     <- 38L
N_AFTER_OUT   <- 36L       # 38 - 2 DO outliers
N_FILTER_ISO  <- nrow(expr_mat)           # 123,790
N_FILTER_SMP  <- ncol(expr_mat)           # 36 (all 6 CTs at this stage)
N_4CT_SMP     <- 26L                       # AT 6 + DD 8 + FB 6 + MV 6
N_DROPPED_SMP <- N_FILTER_SMP - N_4CT_SMP  # 10 (6 DD_ALI + 4 DO)
N_NMD_AS      <- length(nmd_class[["all_samples"]]$nmd)
N_NONNMD_AS   <- length(nmd_class[["all_samples"]]$non_nmd)
N_C2          <- nrow(profiles_c2)         # 3,009
N_C4          <- nrow(profiles_c4)         # 8,323

fmt <- function(x) format(x, big.mark = ",")

# ─────────────────────────────────────────────────────────────────────
# DOT spec — three logical zones:
#   1. Upstream pipeline (single column, top)
#   2. Matched-pair hub (3,009 / 3,009)
#   3. Two subset clusters (left = strict / right = broader), each with
#      cluster-level header as its title and filter nodes inside
#   + Hidden-PTC sub-cascade off Subset 2
#
# Box formatting: each filter node uses a small HTML table inside the
# label so the kept / dropped / why rows align cleanly.
# ─────────────────────────────────────────────────────────────────────
filter_node <- function(what, kept, dropped, why) {
  # 2-column TABLE — left column = label, right column = value.
  # Eliminates the "<B>72</B> PTC+" kerning trap (graphviz collapses the
  # space after </B>); proper grid layout aligns counts visually.
  sprintf(
    '<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
       <TR><TD ALIGN="LEFT" COLSPAN="2"><B>Filter:</B>&nbsp;%s</TD></TR>
       <TR>
         <TD ALIGN="LEFT"><FONT COLOR="#1a5a2a"><B>kept:</B></FONT></TD>
         <TD ALIGN="LEFT"><FONT COLOR="#1a5a2a"><B>%s</B></FONT></TD>
       </TR>
       <TR>
         <TD ALIGN="LEFT"><FONT COLOR="#b22222">dropped:</FONT></TD>
         <TD ALIGN="LEFT"><FONT COLOR="#b22222">%s</FONT></TD>
       </TR>
       <TR>
         <TD ALIGN="LEFT" COLSPAN="2"><FONT POINT-SIZE="14" COLOR="#444">Why: %s</FONT></TD>
       </TR>
     </TABLE>>', what, kept, dropped, why)
}

upstream_node <- function(stage, detail, what_dropped = NULL) {
  rows <- sprintf(
    '<TR><TD ALIGN="LEFT"><B>%s</B></TD></TR><TR><TD ALIGN="LEFT">%s</TD></TR>',
    stage, detail)
  if (!is.null(what_dropped)) {
    rows <- paste0(rows, sprintf(
      '<TR><TD ALIGN="LEFT"><FONT POINT-SIZE="14" COLOR="#444">%s</FONT></TD></TR>',
      what_dropped))
  }
  sprintf('<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">%s</TABLE>>', rows)
}

dot_spec <- sprintf('
digraph cohort_flow {
  graph [rankdir=TB, fontname="Helvetica", nodesep=0.30, ranksep=0.45,
         splines=ortho, compound=true]
  node  [shape=box, style="rounded,filled", fontname="Helvetica",
         fontsize=18, margin="0.25,0.18"]
  edge  [color="#666", penwidth=1.4, fontname="Helvetica",
         fontsize=14, fontcolor="#444"]

  // ─── Starting cohort + upstream classification + pair construction
  U1 [label=%s, fillcolor="#e8eef4", color="#2f4858"]
  U2 [label=%s, fillcolor="#fff7bc", color="#d95f0e"]
  U3 [label=%s, fillcolor="#fff7bc", color="#d95f0e"]

  U1 -> U2 -> U3

  // ─── Hub: matched isoform pairs ──────────────────────────────────
  HUB [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                <TR><TD ALIGN="CENTER"><FONT COLOR="white" POINT-SIZE="22"><B>Matched isoform pairs</B></FONT></TD></TR>
                <TR><TD ALIGN="CENTER"><FONT COLOR="white"><B>3,009 NMD susceptible  ·  3,009 Control</B></FONT></TD></TR>
                <TR><TD ALIGN="CENTER"><FONT COLOR="#cfd8e2" POINT-SIZE="14">one triplet per gene</FONT></TD></TR>
              </TABLE>>,
       fillcolor="#1f3a5f", color="#0a1a30", penwidth=2.0]

  U3 -> HUB

  // ─── Subset 1 cluster (left) ─────────────────────────────────────
  subgraph cluster_S1 {
    label = <<FONT POINT-SIZE="22"><B>GENCODE-restricted subset (n = 190)</B></FONT><BR/>all three transcripts in the triplet are annotated coding transcripts>
    labelloc = "t"
    style = "rounded,filled"
    fillcolor = "#fef0e6"
    color = "#7d4030"
    fontcolor = "#5a2d1e"
    margin = 12

    S1_F1 [label=%s, fillcolor="white", color="#a2604a"]
    S1_F2 [label=%s, fillcolor="white", color="#a2604a"]

    S1_GROUPS [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                       <TR><TD ALIGN="LEFT" COLSPAN="2"><B>Classification:</B>&nbsp;premature stop in the NMD-susceptible transcript?</TD></TR>
                       <TR><TD ALIGN="RIGHT"><B>72</B></TD><TD ALIGN="LEFT">PTC+ (38%%)</TD></TR>
                       <TR><TD ALIGN="RIGHT"><B>118</B></TD><TD ALIGN="LEFT">PTC−</TD></TR>
                       <TR><TD ALIGN="RIGHT"><B>190</B></TD><TD ALIGN="LEFT">Control</TD></TR>
                     </TABLE>>,
               fillcolor="#fff7bc", color="#d95f0e"]

    S1_CONS [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                     <TR><TD ALIGN="LEFT"><B>Used by:</B></TD></TR>
                     <TR><TD ALIGN="LEFT">Figure 3 (Panels D / E / F)</TD></TR>
                     <TR><TD ALIGN="LEFT">Figure 4 (Panels A / B)</TD></TR>
                     <TR><TD ALIGN="LEFT">SF32</TD></TR>
                   </TABLE>>,
             fillcolor="#eeeeee", color="#666"]

    S1_F1 -> S1_F2 -> S1_GROUPS -> S1_CONS [style=dashed, color="#7d4030"]
  }

  // ─── Subset 2 cluster (right) ────────────────────────────────────
  subgraph cluster_S2 {
    label = <<FONT POINT-SIZE="22"><B>Reference-AUG-traceable subset (n = 1,166)</B></FONT><BR/>reference is annotated; NMD / Control can be novel if the reference start codon maps to the comparator>
    labelloc = "t"
    style = "rounded,filled"
    fillcolor = "#e6f1f8"
    color = "#1f5570"
    fontcolor = "#0e2c3e"
    margin = 12

    S2_F1 [label=%s, fillcolor="white", color="#406080"]
    S2_F2 [label=%s, fillcolor="white", color="#406080"]
    S2_F3 [label=%s, fillcolor="white", color="#406080"]

    S2_GROUPS [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                       <TR><TD ALIGN="LEFT" COLSPAN="2"><B>Classification:</B>&nbsp;reference start codon projected to a premature stop?</TD></TR>
                       <TR><TD ALIGN="RIGHT"><B>1,050</B></TD><TD ALIGN="LEFT">PTC+ (90%%)</TD></TR>
                       <TR><TD ALIGN="RIGHT"><B>116</B></TD><TD ALIGN="LEFT">PTC−</TD></TR>
                       <TR><TD ALIGN="RIGHT"><B>1,166</B></TD><TD ALIGN="LEFT">Control</TD></TR>
                     </TABLE>>,
               fillcolor="#fff7bc", color="#d95f0e"]

    S2_CONS [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                     <TR><TD ALIGN="LEFT"><B>Used by:</B></TD></TR>
                     <TR><TD ALIGN="LEFT">Figure 4 (Panels C / D)</TD></TR>
                     <TR><TD ALIGN="LEFT">SF30 / SF31</TD></TR>
                     <TR><TD ALIGN="LEFT">SF35</TD></TR>
                   </TABLE>>,
             fillcolor="#eeeeee", color="#666"]

    S2_F1 -> S2_F2 -> S2_F3 -> S2_GROUPS -> S2_CONS [style=dashed, color="#1f5570"]
  }

  // ─── Hidden-PTC sub-cascade ─────────────────────────────────────
  HPTC [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                <TR><TD ALIGN="LEFT"><B>Occult-PTC subset (n = 492)</B></TD></TR>
                <TR><TD ALIGN="LEFT">reference-anchored analysis flags a premature stop,</TD></TR>
                <TR><TD ALIGN="LEFT">but the TransDecoder2 CDS caller does not</TD></TR>
                <TR><TD ALIGN="LEFT"><FONT COLOR="#b22222">674 / 1,166 dropped</FONT> — no disagreement</TD></TR>
              </TABLE>>,
        fillcolor="#a48bbf", color="#4a3556", fontcolor="#1a0a25"]

  HPTC_CONS [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                     <TR><TD ALIGN="LEFT"><B>Used by:</B></TD></TR>
                     <TR><TD ALIGN="LEFT">SF33 / SF34</TD></TR>
                   </TABLE>>,
             fillcolor="#eeeeee", color="#666"]

  // ─── Direct-from-hub consumers (whole pair set, no further narrowing) ─
  HUB_CONS [label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">
                    <TR><TD ALIGN="LEFT"><B>Used by (whole pair set):</B></TD></TR>
                    <TR><TD ALIGN="LEFT">SF25 / SF26 / SF27</TD></TR>
                    <TR><TD ALIGN="LEFT">SF29</TD></TR>
                  </TABLE>>,
            fillcolor="#eeeeee", color="#666"]

  // ─── Edges from hub into subset clusters and TD2 sub-cascade ─────
  // HUB → cluster arrows: no tailport specifier on either edge (mixed ports
  // like sw vs se on a single source produce visually asymmetric routes;
  // default auto-routing under splines=ortho + symmetric cluster targets
  // gives a balanced fan-out). minlen=2 forces each edge to span 2 ranks
  // for visible length. Both rules enforced by
  // figures/lib/validate_flowchart_dot.R.
  HUB -> S1_F1 [lhead=cluster_S1, color="#7d4030", penwidth=2.0, minlen=2]
  HUB -> S2_F1 [lhead=cluster_S2, color="#1f5570", penwidth=2.0, minlen=2]
  // HUB_CONS placed to the RIGHT of HUB (same rank), not below — the
  // whole-pair-set descriptive analyses are a sibling of the two subset
  // narrowings, so a lateral placement reads more accurately.
  HUB -> HUB_CONS [style=dashed, color="#555", constraint=false]
  { rank = same; HUB; HUB_CONS }

  S2_F3 -> HPTC [color="#4a3556"]
  HPTC -> HPTC_CONS [style=dashed, color="#4a3556"]
}',
  upstream_node("Starting cohort",
                sprintf("AT (n = 6) · DD (n = 8) · FB (n = 6) · MV (n = 6) · %d samples · %s expressed isoforms",
                        N_4CT_SMP, fmt(N_FILTER_ISO))),
  upstream_node("Identify NMD susceptible isoforms",
                sprintf("%s NMD susceptible · %s non-NMD",
                        fmt(N_NMD_AS), fmt(N_NONNMD_AS))),
  upstream_node("Construct paired comparisons",
                sprintf("per gene: dominant NMD-susceptible vs reference · two dominant non-NMD (Control + reference) → %s NMD-susceptible pairs · %s Control pairs (before matching)",
                        fmt(N_C2), fmt(N_C4))),

  # Subset 1 filter nodes — Why = single clause, no in-cell line wrapping
  filter_node("all three transcripts are annotated (in GENCODE)",
              "301 / 301", "2,708 / 2,708",
              "most NMD-susceptible transcripts are novel"),
  filter_node("all three transcripts have a curated coding sequence",
              "190 / 190", "111 / 111",
              "annotated non-coding transcripts (lncRNA, etc.)"),

  # Subset 2 filter nodes
  filter_node("reference is annotated (NMD and Control can be novel)",
              "2,098 / 2,098", "911 / 911",
              "reference itself is a novel isoform"),
  filter_node("the reference start codon maps onto the comparator transcript",
              "1,659 NMD  ·  1,286 Control", "439 NMD  ·  812 Control",
              "reference exon missing from comparator"),
  filter_node("keep only genes where both NMD and Control side passed",
              "1,166 / 1,166", "493 NMD  ·  120 Control",
              "reference passed on one arm but not the other")
)

dot_path  <- file.path(HERE, "figure_s_pair_analysis_flowchart.dot")
html_path <- file.path(HERE, "figure_s_pair_analysis_flowchart.html")

# Static flowchart validator — catches recurring formatting traps before render
source(file.path(NMD_ROOT, "figures/lib/validate_flowchart_dot.R"))
v <- validate_flowchart_dot(dot_spec)
if (length(v$errors) > 0) {
  stop(sprintf("validate_flowchart_dot found %d error(s) — fix before render",
                length(v$errors)))
}

writeLines(dot_spec, dot_path)
cat(sprintf("[dot]  → %s\n", dot_path))

# Explicit container size matching the SVG's natural extent at the chosen
# font sizes — prevents the htmlwidget's small default container from
# scaling the SVG (which makes fonts look microscopic). Width/height in
# pixels; if you bump font sizes or add more nodes, bump these too.
widget <- DiagrammeR::grViz(dot_spec, width = 1800, height = 2200)
htmlwidgets::saveWidget(widget, html_path,
                         selfcontained = TRUE,
                         title = "Pair-analysis cohort flowchart")
cat(sprintf("[html] → %s\n", html_path))

cat("\nView: open figure_s_pair_analysis_flowchart.html in a browser.\n",
    "For PDF: Chrome → File → Print → Save as PDF.\n", sep = "")
