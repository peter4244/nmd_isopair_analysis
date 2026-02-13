#!/bin/bash
# Master pipeline script for detailed event analysis
# Runs all phases of Version 2.0.0

set -e  # Exit on error

RESULTS_DIR="/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions"
CODE_DIR="/Users/petecastaldi/claude_projects/nmd/code"

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║    DETAILED EVENT ANALYSIS PIPELINE - VERSION 2.0.0            ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# Check if detailed event detection is complete
if [ ! -f "$RESULTS_DIR/detailed_event_vectors.rds" ]; then
    echo "ERROR: Detailed event vectors not found!"
    echo "Please run detailed_event_detection.R first"
    exit 1
fi

echo "✓ Detailed event vectors found"
echo ""

# Phase 5: Q1-Q3 Analysis
echo "═══ Phase 5: Running Q1-Q3 Analysis ═══"
Rscript "$CODE_DIR/detailed_event_analysis_q1_q3.R"
echo "✓ Q1-Q3 analysis complete"
echo ""

# Phase 6: Visualizations
echo "═══ Phase 6: Generating Visualizations ═══"
Rscript "$CODE_DIR/detailed_event_visualizations.R"
echo "✓ Visualizations complete"
echo ""

# Summary
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║         DETAILED ANALYSIS PIPELINE COMPLETE                    ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Results saved to: $RESULTS_DIR"
echo ""
echo "Generated files:"
echo "  - q1_event_incidence_by_celltype.tsv"
echo "  - q2_event_cooccurrence_matrix.tsv"
echo "  - q2_positively_correlated_pairs.tsv"
echo "  - q2_negatively_correlated_pairs.tsv"
echo "  - q3_positive_pair_distances.tsv (if applicable)"
echo "  - q3_distance_summary.tsv (if applicable)"
echo "  - fig1_detailed_event_incidence.pdf"
echo "  - fig2_event_proportions.pdf"
echo "  - fig3_complexity_distribution.pdf"
echo "  - fig4_cooccurrence_heatmap.pdf"
echo "  - fig5_mean_events_by_celltype.pdf"
echo ""
