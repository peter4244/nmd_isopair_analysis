#!/usr/bin/bash
# Archive obsolete files from v3.0 junction-based approach
# Keep only v4.0 union model files

echo "=== Cleaning up transition to v4.0 Union Model approach ==="
echo

# Create archive directories
ARCHIVE_CODE="code/archive_v3.0_junction_based"
ARCHIVE_RESULTS="results/isoform_transitions/v3.0_reference_based/archive_v3.0_documentation"

mkdir -p ${ARCHIVE_CODE}
mkdir -p ${ARCHIVE_RESULTS}

echo "Created archive directories"
echo

# ============================================================================
# Archive OLD scripts (junction-based approach)
# ============================================================================

echo "Archiving old junction-based scripts..."

# Old event detection scripts (junction-based)
mv code/detailed_event_detection.R ${ARCHIVE_CODE}/ 2>/dev/null
mv code/detailed_event_detection_v2.1.R ${ARCHIVE_CODE}/ 2>/dev/null
mv code/detailed_event_detection_v2.2_directional.R ${ARCHIVE_CODE}/ 2>/dev/null
mv code/detailed_event_detection_v3.0_reference_based.R ${ARCHIVE_CODE}/ 2>/dev/null

# Old exon structure scripts
mv code/regenerate_exon_structures.R ${ARCHIVE_CODE}/ 2>/dev/null

# Test/development scripts (keep for reference but archive)
mv code/build_union_exon_model.R ${ARCHIVE_CODE}/ 2>/dev/null
mv code/build_union_exon_model_test_suox.R ${ARCHIVE_CODE}/ 2>/dev/null
mv code/detect_events_from_union_model.R ${ARCHIVE_CODE}/ 2>/dev/null

# Failed version (replaced by batched)
mv code/build_union_exon_model_full.R ${ARCHIVE_CODE}/ 2>/dev/null

echo "  Archived junction-based and test scripts"

# ============================================================================
# Archive OLD documentation
# ============================================================================

echo "Archiving old documentation..."

cd results/isoform_transitions/v3.0_reference_based

mv EVENT_VECTOR_FIX_SUMMARY.md ${ARCHIVE_RESULTS}/ 2>/dev/null
mv REFERENCE_BASED_ANALYSIS_V3.0_SUMMARY.md ${ARCHIVE_RESULTS}/ 2>/dev/null
mv v3.0_detection.log ${ARCHIVE_RESULTS}/ 2>/dev/null

cd ../../..

echo "  Archived old documentation files"

# ============================================================================
# Clean up temporary log files
# ============================================================================

echo "Cleaning up temporary files..."

rm -f /tmp/union_model_full.log
rm -f /tmp/union_model_batched.log

echo "  Removed temporary log files"

# ============================================================================
# Summary
# ============================================================================

echo
echo "=== Cleanup Complete ==="
echo
echo "Archived directories:"
echo "  - ${ARCHIVE_CODE}/"
echo "  - ${ARCHIVE_RESULTS}/"
echo
echo "Current v4.0 scripts (code/):"
echo "  - filter_genes_for_union_model.R"
echo "  - regenerate_exon_structures_for_union_model.R"
echo "  - build_union_exon_model_full_batched.R"
echo "  - detect_events_from_union_model_full.R (to be created)"
echo "  - detect_mxe_from_union_model.R"
echo
echo "Current v4.0 data files (results/isoform_transitions/v3.0_reference_based/):"
echo "  - genes_for_union_model.rds/tsv"
echo "  - isoforms_for_union_model.rds/tsv"
echo "  - exon_structures_by_isoform_full.rds"
echo "  - transitions_for_union_model_unfiltered.rds"
echo "  - union_exon_models_full.rds"
echo "  - union_model_statistics.tsv"
echo "  - filtered_genes_full.rds/tsv"
echo
echo "Next step: Run event detection"
echo "  Rscript code/detect_events_from_union_model_full.R"
echo
