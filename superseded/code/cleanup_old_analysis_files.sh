#!/usr/bin/bash
# Remove old junction-based analysis files and test files
# Keep only union model approach files

DIR="results/isoform_transitions/v3.0_reference_based"

echo "=== Cleaning up old analysis files ==="
echo

# Create backup directory
mkdir -p ${DIR}/archive_old_junction_based
echo "Created archive directory"

# Move old junction-based analysis files
echo "Archiving old junction-based analysis files..."
mv ${DIR}/reference_event_vectors_v3.0_all.rds ${DIR}/archive_old_junction_based/
mv ${DIR}/reference_event_vectors_v3.0_filtered.rds ${DIR}/archive_old_junction_based/
mv ${DIR}/reference_event_vectors_v3.0_all_BACKUP.rds ${DIR}/archive_old_junction_based/
mv ${DIR}/reference_event_vectors_v3.0_filtered_BACKUP.rds ${DIR}/archive_old_junction_based/
mv ${DIR}/reference_cooccurrence_*.tsv ${DIR}/archive_old_junction_based/
mv ${DIR}/reference_event_summary_*.tsv ${DIR}/archive_old_junction_based/
mv ${DIR}/reference_coordinated_*.tsv ${DIR}/archive_old_junction_based/
mv ${DIR}/v3.0_topological_*.tsv ${DIR}/archive_old_junction_based/
mv ${DIR}/reference_isoforms_dmso.rds ${DIR}/archive_old_junction_based/
mv ${DIR}/reference_isoforms_dmso.tsv ${DIR}/archive_old_junction_based/

echo "Archived old junction-based files"

# Remove test files
echo "Removing test files..."
rm -f ${DIR}/union_exon_models_test.rds
rm -f ${DIR}/union_model_events_test.rds
rm -f ${DIR}/filtered_genes_test.rds
rm -f ${DIR}/mxe_events_test.rds
rm -f ${DIR}/mxe_events_test.tsv

echo "Removed test files"

echo
echo "=== Cleanup complete ==="
echo
echo "Archived files in: ${DIR}/archive_old_junction_based/"
echo
echo "Current files (union model approach):"
ls -lh ${DIR}/ | grep -v "^d" | grep -v "^total"
