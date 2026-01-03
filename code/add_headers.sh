#!/bin/bash
# Add column headers to isoform proportion files

cd /Users/petecastaldi/claude_projects/nmd/tmp

# Column header line
HEADER="cell_type\ttreatment\tsample\tgene_id\ttranscript_id\traw_count\tgene_total\tisoform_proportion\tn_isoforms"

# Process each isoform proportion file
for file in isoform_proportions_*_2026.1.3.tsv; do
    echo "Processing $file..."

    # Create temporary file with header inserted after comments
    awk -v header="$HEADER" '
        /^#/ { print; next }
        !printed { print header; printed=1 }
        { print }
    ' "$file" > "${file}.tmp"

    # Replace original with updated file
    mv "${file}.tmp" "$file"

    echo "  Added header to $file"
done

echo ""
echo "All files updated with column headers"
