#!/bin/bash
# Create PROVENANCE files for isoform proportion outputs

cd /Users/petecastaldi/claude_projects/nmd/tmp

# Get generation timestamp from first file
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

# Create PROVENANCE file for each isoform proportion file
for file in isoform_proportions_*_2026.1.3.tsv; do
    # Extract cell type and treatment from filename
    basename="${file%.tsv}"
    parts=(${basename//_/ })

    # Get stats from file
    n_samples=$(grep "^# Samples:" "$file" | cut -d: -f2 | tr -d ' ')
    n_genes=$(grep "^# Genes:" "$file" | cut -d: -f2 | tr -d ' ')
    n_transcripts=$(grep "^# Transcripts:" "$file" | cut -d: -f2 | tr -d ' ')

    # Create PROVENANCE filename
    prov_file="PROVENANCE_${basename}.md"

    echo "Creating $prov_file..."

    # Write PROVENANCE file
    cat > "$prov_file" << EOF
# Data Provenance: ${file}

**Generated:** ${TIMESTAMP}
**Script:** code/calculate_isoform_proportions_2026.1.3.Rmd
**Description:** Per-sample isoform proportions with expression-based filtering

## Source Data

- **Raw counts:** data/qdf_raw_counts_2026.1.3.csv
- **Phenotype:** pheno/nmd_pheno_longreadbamids_2025.1.1.csv
- **SQANTI mapping:** results/sqanti_runs/merged_collapsed/isoforms_classification_TPM.txt
- **GTF annotation:** reference_files/gencode.v49.primary_assembly.annotation.chrnamesedited.gtf

## Data Dimensions

- **Samples:** ${n_samples}
- **Genes:** ${n_genes}
- **Transcripts:** ${n_transcripts}

## Filters Applied

1. **Zero-count filter:** Removed isoforms with zero counts across all 38 samples
2. **Gene mapping filter:** Filtered to ENSG-mapped genes only (excluded novelGenes)
3. **Expression filter:**
   - Total counts > 10 per treatment/cell type group
   - Minimum 2 samples with ≥2 counts
   - Transcripts passing in EITHER DMSO OR Smg1i are retained
4. **High-count single-sample exception:** Transcripts with >10 total counts but only 1 sample with ≥2 counts saved separately

## Method

For each gene with multiple isoforms:
\`\`\`
isoform_proportion = transcript_count / gene_total_count
\`\`\`

Where:
- \`transcript_count\`: Raw read count for specific transcript in sample
- \`gene_total_count\`: Sum of all transcript counts for that gene in sample

## Output Format

Tab-separated file with columns:
- \`cell_type\`: Cell type (AT, DD, DD_ALI, DO, FB, MV)
- \`treatment\`: Treatment condition (DMSO or Smg1i)
- \`sample\`: Sample identifier
- \`gene_id\`: Gene identifier (ENSG format)
- \`transcript_id\`: Transcript identifier (PB.x.x or ENST format)
- \`raw_count\`: Raw read count for transcript in sample
- \`gene_total\`: Total read count for gene in sample
- \`isoform_proportion\`: Proportion of transcript expression (0-1)
- \`n_isoforms\`: Total number of isoforms detected for this gene

## Related Files

- **Novel isoforms:** tmp/novel_isoforms_expression_2026.1.3.tsv
- **High-count single-sample transcripts:** tmp/high_count_singlesample_transcripts_2026.1.3.tsv
- **HTML report:** tmp/calculate_isoform_proportions_2026.1.3.html

## Usage Notes

- Isoform proportions sum to 1.0 for each gene within each sample
- Zero proportions indicate transcript not expressed in that sample
- Single-isoform genes have proportion = 1.0
- Use with DGE results from longread_dge/ directory for interpretive analyses
EOF

done

echo ""
echo "Created PROVENANCE files for all isoform proportion outputs"
