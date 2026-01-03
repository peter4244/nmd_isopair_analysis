# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains differential gene expression (DGE) analysis results for an NMD (Nonsense-Mediated Decay) study using lung cell lines. The project includes both long-read (PacBio) and short-read RNA sequencing data across multiple experimental conditions.

## Data Architecture

### Directory Structure

- **pheno/** - Phenotype metadata files linking sample IDs to BAM files, treatments, and experimental conditions
- **longread_dge/** - Long-read DGE results (PacBio IsoSeq data) with transcript-level resolution
- **shortread_dge/** - Short-read DGE results with gene-level quantification
- **code/** - Analysis scripts (currently empty)
- **results/** - Output directory for analysis results

### Experimental Design

**Cell Types:**
- DD (Differentiated Cells)
- DD_ALI (Differentiated Cells - Air Liquid Interface)
- DO (Donor Cells)
- AT2 (Alveolar Type 2 cells)
- FB (Fibroblasts)
- MV (Microvascular cells)

**Treatments:**
- Smg1i (SMG1 inhibitor - blocks NMD pathway)
- DMSO (vehicle control)

### Data File Naming Convention

**Long-read DGE:** `nmd_dge_{celltype}_YYYY.M.D.csv`
**Short-read DGE:** `nmd_sr_dge_{celltype}_YYYY.M.D.csv`
**Phenotype:** `nmd_pheno_longreadbamids_YYYY.M.D.csv`

Cell type codes: `at2`, `dd`, `ddali`, `do`, `doali`, `fb`, `mv`

## Data Formats

### Phenotype Files (pheno/)

Columns:
- `lib.size` - Library size for normalization
- `norm.factors` - Normalization factors
- `id` - Sample donor ID
- `sample` - Full sample identifier (format: {CELLTYPE}{DONOR_ID}_{TREATMENT})
- `treatment` - Smg1i or DMSO
- `ct` - Cell type
- `sample_num` - Sample number
- `bam` - Full path to BAM file on cluster

### Long-read DGE Files (longread_dge/)

PacBio IsoSeq transcript-level analysis with columns:
- `txid` - Transcript ID (Ensembl or PacBio IsoSeq ID like PB.xxxxx.xxx)
- `transcript.gc` - Transcript name from GENCODE
- `hgnc_id.gc` - Gene ID from GENCODE
- `biotype` - Transcript biotype (protein_coding, lncRNA, processed_pseudogene, etc.)
- `transcript.sq` - Transcript ID from SQANTI
- `hgnc_id.sq` - Gene ID from SQANTI
- `logFC` - Log2 fold change (Smg1i vs DMSO)
- `CI.L`, `CI.R` - Confidence interval bounds
- `AveExpr` - Average expression
- `t` - t-statistic
- `P.Value` - Raw p-value
- `adj.P.Val` - Adjusted p-value (FDR)
- `B` - B-statistic (log-odds)
- `hgnc_id` - Final gene identifier

**Important:** Long-read data includes novel isoforms (indicated by "novel" in transcript.sq or PB. prefix) not in reference annotations.

### Short-read DGE Files (shortread_dge/)

Gene-level analysis with columns:
- `ensembl_gene_id_version` - Ensembl gene ID with version
- `hgnc_symbol` - HGNC gene symbol
- `external_gene_name` - Gene name
- `logFC` - Log2 fold change (Smg1i vs DMSO)
- `CI.L`, `CI.R` - Confidence interval bounds
- `AveExpr` - Average expression
- `t` - t-statistic
- `P.Value` - Raw p-value
- `adj.P.Val` - Adjusted p-value (FDR)
- `B` - B-statistic

## Key Analysis Considerations

### Transcript vs Gene Level
- Long-read data provides **transcript-level** resolution, capturing isoform diversity and novel transcripts
- Short-read data provides **gene-level** resolution with better quantification accuracy
- When comparing datasets, aggregate long-read transcripts to gene level or consider isoform-specific effects

### NMD Pathway Effects
- Positive logFC values indicate **upregulation with Smg1i** (transcripts/genes rescued from NMD)
- Negative logFC values indicate **downregulation with Smg1i** (possibly indirect effects)
- NMD targets typically show strong upregulation when SMG1 is inhibited

### Statistical Thresholds
- Conventional thresholds: `|logFC| > 1` and `adj.P.Val < 0.05`
- Highly significant results have `adj.P.Val < 1e-6` and `B > 10`

### Novel Isoforms (Long-read Only)
- Transcripts with `PB.xxxxx` IDs are from PacBio IsoSeq assembly
- "novel" annotations indicate isoforms not in GENCODE reference
- These may represent NMD-specific isoforms or cell-type-specific transcripts

## Working with the Data

### Loading Data
```python
import pandas as pd

# Load phenotype data
pheno = pd.read_csv('nmd/pheno/nmd_pheno_longreadbamids_2025.1.1.csv')

# Load long-read DGE for a specific cell type
lr_dge = pd.read_csv('nmd/longread_dge/nmd_dge_dd_2026.1.2.csv')

# Load short-read DGE
sr_dge = pd.read_csv('nmd/shortread_dge/nmd_sr_dge_dd_2026.1.2.csv')
```

```r
# Load in R
pheno <- read.csv('nmd/pheno/nmd_pheno_longreadbamids_2025.1.1.csv')
lr_dge <- read.csv('nmd/longread_dge/nmd_dge_dd_2026.1.2.csv')
sr_dge <- read.csv('nmd/shortread_dge/nmd_sr_dge_dd_2026.1.2.csv')
```

### Common Analysis Tasks

**Filter for significant genes:**
```python
sig = lr_dge[(lr_dge['adj.P.Val'] < 0.05) & (abs(lr_dge['logFC']) > 1)]
```

**Identify NMD targets (upregulated with Smg1i):**
```python
nmd_targets = lr_dge[(lr_dge['logFC'] > 1) & (lr_dge['adj.P.Val'] < 0.05)]
```

**Find novel isoforms:**
```python
novel = lr_dge[lr_dge['transcript.sq'] == 'novel']
```

**Aggregate transcripts to gene level:**
```python
gene_level = lr_dge.groupby('hgnc_id').agg({
    'logFC': 'mean',
    'adj.P.Val': 'min',
    'txid': 'count'  # number of isoforms per gene
})
```

### Cross-dataset Comparisons

When comparing long-read and short-read results:
1. Aggregate long-read data to gene level using `hgnc_id`
2. Match short-read data using `hgnc_symbol` or `ensembl_gene_id_version`
3. Be aware that long-read data may have lower statistical power but better isoform resolution

## Data Provenance

- Long-read data: PacBio IsoSeq sequencing on Randell Lung Cells
- BAM files located at: `/proj/regeps/regep00/studies/ExternalCellLines/data/longread/mrna/Randell_Lung_Cells_2025/`
- Sequencing runs from July-August 2025
- DGE analysis performed using limma-voom pipeline

## Analysis Workflows

### Isoform Proportion Analysis

**Purpose:** Calculate per-sample isoform proportions to understand transcript diversity and NMD effects

**Scripts:**
- `code/calculate_isoform_proportions_2026.1.3.Rmd` - Calculate isoform proportions with expression-based filtering
- `code/interpret_isoform_patterns_2026.1.3.Rmd` - Interpretive analyses (Q1-Q4)

**Workflow:**

1. **Data Processing:**
   - Load raw transcript counts from `data/qdf_raw_counts_2026.1.3.csv`
   - Map transcripts to genes using SQANTI and GTF annotations
   - Filter to ENSG-mapped genes only (exclude novelGenes)
   - Calculate proportions: `isoform_proportion = transcript_count / gene_total`

2. **Expression Filtering:**
   - Remove zero-count isoforms across all samples
   - Apply expression filter: total_counts > 10 AND ≥2 samples with ≥2 counts per treatment/cell type
   - Transcripts passing in EITHER DMSO OR Smg1i are retained
   - High-count single-sample transcripts (>10 counts, 1 sample) saved separately

3. **Output Files (in tmp/):**
   - **Isoform proportions:** `isoform_proportions_{celltype}_{treatment}_2026.1.3.tsv` (12 files: 6 cell types × 2 treatments)
   - **Novel isoforms:** `novel_isoforms_expression_2026.1.3.tsv` (145K novel transcripts with expression data)
   - **High-count single-sample:** `high_count_singlesample_transcripts_2026.1.3.tsv` (214K transcripts)
   - **HTML reports:** `calculate_isoform_proportions_2026.1.3.html` and `interpret_isoform_patterns_2026.1.3.html`
   - **PROVENANCE files:** `PROVENANCE_isoform_proportions_{celltype}_{treatment}_2026.1.3.md` (12 files)

4. **Interpretive Analyses:**
   - **Q1:** Isoform diversity - NMD inhibition increases detected isoforms per gene
   - **Q2:** Dominant isoforms - Percentage of genes with >50% expression from single isoform
   - **Q3:** NMD target baseline - Upregulated transcripts were predominantly minor isoforms in DMSO
   - **Q4:** Cell-type specificity - Patterns are largely consistent across cell types

**Data Dimensions:**
- ~19,600 genes per cell type
- ~200,000-290,000 transcripts per cell type (after filtering)
- 38 samples total (19 DMSO + 19 Smg1i across 6 cell types)

**Isoform Proportion File Format:**

Tab-separated with header comments and column names:
```
cell_type       treatment       sample          gene_id         transcript_id   raw_count       gene_total      isoform_proportion      n_isoforms
DD              DMSO            DD010R_DMSO     ENSG00000227232 PB.1.1          15.2            78.6            0.193               23
```

**Loading Isoform Proportions:**
```r
library(data.table)

# Load single file
dd_dmso <- fread("tmp/isoform_proportions_dd_dmso_2026.1.3.tsv",
                 skip = 20, data.table = FALSE)  # Skip header + column names

# Load all files
files <- list.files("tmp", pattern = "^isoform_proportions_.*\\.tsv$", full.names = TRUE)
all_props <- map_df(files, ~fread(.x, skip = 20, data.table = FALSE))
```

### Directory Organization

**tmp/** - Working directory for initial analyses
- Isoform proportion files
- Interpretive analysis results
- HTML reports
- PROVENANCE documentation

**final/** - Paper-ready analyses (after review/approval)
- Finalized datasets
- Publication-ready figures
- Supplementary data files

**Workflow:** Analyses are developed in `tmp/`, reviewed, and moved to `final/` when approved for publication.
