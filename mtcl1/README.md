# Gene Isoform Annotation Pipeline

Extract comprehensive isoform-level annotations for genes in the NMD (Nonsense-Mediated Decay) study, combining GENCODE reference annotations with SQANTI PacBio long-read isoforms.

---

## Quick Start

```bash
# 1. Run for a gene of interest (TSV output only)
Rscript gene_isoform_annotation.R MTCL1

# 2. Include GTF and FASTA outputs
Rscript gene_isoform_annotation.R MTCL1 --gtf --fasta

# 3. Skip expression calculation (faster)
Rscript gene_isoform_annotation.R MTCL1 --no-expr

# 4. Use Ensembl ID instead of gene name
Rscript gene_isoform_annotation.R ENSG00000168502 --gtf --fasta
```

**Output:** Tab-separated file with comprehensive isoform annotations (exon structures, CDS positions, UTR lengths, expression levels, and more).

**Performance:** ~30-60 seconds per gene using tabix-indexed GTF files.

**Prerequisite:** Run `bash preprocess_sqanti_tabix.sh` once to create indexed GTF files (~10-15 minutes).

---

## Required Files

Before running the pipeline, ensure these files are accessible and paths are correctly set in the script (lines 76-83):

### Input Data Files

| File Type | Description | Typical Size | Path in Script |
|-----------|-------------|--------------|----------------|
| **GENCODE GTF (unindexed)** | Reference genome annotations (v49) | ~2.9 GB | `GENCODE_GTF_UNINDEXED` |
| **GENCODE GTF (indexed)** | Tabix-indexed GENCODE GTF | ~140 MB | `GENCODE_GTF_INDEXED` |
| **GENCODE GTF index** | Tabix index for GENCODE | ~500 KB | `<GENCODE_GTF_INDEXED>.tbi` |
| **GENCODE FASTA** | Reference transcript sequences | ~135 MB | `GENCODE_FASTA` |
| **SQANTI Classification** | PacBio isoform metadata | ~770 MB | `SQANTI_CLASSIFICATION` |
| **SQANTI GTF (indexed)** | Tabix-indexed PacBio isoforms | ~41 MB | `SQANTI_GTF_INDEXED` |
| **SQANTI GTF index** | Tabix index for SQANTI | ~224 KB | `<SQANTI_GTF_INDEXED>.tbi` |
| **SQANTI FASTA** | PacBio transcript sequences | ~497 MB | `SQANTI_FASTA` |
| **DGEList RDS** | Expression data with sample metadata | ~29 MB | `DGELIST_RDS` |

### Required Preprocessing (One-Time Setup)

**IMPORTANT:** You must run preprocessing once to create the indexed GTF files:

```bash
bash preprocess_sqanti_tabix.sh
```

This script will:
1. Create tabix-indexed version of GENCODE GTF (2.9 GB → 140 MB compressed)
2. Create tabix-indexed version of SQANTI GTF (1.0 GB → 41 MB compressed)
3. Generate tabix indices for fast coordinate lookups

**Time:** ~10-15 minutes (one-time setup)

**Output files:**
- `/path/to/gencode.v49.primary_assembly.annotation.chrnamesedited.sorted.gtf.gz` (+ `.tbi` index)
- `/path/to/sqanti3_corrected.sorted.gtf.gz` (+ `.tbi` index)

---

## Sample Output

### Example: MTCL1 Gene

**Input:**
```bash
Rscript gene_isoform_annotation.R MTCL1 --gtf --fasta
```

**Output Files:**
```
isoform_annotation_MTCL1.tsv    # 58 isoforms with 25 annotation columns
isoform_annotation_MTCL1.gtf    # 640 genomic features
isoform_annotation_MTCL1.fasta  # 58 transcript sequences
```

**TSV Columns (25 total):**
```
isoform_id          isoform_name    source      chromosome  strand
tss                 tes             cds_start   cds_end     transcript_length
n_exons             is_protein_coding           cds_length  protein_length
utr5_length         utr3_length
expr_DMSO_ddali     expr_DMSO_dd    expr_DMSO_doali
expr_DMSO_at2       expr_DMSO_fb    expr_DMSO_mv
exons               junctions       coding_source
```

**Sample Data (first isoform):**
```
ENST00000695635.1  MTCL1-216  GENCODE  18  +  8705271  8832780  8706689  8826232
127510  15  TRUE  5830  1943  984  593
0.18  0.7  NA  0.42  0  0.11
8705271_8706713_+:CDS;8718424_8718648_+:CDS;...
8706713_8718424_+;8718648_8720338_+;...
GENCODE_biotype
```

**Summary Statistics:**
- 21 GENCODE reference isoforms
- 37 SQANTI novel isoforms
- 58 total isoforms annotated
- 41 protein coding, 17 non-coding
- 8-9 isoforms with detectable expression per cell type

---

## Definitions

### Terminology

**DMSO**: Dimethyl sulfoxide, used as vehicle control in NMD inhibitor treatment. DMSO expression represents **baseline expression levels** in untreated cells (before NMD pathway inhibition).

**NMD**: Nonsense-Mediated Decay, a cellular quality control pathway that degrades aberrant transcripts. Inhibiting NMD (with Smg1i) reveals transcripts normally targeted for degradation.

**Isoform**: Alternative transcript variant of a gene, created through alternative splicing, alternative TSS/TES, or other mechanisms.

### Cell Types (Primary Human Lung Cells)

| Code | Full Name | Description |
|------|-----------|-------------|
| **DD** | Basal large airway epithelial cells | Differentiated cells from large airways |
| **DD_ALI** | Basal airway epithelial cells (ALI) | Basal airway cells cultured at air-liquid interface |
| **DO_ALI** | Basal small airway cells (ALI) | Small airway cells cultured at air-liquid interface |
| **AT2** | Alveolar epithelial type 2 cells | Surfactant-producing cells in alveoli |
| **FB** | Lung fibroblasts | Connective tissue cells |
| **MV** | Lung microvascular cells | Endothelial cells from lung microvasculature |

**ALI (Air-Liquid Interface)**: Culture condition where cells are grown at the interface between air and liquid medium, promoting differentiation and mucus production.

### Data Sources

**GENCODE**: Reference human genome annotation consortium providing high-quality gene and transcript annotations.

**SQANTI**: Software for characterizing PacBio long-read transcripts, including novel isoforms not in reference annotations.

---

## What the Pipeline Does

### Overview

The Gene Isoform Annotation Pipeline extracts and integrates comprehensive transcript-level information for a gene of interest, combining:

1. **Reference annotations** from GENCODE (canonical transcripts)
2. **Novel isoforms** from PacBio long-read sequencing (SQANTI)
3. **Structural features** (exons, CDS regions, UTRs, splice junctions)
4. **Expression data** (baseline DMSO levels across 6 lung cell types)

### Step-by-Step Process

#### Step 1: Gene Identification
- Accepts either gene name (e.g., "MTCL1") or Ensembl ID (e.g., "ENSG00000168502")
- Searches GENCODE GTF to find gene ID and confirm gene name
- Extracts gene coordinates (chromosome, start, end) for tabix queries
- Validates gene exists in reference genome

#### Step 2: GENCODE Isoform Extraction
- Queries tabix-indexed GENCODE GTF by gene coordinates (<1 second)
- Parses GTF to extract transcript features
- For each transcript:
  - **Exon coordinates** with functional types (CDS, 5'UTR, 3'UTR, mixed)
  - **Splice junctions** between consecutive exons
  - **TSS/TES** (transcription start/end sites, strand-aware)
  - **CDS boundaries** (genomic start/end of coding sequence)
  - **UTR lengths** calculated from TSS→CDS (5'UTR) and CDS→TES (3'UTR)
  - **Protein coding status** from transcript biotype and CDS presence
  - **Protein length** calculated from CDS length ÷ 3

**Key feature:** Strand-aware calculations ensure correct UTR assignment:
- **Plus strand** (+): TSS < CDS_start < CDS_end < TES
  - 5'UTR = CDS_start - TSS
  - 3'UTR = TES - CDS_end
- **Minus strand** (-): TES < CDS_start < CDS_end < TSS
  - 5'UTR = TSS - CDS_end
  - 3'UTR = CDS_start - TES

#### Step 3: SQANTI Isoform Extraction

**Tabix-indexed coordinate queries (<1 second):**
- Queries preprocessed, coordinate-indexed SQANTI GTF
- Retrieves only features overlapping gene coordinates
- Uses same coordinates as GENCODE for consistent coverage

For each SQANTI isoform:
- Extracts exon structures (simpler than GENCODE, no UTR subtypes)
- Gets coding status from SQANTI classification
- Calculates CDS positions and UTRs if available
- Identifies novel isoforms (PB.* IDs) not in reference

#### Step 4: Expression Calculation (Optional)

**Data Source:** DGEList object containing:
- TMM-normalized counts from edgeR
- Sample metadata (treatment, cell type, replicate info)
- 38 samples across 6 cell types × 2 treatments (DMSO + Smg1i)

**Process:**
1. Load DGEList and calculate normalized CPM (counts per million)
2. Filter samples: Keep only DMSO (control) samples per cell type
3. Calculate mean DMSO expression across biological replicates
4. Match isoforms by transcript ID (handles versioned IDs)
5. Round to 2 significant figures for readability

**Result:** Cell-type-specific baseline expression (6 columns)
- Values represent transcript abundance in untreated cells
- Missing values (NA) indicate isoform not detected in sequencing
- Zero values indicate detected but not expressed

#### Step 5: Output Generation

**TSV (always created):**
- 25 columns of integrated annotation data
- One row per isoform
- Expression values rounded to 2 significant figures
- Ready for downstream analysis in R/Python

**GTF (optional, `--gtf` flag):**
- Contains all genomic features for visualization
- Combines GENCODE + SQANTI features
- Compatible with genome browsers (IGV, UCSC)

**FASTA (optional, `--fasta` flag):**
- Transcript sequences for all isoforms
- Extracts from two sources:
  - GENCODE transcripts → GENCODE FASTA
  - SQANTI transcripts → SQANTI FASTA
- Handles pipe-delimited GENCODE headers automatically
- Useful for sequence analysis, primer design, structure prediction

---

## Data Flow Diagram

```
USER INPUT
    ↓
┌───────────────────────────────────────────┐
│  Gene Name (MTCL1)                        │
│  OR Ensembl ID (ENSG00000168502)          │
│  + Optional Flags (--gtf, --fasta)        │
└───────────────────────────────────────────┘
    ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 1: GENE IDENTIFICATION                                     │
│   • Search GENCODE GTF (grep gene feature)                     │
│   • Resolve gene name ↔ Ensembl ID                              │
│   • Extract gene coordinates: chr18:8705271-8832780            │
│   • Result: ENSG00000168502.19 → "MTCL1"                        │
└─────────────────────────────────────────────────────────────────┘
    ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 2: GENCODE EXTRACTION                                      │
│   Input: Tabix-indexed GENCODE GTF (140 MB compressed)         │
│   Method: tabix query by coordinates (<1 second)               │
│   Query: tabix gencode.gtf.gz chr18:8705271-8832780            │
│   ┌───────────────────────────────────────────────────────┐    │
│   │ For each transcript:                                  │    │
│   │  • Parse exons (with CDS/UTR types)                   │    │
│   │  • Calculate junctions (exon boundaries)              │    │
│   │  • Extract TSS/TES (strand-aware)                     │    │
│   │  • Get CDS start/end positions                        │    │
│   │  • Calculate UTR lengths (TSS↔CDS↔TES)               │    │
│   │  • Determine protein coding status                    │    │
│   └───────────────────────────────────────────────────────┘    │
│   Output: 21 GENCODE isoforms with full annotations            │
└─────────────────────────────────────────────────────────────────┘
    ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 3: SQANTI EXTRACTION                                       │
│   Input: Tabix-indexed SQANTI GTF (41 MB compressed)           │
│          + SQANTI Classification (770 MB)                       │
│   Method: tabix query by coordinates (<1 second)               │
│   Query: tabix sqanti.gtf.gz chr18:8705271-8832780             │
│   ┌───────────────────────────────────────────────────────┐    │
│   │ For each PacBio isoform:                              │    │
│   │  • Extract exon positions from GTF                    │    │
│   │  • Get coding status from SQANTI classification       │    │
│   │  • Calculate CDS positions (if available)             │    │
│   │  • Calculate UTRs from TSS/CDS/TES                    │    │
│   └───────────────────────────────────────────────────────┘    │
│   Output: 37 SQANTI isoforms with annotations                  │
└─────────────────────────────────────────────────────────────────┘
    ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 4: EXPRESSION CALCULATION (if --no-expr NOT set)          │
│   Input: DGEList RDS (29 MB) - 38 samples                      │
│   ┌───────────────────────────────────────────────────────┐    │
│   │ Load normalized counts (TMM/CPM)                      │    │
│   │         ↓                                              │    │
│   │ Filter: treatment=="DMSO" for each cell type          │    │
│   │         ↓                                              │    │
│   │ Calculate: mean across biological replicates          │    │
│   │         ↓                                              │    │
│   │ Match: isoform_id to DGEList rownames                 │    │
│   │         ↓                                              │    │
│   │ Round: 2 significant figures                          │    │
│   └───────────────────────────────────────────────────────┘    │
│   Cell Types:                                                   │
│     DD (n=3)  DD_ALI (n=3)  DO_ALI (n=3)                       │
│     AT2 (n=3)  FB (n=3)     MV (n=3)                           │
│   Output: 6 expression columns added to all isoforms           │
└─────────────────────────────────────────────────────────────────┘
    ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 5: OUTPUT GENERATION                                       │
│                                                                 │
│ ┌─────────────────┐  ┌──────────────────┐  ┌────────────────┐ │
│ │  TSV (always)   │  │  GTF (--gtf)     │  │ FASTA (--fasta)│ │
│ │  ─────────────  │  │  ─────────────   │  │ ────────────── │ │
│ │ 58 isoforms     │  │ 640 features     │  │ 58 sequences   │ │
│ │ 25 columns      │  │ • Genes          │  │ • GENCODE: 21  │ │
│ │ • Structures    │  │ • Transcripts    │  │ • SQANTI: 37   │ │
│ │ • CDS/UTR info  │  │ • Exons          │  │ From 2 sources │ │
│ │ • Expression    │  │ • CDS regions    │  │ Merged output  │ │
│ │ • 21 GENCODE    │  │ For genome       │  │ For sequence   │ │
│ │ • 37 SQANTI     │  │ browsers         │  │ analysis       │ │
│ └─────────────────┘  └──────────────────┘  └────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
    ↓
FINAL OUTPUT: isoform_annotation_<GENE>.tsv (+ .gtf + .fasta)
```

---

## Usage Details

### Command-Line Options

```bash
Rscript gene_isoform_annotation.R <GENE_NAME_OR_ID> [OPTIONS]
```

**Required Argument:**
- `GENE_NAME_OR_ID` - Gene symbol (e.g., `MTCL1`) or Ensembl ID (e.g., `ENSG00000168502`)

**Optional Flags:**
- `--gtf` - Generate GTF file with all genomic features
- `--fasta` - Generate FASTA file with transcript sequences
- `--no-expr` - Skip expression calculation (faster, no DGEList needed)

### Examples

```bash
# Basic usage - TSV output only
Rscript gene_isoform_annotation.R MTCL1

# Include GTF for genome browser visualization
Rscript gene_isoform_annotation.R MTCL1 --gtf

# Include FASTA for sequence analysis
Rscript gene_isoform_annotation.R MTCL1 --fasta

# Both GTF and FASTA
Rscript gene_isoform_annotation.R MTCL1 --gtf --fasta

# Structural annotations only (skip expression)
Rscript gene_isoform_annotation.R MTCL1 --no-expr

# All options combined
Rscript gene_isoform_annotation.R MTCL1 --gtf --fasta --no-expr

# Using Ensembl ID instead of gene name
Rscript gene_isoform_annotation.R ENSG00000168502

# Batch processing multiple genes
for gene in MTCL1 TP53 EGFR CFTR; do
  Rscript gene_isoform_annotation.R $gene --gtf --fasta
done
```

### Performance

**Runtime:** ~30-60 seconds per gene

The pipeline uses tabix-indexed GTF files for fast coordinate-based queries:
- GENCODE extraction: <1 second (via tabix index)
- SQANTI extraction: <1 second (via tabix index)
- Expression calculation: ~15 seconds (DGEList loading and CPM calculation)
- Total: 30-60 seconds per gene

---

## Output File Formats

### TSV File (Primary Output)

**File:** `isoform_annotation_<GENE>.tsv`

Tab-separated table with one row per isoform and 25 columns:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `isoform_id` | String | Transcript identifier | ENST00000695635.1 |
| `isoform_name` | String | Transcript name | MTCL1-216 |
| `source` | String | GENCODE or SQANTI | GENCODE |
| `chromosome` | String | Chromosome location | 18 |
| `strand` | String | Strand | + |
| `tss` | Integer | Transcription start site | 8705271 |
| `tes` | Integer | Transcription end site | 8832780 |
| `transcript_length` | Integer | Total length (bp) | 127510 |
| `n_exons` | Integer | Number of exons | 15 |
| `exons` | String | Exon coords with types | 8705271_8706713_+:CDS;... |
| `junctions` | String | Splice junction coords | 8706713_8718424_+;... |
| `cds_start` | Integer | CDS genomic start | 8706689 |
| `cds_end` | Integer | CDS genomic end | 8826232 |
| `cds_length` | Integer | CDS length (bp) | 5830 |
| `protein_length` | Integer | Protein length (aa) | 1943 |
| `utr5_length` | Integer | 5' UTR length (bp) | 984 |
| `utr3_length` | Integer | 3' UTR length (bp) | 593 |
| `is_protein_coding` | Boolean | Protein coding status | TRUE |
| `coding_source` | String | Source of coding annotation | GENCODE_biotype |
| `expr_DMSO_ddali` | Float | DMSO expr in DD_ALI (CPM) | 0.18 |
| `expr_DMSO_dd` | Float | DMSO expr in DD (CPM) | 0.7 |
| `expr_DMSO_doali` | Float | DMSO expr in DO_ALI (CPM) | NA |
| `expr_DMSO_at2` | Float | DMSO expr in AT2 (CPM) | 0.42 |
| `expr_DMSO_fb` | Float | DMSO expr in FB (CPM) | 0 |
| `expr_DMSO_mv` | Float | DMSO expr in MV (CPM) | 0.11 |

**Exon Format:** `start_end_strand:type`
- Types: `CDS`, `5UTR`, `3UTR`, `mixed`, `exon`
- Multiple exons separated by semicolons
- Example: `8705271_8706713_+:CDS;8718424_8718648_+:CDS`

**Junction Format:** `exon1_end_exon2_start_strand`
- Represents splice sites between consecutive exons
- Example: `8706713_8718424_+;8718648_8720338_+`

**Expression Values:**
- `NA` - Isoform not detected in sequencing
- `0` - Detected but not expressed in DMSO
- `>0` - Expressed; CPM values rounded to 2 sig figs

### GTF File (Optional)

**File:** `isoform_annotation_<GENE>.gtf`

Standard GTF format compatible with genome browsers:
- Gene, transcript, exon, CDS, and UTR features
- Both GENCODE and SQANTI isoforms
- Can be loaded into IGV, UCSC Genome Browser, etc.

### FASTA File (Optional)

**File:** `isoform_annotation_<GENE>.fasta`

Standard FASTA format with transcript sequences:
```
>ENST00000695635.1
ATGGCGGCGGCG...
>PB.14156.3
ATGGAGGAGCTG...
```
- Contains all isoforms (GENCODE + SQANTI)
- Headers are simple transcript IDs
- Extracted from two sources (GENCODE and SQANTI FASTA files)

---

## Installation & Dependencies

### System Requirements

- **OS**: Linux or macOS
- **RAM**: 8 GB minimum, 16 GB recommended
- **Disk**: ~5 GB for reference files

### Required Software

#### R (version ≥ 4.0)
Install from https://cran.r-project.org/

#### R Packages
```r
# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "rtracklayer",      # GTF/GFF parsing
  "GenomicRanges",    # Genomic feature handling
  "edgeR",            # DGEList handling, CPM calculation
  "Biostrings"        # FASTA file handling (only needed for --fasta)
))

# Install CRAN packages
install.packages(c(
  "dplyr",            # Data manipulation
  "tidyr",            # Data tidying
  "readr"             # Fast file reading
))
```

#### tabix and bgzip (part of htslib)

**Required for preprocessing and running the pipeline.**

**Via Conda** (recommended):
```bash
conda install -c bioconda htslib
```

**Via Homebrew** (macOS):
```bash
brew install htslib
```

**Via apt** (Ubuntu/Debian):
```bash
sudo apt-get install tabix
```

**Purpose:** Enables fast coordinate-based GTF queries (<1 second) for both GENCODE and SQANTI files.

**Verification:**
```bash
which tabix && which bgzip && echo "Tools installed!" || echo "Need to install htslib"
```

---

## Troubleshooting

### Common Issues

**Error: "Gene not found in GENCODE"**
```
Solution:
- Check gene name spelling (case-sensitive)
- Try Ensembl ID: Rscript gene_isoform_annotation.R ENSG00000168502
- Verify gene is in GENCODE v49
```

**Error: "DGEList RDS not found" or "bamid column not found"**
```
Solution:
- If you don't have expression data, use --no-expr flag
- Rscript gene_isoform_annotation.R MTCL1 --no-expr
```

**Error: "Missing required files" or "tabix index not found"**
```
Solution:
- Run preprocessing first: bash preprocess_sqanti_tabix.sh
- This creates the required tabix-indexed GTF files
- Takes 10-15 minutes but only needed once
```

**Error: "tabix command not found" or "bgzip command not found"**
```
Solution:
- Install htslib: conda install -c bioconda htslib
- OR continue without preprocessing (grep mode works fine)
```

**Warning: "N isoforms not found in DGEList"**
```
This is normal:
- Newer GENCODE versions may not be in DGE analysis
- Novel SQANTI isoforms may not be quantified
- Affected isoforms will have NA expression values
```

**Error: Package 'Biostrings' not found (with --fasta)**
```
Solution:
- Install Biostrings: BiocManager::install("Biostrings")
- OR run without --fasta flag
```

### Performance Issues

**Script taking longer than expected?**
1. Check if tabix preprocessing was done: `ls sqanti3_corrected.sorted.gtf.gz.tbi`
2. Monitor system resources: `top` or `htop`
3. Large genes (>50 isoforms) take longer
4. Expression calculation adds ~15-30 seconds - use `--no-expr` if not needed

**Out of memory errors?**
1. Close other applications
2. Use `--no-expr` to reduce memory usage
3. Process genes one at a time instead of batch mode

---

## Citation

If you use this pipeline in your research, please cite:

- **GENCODE**: Frankish et al. (2021). GENCODE 2021. *Nucleic Acids Research*, 49(D1), D916-D923.
- **SQANTI3**: Tardaguila et al. (2018). SQANTI: extensive characterization of long-read transcript sequences. *Genome Research*, 28(3), 396-411.
- **edgeR**: Robinson et al. (2010). edgeR: a Bioconductor package for differential expression analysis. *Bioinformatics*, 26(1), 139-140.

---

## Version History

- **v1.0** (2026-02-14)
  - Initial production release
  - Command-line interface with gene name or Ensembl ID
  - Optional GTF and FASTA outputs via flags
  - Tabix optimization with automatic fallback
  - DMSO-specific expression from DGEList
  - CDS positions and strand-aware UTR calculations
  - Expression values rounded to 2 significant figures
  - Comprehensive error checking and diagnostics

---

## Support

For questions or issues:
1. Check the Troubleshooting section above
2. Verify file paths in configuration (lines 66-72 of the script)
3. Ensure all required files are accessible
4. Review diagnostic messages in output

---

## License

This pipeline is part of the NMD research project. Please cite appropriately if used in publications.
