# Isoform Transition Detection Methods v2.0

## Overview

This document describes a complete system for detecting and representing splicing differences between isoforms, designed to enable full reconstruction of one isoform from another.

---

## 1. Union Exon Creation

### 1.1 Inputs

**Input 1: GTF File**
- Standard GTF format (GTF2.2 or GFF3)
- Required features: `exon` records
- Required attributes:
  - `gene_id` - unique gene identifier
  - `transcript_id` - unique transcript identifier
  - `exon_number` - exon numbering (or inferred from coordinates)
- Coordinates: 1-based, inclusive (standard GTF)
- All isoforms for genes of interest must be present

**Input 2: Dominant Isoform Specification**
- Format: Tab-separated file with columns:
  - `gene_id` - matches GTF gene_id
  - `dominant_transcript_id` - the reference isoform for this gene
- One dominant isoform per gene
- All comparisons within a gene are made relative to this dominant isoform

### 1.2 Union Exon Definition

A **union exon** is a contiguous genomic region that represents exonic sequence in at least one isoform of a gene. Union exons are created by:

1. For each gene, extract all exons from all isoforms
2. Merge overlapping exons to create non-overlapping genomic intervals
3. Each merged interval is a union exon

**Properties:**
- Union exons within a gene do not overlap
- Union exons are strand-specific
- A union exon may be present in one, multiple, or all isoforms of a gene
- Union exons can be shared across different genes if genes overlap

### 1.3 Relationship Between Union Exons and Splicing Events

**Important:** Union exons are created by merging overlapping exonic regions across **all** isoforms of a gene. However, when comparing **two specific isoforms**, a single splicing event may span **multiple union exons**.

**Examples:**

*Intron Retention:*
- Isoform A: single exon 1000-2000
- Isoform B: three exons 1000-1200, 1500-1700, 1900-2000
- Union exons: ue_1 (1000-1200), ue_2 (1500-1700), ue_3 (1900-2000)
- IR event in A spans: **ue_1, ue_2, ue_3** (all three union exons)

*Alternative Splice Sites:*
- Isoform A: exon 1000-2000
- Isoform B: exon 1000-1800
- Isoform C: exon 1200-2000
- Union exons: ue_1 (1000-1200), ue_2 (1200-1800), ue_3 (1800-2000)
- A5SS between A and B involves: **ue_2, ue_3** (donor differs)

**Implication:** Event detection must record **all affected union exons** for each splicing event to enable accurate reconstruction and provide a complete description of the splicing difference.

### 1.4 Union Exon ID Assignment

Union exons are assigned globally unique sequential identifiers:
- Format: `ue_N` where N is a sequential integer
- IDs assigned consecutively during union exon creation across all genes
- IDs are stable (do not change after genomic sorting)
- IDs are independent of gene membership (shared exons have one ID)

### 1.5 Output Files

#### 1.5.1 Union Exon File

**Format:** Tab-separated file, sorted for tabix indexing. 1-based coordinates.

**Columns:**
```
chr    start    end    union_exon_id    strand
```

**Field descriptions:**
- `chr` - chromosome/contig name (matches GTF seqname)
- `start` - 1-based start coordinate (inclusive)
- `end` - 1-based end coordinate (inclusive)
- `union_exon_id` - unique identifier (ue_1, ue_2, ue_3, ...)
- `strand` - genomic strand (+ or -)

**File properties:**
- Sorted by chr (lexicographic), then start (numeric)
- Indexed with tabix for efficient genomic range queries
- Header line (starts with #) describes columns
- Union exons are globally unique by coordinates
- **All coordinates are 1-based (matches GTF convention)**

**Example:**
```
#chr    start    end      union_exon_id    strand
chr1    1000     1200     ue_1            +
chr1    1500     1700     ue_2            +
chr1    2000     2200     ue_3            +
chr2    5000     5500     ue_4            -
```

#### 1.5.2 Junction File

**Format:** Tab-separated file, sorted for tabix indexing. 1-based coordinates.

**Columns:**
```
junction_id    chr    start    end    5_prime    3_prime    strand    type
```

**Field descriptions:**
- `junction_id` - unique identifier (j_1, j_2, j_3, ...)
- `chr` - chromosome/contig name
- `start` - 1-based min(5_prime, 3_prime) (for tabix indexing)
- `end` - 1-based max(5_prime, 3_prime) (for tabix indexing)
- `5_prime` - 1-based position of 5' splice site (donor)
  - Plus strand: end coordinate of upstream exon
  - Minus strand: start coordinate of downstream exon
- `3_prime` - 1-based position of 3' splice site (acceptor)
  - Plus strand: start coordinate of downstream exon
  - Minus strand: end coordinate of upstream exon
- `strand` - genomic strand (+ or -)
  - Can be inferred from coordinates (5_prime < 3_prime = plus strand)
  - Included explicitly for clarity and backsplice junctions
- `type` - junction type (optional)
  - `linear` - canonical junction (5' → 3' in transcription direction)
  - `circular` - backsplice junction (circular RNA)
  - Empty for standard analysis

**File properties:**
- Sorted by chr, then 5_prime position
- Indexed with tabix for efficient queries
- Junctions are globally unique by coordinates
- **All coordinates are 1-based**
- Union exons connected by a junction can be inferred by matching junction coordinates to union exon boundaries

**Example:**
```
#junction_id    chr     start    end     5_prime    3_prime    strand    type
j_1            chr1    1200     1500    1200       1500       +
j_2            chr1    1700     2000    1700       2000       +
j_3            chr2    5000     5500    5500       5000       -
j_4            chr1    1000     2200    2200       1000       +         circular
```

**Note on circular junctions:** For backsplice junctions, 3_prime position comes before 5_prime in genomic coordinates, breaking the normal strand inference rule.

### 1.6 Union Exon Creation Algorithm

**For each gene:**
1. Extract all exon coordinates from all isoforms
2. Sort exons by genomic position (start coordinate)
3. Merge overlapping/adjacent exons:
   - If exons overlap by ≥1bp, merge them
   - Create single union exon spanning merged region
4. Assign next sequential union_exon_id to each merged region
5. Record: chr, start, end (0-based), union_exon_id, strand

**After processing all genes:**
1. Combine all union exons from all genes
2. Sort by chr, start (for tabix)
3. Write union exon file
4. Extract all junctions between consecutive exons
5. Assign sequential junction_ids
6. Sort junctions by chr, 5_prime
7. Write junction file
8. Index both files with tabix

**Note on coordinates:**
- All coordinates remain 1-based inclusive (GTF convention)
- No coordinate conversion needed

---

## 2. Event Detection

### 2.1 Workflow Organization

Event detection is organized by chromosome for computational efficiency:

1. **Process genes chromosome-by-chromosome**
2. **Preload union exons** for current chromosome into memory
3. **For each gene:**
   - Compare dominant isoform to comparator isoform(s)
   - Detect splicing events
   - Map events to genomic coordinates
   - Identify missing terminal exons

### 2.2 Event File Format

**Format:** Tab-separated file with detected splicing events

**Columns:**
```
gene_id                    Gene identifier
isoform_dominant           Dominant isoform ID (reference)
isoform_comparator         Comparison isoform ID
event_type                 A5SS, A3SS, Partial_IR_5, Partial_IR_3, IR, SE, Alt_TSS, Alt_TES
direction                  GAIN or LOSS
chr                        Chromosome
5_prime                    1-based position of biological 5' boundary
3_prime                    1-based position of biological 3' boundary
strand                     + or -
bp_diff                    Size difference in bp (for boundary events)
missing_terminal_exons     Comma-separated exonic coordinate ranges
```

### 2.3 Event Type Specifications

**5_prime and 3_prime usage by event type:**

| Event Type | 5_prime | 3_prime | Notes |
|------------|---------|---------|-------|
| A5SS | Donor position | Shared acceptor | Donor differs, **<100bp** |
| A3SS | Shared donor | Acceptor position | Acceptor differs, **<100bp** |
| Partial_IR_5 | Donor position | Shared acceptor | **≥100bp** at donor |
| Partial_IR_3 | Shared donor | Acceptor position | **≥100bp** at acceptor |
| IR | 5' boundary | 3' boundary | **Boundaries of IR event (not always exact intron boundaries)** |
| SE | Skipped exon start | Skipped exon end | Coordinates of missing exon |
| Alt_TSS | **5' most TSS** | **3' most TSS** | Use GAIN/LOSS + strand to infer which is dominant vs comparator |
| Alt_TES | **5' most TES** | **3' most TES** | Use GAIN/LOSS + strand to infer which is dominant vs comparator |

### 2.4 Direction Interpretation

**GAIN:** Comparator isoform has more sequence than dominant
**LOSS:** Comparator isoform has less sequence than dominant (dominant is longer)

**All event types include direction:**
- Boundary events (A5SS, A3SS, Partial_IR, IR) - which isoform extends further
- Terminal events (Alt_TSS, Alt_TES) - which isoform has longer terminal region
- SE events - marked as "-" (no direction)

### 2.5 Alt_TSS and Alt_TES Coordinate Interpretation

For terminal boundary events, the 5_prime and 3_prime positions represent the transcriptionally 5' and 3' boundaries, with one from the dominant and one from the comparator.

**Inference logic:**
- **Plus strand:** 5_prime = lower genomic coordinate, 3_prime = higher genomic coordinate
- **Minus strand:** 5_prime = higher genomic coordinate, 3_prime = lower genomic coordinate
- **LOSS direction:** Dominant extends beyond comparator → missing_terminal_exons populated
- **GAIN direction:** Comparator extends beyond dominant → missing_terminal_exons empty

**Examples:**

*Plus strand, Alt_TSS with LOSS (dominant has longer 5' end):*
```
Event: Alt_TSS  LOSS  chr1  1000  1100  +  100  missing: "1000-1100"
  → 5_prime (1000) = dominant TSS (extends further 5')
  → 3_prime (1100) = comparator TSS
  → Comparator is missing exonic region 1000-1100
```

*Minus strand, Alt_TES with LOSS (dominant has longer 3' end):*
```
Event: Alt_TES  LOSS  chr1  2100  2000  -  100  missing: "2000-2100"
  → 5_prime (2100) = comparator TES
  → 3_prime (2000) = dominant TES (extends further 3')
  → Comparator is missing exonic region 2000-2100
```

### 2.6 Missing Terminal Exons

The `missing_terminal_exons` field captures exonic regions present in the dominant isoform but absent from the comparator due to Alt_TSS or Alt_TES events.

**Format:** Comma-separated coordinate ranges
- Example: `"1000-1200,1500-1700,2000-2150"`

**Key properties:**
- Represents **exonic regions only** (not introns)
- May include **multiple non-contiguous regions** separated by introns
- Each coordinate range is a contiguous exonic segment
- Empty for GAIN direction (comparator is longer)
- Empty for non-terminal events (A5SS, A3SS, IR, SE, etc.)

**Example - Multiple missing exons:**
```
Dominant:   [==exon1==] --intron-- [==exon2==] --intron-- [==exon3==]...
            1000-1200               1500-1700               2000-2200

Comparator:                                     [====exon1====]...
                                                2150-2400

Missing terminal exons: "1000-1200,1500-1700,2000-2150"
  → Entire exon 1 (1000-1200)
  → Entire exon 2 (1500-1700)
  → Partial exon 3 (2000-2150)
  → Introns NOT included (1200-1500, 1700-2000)
```

**Reconstruction usage:**
To reconstruct the dominant isoform from the comparator, these exonic regions must be added back at the appropriate terminal position.

---

## 3. Reconstruction and Validation
*(To be documented)*
