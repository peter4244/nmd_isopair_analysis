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
*(To be documented)*

## 3. Reconstruction and Validation
*(To be documented)*
