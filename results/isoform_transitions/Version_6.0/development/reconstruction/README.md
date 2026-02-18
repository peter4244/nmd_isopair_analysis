# Event Detection Workflow

## Purpose

The Event Detection workflow identifies and validates splicing events by comparing isoform pairs within genes. It uses a reconstruction-based validation approach: if we can correctly reconstruct the dominant isoform from the comparator isoform plus detected events, we know the event detection is correct.

This workflow is designed as a modular, standalone component that can integrate into larger isoform analysis pipelines without modification.

## Workflow Overview

```
Input GTF + Dominant Mapping
         ↓
Generate Isoform Pairs (dominant vs non-dominant)
         ↓
Create Atomic Union Exons
         ↓
Detect Splicing Events
         ↓
Validate via Reconstruction
```

## Input Data Requirements

### 1. GTF File

**Format**: GTF2 format with exon features

**Structure**: GTF files have 9 tab-separated columns. The 9th column (called "Attributes") contains semicolon-separated key-value pairs.

**Required Attributes** (within column 9):
- `gene_id`: Gene identifier (quoted string)
- `transcript_id`: Transcript/isoform identifier (quoted string)

**Note on Exon Numbering**: The `exon_number` attribute is not required. We compute it internally. The exon numbers in the GTF are validated against our own computed numbers.

**Required Columns** (all 9 tab-separated fields):
1. Chromosome/scaffold name
2. Source (e.g., "GENCODE", "PacBio", "test")
3. Feature type (must include "exon")
4. Start position (1-based, inclusive)
5. End position (1-based, inclusive)
6. Score (can be ".")
7. Strand (+ or -)
8. Frame (can be ".")
9. Attributes (semicolon-separated key-value pairs)

**Example**:
```
chr1    GENCODE    exon    1000    1200    .    +    .    gene_id "GENE1"; transcript_id "GENE1.1"; exon_number "1";
chr1    GENCODE    exon    1500    1700    .    +    .    gene_id "GENE1"; transcript_id "GENE1.1"; exon_number "2";
chr1    GENCODE    exon    1000    1200    .    +    .    gene_id "GENE1"; transcript_id "GENE1.2"; exon_number "1";
```

**Assumptions**:
- Coordinates are 1-based, inclusive (standard GTF format)
- All exons for a given transcript_id belong to the same gene_id
- All exons for a given gene_id are on the same chromosome and strand
- Exon coordinates are non-overlapping within a transcript

**Notes**:
- Comment lines starting with `#` are ignored
- GTF files with or without explicit transcript/gene features are supported through initial processing

### 2. Dominant Isoform Mapping

**Format**: Tab-separated values (TSV) file

**Required Columns**:
1. `gene_id`: Gene identifier (must match GTF gene_ids)
2. `dominant_isoform_id`: Transcript identifier designated as dominant (must match GTF transcript_ids)

**Example**:
```
gene_id             dominant_isoform_id
GENE1               GENE1.1
GENE2               GENE2.3
```

**Assumptions**:
- Exactly one dominant isoform per gene
- All gene_ids in the mapping exist in the GTF
- All dominant_isoform_ids exist in the GTF
- Genes with only one isoform in the GTF will not generate pairs (skipped with warning)

**Notes**:
- The definition of "dominant" is external to this workflow
- For the NMD project: dominant isoforms are defined by Differential Isoform Expression (DIE) analysis
- For testing: dominant isoforms may be randomly assigned

## Output Data

### 1. Isoform Pairs File
**Format**: TSV with columns: gene_id, isoform_A (dominant), isoform_B (comparator), chr, strand

### 2. Atomic Union Exons File
**Format**: Bgzip-compressed, tabix-indexed TSV
**Columns**: chr, start, end, union_exon_id, strand, gene_id
**Files**: `*.tsv.gz` (bgzipped data), `*.tsv.gz.tbi` (tabix index)
**Note**: Genes can overlap genomically. The same coordinates may appear in multiple union exons with different gene_ids. Union exon IDs are unique within a gene but not globally.

### 3. Detected Events File
**Format**: RDS file containing `event_results` list with elements:

- **events**: Data frame with columns:
  - gene_id
  - dominant_transcript_id
  - comparator_transcript_id
  - event_type
  - direction (GAIN/LOSS)
  - chr
  - strand
  - event_coordinates (comma-separated genomic positions corresponding to union exon boundaries)

- **TSS_tolerance**: Numeric value (bp tolerance for transcription start site differences)
- **TES_tolerance**: Numeric value (bp tolerance for transcription end site differences)
- **source_gtfs**: Character vector of input GTF file path(s)

### 4. Validation Results File
**Format**: TSV with columns: gene_id, dominant_isoform, comparator_isoform, status (PASS/FAIL), reason

## Creation of Atomic Union Exons

**Purpose**: Atomic union exons provide a common coordinate system for comparing isoforms within a gene. By splitting all exons at every boundary point, we create non-overlapping genomic segments that serve as the basis for event detection and reconstruction.

**Algorithm**:
1. For each gene independently:
   - Collect all exon boundaries (start and end positions) for all isoforms in the gene
   - Sort boundaries and identify unique positions
   - Create segments by splitting the genomic space between consecutive boundaries
   - Filter to segments that are actually covered by at least one exon
2. Assign global union exon IDs across all genes

**Output**: Each atomic union exon represents a contiguous genomic region that is not split by any exon boundary within its gene. These segments are indexed with tabix for efficient coordinate-based queries.

**Overlapping Genes**: When genes overlap genomically, the same genomic coordinates may appear in multiple union exons with different gene_ids and different union_exon_ids. Each gene's union exons are computed independently based on its own isoforms' boundaries.

**Script**: `build_atomic_union_exons.R`

**Key Property**: Union exon boundaries correspond exactly to splice sites and terminal boundaries (TSS/TES), ensuring all event coordinates can be expressed as union exon boundaries.

## Event Detection

**Core Event Types:**

- **SE** (Skipped Exon): Exon present in one isoform but absent in the other, with flanking exons shared between isoforms
  - Function: `detect_se()`

- **A5SS** (Alternative 5' Splice Site): Different 5' donor boundary, <100bp difference
  - Function: `detect_shared_boundary_event()` (exact or overlap-based)

- **A3SS** (Alternative 3' Splice Site): Different 3' acceptor boundary, <100bp difference
  - Function: `detect_shared_boundary_event()` (exact or overlap-based)

- **Alt_TSS** (Alternative Transcription Start Site): Different transcription start site positions (first exon) with >20bp difference
  - Function: `detect_tss_change()`

- **Alt_TES** (Alternative Transcription End Site): Different transcription end site positions (last exon) with >20bp difference
  - Function: `detect_tes_change()`

- **IR** (Intron Retention): Exon in one isoform spans multiple (≥2) consecutive exons in the other isoform
  - Function: `detect_ir_simple()` and `detect_shared_boundary_event()` (for boundary-spanning cases)

- **Partial_IR** (Partial Intron Retention): One boundary similar, other differs by ≥100bp
  - Subtypes: Partial_IR_5 (5' donor differs), Partial_IR_3 (3' acceptor differs)
  - Function: `detect_shared_boundary_event()`

**Detection Modes:**

1. **Exact Boundary Match**: Events are detected when one exon boundary matches exactly between isoforms while the other differs. Used as the primary detection method.

2. **Overlap-Based Detection**: Fallback method when no exact boundary match exists but exons overlap. This relaxed approach:
   - Checks if exons overlap by ≥1bp
   - Compares 5' and 3' boundaries in a strand-aware manner
   - Determines if boundary differences fall within an exon or extend beyond into flanking regions
   - Classifies based on distance: <100bp → splice site variation (A5SS/A3SS), ≥100bp → retention (Partial_IR/IR)
   - Implemented in: `check_boundary_within_exon()` helper function
   - Skips monoexonic vs multi-exonic comparisons (handled by IR detection)

**Event Detection Flow:**

```
For each isoform pair (dominant vs comparator):

STEP 1: Terminal Boundary Detection
├─ Compare first exons → detect_tss_change()
│  └─ If TSS differs by >20bp → Alt_TSS event
└─ Compare last exons → detect_tes_change()
   └─ If TES differs by >20bp → Alt_TES event

STEP 2: Intron Retention Detection
├─ Check each comparator exon vs all dominant exons → detect_ir_simple()
│  └─ If overlaps (≥1bp) with ≥2 dominant exons → IR event (GAIN direction)
├─ Check each dominant exon vs all comparator exons → detect_ir_simple()
│  └─ If overlaps (≥1bp) with ≥2 comparator exons → IR event (LOSS direction)
└─ Track pairwise (comparator, dominant) exon combinations involved in IR
   (Each IR event involves ≥3 exons: 1 spanning, ≥2 spanned)

STEP 3: Boundary Event Detection (for overlapping exon pairs not in IR)
└─ For each (comparator, dominant) exon pair with genomic overlap (≥1bp):
   └─ Skip if this pair is involved in an IR event (from Step 2)
   └─ Call detect_shared_boundary_event()
      │
      ├─ MODE A: Exact Boundary Match (tried first)
      │  ├─ Check if acceptor coordinates exactly equal
      │  ├─ Check if donor coordinates exactly equal
      │  └─ If one shared, one differs:
      │     ├─ Difference <100bp → A5SS or A3SS
      │     └─ Difference ≥100bp → Partial_IR_5 or Partial_IR_3
      │
      └─ MODE B: Overlap-Based Detection (fallback if Mode A finds nothing)
         ├─ Verify exons overlap by ≥1bp
         ├─ Compare boundaries in strand-aware manner
         └─ Call check_boundary_within_exon() for each boundary:
            ├─ Case A: Comparator boundary within dominant exon
            │  ├─ <100bp diff → A5SS/A3SS (LOSS direction)
            │  └─ ≥100bp diff → Partial_IR (LOSS direction)
            └─ Case B: Comparator boundary extends beyond dominant
               ├─ Check if extends into flanking exon
               │  ├─ <100bp → A5SS/A3SS (GAIN direction)
               │  └─ ≥100bp → IR (GAIN direction)
               └─ Doesn't reach flanking
                  ├─ <100bp → A5SS/A3SS (GAIN direction)
                  └─ ≥100bp → Partial_IR (GAIN direction)
```

**Technical Note: IR Exon Pair Tracking**

In Step 2, IR events are identified when one exon overlaps with multiple exons. Each IR event involves at least 3 exons total: 1 spanning exon and ≥2 spanned exons. To prevent these exons from being re-analyzed in Step 3 (boundary detection), all pairwise combinations of (comparator_exon, dominant_exon) involved in the IR are tracked.

**Why pairs?** Each tracked pair consists of exactly one comparator exon index and one dominant exon index, both involved in the same IR event. In Step 3, before performing boundary analysis on any (comparator, dominant) exon pair, the code checks if that pair is in the IR tracking list and skips it if present.

**Order consistency:** The tracking always stores pairs in `(comparator_index, dominant_index)` format, regardless of IR direction:
- GAIN direction (comparator spans dominant): stores `(comp=i, dom=j)` for each overlapping dominant exon j
- LOSS direction (dominant spans comparator): stores `(comp=j, dom=i)` for each overlapping comparator exon j

This consistent ordering ensures Step 3 correctly identifies and skips all exon combinations already classified as IR, preventing redundant or conflicting event assignments.

**Script**: `detect_and_save_events.R` (sources `scripts/event_detection_functions.R`)

## Validation with Reconstruction

**Purpose**: Reconstruction serves as validation of event detection correctness. If we can accurately reconstruct the dominant isoform from the comparator isoform plus the detected events, this confirms that our event detection captured all meaningful differences between the isoforms.

**Reconstruction Algorithm**:
1. Start with the comparator isoform exon structure
2. For each detected event (in order):
   - **LOSS events**: Add regions to comparator (dominant has MORE)
     - Alt_TSS/Alt_TES: Extend terminal exons using missing_terminal_exons coordinates
     - A5SS/A3SS/Partial_IR: Add sequence at splice boundaries
     - SE: Add skipped exons
     - IR: Split retained introns
   - **GAIN events**: Remove regions from comparator (dominant has LESS)
     - Alt_TSS/Alt_TES: Trim terminal exons
     - A5SS/A3SS/Partial_IR: Remove sequence at splice boundaries
     - SE: Remove extra exons
3. Match reconstructed exons to union exon boundaries
4. Merge adjacent segments: Combine adjacent or overlapping exon segments into continuous exons. This ensures that when sequence is added to extend an exon (e.g., Alt_TSS extending the first exon), the extension and original exon are joined into a single exon rather than treated as separate exons.
   - Function: `merge_adjacent_exons()`
5. Compare reconstructed structure to expected dominant isoform

**Validation Criteria**:
- **Exact match**: Reconstructed exon structure matches expected structure exactly
- **TSS tolerance**: Transcription start site differences within configured tolerance are acceptable (reported in reason field)
- **TES tolerance**: Transcription end site differences within configured tolerance are acceptable (reported in reason field)

**Scripts**:
- `reconstruction_functions.R` - Core reconstruction logic (`reconstruct_dominant_v2()`, `merge_adjacent_exons()`)
- `reconstruct_dominant_isoforms.R` - Main reconstruction pipeline
- `test_reconstruction_v2.R` - Validation testing script that compares reconstructed isoforms to expected structures and reports pass/fail status

## Technical Notes

- **GAIN/LOSS Semantics**: Event directions (GAIN/LOSS) are defined from the comparator's perspective:
  - LOSS = comparator LOST sequence (dominant has more exonic sequence)
  - GAIN = comparator GAINED sequence (dominant has less exonic sequence)
- **Union Exon Approach**: Uses atomic (non-overlapping) union exon segments
- **Strand-Aware**: All event detection and reconstruction is strand-aware
