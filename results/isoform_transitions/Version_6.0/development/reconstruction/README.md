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
**Format**: Tab-separated values (TSV)

**Columns**:
- `gene_id`
- `dominant_transcript_id`
- `comparator_transcript_id`
- `event_type`: SE, A5SS, A3SS, Partial_IR_5, Partial_IR_3, IR, Alt_TSS, Alt_TES
- `direction`: GAIN or LOSS (from comparator's perspective — see GAIN/LOSS Semantics below)
- `chr`
- `five_prime`: 5'-most coordinate of the event (genomic; not strand-corrected)
- `three_prime`: 3'-most coordinate of the event (genomic; not strand-corrected)
- `strand`
- `bp_diff`: size of the event in bp (NA for IR)
- `missing_terminal_exons`: for Alt_TSS/Alt_TES, comma-separated coordinate ranges of the
  longer isoform's terminal extension (e.g., `"1000-1200,1500-1700"`); empty for other event types
- `n_terminal_regions`: for Alt_TSS/Alt_TES, count of coordinate ranges in missing_terminal_exons;
  1 = simple boundary shift on a shared exon, >1 = multiple exons or partial+whole exons differ
- `missing_internal_exons`: comma-separated coordinate ranges of internal exons present in one
  isoform but absent in the other (within the overlap span); populated for both GAIN and LOSS;
  diagnostic only — not consumed by reconstruction
- `ir_split_exons`: for IR events, comma-separated coordinate ranges of the spanned exons

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

- **SE** (Skipped Exon): Exon present in one isoform but absent in the other, with both flanking
  exons (upstream and downstream) present in the other isoform (strict flanking requirement)

- **A5SS** (Alternative 5' Splice Site): Different 5' donor boundary, <100bp difference

- **A3SS** (Alternative 3' Splice Site): Different 3' acceptor boundary, <100bp difference

- **IR** (Intron Retention): Exon in one isoform spans multiple (≥2) consecutive exons in the other

- **Partial_IR_5** (Partial Intron Retention, 5' boundary): One boundary similar, 5' boundary
  differs by ≥100bp

- **Partial_IR_3** (Partial Intron Retention, 3' boundary): One boundary similar, 3' boundary
  differs by ≥100bp

- **Alt_TSS** (Alternative Transcription Start Site): Different biological first exon TSS position;
  includes `missing_terminal_exons` coordinate payload and `n_terminal_regions` count

- **Alt_TES** (Alternative Transcription End Site): Different biological last exon TES position;
  includes `missing_terminal_exons` coordinate payload and `n_terminal_regions` count

**Event Detection Flow:**

```
For each isoform pair (dominant vs comparator):

STEP 1: Boundary Determination
├─ Order exons TSS→TES (biological order, strand-aware)
├─ Compute TSS and TES positions for both isoforms
└─ Define overlap region for use in Steps 2 and 3

STEP 2: Within-Boundary Event Detection
│  (Skipped entirely for non-overlapping isoforms)
│
├─ 2a: IR Detection
│  ├─ Check each comparator exon vs all dominant exons → detect_ir_simple()
│  │  └─ If overlaps (≥1bp) with ≥2 dominant exons → IR GAIN event
│  ├─ Check each dominant exon vs all comparator exons → detect_ir_simple()
│  │  └─ If overlaps (≥1bp) with ≥2 comparator exons → IR LOSS event
│  └─ Track all (comparator_idx, dominant_idx) exon pairs involved in IR
│     to prevent re-analysis in the next step
│
├─ 2b: Boundary Event Detection (A5SS / A3SS / Partial_IR)
│  └─ For each (comparator, dominant) exon pair with genomic overlap (≥1bp):
│     └─ Skip if this pair is involved in an IR event (from 2a)
│     └─ Call detect_shared_boundary_event()
│        │
│        ├─ MODE A: Exact Boundary Match (tried first)
│        │  ├─ Check if acceptor coordinates exactly equal
│        │  ├─ Check if donor coordinates exactly equal
│        │  └─ If one shared, one differs:
│        │     ├─ Difference <100bp → A5SS or A3SS
│        │     └─ Difference ≥100bp → Partial_IR_5 or Partial_IR_3
│        │
│        └─ MODE B: Overlap-Based Detection (fallback if Mode A finds nothing)
│           ├─ Verify exons overlap by ≥1bp
│           ├─ Compare boundaries in strand-aware manner
│           └─ Call check_boundary_within_exon() for each boundary:
│              ├─ Case A: Comparator boundary within dominant exon
│              │  ├─ <100bp diff → A5SS/A3SS (LOSS direction)
│              │  └─ ≥100bp diff → Partial_IR (LOSS direction)
│              └─ Case B: Comparator boundary extends beyond dominant
│                 ├─ Check if extends into flanking exon
│                 │  ├─ <100bp → A5SS/A3SS (GAIN direction)
│                 │  └─ ≥100bp → IR (GAIN direction)
│                 └─ Doesn't reach flanking exon
│                    ├─ <100bp → A5SS/A3SS (GAIN direction)
│                    └─ ≥100bp → Partial_IR (GAIN direction)
│
├─ 2c: SE Detection (Skipped Exon — strict flanking required)
│  └─ For each dominant exon i within overlap region:
│     ├─ Skip if already involved in IR or boundary events
│     ├─ Require: dom[i-1] overlaps a comp exon AND dom[i+1] overlaps a comp exon
│     └─ If dom[i] has no overlap with any comp exon → SE LOSS event
│  └─ For each comparator exon j within overlap region (symmetric):
│     ├─ Skip if already involved in IR or boundary events
│     ├─ Require: comp[j-1] overlaps a dom exon AND comp[j+1] overlaps a dom exon
│     └─ If comp[j] has no overlap with any dom exon → SE GAIN event
│
└─ 2d: Missing Internal Exons (diagnostic, both directions)
   ├─ LOSS: dominant exons within comparator span with no comparator overlap
   └─ GAIN: comparator exons within dominant span with no dominant overlap
   (Recorded in missing_internal_exons field; not consumed by reconstruction)

STEP 3: Terminal Event Detection
├─ Alt_TSS: Compare biological first exons (strand-aware)
│  ├─ If TSS positions differ → Alt_TSS event
│  ├─ direction: LOSS if dominant extends further 5', GAIN if comparator extends further 5'
│  ├─ missing_terminal_exons: exonic coordinate ranges of the longer isoform's 5' extension,
│  │  computed by walking the longer isoform's exons from its TSS to the shorter isoform's TSS
│  └─ n_terminal_regions: count of coordinate ranges in missing_terminal_exons
└─ Alt_TES: Compare biological last exons (strand-aware)
   ├─ If TES positions differ → Alt_TES event
   ├─ direction: LOSS if dominant extends further 3', GAIN if comparator extends further 3'
   ├─ missing_terminal_exons: exonic coordinate ranges of the longer isoform's 3' extension,
   │  computed by walking the longer isoform's exons from its TES to the shorter isoform's TES
   └─ n_terminal_regions: count of coordinate ranges in missing_terminal_exons

NOTE on non-overlapping isoforms:
  When the two isoforms share no exon-level overlap, Step 2 is skipped entirely.
  Step 3 still fires and correctly captures the complete structural difference via
  Alt_TSS + Alt_TES events. missing_terminal_exons in each event will contain all
  exonic coordinate ranges of the longer isoform's terminal extension (which may
  span multiple complete exons).
```

**Technical Note: IR Exon Pair Tracking**

In Step 2a, IR events are identified when one exon overlaps with multiple exons. Each IR event involves at least 3 exons total: 1 spanning exon and ≥2 spanned exons. To prevent these exons from being re-analyzed in Step 2b (boundary detection), all pairwise combinations of (comparator_exon, dominant_exon) involved in the IR are tracked.

**Why pairs?** Each tracked pair consists of exactly one comparator exon index and one dominant exon index, both involved in the same IR event. In Step 2b, before performing boundary analysis on any (comparator, dominant) exon pair, the code checks if that pair is in the IR tracking list and skips it if present.

**Order consistency:** The tracking always stores pairs in `(comparator_index, dominant_index)` format, regardless of IR direction:
- GAIN direction (comparator spans dominant): stores `(comp=i, dom=j)` for each overlapping dominant exon j
- LOSS direction (dominant spans comparator): stores `(comp=j, dom=i)` for each overlapping comparator exon j

This consistent ordering ensures Step 2b correctly identifies and skips all exon combinations already classified as IR, preventing redundant or conflicting event assignments.

**Technical Note: SE Strict Flanking**

SE detection (Step 2c) requires that both flanking exons (the exons immediately upstream and downstream of the candidate SE) are present and overlapping in the other isoform. An internal exon without both flanking exons present is recorded in `missing_internal_exons` instead, which serves as a diagnostic field indicating structural differences that do not fit the SE pattern. This distinction is important because an internal exon absent in the other isoform but without both flanking exons present may indicate a more complex structural rearrangement.

**Script**: `detect_and_save_events.R` (sources local `event_detection_functions.R`)

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
     - Alt_TSS/Alt_TES: Trim terminal exons using missing_terminal_exons coordinates
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
- `reconstruction_functions.R` — Core reconstruction logic (`reconstruct_dominant_v2()`, `merge_adjacent_exons()`)
- `reconstruct_dominant_isoforms.R` — Main reconstruction pipeline
- `verify_reconstruction.R` — Validation script that compares reconstructed isoforms to expected structures and reports pass/fail status

## Technical Notes

- **GAIN/LOSS Semantics**: Event directions (GAIN/LOSS) are defined from the comparator's perspective:
  - LOSS = comparator LOST sequence (dominant has more exonic sequence)
  - GAIN = comparator GAINED sequence (dominant has less exonic sequence)
- **Biological Exon Order**: Exons are ordered TSS→TES throughout detection and reconstruction.
  Plus strand: ascending genomic coordinate. Minus strand: descending genomic coordinate.
- **Union Exon Approach**: Uses atomic (non-overlapping) union exon segments
- **Strand-Aware**: All event detection and reconstruction is strand-aware
