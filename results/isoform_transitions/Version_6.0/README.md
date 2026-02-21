# Isoform Transition Detection Pipeline - Version 6.0

Complete system for detecting splicing differences between isoform pairs, validating detection accuracy via reconstruction, and analyzing splicing choice profiles across NMD conditions.

## Overview

This pipeline:
1. **Prepares data** from DGEList objects + GFF annotations (nmd/01, core/02-05, nmd/06)
2. **Classifies isoforms and generates comparison pairs** for NMD analysis (nmd/07)
3. **Detects splicing events** between isoform pairs using hierarchical event detection (core/08)
4. **Validates detection** by reconstructing the dominant isoform from comparator + events (core/08 `--reconstruction_check` or core/09)
5. **Analyzes splicing patterns** across cell types and NMD conditions (nmd/09-14)

**Validation Status:**
- Curated test suite (synthetic + real failures): 126/126 tests (100%)
- Real data (GENCODE): 4,258/4,274 pairs (99.6%) — 0 FAILs, 16 ERRORs (missing comparator data)

## Directory Structure

```
Version_6.0/
├── README.md                          # This file
├── ANALYSIS_PLAN.md                   # Full analysis framework (Sections 1-6)
├── METHODS.md                         # Algorithmic documentation
├── DATA_FLOW.md                       # Data pipeline documentation
├── MIGRATION_PLAN.md                  # Dev → scripts migration plan
│
├── scripts/
│   ├── core/                          # Generic, reusable pipeline scripts
│   │   ├── 02_extract_isoform_structures.R   # Parse GENCODE/SQANTI GFF to exon structures
│   │   ├── 03_build_union_exons.R            # Build atomic union exon models per gene
│   │   ├── 04_extract_cds_annotations.R      # Extract CDS start/stop coordinates
│   │   ├── 05_annotate_region_types.R        # Classify union exons as 5'UTR/CDS/3'UTR
│   │   ├── 08_extract_splicing_profiles.R    # Event detection + splicing choice profiles
│   │   │   # --reconstruction_check: on-the-fly reconstruction verification
│   │   │   # --pairs-file <path>: explicit contrast pairs (TSV)
│   │   │   # --output <path>: custom output path
│   │   │   # --test N: limit to first N genes
│   │   ├── 09_validate_reconstruction.R      # Standalone reconstruction validation
│   │   ├── event_detection_functions.R       # Shared: event detection + thresholds
│   │   ├── reconstruction_functions.R        # Shared: reconstruct_dominant_v2 + verify_transcript
│   │   └── visualization_functions.R         # Shared: isoform pair plotting + visual validation
│   │
│   ├── nmd/                           # NMD study-specific scripts
│   │   ├── 01_prepare_expression_data.R      # Load DGEList, filter major isoforms, identify dominants
│   │   ├── 06_filter_to_analysis_subset.R    # Apply expression filters (filterByExpr)
│   │   ├── 07_classify_and_pair.R            # NMD classification + comparison pair generation (C1-C4)
│   │   ├── 09_analyze_complexity_relationship.R
│   │   ├── 10_analyze_cooccurrence.R
│   │   ├── 11_analyze_spatial_patterns.R
│   │   ├── 12_analyze_functional_context.R
│   │   ├── 13_analyze_patterns.R
│   │   └── 14_generate_report.Rmd
│   │
│   ├── tests/                         # Validation test suite (126/126 PASS)
│   │   ├── run_tests.R                       # Test runner (synthetic + real-data cases)
│   │   ├── extract_failure_cases.R
│   │   └── visualize_results.R
│   │
│   ├── dev/                           # Development and diagnostic scripts
│   │   ├── visualize_comparisons.R           # CLI wrapper for visual validation
│   │   ├── summarize_comparisons.R           # Profile summary by comparison type
│   │   ├── diagnose_failures.R
│   │   ├── test_08_on_gencode.R
│   │   └── analyze_08_failures.R
│   │
│   └── archive/                       # Archived earlier versions
│
├── data/                              # Pipeline output data (.rds files)
├── results/                           # Analysis results
├── figures/                           # Generated figures
├── logs/                              # Execution logs
│
├── run_pipeline_03_to_06.sh           # Pipeline runner: data prep
├── run_analysis_07_to_12.sh           # Pipeline runner: analysis
└── check_pipeline_status.sh           # Status checker
```

## Event Detection Algorithm

### Hierarchical Detection

Events are detected in a strict order to prevent double-counting:

```
For each isoform pair (dominant vs comparator):

STEP 1: Boundary Determination
├─ Order exons TSS→TES (biological order, strand-aware)
├─ Compute TSS and TES positions for both isoforms
└─ Define overlap region

STEP 2: Within-Boundary Event Detection
│  (Skipped entirely for non-overlapping isoforms)
│
├─ 2a: IR Detection (FIRST — prevents re-analysis in later steps)
│  ├─ Check each comp exon vs all dom exons → IR GAIN if spans ≥2
│  ├─ Check each dom exon vs all comp exons → IR LOSS if spans ≥2
│  ├─ Emit IR_diff events when retained exon boundaries differ
│  └─ Track all exon pairs involved in IR
│
├─ 2b: Boundary Event Detection (A5SS / A3SS / Partial_IR)
│  └─ For overlapping exon pairs NOT involved in IR:
│     ├─ Boundary difference → A5SS or A3SS or Partial_IR
│     ├─ Asymmetric terminal pairs → emit second_event for both boundaries
│     └─ Dual-boundary internal pairs → decompose into two events
│
├─ 2c: SE Detection (strict flanking required)
│  └─ Internal exons with no overlap AND both flanks present
│
└─ 2d: Missing Internal Exons
   ├─ Dom exons within comp span with no comp overlap → Missing_Internal LOSS
   └─ Comp exons within dom span with no dom overlap → Missing_Internal GAIN

STEP 3: Terminal Event Detection
├─ Alt_TSS: Compare biological first exons
│  ├─ Triggered when TSS differs > 20bp OR first exons don't overlap
│  └─ missing_terminal_exons: coordinate ranges of longer isoform's 5' extension
└─ Alt_TES: Compare biological last exons
   ├─ Triggered when TES differs > 20bp OR last exons don't overlap
   └─ missing_terminal_exons: coordinate ranges of longer isoform's 3' extension
```

### Event Types

| Event Type | Description |
|------------|-------------|
| **SE** | Skipped Exon — exon absent in other isoform, both flanking exons present |
| **Missing_Internal** | Internal exon absent, flanking condition NOT met |
| **A5SS** | Alternative 5' splice site — shared acceptor, different donor (< 100bp) |
| **A3SS** | Alternative 3' splice site — shared donor, different acceptor (< 100bp) |
| **Partial_IR_5** | Partial intron retention at 5' boundary (≥ 100bp) |
| **Partial_IR_3** | Partial intron retention at 3' boundary (≥ 100bp) |
| **IR** | Intron retention — one exon spans ≥2 exons in the other isoform |
| **IR_diff_5/3/5_3** | IR with boundary mismatches at 5', 3', or both ends |
| **Alt_TSS** | Alternative transcription start site (> 20bp difference OR non-overlapping first exons) |
| **Alt_TES** | Alternative transcription end site (> 20bp difference OR non-overlapping last exons) |

### GAIN/LOSS Semantics

From the **comparator's** perspective:
- **LOSS** = comparator lost sequence (dominant has more) → reconstruction **ADDs** to comparator
- **GAIN** = comparator gained sequence (dominant has less) → reconstruction **REMOVEs** from comparator

### Dominance

In the development/validation workflow, dominance is determined by total exonic sequence length (dominant = more exonic bp). In the main pipeline, dominance is determined by mean CPM expression across samples.

## Reconstruction Validation

The reconstruction system validates event detection accuracy by rebuilding the dominant isoform from the comparator plus detected events.

### Two-Phase Reconstruction

**Phase 1 — Internal Events (LOSS before GAIN, within shared boundaries):**
- IR: Direct coordinate reconstruction — LOSS merges split exons, GAIN splits using `ir_split_exons`
- SE/Missing_Internal: Add or remove exons by coordinate range
- A5SS/A3SS/Partial_IR: Algebraic boundary adjustment from five_prime/three_prime
- Pre-merge orphan removal: terminal orphans removed before merge step
- `merge_adjacent_exons()` combines touching/overlapping segments

**Phase 2 — Terminal Events (LOSS before GAIN):**
- LOSS: Add missing exon ranges from `missing_terminal_exons`
- GAIN: Use `five_prime` (dominant's terminal boundary) to truncate/remove exons beyond the dominant's extent
- Remove any remaining orphan terminal exons

**Key Principle:** Reconstruction uses ONLY event type + associated coordinates. No union exon lookups are performed — all event types (including IR) reconstruct directly from event coordinate fields.

## Events File Format

Tab-separated, 15 columns:

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier |
| `dominant_transcript_id` | Dominant isoform |
| `comparator_transcript_id` | Comparator isoform |
| `event_type` | Event classification |
| `direction` | GAIN or LOSS |
| `chr`, `strand` | Genomic location |
| `five_prime`, `three_prime` | Event boundary coordinates |
| `bp_diff` | Base pair difference |
| `missing_terminal_exons` | Coordinate ranges for Alt_TSS/Alt_TES |
| `orphan_terminal_exons` | Comparator exons to remove |
| `ir_split_exons` | Exon coordinate ranges for IR reconstruction |
| `dom_junctions`, `comp_junctions` | Splice junctions (informational) |

## Dependencies

**R packages:** tidyverse, ggplot2, patchwork, edgeR, rtracklayer
**System tools:** bgzip, tabix

## Comparison Framework (NMD-Specific)

Four comparison types, each run across 7 sample sets (all_samples + 6 cell types):

| | Dominant | Comparator | Pairs/gene |
|---|---|---|---|
| **C1** | Dominant non-NMD (DMSO) | Dominant NMD-sensitive (Smg1i) | 1 |
| **C2** | Top non-NMD by CPM (DMSO) | Top NMD-sensitive by CPM (Smg1i) | 1 |
| **C3** | Dominant non-NMD (DMSO) | All other non-NMD | Multiple |
| **C4** | Dominant non-NMD (DMSO) | Next-best non-NMD by CPM | 1 |

Pairs are deduplicated across all 28 sets before running event detection once via core/08.

## Contact

Project: NMD Lung Cells 2026
