# Isoform Transition Detection Pipeline - Version 6.0

Complete system for detecting splicing differences between isoform pairs, validating detection accuracy via reconstruction, and analyzing splicing choice profiles across NMD conditions.

## Overview

This pipeline:
1. **Prepares data** from DGEList objects + GFF annotations (nmd/01, core/02-05, nmd/06)
2. **Classifies isoforms and generates comparison pairs** for NMD analysis (nmd/07)
3. **Detects splicing events** between isoform pairs using hierarchical event detection (core/08)
4. **Validates detection** by reconstructing the dominant isoform from comparator + events (core/08 `--reconstruction_check` or core/09)
5. **Analyzes splicing patterns** — complexity, co-occurrence, spatial, functional, patterns (nmd/09-13)
6. **Runs per-comparison downstream** — filters profiles to each comparison×run and runs 09-13 (nmd/14)
7. **Compares NMD vs baseline** — paired/unpaired tests, meta-analysis across cell types (nmd/15)
8. **PTC analysis** — premature termination codon prevalence, NMD association, frame disruption, rMATS concordance (ptc/01-05, 11)

All pipeline scripts accept `--source oarfish` (default) or `--source isocall`. The **isocall** analysis is the primary/active analysis. Oarfish results have been archived.

**Validation Status (isocall):**
- Curated test suite (synthetic + real failures): 126/126 tests (100%)
- GENCODE validation: 4,274/4,274 pairs (100%)
- Reconstruction verification: 43,513/43,513 pairs (100%)

## Directory Structure

```
Version_6.0/
├── README.md                          # This file
├── ISOCALL_RESULTS_SUMMARY.md         # Full results summary (isocall analysis)
├── METHODS.md                         # Algorithmic documentation
│
├── scripts/
│   ├── core/                          # Generic, reusable pipeline scripts
│   │   ├── 02_extract_isoform_structures.R   # Parse GFF to exon structures
│   │   ├── 03_build_union_exons.R            # Build atomic union exon models per gene
│   │   ├── 04_extract_cds_annotations.R      # Extract CDS start/stop coordinates
│   │   ├── 05_annotate_region_types.R        # Classify union exons as 5'UTR/CDS/3'UTR
│   │   ├── 08_extract_splicing_profiles.R    # Event detection + splicing choice profiles
│   │   ├── 09_validate_reconstruction.R      # Standalone reconstruction validation
│   │   ├── event_detection_functions.R       # Shared: event detection + thresholds
│   │   ├── reconstruction_functions.R        # Shared: reconstruct_dominant_v2 + verify_transcript
│   │   └── visualization_functions.R         # Shared: isoform pair plotting
│   │
│   ├── nmd/                           # NMD study-specific scripts
│   │   ├── 01_prepare_expression_data.R      # Load DGEList, filter major isoforms, identify dominants
│   │   ├── 06_filter_to_analysis_subset.R    # Apply expression filters (filterByExpr)
│   │   ├── 07_classify_and_pair.R            # NMD classification + comparison pair generation (C1-C4)
│   │   ├── 09_analyze_complexity_relationship.R  # Complexity vs event count
│   │   ├── 10_analyze_cooccurrence.R             # Event co-occurrence
│   │   ├── 11_analyze_spatial_patterns.R         # Positional bias, topology, proximity
│   │   ├── 12_analyze_functional_context.R       # Regional distribution, ORF impact
│   │   ├── 13_analyze_patterns.R                 # Profile type classification
│   │   ├── 14_run_per_comparison.R               # Per-comparison downstream runner
│   │   └── 15_compare_across_comparisons.R       # Cross-comparison statistical framework
│   │
│   ├── tests/                         # Validation test suite (126/126 PASS)
│   │   ├── run_tests.R
│   │   ├── extract_failure_cases.R
│   │   └── visualize_results.R
│   │
│   ├── dev/                           # Development and diagnostic scripts
│   │   ├── visualize_comparisons.R
│   │   ├── summarize_comparisons.R
│   │   ├── diagnose_failures.R
│   │   ├── test_08_on_gencode.R
│   │   └── analyze_08_failures.R
│   │
│   └── archive/                       # Archived earlier versions
│
├── results/
│   └── ptc/                           # PTC analysis
│       ├── scripts/                   # PTC analysis scripts (01-11)
│       │   ├── 01_compute_ptc_status.R           # Compute PTC status for coding isoforms
│       │   ├── 02_ptc_nmd_association.R          # PTC → NMD association tests
│       │   ├── 03_ptc_logfc_distributions.R      # logFC distributions by PTC status
│       │   ├── 04_ptc_signal_exploration.R       # Systematic 10-angle signal search
│       │   ├── 05_ptc_in_comparisons.R           # PTC prevalence + frame disruption in pairs
│       │   └── 11_rmats_longread_concordance.R   # rMATS–long-read concordance
│       ├── results/isocall/           # Isocall PTC results (TSVs, RDS)
│       ├── figures/isocall/           # Isocall PTC figures (PDFs)
│       ├── METHODS.md                 # PTC methodology
│       └── RESULTS.md                 # PTC results narrative
│
├── data/
│   └── isocall/                       # Active isocall pipeline data (.rds files)
│
├── comparisons/
│   └── isocall_nonNMD_0.50/           # Active isocall comparison results
│       ├── C1/, C2/, C4/             # Per-comparison pair files + per-run results
│       ├── deduplicated/             # Pooled dedup profiles + downstream results
│       └── cross_comparison/         # NMD vs baseline statistical framework
│
└── archive/                           # Deprecated oarfish results (pending deletion)
    ├── data/                          # Oarfish pipeline data
    ├── nonNMD_0.50/, nonNMD_0.95/    # Oarfish comparison results
    ├── deduplicated/                  # Oarfish pooled dedup
    ├── results/                       # Oarfish PTC + pipeline results
    ├── figures/, logs/, testing/      # Oarfish auxiliary
    ├── development/, misc_reports/    # Dev artifacts
    └── RESULTS_SUMMARY.md, etc.      # Oarfish-era documentation
```

## Source Parameterization

All scripts accept `--source oarfish` (default) or `--source isocall`:

| Component | Oarfish | Isocall |
|-----------|---------|---------|
| Data dir | `data/` | `data/isocall/` |
| DE dir | `longread_dge/` | `isocall_dge/` |
| DE ID column | `txid` | `transcript_id` |
| Novel ID format | `PB.NNNN.NNN` | `ENSG*.novelN` |
| Gene ID format | Unversioned ENSG | Versioned ENSG |
| Comparison prefix | `nonNMD_` | `isocall_nonNMD_` |
| Output dir (PTC) | `results/ptc/results/` | `results/ptc/results/isocall/` |

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

## Comparison Framework (NMD-Specific)

Four comparison types, each run across 7 sample sets (all_samples + 6 cell types):

| | Dominant | Comparator | Pairs/gene |
|---|---|---|---|
| **C1** | Dominant non-NMD (DMSO) | Dominant NMD-sensitive (Smg1i) | 1 |
| **C2** | Top non-NMD by CPM (DMSO) | Top NMD-sensitive by CPM (Smg1i) | 1 |
| **C3** | Dominant non-NMD (DMSO) | All other non-NMD | Multiple |
| **C4** | Dominant non-NMD (DMSO) | Next-best non-NMD by CPM | 1 |

Pairs are deduplicated across all 28 sets before running event detection once via core/08.

**Downstream scope (Scripts 14-15):** C3 is excluded from cross-comparison analysis because C4 is a strict subset of C3, creating pseudo-replication. DO cell type is excluded from downstream analysis due to insufficient NMD pairs (C1: 1, C2: 6 at 0.50 threshold), though DO still contributes to all_samples isoform classification in Script 07.

## Dependencies

**R packages:** tidyverse, ggplot2, patchwork, edgeR, rtracklayer, metafor
**System tools:** bgzip, tabix

## Contact

Project: NMD Lung Cells 2026
