# Event Detection and Reconstruction Workflow

## Purpose

This workflow detects splicing events between isoform pairs and validates detection accuracy through reconstruction. Starting from a comparator isoform plus detected events, the system reconstructs the dominant isoform. If the reconstruction exactly matches the original dominant, the event detection was correct.

**Current Status (2026-02-20):**
- Curated test suite (synthetic + real failures): **126/126 tests passing (100%)**
- Real data (GENCODE): **4,258/4,274 = 99.6%** — 0 FAILs, 16 ERRORs (all "No comparator exons")

## Architecture

This is a self-contained development environment. All scripts source dependencies locally or from `../../scripts/`. The workflow has been validated on both synthetic and real (GENCODE) data.

### Key Design Principle

**Reconstruction uses ONLY event type + associated coordinates.** All event types — including IR — reconstruct directly from event coordinate fields (`five_prime`, `three_prime`, `missing_terminal_exons`, `ir_split_exons`). No union exon lookups are performed during reconstruction. This eliminates dependencies on UE coverage and avoids gene_id mismatch issues between PacBio (PB.xxxx) and GENCODE (ENSG) identifiers.

## Workflow Overview

```
Input GTF + Pairs File
         ↓
  1. Create Atomic Union Exons        → union_exons.tsv.gz
         ↓
  2. Detect Splicing Events            → events.tsv
     (determines dominance internally)
         ↓
  3. Extract Comparator GTF            → comparator.gtf
         ↓
  4. Reconstruct Dominant Isoforms     → reconstructed.gtf
         ↓
  5. Verify Reconstruction             → verification.tsv
```

**Key Insight:** Dominance is NOT pre-specified. The `detect_and_save_events.R` script determines which isoform is dominant based on total exonic sequence length (dominant = more exonic bp).

## Input Data Requirements

### 1. GTF File

Standard GTF2 format with exon features. Required attributes: `gene_id`, `transcript_id`. Coordinates are 1-based, inclusive.

### 2. Pairs File

Tab-separated with columns: `gene_id`, `isoform_A`, `isoform_B`, `chr`, `strand`

No dominance specified — just two isoforms to compare per gene. Dominance is determined internally.

## Event Types

### Core Event Types

| Event Type | Description | Classification Threshold |
|------------|-------------|------------------------|
| **SE** | Skipped Exon — exon in one isoform absent in the other, with both flanking exons present | Strict flanking requirement |
| **Missing_Internal** | Internal exon absent in other isoform, flanking condition NOT met | Weaker than SE |
| **A5SS** | Alternative 5' splice site — shared acceptor, different donor | < 100bp difference |
| **A3SS** | Alternative 3' splice site — shared donor, different acceptor | < 100bp difference |
| **Partial_IR_5** | Partial intron retention at 5' boundary | ≥ 100bp, one shared boundary |
| **Partial_IR_3** | Partial intron retention at 3' boundary | ≥ 100bp, one shared boundary |
| **IR** | Intron retention — one exon spans ≥2 exons in the other isoform | ≥2 overlapping exons |
| **IR_diff_5** | IR where retained exon 5' boundary differs from spanned exons | Emitted alongside IR |
| **IR_diff_3** | IR where retained exon 3' boundary differs from spanned exons | Emitted alongside IR |
| **IR_diff_5_3** | IR where both boundaries differ | Emitted alongside IR |
| **Alt_TSS** | Alternative transcription start site | > 20bp TSS difference |
| **Alt_TES** | Alternative transcription end site | > 20bp TES difference |

### GAIN/LOSS Semantics (CRITICAL)

**GAIN and LOSS are from the COMPARATOR'S perspective:**

- **LOSS** = Comparator LOST sequence (dominant has MORE)
  - Reconstruction action: **ADD** regions to comparator
- **GAIN** = Comparator GAINED sequence (dominant has LESS)
  - Reconstruction action: **REMOVE** regions from comparator

## Event Detection Algorithm

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
│     ├─ Exact boundary match → A5SS or A3SS or Partial_IR
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
│  └─ missing_terminal_exons: coordinate ranges of longer isoform's 5' extension
└─ Alt_TES: Compare biological last exons
   └─ missing_terminal_exons: coordinate ranges of longer isoform's 3' extension
```

### Key Detection Features

- **Dual-boundary decomposition**: Internal exon pairs where both boundaries differ are decomposed into two independent events (e.g., A5SS + A3SS, or Partial_IR_5 + Partial_IR_3)
- **Asymmetric terminal handling**: When one exon is terminal but its pair is not, both boundaries are emitted (the terminal-facing boundary and the splice-site-facing boundary)
- **Junction tracking**: All events include `dom_junctions` and `comp_junctions` fields
- **Non-overlapping isoforms**: Step 2 is skipped; Alt_TSS + Alt_TES capture the complete structural difference via `missing_terminal_exons`

## Reconstruction Algorithm

### Two-Phase Approach

**Phase 1: Internal Events (IR, SE, A5SS/A3SS, Partial_IR, Missing_Internal)**
Events are sorted LOSS-before-GAIN, then applied:
1. IR LOSS → Merge comparator's split exons into one retained-intron exon using event coordinates
2. IR GAIN → Remove comparator's retained-intron exon, add split exons from `ir_split_exons` field
3. SE/Missing_Internal LOSS → Add exons by coordinate range
4. SE/Missing_Internal GAIN → Remove exons by coordinate range
5. A5SS/A3SS/Partial_IR LOSS → Extend exon boundaries algebraically
6. A5SS/A3SS/Partial_IR GAIN → Shrink exon boundaries algebraically

**Pre-merge orphan removal:** Before merging, terminal event orphans (comparator exons at 5'/3' ends that don't overlap the dominant) are removed. This prevents `merge_adjacent_exons` from erroneously consolidating orphan exons with internally-modified adjacent exons.

**Post-Phase 1:**
- `merge_adjacent_exons()` combines touching/overlapping exon segments

**Phase 2: Terminal Events (Alt_TSS, Alt_TES)**
Applied after merge to avoid interference with internal event boundaries:
1. Alt_TSS LOSS → Add missing terminal exons at 5' end
2. Alt_TES LOSS → Add missing terminal exons at 3' end
3. Alt_TSS GAIN → Remove terminal exons/trim boundary
4. Alt_TES GAIN → Remove terminal exons/trim boundary
5. Remove any remaining orphan terminal exons

**Final:** Exons are sorted in biological order

### Core Functions (reconstruction_functions.R)

| Function | Purpose |
|----------|---------|
| `reconstruct_dominant_v2()` | Main reconstruction pipeline (two-phase + pre-merge orphan removal) |
| `apply_event_union_based()` | Routes events to appropriate handler (all use direct coordinates) |
| `modify_exon_boundary()` | Algebraic boundary adjustment for A5SS/A3SS/Partial_IR |
| `modify_terminal_exon()` | Handles Alt_TSS/Alt_TES using coordinate ranges |
| `add_union_exons()` | Adds exons for SE/Missing_Internal LOSS events |
| `remove_union_exons()` | Removes exons for SE/Missing_Internal GAIN events |
| `merge_adjacent_exons()` | Combines adjacent/overlapping exon segments |

### Boundary Modification Logic

`modify_exon_boundary()` uses algebraic derivation from event coordinates — no union exon lookups:
- Identifies the target exon (overlapping five_prime/three_prime)
- Computes new boundary from `min()`/`max()` of event coordinates ± 1
- Strand-aware: accounts for donor/acceptor semantics on plus vs minus strand

### Terminal Modification Logic

`modify_terminal_exon()` uses event coordinates directly:
- **LOSS**: Add missing exon ranges from `missing_terminal_exons` field
- **GAIN**: Use `five_prime` (dominant's terminal boundary) to truncate/remove exons.
  Cut direction is determined by `event_type + strand`:
  Alt_TES+/Alt_TSS- → remove high side; Alt_TSS+/Alt_TES- → remove low side.

## Events File Format

Tab-separated with columns:
- `gene_id`, `dominant_transcript_id`, `comparator_transcript_id`
- `event_type`, `direction`
- `chr`, `five_prime`, `three_prime`, `strand`
- `bp_diff` (NA for IR)
- `missing_terminal_exons` (coordinate ranges for Alt_TSS/Alt_TES)
- `orphan_terminal_exons` (comparator exons to remove)
- `ir_split_exons` (for IR events)
- `dom_junctions`, `comp_junctions`

## Reconstructed GTF Format

Transcript IDs use `DOM::COMP` format: `"ENST00000497506.5::ENST00000412894.5"`

## Validation

### Verification Criteria
- **Exact match**: All exon coordinates identical
- **TSS tolerance**: ±20bp at transcription start site
- **TES tolerance**: ±20bp at transcription end site

### Current Results

**Curated Test Suite**: 126/126 (100%) — 44 synthetic + 82 real-world cases including all IR subtypes, both GTF-built and production union exons

**Real Data (GENCODE)**: 4,258/4,274 (99.6%)
- 0 FAILs (exact coordinate mismatches)
- 16 ERRORs: all "No comparator exons" (comparator transcript not in comparator.gtf)

## Scripts

### Pipeline Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| `prepare_test_data.R` | Validate GTF, randomly select dominants | `Rscript prepare_test_data.R <input.gtf> <output_dominant_mapping.tsv>` |
| `build_atomic_union_exons.R` | Create tabix-indexed atomic union exons | `Rscript build_atomic_union_exons.R <input.gtf> <output_file>` |
| `generate_pairs_from_dominant.R` | Generate pairs from dominant mapping | `Rscript generate_pairs_from_dominant.R <gtf> <dominant.tsv> <pairs.tsv>` |
| `detect_and_save_events.R` | Detect events, determine dominance | `Rscript detect_and_save_events.R <gtf> <pairs.tsv> <events.tsv>` |
| `reconstruct_dominant_isoforms.R` | Reconstruct dominant from comparator + events | `Rscript reconstruct_dominant_isoforms.R <comp.gtf> <events.tsv> <ue.tsv.gz> <out.gtf> <log.tsv>` |
| `verify_reconstruction.R` | Compare reconstructed vs original | `Rscript verify_reconstruction.R <original.gtf> <reconstructed.gtf> <events.tsv> <verification.tsv>` |

### Function Libraries

| Script | Purpose |
|--------|---------|
| `event_detection_functions.R` | All event detection logic (sourced by `detect_and_save_events.R`) |
| `reconstruction_functions.R` | All reconstruction logic (sourced by `reconstruct_dominant_isoforms.R`) |
| `visualization_functions.R` | Gene structure plotting |

### Visualization Scripts

| Script | Purpose |
|--------|---------|
| `visualize_verification.R` | Visualize verification failures |
| `visualize_all_three.R` | Plot dominant, reconstructed, and comparator side-by-side |
| `visualize_failures.R` | Debug failing test cases |

## Complete Validation Workflow

```bash
cd development/reconstruction/

# Step 1: Validate GTF and select dominants
Rscript prepare_test_data.R <input.gtf> <data_dir>/dominant_mapping.tsv

# Step 2: Build atomic union exons
Rscript build_atomic_union_exons.R <input.gtf> <data_dir>/union_exons.tsv

# Step 3: Generate pairs
Rscript generate_pairs_from_dominant.R <input.gtf> <data_dir>/dominant_mapping.tsv <data_dir>/pairs.tsv

# Step 4: Detect events
Rscript detect_and_save_events.R <input.gtf> <data_dir>/pairs.tsv <data_dir>/events.tsv

# Step 5: Extract comparator GTF
tail -n +2 <data_dir>/events.tsv | cut -f3 | sort -u > /tmp/comparator_ids.txt
grep -Ff /tmp/comparator_ids.txt <input.gtf> > <data_dir>/comparator.gtf

# Step 6: Reconstruct
Rscript reconstruct_dominant_isoforms.R \
  <data_dir>/comparator.gtf <data_dir>/events.tsv \
  <data_dir>/union_exons.tsv.gz <data_dir>/reconstructed.gtf \
  <data_dir>/reconstruction_log.tsv

# Step 7: Verify
Rscript verify_reconstruction.R \
  <input.gtf> <data_dir>/reconstructed.gtf \
  <data_dir>/events.tsv <data_dir>/verification.tsv

# Step 8 (Optional): Visualize
Rscript visualize_all_three.R \
  <input.gtf> <data_dir>/reconstructed.gtf \
  <data_dir>/comparator.gtf <data_dir>/events.tsv \
  <data_dir>/three_isoform_viz.pdf
```

## Working Directories

- `synthetic_data/` — Synthetic test data (44 test cases, 100% pass)
- `real_data/` — Real GENCODE data testing (4,274 pairs, 97.2% pass)

## Detection Thresholds

| Constant | Value | Purpose |
|----------|-------|---------|
| `TSS_TOLERANCE` | 20 bp | Minimum TSS difference for Alt_TSS event |
| `TES_TOLERANCE` | 20 bp | Minimum TES difference for Alt_TES event |
| `SPLICE_SITE_THRESHOLD` | 100 bp | < 100bp → A5SS/A3SS; ≥ 100bp → Partial_IR |

## Dependencies

**R packages:** tidyverse (dplyr, readr, tidyr, stringr, purrr), ggplot2, patchwork

**System tools:** bgzip, tabix
