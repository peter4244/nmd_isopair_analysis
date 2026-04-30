# Isoform Choice Analysis: Computational Methods

**Version:** 6.0
**Date:** 2026-02-21
**Purpose:** Detailed algorithmic documentation for splicing event detection and isoform comparison

---

## Table of Contents

1. [Overview](#overview)
2. [Union Exon Construction](#union-exon-construction)
3. [Event Detection Algorithms](#event-detection-algorithms)
4. [Terminal Boundary Rules](#terminal-boundary-rules)
5. [Validation Framework](#validation-framework)
6. [Statistical Methods](#statistical-methods)
7. [Reading Frame and PTC Analysis](#reading-frame-and-ptc-analysis)
8. [Cross-Comparison Statistical Framework](#cross-comparison-statistical-framework)
9. [Implementation Details](#implementation-details)
10. [Publication Report Methods (mashr Classification)](#publication-report-methods-mashr-classification)

---

## Overview

### Biological Framework

When cells produce alternative isoforms instead of the dominant isoform, they make deliberate splicing choices. We characterize these choices by comparing each non-dominant isoform to the gene's dominant isoform, identifying:

1. **Transcript boundary changes** (Alt_TSS, Alt_TES)
2. **Internal exon inclusion/exclusion** (SE, Missing_Internal)
3. **Exon boundary modifications** (A5SS, A3SS, Partial_IR_5, Partial_IR_3)
4. **Intron retention** (IR, IR_diff_5, IR_diff_3, IR_diff_5_3)

Each event carries a **direction** (GAIN or LOSS) from the comparator's perspective:
- **LOSS** = comparator lost sequence relative to dominant (dominant has more). Reconstruction: ADD regions to comparator.
- **GAIN** = comparator gained sequence relative to dominant (dominant has less). Reconstruction: REMOVE regions from comparator.

### Computational Approach

Our analysis proceeds in three phases:

**Phase 1: Data Preparation**
- Extract isoform structures from GFF3 files (GENCODE + SQANTI)
- Construct union exon models per gene
- Annotate region types (5'UTR, CDS, 3'UTR)
- Identify dominant isoforms

**Phase 2: Classification, Pair Generation, and Event Detection**
- Classify isoforms as NMD-sensitive or non-NMD based on DE results (nmd/07)
- Generate comparison pairs across 4 comparison types x 7 sample sets (nmd/07)
- Compare each comparator to dominant isoform (core/08)
- Detect all splicing events using a hierarchical pipeline: IR -> boundary shifts -> exon skipping -> terminal events
- Reconstruct dominant isoform from comparator + events and verify coordinates match
- Generate splicing choice profiles
- Shared libraries: `scripts/core/event_detection_functions.R`, `scripts/core/reconstruction_functions.R`

**Phase 3: Downstream Statistical Analysis**
- Complexity relationship (nmd/09)
- Co-occurrence analysis (nmd/10)
- Spatial organization (nmd/11)
- Functional context (nmd/12)
- Pattern classification (nmd/13)

**Phase 4: Cross-Comparison Analysis**
- Per-comparison downstream runner (nmd/14): filters profiles to each comparison×run, runs Scripts 09-13
- Cross-comparison statistical framework (nmd/15): NMD vs baseline comparison, meta-analysis, sensitivity analyses

---

## Union Exon Construction

### Algorithm Overview

For each gene, we construct a **union exon model** that represents all possible exonic regions used by any isoform of that gene. This creates a common coordinate system for comparing isoforms.

### Definition

An **atomic union exon** is a minimal genomic segment such that:
1. Every isoform exon boundary is also a union exon boundary
2. Each atomic UE is fully contained within at least one isoform exon
3. No atomic UE can be further split without violating these properties

### Construction Algorithm (Atomic Boundary-Type Splitting)

```
Input: Set of isoforms I₁, I₂, ..., Iₙ for gene G
       Each isoform Iᵢ has exons E₁ᵢ, E₂ᵢ, ..., Eₘᵢ
       Each exon Eⱼᵢ has coordinates [start, end] (inclusive)

Output: Set of atomic union exons U₁, U₂, ..., Uₖ
        Mapping M: (isoform, union_exon) → present/absent

Step 1: Collect all unique boundaries
  B = sorted unique set of all exon starts and ends across all isoforms

Step 2: Classify each boundary
  For each b in B:
    is_start = b appears as an exon_start in any isoform
    is_end = b appears as an exon_end in any isoform
    type(b) = "BOTH" if is_start AND is_end
              "START" if is_start only
              "END" if is_end only

Step 3: Split into atomic segments
  current_start = B[1]
  segments = []

  For each internal boundary b (B[2] through B[n-1]):
    If type(b) == "START":
      # Split before: exon starting here needs its own boundary
      segments += [current_start, b-1]
      current_start = b

    Else if type(b) == "END":
      # Split after: exon ending here gets included
      segments += [current_start, b]
      current_start = b + 1

    Else if type(b) == "BOTH":
      # Position is both an exon-end and exon-start (different isoforms)
      # Must produce THREE segments so both can claim position b
      segments += [current_start, b-1]
      segments += [b, b]           # Single-base segment
      current_start = b + 1

  # Final segment
  segments += [current_start, B[n]]

  # Filter: keep only segments where start <= end
  # Filter: keep only segments covered by at least one exon

Step 4: Create isoform-union exon mapping
  For each isoform Iᵢ:
    For each atomic UE Uₖ:
      # Strict containment: UE must be fully within isoform exon
      M[Iᵢ, Uₖ] = any exon Eⱼᵢ where
        Eⱼᵢ.start <= Uₖ.start AND Eⱼᵢ.end >= Uₖ.end
```

### Properties

1. **Atomicity**: Every isoform exon boundary is also a union exon boundary
2. **Completeness**: Every exonic position belongs to exactly one atomic UE
3. **Non-overlapping**: Atomic UEs are disjoint
4. **Strict containment**: Every isoform exon decomposes into a set of complete atomic UEs

### Example (BOTH Boundary Case)

```
Isoform A: Exons [100-300], [500-700]
Isoform B: Exons [100-300], [300-600]

Boundaries: 100(START), 300(BOTH), 500(START), 600(END), 700(END)

Atomic segments:
  [100-299]: covered by A and B → UE1
  [300-300]: covered by A and B → UE2  (single-base, BOTH boundary)
  [301-499]: covered by B only → UE3
  [500-600]: covered by A and B → UE4
  [601-700]: covered by A only → UE5

Mapping:
  A: [UE1, UE2, UE4, UE5]         (exon [100-300] → UE1+UE2, exon [500-700] → UE4+UE5)
  B: [UE1, UE2, UE3, UE4]  (exon [100-300] → UE1+UE2, exon [300-600] → UE2+UE3+UE4)
```

---

## Event Detection Algorithms

### Overview

All events are detected by comparing a comparator (non-dominant) isoform to the dominant isoform. We use the union exon framework to identify differences in exon usage and boundaries.

**Detection proceeds in biological priority order:**

| Step | What it detects | Event types |
|------|----------------|-------------|
| 2a | Intron retention | IR, IR_diff_5, IR_diff_3, IR_diff_5_3 |
| 2b | Splice-site boundary shifts | A5SS, A3SS, Partial_IR_5, Partial_IR_3 |
| 2c | Exon skipping / missing exons | SE, Missing_Internal |
| gap | Non-overlapping dom exons in comp span | Missing_Internal |
| 3 | Terminal boundary differences | Alt_TSS, Alt_TES |

**Event Direction (GAIN/LOSS)**:

All events carry a direction from the comparator's perspective:
- **LOSS**: Comparator lost sequence (dominant has more). Reconstruction: ADD regions to comparator.
- **GAIN**: Comparator gained sequence (dominant has less). Reconstruction: REMOVE regions from comparator.

**Junction Tracking**:

Events include `dom_junctions` and `comp_junctions` fields recording the splice junctions involved in each event.

### Coordinate System Conventions

**Splice Site Terminology** (biological standard):
- **5' splice site (donor)**: Where the exon ends and intron begins
  - Plus strand: exon_end
  - Minus strand: exon_start
- **3' splice site (acceptor)**: Where intron ends and exon begins
  - Plus strand: exon_start
  - Minus strand: exon_end

**Strand-Specific Transcription Boundaries**:
- **TSS (Transcription Start Site)**: 5' end of transcript
  - Plus strand: exon_start of first exon
  - Minus strand: exon_end of first exon
- **TES (Transcription End Site)**: 3' end of transcript
  - Plus strand: exon_end of last exon
  - Minus strand: exon_start of last exon

---

### 1. Alternative TSS Detection

**Definition**: Transcription start sites differ between dominant and non-dominant isoforms by more than a tolerance threshold.

**Algorithm**:

```
Input: first_exon_dom (first exon of dominant isoform)
       first_exon_non_dom (first exon of non-dominant isoform)
       strand (+ or -)
       tolerance (default: 20 bp)

Output: tss_changed (TRUE/FALSE)

If strand == "+":
  tss_dom = first_exon_dom.exon_start
  tss_non_dom = first_exon_non_dom.exon_start
Else if strand == "-":
  tss_dom = first_exon_dom.exon_end
  tss_non_dom = first_exon_non_dom.exon_end

tss_diff = |tss_dom - tss_non_dom|

# Primary check: coordinate distance
If tss_diff > tolerance:
  Return TRUE

# Secondary check: even if within tolerance, detect structural change
# when first exons don't overlap — indicates an extra terminal exon,
# not just a minor TSS coordinate shift
first_overlap = (first_exon_dom.exon_start <= first_exon_non_dom.exon_end) AND
                (first_exon_dom.exon_end >= first_exon_non_dom.exon_start)

If NOT first_overlap:
  Return TRUE

Return FALSE
```

**Rationale**: Small coordinate differences (≤20bp) may represent annotation noise or minor isoform variation; we focus on biologically meaningful changes >20bp. However, when the first exons don't overlap at all, this indicates a true structural change (e.g., an extra terminal exon in one isoform) regardless of the TSS coordinate distance. This overlap check catches cases where a short extra exon (e.g., 9-19bp) creates a TSS difference within tolerance but represents a genuine alternative first exon.

**Biological Interpretation**:
- Alternative promoter usage
- Different transcription initiation
- 5'UTR length changes

---

### 2. Alternative TES Detection

**Definition**: Transcription end sites differ between dominant and non-dominant isoforms by more than a tolerance threshold.

**Algorithm**:

```
Input: last_exon_dom (last exon of dominant isoform)
       last_exon_non_dom (last exon of non-dominant isoform)
       strand (+ or -)
       tolerance (default: 20 bp)

Output: tes_changed (TRUE/FALSE)

If strand == "+":
  tes_dom = last_exon_dom.exon_end
  tes_non_dom = last_exon_non_dom.exon_end
Else if strand == "-":
  tes_dom = last_exon_dom.exon_start
  tes_non_dom = last_exon_non_dom.exon_start

tes_diff = |tes_dom - tes_non_dom|

# Primary check: coordinate distance
If tes_diff > tolerance:
  Return TRUE

# Secondary check: non-overlapping last exons indicate structural change
# regardless of coordinate distance
last_overlap = (last_exon_dom.exon_start <= last_exon_non_dom.exon_end) AND
               (last_exon_dom.exon_end >= last_exon_non_dom.exon_start)

If NOT last_overlap:
  Return TRUE

Return FALSE
```

**Biological Interpretation**:
- Alternative polyadenylation
- Different transcription termination
- 3'UTR length changes

---

### 3. Alternative 5' Splice Site (A5SS) Detection

**Definition**: Two exons share their acceptor (3' splice site) but have different donors (5' splice sites), with difference <100bp.

**Algorithm**:

```
Input: exon_dom (exon from dominant isoform)
       exon_non_dom (exon from non-dominant isoform)
       strand (+ or -)

Output: detected (TRUE/FALSE), bp_diff (integer)

# Determine acceptor and donor based on strand
If strand == "+":
  acceptor_dom = exon_dom.exon_start
  acceptor_non_dom = exon_non_dom.exon_start
  donor_dom = exon_dom.exon_end
  donor_non_dom = exon_non_dom.exon_end
Else if strand == "-":
  acceptor_dom = exon_dom.exon_end
  acceptor_non_dom = exon_non_dom.exon_end
  donor_dom = exon_dom.exon_start
  donor_non_dom = exon_non_dom.exon_start

# Check if acceptors are shared and donors differ
shares_acceptor = (acceptor_dom == acceptor_non_dom)
differs_donor = (donor_dom != donor_non_dom)

If differs_donor:
  bp_diff = |donor_dom - donor_non_dom|
Else:
  bp_diff = 0

# A5SS: shared acceptor, different donor, <100bp difference
detected = shares_acceptor AND differs_donor AND bp_diff < 100

Return (detected, bp_diff)
```

**Biological Interpretation**:
- Partial exon inclusion/exclusion
- Affects protein sequence (may maintain or disrupt reading frame)
- Changes exon length on the 3' end (donor side)

**100bp Threshold**: Distinguishes A5SS from Partial_IR. Differences ≥100bp are classified as Partial_IR.

---

### 4. Alternative 3' Splice Site (A3SS) Detection

**Definition**: Two exons share their donor (5' splice site) but have different acceptors (3' splice sites), with difference <100bp.

**Algorithm**:

```
Input: exon_dom (exon from dominant isoform)
       exon_non_dom (exon from non-dominant isoform)
       strand (+ or -)

Output: detected (TRUE/FALSE), bp_diff (integer)

# Determine acceptor and donor based on strand
If strand == "+":
  acceptor_dom = exon_dom.exon_start
  acceptor_non_dom = exon_non_dom.exon_start
  donor_dom = exon_dom.exon_end
  donor_non_dom = exon_non_dom.exon_end
Else if strand == "-":
  acceptor_dom = exon_dom.exon_end
  acceptor_non_dom = exon_non_dom.exon_end
  donor_dom = exon_dom.exon_start
  donor_non_dom = exon_non_dom.exon_start

# Check if donors are shared and acceptors differ
shares_donor = (donor_dom == donor_non_dom)
differs_acceptor = (acceptor_dom != acceptor_non_dom)

If differs_acceptor:
  bp_diff = |acceptor_dom - acceptor_non_dom|
Else:
  bp_diff = 0

# A3SS: shared donor, different acceptor, <100bp difference
detected = shares_donor AND differs_acceptor AND bp_diff < 100

Return (detected, bp_diff)
```

**Biological Interpretation**:
- Partial exon inclusion/exclusion
- Affects protein sequence (may maintain or disrupt reading frame)
- Changes exon length on the 5' end (acceptor side)

---

### 5. Partial Intron Retention (Partial_IR_5 / Partial_IR_3) Detection

**Definition**: Two exons share one boundary but differ by ≥100bp on the other boundary, suggesting partial intron retention or extension. Partial IR is now split into two subtypes:
- **Partial_IR_5**: Retention/extension on the 5' (donor) side of the exon
- **Partial_IR_3**: Retention/extension on the 3' (acceptor) side of the exon

**Algorithm**:

```
Input: exon_dom (exon from dominant isoform)
       exon_non_dom (exon from non-dominant isoform)
       strand (+ or -)

Output: detected (TRUE/FALSE), bp_added (integer), side (5' or 3')

# Check if exons share one boundary
shares_start = (exon_dom.exon_start == exon_non_dom.exon_start)
shares_end = (exon_dom.exon_end == exon_non_dom.exon_end)

# Must share exactly one boundary
If NOT (shares_start XOR shares_end):
  Return (FALSE, 0, NA)

# Calculate length difference
len_dom = exon_dom.exon_end - exon_dom.exon_start + 1
len_non_dom = exon_non_dom.exon_end - exon_non_dom.exon_start + 1
len_diff = |len_dom - len_non_dom|

# Partial_IR if difference ≥100bp
detected = (len_diff >= 100)

If NOT detected:
  Return (FALSE, 0, NA)

# Determine which side is affected (5' or 3' in transcript direction)
# Plus strand: shares_start → 3' side variable, shares_end → 5' side variable
# Minus strand: shares_start → 5' side variable, shares_end → 3' side variable
If strand == "+":
  side = "3'" if shares_start else "5'"
Else:
  side = "5'" if shares_start else "3'"

Return (detected, len_diff, side)
```

**Biological Interpretation**:
- Partial intron retention
- Exon extension into adjacent intron
- May introduce premature termination codons (PTCs)

---

### 6. Intron Retention (IR) Detection

**Definition**: An exon in one isoform spans two or more exons in the other isoform, suggesting a retained intron.

**Algorithm**:

```
Input: exon (exon to check for IR)
       other_exons (list of exons from other isoform)

Output: is_ir (TRUE/FALSE)

# Find exons in other_exons that overlap with exon
overlapping = []
For each other_exon in other_exons:
  If exon_overlaps(exon, other_exon):
    Append other_exon to overlapping

# IR if exon spans 2 or more exons in the other isoform
is_ir = (length(overlapping) >= 2)

Return is_ir

# Helper function: exon_overlaps
Function exon_overlaps(exon1, exon2):
  # Check if exon1 and exon2 have any overlapping coordinates
  overlap_start = max(exon1.exon_start, exon2.exon_start)
  overlap_end = min(exon1.exon_end, exon2.exon_end)

  If overlap_end >= overlap_start:
    Return TRUE  # They overlap
  Else:
    Return FALSE  # They don't overlap
```

**Detection Strategy**: Check both directions:
1. For each exon in non-dominant, check if it spans multiple dominant exons
2. For each exon in dominant, check if it spans multiple non-dominant exons

**Biological Interpretation**:
- Complete intron retention
- Often introduces PTCs → triggers NMD
- May represent splicing errors or regulated NMD-targeting

**IR Subtypes (IR_diff)**: When the retained-intron exon's boundaries don't exactly match the outermost split exons, the event is classified as:
- **IR_diff_5**: 5' boundary mismatch (retained exon starts at different position than outermost split exon)
- **IR_diff_3**: 3' boundary mismatch
- **IR_diff_5_3**: Both boundaries mismatch

These capture cases where intron retention co-occurs with a boundary shift.

---

### 7. Skipped Exon (SE) Detection

**Definition**: An exon present in one isoform is absent in the other isoform, AND the flanking exons are "comparable" (share exact coordinates or overlap).

**Key Insight**: SE is distinct from other missing exons because the flanking context is preserved, indicating a deliberate cassette exon skipping event rather than alternative promoter usage or structural rearrangement.

**Critical Requirement**: SE requires splicing events on both sides of the skipped exon. Therefore:
- **Both isoforms must have ≥2 exons** to provide flanking context
- **Monoexonic isoforms cannot have SE** (no internal exons with flanking splices)
- If either isoform is monoexonic, missing exons are due to Alt_TSS/Alt_TES or structural differences, not cassette exon skipping

**Algorithm** (Using Isoform Exon Lists):

```
Input: struct_a (list of exons in isoform A, ordered 5' → 3')
       struct_b (list of exons in isoform B, ordered 5' → 3')

Output: n_se (count of skipped exon events)

# CRITICAL: SE requires splicing on both sides
# If either isoform is monoexonic, SE = 0
If length(struct_a) < 2 OR length(struct_b) < 2:
  Return 0  # Cannot have SE without flanking exons

# Helper function: check if two exons are "comparable"
Function exons_comparable(exon1, exon2):
  # Exact match
  If exon1.exon_start == exon2.exon_start AND exon1.exon_end == exon2.exon_end:
    Return TRUE

  # OR overlap
  If exon_overlaps(exon1, exon2):
    Return TRUE

  Return FALSE

# Check exons in A that are missing in B
n_se = 0

For i = 1 to length(struct_a):
  exon_a = struct_a[i]

  # Is this exon present in B (comparable)?
  in_b = FALSE
  For each exon_b in struct_b:
    If exons_comparable(exon_a, exon_b):
      in_b = TRUE
      break

  # If not in B, check if it's an internal exon with comparable flanks
  If NOT in_b AND i > 1 AND i < length(struct_a):
    prev_exon = struct_a[i-1]
    next_exon = struct_a[i+1]

    # Check if both flanking exons are comparable in B
    prev_comparable = FALSE
    next_comparable = FALSE

    For each exon_b in struct_b:
      If exons_comparable(prev_exon, exon_b):
        prev_comparable = TRUE
      If exons_comparable(next_exon, exon_b):
        next_comparable = TRUE

    # SE if both flanks are comparable
    If prev_comparable AND next_comparable:
      n_se = n_se + 1

# Check exons in B that are missing in A (symmetric)
For i = 1 to length(struct_b):
  exon_b = struct_b[i]

  in_a = FALSE
  For each exon_a in struct_a:
    If exons_comparable(exon_b, exon_a):
      in_a = TRUE
      break

  If NOT in_a AND i > 1 AND i < length(struct_b):
    prev_exon = struct_b[i-1]
    next_exon = struct_b[i+1]

    prev_comparable = FALSE
    next_comparable = FALSE

    For each exon_a in struct_a:
      If exons_comparable(prev_exon, exon_a):
        prev_comparable = TRUE
      If exons_comparable(next_exon, exon_a):
        next_comparable = TRUE

    If prev_comparable AND next_comparable:
      n_se = n_se + 1

Return n_se
```

**Algorithm** (Using Union Exon Framework):

```
Input: comparison (dataframe of union exons with status: "shared", "dominant_only", "non_dominant_only")

Output: n_se (count of skipped exon events)

n_se = 0

For i = 1 to nrow(comparison):
  exon_status = comparison[i].exon_status

  # Skip if this union exon is shared or in neither isoform
  If exon_status in ["shared", "neither"]:
    Continue

  # Check if flanking union exons are shared (comparable)
  has_prev = (i > 1)
  has_next = (i < nrow(comparison))

  prev_shared = has_prev AND comparison[i-1].exon_status == "shared"
  next_shared = has_next AND comparison[i+1].exon_status == "shared"

  # SE if both flanking union exons are shared
  If prev_shared AND next_shared:
    n_se = n_se + 1

Return n_se
```

**Comparable Flanking Exons**:
- **Exact match**: Flanking exons have identical coordinates in both isoforms
- **OR Overlap**: Flanking exons overlap (e.g., due to A5SS/A3SS), indicating they correspond to the same exonic region despite coordinate differences

**Distinction from Other Missing Exons**:
- **SE**: Missing exon with comparable flanks → cassette exon skipping
- **Missing (not SE)**: Missing exon without comparable flanks → likely due to alternative promoter, structural rearrangement, or different terminal exons

**Biological Interpretation**:
- Cassette exon: discrete functional unit that can be included or excluded
- Often regulated in a tissue-specific or developmental manner
- May affect protein domain composition

---

### 8. Missing_Internal Detection

**Definition**: An exon present in only one isoform (within the other's genomic span) where the strict flanking condition for SE is not met.

When a comparator-only or dominant-only exon is found within the other isoform's span, we check if both flanking exons overlap the other isoform. If yes, the event is classified as SE. If no, it is classified as Missing_Internal.

Missing_Internal is also emitted for non-overlapping isoforms: when two isoforms share no exonic overlap, dominant exons within the comparator's genomic span are emitted as Missing_Internal LOSS events.

**Distinction from SE**:
- **SE**: Missing exon with comparable flanks on both sides (cassette exon skipping)
- **Missing_Internal**: Missing exon where at least one flanking exon does not overlap the other isoform (structural context differs)

---

## Terminal Boundary Rules

### Critical Principle

**Transcript boundaries (TSS/TES) are NOT splice sites.** They cannot participate in splice site events (A5SS/A3SS) because there is no splicing occurring at these positions.

### Hierarchical Detection Logic

When comparing terminal exons between isoforms, we apply a hierarchical approach:

#### Rule 1: Never Check A5SS on Last Exons

**Rationale**: On the last exon, the donor position (exon_end on plus, exon_start on minus) is the TES, not a splice site. There is no downstream exon to splice to.

**Implementation**:
```
If comparing last exons:
  Skip A5SS detection
  # A5SS requires a donor site that splices to the next exon
  # Last exon has no next exon → no donor splice site
```

**Example (Plus Strand)**:
```
Last exon A: [1000-1500]  (TES at 1500)
Last exon B: [1000-1600]  (TES at 1600)

These exons:
- Share acceptor at 1000 ✓
- Differ at position 1500 vs 1600 (50bp)
- But 1500 and 1600 are TES positions, NOT donor sites

Correct classification: Alt_TES (50bp > 20bp tolerance)
Incorrect classification: A5SS (this would be wrong!)
```

#### Rule 2: Never Check A3SS on First Exons

**Rationale**: On the first exon, the acceptor position (exon_start on plus, exon_end on minus) is the TSS, not a splice site. There is no upstream exon to splice from.

**Implementation**:
```
If comparing first exons:
  Skip A3SS detection
  # A3SS requires an acceptor site that receives splice from previous exon
  # First exon has no previous exon → no acceptor splice site
```

#### Rule 3: A5SS on First Exons - Requires Overlap Only

**Rationale**: On first exons, the acceptor position is the TSS, which may differ between isoforms (Alt_TSS). A5SS detection on first exons only requires that the exons overlap by ≥50%, indicating they represent the same exonic region with different donor sites. The TSS (acceptor) positions do NOT need to be identical.

**Implementation**:
```
If comparing first exons:
  If first exons overlap by ≥50%:
    # Overlapping first exons → same exonic region
    # Check A5SS at donor side (exon_end on plus, exon_start on minus)
    # TSS can differ (Alt_TSS) - we only require overlap
    Check A5SS detection
  Else:
    # No overlap → different exonic regions
    Skip A5SS detection
```

**Example (Plus Strand)**:
```
Case 1: Overlap ≥50%, A5SS possible (identical TSS)
  First exon A: [1000-1200]  (TSS at 1000)
  First exon B: [1000-1300]  (TSS at 1000)

  - Overlap: 100% ✓
  - Differ at donor: 1200 vs 1300 (100bp)

  Classification: A5SS detected

Case 2: Overlap ≥50%, A5SS possible (different TSS)
  First exon A: [1000-1200]  (TSS at 1000, length 200bp)
  First exon B: [1050-1300]  (TSS at 1050, length 250bp)

  - Different TSS (1000 vs 1050) → Alt_TSS detected
  - Overlap: 1050-1200 = 150bp
  - Overlap % = 150/200 = 75% ≥ 50% ✓
  - Differ at donor: 1200 vs 1300 (100bp)

  Classification: Both Alt_TSS AND A5SS detected

Case 3: No overlap, cannot have A5SS
  First exon A: [1000-1100]  (TSS at 1000)
  First exon B: [1200-1400]  (TSS at 1200)

  - No overlap (0%) → different exonic regions
  - Cannot check A5SS
```

#### Rule 4: A3SS on Last Exons - Requires Overlap Only

**Rationale**: On last exons, the donor position is the TES, which may differ between isoforms (Alt_TES). A3SS detection on last exons only requires that the exons overlap by ≥50%, indicating they represent the same exonic region with different acceptor sites. The TES (donor) positions do NOT need to be identical.

**Implementation**:
```
If comparing last exons:
  If last exons overlap by ≥50%:
    # Overlapping last exons → same exonic region
    # Check A3SS at acceptor side (exon_start on plus, exon_end on minus)
    # TES can differ (Alt_TES) - we only require overlap
    Check A3SS detection
  Else:
    # No overlap → different exonic regions
    Skip A3SS detection
```

### Overlap Calculation

**Purpose**: Determine if two exons correspond to the same exonic region despite coordinate differences.

**Algorithm**:
```
Function calculate_overlap(exon1, exon2):
  # Calculate overlapping region
  overlap_start = max(exon1.exon_start, exon2.exon_start)
  overlap_end = min(exon1.exon_end, exon2.exon_end)
  overlap_length = max(0, overlap_end - overlap_start + 1)

  # Calculate exon lengths
  len1 = exon1.exon_end - exon1.exon_start + 1
  len2 = exon2.exon_end - exon2.exon_start + 1

  # Overlap percentage relative to shorter exon
  overlap_pct = overlap_length / min(len1, len2)

  Return (overlap_pct >= 0.5)
```

**50% Threshold**: We require ≥50% overlap of the shorter exon to consider exons "comparable". This ensures the exons correspond to the same genomic region.

### Summary Table

| Exon Type | Can Have Alt_TSS | Can Have Alt_TES | Can Have A5SS | Can Have A3SS |
|-----------|------------------|------------------|---------------|---------------|
| First exon (plus) | Yes (if >20bp diff) | N/A | Only if TSS shared + overlap ≥50% | **Never** (TSS is not acceptor) |
| Last exon (plus) | N/A | Yes (if >20bp diff) | **Never** (TES is not donor) | Only if TES shared + overlap ≥50% |
| Internal exon | N/A | N/A | Yes (if shared acceptor) | Yes (if shared donor) |
| First exon (minus) | Yes (if >20bp diff) | N/A | Only if TSS shared + overlap ≥50% | **Never** (TSS is not acceptor) |
| Last exon (minus) | N/A | Yes (if >20bp diff) | **Never** (TES is not donor) | Only if TES shared + overlap ≥50% |

**Key**:
- Plus strand: TSS = exon_start, TES = exon_end, acceptor = exon_start, donor = exon_end
- Minus strand: TSS = exon_end, TES = exon_start, acceptor = exon_end, donor = exon_start

### Asymmetric Terminal Exon Handling

When one exon in a pair is terminal but its pair exon is not (asymmetric case), `detect_shared_boundary_event` returns a `second_event` for the splice-site-facing boundary. This handles cases where one isoform's terminal exon aligns with an internal exon in the other isoform, requiring both a terminal boundary event and a splice-site event to fully describe the difference.

The asymmetric condition fires when: `xor(is_first_exon_dom, is_first_exon_comp)` or `xor(is_last_exon_dom, is_last_exon_comp)`, or when both first AND last flags fire simultaneously.

### Dual-Boundary Decomposition

Internal exon pairs where both boundaries differ are decomposed into two independent events (one per boundary). For example, if the 5' boundary differs as A5SS and the 3' boundary differs as A3SS, two separate events are emitted rather than a single ambiguous event. This decomposition applies to all boundary event types (A5SS, A3SS, Partial_IR_5, Partial_IR_3).

---

## Validation Framework

### Synthetic Test Data

We validate all event detection algorithms using synthetic test genes with known ground truth.

#### Test Data Structure

**Components**:
1. **GTF file**: Contains synthetic isoform structures with carefully designed events
2. **Annotations file** (TSV): Describes expected events for each test case
3. **Validation script**: Parses GTF, runs event detection, compares to expected

**Test Coverage**:
- Alt_TSS / Alt_TES (both strands)
- A5SS / A3SS (both strands, on appropriate exon types)
- Partial_IR / IR
- SE (exact flanking matches and overlapping flanks)
- Multi-event cases (e.g., SE + A3SS)
- Edge cases (TSS within tolerance, monoexonic IR)

#### Validation Results

**Version 6.0 Validation**: 44/44 synthetic tests passing (100%)

**Test Categories**:
1. **Basic events** (1 per type): SE, A5SS, A3SS, IR, Alt_TSS, Alt_TES
2. **Strand coverage**: Each event type tested on both plus and minus strands
3. **Terminal boundary rules**: Tests confirming A5SS/A3SS cannot occur at TSS/TES
4. **Overlap-based SE**: Tests SE detection with exact and overlapping flanks
5. **Multi-event**: Tests detecting multiple co-occurring events
6. **Partial_IR events**: Tests for Partial_IR_5 and Partial_IR_3 on internal, first, and last exons
7. **Dual-boundary decomposition**: Tests for internal exon pairs with both boundaries differing
8. **Missing_Internal events**: Tests for exons missing without comparable flanks
9. **IR_diff subtypes**: Tests for IR_diff_5, IR_diff_3, IR_diff_5_3
10. **Non-overlapping isoforms**: Tests for isoform pairs sharing no exonic overlap
11. **Real-world complex cases**: Including minus strand genes with Alt_TSS, SE, and A3SS combinations

**Reconstruction Validation**:
- Reconstruct dominant isoform from comparator + detected events, verify all exon coordinates match
- Curated test suite (synthetic + real failures): **126/126 = 100%** (all event types, both GTF and production UEs)
- GENCODE test: **4258/4274 = 99.6%** (0 FAILs, 16 ERRORs from "No comparator exons")

**Key Validation Case**:
- **TEST_RealCase19_Minus_AltTSS_SE_A3SS**: Based on PB.19746, a real PacBio gene with 15 vs 12 exons on minus strand
  - Tests proper TSS detection on minus strand (TSS at high coordinates)
  - Tests skipped exon detection (2 missing exons)
  - Tests A3SS classification on minus strand internal exons
  - Expected: Alt_TSS, SE, A3SS
  - Detected: Alt_TSS, SE, A3SS ✓

#### Example Test Cases

**TEST_SE_Basic** (Skipped Exon, Exact Flanks):
```
Isoform T1: [1000-1200], [1500-1700], [2000-2200]
Isoform T2: [1000-1200], [2000-2200]

Skipped: [1500-1700]
Flanks: [1000-1200] (exact match), [2000-2200] (exact match)
Expected: SE
Detected: SE ✓
```

**TEST_MultiEvent** (SE + A3SS, Overlapping Flanks):
```
Isoform T1: [19000-19200], [19400-19600], [19800-20000]
Isoform T2: [19000-19200], [19850-20000]

Skipped: [19400-19600]
Flanks: [19000-19200] (exact), [19800-20000] vs [19850-20000] (overlap, A3SS)
Expected: SE, A3SS
Detected: SE, A3SS ✓
```

**TEST_A5SS_Plus** (A5SS on First Exon, Shared TSS):
```
Isoform T1: [35000-35200], [35400-35600]
Isoform T2: [35000-35250], [35400-35600]

First exons: both start at 35000 (TSS shared)
First exons: 100% overlap
Donors differ: 35200 vs 35250 (50bp)
Expected: A5SS
Detected: A5SS ✓
```

**TEST_A3SS_Basic** (Cannot be A5SS on Last Exon):
```
Isoform T1: [3000-3200], [3500-3700]
Isoform T2: [3000-3200], [3550-3700]

Last exons: share donor (TES) at 3700
Last exons: differ at acceptor (3500 vs 3550, 50bp)
A5SS? NO - exon_end is TES, not a donor site
A3SS? YES - exons share donor at 3700, differ at acceptor
Expected: A3SS
Detected: A3SS ✓
```

---

## Statistical Methods

### Co-occurrence Analysis

**Two-Level Approach**:

#### Level 1: Crude Association (Contingency Tables)

**Method**: 2×2 contingency table with Fisher's exact test

```
                Event B Present | Event B Absent
Event A Present      n_both          n_A_only
Event A Absent       n_B_only        n_neither

Odds Ratio = (n_both × n_neither) / (n_A_only × n_B_only)
```

**Purpose**: Rapid screening for potential associations

#### Level 2: Structure-Controlled Association (Regression)

**Method**: Logistic regression with LASSO/ridge regularization

```
Model: has_event_B ~ has_event_A + union_exon_composition + complexity
       + cluster(gene_id)  # Robust standard errors
```

**Purpose**: Control for confounding by isoform complexity and exon composition

**Interpretation**:
- Crude OR ≈ Adjusted OR → minimal confounding
- Crude OR > Adjusted OR → confounded by structure
- Adjusted OR significant → true biological co-occurrence

### Spatial Analysis

**Positional Bias**: Kolmogorov-Smirnov test vs uniform distribution
**Topology Enrichment**: Chi-square test (observed vs expected from simulations)
**Proximity**: Permutation test with 10,000 permutations

### Regional Distribution

**Null Model**: Events distributed proportional to region size (5'UTR, CDS, 3'UTR)

**Test**: Multinomial goodness-of-fit or per-region binomial tests

**Effect Size**: Enrichment ratio = (observed %) / (expected %)

### Multiple Testing Correction

**Strategy**: Benjamini-Hochberg FDR correction within each category
- Category 1: All co-occurrence tests
- Category 2: All spatial tests
- Category 3: All regional tests

**Threshold**: FDR q < 0.05

---

## Reading Frame and PTC Analysis

### Splice Group Assignment

Before evaluating reading frame effects, splicing events are grouped into **splice groups** — sets of events that arise from a single splicing decision (one donor-acceptor choice). Events are grouped using a union-find algorithm based on junction sharing: two events that share at least one splice junction in either `ref_junctions` or `comp_junctions` are assigned to the same group. Events linked transitively (A shares a junction with B, B shares with C) form a single group.

**Rationale**: A multi-exon skip, for example, produces multiple Missing_Internal events (one per skipped exon), but these all reside within a single intron of the comparator and arise from a single splicing decision. Evaluating each skipped exon independently for frameshift potential is misleading — a 4 bp exon and a 2 bp exon skipped together represent a single 6 bp (in-frame) event, not two independent frameshifts that happen to compensate.

**Dual counting**: The analysis tracks both levels of granularity:
- `n_cds_events`: individual exonic changes (how many distinct regions of coding sequence differ)
- `n_splice_groups`: independent splice-site decisions (how many donor-acceptor choices produced those changes)

Frameshift and compensatory evaluation operates at the splice-group level. Event frequency and structural complexity analyses can use either level as appropriate.

### Frame Walk Algorithm

The frame walk (`analyzeFrameWalk()`) traces the cumulative effect of splicing events on the reading frame along each isoform pair. For each pair, CDS-overlapping events are sorted in 5'→3' order (ascending genomic coordinate on plus strand, descending on minus strand) and assigned to splice groups via junction sharing (see above).

The walk proceeds event-by-event, updating the cumulative frame offset for each event. However, frameshift and compensatory status are evaluated at **splice group boundaries** — only after all events in a group have been processed. The group's **net signed CDS change** (sum of signed changes of all member events) determines whether the group shifts the reading frame:

- **Frameshift**: Group net CDS change mod 3 ≠ 0, and the cumulative offset transitions from 0 to non-zero
- **Compensatory**: Group net CDS change restores the cumulative offset to 0
- **Frame-preserving**: Group net CDS change mod 3 = 0 (no frame effect)

**Resolution**: A pair is "resolved" if the cumulative frame offset returns to 0 by the end of the CDS (all frameshifts are compensated). "Unresolved" pairs carry a persistent reading frame shift through the remainder of the CDS, typically producing a premature termination codon in the new frame.

**Reference frame**: The frame walk uses the **reference (dominant) isoform's** CDS as the coordinate system for determining which events overlap coding sequence. This avoids the circular problem of using each isoform's own CDS annotation.

### PTC Detection

PTCs are identified using the canonical 50-nucleotide rule: a stop codon is classified as premature if it is located >50 nt upstream of the last exon-exon junction (EJC) in mRNA coordinates. Additional features:

- `n_downstream_ejcs`: number of EJCs downstream of the stop codon
- `stop_in_last_exon`: whether the stop codon resides in the final exon (expected for normal termination)

Single-exon transcripts are classified as non-PTC (no EJCs possible). CDS annotations come from GENCODE (ENST isoforms) and SQANTI ORF predictions (novel PacBio isoforms).

### Frameshift → PTC Funnel

The funnel traces coding pairs through sequential stages to quantify the fate of frameshifts:

1. **n_total**: All analyzed pairs (coding and non-coding)
2. **n_coding**: Pairs where the reference isoform is coding (`frame_resolved` not NA)
3. **n_frameshift**: Pairs with ≥1 frameshift splice group (group net CDS change not divisible by 3)
4. **n_unresolved**: Frameshifts not compensated by a downstream splice group
5. **PTC breakdown of unresolved pairs**:
   - **n_ptc**: Comparator has a PTC (stop >50 nt upstream of last EJC)
   - **n_no_ptc**: Comparator has no PTC (stop in or near last exon)
   - **n_ptc_unknown**: Comparator lacks CDS annotation (no PTC call possible); excluded from PTC/no-PTC breakdown
6. **n_stop_last_exon**: Among no-PTC unresolved, stop codon resides in the last exon
7. **n_res_ptc**: Resolved frameshifts where comparator nevertheless has a PTC
8. **n_nofs_ptc**: No-frameshift pairs where comparator has a PTC

**Validation**: Funnel arithmetic is verified: `n_ptc + n_no_ptc + n_ptc_unknown = n_unresolved`.

The funnel is computed for both NMD (C2) and Control (C4) comparisons across all cell types. The key biological insight is the asymmetry: NMD unresolved frameshifts are enriched for PTCs (~46%), while Control unresolved frameshifts overwhelmingly lack PTCs (~97% no-PTC) with stop codons in the last exon.

### Novel Protein Sequence Quantification

Control (C4) pairs with unresolved frameshifts and no PTC represent productive alternative splicing: the frameshift rewrites the C-terminal protein sequence without triggering NMD. We quantify how much of the protein is novel (translated in a non-original reading frame).

#### Frameshift Boundary Determination

The first frameshift event per pair determines the boundary between conserved and novel protein sequence. Events are already sorted 5'→3' by `analyzeFrameWalk()`. The boundary position is:

- Plus strand: `genomic_end` of the first frameshift event
- Minus strand: `genomic_start` of the first frameshift event

This follows the same logic as Isopair's `last_shift_pos` computation. For pairs with multiple frameshifts and partial compensation, the first frameshift boundary is the correct choice — everything downstream is in a non-original reading frame.

#### Conserved CDS Base Pair Counting

The `count_conserved_cds_bp()` function walks the **comparator's** exon structure 5'→3', intersecting each exon with the comparator's CDS range, and accumulates CDS bp upstream of the frameshift boundary:

1. Get the comparator's CDS genomic range (`cds_start`, `cds_stop`) and exon coordinates
2. Order exons 5'→3' (ascending for plus strand, descending for minus strand)
3. For each exon, compute the CDS-overlapping segment: `max(exon_start, cds_start)` to `min(exon_end, cds_stop)`
4. Walk exons, accumulating CDS bp until the frameshift boundary is reached:
   - If the boundary is past the current exon's CDS segment: add the full segment
   - If the boundary falls within the segment: add the partial bp and stop
   - If the boundary is before the segment: stop (all remaining CDS is novel)
5. Cap conserved bp at `orf_length`

**Note**: The frame walk identifies frameshifts based on the *reference* CDS overlap, but conserved bp are measured using the *comparator's* CDS and exon structure. This is correct because we want the comparator's actual protein product. Comparators in scope are all coding (filtered by merge with PTC table).

#### Novel Protein Metrics

Per pair:
- `novel_bp = max(0, orf_length - conserved_bp)`
- `novel_aa = novel_bp / 3`
- `pct_novel = 100 × novel_aa / total_aa`

Summarized per cell type as median and IQR of `novel_aa` and `pct_novel`.

#### Edge Cases

- Comparator not in CDS table or structures → excluded (NA)
- Frameshift boundary outside comparator CDS range → `conserved_bp = orf_length`, `novel_bp = 0`
- Single-exon comparator → works naturally (one exon to walk)

#### Caching

Per-pair exon walking is cached via `cached_compute()` with key `"novel_protein_c4_{ct_label}"` since the computation is moderately expensive.

---

## Cross-Comparison Statistical Framework

### Overview

After characterizing splicing patterns within the pooled dataset, we compare NMD-triggering transitions to baseline splicing variation to determine whether NMD-associated isoform differences are qualitatively distinct from normal splicing diversity.

### Comparison Scope

**Comparisons analyzed:**
- **C1** (dominant NMD): Dominant non-NMD (DMSO) vs dominant NMD-sensitive (Smg1i)
- **C2** (top-CPM NMD): Top non-NMD by CPM vs top NMD-sensitive by CPM
- **C4** (baseline): Dominant non-NMD vs next-best non-NMD by CPM

**C3 excluded:** C4 pairs are a strict subset of C3. Using both introduces pseudo-replication because C3's multiple pairs per gene inflate sample sizes, and C3 and C4 are not statistically independent.

**Cell types analyzed (paper scope, 2026-04-29):** AT2, DD, FB, MV (4 non-ALI cell types + all_samples aggregate). DD_ALI, DO_ALI, and DO are excluded from the paper.

### Per-Comparison Filtering (Script 14)

For each comparison×run combination, the deduplicated splicing choice profiles are filtered to include only pairs present in that comparison's pairs file. The join matches on `(gene_id, dominant_isoform_id, non_dominant_isoform_id = comparator_isoform_id)`. Combinations with fewer than 50 matched profiles are skipped. Scripts 09-13 are then run on each filtered profile set.

### Statistical Design: Shared Isoform Structure

A key feature of our comparison framework is that C1/C2 and C4 share the same dominant isoform within each gene. This creates a natural pairing: for genes that appear in both an NMD comparison (C1 or C2) and the baseline (C4), we can perform paired tests that control for gene-level confounds. We use both paired and unpaired analyses:

- **Paired analysis** (gene overlap set, minimum 50 genes): controls for gene-specific effects since both NMD and baseline transitions originate from the same dominant isoform
- **Unpaired analysis** (full sets): uses all available data for maximum power

### Phase 1: NMD vs Baseline Tests

For each NMD comparison (C1, C2) against baseline (C4), within each run:

#### Test 1: Event Complexity

**Unpaired:** Wilcoxon rank-sum test comparing event counts (n_events) between NMD and baseline profiles. Effect size: Cliff's delta (nonparametric, range -1 to +1).

**Paired:** Wilcoxon signed-rank test on matched gene pairs + sign test (proportion of genes where NMD pair has more events than baseline pair).

#### Test 2: Profile Type Distribution

**Unpaired:** Chi-square test (or Monte Carlo simulation when expected counts < 5) comparing profile type frequencies (Terminal-only, Boundary-only, Inclusion-only, Partial_Retention-only, Full_Retention-only, Combined) between NMD and baseline. Effect size: Cramer's V.

**Paired:** Profile type transition matrix (proportion of genes changing profile type between NMD and baseline).

#### Test 3: Event Type Prevalence

**Unpaired:** Two-proportion z-test for each of 8 event types (Alt_TSS, Alt_TES, A5SS, A3SS, Partial_IR, IR, SE, Missing_Internal). Effect size: odds ratio with Haldane correction (+0.5 to all cells when any cell is zero). FDR correction across 8 tests within each comparison×run.

**Paired:** McNemar's test per event type on discordant gene pairs (exact binomial when discordant count < 10).

#### Test 4: Regional Enrichment Bootstrap

For comparisons where event_regions data is available (from Script 12):
- Resample profiles (not individual events) with replacement to preserve within-profile event correlation
- Recompute enrichment ratios (observed proportion / expected proportion based on region sizes) for each event type × region type
- 1,000 bootstrap iterations; report 95% CI for the NMD-minus-baseline enrichment difference
- Significance: CI excludes zero

### Phase 2: Random-Effects Meta-Analysis

Meta-analyze key effect sizes across the 5 cell types (AT2, DD, DD_ALI, FB, MV), excluding all_samples (not independent of individual cell types) and DO (excluded from downstream).

**Method:** Random-effects model (REML estimator) via the metafor R package. Random effects are preferred over fixed effects because cell types are biologically distinct populations, not replicate samples from a single population.

**Metrics meta-analyzed:**

1. **Cliff's delta** (event complexity): approximate SE from `sqrt((1 - delta^2) / (n_eff - 1))` where n_eff is the harmonic mean of group sizes
2. **Cramer's V** (profile type shift): meta-analyzed on raw scale with delta-method SE `sqrt(2/N)`
3. **Log-odds ratio** per event type: standard SE from `sqrt(1/a + 1/b + 1/c + 1/d)` with Haldane correction

**Heterogeneity:** I² statistic (proportion of variance due to between-study heterogeneity) and Cochran's Q test.

**Output:** Forest plots for each metric; pooled estimates with 95% CI.

### Phase 3: Sensitivity Analyses

#### 3.1: C1 vs C2 Concordance

Tests whether the two NMD comparison definitions yield similar results. For each run, computes Pearson correlation between C1 and C2 event prevalence differences (NMD - baseline). High concordance supports the robustness of findings.

#### 3.2: Complexity Confound Check

Tests whether NMD comparator isoforms are structurally different from baseline comparator isoforms. Compares n_exons, transcript length, and n_junctions of the comparator (non-dominant) isoform between NMD and baseline groups. If NMD comparators are systematically more/less complex, this confounds event count comparisons.

#### 3.3: Event Direction Analysis

Compares the proportion of LOSS vs GAIN events between NMD and baseline transitions. A shift toward more LOSS events in NMD would indicate that NMD-sensitive isoforms are more often "reduced" versions of the dominant (missing exons), while more GAIN events would indicate they carry additional sequence (retained introns, extra exons).

### Multiple Testing Correction

- **Within test families:** FDR (Benjamini-Hochberg) across 8 event types within each comparison×run×analysis_type
- **Across test families:** Results reported per family without cross-family adjustment (different biological questions)
- **C1 results flagged:** Lower confidence where sample sizes are small (< 100 profiles)

---

## Implementation Details

### GTF Parsing and Data Preparation

#### Critical Implementation Details

**GTF Feature Type Filtering**:

The GTF/GFF3 format includes multiple feature types in column 3:
- `transcript`: Full transcript span (single line per transcript)
- `exon`: Individual exons
- `CDS`: Coding sequence regions
- `UTR`: Untranslated regions

**CRITICAL**: Only `exon` features should be parsed for exon structure analysis. Including `transcript` features causes incorrect exon counts and structure misidentification.

**Implementation**:
```r
parse_gtf <- function(file) {
  for (line in lines) {
    parts <- strsplit(line, '\t')[[1]]

    # Only parse exon features, skip transcript/gene/CDS features
    if (parts[3] != "exon") {
      next
    }

    # ... parse exon attributes
  }
}
```

**Transcriptional Exon Ordering**:

**Issue**: GTF files may not include `exon_number` attributes (e.g., SQANTI output), causing exons to be ordered by genomic coordinates rather than transcriptional order. This breaks TSS/TES detection on minus strand genes.

**Solution**: Automatically assign transcriptional exon numbers when missing:

```r
# Assign transcriptional exon numbers if missing
gtf_data <- gtf_data %>%
  group_by(transcript_id, strand) %>%
  mutate(
    exon_number = if (all(is.na(exon_number))) {
      if (first(strand) == "+") {
        rank(start, ties.method = "first")      # Ascending for plus strand
      } else {
        rank(desc(start), ties.method = "first")  # Descending for minus strand
      }
    } else {
      exon_number
    }
  ) %>%
  ungroup()
```

**Transcriptional Order Convention**:
- **Exon 1** = TSS exon (first exon in transcriptional 5' → 3' order)
  - Plus strand: exon 1 = lowest genomic coordinates
  - Minus strand: exon 1 = highest genomic coordinates
- **Exon N** = TES exon (last exon in transcriptional order)
  - Plus strand: exon N = highest genomic coordinates
  - Minus strand: exon N = lowest genomic coordinates

This ensures that:
- `exons[1, ]` always refers to the TSS exon (for both strands)
- `exons[nrow(exons), ]` always refers to the TES exon (for both strands)
- TSS/TES detection functions work correctly without strand-specific logic

**Impact**: This fix resolved incorrect TSS/TES calls on minus strand genes in real data validation (tested on 25 PacBio gene pairs, 100% accuracy achieved).

---

### Data Structures

#### Isoform Structure Table
```
Columns: isoform_id, gene_id, seqnames, strand, exon_number, exon_start, exon_end
Format: One row per exon
Ordering: Exons sorted by exon_number (1, 2, 3, ...)
```

#### Union Exon Table
```
Columns: gene_id, union_exon_id, union_exon_index, seqnames, start, end, strand
Format: One row per union exon
Ordering: union_exon_index reflects 5' → 3' transcript order
```

#### Isoform-Union Exon Mapping
```
Columns: isoform_id, gene_id, union_exon_id, union_exon_index, present (1/0), region_type
Format: One row per (isoform, union_exon) pair
Region Types: "5'UTR", "CDS", "3'UTR", "contains_orf_start", "contains_orf_stop", "contains_orf_start_stop", "non_coding", "unknown"
```

#### Splicing Choice Profile
```
Columns: gene_id, dominant_isoform_id, non_dominant_isoform_id,
         n_union_exons_total, n_exons_shared, n_exons_dominant_only, n_exons_non_dominant_only,
         tss_changed, tes_changed,
         n_a5ss, n_a3ss, n_partial_ir, n_ir, n_se,
         complexity metrics, spatial metrics, ...
Format: One row per (dominant, non-dominant) comparison
```

### Computational Performance

**Union Exon Construction**: O(n log n) per gene, where n = total exons across all isoforms
**Event Detection**: O(m²) per gene, where m = exons per isoform (typically small)
**Total Runtime**: ~10-20 minutes for ~23,000 genes (human transcriptome)

### Quality Control

**Validation Checkpoints**:
1. After union exon construction: Check coverage (all exonic bases assigned)
2. After event detection: Run synthetic validation (100% pass required)
3. After profile generation: Check distribution of event frequencies

**Error Handling**:
- Skip genes with incomplete annotations
- Flag isoforms with unknown strand
- Report genes with unexpected union exon patterns

---

## Software Implementation

**Language**: R (version 4.3+)

**Key Packages**:
- `tidyverse`: Data manipulation and visualization
- `GenomicRanges`: Genomic coordinate operations (optional)
- `glmnet`: LASSO/ridge regression for co-occurrence analysis
- `metafor`: Random-effects meta-analysis (Phase 2 cross-comparison)

**Scripts** (Version 6.0) — organized into `scripts/core/` (generic) and `scripts/nmd/` (NMD study-specific):

**NMD-specific (`scripts/nmd/`):**
1. `nmd/01_prepare_expression_data.R`: Load DGEList, calculate CPM, identify dominant isoforms
6. `nmd/06_filter_to_analysis_subset.R`: Apply expression (filterByExpr) and gene category filters
7. `nmd/07_classify_and_pair.R`: NMD/non-NMD classification + comparison pair generation (C1-C4, 4 comparisons x 7 runs = 28 sets, with deduplication)
9. `nmd/09_analyze_complexity_relationship.R`: Complexity vs event count, quartile binning
10. `nmd/10_analyze_cooccurrence.R`: Event co-occurrence (Fisher's exact, LASSO-controlled)
11. `nmd/11_analyze_spatial_patterns.R`: Positional bias, topology enrichment, proximity
12. `nmd/12_analyze_functional_context.R`: Regional distribution, ORF boundary susceptibility, ORF impact; saves event_regions.rds for cross-comparison bootstrap
13. `nmd/13_analyze_patterns.R`: Profile type classification, pattern frequencies by complexity
14. `nmd/14_run_per_comparison.R`: Filter profiles per comparison×run, run Scripts 09-13 on each set (C1, C2, C4 × 6 runs; excludes C3 and DO)
15. `nmd/15_compare_across_comparisons.R`: Cross-comparison statistical framework — NMD vs baseline tests (paired + unpaired), random-effects meta-analysis, sensitivity analyses

**Core pipeline (`scripts/core/`):**
2. `core/02_extract_isoform_structures.R`: Parse GFF files to exon structures
3. `core/03_build_union_exons.R`: Construct atomic union exon models
4. `core/04_extract_cds_annotations.R`: Extract CDS coordinates
5. `core/05_annotate_region_types.R`: Classify union exons by region
8. `core/08_extract_splicing_profiles.R`: Hierarchical event detection, generate profiles
   - `--reconstruction_check`: On-the-fly reconstruction verification (eliminates need for separate core/09 run)
   - `--pairs-file <path>`: Specify explicit contrast pairs (TSV: gene_id, dominant_isoform_id, comparator_isoform_id)
   - `--output <path>`: Custom output path for splicing profiles RDS
   - `--test N`: Limit to first N genes
9. `core/09_validate_reconstruction.R`: Standalone batch reconstruction validation

**Shared Libraries** (`scripts/core/`):
- `event_detection_functions.R`: Core event detection functions (19 functions), detection thresholds (TSS_TOLERANCE, TES_TOLERANCE, SPLICE_SITE_THRESHOLD)
- `reconstruction_functions.R`: Reconstruction from comparator + events (`reconstruct_dominant_v2`), verification (`verify_transcript`)
- `visualization_functions.R`: Isoform pair plotting, event annotation brackets, reconstruction mismatch highlighting

**Testing** (`scripts/tests/`):
- `run_tests.R`: Automated test runner (44 synthetic + 80 real-data + 2 dedup = 126 cases, 100% PASS)
- `extract_failure_cases.R`: Extract test cases from pipeline failures

**Documentation**:
- `SPLICE_FUNCTION_CATALOG.md`: Comprehensive function catalog linking to splicing biology
- `EVENT_DETECTION_FLOW.md`: Mermaid flowchart of data flow through event detection

---

## Publication Report Methods (mashr Classification)

### Sample Exclusion

DO donor 029T was identified as a PCA outlier (distance 252 from DO centroid,
>3× the next most distant DO sample) and excluded along with its treatment
partner, reducing the dataset from 38 to 36 samples (4 DO samples from
donors 001V and 027U).

### NMD Classification (mashr)

Isoforms are classified using mashr (multivariate adaptive shrinkage)
differential isoform expression results:
- **NMD-sensitive**: `nmd_responsive == TRUE` (pre-computed by mashr;
  lfsr < 0.05 & posterior logFC > 0 with mashr-shrunken effect sizes)
- **Non-NMD**: `adj.P.Val > 0.30`
- **all_samples**: NMD = union across AT, DD, FB, MV cell types;
  non-NMD = intersection across the same 4 cell types.

The non-NMD threshold of 0.30 was selected (2026-04-29) to restore the
historical proportion of non-NMD isoforms (~52–60% of expressed isoforms
across the four cell types) after the mashr model was refit on the
4-cell-type paper scope. The earlier 0.50 threshold, while appropriate
when the mashr model was fit on six cell types, became overly restrictive
under the 4-CT refit because the per-CT limma fits feeding mashr were
rerun on the smaller subset and produced tighter adj.P.Val distributions.
At 0.50 only 7–13% of isoforms per CT qualified as non-NMD; at 0.30 the
per-CT proportions are AT 60.1%, DD 52.4%, FB 59.6%, MV 58.1%.

Paper scope (2026-04-29): the analysis is restricted to four non-ALI
primary lung cell types — AT, DD, FB, MV. DD_ALI, DO_ALI, and DO are
not used.

### Pair Construction and Gene-Matching

- **C4 (Control)**: Top two non-NMD isoforms by DMSO CPM per gene
  (`generatePairsExpression(..., method = "top_two")`)
- **C2 (NMD)**: Same C4 reference paired with top NMD isoform by Smg1i CPM
  (via `identifyDominantIsoforms()` + inner_join on `gene_id`)
- **Gene-matching**: C2 and C4 restricted to shared
  `(gene_id, reference_isoform_id)` pairs, ensuring the same genes and
  reference isoforms appear in both NMD and Control analyses

### Transcriptional Diversity Analysis (Section 1)

For genes in the gene-matched pair set, isoform diversity is characterized
using all expressed isoforms (coding and non-coding) from `structures.rds`
and `expression_data.rds`. Metrics include:

- **Isoform count distribution**: Number of expressed isoforms per gene.
- **TSS/TES diversity**: Number of distinct strand-aware TSS and TES
  positions per gene (multi-isoform genes only).
- **TSS-TES independence**: Within-gene Spearman correlation between TSS
  and TES positions across isoforms (genes with ≥3 isoforms). Combination
  usage: percentage of possible TSS×TES pairs actually observed.
- **TSS-CDS and TES-CDS coupling**: Per-gene Spearman correlation between
  TSS and CDS 5' boundary (strand-aware), and TES vs CDS 3' boundary
  (coding isoforms only, genes with ≥3 coding isoforms).
- **CDS start diversity**: Number of distinct CDS start positions per gene,
  and whether starts fall in different reading frames (positions differing
  by a non-multiple of 3).
- **Expression concentration**: Fraction of gene-level DMSO CPM contributed
  by the dominant isoform (highest mean CPM across DMSO samples).

### Combined Prediction Model: TD2 vs Reference CDS (Section 3d)

#### Overview

Three model sets compare CDS annotation strategies for predicting NMD status
from isoform structural features. All models are trained and evaluated on the
same matched population: the intersection of mashr-classified NMD and non-NMD
coding isoforms with complete features from both CDS sources. Holdout
chromosomes 1, 3, 5, 7 are never used for feature selection or model tuning.

Source scripts: `05l_unified_model.R` (TD2 features), `05t_ref_cds_features.R`
(reference-CDS features), `05u_paralog_annotation.R` (paralog annotation),
`05v_model_comparison.R` (model fitting) → `analysis_cache/model_comparison.rds`.

#### Reference-CDS Feature Construction

For each gene, the dominant non-NMD coding isoform is identified as the
non-NMD isoform with the highest mean DMSO CPM. This isoform's CDS start
(ATG) position serves as the reference reading frame anchor for all other
isoforms in the gene.

**Algorithm** (implemented in `05t_ref_cds_features.R`):

1. **Dominant isoform identification**: For each gene, identify all non-NMD
   coding isoforms present in the CPM matrix. Compute mean CPM across all
   DMSO samples. Select the isoform with the highest mean DMSO CPM as the
   reference.

2. **ATG tracing**: For each target isoform, check if the reference's
   strand-aware CDS start ATG is exonic (all 3 nucleotides) in the target's
   exon structure. If exonic, map the ATG to transcript-space and verify
   the codon is ATG in the target's sequence.

3. **ORF walking**: From the mapped ATG, walk the target's transcript
   sequence in 3-nucleotide steps until the first in-frame stop codon.
   Record the ORF length and stop codon position.

4. **Downstream EJC count**: Compute exon-exon junction positions in the
   target's transcript-space. Count junctions downstream of the stop codon's
   last nucleotide. Truncate at 5 (values ≥5 set to 5).

5. **5'UTR features**: Construct synthetic CDS metadata using the traced ATG
   and stop codon positions (mapped back to genomic coordinates via
   `transcript_to_genomic()`). Call `Isopair::scan5UtrFeatures()` with the
   target's exon structure, synthetic CDS, and transcript sequence. This
   produces ATG density, strong Kozak count, orphan ATG count, uORF counts
   (overlapping, in-frame, out-of-frame), longest uORF, ORF coverage, and
   stop density — all computed relative to the reference CDS boundaries.

6. **3'UTR length**: Exonic basepairs from the traced stop codon to the 3'
   transcript end, log-transformed (log1p). Uses the same strand-aware
   formula as the TD2 3'UTR computation.

**Coordinate conventions**: `cds_start < cds_stop` (genomic order). For +
strand: ATG at `cds_start`, stop codon 3' end at `cds_stop`. For - strand:
stop codon 3' end at `cds_start`, ATG at `cds_stop`. The
`transcript_to_genomic()` function handles strand-aware mapping of the stop
codon's last transcript-space position back to the correct genomic coordinate.

**Edge cases**: Isoforms where the reference ATG is not exonic receive NA for
all reference-CDS features. Isoforms that ARE the reference receive consistent
features (self-reference validation). Genes with no non-NMD coding isoform
have no reference and receive NA.

Output: `analysis_cache/ref_cds_features_all.rds`.

#### Paralog Annotation and Test Set Filtering

Genes with high-similarity paralogs on opposite sides of the train/test
chromosome split pose a data leakage risk. Paralog pairs are queried from
Ensembl via biomaRt (`ensembl_gene_id`, `hsapiens_paralog_ensembl_gene`,
`hsapiens_paralog_perc_id`). Three filters are applied sequentially:

1. **Both expressed**: Gene and paralog must both be in the classified
   expression dataset.
2. **High similarity**: Paralog protein sequence identity ≥ 80%.
3. **Cross-split**: Gene and paralog on opposite sides of the holdout
   chromosome boundary (one on chr 1/3/5/7, one on a training chromosome).

Only holdout-side genes from cross-split pairs are removed from the test set.
Training data is not modified.

Source script: `05u_paralog_annotation.R` →
`analysis_cache/paralog_genes.rds`.

#### TD2 Feature Construction

TD2 features are computed from each isoform's TransDecoder2-predicted CDS
boundaries, as in the original unified model:

- **Downstream EJC count**: From `ptc.rds` (`n_downstream_ejcs`, truncated
  at 5).
- **5'UTR features**: From `utr5_features_all.rds` (ATG density/count,
  strong Kozak, orphan ATG, uORF counts, ORF coverage, stop density).
- **3'UTR length**: Exonic bp from `cds_stop` to 3' transcript end (log1p).

Source script: `05l_unified_model.R` → `analysis_cache/unified_model.rds`.

#### Model Sets

All three model sets use elastic net logistic regression (glmnet, alpha = 0.5,
10-fold cross-validation, `type.measure = "auc"`). Step 1 (single feature)
uses standard logistic regression. Features are centered and scaled using
training set statistics; the same scaling parameters are applied to the test
set. `set.seed(42)` is called before each `cv.glmnet` fit.

**Model Set 1 — TD2 CDS** (progressive):
- Step 1: `downstream_ejc` (1 feature)
- Step 2: + `atg_density`, `atg_count`, `atg_strong_kozak` (4 features)
- Step 3: + `uorf_count_overlapping`, `uorf_longest_nt`, `uorf_count_inframe`,
  `uorf_count_outframe`, `utr5_orf_coverage`, `stop_density`,
  `atg_orphan_count` (11 features)
- Step 4: + `log_utr3_length` (12 features)

**Model Set 2 — Reference CDS** (progressive, same structure):
- Steps 1–4 use the `ref_` prefixed equivalents of the same features.

**Model Set 3 — Combined**:
- All 24 features (12 TD2 + 12 reference CDS) provided to a single elastic
  net. The L1 penalty selects which features survive, revealing complementarity
  vs redundancy between CDS sources.

AUC is evaluated on the paralog-free holdout test set for all models.

#### SHAP-Based Isoform Clustering

SHAP (SHapley Additive exPlanations) values decompose each isoform's
predicted log-odds into per-feature contributions. For a linear model
(elastic net logistic regression), SHAP values have a closed-form solution:

    SHAP_j(i) = beta_j × x_j,scaled(i)

where `beta_j` is the standardized coefficient at `lambda.min` and
`x_j,scaled(i)` is the centered/scaled feature value for isoform i. Features
with zero coefficients have zero SHAP values.

**All-isoform clustering**: K-means clustering is applied to the SHAP matrix
(non-zero features only) of the full matched population. Optimal k is selected
by silhouette score (sampled to 5,000 isoforms for computational efficiency)
over k = 2..8.

**NMD-only clustering**: The same procedure is repeated on NMD isoforms only,
to identify mechanistically distinct subpopulations within the NMD class.
Cluster profiles are characterized by mean SHAP values, downstream EJC
patterns from both CDS sources, and predicted probabilities.

Visualization: cluster-mean SHAP heatmaps and beeswarm summary plots.

#### Dose-Response Analysis

The combined model's predicted NMD probability is tested for correlation with
actual NMD response magnitude. Raw limma logFC values (per-cell-type
unshrunk estimates from 6 cell types: AT, DD, DD_ALI, DO, FB, MV) are
averaged per isoform to produce a mean logFC. Spearman correlation between
predicted probability and mean logFC is computed for holdout NMD isoforms.

### uORF ATG Position Sharing (Section 3)

For genes with PTC-negative NMD comparators that have overlapping uORFs,
all expressed coding isoforms of those genes are examined to determine
whether the uORF ATG position also serves as a CDS start in other isoforms.
CDS start positions are extracted from `cds.rds` (strand-aware: `cds_start`
on + strand, `cds_stop` on - strand). ATG positions shared between NMD and
non-NMD isoforms are identified by matching CDS 5' boundaries across all
coding isoforms of the gene.

### PTC-Causing Event Attribution

For each PTC+ comparator, the specific splice event that caused the PTC is
identified through two mechanisms:

1. **Frameshift PTCs**: The frameshift-causing event identified by
   `analyzeFrameWalk()` — the splice event that shifted the reading frame,
   leading to a premature stop codon in the new frame.
2. **Non-frameshift (in-frame stop) PTCs**: The event whose genomic
   coordinates contain the strand-aware PTC position. Since the isoform is
   in-frame, the stop codon must reside within the sequence that differs
   between reference and comparator.

**Strand-aware PTC position**: The genomic coordinate of the stop codon
depends on the gene's strand. For + strand genes, the stop codon is at
`cds_stop` (the larger genomic coordinate). For - strand genes, the stop
codon is at `cds_start` (the smaller genomic coordinate, which corresponds
to the 3' end in transcription direction). This mirrors the strand handling
in `Isopair::computePtcStatus`.

Pairs where the PTC position does not fall within any event's coordinates
are labeled "unresolved" and excluded from the Sankey diagram and enrichment
analyses. Similarly, pairs classified as frameshift by `compareIsoformFrames()`
but lacking an `is_frameshift` event in `analyzeFrameWalk()` are labeled
unresolved (classification disagreement). A diagnostic table characterizes
all unresolved cases by mechanism type and distance to the nearest event
boundary.

This attribution enables a comprehensive Sankey diagram tracing PTC-causing
splice events to their mechanism (frameshift vs in-frame stop).

**Implementation**: The attribution logic is implemented in
`attribute_ptc_events()` (defined in `analysis_functions.R`), called by both
the report (Section 2b, for PTC-positive pairs) and the precompute script
(`05r_ref_atg_analysis.R`, for reclassified PTC-negative pairs). The function
takes configurable parameters: `atg_genomic_pos` (NULL for SQANTI CDS
attribution, populated for reference-ATG attribution) and `is_frameshift_vec`
(pre-computed from `compareIsoformFrames()` frame categories, or inferred from
frame walk events).

**Split-codon splice junctions**: A stop codon can span a splice junction,
with 2 nucleotides from one exon and 1 from the adjacent exon. The event that
creates the junction is technically just outside the stop codon's first-
nucleotide position. To handle this, coordinate containment checks use a ±2 bp
buffer (`stop_g >= ev_min - 2 & stop_g <= ev_max + 2`), and the ATG-to-stop
region filter uses overlap semantics with a 3 bp extension rather than strict
containment.

### Reference-ATG Tracing (PTC-Negative Reclassification)

For PTC-negative NMD pairs where both isoforms are coding, we test whether the
comparator is effectively PTC-positive when analyzed from the reference isoform's
CDS start position rather than the TransDecoder2-predicted CDS.

**Motivation**: TransDecoder2 (TD2) selects the CDS with the highest PSAURON
score and longest significant ORF. When a splice event creates a frameshift that
truncates the reference reading frame, TD2 may select a longer alternative ORF
whose stop codon has no downstream EJCs — making the isoform appear PTC-negative.
For ENST isoforms, TD2 independently predicts the CDS (92.2% agreement with
GENCODE at the CDS 5' position). For novel isoforms (76.8% of NMD comparators),
the CDS is entirely TD2-predicted. CDS prediction provenance verified from
SQANTI3 source code: `helpers.py:predictORF()` calls `run_td2()` (not
GeneMarkS-T).

**Algorithm** (implemented in `05r_ref_atg_analysis.R`):

1. **ATG availability check**: For each gene-matched pair, determine the
   reference isoform's strand-aware CDS start position (+ strand: `cds_start`;
   - strand: `cds_stop`). Check whether the full 3-nucleotide ATG codon is
   exonic in the comparator's exon structure.

2. **ORF tracing**: If the ATG is available, map it to the comparator's
   transcript-space position. Walk the comparator's transcript sequence from
   the ATG in 3-nucleotide steps until the first in-frame stop codon
   (TAA/TAG/TGA). Record the ORF length.

3. **Downstream EJC counting**: Compute exon-exon junction positions in the
   comparator's transcript-space. Count junctions >50 nt downstream of the
   stop codon (the NMD rule).

4. **Classification**: Pairs are classified as:
   - `effectively_ptc`: ORF is shorter than the reference AND ≥1 downstream EJC
   - `truncated_no_ejc`: ORF is shorter but no downstream EJC
   - `ref_atg_lost`: Reference ATG codon not exonic in comparator
   - `no_downstream_ejc`: ORF same or longer, no downstream EJC
   - `same_or_longer_with_ejc`: ORF same or longer with downstream EJCs

5. **PTC-causing event attribution**: For `effectively_ptc` pairs, the
   PTC-causing splice event is identified using the same `attribute_ptc_events()`
   function as for the original PTC-positive pairs (see "PTC-Causing Event
   Attribution" above). The `atg_genomic_pos` parameter restricts attribution
   to events overlapping the ATG-to-premature-stop genomic region.

**Output**: `analysis_cache/ref_atg_analysis.rds` containing per-pair
classification and attribution results for both C2 (NMD) and C4 (Control) pairs.

### 3'UTR Length Analysis with PTC Correction

Standard 3'UTR length measurement (`utr3_bp_comp` from
`quantifyPairDivergence()`) uses each isoform's own `cds_stop`. For PTC+
isoforms, this is the premature stop, so the measured 3'UTR is inflated
(includes former CDS sequence downstream of the PTC).

To distinguish genuine 3'UTR length from PTC-induced inflation, PTC+
comparators are split by stop codon mechanism:

- **Same-stop PTC+**: Comparator and reference share the same `cds_stop`.
  The PTC arises because a 3'UTR splicing event repositioned the last EJC
  downstream. The 3'UTR length measurement is correct and unaffected.
- **Diff-stop PTC+**: Comparator has a premature stop (CDS frameshift or
  in-frame stop). The measured 3'UTR is inflated. A corrected
  "reference-based 3'UTR" is computed: exonic bp downstream of the
  reference's `cds_stop` in the comparator's exon structure. This is only
  valid when the reference stop position falls within the comparator's exons
  (~86% of cases).

### 5'UTR Feature Analysis and uORF Detection

#### Phase 0: Data-Driven Feature Scan

A comprehensive scan of 28 5'UTR sequence features is performed for all coding
comparator and reference isoforms using `Isopair::scan5UtrFeatures()`. Features
include ATG and stop codon counts (total and stratified by reading frame relative
to the main CDS ATG), ORF counts (in-frame, out-of-frame, overlapping), density
measures (per 100 bp), Kozak initiation context, and spatial features. The
fraction of 5'UTR occupied by uORFs (`pct_utr5_in_orfs`) is computed as the
union of all uORF intervals within the 5'UTR, with no double-counting of
overlapping ORF regions.

Feature importance is assessed by elastic net logistic regression (glmnet,
alpha=0.5, 10-fold cross-validation) predicting PTC-negative NMD status vs
Control. Three models are fit: Model A (PTC_neg vs Control, primary), Model B
(PTC_pos vs Control, specificity control), Model C (PTC_neg vs PTC_pos).
Cross-validated AUC is computed with `type.measure="auc"`.

#### uORF Detection

1. **5'UTR extraction**: For each isoform, the 5'UTR boundaries are computed
   from genomic CDS coordinates + exon structure (strand-aware). Validated by
   checking for ATG at the computed CDS start position (98.6% pass rate).
2. **uORF detection**: All AUG positions in the 5'UTR sequence (from SQANTI
   corrected transcript FASTA) are identified via `Isopair::detectUorfs()`.
   For each ATG, the reading frame is walked through the **full transcript
   sequence** until the first in-frame stop codon (TAA, TAG, TGA). Minimum
   uORF size: 3 codons (9 nt). Walking through the full transcript (not just
   the 5'UTR) is necessary to detect overlapping uORFs whose stop codon falls
   in or past the main CDS.
3. **Overlapping uORFs**: uORFs whose stop codon falls at or beyond the main
   CDS start are classified as overlapping — a known NMD trigger.
4. **Validation against ORFik**: Non-overlapping uORF counts are compared to
   ORFik `findORFs(startCodon="ATG", longestORF=FALSE, minimumLength=0)`
   applied to 5'UTR-only sequences. ORFik's `minimumLength` counts body
   codons excluding start and stop; `minimumLength=0` matches our 3-codon
   total minimum. ORFik cannot detect overlapping uORFs from 5'UTR-only
   input because their stop codons fall outside the provided sequence.
5. **Groups compared**: PTC-negative NMD comparators, PTC-positive NMD
   comparators, and Control (C4) comparators, all restricted to coding pairs.

#### uORF Set Operations

For each isoform pair, uORFs in the reference and comparator are compared by
genomic ATG position (`Isopair::compareUorfs()`). uORFs sharing the same ATG
position are classified as "shared"; those present only in the comparator are
"gained"; those only in the reference are "lost". Gained and lost uORFs are
stratified by overlapping status.

#### Splice Event Co-Occurrence

For pairs where the comparator gains overlapping uORFs, the prevalence of each
splice event type (from `Isopair::buildProfiles()`) is tabulated and compared to
pairs without gained overlapping uORFs using Fisher's exact test. This describes
the structural context of overlapping uORF gain without claiming individual
causal attribution.

#### NMD Response Correlation

The elastic net Model A predicted probability (from 5'UTR features) is tested
for correlation with NMD response magnitude (mashr posterior mean logFC). Each
isoform has 6 cell-type-specific logFC values (from `ashr::get_pm()` on the
mashr model). A linear mixed-effects model (`lme4::lmer()`) is fit within each
group: `logFC ~ pred_nmd_prob + (1 | isoform_id)`, with the isoform random
intercept accounting for repeated cell-type measurements.

### 3'UTR Splicing Enrichment

For each comparator isoform with CDS annotation, splice events in the 3'UTR
are identified by checking whether any event's genomic coordinates fall
downstream of the strand-aware stop codon position (`cds_stop` on + strand,
`cds_start` on - strand) in the transcription direction.
The rate of 3'UTR splicing is compared across PTC-negative NMD, PTC-positive
NMD, and Control groups using Fisher's exact test.

### Cumulative NMD Mechanism Attribution

A three-tier framework summarizes the evidence for NMD mechanisms across all
coding gene-matched NMD pairs:

1. **Tier 1 — PTC-mediated (mechanistic):** Pairs where the comparator has a
   PTC (stop >50 nt upstream of last EJC). The specific causal splice event
   is identified as described in "PTC-Causing Event Attribution" above.
   Subdivided by stop codon origin:
   - *Diff-stop frameshift*: Splice event shifts reading frame → new stop codon
   - *Diff-stop in-frame stop*: Splice event introduces in-frame stop codon
   - *Same-stop / 3'UTR splice*: Comparator and reference share the same stop
     codon. A 3'UTR splicing event repositions an EJC >50 nt downstream of the
     shared stop, triggering NMD via the canonical EJC model.

2. **Tier 2 — uORF-associated (associative):** PTC-negative NMD comparators
   with detectable 5'UTR uORFs. Group-level enrichment is significant
   (Fisher's exact test), but individual causal attribution is not possible
   from sequence data alone. Overlapping uORFs (stop codon at or past main
   CDS start) are distinguished as a stronger NMD trigger.

3. **Tier 3 — Unexplained:** PTC-negative pairs without detectable uORFs.
   Possible mechanisms include long 3'UTR-mediated NMD, EJC-independent UPF1
   recruitment, or structural features not captured by this analysis.

### Protein-Level Consequence Analysis (Goal 5)

Goal 5 uses **all non-NMD C4 isoform comparisons** (unmatched — every gene with
≥2 non-NMD isoforms), deduplicated across cell types. This broader scope maximizes
statistical power for domain and mass-spec analyses, unlike Goals 1-4 which use
gene-matched NMD vs Control pairs.

#### Splice-to-Protein Mapping

`Isopair::mapSpliceToProtein(profiles, cds, structures, buffer_aa = 5)` maps
each genomic splice event to protein coordinates:

1. For each CDS-overlapping event, the genomic coordinates (`five_prime`,
   `three_prime`) are intersected with the isoform's CDS boundaries.
2. The overlapping genomic region is converted to CDS-relative nucleotide
   positions, then to amino acid positions (`protein_start_aa`,
   `protein_end_aa`).
3. A ±5 aa buffer (`buffer_aa`) is added to capture junction-spanning
   peptides, stored as `protein_start_aa_buffered` / `protein_end_aa_buffered`.
4. Events entirely in UTR regions are excluded.

Output: one row per CDS-overlapping event, with protein coordinates and a
flag (`event_affects_frame`) indicating whether the event causes a frameshift.

#### Domain Detection (hmmscan)

Protein domains are annotated using HMMER hmmscan against Pfam-A:

1. Protein sequences for all coding C4 isoforms are extracted from SQANTI
   corrected protein FASTA and written to `c4_coding_proteins.faa`.
2. hmmscan is run against the Pfam-A HMM database with default thresholds.
3. The domtblout output is parsed, extracting: domain name, Pfam accession,
   isoform ID, envelope coordinates (`env_from`, `env_to` in amino acid
   positions), and per-domain E-value.
4. Domains are mapped per-isoform (not per-gene), enabling isoform-specific
   domain annotation. This replaces the earlier biomaRt approach which only
   provided gene-level canonical domain positions.

#### Domain Effect Classification

For each CDS-overlapping splice event, domain overlap is assessed:

1. **Reference domains**: All Pfam domains on the reference isoform are
   checked for overlap with the event's protein region. Overlapping domains
   are classified as `disrupted_lost` — the splice event disrupts or removes
   this domain in the comparator.
2. **Comparator domains (GAIN events only)**: For GAIN-direction events,
   domains on the comparator isoform overlapping the gained region are
   classified as `gained_intact` — the comparator has a domain in a region
   absent from the reference.

#### Domain Enrichment Testing

Over-representation of specific Pfam domains among disrupted domains is
tested using Fisher's exact test (`Isopair::testDomainEnrichment()`):

- **Numerator**: Number of comparisons where this domain is disrupted.
- **Denominator**: Number of reference isoforms carrying this domain
  (the baseline prevalence).
- **2×2 table**: (disrupted / not-disrupted) × (this domain / all other domains).
- **Multiple testing**: Benjamini-Hochberg FDR correction across all
  tested domains. Significance threshold: FDR < 0.05.

#### Directional Domain Enrichment

Domain effects are further stratified by splice event direction (GAIN vs LOSS):

1. **Event-level gating**: Each CDS-overlapping event is classified as
   overlapping ≥1 Pfam domain or not. A Fisher's exact test compares domain
   overlap rates between GAIN and LOSS events (2×2: direction × domain overlap).

2. **Gained domain enrichment** (Table 9b): Analogous to the disruption
   enrichment test above, but using `gained_intact` effects as the numerator.
   The denominator uses comparator isoform domain prevalence (`n_comp`),
   computed directly from `parsed_doms` (per-isoform hmmscan annotations).
   This is the correct baseline: gained domains are carried by comparator
   isoforms, so the prevalence should reflect how many comparators have each
   domain.

3. **Disrupted vs gained comparison** (Table 9c): For each domain with ≥5
   affected pairs (disrupted + gained), a one-sided binomial test assesses
   whether the gained fraction exceeds the overall baseline gained rate.
   The baseline is computed as the fraction of all domain-affected pairs
   with any gained-intact effect. BH correction is applied across all tested
   domains. Domains with significantly elevated gained fractions represent
   protein functions that alternative splicing preferentially introduces
   rather than ablates.

#### Cell-Type Domain Analysis

Domain overlap rates are stratified by cell type using `ct_membership`,
which maps each deduplicated isoform pair to the cell types in which it
was observed. Because pairs are not exclusive to a single cell type (a pair
present in 3 CTs contributes to all 3), per-CT counts are not independent.

For each cell type, we report: number of coding pairs, pairs with CDS events,
pairs with ≥1 domain affected, and the disrupted/gained breakdown. A
Cochran-Mantel-Haenszel (CMH) test assesses the direction × domain overlap
association stratified by cell type. The CMH assumes independent strata; the
shared-pair non-independence means the p-value should be interpreted
conservatively.

Per-domain cell-type heterogeneity is tested using chi-squared tests with
Monte Carlo simulation (B = 10,000). For each Pfam domain with ≥10 total
affected pairs across ≥2 cell types, a 2×5 contingency table is constructed:
rows = (domain affected / not affected), columns = 5 cell types. The
denominator per cell type is all pairs with CDS-overlapping events. The null
hypothesis is that the rate of domain X being affected is the same across all
cell types. Standardized residuals identify which cell types drive
heterogeneity (positive = enriched, negative = depleted). BH correction is
applied across all tested domains. The non-independence caveat from shared
pairs applies; anti-conservative p-values should be interpreted as hypothesis-
generating.

#### Mass-Spectrometry Integration with Domain Results

Domain enrichment tables (Tables 9, 9b) are annotated with mass-spectrometry
support columns: for each domain, `n_ms_supported` counts how many
disrupted/gained pairs also have ≥1 PeptideAtlas-confirmed splice-specific
peptide, and `pct_ms` gives the percentage.

A domain × mass-spec cross-tabulation (Table 10) tests whether
domain-affecting splice events are more or less likely to have mass-spec
support than non-domain-affecting events, using a 2×2 Fisher's exact test
at the pair level (universe = pairs with CDS-overlapping events).

#### PeptideAtlas Mass Spectrometry Validation

Splice-specific peptides are validated against the Human PeptideAtlas
(2024-01 build, ~1.5M peptides):

1. **Peptide generation**: For each comparison with CDS-overlapping events,
   the event-affected protein region (with ±5 aa buffer) is extracted from
   the comparator's protein sequence. Tryptic digestion (cleave after K/R,
   not before P) produces peptides of 7-30 amino acids.
2. **Splice-specific filtering**: Only peptides unique to the comparator
   (not present in tryptic digest of the reference protein) are retained.
3. **Proteotypic filtering**: Each peptide is checked against a pre-computed
   lookup of how many isoforms in the full proteome contain it. Proteotypic
   peptides (found in exactly 1 isoform) are flagged — these unambiguously
   identify the alternative isoform.
4. **PeptideAtlas matching**: Exact sequence match against PeptideAtlas.
   Quality metrics retained: `n_observations` (spectral count),
   `n_samples` (independent samples), `best_probability` (identification
   confidence).
5. **Summary**: Per-comparison metrics include total splice-specific
   peptides, PeptideAtlas hits, proteotypic hits, and best quality scores.

---

## References

**Splice Site Nomenclature**:
- Breathnach, R., & Chambon, P. (1981). Organization and expression of eucaryotic split genes coding for proteins. *Annual Review of Biochemistry*, 50(1), 349-383.
- Mount, S. M. (1982). A catalogue of splice junction sequences. *Nucleic Acids Research*, 10(2), 459-472.

**Event Classification**:
- Katz, Y., Wang, E. T., Airoldi, E. M., & Burge, C. B. (2010). Analysis and design of RNA sequencing experiments for identifying isoform regulation. *Nature Methods*, 7(12), 1009-1015.
- Shen, S., et al. (2014). rMATS: Robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. *PNAS*, 111(51), E5593-E5601.

**NMD Rules**:
- Nagy, E., & Maquat, L. E. (1998). A rule for termination-codon position within intron-containing genes. *Genes & Development*, 12(3), 665-676.

---

**Document Version**: 2.8
**Last Updated**: 2026-03-19
**Status**: Complete and validated 2026-03-20 — added methods for transcriptional diversity (Section 1), unified model (Section 2a), dose-response (Section 2d), uORF ATG position sharing (Section 3), reference-ATG tracing, PTC attribution implementation details
