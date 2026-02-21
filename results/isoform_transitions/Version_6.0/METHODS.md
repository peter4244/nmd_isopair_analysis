# Isoform Choice Analysis: Computational Methods

**Version:** 6.0
**Date:** 2026-02-20
**Purpose:** Detailed algorithmic documentation for splicing event detection and isoform comparison

---

## Table of Contents

1. [Overview](#overview)
2. [Union Exon Construction](#union-exon-construction)
3. [Event Detection Algorithms](#event-detection-algorithms)
4. [Terminal Boundary Rules](#terminal-boundary-rules)
5. [Validation Framework](#validation-framework)
6. [Statistical Methods](#statistical-methods)
7. [Implementation Details](#implementation-details)

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

**Phase 2: Event Detection and Reconstruction Validation**
- Compare each non-dominant (comparator) to dominant isoform
- Detect all splicing events using a hierarchical pipeline: IR -> boundary shifts -> exon skipping -> terminal events
- Reconstruct dominant isoform from comparator + events and verify coordinates match
- Generate splicing choice profiles
- Shared libraries: `event_detection_functions.R`, `reconstruction_functions.R`

**Phase 3: Statistical Analysis**
- Co-occurrence analysis
- Spatial organization
- Functional context
- Pattern classification

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
Region Types: "5'UTR", "CDS", "3'UTR", "contains_orf_start", "contains_orf_stop", "non_coding", "unknown"
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

**Scripts** (Version 6.0):
1. `01_prepare_dge_data_v2.R`: Load expression data, identify dominant isoforms
2. `02_extract_isoform_structures.R`: Parse GFF files
3. `03_build_union_exons.R`: Construct atomic union exon models
4. `04_extract_cds_annotations.R`: Extract CDS coordinates
5. `05_annotate_region_types.R`: Classify union exons by region
6. `06_filter_to_analysis_subset.R`: Apply expression and gene category filters
7. `07_extract_splicing_profiles.R`: Hierarchical event detection, generate profiles
   - `--reconstruction_check`: On-the-fly reconstruction verification (eliminates need for separate Script 08 run)
   - `--pairs-file <path>`: Specify explicit contrast pairs (TSV: gene_id, dominant_isoform_id, comparator_isoform_id)
   - `--test N`: Limit to first N genes
8. `08_validate_reconstruction.R`: Standalone reconstruction validation (batch verification of Script 07 output)
9-13: Downstream statistical analyses

**Shared Libraries**:
- `event_detection_functions.R`: Core event detection functions (19 functions), detection thresholds (TSS_TOLERANCE, TES_TOLERANCE, SPLICE_SITE_THRESHOLD)
- `reconstruction_functions.R`: Reconstruction from comparator + events (`reconstruct_dominant_v2`), verification (`verify_transcript`)

**Documentation**:
- `SPLICE_FUNCTION_CATALOG.md`: Comprehensive function catalog linking to splicing biology
- `EVENT_DETECTION_FLOW.md`: Mermaid flowchart of data flow through event detection

**Validation**:
- `testing/validate_synthetic_simple.R`: Automated validation on synthetic data (44/44)
- `testing/synthetic/TestData/`: Synthetic test genes and annotations
- `testing/VALIDATED_REAL_EXAMPLES.md`: Registry of manually validated real examples

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

**Document Version**: 2.1
**Last Updated**: 2026-02-20
**Status**: Complete and validated
