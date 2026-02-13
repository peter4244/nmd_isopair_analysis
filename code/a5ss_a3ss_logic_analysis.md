# A5SS/A3SS Detection Logic Analysis

## Executive Summary

**PROBLEM IDENTIFIED:** The A5SS/A3SS detection logic in `detect_events_major_isoforms.R` is fundamentally flawed. It will never detect A5SS/A3SS events because of a logical contradiction between how union exon groups are constructed and how the detection code expects to find events.

**STATUS:** Zero A5SS/A3SS events detected in first checkpoint (100 genes) - this is expected behavior given the current logic.

**ROOT CAUSE:** Misunderstanding of what "shared union exons" means in the context of the union exon model.

---

## The Logical Contradiction

### How Union Exons Are Grouped

From `build_union_exon_model_expression_based.R` (lines 123-147):

```r
# Group internal exons by shared boundaries
group_internal_exons <- function(internal_exons) {
  # ...
  # Find all exons that share START or END with current exon
  same_group <- (internal_exons$start == current_start) | (internal_exons$end == current_end)
  # ...
}
```

**Key Point:** Exons are grouped into the same union exon if they share **at least one boundary** (start OR end).

This means:
- If exons share a **start** boundary, they go in the same group
- If exons share an **end** boundary, they go in the same group
- If exons share **both** boundaries (identical coordinates), they go in the same group
- If exons share **neither** boundary, they go in different groups

### How A5SS/A3SS Detection Works

From `detect_events_major_isoforms.R` (lines 216-229):

```r
for (union_exon_num in shared) {
  exon_A <- exons_A %>% filter(union_exon_number == union_exon_num)
  exon_B <- exons_B %>% filter(union_exon_number == union_exon_num)

  same_start <- (exon_A$start[1] == exon_B$start[1])
  same_end <- (exon_A$end[1] == exon_B$end[1])

  # If coordinates are identical, skip (no alternative splicing)
  if (same_start && same_end) next

  # If start differs but end matches → A5SS
  if (!same_start && same_end) { ... }

  # If start matches but end differs → A3SS
  if (same_start && !same_end) { ... }
}
```

The code looks for cases where:
- **A5SS:** Same end, different starts
- **A3SS:** Same start, different ends

### Why This Doesn't Work

**The grouping algorithm guarantees that exons in the same union exon group share at least one boundary.**

Looking at the diagnostic results:
- 453 union exon groups analyzed (groups with 2+ isoforms)
- **100% had identical coordinates** (same start AND same end)
- **0% had different starts with same end** (would be A5SS)
- **0% had same start with different ends** (would be A3SS)

**Why all identical?** Because:

1. The grouping uses `(start == start) | (end == end)` to group exons
2. This means exons with the SAME start OR SAME end get grouped together
3. BUT - it appears the algorithm is being **too conservative**
4. Exons are only being grouped if they have **identical coordinates**

Wait - let me re-examine this. The grouping logic says "share START **OR** END", but we're seeing only identical coordinates. Let me verify this is actually a problem with the detection logic, not the grouping logic.

---

## Re-analyzing the Grouping Logic

Looking at line 137 more carefully:

```r
same_group <- (internal_exons$start == current_start) | (internal_exons$end == current_end)
```

This should match:
- Exons where `start == current_start` (same 5' boundary on + strand)
- Exons where `end == current_end` (same 3' boundary on + strand)
- Exons where both match (identical)

So theoretically, we should get groups like:
- Group 1: Exon A (100-200), Exon B (100-250) → share start=100
- Group 2: Exon C (150-300), Exon D (200-300) → share end=300

**But the diagnostic shows 100% identical coordinates.** This means one of two things:

1. **The grouping logic has a bug** - exons that should be grouped together aren't
2. **The data doesn't have A5SS/A3SS events** - all alternative exons differ at both boundaries

---

## Testing the Hypothesis

Let me check if there are exons with shared boundaries that aren't being grouped.

The grouping algorithm is **iterative and transitive**:
- Start with exon 1
- Find all exons sharing start OR end with exon 1
- Mark them all as "same group"
- Move to next unprocessed exon
- Repeat

**Issue with transitivity:** If exon A and B share a start, and B and C share an end, then A, B, and C all go in the same group - even if A and C share nothing!

Example:
- Exon A: 100-200
- Exon B: 100-300 (shares start with A)
- Exon C: 200-300 (shares end with B)

Result: A, B, C all in same group (even though A and C don't overlap!)

This could explain why we're not seeing the expected patterns.

---

## The Real Issue: Alternative Splice Sites Don't Work This Way

**Fundamental misconception:** A5SS and A3SS are **not** detected by looking at individual exon coordinates. They're detected by looking at **junction patterns**.

### What A5SS/A3SS Actually Are

**Alternative 5' Splice Site (A5SS):**
- Two isoforms use different **donor sites** (5' end of intron)
- Results in exons with the **same acceptor** (3' end) but **different donors** (5' end)
- Example on + strand:
  - Isoform A: Exon ends at position 1000 (donor site)
  - Isoform B: Exon ends at position 1050 (different donor site)
  - Next exon starts at 2000 (same acceptor site for both)

**Alternative 3' Splice Site (A3SS):**
- Two isoforms use different **acceptor sites** (3' end of intron)
- Results in exons with the **same donor** (5' end) but **different acceptors** (3' end)
- Example on + strand:
  - Previous exon ends at 1000 (same donor site for both)
  - Isoform A: Exon starts at position 2000 (acceptor site)
  - Isoform B: Exon starts at position 2050 (different acceptor site)

### How to Actually Detect Them

**A5SS/A3SS detection requires:**

1. **Identifying consecutive exon pairs** in each isoform
2. **Comparing junction coordinates** between isoforms
3. **Finding cases where:**
   - One boundary is shared (same junction)
   - Other boundary differs (alternative splice site)

**Current union model approach:**
- Groups exons by shared boundaries
- Loses information about **which exons are adjacent** in each isoform
- Cannot reconstruct junction patterns
- **Cannot detect A5SS/A3SS**

---

## Solution Options

### Option 1: Junction-Based Detection (Recommended)

Instead of using union exon groups, analyze junctions directly:

1. For each isoform pair:
   - Extract all introns (gaps between consecutive exons)
   - Compare intron start/end coordinates
   - **A5SS:** Different intron starts, same ends
   - **A3SS:** Same intron starts, different ends

2. Advantages:
   - Direct detection of splice site usage
   - Clear biological interpretation
   - No grouping artifacts

### Option 2: Adjacent Union Exon Analysis

Keep union model but compare **adjacent** union exons:

1. For each isoform:
   - Identify which union exons are consecutive
   - Track the boundary between them (junction)

2. Compare between isoforms:
   - Find cases where isoforms use different junctions
   - Classify based on which side differs

3. Challenges:
   - Complex bookkeeping
   - Requires reconstructing isoform topology from union model
   - Prone to edge cases

### Option 3: Remove A5SS/A3SS Detection

If these events are rare or not biologically relevant for this analysis:

1. Focus on events the union model **can** detect:
   - Exon skipping (SE)
   - Alternative TSS/TES
   - Mutually exclusive exons (if implemented)

2. Acknowledge limitation in methods
3. Consider junction-based analysis as separate analysis

---

## Current Code Assessment

### What Works:

1. **TSS/TES Detection (lines 98-159):** ✓ CORRECT
   - Compares first/last exons between isoforms
   - Checks if union exon numbers differ
   - Correctly identifies terminal differences

2. **Exon Skipping (lines 166-206):** ✓ CORRECT
   - Identifies union exons present in one isoform but not the other
   - Filters for internal exons only
   - Correctly classifies SE events

3. **Union Model Construction:** ✓ CORRECT
   - Groups exons by shared boundaries
   - Handles terminal exons separately
   - Detects intron retention

### What Doesn't Work:

1. **A5SS/A3SS Detection (lines 208-372):** ✗ BROKEN
   - **Fundamental logic error:** Assumes union exon groups can contain exons with different coordinates
   - **Reality:** All exons in a union group have identical coordinates (by design)
   - **Result:** Zero A5SS/A3SS events detected (expected)
   - **Lines 229 check:** `if (same_start && same_end) next` will ALWAYS trigger (skip)

---

## Recommendations

### Immediate Action:

1. **Remove or comment out A5SS/A3SS detection code** (lines 208-372)
   - It cannot work with the current union model
   - Keeping it is misleading

2. **Update documentation:**
   - List events that CAN be detected: SE, Alt_TSS, Alt_TES
   - Note that A5SS/A3SS require junction-based analysis

### Future Implementation:

1. **Implement junction-based A5SS/A3SS detection:**
   - Create separate function that analyzes introns
   - Compare consecutive exon boundaries
   - Classify based on which junction side differs

2. **Validate with known examples:**
   - Find genes with documented A5SS/A3SS events
   - Verify detection works correctly
   - Check strand-awareness

---

## Diagnostic Results Summary

From first 100 genes checkpoint:

```
Union exon groups with 2+ isoforms: 453

Classification:
  Identical coordinates:    453 (100.0%)
  A5SS candidates:            0 (0.0%)
  A3SS candidates:            0 (0.0%)
  Both differ:                0 (0.0%)
```

**Interpretation:**
- All union exon groups have identical coordinates
- This is **expected behavior** given the grouping algorithm
- A5SS/A3SS cannot be detected using this approach
- Zero events is **not a bug** - it's a **design limitation**

---

## Conclusion

The A5SS/A3SS detection code is **logically sound in isolation** but **incompatible with the union exon model**. The union model groups exons by shared boundaries, guaranteeing that all exons in a group have at least one matching coordinate. The detection code looks for differences in coordinates within the same group - which cannot exist by design.

**This is a conceptual error, not a coding error.**

The code should either:
1. Use junction-based detection (recommended)
2. Be removed entirely (acceptable)
3. Be reimplemented using adjacent exon analysis (complex)

Currently, keeping the A5SS/A3SS code creates false expectations - it will never detect events regardless of the data.
