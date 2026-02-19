# Reconstruction Validation Results - Issues Found

## Summary

**Verification Results:** 1/37 PASS (2.7% success rate)

The reconstruction pipeline successfully runs but produces incorrect results for 36/37 test cases. This indicates systematic issues in the reconstruction logic and/or event detection.

---

## Root Cause Analysis

### Issue #1: Missing Direction Information ⚠️ CRITICAL

**Problem:** Many events have `direction = "-"` instead of "GAIN" or "LOSS"

**Example:**
```
Event: A3SS  direction: "-"  bp_diff: 50
Original T1 exon: 25500-25700
Comparator T2 exon: 25550-25700
Reconstructed: 25550-25700 (WRONG - should be 25500-25700)
```

**Expected:**
- T2 is 50bp shorter (comparator lost sequence)
- Direction should be: **LOSS**
- Reconstruction: Extend acceptor by 50bp → 25550 - 50 = 25500 ✓

**Root Cause:**
The `detect_shared_boundary_event()` function only returns direction for overlap-based detection. For exact shared boundary matches (standard A5SS/A3SS), it returns `direction = NULL`, which gets recorded as "-".

**Impact:** Affects almost all boundary events (A5SS, A3SS, Partial_IR)

---

### Issue #2: Terminal Event Reconstruction Failures

**Problem:** Alt_TSS and Alt_TES events show exon count mismatches

**Examples:**
```
TEST_Alt_TSS_Basic: Expected 2 exons, reconstructed 3 exons
TEST_Alt_TES_NoOverlap_Plus: Expected 2 exons, reconstructed 3 exons
```

**Potential Causes:**
1. `missing_terminal_exons` field is not populated correctly
2. Terminal event coordinates (5_prime/3_prime) are ambiguous
3. Reconstruction logic for Alt_TSS/Alt_TES is incorrect

---

### Issue #3: IR Reconstruction Failures

**Problem:** IR events show exon count mismatches

**Examples:**
```
TEST_IR_Basic: Expected 2 exons, reconstructed 1 exon
TEST_IR_Minus: Expected 2 exons, reconstructed 1 exon
```

**Analysis:**
- IR with GAIN direction: Should split comparator's long exon → not happening
- IR with LOSS direction: Should merge comparator's exons → may be over-merging

**Likely cause:** Union exon lookup is not working correctly for IR splitting

---

### Issue #4: Coordinate Mismatches

**Pattern:** Many events show small coordinate differences (50-150bp)

**Examples:**
```
A5SS: Original [35000-35200], Reconstructed [35000-35250]  (50bp diff at donor)
Partial_IR_5: Original [40500-40650], Reconstructed [40500-40800]  (150bp diff)
```

**Root Cause:** Direction is missing or incorrectly interpreted, leading to:
- Extending when should shorten
- Shortening when should extend
- Using wrong sign for coordinate adjustment

---

## Detailed Failure Categories

### Category 1: Splice Site Events (A5SS, A3SS) - 18 failures
**Symptom:** Coordinate mismatches at donor or acceptor
**Cause:** Missing direction → cannot determine extend vs shorten

### Category 2: Partial IR Events - 8 failures
**Symptom:** Coordinate mismatches (usually larger, ≥100bp)
**Cause:** Same as splice site events + larger bp_diff

### Category 3: IR Events - 4 failures
**Symptom:** Wrong exon count
**Cause:** Splitting/merging logic not working

### Category 4: Terminal Events (Alt_TSS, Alt_TES) - 6 failures
**Symptom:** Wrong exon count or coordinate mismatch
**Cause:** missing_terminal_exons not populated or incorrectly used

---

## Required Fixes

### Fix #1: Add Direction Computation for All Events

**In `detect_and_save_events.R`:**

For A5SS/A3SS/Partial_IR events where direction is NULL:

```r
# Compute direction based on exon lengths
if (is.null(event_result$direction) || event_result$direction == "") {
  # Compare exon lengths to determine direction
  if (event_result$event_type %in% c("A5SS", "Partial_IR_5")) {
    # Donor differs - compare exon ends (plus) or starts (minus)
    if (strand == "+") {
      if (dom_exon$exon_end > comp_exon$exon_end) {
        direction <- "LOSS"  # Comparator is shorter
      } else {
        direction <- "GAIN"  # Comparator is longer
      }
    } else {
      if (dom_exon$exon_start < comp_exon$exon_start) {
        direction <- "LOSS"
      } else {
        direction <- "GAIN"
      }
    }
  } else if (event_result$event_type %in% c("A3SS", "Partial_IR_3")) {
    # Acceptor differs - compare exon starts (plus) or ends (minus)
    if (strand == "+") {
      if (dom_exon$exon_start < comp_exon$exon_start) {
        direction <- "LOSS"
      } else {
        direction <- "GAIN"
      }
    } else {
      if (dom_exon$exon_end > comp_exon$exon_end) {
        direction <- "LOSS"
      } else {
        direction <- "GAIN"
      }
    }
  }
}
```

### Fix #2: Improve Terminal Event Detection

**For Alt_TSS/Alt_TES events:**

1. Properly compute missing_terminal_exons:
   ```r
   # Walk through all 5' exons in dominant
   # Identify which ones are missing in comparator
   # Record as comma-separated ranges
   ```

2. Record correct 5_prime/3_prime coordinates:
   ```r
   # 5_prime = dominant TSS/TES
   # 3_prime = comparator TSS/TES
   # Or vice versa depending on direction
   ```

### Fix #3: Fix IR Reconstruction

**For IR splitting (GAIN direction):**

1. Ensure union exons are correctly passed to reconstruction
2. Filter union exons to the correct chromosome and region
3. Validate that union exons cover the IR region

### Fix #4: Improve Event Coordinate Recording

**Store exact boundary coordinates:**

For each event, record:
- Which exon (by position: first, last, internal)
- Which boundary (donor vs acceptor)
- Dominant boundary coordinate
- Comparator boundary coordinate

This removes ambiguity about which coordinate belongs to which isoform.

---

## Testing Strategy

### Phase 1: Fix Direction (Highest Priority)
1. Modify detect_and_save_events.R to compute direction for all events
2. Regenerate events file
3. Run reconstruction
4. Expected improvement: ~50% of cases should pass

### Phase 2: Fix Terminal Events
1. Improve missing_terminal_exons calculation
2. Fix Alt_TSS/Alt_TES reconstruction logic
3. Expected improvement: +20% pass rate

### Phase 3: Fix IR Events
1. Debug union exon usage in IR reconstruction
2. Test IR splitting and merging separately
3. Expected improvement: +10% pass rate

### Phase 4: Edge Cases
1. Address remaining coordinate mismatches
2. Handle complex multi-event cases
3. Target: 100% pass rate

---

## Next Steps

1. **Immediate:** Fix direction computation in event detection
2. **Short term:** Improve terminal event handling
3. **Medium term:** Debug IR reconstruction
4. **Final:** Achieve 100% verification pass rate

Once all fixes are implemented, we should see perfect reconstruction of all 37 synthetic test cases.

---

## Success Criteria

✅ **Goal:** 37/37 transcripts pass verification (100%)
✅ **Metric:** All exon coordinates match exactly
✅ **Evidence:** Event detection + reconstruction can reproduce original isoforms

This will validate that our event detection system is complete and accurate.
