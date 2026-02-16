# Isoform Reconstruction and Verification System

## Overview

This system validates that our event detection is complete and accurate by:
1. Starting with **comparator isoforms**
2. Applying detected **events** to reconstruct **dominant isoforms**
3. Verifying reconstructed isoforms **exactly match** original dominant isoforms

---

## Workflow

```
Original GTF (dominant + comparator)
         ↓
   Event Detection
         ↓
   Events File (with transcript IDs)
         ↓
   Extract Comparator GTF
         ↓
   Reconstruction Script (apply events)
         ↓
   Reconstructed Dominant GTF
         ↓
   Verification (compare to original)
```

---

## 1. Event File Requirements

The event detection output must include transcript identifiers to link events to specific isoform pairs.

### Required Columns

```
gene_id                    Gene identifier
dominant_transcript_id     Dominant isoform transcript ID
comparator_transcript_id   Comparator isoform transcript ID
event_type                 A5SS, A3SS, Partial_IR_5, Partial_IR_3, IR, SE, Alt_TSS, Alt_TES
direction                  GAIN or LOSS
chr                        Chromosome
5_prime                    1-based position of biological 5' boundary
3_prime                    1-based position of biological 3' boundary
strand                     + or -
bp_diff                    Size difference in bp (for boundary events)
missing_terminal_exons     Comma-separated exonic coordinate ranges
```

**Critical:** Every event must be annotated with both `dominant_transcript_id` and `comparator_transcript_id`.

---

## 2. Comparator GTF Extraction

**Purpose:** Create a GTF file containing only comparator isoforms.

**Input:**
- Original GTF (all isoforms)
- Events file (to identify comparator transcript IDs)

**Output:**
- Comparator-only GTF

**Algorithm:**
```r
# Extract unique comparator IDs from events file
comparator_ids <- unique(events$comparator_transcript_id)

# Filter GTF to only these transcripts
comparator_gtf <- original_gtf %>%
  filter(transcript_id %in% comparator_ids)

# Write comparator GTF
write_gtf(comparator_gtf, "comparator_only.gtf")
```

---

## 3. Reconstruction Rules

For each event type, define the inverse operation to reconstruct the dominant isoform from the comparator.

### 3.1 Alternative Splice Sites (A5SS, A3SS)

**Event:** Boundary differs by `bp_diff` bases

**Reconstruction:**
- **Direction = LOSS**: Comparator is shorter → extend boundary by `bp_diff`
- **Direction = GAIN**: Comparator is longer → shorten boundary by `bp_diff`

**A5SS (donor differs):**
```
# Plus strand: donor = exon end
# Minus strand: donor = exon start

if (direction == "LOSS") {
  # Comparator lost sequence, extend to reconstruct dominant
  if (strand == "+") {
    exon_end = exon_end + bp_diff  # Move donor downstream
  } else {
    exon_start = exon_start - bp_diff  # Move donor upstream
  }
} else if (direction == "GAIN") {
  # Comparator gained sequence, shorten to reconstruct dominant
  if (strand == "+") {
    exon_end = exon_end - bp_diff
  } else {
    exon_start = exon_start + bp_diff
  }
}
```

**A3SS (acceptor differs):**
```
# Plus strand: acceptor = exon start
# Minus strand: acceptor = exon end

if (direction == "LOSS") {
  if (strand == "+") {
    exon_start = exon_start - bp_diff  # Move acceptor upstream
  } else {
    exon_end = exon_end + bp_diff  # Move acceptor downstream
  }
} else if (direction == "GAIN") {
  if (strand == "+") {
    exon_start = exon_start + bp_diff
  } else {
    exon_end = exon_end - bp_diff
  }
}
```

### 3.2 Partial Intron Retention (Partial_IR_5, Partial_IR_3)

**Event:** One boundary extends by `≥100bp` into intron

**Reconstruction:** Same as A5SS/A3SS but with larger `bp_diff`

### 3.3 Intron Retention (IR)

**Event:** Comparator has monoexonic region spanning multiple dominant exons

**Reconstruction:**
```
if (direction == "LOSS") {
  # Dominant has multiple exons, comparator has one
  # Split the comparator's long exon at the junction positions
  # Use 5_prime and 3_prime to identify junction boundaries

  # This is complex - need to query union exons within the IR region
  # to determine where to split the exon

} else if (direction == "GAIN") {
  # Comparator has intron retention, dominant does not
  # Merge adjacent dominant exons to create retained intron
  # (This case is rare - usually we compare dominant → comparator)
}
```

**Note:** IR reconstruction requires union exon structure to identify split positions.

### 3.4 Skipped Exon (SE)

**Event:** Exon present in one isoform, absent in the other

**Reconstruction:**
```
# SE events store the coordinates of the missing exon in 5_prime/3_prime
# Direction is always "-" (no clear gain/loss)

# Check which isoform has the exon by examining context
# If comparator lacks the exon, insert it between flanking exons
# Use 5_prime and 3_prime as the exon boundaries

# Add new exon:
new_exon = tibble(
  exon_start = min(5_prime, 3_prime),
  exon_end = max(5_prime, 3_prime)
)

# Insert between appropriate flanking exons
```

### 3.5 Alternative TSS (Alt_TSS)

**Event:** Transcription start site differs

**Reconstruction:**
```
if (direction == "LOSS") {
  # Dominant has longer 5' end
  # missing_terminal_exons contains regions to add back
  # Parse: "1000-1200,1500-1700,2000-2150"

  terminal_regions <- parse_coordinate_ranges(missing_terminal_exons)

  # Add these exonic regions to the 5' end of comparator
  for (region in terminal_regions) {
    add_exon(start = region$start, end = region$end)
  }

} else if (direction == "GAIN") {
  # Comparator has longer 5' end, dominant is shorter
  # Trim 5' end to match dominant TSS

  # Dominant TSS is stored in one of 5_prime or 3_prime
  # Use strand to determine which
  if (strand == "+") {
    dominant_tss = min(5_prime, 3_prime)  # Lower coordinate
    # Trim first exon start to dominant_tss
  } else {
    dominant_tss = max(5_prime, 3_prime)  # Higher coordinate
    # Trim first exon end to dominant_tss
  }
}
```

### 3.6 Alternative TES (Alt_TES)

**Event:** Transcription end site differs

**Reconstruction:** Similar to Alt_TSS but at 3' end

```
if (direction == "LOSS") {
  # Dominant has longer 3' end
  # Add back missing_terminal_exons at 3' end

} else if (direction == "GAIN") {
  # Comparator has longer 3' end
  # Trim 3' end to match dominant TES
}
```

---

## 4. Reconstruction Algorithm

**Inputs:**
1. Comparator GTF (exons for comparator isoform)
2. Events file (filtered to this gene, this comparator)
3. Union exons file (for IR reconstruction)

**Output:**
- Reconstructed dominant GTF

**Algorithm:**

```r
reconstruct_dominant <- function(comparator_exons, events, union_exons) {

  # Start with comparator structure
  reconstructed <- comparator_exons

  # Sort events by priority:
  # 1. Terminal events first (Alt_TSS, Alt_TES)
  # 2. Internal events (A5SS, A3SS, Partial_IR, IR, SE)
  events <- events %>%
    arrange(
      case_when(
        event_type %in% c("Alt_TSS", "Alt_TES") ~ 1,
        TRUE ~ 2
      )
    )

  # Apply each event
  for (i in seq_len(nrow(events))) {
    event <- events[i, ]

    if (event$event_type == "A5SS") {
      reconstructed <- apply_a5ss(reconstructed, event)

    } else if (event$event_type == "A3SS") {
      reconstructed <- apply_a3ss(reconstructed, event)

    } else if (event$event_type == "Partial_IR_5") {
      reconstructed <- apply_partial_ir_5(reconstructed, event)

    } else if (event$event_type == "Partial_IR_3") {
      reconstructed <- apply_partial_ir_3(reconstructed, event)

    } else if (event$event_type == "IR") {
      reconstructed <- apply_ir(reconstructed, event, union_exons)

    } else if (event$event_type == "SE") {
      reconstructed <- apply_se(reconstructed, event)

    } else if (event$event_type == "Alt_TSS") {
      reconstructed <- apply_alt_tss(reconstructed, event)

    } else if (event$event_type == "Alt_TES") {
      reconstructed <- apply_alt_tes(reconstructed, event)
    }
  }

  # Return reconstructed exon structure
  return(reconstructed)
}
```

---

## 5. Verification Algorithm

**Inputs:**
1. Reconstructed dominant GTF
2. Original dominant GTF

**Output:**
- Verification report (PASS/FAIL per transcript)

**Checks:**
1. **Exon count match**: Same number of exons
2. **Coordinate match**: All exon start/end coordinates identical
3. **Strand match**: Same strand
4. **Order match**: Exons in same biological order

**Algorithm:**

```r
verify_reconstruction <- function(original_exons, reconstructed_exons) {

  # Check exon count
  if (nrow(original_exons) != nrow(reconstructed_exons)) {
    return(list(
      pass = FALSE,
      reason = sprintf("Exon count mismatch: %d original vs %d reconstructed",
                      nrow(original_exons), nrow(reconstructed_exons))
    ))
  }

  # Order both by genomic position
  original_sorted <- original_exons %>% arrange(exon_start, exon_end)
  reconstructed_sorted <- reconstructed_exons %>% arrange(exon_start, exon_end)

  # Compare each exon
  for (i in seq_len(nrow(original_sorted))) {
    orig <- original_sorted[i, ]
    recon <- reconstructed_sorted[i, ]

    if (orig$exon_start != recon$exon_start || orig$exon_end != recon$exon_end) {
      return(list(
        pass = FALSE,
        reason = sprintf("Exon %d mismatch: original [%d-%d] vs reconstructed [%d-%d]",
                        i, orig$exon_start, orig$exon_end,
                        recon$exon_start, recon$exon_end)
      ))
    }
  }

  # All checks passed
  return(list(pass = TRUE, reason = "Perfect match"))
}
```

---

## 6. Implementation Plan

### Phase 1: Update Event Detection Output
**Script:** `validate_synthetic_simple.R` (or new `detect_and_save_events.R`)

**Tasks:**
- [ ] Add `dominant_transcript_id` and `comparator_transcript_id` to event records
- [ ] Write events to TSV file with all required columns
- [ ] Test on synthetic data

### Phase 2: Comparator Extraction
**Script:** `extract_comparator_gtf.R`

**Tasks:**
- [ ] Read events file to get comparator IDs
- [ ] Filter GTF to comparator transcripts only
- [ ] Write comparator GTF

### Phase 3: Reconstruction Functions
**Script:** `reconstruction_functions.R`

**Tasks:**
- [ ] Implement `apply_a5ss()`
- [ ] Implement `apply_a3ss()`
- [ ] Implement `apply_partial_ir_5()` and `apply_partial_ir_3()`
- [ ] Implement `apply_ir()` (requires union exon lookup)
- [ ] Implement `apply_se()`
- [ ] Implement `apply_alt_tss()`
- [ ] Implement `apply_alt_tes()`

### Phase 4: Main Reconstruction Script
**Script:** `reconstruct_dominant_isoforms.R`

**Tasks:**
- [ ] Load comparator GTF
- [ ] Load events file
- [ ] Load union exons
- [ ] For each gene/transcript pair:
  - [ ] Apply all events
  - [ ] Reconstruct dominant
- [ ] Write reconstructed GTF

### Phase 5: Verification Script
**Script:** `verify_reconstruction.R`

**Tasks:**
- [ ] Load original dominant GTF
- [ ] Load reconstructed dominant GTF
- [ ] Compare each transcript
- [ ] Generate verification report
- [ ] Report statistics (% perfect matches)

---

## 7. Success Criteria

**Perfect Reconstruction:** All synthetic test cases should achieve 100% exact coordinate match.

**Expected Results:**
- 38 synthetic test transcript pairs
- 38/38 (100%) perfect reconstructions
- All exon coordinates match exactly

**If reconstruction fails:**
- Indicates missing events or incorrect event classification
- Indicates bugs in event detection logic
- Requires fixes to event detection functions

---

## 8. Future Extensions

Once validation on synthetic data passes:

1. **Real data validation**: Test on sampled PacBio isoform pairs
2. **Round-trip testing**: Reconstruct dominant → detect events → reconstruct comparator
3. **Union exon mapping**: Add union_exon_ids to events for explicit reconstruction
4. **Multiple comparators**: Handle genes with >2 isoforms (star topology)

---

## 9. File Structure

```
Version_6.0/
├── testing/
│   ├── reconstruction/
│   │   ├── extract_comparator_gtf.R
│   │   ├── reconstruction_functions.R
│   │   ├── reconstruct_dominant_isoforms.R
│   │   ├── verify_reconstruction.R
│   │   └── run_reconstruction_test.R  # Main test runner
│   └── synthetic/
│       ├── base_events.gtf            # Original (both isoforms)
│       ├── base_events_comparator.gtf # Extracted comparators
│       ├── base_events_events.tsv     # Detected events with IDs
│       └── base_events_reconstructed.gtf  # Reconstructed dominants
```

---

## Next Steps

1. **Update validation script** to save events with transcript IDs
2. **Implement comparator extraction** script
3. **Build reconstruction functions** one event type at a time
4. **Test incrementally** on simple cases first (A5SS, A3SS only)
5. **Add complex events** (IR, SE, terminal events)
6. **Run full verification** on all 38 synthetic tests

