# Reconstruction System Implementation Progress

## Status: Phases 1-3 Complete ✅

---

## Phase 1: Event Detection ✅ COMPLETE

**Script:** `detect_and_save_events.R`

**Functionality:**
- Detects all splicing events between isoform pairs
- Saves complete event information with transcript IDs

**Output:** `synthetic/base_events_events.tsv`
- 69 events detected across 47 test cases
- All 11 required columns present
- Transcript IDs link events to specific isoform pairs

**Key Features:**
- Terminal boundary detection (Alt_TSS, Alt_TES)
- Internal event detection (A5SS, A3SS, Partial_IR, IR, SE)
- Direction tracking (GAIN/LOSS)
- Missing terminal exons recorded

---

## Phase 2: Comparator Extraction ✅ COMPLETE

**Script:** `extract_comparator_gtf.R`

**Functionality:**
- Filters original GTF to only comparator isoforms
- Uses events file to identify comparator transcript IDs

**Output:** `synthetic/base_events_comparator.gtf`
- 80 exons from 37 comparator transcripts
- Standard GTF format
- Ready for reconstruction

**Usage:**
```bash
Rscript extract_comparator_gtf.R \
  original.gtf \
  events.tsv \
  comparator_only.gtf
```

---

## Phase 3: Reconstruction Functions ✅ COMPLETE

**Script:** `reconstruction_functions.R`

**Functionality:**
Complete library of reconstruction functions for all event types.

### Implemented Functions

#### Boundary Event Functions
- **`apply_a5ss()`** - Alternative 5' splice site (donor)
  - Adjusts exon end (plus) or start (minus)
  - Uses direction (GAIN/LOSS) to determine extension vs shortening

- **`apply_a3ss()`** - Alternative 3' splice site (acceptor)
  - Adjusts exon start (plus) or end (minus)
  - Direction-aware boundary modification

- **`apply_partial_ir_5()`** - Partial intron retention at donor
  - Same logic as A5SS but with larger bp_diff (≥100bp)

- **`apply_partial_ir_3()`** - Partial intron retention at acceptor
  - Same logic as A3SS but with larger bp_diff (≥100bp)

#### Complex Event Functions
- **`apply_ir()`** - Intron retention
  - **GAIN direction**: Split comparator's long exon into multiple exons
    - Uses union exons to determine split points
  - **LOSS direction**: Merge comparator's exons into retained intron

- **`apply_se()`** - Skipped exon
  - Inserts missing exon between flanking exons
  - Uses event coordinates for exon boundaries

#### Terminal Event Functions
- **`apply_alt_tss()`** - Alternative transcription start site
  - **LOSS**: Add missing_terminal_exons to 5' end
  - **GAIN**: Trim 5' end to dominant TSS

- **`apply_alt_tes()`** - Alternative transcription end site
  - **LOSS**: Add missing_terminal_exons to 3' end
  - **GAIN**: Trim 3' end to dominant TES

#### Main Reconstruction
- **`reconstruct_dominant()`** - Apply all events in order
  - Prioritizes terminal events first
  - Applies internal events second
  - Skips Dual_boundary (decomposed into component events)
  - Error handling for robustness

### Helper Functions
- `find_exon_with_boundary()` - Find exon by boundary coordinate
- `find_exon_containing()` - Find exon overlapping coordinate
- `parse_coordinate_ranges()` - Parse missing_terminal_exons string

### Reconstruction Logic

**Event Priority:**
1. **Alt_TSS** (defines 5' boundary)
2. **Alt_TES** (defines 3' boundary)
3. **IR** (structural changes)
4. **SE** (exon insertion)
5. **A5SS, A3SS, Partial_IR** (boundary adjustments)

**Coordinate System:**
- All coordinates are 1-based (GTF convention)
- Strand-aware logic for all functions
- Maintains exon ordering after modification

---

## Generated Files

```
testing/
├── reconstruction/
│   ├── detect_and_save_events.R        ✅ Phase 1
│   ├── extract_comparator_gtf.R        ✅ Phase 2
│   └── reconstruction_functions.R      ✅ Phase 3
└── synthetic/
    ├── base_events_events.tsv          ✅ Events with transcript IDs
    └── base_events_comparator.gtf      ✅ Comparator-only GTF
```

---

## Next Steps: Phases 4-5

### Phase 4: Main Reconstruction Pipeline
**Script:** `reconstruct_dominant_isoforms.R`

**Tasks:**
- [ ] Load comparator GTF
- [ ] Load events file
- [ ] Load union exons (for IR reconstruction)
- [ ] For each gene/transcript pair:
  - [ ] Extract comparator exons
  - [ ] Filter events to this pair
  - [ ] Call `reconstruct_dominant()`
  - [ ] Write reconstructed exon to output GTF
- [ ] Generate summary statistics

**Expected Output:**
- Reconstructed dominant GTF
- Reconstruction log (events applied per transcript)

### Phase 5: Verification
**Script:** `verify_reconstruction.R`

**Tasks:**
- [ ] Load original dominant GTF (ground truth)
- [ ] Load reconstructed dominant GTF
- [ ] For each transcript:
  - [ ] Compare exon count
  - [ ] Compare all exon coordinates (exact match)
  - [ ] Record PASS/FAIL
- [ ] Generate verification report

**Success Criteria:**
- 100% exact coordinate matches for all synthetic test cases
- All 37 reconstructed transcripts identical to originals

---

## Design Decisions

### 1. Event Priority
Terminal events (Alt_TSS/Alt_TES) are applied first because they define the transcript boundaries. Internal events are then applied within those boundaries.

### 2. Union Exon Usage
IR reconstruction requires union exons to determine where to split merged exons. This ensures accurate reconstruction of multi-exon regions.

### 3. Error Handling
Each reconstruction function is wrapped in tryCatch to prevent single event failures from breaking the entire pipeline.

### 4. Missing Terminal Exons Format
Uses comma-separated coordinate ranges: `"1000-1200,1500-1700"`
- Handles non-contiguous exonic regions
- Easy to parse and reconstruct

### 5. Direction Semantics
- **LOSS**: Dominant is longer → extend comparator
- **GAIN**: Comparator is longer → shorten comparator

This consistent interpretation makes reconstruction logic straightforward.

---

## Testing Strategy

1. **Start simple**: Test A5SS and A3SS only (easiest cases)
2. **Add complexity**: Test Partial_IR, then IR
3. **Terminal events**: Test Alt_TSS and Alt_TES
4. **Complex cases**: Test real PacBio cases with multiple events
5. **Full validation**: Run all 37 synthetic test cases

**Validation approach:**
- Compare each reconstructed exon coordinate exactly
- No tolerance - must match perfectly
- Report first mismatch for debugging

---

## Implementation Notes

### Challenges Addressed

**Challenge 1: Finding the correct exon to modify**
- Solution: `find_exon_with_boundary()` and `find_exon_containing()`
- Tries exact boundary match first, then overlap match

**Challenge 2: IR reconstruction without union exons**
- Partial solution: Can merge exons (LOSS direction) without union exons
- Full solution: Need union exons for splitting (GAIN direction)

**Challenge 3: Terminal event coordinate interpretation**
- Solution: Use strand to determine which coordinate is TSS/TES
- Plus: lower = TSS, higher = TES
- Minus: higher = TSS, lower = TES

**Challenge 4: Exon insertion order for SE**
- Solution: Find insertion point based on genomic coordinates
- Insert after exon that ends before SE start

---

## Code Quality

**Features:**
- ✅ Comprehensive documentation (roxygen comments)
- ✅ Error handling (tryCatch wrappers)
- ✅ Helper functions for code reuse
- ✅ Strand-aware logic throughout
- ✅ Direction-aware reconstruction (GAIN/LOSS)

**Ready for:**
- Unit testing
- Integration with main pipeline
- Real data validation

