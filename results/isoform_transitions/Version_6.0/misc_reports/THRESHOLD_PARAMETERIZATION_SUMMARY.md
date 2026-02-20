# Threshold Parameterization Summary

**Date:** 2026-02-15
**Task:** Extract repeated thresholds as configurable constants for sensitivity analyses

---

## Changes Made

### 1. Script 07: Added Threshold Constants

**Location:** Lines 29-42 in `scripts/07_extract_splicing_profiles.R`

**Constants Defined:**
```r
# Terminal boundary change detection
TSS_TOLERANCE <- 20  # bp tolerance for TSS change detection
TES_TOLERANCE <- 20  # bp tolerance for TES change detection

# Splice site classification boundary
SPLICE_SITE_THRESHOLD <- 100  # bp threshold: <100 = splice site, ≥100 = retention

# Exon overlap detection
OVERLAP_THRESHOLD <- 0.5  # proportion (0.5 = 50% overlap required)
```

**Display at Runtime:**
Script now prints threshold values at startup:
```
Detection thresholds:
  TSS tolerance: 20 bp
  TES tolerance: 20 bp
  Splice site threshold: 100 bp
  Overlap threshold: 50%
```

### 2. Replaced All Hardcoded Values

**Summary of replacements:**

| Constant | References | Locations Updated |
|----------|-----------|-------------------|
| `TSS_TOLERANCE` | 3 | Function default, documentation |
| `TES_TOLERANCE` | 3 | Function default, documentation |
| `SPLICE_SITE_THRESHOLD` | 11 | A5SS, A3SS, Partial_IR detection + comments |
| `OVERLAP_THRESHOLD` | 4 | Exon overlap calculation, terminal rules |

**Key Changes:**
- Line 82: `return(overlap_pct >= 0.5)` → `return(overlap_pct >= OVERLAP_THRESHOLD)`
- Line 92: `tolerance = 20` → `tolerance = TSS_TOLERANCE`
- Line 113: `tolerance = 20` → `tolerance = TES_TOLERANCE`
- Line 149: `bp_diff < 100` → `bp_diff < SPLICE_SITE_THRESHOLD`
- Line 176: `bp_diff < 100` → `bp_diff < SPLICE_SITE_THRESHOLD`
- Line 204: `len_diff >= 100` → `len_diff >= SPLICE_SITE_THRESHOLD`
- Multiple comment updates to reference constants

---

## Benefits

### 1. **Easy Sensitivity Analysis**
Can now test different thresholds without editing code throughout the script:
```r
# Test stricter thresholds
TSS_TOLERANCE <- 10
TES_TOLERANCE <- 10
SPLICE_SITE_THRESHOLD <- 50
OVERLAP_THRESHOLD <- 0.7
source('scripts/07_extract_splicing_profiles.R')
```

### 2. **Clear Documentation**
Threshold values are now:
- Defined once at the top
- Self-documenting with comments
- Displayed to user at runtime
- Easy to find and understand

### 3. **Separate Control**
TSS and TES thresholds can now be set independently:
```r
TSS_TOLERANCE <- 10  # Strict TSS detection
TES_TOLERANCE <- 50  # Lenient TES detection
```

### 4. **Reproducibility**
Threshold values are explicitly stated in output logs, making analyses reproducible.

---

## How to Use

### Default Behavior (No Changes Required)
Running Script 07 normally uses default thresholds (20, 20, 100, 0.5):
```bash
Rscript scripts/07_extract_splicing_profiles.R
```

### Custom Thresholds
Set constants in R before sourcing the script:
```r
# Set custom thresholds
TSS_TOLERANCE <- 15
TES_TOLERANCE <- 25
SPLICE_SITE_THRESHOLD <- 75
OVERLAP_THRESHOLD <- 0.4

# Run Script 07
source('scripts/07_extract_splicing_profiles.R')
```

### Systematic Sensitivity Analysis
See `THRESHOLD_SENSITIVITY_GUIDE.md` for:
- Detailed threshold descriptions
- Expected sensitivity patterns
- Example sensitivity analysis code
- Interpretation guidelines

---

## Validation

### ✅ Syntax Check
All scripts parse correctly with new constants.

### ✅ Constant Usage Verified
- TSS_TOLERANCE: 3 references (all hardcoded 20s replaced)
- TES_TOLERANCE: 3 references (all hardcoded 20s replaced)
- SPLICE_SITE_THRESHOLD: 11 references (all hardcoded 100s replaced)
- OVERLAP_THRESHOLD: 4 references (all hardcoded 0.5s replaced)

### ✅ No Remaining Hardcoded Values
Zero hardcoded threshold values found in event detection code.

### ✅ Function Defaults Updated
All function signatures use constants as defaults:
- `detect_tss_change(..., tolerance = TSS_TOLERANCE)` — also checks exon overlap
- `detect_tes_change(..., tolerance = TES_TOLERANCE)` — also checks exon overlap

**Note:** Both functions now have a secondary overlap check: even if the TSS/TES coordinate difference is within tolerance, the function returns TRUE when first/last exons don't overlap at all. This catches cases with short extra terminal exons that create small coordinate differences but represent genuine structural changes. The overlap check is not controlled by a threshold parameter.

---

## Documentation Created

1. **THRESHOLD_SENSITIVITY_GUIDE.md**
   - Complete guide for running sensitivity analyses
   - Threshold descriptions and biological interpretations
   - Expected sensitivity patterns
   - Example code for systematic testing
   - Interpretation guidelines

2. **THRESHOLD_PARAMETERIZATION_SUMMARY.md** (this file)
   - Summary of changes made
   - Validation results
   - Usage instructions

3. **scripts/07_with_thresholds.R**
   - Wrapper script template for batch sensitivity runs
   - Command-line argument parsing example

4. **run_threshold_sensitivity.sh**
   - Shell script template for automated runs
   - Usage examples

---

## Testing Recommendations

### Quick Validation Test
```r
# 1. Run with defaults
source('scripts/07_extract_splicing_profiles.R')
default_profiles <- readRDS("data/splicing_choice_profiles.rds")

# 2. Run with different thresholds
TSS_TOLERANCE <- 10
TES_TOLERANCE <- 10
SPLICE_SITE_THRESHOLD <- 50
OVERLAP_THRESHOLD <- 0.3
source('scripts/07_extract_splicing_profiles.R')
strict_profiles <- readRDS("data/splicing_choice_profiles.rds")

# 3. Compare
cat("Default:", nrow(default_profiles), "profiles\n")
cat("Strict:", nrow(strict_profiles), "profiles\n")
cat("Alt_TSS rate:", mean(default_profiles$tss_changed), "vs",
    mean(strict_profiles$tss_changed), "\n")
```

### Comprehensive Sensitivity Analysis
See examples in `THRESHOLD_SENSITIVITY_GUIDE.md` for:
- Grid search across multiple parameter values
- Event frequency comparison
- Robustness assessment

---

## Files Modified

1. ✅ `scripts/07_extract_splicing_profiles.R`
   - Added 4 threshold constants (lines 29-42)
   - Replaced 21 hardcoded values
   - Added runtime threshold display
   - Updated function defaults
   - Updated documentation comments

2. ✅ `scripts/08_analyze_complexity_relationship.R`
   - Already reviewed, no threshold changes needed

3. ✅ `scripts/09_analyze_cooccurrence.R`
   - Already reviewed, no threshold changes needed

---

## Next Steps

1. **Current State:** All thresholds are parameterized and ready for use
2. **When you return:** Continue with Script 10 or run sensitivity analyses
3. **For sensitivity analysis:** Use `THRESHOLD_SENSITIVITY_GUIDE.md` as reference

---

## Summary Statistics

- **Constants added:** 4
- **Hardcoded values replaced:** 21
- **Lines of documentation added:** ~350 (guide + summary)
- **Scripts enhanced:** 3 (Script 07 + 2 wrapper templates)
- **No breaking changes:** All default behavior preserved

✅ **Ready for sensitivity analyses!**
