# Threshold Sensitivity Analysis Guide

**Purpose:** Run Script 07 with different threshold parameters to understand how event detection sensitivity changes with parameter choices.

---

## Configurable Thresholds

Script 07 now uses clearly defined constants (at the top of the script) that can be modified for sensitivity analyses:

### 1. **TSS_TOLERANCE** (default: 20 bp)
- **What it controls:** Minimum coordinate difference to call an alternative TSS event
- **When Alt_TSS is detected:** When first exons differ at TSS by MORE than this threshold, OR when first exons don't overlap at all (regardless of coordinate distance)
- **Note:** The overlap check is not affected by this threshold — non-overlapping first exons always trigger Alt_TSS detection
- **Interpretation:**
  - **Smaller value (e.g., 10 bp)**: More sensitive - detects even small TSS shifts
  - **Larger value (e.g., 50 bp)**: More stringent - only calls substantial TSS changes (but non-overlapping exons still detected)

### 2. **TES_TOLERANCE** (default: 20 bp)
- **What it controls:** Minimum coordinate difference to call an alternative TES event
- **When Alt_TES is detected:** When last exons differ at TES by MORE than this threshold, OR when last exons don't overlap at all (regardless of coordinate distance)
- **Note:** Same overlap check as TSS — non-overlapping last exons always trigger Alt_TES detection
- **Interpretation:**
  - **Smaller value (e.g., 10 bp)**: Detects subtle polyadenylation site changes
  - **Larger value (e.g., 50 bp)**: Only calls major TES shifts (but non-overlapping exons still detected)

### 3. **SPLICE_SITE_THRESHOLD** (default: 100 bp)
- **What it controls:** Boundary between splice site variation vs. intron retention
- **Classification rules:**
  - **Difference < threshold**: Classified as A5SS or A3SS (splice site choice)
  - **Difference ≥ threshold**: Classified as Partial_IR (intron retention)
- **Interpretation:**
  - **Smaller value (e.g., 50 bp)**: More events classified as Partial_IR
  - **Larger value (e.g., 200 bp)**: More events classified as A5SS/A3SS

### 4. **OVERLAP_THRESHOLD** (default: 0.5 = 50%)
- **What it controls:** Minimum overlap to consider two exons as "the same exonic region"
- **Where it's used:**
  - Determining if first/last exons overlap for terminal splice site detection
  - Skipped exon (SE) detection - identifying comparable flanking exons
- **Interpretation:**
  - **Smaller value (e.g., 0.3 = 30%)**: More lenient - allows larger coordinate differences
  - **Larger value (e.g., 0.7 = 70%)**: More stringent - requires closer coordinate match

---

## How to Run Sensitivity Analyses

### Method 1: Manual in R Console (Recommended)

The simplest approach for one-off sensitivity tests:

```r
# Set thresholds BEFORE sourcing Script 07
TSS_TOLERANCE <- 10          # Stricter TSS detection
TES_TOLERANCE <- 10          # Stricter TES detection
SPLICE_SITE_THRESHOLD <- 50  # Lower boundary for retention
OVERLAP_THRESHOLD <- 0.3     # More lenient overlap

# Run Script 07 with these thresholds
source('scripts/07_extract_splicing_profiles.R')

# Results will be saved to:
#   data/splicing_choice_profiles.rds
# (Rename this file to preserve it before running again)
```

### Method 2: Modify Script 07 Directly

For temporary analyses, edit the constant definitions at the top of Script 07:

```r
# Lines 29-42 in Script 07
TSS_TOLERANCE <- 10   # Change from default 20
TES_TOLERANCE <- 10   # Change from default 20
# ... etc
```

Then run normally:
```bash
Rscript scripts/07_extract_splicing_profiles.R
```

### Method 3: Systematic Sensitivity Analysis

For testing multiple threshold combinations:

```r
# sensitivity_analysis.R
library(tidyverse)

# Define parameter grid
param_grid <- expand.grid(
  tss_tol = c(10, 20, 30),
  tes_tol = c(10, 20, 30),
  splice_tol = c(50, 100, 150),
  overlap_tol = c(0.3, 0.5, 0.7)
)

results_summary <- list()

for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]

  cat(sprintf("\n=== Run %d/%d ===\n", i, nrow(param_grid)))
  cat(sprintf("TSS=%d, TES=%d, Splice=%d, Overlap=%.1f\n",
              params$tss_tol, params$tes_tol,
              params$splice_tol, params$overlap_tol))

  # Set thresholds
  TSS_TOLERANCE <- params$tss_tol
  TES_TOLERANCE <- params$tes_tol
  SPLICE_SITE_THRESHOLD <- params$splice_tol
  OVERLAP_THRESHOLD <- params$overlap_tol

  # Source Script 07
  source('scripts/07_extract_splicing_profiles.R')

  # Load results
  profiles <- readRDS("data/splicing_choice_profiles.rds")

  # Summarize event frequencies
  results_summary[[i]] <- tibble(
    run_id = i,
    tss_tol = params$tss_tol,
    tes_tol = params$tes_tol,
    splice_tol = params$splice_tol,
    overlap_tol = params$overlap_tol,
    n_profiles = nrow(profiles),
    pct_alt_tss = mean(profiles$tss_changed),
    pct_alt_tes = mean(profiles$tes_changed),
    pct_a5ss = mean(profiles$n_a5ss > 0),
    pct_a3ss = mean(profiles$n_a3ss > 0),
    pct_partial_ir = mean(profiles$n_partial_ir > 0),
    pct_ir = mean(profiles$n_ir > 0),
    pct_se = mean(profiles$n_se > 0)
  )

  # Save this run's results
  saveRDS(profiles, sprintf("data/profiles_run%03d.rds", i))
}

# Compile results
sensitivity_results <- bind_rows(results_summary)
write_csv(sensitivity_results, "results/threshold_sensitivity_summary.csv")

# Analyze sensitivity
cat("\n=== SENSITIVITY ANALYSIS RESULTS ===\n")
print(sensitivity_results)
```

---

## Expected Sensitivity Patterns

### TSS/TES Tolerance Effects

| Threshold | Effect on Alt_TSS/TES Detection |
|-----------|----------------------------------|
| 5-10 bp   | Very sensitive - captures minor isoform TSS/TES shifts |
| 20 bp ✓   | **Default** - balanced between noise and signal |
| 50+ bp    | Conservative - only major promoter/polyA site changes |

**Recommendation:** Start with 20bp. Use 10bp if interested in subtle regulatory changes, 50bp if focused on major structural differences.

### Splice Site Threshold Effects

| Threshold | A5SS/A3SS Events | Partial_IR Events |
|-----------|------------------|-------------------|
| 50 bp     | Fewer (more stringent) | More (includes moderate size differences) |
| 100 bp ✓  | **Default** - biological convention | Standard definition |
| 200 bp    | More (includes larger changes) | Fewer (only large retention) |

**Recommendation:** 100bp is standard in splice site literature. Only change if you have specific biological reasons (e.g., studying microexons: use 50bp).

### Overlap Threshold Effects

| Threshold | SE Detection | Terminal Splice Site Detection |
|-----------|--------------|-------------------------------|
| 0.3 (30%) | More SE events (lenient flanking) | More terminal A5SS/A3SS calls |
| 0.5 (50%) ✓ | **Default** - balanced | Standard overlap definition |
| 0.7 (70%) | Fewer SE events (stringent) | Fewer terminal splice site calls |

**Recommendation:** 50% is standard for overlap detection. Use 30% if you suspect splice site variation at exon boundaries.

---

## Interpreting Sensitivity Results

### High Sensitivity (Expected)
If changing a threshold by ±50% causes ±50% change in event counts:
- **TSS/TES tolerance**: Expected - many isoforms have TSS/TES shifts near the threshold
- **Interpretation:** Threshold is in a "decision boundary" region - many events are borderline

### Low Sensitivity (Robust)
If changing threshold by ±50% causes <10% change in event counts:
- **Good news:** Detection is robust to parameter choice
- **Interpretation:** Events have clear, unambiguous characteristics

### Recommended Checks

1. **Rerun with stricter thresholds** (TSS=10, TES=10, Splice=50, Overlap=0.7)
   - If results are very similar: detection is robust ✓
   - If results change substantially: thresholds matter for your dataset

2. **Rerun with more lenient thresholds** (TSS=50, TES=50, Splice=200, Overlap=0.3)
   - Compare event proportions to defaults
   - Document which thresholds cause largest changes

3. **Report in methods:**
   ```
   "Event detection used thresholds of 20bp for TSS/TES changes,
   100bp for splice site classification, and 50% for exon overlap.
   Sensitivity analysis with ±50% threshold variation showed
   <15% change in event counts, indicating robust detection."
   ```

---

## Files Modified

1. **scripts/07_extract_splicing_profiles.R**
   - Added threshold constants at top (lines 29-42)
   - Replaced all hardcoded values with constants
   - Added threshold display at script start

2. **THRESHOLD_SENSITIVITY_GUIDE.md** (this file)
   - Complete documentation for sensitivity analyses

---

## Quick Start Example

```r
# 1. Test stricter detection
TSS_TOLERANCE <- 10
TES_TOLERANCE <- 10
SPLICE_SITE_THRESHOLD <- 50
OVERLAP_THRESHOLD <- 0.7
source('scripts/07_extract_splicing_profiles.R')
profiles_strict <- readRDS("data/splicing_choice_profiles.rds")

# 2. Test lenient detection
TSS_TOLERANCE <- 50
TES_TOLERANCE <- 50
SPLICE_SITE_THRESHOLD <- 200
OVERLAP_THRESHOLD <- 0.3
source('scripts/07_extract_splicing_profiles.R')
profiles_lenient <- readRDS("data/splicing_choice_profiles.rds")

# 3. Compare
cat("Strict:  ", nrow(profiles_strict), "profiles\n")
cat("Lenient: ", nrow(profiles_lenient), "profiles\n")

cat("\nEvent frequencies:\n")
cat("Alt_TSS: ", mean(profiles_strict$tss_changed), "vs",
    mean(profiles_lenient$tss_changed), "\n")
# ... etc for all event types
```

---

## Contact

For questions about threshold choices or sensitivity analysis results, consult with Pete Castaldi or refer to the METHODS.md document for detailed algorithm descriptions.
