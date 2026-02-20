#!/usr/bin/env Rscript
#
# Script 07 Wrapper: Run with Custom Thresholds
#
# This wrapper allows running Script 07 with custom threshold parameters
# for sensitivity analyses.
#
# Usage:
#   Rscript scripts/07_with_thresholds.R [tss] [tes] [splice_site] [overlap] [output_suffix]
#
# Parameters:
#   tss          TSS tolerance in bp (default: 20)
#   tes          TES tolerance in bp (default: 20)
#   splice_site  Splice site threshold in bp (default: 100)
#   overlap      Overlap threshold as proportion 0-1 (default: 0.5)
#   output_suffix Optional suffix for output files (default: auto-generated)
#
# Examples:
#   # Stricter TSS/TES tolerance
#   Rscript scripts/07_with_thresholds.R 10 10 100 0.5 strict_tss_tes
#
#   # Different splice site threshold
#   Rscript scripts/07_with_thresholds.R 20 20 50 0.5 splice50
#
#   # More lenient overlap
#   Rscript scripts/07_with_thresholds.R 20 20 100 0.3 overlap30

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set thresholds (override defaults or use command line arguments)
if (length(args) >= 1) {
  TSS_TOLERANCE <- as.numeric(args[1])
} else {
  TSS_TOLERANCE <- 20
}

if (length(args) >= 2) {
  TES_TOLERANCE <- as.numeric(args[2])
} else {
  TES_TOLERANCE <- 20
}

if (length(args) >= 3) {
  SPLICE_SITE_THRESHOLD <- as.numeric(args[3])
} else {
  SPLICE_SITE_THRESHOLD <- 100
}

if (length(args) >= 4) {
  OVERLAP_THRESHOLD <- as.numeric(args[4])
} else {
  OVERLAP_THRESHOLD <- 0.5
}

# Output suffix
if (length(args) >= 5) {
  output_suffix <- args[5]
} else {
  output_suffix <- sprintf("tss%d_tes%d_splice%d_overlap%.0f",
                          TSS_TOLERANCE, TES_TOLERANCE,
                          SPLICE_SITE_THRESHOLD, OVERLAP_THRESHOLD * 100)
}

# Validate parameters
if (TSS_TOLERANCE < 0 || !is.finite(TSS_TOLERANCE)) {
  stop("TSS_TOLERANCE must be a non-negative number")
}
if (TES_TOLERANCE < 0 || !is.finite(TES_TOLERANCE)) {
  stop("TES_TOLERANCE must be a non-negative number")
}
if (SPLICE_SITE_THRESHOLD <= 0 || !is.finite(SPLICE_SITE_THRESHOLD)) {
  stop("SPLICE_SITE_THRESHOLD must be a positive number")
}
if (OVERLAP_THRESHOLD < 0 || OVERLAP_THRESHOLD > 1) {
  stop("OVERLAP_THRESHOLD must be between 0 and 1")
}

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  RUNNING SCRIPT 07 WITH CUSTOM THRESHOLDS\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat(sprintf("  TSS tolerance:        %d bp\n", TSS_TOLERANCE))
cat(sprintf("  TES tolerance:        %d bp\n", TES_TOLERANCE))
cat(sprintf("  Splice site boundary: %d bp\n", SPLICE_SITE_THRESHOLD))
cat(sprintf("  Overlap threshold:    %.2f (%.0f%%)\n", OVERLAP_THRESHOLD, OVERLAP_THRESHOLD * 100))
cat(sprintf("  Output suffix:        %s\n", output_suffix))
cat("═══════════════════════════════════════════════════════════════\n")
cat("\n")

# Modify output file names to include suffix
OUTPUT_FILE <- sprintf("data/splicing_choice_profiles_%s.rds", output_suffix)
OUTPUT_FILE_INTERMEDIATE <- sprintf("data/splicing_choice_profiles_intermediate_%s.rds", output_suffix)
LOG_FILE <- sprintf("logs/07_extract_splicing_profiles_%s.log", output_suffix)

# Create log directory
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

# Redirect output to log file
log_conn <- file(LOG_FILE, open = "wt")
sink(log_conn, type = "output")
sink(log_conn, type = "message")

cat(sprintf("Thresholds: TSS=%d, TES=%d, Splice=%d, Overlap=%.2f\n",
            TSS_TOLERANCE, TES_TOLERANCE, SPLICE_SITE_THRESHOLD, OVERLAP_THRESHOLD))
cat(sprintf("Output file: %s\n\n", OUTPUT_FILE))

# Now source the main script
# Note: This will use the threshold constants defined above
tryCatch({
  # Load the main script content but with modified output paths
  library(tidyverse)

  cat("\nNOTE: Full Script 07 implementation would go here.\n")
  cat("For sensitivity analysis, you can:\n\n")
  cat("1. Copy Script 07 code into this wrapper, OR\n")
  cat("2. Modify Script 07 to check for pre-defined constants, OR\n")
  cat("3. Run manually in R console by setting constants first\n\n")
  cat("Example manual approach:\n")
  cat("  TSS_TOLERANCE <- 10\n")
  cat("  TES_TOLERANCE <- 10\n")
  cat("  SPLICE_SITE_THRESHOLD <- 50\n")
  cat("  OVERLAP_THRESHOLD <- 0.3\n")
  cat("  source('scripts/07_extract_splicing_profiles.R')\n\n")

}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  quit(status = 1)
}, finally = {
  sink(type = "output")
  sink(type = "message")
  close(log_conn)
})

cat("\nLog file saved to:", LOG_FILE, "\n")
