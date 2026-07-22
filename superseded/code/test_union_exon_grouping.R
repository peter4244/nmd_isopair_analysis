#!/usr/bin/env Rscript
# Test Union Exon Grouping Functions
# Verify coordinate-based grouping before full rebuild

library(tidyverse)

TSS_TES_TOLERANCE <- 20

# Source the grouping functions
source("code/build_union_exon_model_expression_based.R")

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("   TESTING UNION EXON GROUPING FUNCTIONS\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("\n")

# ============================================================================
# Test 1: group_first_exons() - same start
# ============================================================================
cat("TEST 1: First exons with same start\n")
cat("───────────────────────────────────────\n")

first_exons_test1 <- tibble(
  isoform_id = c("A", "B"),
  start = c(1000, 1000),
  end = c(1200, 1500),
  is_first = c(TRUE, TRUE),
  is_last = c(FALSE, FALSE),
  transcript_exon_number = c(1, 1)
)

result1 <- group_first_exons(first_exons_test1, tolerance = 20)
cat(sprintf("Input: 2 exons with same start (1000)\n"))
cat(sprintf("Expected: 1 group\n"))
cat(sprintf("Result: %d groups\n", length(result1)))
if (length(result1) == 1 && nrow(result1[[1]]) == 2) {
  cat("✓ PASS\n\n")
} else {
  cat("✗ FAIL\n\n")
  print(result1)
}

# ============================================================================
# Test 2: group_first_exons() - different starts (beyond tolerance)
# ============================================================================
cat("TEST 2: First exons with different starts (beyond tolerance)\n")
cat("─────────────────────────────────────────────────────────────\n")

first_exons_test2 <- tibble(
  isoform_id = c("A", "B"),
  start = c(1000, 1100),
  end = c(1200, 1500),
  is_first = c(TRUE, TRUE),
  is_last = c(FALSE, FALSE),
  transcript_exon_number = c(1, 1)
)

result2 <- group_first_exons(first_exons_test2, tolerance = 20)
cat(sprintf("Input: 2 exons, starts differ by 100bp (tolerance=20bp)\n"))
cat(sprintf("Expected: 2 groups\n"))
cat(sprintf("Result: %d groups\n", length(result2)))
if (length(result2) == 2) {
  cat("✓ PASS\n\n")
} else {
  cat("✗ FAIL\n\n")
  print(result2)
}

# ============================================================================
# Test 3: group_first_exons() - within tolerance
# ============================================================================
cat("TEST 3: First exons with starts within tolerance\n")
cat("────────────────────────────────────────────────\n")

first_exons_test3 <- tibble(
  isoform_id = c("A", "B"),
  start = c(1000, 1015),
  end = c(1200, 1500),
  is_first = c(TRUE, TRUE),
  is_last = c(FALSE, FALSE),
  transcript_exon_number = c(1, 1)
)

result3 <- group_first_exons(first_exons_test3, tolerance = 20)
cat(sprintf("Input: 2 exons, starts differ by 15bp (tolerance=20bp)\n"))
cat(sprintf("Expected: 1 group\n"))
cat(sprintf("Result: %d groups\n", length(result3)))
if (length(result3) == 1 && nrow(result3[[1]]) == 2) {
  cat("✓ PASS\n\n")
} else {
  cat("✗ FAIL\n\n")
  print(result3)
}

# ============================================================================
# Test 4: group_monoexonic_exons() - exact match
# ============================================================================
cat("TEST 4: Monoexonic exons with exact coordinates\n")
cat("───────────────────────────────────────────────\n")

monoexonic_test4 <- tibble(
  isoform_id = c("A", "B"),
  start = c(1000, 1000),
  end = c(1500, 1500),
  is_first = c(TRUE, TRUE),
  is_last = c(TRUE, TRUE),
  transcript_exon_number = c(1, 1)
)

result4 <- group_monoexonic_exons(monoexonic_test4)
cat(sprintf("Input: 2 monoexonic with same coords (1000-1500)\n"))
cat(sprintf("Expected: 1 group\n"))
cat(sprintf("Result: %d groups\n", length(result4)))
if (length(result4) == 1 && nrow(result4[[1]]) == 2) {
  cat("✓ PASS\n\n")
} else {
  cat("✗ FAIL\n\n")
  print(result4)
}

# ============================================================================
# Test 5: group_monoexonic_exons() - different coordinates
# ============================================================================
cat("TEST 5: Monoexonic exons with different coordinates\n")
cat("───────────────────────────────────────────────────\n")

monoexonic_test5 <- tibble(
  isoform_id = c("A", "B"),
  start = c(1000, 1000),
  end = c(1500, 1600),
  is_first = c(TRUE, TRUE),
  is_last = c(TRUE, TRUE),
  transcript_exon_number = c(1, 1)
)

result5 <- group_monoexonic_exons(monoexonic_test5)
cat(sprintf("Input: 2 monoexonic, different ends (1500 vs 1600)\n"))
cat(sprintf("Expected: 2 groups\n"))
cat(sprintf("Result: %d groups\n", length(result5)))
if (length(result5) == 2) {
  cat("✓ PASS\n\n")
} else {
  cat("✗ FAIL\n\n")
  print(result5)
}

# ============================================================================
# Test 6: group_last_exons() - same end
# ============================================================================
cat("TEST 6: Last exons with same end\n")
cat("─────────────────────────────────\n")

last_exons_test6 <- tibble(
  isoform_id = c("A", "B"),
  start = c(1200, 1300),
  end = c(1500, 1500),
  is_first = c(FALSE, FALSE),
  is_last = c(TRUE, TRUE),
  transcript_exon_number = c(5, 3)
)

result6 <- group_last_exons(last_exons_test6, tolerance = 20)
cat(sprintf("Input: 2 exons with same end (1500)\n"))
cat(sprintf("Expected: 1 group\n"))
cat(sprintf("Result: %d groups\n", length(result6)))
if (length(result6) == 1 && nrow(result6[[1]]) == 2) {
  cat("✓ PASS\n\n")
} else {
  cat("✗ FAIL\n\n")
  print(result6)
}

# ============================================================================
# Summary
# ============================================================================
cat("════════════════════════════════════════════════════════════════\n")
cat("   ALL TESTS COMPLETED\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("\n")
