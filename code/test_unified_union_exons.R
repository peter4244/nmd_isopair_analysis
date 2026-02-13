#!/usr/bin/env Rscript
# Test Unified Union Exon Grouping
# Verifies position-based grouping with structural metadata

library(tidyverse)

# Source the build functions
source("code/build_union_exon_model_expression_based.R")

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   TESTING UNIFIED UNION EXON APPROACH                         ║\n")
cat("║   Position-Based Grouping with Structural Metadata            ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ============================================================================
# Test 1: A5SS - Exons sharing 3' end (different starts)
# ============================================================================
cat("TEST 1: A5SS - Exons with same end, different starts\n")
cat("──────────────────────────────────────────────────────────────\n")

test1_exons <- tibble(
  isoform_id = c("A", "B"),
  start = c(100, 150),
  end = c(200, 200),
  is_first = c(FALSE, FALSE),
  is_last = c(FALSE, FALSE),
  transcript_exon_number = c(2, 2)
)

result1 <- group_exons_by_boundaries(test1_exons)
result1_meta <- add_structural_metadata(result1[[1]])

cat("Input: 2 exons (100-200, 150-200) sharing end 200\n")
cat(sprintf("Expected: 1 union exon group containing both\n"))
cat(sprintf("Result: %d groups\n", length(result1)))

if (length(result1) == 1 && nrow(result1[[1]]) == 2) {
  cat("✓ PASS: Grouped correctly\n")
  cat(sprintf("  Metadata: first_in = %s, last_in = %s\n",
              ifelse(is.null(result1_meta$first_in_isoforms[[1]]), "NULL", "set"),
              ifelse(is.null(result1_meta$last_in_isoforms[[1]]), "NULL", "set")))
} else {
  cat("✗ FAIL\n")
}
cat("\n")

# ============================================================================
# Test 2: A3SS - Exons sharing 5' start (different ends)
# ============================================================================
cat("TEST 2: A3SS - Exons with same start, different ends\n")
cat("──────────────────────────────────────────────────────────────\n")

test2_exons <- tibble(
  isoform_id = c("A", "B"),
  start = c(100, 100),
  end = c(200, 250),
  is_first = c(FALSE, FALSE),
  is_last = c(FALSE, FALSE),
  transcript_exon_number = c(2, 2)
)

result2 <- group_exons_by_boundaries(test2_exons)

cat("Input: 2 exons (100-200, 100-250) sharing start 100\n")
cat(sprintf("Expected: 1 union exon group\n"))
cat(sprintf("Result: %d groups\n", length(result2)))

if (length(result2) == 1 && nrow(result2[[1]]) == 2) {
  cat("✓ PASS\n\n")
} else {
  cat("✗ FAIL\n\n")
}

# ============================================================================
# Test 3: Monoexonic + First-only sharing start
# ============================================================================
cat("TEST 3: Monoexonic + First-only sharing start boundary\n")
cat("──────────────────────────────────────────────────────────────\n")

test3_exons <- tibble(
  isoform_id = c("A", "B"),
  start = c(1000, 1000),
  end = c(1500, 1200),
  is_first = c(TRUE, TRUE),
  is_last = c(TRUE, FALSE),
  transcript_exon_number = c(1, 1)
)

result3 <- group_exons_by_boundaries(test3_exons)
result3_meta <- add_structural_metadata(result3[[1]])

cat("Input: Monoexonic (1000-1500) + First-only (1000-1200)\n")
cat(sprintf("Expected: 1 union exon group (shared start)\n"))
cat(sprintf("Result: %d groups\n", length(result3)))

if (length(result3) == 1 && nrow(result3[[1]]) == 2) {
  cat("✓ PASS: Grouped together\n")
  first_in <- result3_meta$first_in_isoforms[[1]]
  last_in <- result3_meta$last_in_isoforms[[1]]
  cat(sprintf("  First in: %s\n", paste(first_in, collapse=", ")))
  cat(sprintf("  Last in: %s\n", ifelse(is.null(last_in), "NULL", paste(last_in, collapse=", "))))

  if (length(first_in) == 2 && "A" %in% first_in && "B" %in% first_in) {
    cat("  ✓ Both marked as first\n")
  }
  if (!is.null(last_in) && "A" %in% last_in && length(last_in) == 1) {
    cat("  ✓ Only A marked as last\n")
  }
} else {
  cat("✗ FAIL\n")
}
cat("\n")

# ============================================================================
# Test 4: Cross-type coordinates (first vs internal)
# ============================================================================
cat("TEST 4: Same coordinates, different structural roles\n")
cat("──────────────────────────────────────────────────────────────\n")

test4_exons <- tibble(
  isoform_id = c("A", "B"),
  start = c(1000, 1000),
  end = c(1200, 1200),
  is_first = c(TRUE, FALSE),   # A: first exon, B: internal exon
  is_last = c(FALSE, FALSE),
  transcript_exon_number = c(1, 3)
)

result4 <- group_exons_by_boundaries(test4_exons)
result4_meta <- add_structural_metadata(result4[[1]])

cat("Input: Same coords (1000-1200), A=first, B=internal\n")
cat(sprintf("Expected: 1 union exon group (position-based)\n"))
cat(sprintf("Result: %d groups\n", length(result4)))

if (length(result4) == 1 && nrow(result4[[1]]) == 2) {
  cat("✓ PASS: Grouped together despite different roles\n")
  first_in <- result4_meta$first_in_isoforms[[1]]
  cat(sprintf("  First in: %s (only A should be marked)\n", paste(first_in, collapse=", ")))

  if (length(first_in) == 1 && first_in == "A") {
    cat("  ✓ Metadata correctly tracks structural role\n")
  }
} else {
  cat("✗ FAIL\n")
}
cat("\n")

# ============================================================================
# Test 5: No shared boundaries (should be separate groups)
# ============================================================================
cat("TEST 5: Exons with no shared boundaries\n")
cat("──────────────────────────────────────────────────────────────\n")

test5_exons <- tibble(
  isoform_id = c("A", "B"),
  start = c(100, 300),
  end = c(200, 400),
  is_first = c(FALSE, FALSE),
  is_last = c(FALSE, FALSE),
  transcript_exon_number = c(2, 3)
)

result5 <- group_exons_by_boundaries(test5_exons)

cat("Input: 2 exons (100-200, 300-400) no shared boundaries\n")
cat(sprintf("Expected: 2 separate union exon groups\n"))
cat(sprintf("Result: %d groups\n", length(result5)))

if (length(result5) == 2) {
  cat("✓ PASS\n\n")
} else {
  cat("✗ FAIL\n\n")
}

# ============================================================================
# Test 6: Complex case - multiple exons sharing boundaries
# ============================================================================
cat("TEST 6: Multiple exons with cascading shared boundaries\n")
cat("──────────────────────────────────────────────────────────────\n")

test6_exons <- tibble(
  isoform_id = c("A", "B", "C", "D"),
  start = c(100, 100, 100, 150),
  end = c(200, 250, 300, 200),
  is_first = c(FALSE, FALSE, FALSE, FALSE),
  is_last = c(FALSE, FALSE, FALSE, FALSE),
  transcript_exon_number = c(2, 2, 2, 2)
)

result6 <- group_exons_by_boundaries(test6_exons)

cat("Input: 4 exons - A(100-200), B(100-250), C(100-300), D(150-200)\n")
cat("  A,B,C share start 100\n")
cat("  A,D share end 200\n")
cat(sprintf("Expected: 1 union exon group (all connected via shared boundaries)\n"))
cat(sprintf("Result: %d groups\n", length(result6)))

if (length(result6) == 1 && nrow(result6[[1]]) == 4) {
  cat("✓ PASS: All grouped together via transitive boundary sharing\n\n")
} else {
  cat("✗ FAIL\n")
  print(result6)
  cat("\n")
}

# ============================================================================
# Test 7: Verify no cascading duplicates
# ============================================================================
cat("TEST 7: Prevent cascading duplicates (processed flag test)\n")
cat("──────────────────────────────────────────────────────────────\n")

test7_exons <- tibble(
  isoform_id = c("A", "B", "C"),
  start = c(1000, 1010, 1020),
  end = c(1100, 1100, 1100),
  is_first = c(FALSE, FALSE, FALSE),
  is_last = c(FALSE, FALSE, FALSE),
  transcript_exon_number = c(2, 2, 2)
)

result7 <- group_exons_by_boundaries(test7_exons)

cat("Input: 3 exons all sharing end 1100, different starts (1000, 1010, 1020)\n")
cat(sprintf("Expected: 1 union exon group (all share end)\n"))
cat(sprintf("Result: %d groups\n", length(result7)))

if (length(result7) == 1 && nrow(result7[[1]]) == 3) {
  cat("✓ PASS: No duplicates, all in one group\n")

  # Check each exon appears exactly once
  exon_counts <- result7[[1]] %>%
    group_by(isoform_id, start, end) %>%
    summarise(n = n(), .groups = 'drop')

  if (all(exon_counts$n == 1)) {
    cat("  ✓ Each exon appears exactly once\n")
  } else {
    cat("  ✗ FAIL: Found duplicate exons within group\n")
    print(exon_counts)
  }
} else {
  cat("✗ FAIL\n")
}
cat("\n")

# ============================================================================
# Summary
# ============================================================================
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   SYNTHETIC TESTS COMPLETED                                    ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")
