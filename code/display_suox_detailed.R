#!/usr/bin/env Rscript
# Detailed visualization of SUOX to understand union exon model and event detection

library(tidyverse)

output_dir <- "results/isoform_transitions/v4.0_reference_based"

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   DETAILED SUOX EXAMPLE - Union Exon Model                   ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Load data
union_models <- readRDS(file.path(output_dir, "union_exon_models_full.rds"))
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))
test_results <- readRDS(file.path(output_dir, "test_results_10genes.rds"))

# Get SUOX union model
suox_model <- union_models[["SUOX"]]
suox_strand <- "+"

cat("Gene: SUOX (strand: +)\n")
cat(sprintf("Number of union exons: %d\n", suox_model$n_union_exons))
cat(sprintf("Number of isoforms: %d\n\n", suox_model$n_isoforms))

# Display union model structure
cat("═══════════════════════════════════════════════════════════════════\n")
cat("UNION EXON MODEL STRUCTURE\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

for (i in seq_along(suox_model$union_exons)) {
  union_exon <- suox_model$union_exons[[i]]

  cat(sprintf("─── Union Exon %d (type: %s) ───\n",
              union_exon$exon_number, union_exon$exon_type))

  variants <- union_exon$variants
  cat(sprintf("  Variants: %d\n", nrow(variants)))

  for (j in seq_len(nrow(variants))) {
    var <- variants[j, ]
    flags <- c()
    if (var$is_first) flags <- c(flags, "FIRST")
    if (var$is_last) flags <- c(flags, "LAST")
    flag_str <- if (length(flags) > 0) paste0(" [", paste(flags, collapse=", "), "]") else ""

    cat(sprintf("    %s: %d-%d (isoform_exon: %d)%s\n",
                var$isoform_id,
                var$start, var$end,
                var$exon_index,
                flag_str))
  }
  cat("\n")
}

# Now let's look at a specific transition
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SPECIFIC TRANSITION EXAMPLE\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

# Get the transition from test results
suox_transitions <- test_results %>%
  filter(gene_id == "SUOX", n_Alt_TSS > 0)

if (nrow(suox_transitions) > 0) {
  trans <- suox_transitions[1, ]

  iso_A <- trans$isoform_A
  iso_B <- trans$isoform_B

  cat(sprintf("Transition: %s → %s\n\n", iso_A, iso_B))

  # Get exon structures from exon_structures file
  exons_A_full <- exon_structures %>%
    filter(isoform_id == iso_A)

  exons_B_full <- exon_structures %>%
    filter(isoform_id == iso_B)

  # Display isoform A with union exon mapping
  cat(sprintf("─── Isoform A: %s ───\n", iso_A))
  if (nrow(exons_A_full) > 0) {
    starts_A <- exons_A_full$exon_starts[[1]]
    ends_A <- exons_A_full$exon_ends[[1]]

    cat("Isoform_Pos | Coordinates        | Union_Exon | Flags\n")
    cat("------------|--------------------|-----------|---------\n")

    for (i in seq_along(starts_A)) {
      # Find which union exon this maps to
      union_exon_num <- NA
      union_exon_type <- ""
      flags <- c()

      if (i == 1) flags <- c(flags, "FIRST")
      if (i == length(starts_A)) flags <- c(flags, "LAST")

      # Search union model for this exon
      for (ue_idx in seq_along(suox_model$union_exons)) {
        ue <- suox_model$union_exons[[ue_idx]]
        matching_var <- ue$variants %>%
          filter(isoform_id == iso_A,
                 start == starts_A[i],
                 end == ends_A[i])

        if (nrow(matching_var) > 0) {
          union_exon_num <- ue$exon_number
          union_exon_type <- ue$exon_type
          break
        }
      }

      flag_str <- if (length(flags) > 0) paste(flags, collapse=", ") else ""

      cat(sprintf("Exon %-6d | %9d - %9d | %-10s | %s\n",
                  i, starts_A[i], ends_A[i],
                  ifelse(is.na(union_exon_num), "???",
                        sprintf("#%d (%s)", union_exon_num, union_exon_type)),
                  flag_str))
    }
  }

  cat("\n")

  # Display isoform B with union exon mapping
  cat(sprintf("─── Isoform B: %s ───\n", iso_B))
  if (nrow(exons_B_full) > 0) {
    starts_B <- exons_B_full$exon_starts[[1]]
    ends_B <- exons_B_full$exon_ends[[1]]

    cat("Isoform_Pos | Coordinates        | Union_Exon | Flags\n")
    cat("------------|--------------------|-----------|---------\n")

    for (i in seq_along(starts_B)) {
      # Find which union exon this maps to
      union_exon_num <- NA
      union_exon_type <- ""
      flags <- c()

      if (i == 1) flags <- c(flags, "FIRST")
      if (i == length(starts_B)) flags <- c(flags, "LAST")

      # Search union model for this exon
      for (ue_idx in seq_along(suox_model$union_exons)) {
        ue <- suox_model$union_exons[[ue_idx]]
        matching_var <- ue$variants %>%
          filter(isoform_id == iso_B,
                 start == starts_B[i],
                 end == ends_B[i])

        if (nrow(matching_var) > 0) {
          union_exon_num <- ue$exon_number
          union_exon_type <- ue$exon_type
          break
        }
      }

      flag_str <- if (length(flags) > 0) paste(flags, collapse=", ") else ""

      cat(sprintf("Exon %-6d | %9d - %9d | %-10s | %s\n",
                  i, starts_B[i], ends_B[i],
                  ifelse(is.na(union_exon_num), "???",
                        sprintf("#%d (%s)", union_exon_num, union_exon_type)),
                  flag_str))
    }
  }

  cat("\n")

  # Display detected events
  cat("─── Detected Events ───\n")
  events <- trans$event_vector[[1]]

  cat("\nEvent_Type  | Union_Exon | Direction | Details\n")
  cat("------------|------------|-----------|---------------------------\n")

  for (i in seq_len(nrow(events))) {
    evt <- events[i, ]
    exon_str <- ifelse(is.na(evt$exon_number), "N/A", as.character(evt$exon_number))
    cat(sprintf("%-11s | %-10s | %-9s | %s\n",
                evt$event_type,
                exon_str,
                evt$direction,
                evt$detail))
  }

  cat("\n")

  # Explain the Alt_TSS
  cat("═══════════════════════════════════════════════════════════════════\n")
  cat("INTERPRETATION\n")
  cat("═══════════════════════════════════════════════════════════════════\n\n")

  cat("Alt_TSS Detection:\n")
  cat(sprintf("  Isoform A first exon (position 1): %d - %d\n", starts_A[1], ends_A[1]))
  cat(sprintf("  Isoform B first exon (position 1): %d - %d\n", starts_B[1], ends_B[1]))
  cat(sprintf("  TSS difference: %d bp (start coordinates on plus strand)\n",
              abs(starts_A[1] - starts_B[1])))
  cat(sprintf("  Threshold: >20bp → Alt_TSS detected ✓\n\n"))

  # Check if these first exons are in same or different union exons
  union_A_first <- NA
  union_B_first <- NA

  for (ue_idx in seq_along(suox_model$union_exons)) {
    ue <- suox_model$union_exons[[ue_idx]]

    if (any(ue$variants$isoform_id == iso_A &
            ue$variants$start == starts_A[1] &
            ue$variants$end == ends_A[1])) {
      union_A_first <- ue$exon_number
    }

    if (any(ue$variants$isoform_id == iso_B &
            ue$variants$start == starts_B[1] &
            ue$variants$end == ends_B[1])) {
      union_B_first <- ue$exon_number
    }
  }

  cat("Union Exon Grouping:\n")
  if (!is.na(union_A_first) && !is.na(union_B_first)) {
    if (union_A_first == union_B_first) {
      cat(sprintf("  Both first exons are in Union Exon #%d\n", union_A_first))
      cat("  → In v3.0, this WOULD be detected (grouped together)\n")
      cat("  → In v4.0, detected via EXPLICIT check (no longer dependent on grouping) ✓\n")
    } else {
      cat(sprintf("  Isoform A first exon: Union Exon #%d\n", union_A_first))
      cat(sprintf("  Isoform B first exon: Union Exon #%d\n", union_B_first))
      cat("  → Different union exons (too far apart to group)\n")
      cat("  → In v3.0, this would be MISSED (called as SE) ✗\n")
      cat("  → In v4.0, detected via EXPLICIT check ✓\n")
    }
  }
}

cat("\n✓ Detailed display complete!\n")
