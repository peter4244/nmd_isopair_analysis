#!/usr/bin/env Rscript
# Fix Event Vectors v3.0 - Remove Invalid Event Classifications
# Filters out SE and A3SS events at terminal positions

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   FIX EVENT VECTORS v3.0 - Remove Invalid Classifications     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Configuration
input_file <- "results/isoform_transitions/v3.0_reference_based/reference_event_vectors_v3.0_filtered.rds"
output_file <- "results/isoform_transitions/v3.0_reference_based/reference_event_vectors_v3.0_filtered_FIXED.rds"

# Load data
cat("Loading v3.0 event vectors...\n")
data <- readRDS(input_file)
cat("  Loaded:", nrow(data), "transitions\n\n")

# Calculate n_exons for each transition
cat("Calculating n_exons...\n")
data <- data %>%
  mutate(n_exons = map_int(event_vector, function(ev) {
    if (is.null(ev) || nrow(ev) == 0) return(NA_integer_)

    # If Alt_TES exists, use its exon_number (actual last exon)
    tes_events <- ev %>% filter(event_type == "Alt_TES")
    if (nrow(tes_events) > 0) {
      return(as.integer(tes_events$exon_number[1]))
    }

    # Otherwise, use max junction index + 1
    junction_events <- ev %>% filter(event_type %in% c("A5SS", "A3SS", "CONST", "SE"))
    if (nrow(junction_events) > 0) {
      max_junction <- max(junction_events$exon_number, na.rm = TRUE)
      if (is.infinite(max_junction)) return(NA_integer_)
      return(as.integer(max_junction + 1))
    }

    return(1L)
  })) %>%
  filter(!is.na(n_exons))

cat("  Calculated n_exons for:", nrow(data), "transitions\n\n")

# Count invalid events before filtering
cat("Counting invalid events BEFORE filtering...\n")
invalid_before <- data %>%
  mutate(
    invalid_se_e1 = map2_int(event_vector, n_exons, function(ev, n) {
      if (is.null(ev) || nrow(ev) == 0) return(0L)
      sum(ev$event_type == "SE" & ev$exon_number == 1)
    }),
    invalid_se_en = map2_int(event_vector, n_exons, function(ev, n) {
      if (is.null(ev) || nrow(ev) == 0) return(0L)
      n <- as.integer(n)
      sum(ev$event_type == "SE" & (ev$exon_number == n | ev$exon_number > n))
    }),
    invalid_a3ss_e1 = map_int(event_vector, function(ev) {
      if (is.null(ev) || nrow(ev) == 0) return(0L)
      sum(ev$event_type == "A3SS" & ev$exon_number == 1)
    })
  )

cat("  SE at E1 (junction 1):", sum(invalid_before$invalid_se_e1), "\n")
cat("  SE at En (last junction):", sum(invalid_before$invalid_se_en), "\n")
cat("  A3SS at E1:", sum(invalid_before$invalid_a3ss_e1), "\n")
cat("  Total invalid events:", sum(invalid_before$invalid_se_e1) +
                              sum(invalid_before$invalid_se_en) +
                              sum(invalid_before$invalid_a3ss_e1), "\n\n")

# Fix event vectors
cat("Filtering invalid events...\n")
data_fixed <- data %>%
  mutate(
    # Filter event vectors - ALWAYS remove invalid positions
    event_vector = map2(event_vector, n_exons, function(ev, n) {
      if (is.null(ev) || nrow(ev) == 0) return(ev)

      # Convert n to integer for proper comparison
      n <- as.integer(n)

      ev %>%
        # SE events: exon_number is ACTUAL EXON INDEX (not junction index)
        # Remove SE at exon 1 (first exon - terminal)
        filter(!(event_type == "SE" & exon_number == 1)) %>%
        # Remove SE at exon n (last exon - terminal)
        filter(!(event_type == "SE" & exon_number == n)) %>%
        # Remove SE at exons beyond n (invalid data)
        filter(!(event_type == "SE" & exon_number > n)) %>%
        # A3SS events: exon_number is JUNCTION INDEX
        # Remove A3SS at junction 1 (no upstream boundary at first exon)
        filter(!(event_type == "A3SS" & exon_number == 1))
    })
  )

cat("  Filtering complete\n\n")

# Recalculate event counts
cat("Recalculating event counts...\n")
data_fixed <- data_fixed %>%
  mutate(
    n_alt_tss = map_int(event_vector, ~sum(.x$event_type == "Alt_TSS")),
    n_alt_tes = map_int(event_vector, ~sum(.x$event_type == "Alt_TES")),
    n_se = map_int(event_vector, ~sum(.x$event_type == "SE")),
    n_a5ss = map_int(event_vector, ~sum(.x$event_type == "A5SS")),
    n_a3ss = map_int(event_vector, ~sum(.x$event_type == "A3SS")),
    n_ir = map_int(event_vector, ~sum(.x$event_type == "IR")),
    n_constitutive = map_int(event_vector, ~sum(.x$event_type == "CONST"))
  ) %>%
  select(-n_exons)  # Remove temporary n_exons column

# Validate fix
cat("Validating fix...\n")
validation <- data_fixed %>%
  mutate(n_exons = map_int(event_vector, function(ev) {
    if (is.null(ev) || nrow(ev) == 0) return(NA_integer_)
    tes_events <- ev %>% filter(event_type == "Alt_TES")
    if (nrow(tes_events) > 0) return(as.integer(tes_events$exon_number[1]))
    junction_events <- ev %>% filter(event_type %in% c("A5SS", "A3SS", "CONST", "SE"))
    if (nrow(junction_events) > 0) {
      max_junction <- max(junction_events$exon_number, na.rm = TRUE)
      if (is.infinite(max_junction)) return(NA_integer_)
      return(as.integer(max_junction + 1))
    }
    return(1L)
  })) %>%
  mutate(
    invalid_se_e1 = map2_int(event_vector, n_exons, function(ev, n) {
      if (is.null(ev) || nrow(ev) == 0) return(0L)
      sum(ev$event_type == "SE" & ev$exon_number == 1)
    }),
    invalid_se_en = map2_int(event_vector, n_exons, function(ev, n) {
      if (is.null(ev) || nrow(ev) == 0) return(0L)
      n <- as.integer(n)
      sum(ev$event_type == "SE" & (ev$exon_number == n | ev$exon_number > n))
    }),
    invalid_a3ss_e1 = map_int(event_vector, function(ev) {
      if (is.null(ev) || nrow(ev) == 0) return(0L)
      sum(ev$event_type == "A3SS" & ev$exon_number == 1)
    })
  )

cat("  SE at E1 (junction 1):", sum(validation$invalid_se_e1), "✓\n")
cat("  SE at En (last junction):", sum(validation$invalid_se_en), "✓\n")
cat("  A3SS at E1:", sum(validation$invalid_a3ss_e1), "✓\n")
cat("  Total invalid events:", sum(validation$invalid_se_e1) +
                              sum(validation$invalid_se_en) +
                              sum(validation$invalid_a3ss_e1), "✓\n\n")

if (sum(validation$invalid_se_e1) + sum(validation$invalid_se_en) + sum(validation$invalid_a3ss_e1) == 0) {
  cat("✓ SUCCESS: All invalid events removed!\n\n")
} else {
  cat("⚠ WARNING: Some invalid events remain!\n\n")
}

# Save fixed data
cat("Saving fixed event vectors...\n")
saveRDS(data_fixed, output_file)
cat("  Saved to:", output_file, "\n\n")

# Summary
cat("═══ Summary ═══\n\n")
cat("Events removed:\n")
cat("  SE at junction 1:", sum(invalid_before$invalid_se_e1), "\n")
cat("  SE at last junction:", sum(invalid_before$invalid_se_en), "\n")
cat("  A3SS at E1:", sum(invalid_before$invalid_a3ss_e1), "\n")
cat("  TOTAL:", sum(invalid_before$invalid_se_e1) +
                sum(invalid_before$invalid_se_en) +
                sum(invalid_before$invalid_a3ss_e1), "\n\n")

cat("Event counts AFTER fix:\n")
cat("  Alt_TSS:", sum(data_fixed$n_alt_tss), "\n")
cat("  Alt_TES:", sum(data_fixed$n_alt_tes), "\n")
cat("  SE:", sum(data_fixed$n_se), "\n")
cat("  A5SS:", sum(data_fixed$n_a5ss), "\n")
cat("  A3SS:", sum(data_fixed$n_a3ss), "\n")
cat("  IR:", sum(data_fixed$n_ir), "\n")
cat("  CONST:", sum(data_fixed$n_constitutive), "\n\n")

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   FIX COMPLETE - Event vectors cleaned and validated          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")
