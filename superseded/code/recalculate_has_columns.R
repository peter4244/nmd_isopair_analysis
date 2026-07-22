#!/usr/bin/env Rscript
# Recalculate has_* indicator columns after event vector fix

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   RECALCULATE has_* INDICATOR COLUMNS                         ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

output_dir <- "results/isoform_transitions/v3.0_reference_based"

# ============================================================================
# Fix FILTERED version
# ============================================================================

cat("═══ Processing FILTERED Version ═══\n\n")

filtered <- readRDS(file.path(output_dir, "reference_event_vectors_v3.0_filtered.rds"))
cat("Loaded:", nrow(filtered), "transitions\n\n")

# Recalculate has_* columns from actual event vectors
cat("Recalculating has_* indicator columns...\n")
filtered_fixed <- filtered %>%
  mutate(
    # Alt_TSS
    has_Alt_TSS_gain = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "Alt_TSS" & ev$direction == "gain")
    }),
    has_Alt_TSS_loss = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "Alt_TSS" & ev$direction == "loss")
    }),

    # Alt_TES
    has_Alt_TES_gain = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "Alt_TES" & ev$direction == "gain")
    }),
    has_Alt_TES_loss = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "Alt_TES" & ev$direction == "loss")
    }),

    # SE
    has_SE_gain = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "SE" & ev$direction == "gain")
    }),
    has_SE_loss = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "SE" & ev$direction == "loss")
    }),

    # A5SS
    has_A5SS_gain = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "A5SS" & ev$direction == "gain")
    }),
    has_A5SS_loss = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "A5SS" & ev$direction == "loss")
    }),

    # A3SS
    has_A3SS_gain = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "A3SS" & ev$direction == "gain")
    }),
    has_A3SS_loss = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "A3SS" & ev$direction == "loss")
    }),

    # IR
    has_IR_gain = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "IR" & ev$direction == "gain")
    }),
    has_IR_loss = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "IR" & ev$direction == "loss")
    })
  )

# Compare old vs new counts
comparison <- tibble(
  event_type = c("Alt_TSS_gain", "Alt_TSS_loss", "Alt_TES_gain", "Alt_TES_loss",
                 "SE_gain", "SE_loss", "A5SS_gain", "A5SS_loss",
                 "A3SS_gain", "A3SS_loss", "IR_gain", "IR_loss"),
  old_count = c(
    sum(filtered$has_Alt_TSS_gain), sum(filtered$has_Alt_TSS_loss),
    sum(filtered$has_Alt_TES_gain), sum(filtered$has_Alt_TES_loss),
    sum(filtered$has_SE_gain), sum(filtered$has_SE_loss),
    sum(filtered$has_A5SS_gain), sum(filtered$has_A5SS_loss),
    sum(filtered$has_A3SS_gain), sum(filtered$has_A3SS_loss),
    sum(filtered$has_IR_gain), sum(filtered$has_IR_loss)
  ),
  new_count = c(
    sum(filtered_fixed$has_Alt_TSS_gain), sum(filtered_fixed$has_Alt_TSS_loss),
    sum(filtered_fixed$has_Alt_TES_gain), sum(filtered_fixed$has_Alt_TES_loss),
    sum(filtered_fixed$has_SE_gain), sum(filtered_fixed$has_SE_loss),
    sum(filtered_fixed$has_A5SS_gain), sum(filtered_fixed$has_A5SS_loss),
    sum(filtered_fixed$has_A3SS_gain), sum(filtered_fixed$has_A3SS_loss),
    sum(filtered_fixed$has_IR_gain), sum(filtered_fixed$has_IR_loss)
  )
) %>%
  mutate(difference = new_count - old_count)

cat("\nComparison of has_* column counts:\n")
print(comparison %>% filter(difference != 0))

# Save
saveRDS(filtered_fixed, file.path(output_dir, "reference_event_vectors_v3.0_filtered.rds"))
cat("\n✓ Saved updated FILTERED file\n\n")

# ============================================================================
# Fix ALL version
# ============================================================================

cat("═══ Processing ALL Version ═══\n\n")

all_data <- readRDS(file.path(output_dir, "reference_event_vectors_v3.0_all.rds"))
cat("Loaded:", nrow(all_data), "transitions\n\n")

cat("Recalculating has_* indicator columns...\n")
all_fixed <- all_data %>%
  mutate(
    has_Alt_TSS_gain = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "Alt_TSS" & ev$direction == "gain")
    }),
    has_Alt_TSS_loss = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "Alt_TSS" & ev$direction == "loss")
    }),
    has_Alt_TES_gain = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "Alt_TES" & ev$direction == "gain")
    }),
    has_Alt_TES_loss = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "Alt_TES" & ev$direction == "loss")
    }),
    has_SE_gain = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "SE" & ev$direction == "gain")
    }),
    has_SE_loss = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "SE" & ev$direction == "loss")
    }),
    has_A5SS_gain = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "A5SS" & ev$direction == "gain")
    }),
    has_A5SS_loss = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "A5SS" & ev$direction == "loss")
    }),
    has_A3SS_gain = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "A3SS" & ev$direction == "gain")
    }),
    has_A3SS_loss = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "A3SS" & ev$direction == "loss")
    }),
    has_IR_gain = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "IR" & ev$direction == "gain")
    }),
    has_IR_loss = map_lgl(event_vector, function(ev) {
      any(ev$event_type == "IR" & ev$direction == "loss")
    })
  )

# Save
saveRDS(all_fixed, file.path(output_dir, "reference_event_vectors_v3.0_all.rds"))
cat("✓ Saved updated ALL file\n\n")

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   has_* COLUMNS RECALCULATED                                  ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")
