#!/usr/bin/env Rscript
# Investigate why ENST00000886413.1's first exon is missing from SUOX union model

library(tidyverse)

output_dir <- "results/isoform_transitions/v4.0_reference_based"

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   INVESTIGATING MISSING FIRST EXON                           ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Load data
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))
union_models <- readRDS(file.path(output_dir, "union_exon_models_full.rds"))
isoforms_for_union <- readRDS(file.path(output_dir, "isoforms_for_union_model.rds"))

target_isoform <- "ENST00000886413.1"
gene <- "SUOX"

cat("═══ Step 1: Check if isoform exists in input data ═══\n\n")

# Check exon_structures
exon_data <- exon_structures %>%
  filter(isoform_id == target_isoform)

if (nrow(exon_data) > 0) {
  cat(sprintf("✓ Found %s in exon_structures\n", target_isoform))
  cat(sprintf("  Gene: %s\n", exon_data$gene_id))
  cat(sprintf("  Strand: %s\n", exon_data$strand))
  cat(sprintf("  Number of exons: %d\n\n", exon_data$n_exons))

  cat("  Exon coordinates:\n")
  starts <- exon_data$exon_starts[[1]]
  ends <- exon_data$exon_ends[[1]]
  for (i in seq_along(starts)) {
    cat(sprintf("    Exon %d: %d - %d\n", i, starts[i], ends[i]))
  }
  cat("\n")
} else {
  cat(sprintf("✗ %s NOT found in exon_structures\n\n", target_isoform))
}

# Check isoforms_for_union
in_union_list <- isoforms_for_union %>%
  filter(isoform_id == target_isoform)

if (nrow(in_union_list) > 0) {
  cat(sprintf("✓ Found %s in isoforms_for_union_model\n", target_isoform))
  cat(sprintf("  Gene: %s\n\n", in_union_list$gene_id))
} else {
  cat(sprintf("✗ %s NOT found in isoforms_for_union_model\n\n", target_isoform))
}

cat("═══ Step 2: Check SUOX union model construction ═══\n\n")

# Get all SUOX isoforms
suox_isoforms <- isoforms_for_union %>%
  filter(gene_id == gene)

cat(sprintf("SUOX has %d isoforms in input list:\n", nrow(suox_isoforms)))
for (i in seq_len(nrow(suox_isoforms))) {
  cat(sprintf("  %d. %s\n", i, suox_isoforms$isoform_id[i]))
}
cat("\n")

# Check if all these isoforms have exon structures
cat("Checking exon structures for all SUOX isoforms:\n")
for (iso in suox_isoforms$isoform_id) {
  exon_info <- exon_structures %>%
    filter(isoform_id == iso)

  if (nrow(exon_info) > 0) {
    starts <- exon_info$exon_starts[[1]]
    ends <- exon_info$exon_ends[[1]]
    cat(sprintf("  %s: %d exons\n", iso, length(starts)))
    cat(sprintf("    First: %d - %d\n", starts[1], ends[1]))
    cat(sprintf("    Last:  %d - %d\n", starts[length(starts)], ends[length(ends)]))
  } else {
    cat(sprintf("  %s: NO EXON DATA\n", iso))
  }
}
cat("\n")

cat("═══ Step 3: Reconstruct union model for SUOX ═══\n\n")

# Let's manually trace through union model construction for first exons
cat("Extracting all first exons for SUOX isoforms:\n\n")

all_first_exons <- map_dfr(suox_isoforms$isoform_id, function(iso) {
  exon_info <- exon_structures %>%
    filter(isoform_id == iso)

  if (nrow(exon_info) == 0) return(NULL)

  starts <- exon_info$exon_starts[[1]]
  ends <- exon_info$exon_ends[[1]]

  tibble(
    isoform_id = iso,
    exon_index = 1,
    start = starts[1],
    end = ends[1],
    is_first = TRUE,
    is_last = (length(starts) == 1)
  )
})

cat("First exons:\n")
print(all_first_exons %>% arrange(start), n = 100)
cat("\n")

# Check grouping logic
cat("═══ Step 4: Understanding first exon grouping ═══\n\n")

cat("First exon grouping uses TSS_TES_TOLERANCE = 20bp\n")
cat("Grouping criteria: first exons within 20bp at TSS (start on plus)\n\n")

# Sort by start and check which ones would be grouped
sorted_first <- all_first_exons %>%
  arrange(start) %>%
  mutate(
    start_diff_from_prev = start - lag(start, default = start[1]),
    would_group_with_prev = start_diff_from_prev <= 20
  )

cat("First exons sorted by start coordinate:\n")
print(sorted_first %>% select(isoform_id, start, end, start_diff_from_prev, would_group_with_prev), n = 100)
cat("\n")

# Identify groups
cat("═══ Step 5: Expected first exon groups ═══\n\n")

# Manual grouping based on 20bp tolerance
groups <- list()
current_group <- 1
groups[[1]] <- c(sorted_first$isoform_id[1])

for (i in 2:nrow(sorted_first)) {
  if (sorted_first$would_group_with_prev[i]) {
    # Add to current group
    groups[[current_group]] <- c(groups[[current_group]], sorted_first$isoform_id[i])
  } else {
    # Start new group
    current_group <- current_group + 1
    groups[[current_group]] <- c(sorted_first$isoform_id[i])
  }
}

cat(sprintf("Expected number of first exon groups: %d\n\n", length(groups)))

for (i in seq_along(groups)) {
  cat(sprintf("Group %d: %d isoforms\n", i, length(groups[[i]])))
  for (iso in groups[[i]]) {
    exon_info <- filter(all_first_exons, isoform_id == iso)
    cat(sprintf("  %s: %d - %d\n", iso, exon_info$start, exon_info$end))
  }
  cat("\n")
}

# Compare to actual union model
cat("═══ Step 6: Compare to actual union model ═══\n\n")

suox_union <- union_models[["SUOX"]]

first_union_exons <- Filter(function(ue) ue$exon_type == "first", suox_union$union_exons)

cat(sprintf("Actual union model has %d first exon groups\n\n", length(first_union_exons)))

for (i in seq_along(first_union_exons)) {
  ue <- first_union_exons[[i]]
  cat(sprintf("Union Exon %d (type: %s): %d variants\n",
              ue$exon_number, ue$exon_type, nrow(ue$variants)))
  for (j in seq_len(nrow(ue$variants))) {
    var <- ue$variants[j, ]
    cat(sprintf("  %s: %d - %d\n", var$isoform_id, var$start, var$end))
  }
  cat("\n")
}

# Check if target isoform is in ANY union exon
cat("═══ Step 7: Where is ENST00000886413.1 in union model? ═══\n\n")

found_in_union <- FALSE
for (i in seq_along(suox_union$union_exons)) {
  ue <- suox_union$union_exons[[i]]
  if (any(ue$variants$isoform_id == target_isoform)) {
    found_in_union <- TRUE
    cat(sprintf("Found in Union Exon %d (type: %s, isoform_exon: %s)\n",
                ue$exon_number, ue$exon_type,
                ue$variants %>% filter(isoform_id == target_isoform) %>% pull(exon_index)))
  }
}

if (!found_in_union) {
  cat(sprintf("✗ %s NOT found in ANY union exon!\n", target_isoform))
  cat("This isoform was completely excluded from the union model!\n")
}

cat("\n✓ Investigation complete\n")
