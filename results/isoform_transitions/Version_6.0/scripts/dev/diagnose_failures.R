#!/usr/bin/env Rscript
################################################################################
# Diagnose first 20 reconstruction failures
################################################################################

library(tidyverse)

base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"
setwd(base_dir)

source("scripts/core/reconstruction_functions.R")
source("scripts/core/event_detection_functions.R")

# Load data
verification <- readRDS("data/reconstruction_verification.rds")
profiles <- readRDS("data/splicing_choice_profiles.rds")
structures_compact <- readRDS("data/isoform_structures_filtered.rds")
union_exons_raw <- readRDS("data/union_exons_filtered.rds")

# Expand isoform structures
isoform_structures <- structures_compact %>%
  rowwise() %>%
  mutate(
    exon_data = list(tibble(
      exon_number = seq_along(exon_starts),
      exon_start = exon_starts,
      exon_end = exon_ends
    ))
  ) %>%
  unnest(exon_data) %>%
  select(isoform_id, gene_id, seqnames, strand, exon_number, exon_start, exon_end) %>%
  ungroup()

# Rename union exons
union_exons <- union_exons_raw %>%
  rename(start = union_exon_start, end = union_exon_end, chr = seqnames)

# Extract all events
all_events <- profiles %>%
  filter(n_events > 0) %>%
  select(dominant_isoform_id, non_dominant_isoform_id, detailed_events) %>%
  unnest(detailed_events)

# Get failures
failures <- verification %>% filter(status == "FAIL")
cat(sprintf("Total failures: %d\n\n", nrow(failures)))

# Analyze first 20
n_analyze <- min(20, nrow(failures))

for (idx in 1:n_analyze) {
  f <- failures[idx, ]
  dom_id <- f$dominant_isoform_id
  comp_id <- f$comparator_isoform_id
  gene <- f$gene_id

  cat(sprintf("════════════════════════════════════════════════════════════════\n"))
  cat(sprintf("FAILURE %d: %s\n", idx, f$reason))
  cat(sprintf("  Gene: %s\n", gene))
  cat(sprintf("  Dominant: %s\n", dom_id))
  cat(sprintf("  Comparator: %s\n", comp_id))

  # Get exon structures
  dom_exons <- isoform_structures %>%
    filter(isoform_id == dom_id) %>%
    arrange(exon_start)
  comp_exons <- isoform_structures %>%
    filter(isoform_id == comp_id) %>%
    arrange(exon_start)

  strand <- dom_exons$strand[1]

  cat(sprintf("  Strand: %s\n", strand))
  cat(sprintf("  Dom exons: %d  [%d - %d]\n",
              nrow(dom_exons), min(dom_exons$exon_start), max(dom_exons$exon_end)))
  cat(sprintf("  Comp exons: %d  [%d - %d]\n",
              nrow(comp_exons), min(comp_exons$exon_start), max(comp_exons$exon_end)))

  # Get events for this pair
  pair_events <- all_events %>%
    filter(dominant_transcript_id == dom_id,
           comparator_transcript_id == comp_id)

  cat(sprintf("  Events: %d\n", nrow(pair_events)))

  if (nrow(pair_events) > 0) {
    for (e in 1:nrow(pair_events)) {
      ev <- pair_events[e, ]
      cat(sprintf("    [%d] %s %s: %d-%d (bp_diff=%d)\n",
                  e, ev$event_type, ev$direction,
                  ev$five_prime, ev$three_prime,
                  ev$bp_diff))
      # Show ir_split_exons if present
      if (!is.null(ev$ir_split_exons) && !is.na(ev$ir_split_exons) &&
          is.character(ev$ir_split_exons) && nchar(ev$ir_split_exons) > 0) {
        cat(sprintf("         ir_split_exons: %s\n", ev$ir_split_exons))
      }
      if (!is.null(ev$missing_terminal_exons) && !is.na(ev$missing_terminal_exons) &&
          is.character(ev$missing_terminal_exons) && nchar(ev$missing_terminal_exons) > 0) {
        cat(sprintf("         missing_terminal: %s\n", ev$missing_terminal_exons))
      }
      if (!is.null(ev$orphan_terminal_exons) && !is.na(ev$orphan_terminal_exons) &&
          is.character(ev$orphan_terminal_exons) && nchar(ev$orphan_terminal_exons) > 0) {
        cat(sprintf("         orphan_terminal: %s\n", ev$orphan_terminal_exons))
      }
    }
  }

  # Reconstruct and show result
  comp_exons_for_recon <- isoform_structures %>%
    filter(isoform_id == comp_id) %>%
    rename(chr = seqnames, transcript_id = isoform_id) %>%
    select(chr, exon_start, exon_end, strand, gene_id, transcript_id)

  gene_ue <- union_exons %>% filter(gene_id == gene)

  reconstructed <- tryCatch({
    reconstruct_dominant_v2(comp_exons_for_recon, pair_events, gene_ue)
  }, error = function(e) {
    cat(sprintf("  RECONSTRUCTION ERROR: %s\n", e$message))
    tibble()
  })

  if (nrow(reconstructed) > 0) {
    recon_sorted <- reconstructed %>% arrange(exon_start)
    dom_sorted <- dom_exons %>% arrange(exon_start)

    cat(sprintf("\n  Reconstructed exons: %d\n", nrow(recon_sorted)))

    # Side-by-side comparison
    max_rows <- max(nrow(dom_sorted), nrow(recon_sorted))
    cat(sprintf("  %-5s %-25s  %-25s  %s\n", "Exon", "Original Dominant", "Reconstructed", "Match?"))
    cat(sprintf("  %-5s %-25s  %-25s  %s\n", "----", "-------------------", "-------------------", "------"))

    for (r in 1:max_rows) {
      if (r <= nrow(dom_sorted)) {
        orig_str <- sprintf("[%d - %d]", dom_sorted$exon_start[r], dom_sorted$exon_end[r])
      } else {
        orig_str <- "(missing)"
      }
      if (r <= nrow(recon_sorted)) {
        recon_str <- sprintf("[%d - %d]", recon_sorted$exon_start[r], recon_sorted$exon_end[r])
      } else {
        recon_str <- "(missing)"
      }

      match_str <- ""
      if (r <= nrow(dom_sorted) && r <= nrow(recon_sorted)) {
        s_diff <- abs(dom_sorted$exon_start[r] - recon_sorted$exon_start[r])
        e_diff <- abs(dom_sorted$exon_end[r] - recon_sorted$exon_end[r])
        if (s_diff == 0 && e_diff == 0) {
          match_str <- "OK"
        } else {
          match_str <- sprintf("DIFF s=%d e=%d", s_diff, e_diff)
        }
      } else {
        match_str <- "EXTRA/MISSING"
      }

      cat(sprintf("  %-5d %-25s  %-25s  %s\n", r, orig_str, recon_str, match_str))
    }
  }

  # Classify failure type
  exon_diff <- f$n_exons_reconstructed - f$n_exons_original
  is_dom_novel <- grepl("^PB\\.", dom_id)
  is_comp_novel <- grepl("^PB\\.", comp_id)

  has_ir <- any(grepl("^IR", pair_events$event_type))
  has_ir_diff <- any(grepl("^IR_diff", pair_events$event_type))
  has_se <- any(pair_events$event_type == "SE")
  has_missing_internal <- any(pair_events$event_type == "Missing_Internal")
  has_terminal <- any(pair_events$event_type %in% c("Alt_TSS", "Alt_TES"))

  event_types_present <- paste(unique(pair_events$event_type), collapse=", ")

  cat(sprintf("\n  Classification:\n"))
  cat(sprintf("    Exon count diff: %+d (recon - orig)\n", exon_diff))
  cat(sprintf("    Dom is PacBio: %s\n", is_dom_novel))
  cat(sprintf("    Comp is PacBio: %s\n", is_comp_novel))
  cat(sprintf("    Has IR: %s | IR_diff: %s | SE: %s | Missing_Internal: %s | Terminal: %s\n",
              has_ir, has_ir_diff, has_se, has_missing_internal, has_terminal))
  cat(sprintf("    All event types: %s\n", event_types_present))
  cat("\n")
}

# Summary classification of ALL failures
cat("\n\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("AGGREGATE FAILURE CLASSIFICATION (all 927 failures)\n")
cat("════════════════════════════════════════════════════════════════\n\n")

failure_details <- failures %>%
  left_join(
    all_events %>%
      group_by(dominant_transcript_id, comparator_transcript_id) %>%
      summarise(
        n_events = n(),
        event_types = paste(unique(event_type), collapse=","),
        has_IR = any(grepl("^IR$", event_type)),
        has_IR_diff = any(grepl("^IR_diff", event_type)),
        has_SE = any(event_type == "SE"),
        has_Missing_Internal = any(event_type == "Missing_Internal"),
        has_terminal = any(event_type %in% c("Alt_TSS", "Alt_TES")),
        has_A5SS = any(event_type == "A5SS"),
        has_A3SS = any(event_type == "A3SS"),
        has_Partial_IR = any(grepl("Partial_IR", event_type)),
        .groups = "drop"
      ),
    by = c("dominant_isoform_id" = "dominant_transcript_id",
           "comparator_isoform_id" = "comparator_transcript_id")
  ) %>%
  mutate(
    exon_diff = n_exons_reconstructed - n_exons_original,
    dom_is_pb = grepl("^PB\\.", dominant_isoform_id),
    comp_is_pb = grepl("^PB\\.", comparator_isoform_id),
    direction = case_when(
      exon_diff > 0 ~ "too_many_exons",
      exon_diff < 0 ~ "too_few_exons",
      TRUE ~ "same_count_wrong_coords"
    )
  )

# Direction of mismatch
cat("Exon count direction:\n")
print(failure_details %>% count(direction, sort = TRUE))

cat("\nDom is PacBio vs GENCODE:\n")
print(failure_details %>% count(dom_is_pb, sort = TRUE))

cat("\nComp is PacBio vs GENCODE:\n")
print(failure_details %>% count(comp_is_pb, sort = TRUE))

cat("\nEvent type prevalence in failures:\n")
event_prev <- tibble(
  event_type = c("IR", "IR_diff", "SE", "Missing_Internal", "Alt_TSS/TES", "A5SS", "A3SS", "Partial_IR"),
  n_failures = c(
    sum(failure_details$has_IR, na.rm = TRUE),
    sum(failure_details$has_IR_diff, na.rm = TRUE),
    sum(failure_details$has_SE, na.rm = TRUE),
    sum(failure_details$has_Missing_Internal, na.rm = TRUE),
    sum(failure_details$has_terminal, na.rm = TRUE),
    sum(failure_details$has_A5SS, na.rm = TRUE),
    sum(failure_details$has_A3SS, na.rm = TRUE),
    sum(failure_details$has_Partial_IR, na.rm = TRUE)
  ),
  pct = round(100 * n_failures / nrow(failure_details), 1)
)
print(event_prev %>% arrange(desc(n_failures)))

cat("\nEvent type combinations (top 15):\n")
combo_counts <- failure_details %>%
  count(event_types, sort = TRUE) %>%
  head(15)
print(combo_counts, n = 15)

cat("\nExon count diff distribution:\n")
diff_dist <- failure_details %>%
  count(exon_diff, sort = TRUE) %>%
  head(15)
print(diff_dist, n = 15)

# Cross-tab: direction × has_IR
cat("\nDirection × has_IR:\n")
print(failure_details %>% count(direction, has_IR, sort = TRUE))

# How many failures have NO events at all?
no_events <- failure_details %>% filter(is.na(n_events))
cat(sprintf("\nFailures with NO events detected: %d\n", nrow(no_events)))

cat("\n════════════════════════════════════════════════════════════════\n")
cat("Done.\n")
