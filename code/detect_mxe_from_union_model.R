#!/usr/bin/env Rscript
# Mutually Exclusive Exon (MXE) Detection from Union Model
# Identifies exons at same structural position that never co-occur

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   MXE DETECTION FROM UNION EXON MODEL                        ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Configuration
output_dir <- "results/isoform_transitions/v3.0_reference_based"

# ============================================================================
# Helper Functions
# ============================================================================

#' Build isoform structure with flanking exons
#' @param union_exons Union exon model
#' @param isoform_id Target isoform
#' @return Data frame with exon numbers and flanking exon info
build_isoform_structure <- function(union_exons, isoform_id) {

  exon_records <- map_dfr(union_exons, function(union_exon) {
    variants_in_iso <- union_exon$variants %>%
      filter(isoform_id == !!isoform_id)

    if (nrow(variants_in_iso) == 0) {
      return(NULL)
    }

    variant <- variants_in_iso[1, ]

    tibble(
      exon_number = union_exon$exon_number,
      start = variant$start,
      end = variant$end,
      is_first = variant$is_first,
      is_last = variant$is_last
    )
  })

  if (nrow(exon_records) == 0) return(NULL)

  # Add flanking exon numbers
  exon_records <- exon_records %>%
    arrange(exon_number) %>%
    mutate(
      upstream_exon = lag(exon_number),
      downstream_exon = lead(exon_number)
    )

  exon_records
}

#' Detect MXE events between two isoforms
#' @param union_exons Union exon model
#' @param isoform_A_id First isoform
#' @param isoform_B_id Second isoform (reference)
#' @return MXE event records
detect_mxe_in_transition <- function(union_exons, isoform_A_id, isoform_B_id) {

  # Build structures for both isoforms
  struct_A <- build_isoform_structure(union_exons, isoform_A_id)
  struct_B <- build_isoform_structure(union_exons, isoform_B_id)

  if (is.null(struct_A) || is.null(struct_B)) {
    return(NULL)
  }

  mxe_events <- list()

  # Look for exons with same flanking exons but different union numbers
  # This indicates mutually exclusive exons

  for (i in seq_len(nrow(struct_A))) {
    exon_A <- struct_A[i, ]

    # Skip first and last exons (can't have MXE)
    if (exon_A$is_first || exon_A$is_last) next
    if (is.na(exon_A$upstream_exon) || is.na(exon_A$downstream_exon)) next

    # Find exons in B with same flanking exons but different union number
    for (j in seq_len(nrow(struct_B))) {
      exon_B <- struct_B[j, ]

      # Skip if same exon number (not mutually exclusive)
      if (exon_A$exon_number == exon_B$exon_number) next

      # Skip first and last
      if (exon_B$is_first || exon_B$is_last) next
      if (is.na(exon_B$upstream_exon) || is.na(exon_B$downstream_exon)) next

      # Check for same flanking exons
      same_upstream <- exon_A$upstream_exon == exon_B$upstream_exon
      same_downstream <- exon_A$downstream_exon == exon_B$downstream_exon

      if (same_upstream && same_downstream) {
        # Found MXE!
        # Exon A and Exon B have same flanking exons but different positions

        mxe_events[[length(mxe_events) + 1]] <- tibble(
          event_type = "MXE",
          exon_A_number = exon_A$exon_number,
          exon_B_number = exon_B$exon_number,
          upstream_exon = exon_A$upstream_exon,
          downstream_exon = exon_A$downstream_exon,
          direction = "mutual",
          detail = sprintf(
            "MXE: Exon %d (in alt) vs Exon %d (in ref), flanked by %d and %d",
            exon_A$exon_number, exon_B$exon_number,
            exon_A$upstream_exon, exon_A$downstream_exon
          )
        )
      }
    }
  }

  if (length(mxe_events) == 0) return(NULL)

  bind_rows(mxe_events) %>% distinct()
}

#' Verify MXE across all gene isoforms
#' @param union_exons Union exon model
#' @param exon_A_num First exon number
#' @param exon_B_num Second exon number
#' @param all_isoforms List of all isoform IDs for this gene
#' @return TRUE if exons A and B never co-occur
verify_mutual_exclusivity <- function(union_exons, exon_A_num, exon_B_num, all_isoforms) {

  for (iso_id in all_isoforms) {
    # Get this isoform's exons
    iso_exons <- map_dbl(union_exons, function(ue) {
      variants <- ue$variants %>% filter(isoform_id == iso_id)
      if (nrow(variants) > 0) return(ue$exon_number) else return(NA_real_)
    })
    iso_exons <- iso_exons[!is.na(iso_exons)]

    # Check if both exons present
    has_A <- exon_A_num %in% iso_exons
    has_B <- exon_B_num %in% iso_exons

    if (has_A && has_B) {
      # Both exons present in this isoform - NOT mutually exclusive
      return(FALSE)
    }
  }

  # Checked all isoforms, never found both together
  return(TRUE)
}

# ============================================================================
# Main Processing
# ============================================================================

cat("Loading data...\n")

# Load union exon models
union_models <- readRDS(file.path(output_dir, "union_exon_models_test.rds"))
cat("  Loaded:", length(union_models), "genes with union models\n")

# Load major isoforms data
major_isoforms <- readRDS(file.path(output_dir, "reference_event_vectors_v3.0_filtered.rds"))
cat("  Loaded:", nrow(major_isoforms), "isoform transitions\n\n")

cat("Detecting MXE events...\n\n")

# Store MXE results
all_mxe_results <- list()
mxe_idx <- 1

for (gene_id in names(union_models)) {
  model <- union_models[[gene_id]]

  cat(sprintf("[%d/%d] Processing: %s\n",
              which(names(union_models) == gene_id),
              length(union_models),
              gene_id))

  # Get all transitions for this gene
  gene_transitions <- major_isoforms %>% filter(gene_id == !!gene_id)

  if (nrow(gene_transitions) == 0) {
    cat("  No transitions found, skipping\n\n")
    next
  }

  # Get unique isoforms
  all_isoforms <- unique(c(gene_transitions$isoform_ref, gene_transitions$isoform_alt))

  cat(sprintf("  Isoforms: %d, Transitions: %d\n",
              length(all_isoforms), nrow(gene_transitions)))

  n_mxe_found <- 0

  # Process each transition
  for (i in seq_len(nrow(gene_transitions))) {
    trans <- gene_transitions[i, ]

    # Detect MXE between isoform_alt (A) and isoform_ref (B)
    mxe_events <- detect_mxe_in_transition(
      model$union_exons,
      trans$isoform_alt,
      trans$isoform_ref
    )

    if (!is.null(mxe_events) && nrow(mxe_events) > 0) {

      # Verify each MXE event across all gene isoforms
      verified_mxe <- list()

      for (j in seq_len(nrow(mxe_events))) {
        mxe <- mxe_events[j, ]

        is_exclusive <- verify_mutual_exclusivity(
          model$union_exons,
          mxe$exon_A_number,
          mxe$exon_B_number,
          all_isoforms
        )

        if (is_exclusive) {
          verified_mxe[[length(verified_mxe) + 1]] <- mxe %>%
            mutate(verified = TRUE)
          n_mxe_found <- n_mxe_found + 1
        }
      }

      if (length(verified_mxe) > 0) {
        all_mxe_results[[mxe_idx]] <- bind_rows(verified_mxe) %>%
          mutate(
            gene_id = gene_id,
            isoform_ref = trans$isoform_ref,
            isoform_alt = trans$isoform_alt,
            transition_id = paste(trans$isoform_ref, trans$isoform_alt, sep = "->")
          )
        mxe_idx <- mxe_idx + 1
      }
    }
  }

  cat(sprintf("  MXE events found: %d\n\n", n_mxe_found))
}

# Combine results
if (length(all_mxe_results) > 0) {
  mxe_results <- bind_rows(all_mxe_results)

  cat("\n═══ MXE Detection Summary ═══\n\n")
  cat("Total MXE events found:", nrow(mxe_results), "\n")
  cat("Genes with MXE:", length(unique(mxe_results$gene_id)), "\n")
  cat("Transitions with MXE:", length(unique(mxe_results$transition_id)), "\n\n")

  # Summary by gene
  mxe_by_gene <- mxe_results %>%
    group_by(gene_id) %>%
    summarise(
      n_mxe_events = n(),
      n_transitions_with_mxe = n_distinct(transition_id),
      .groups = "drop"
    ) %>%
    arrange(desc(n_mxe_events))

  cat("MXE events per gene:\n")
  print(mxe_by_gene)
  cat("\n")

  # Show examples
  cat("═══ Example MXE Events ═══\n\n")
  for (i in 1:min(5, nrow(mxe_results))) {
    cat(sprintf("Example %d:\n", i))
    cat("  Gene:", mxe_results$gene_id[i], "\n")
    cat("  Transition:", mxe_results$transition_id[i], "\n")
    cat("  ", mxe_results$detail[i], "\n\n")
  }

  # Save results
  cat("Saving results...\n")
  saveRDS(mxe_results, file.path(output_dir, "mxe_events_test.rds"))
  write_tsv(mxe_results, file.path(output_dir, "mxe_events_test.tsv"))
  cat("  Saved: mxe_events_test.rds\n")
  cat("  Saved: mxe_events_test.tsv\n\n")

} else {
  cat("\n═══ No MXE Events Found ═══\n\n")
  cat("No mutually exclusive exons detected in test dataset.\n")
  cat("This could indicate:\n")
  cat("  1. Test genes don't have MXE patterns\n")
  cat("  2. Union model groups MXE exons together (shared boundaries)\n")
  cat("  3. MXE is rare in these particular genes\n\n")
}

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   MXE DETECTION COMPLETE                                      ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")
