#!/usr/bin/env Rscript

################################################################################
# Script 04: Extract CDS Annotations
################################################################################
#
# Purpose:
#   Extract CDS (Coding Sequence) coordinates and metadata for all major isoforms.
#   This includes start codon, stop codon positions, and ORF length calculations.
#
# Input:
#   --source oarfish (default):
#     - data/expression_data.rds
#     - reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz
#     - reference_files/sqanti3_corrected.cds.gff3 (PacBio CDS)
#   --source isocall:
#     - data/isocall/expression_data.rds
#     - reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz (ENST CDS)
#     - SQANTI CDS output for novel isoforms (path TBD)
#
# Output:
#   - data/isoform_cds_metadata.rds (oarfish) or data/isocall/isoform_cds_metadata.rds
#
# Structure:
#   isoform_id, coding_status (coding/non_coding/unknown),
#   cds_start, cds_stop, orf_length
#
################################################################################

library(tidyverse)
library(rtracklayer)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse --source (default: oarfish)
source_type <- "oarfish"
if ("--source" %in% args) {
  src_idx <- which(args == "--source")
  if (src_idx < length(args)) {
    source_type <- args[src_idx + 1]
    if (!source_type %in% c("oarfish", "isocall")) {
      stop("--source must be 'oarfish' or 'isocall'")
    }
  }
}

# Set paths based on source
data_dir <- if (source_type == "isocall") "data/isocall" else "data"
output_file <- file.path(data_dir, "isoform_cds_metadata.rds")

# ==============================================================================
# Input Validation Helpers
# ==============================================================================

validate_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required input file not found: %s\nHave you run the prerequisite scripts?", path))
  }
}

validate_columns <- function(df, required_cols, name) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing required columns: %s", name, paste(missing, collapse = ", ")))
  }
}

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat(sprintf("║   STEP 4: Extract CDS Annotations [source: %s]%s║\n",
            source_type, strrep(" ", 15 - nchar(source_type))))
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ==============================================================================
# 0. Validate Inputs
# ==============================================================================

cat("Validating inputs...\n")
validate_file(file.path(data_dir, "expression_data.rds"))
validate_file("/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz")
if (source_type == "oarfish") {
  validate_file("/Users/petecastaldi/claude_projects/nmd/reference_files/sqanti3_corrected.cds.gff3")
}
# Note: isocall SQANTI CDS path validated below (TBD)
cat("  All required input files found.\n\n")

# ==============================================================================
# 1. Load Major Isoforms List
# ==============================================================================

cat("Loading major isoforms list...\n")
expression_data <- readRDS(file.path(data_dir, "expression_data.rds"))
validate_columns(expression_data, c("isoform_id"), "expression_data")
cat("  Required columns verified.\n")

major_isoforms <- unique(expression_data$isoform_id)
cat(sprintf("  Loaded %d major isoforms\n", length(major_isoforms)))

if (source_type == "isocall") {
  # Isocall: ENST (versioned) + novel (ENSG*.novelN)
  enst_isoforms <- major_isoforms[str_detect(major_isoforms, "^ENST")]
  novel_isoforms <- major_isoforms[str_detect(major_isoforms, "\\.novel\\d+$")]

  # For GENCODE CDS matching: need to strip version from isocall ENST IDs
  # (GENCODE GFF3 has unversioned ENST, isocall has versioned ENST)
  enst_no_version <- str_replace(enst_isoforms, "\\.\\d+$", "")

  cat(sprintf("  ENST isoforms: %d\n", length(enst_isoforms)))
  cat(sprintf("  Novel isoforms: %d\n", length(novel_isoforms)))
} else {
  # Oarfish: ENST + PB.
  enst_isoforms <- major_isoforms[str_detect(major_isoforms, "^ENST")]
  pb_isoforms <- major_isoforms[str_detect(major_isoforms, "^PB\\.")]

  # Strip version numbers for GENCODE only
  # IMPORTANT: PacBio IDs (PB.gene.isoform) are NOT versioned
  enst_no_version <- str_replace(enst_isoforms, "\\.\\d+$", "")
  pb_no_version <- pb_isoforms

  cat(sprintf("  GENCODE isoforms: %d\n", length(enst_isoforms)))
  cat(sprintf("  PacBio isoforms: %d\n", length(pb_isoforms)))
}

# ==============================================================================
# 2. Extract CDS from GENCODE
# ==============================================================================

cat("\nExtracting GENCODE CDS features...\n")
cat("  Loading GENCODE GFF3 (this may take several minutes)...\n")

gencode_gff <- import("/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz")

cat("  Filtering to CDS and codon features...\n")

# Filter to CDS and start_codon/stop_codon features for major isoforms only
gencode_cds <- gencode_gff %>%
  as_tibble() %>%
  filter(
    type %in% c("CDS", "start_codon", "stop_codon")
  ) %>%
  mutate(transcript_id_no_version = str_replace(transcript_id, "\\.\\d+$", "")) %>%
  filter(transcript_id_no_version %in% enst_no_version)

cat(sprintf("  Found %d CDS/codon features for %d isoforms\n",
            nrow(gencode_cds), n_distinct(gencode_cds$transcript_id)))

# Aggregate per isoform (using versionless IDs for grouping)
cat("  Calculating CDS metadata per isoform...\n")
gencode_cds_summary <- gencode_cds %>%
  group_by(transcript_id_no_version) %>%
  summarise(
    coding_status = if_else(any(type == "CDS"), "coding", "non_coding"),
    strand = dplyr::first(strand),

    # STRAND-AWARE CDS coordinate extraction:
    # cds_start = genomic position of ORF START (biological 5' end of CDS)
    # cds_stop = genomic position of ORF STOP (biological 3' end of CDS)
    #
    # + strand: start_codon at lower coords, stop_codon at higher coords
    # - strand: start_codon at higher coords, stop_codon at lower coords
    #
    # We store as min/max genomic coords for easier range checking:
    # cds_min = smallest genomic coordinate of CDS
    # cds_max = largest genomic coordinate of CDS

    cds_min = if_else(
      any(type == "start_codon") & any(type == "stop_codon"),
      pmin(
        min(start[type == "start_codon"]),
        min(start[type == "stop_codon"])
      ),
      if_else(
        any(type == "CDS"),
        min(start[type == "CDS"]),
        NA_integer_
      )
    ),
    cds_max = if_else(
      any(type == "start_codon") & any(type == "stop_codon"),
      pmax(
        max(end[type == "start_codon"]),
        max(end[type == "stop_codon"])
      ),
      if_else(
        any(type == "CDS"),
        max(end[type == "CDS"]),
        NA_integer_
      )
    ),

    # Calculate ORF length from CDS features
    orf_length = if_else(
      any(type == "CDS"),
      as.integer(sum(end[type == "CDS"] - start[type == "CDS"] + 1)),
      NA_integer_
    ),
    .groups = "drop"
  ) %>%
  select(-strand) %>%
  dplyr::rename(
    isoform_id_no_version = transcript_id_no_version,
    cds_start = cds_min,
    cds_stop = cds_max
  )

# Map back to versioned IDs
isoform_version_map <- tibble(
  isoform_id_no_version = enst_no_version,
  isoform_id = enst_isoforms
)

gencode_cds_summary <- gencode_cds_summary %>%
  left_join(isoform_version_map, by = "isoform_id_no_version") %>%
  select(isoform_id, coding_status, cds_start, cds_stop, orf_length)

cat(sprintf("  GENCODE CDS metadata: %d isoforms\n", nrow(gencode_cds_summary)))
cat(sprintf("    Coding: %d\n", sum(gencode_cds_summary$coding_status == "coding")))
cat(sprintf("    Non-coding: %d\n", sum(gencode_cds_summary$coding_status == "non_coding")))

# ==============================================================================
# 3. Extract CDS for Non-ENST Isoforms (PacBio or Novel)
# ==============================================================================

if (source_type == "isocall") {
  # --- Isocall novel isoforms: SQANTI CDS (TBD) ---
  if (length(novel_isoforms) > 0) {
    cat(sprintf("\nNovel isoforms requiring SQANTI CDS: %d\n", length(novel_isoforms)))

    isocall_sqanti_cds_file <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.cds.gff3"

    if (file.exists(isocall_sqanti_cds_file)) {
      cat("  Loading isocall SQANTI CDS GFF...\n")

      tryCatch({
        sqanti_gff <- import(isocall_sqanti_cds_file, format = "gtf")

        sqanti_cds <- sqanti_gff %>%
          as_tibble() %>%
          filter(
            type %in% c("CDS", "start_codon", "stop_codon"),
            transcript_id %in% novel_isoforms
          )

        cat(sprintf("  Found %d CDS/codon features for %d novel isoforms\n",
                    nrow(sqanti_cds), n_distinct(sqanti_cds$transcript_id)))

        novel_cds_summary <- sqanti_cds %>%
          group_by(transcript_id) %>%
          summarise(
            coding_status = if_else(any(type == "CDS"), "coding", "non_coding"),
            cds_start = if_else(
              any(type == "CDS"),
              min(start[type == "CDS"]),
              NA_integer_
            ),
            cds_stop = if_else(
              any(type == "CDS"),
              max(end[type == "CDS"]),
              NA_integer_
            ),
            orf_length = if_else(
              any(type == "CDS"),
              as.integer(sum(end[type == "CDS"] - start[type == "CDS"] + 1)),
              NA_integer_
            ),
            .groups = "drop"
          ) %>%
          dplyr::rename(isoform_id = transcript_id) %>%
          select(isoform_id, coding_status, cds_start, cds_stop, orf_length)

        cat(sprintf("  Novel CDS metadata: %d isoforms\n", nrow(novel_cds_summary)))
      },
      error = function(e) {
        cat(sprintf("  WARNING: Failed to load SQANTI CDS: %s\n", conditionMessage(e)))
        cat("  All novel isoforms will be marked as 'unknown'\n")
        novel_cds_summary <<- tibble(
          isoform_id = character(), coding_status = character(),
          cds_start = integer(), cds_stop = integer(), orf_length = integer()
        )
      })
    } else {
      cat(sprintf("  SQANTI CDS file not yet available: %s\n", isocall_sqanti_cds_file))
      cat("  All novel isoforms will be marked as 'unknown'\n")
      novel_cds_summary <- tibble(
        isoform_id = character(), coding_status = character(),
        cds_start = integer(), cds_stop = integer(), orf_length = integer()
      )
    }
  } else {
    novel_cds_summary <- tibble(
      isoform_id = character(), coding_status = character(),
      cds_start = integer(), cds_stop = integer(), orf_length = integer()
    )
  }
  # Use novel_cds_summary as the "non-ENST" CDS (analogous to pacbio_cds_summary)
  pacbio_cds_summary <- novel_cds_summary

} else {
  # --- Oarfish: PacBio CDS from SQANTI ---
  if (length(pb_isoforms) > 0) {
    cat("\nExtracting PacBio CDS features...\n")
    cat("  Loading SQANTI CDS GFF (this may take a moment)...\n")

    tryCatch({
      # Note: Despite .gff3 extension, this file is in GTF format (uses quotes)
      pacbio_gff <- import("/Users/petecastaldi/claude_projects/nmd/reference_files/sqanti3_corrected.cds.gff3",
                           format = "gtf")

      pacbio_cds <- pacbio_gff %>%
        as_tibble() %>%
        filter(
          type %in% c("CDS", "start_codon", "stop_codon"),
          transcript_id %in% pb_isoforms
        )

      cat(sprintf("  Found %d CDS/codon features for %d PacBio isoforms\n",
                  nrow(pacbio_cds), n_distinct(pacbio_cds$transcript_id)))

      pacbio_cds_summary <- pacbio_cds %>%
        group_by(transcript_id) %>%
        summarise(
          coding_status = if_else(any(type == "CDS"), "coding", "non_coding"),
          cds_start = if_else(
            any(type == "CDS"),
            min(start[type == "CDS"]),
            NA_integer_
          ),
          cds_stop = if_else(
            any(type == "CDS"),
            max(end[type == "CDS"]),
            NA_integer_
          ),
          orf_length = if_else(
            any(type == "CDS"),
            as.integer(sum(end[type == "CDS"] - start[type == "CDS"] + 1)),
            NA_integer_
          ),
          .groups = "drop"
        ) %>%
        dplyr::rename(isoform_id = transcript_id) %>%
        select(isoform_id, coding_status, cds_start, cds_stop, orf_length)

      cat(sprintf("  PacBio CDS metadata: %d isoforms\n", nrow(pacbio_cds_summary)))
    },
    error = function(e) {
      cat("  WARNING: Failed to load PacBio GFF, skipping PacBio CDS extraction\n")
      pacbio_cds_summary <<- tibble(
        isoform_id = character(), coding_status = character(),
        cds_start = integer(), cds_stop = integer(), orf_length = integer()
      )
    })
  } else {
    cat("\nNo PacBio isoforms to extract (skipping)\n")
    pacbio_cds_summary <- tibble(
      isoform_id = character(), coding_status = character(),
      cds_start = integer(), cds_stop = integer(), orf_length = integer()
    )
  }
}

# ==============================================================================
# 4. Combine and Handle Missing Isoforms
# ==============================================================================

cat("\nCombining CDS metadata...\n")

# Combine GENCODE and PacBio
cds_metadata <- bind_rows(
  gencode_cds_summary,
  pacbio_cds_summary
)

# Add missing isoforms as "unknown"
missing_isoforms <- setdiff(major_isoforms, cds_metadata$isoform_id)

if (length(missing_isoforms) > 0) {
  cat(sprintf("  Adding %d isoforms with unknown coding status\n",
              length(missing_isoforms)))

  missing_df <- tibble(
    isoform_id = missing_isoforms,
    coding_status = "unknown",
    cds_start = NA_integer_,
    cds_stop = NA_integer_,
    orf_length = NA_integer_
  )

  cds_metadata <- bind_rows(cds_metadata, missing_df)
}

# ==============================================================================
# 5. Summary Statistics
# ==============================================================================

cat("\nCDS metadata statistics:\n")
cat(sprintf("  Total isoforms: %d\n", nrow(cds_metadata)))
cat(sprintf("  Coding: %d (%.1f%%)\n",
            sum(cds_metadata$coding_status == "coding"),
            100 * mean(cds_metadata$coding_status == "coding")))
cat(sprintf("  Non-coding: %d (%.1f%%)\n",
            sum(cds_metadata$coding_status == "non_coding"),
            100 * mean(cds_metadata$coding_status == "non_coding")))
cat(sprintf("  Unknown: %d (%.1f%%)\n",
            sum(cds_metadata$coding_status == "unknown"),
            100 * mean(cds_metadata$coding_status == "unknown")))

# ORF length distribution (for coding isoforms)
coding_isoforms <- cds_metadata %>%
  filter(coding_status == "coding", !is.na(orf_length))

if (nrow(coding_isoforms) > 0) {
  cat(sprintf("\n  ORF length (coding isoforms):\n"))
  cat(sprintf("    Mean: %.0f bp\n", mean(coding_isoforms$orf_length)))
  cat(sprintf("    Median: %.0f bp\n", median(coding_isoforms$orf_length)))
  cat(sprintf("    Range: %d - %d bp\n",
              min(coding_isoforms$orf_length),
              max(coding_isoforms$orf_length)))
}

# ==============================================================================
# 6. Save Output
# ==============================================================================

cat("\nSaving CDS metadata...\n")
saveRDS(cds_metadata, output_file)
cat(sprintf("  ✓ %s\n", output_file))

cat("\n✓ Step 4 complete\n")
cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Isoforms annotated: %d\n", nrow(cds_metadata)))
cat(sprintf("  Coding: %d\n", sum(cds_metadata$coding_status == "coding")))
cat(sprintf("  Non-coding: %d\n", sum(cds_metadata$coding_status == "non_coding")))
cat(sprintf("  Unknown: %d\n", sum(cds_metadata$coding_status == "unknown")))
if (nrow(coding_isoforms) > 0) {
  cat(sprintf("Mean ORF length: %.0f bp\n", mean(coding_isoforms$orf_length)))
}
cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")
