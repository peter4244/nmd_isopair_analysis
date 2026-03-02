#!/usr/bin/env Rscript

################################################################################
# Script 03: Build Atomic Union Exon Models
################################################################################
#
# Purpose:
#   Create atomic (non-overlapping) union exon segments by splitting at all
#   unique exon boundaries within each gene. Each atomic segment represents the
#   smallest unit that is either fully included or fully excluded in any given
#   isoform. This provides a common coordinate system for comparing isoform
#   structures and is required by the reconstruction pipeline.
#
# Algorithm:
#   For each gene, collect all unique exon start/end coordinates. Classify each
#   boundary as START, END, or BOTH. Split into atomic segments using
#   boundary-type-aware logic, then filter to segments covered by at least one
#   exon. This produces fine-grained segments at every exon boundary transition.
#
# Performance:
#   - Pre-splits isoform_structures by gene (avoids repeated filter on full df)
#   - Vectorized coverage check (no rowwise())
#   - Vectorized containment join for isoform-UE mapping
#   - Writes results incrementally to TSV (constant memory)
#   - Produces both RDS and tabix-indexed TSV outputs
#
# Input:
#   - data/isoform_structures.rds
#
# Output:
#   - data/union_exons.tsv.gz          (tabix-indexed, for random gene-level access)
#   - data/union_exons.tsv.gz.tbi      (tabix index)
#   - data/union_exons.rds             (R native, for bulk loading)
#   - data/isoform_union_mapping.rds   (which union exons each isoform uses)
#
# Usage:
#   Rscript scripts/03_build_union_exons.R              # full run
#   Rscript scripts/03_build_union_exons.R --test 1000  # first 1000 genes only
#
################################################################################

library(tidyverse)

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

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
n_genes_limit <- NULL
if ("--test" %in% args) {
  test_idx <- which(args == "--test")
  if (test_idx < length(args)) {
    n_genes_limit <- as.integer(args[test_idx + 1])
  } else {
    n_genes_limit <- 1000L  # default test size
  }
}

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

# Set data directory based on source
data_dir <- if (source_type == "isocall") "data/isocall" else "data"

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   STEP 3: Build Atomic Union Exon Models                     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

if (!is.null(n_genes_limit)) {
  cat(sprintf("*** TEST MODE: Processing first %d genes only ***\n\n", n_genes_limit))
}

# ==============================================================================
# 0. Validate Inputs
# ==============================================================================

cat("Validating inputs...\n")
validate_file(file.path(data_dir, "isoform_structures.rds"))
cat("  All required input files found.\n\n")

# ==============================================================================
# 1. Load Isoform Structures
# ==============================================================================

cat("Loading isoform structures...\n")
isoform_structures <- readRDS(file.path(data_dir, "isoform_structures.rds"))
validate_columns(isoform_structures,
                 c("gene_id", "isoform_id", "seqnames", "strand", "exon_starts", "exon_ends"),
                 "isoform_structures")
cat("  Required columns verified.\n")

cat(sprintf("  Loaded %d isoforms from %d genes\n",
            nrow(isoform_structures),
            n_distinct(isoform_structures$gene_id)))

# ==============================================================================
# 2. Pre-split by gene (avoids repeated filter on full 500K-row tibble)
# ==============================================================================

cat("\nPre-splitting isoform structures by gene...\n")
gene_list <- split(isoform_structures, isoform_structures$gene_id)
genes <- names(gene_list)

if (!is.null(n_genes_limit)) {
  genes <- genes[seq_len(min(n_genes_limit, length(genes)))]
  gene_list <- gene_list[genes]
  cat(sprintf("  %d genes selected (test mode, %d total available)\n",
              length(genes), n_distinct(isoform_structures$gene_id)))
} else {
  cat(sprintf("  %d genes to process\n", length(genes)))
}

# ==============================================================================
# 3. Build Atomic Union Exons Per Gene (incremental output)
# ==============================================================================

cat("\nBuilding atomic union exon models...\n")

# Temporary TSV files for incremental writing
ue_tmp_file <- file.path(data_dir, "union_exons_tmp.tsv")
map_tmp_file <- file.path(data_dir, "isoform_union_mapping_tmp.tsv")

# Write headers
ue_header <- "seqnames\tunion_exon_start\tunion_exon_end\tunion_exon_id\tstrand\tgene_id\tunion_exon_number\tunion_exon_length"
map_header <- "union_exon_id\tunion_exon_number\tunion_exon_start\tunion_exon_end\tunion_exon_length\tgene_id\tseqnames\tstrand\tisoform_id\texon_number_in_isoform\tisoform_exon_start\tisoform_exon_end"
writeLines(ue_header, ue_tmp_file)
writeLines(map_header, map_tmp_file)

# Open file connections for appending
ue_con <- file(ue_tmp_file, open = "a")
map_con <- file(map_tmp_file, open = "a")

pb <- txtProgressBar(min = 0, max = length(genes), style = 3)

total_ues <- 0L
total_mappings <- 0L
genes_processed <- 0L

for (i in seq_along(genes)) {
  gene <- genes[i]
  gene_isoforms <- gene_list[[gene]]

  if (nrow(gene_isoforms) == 0) next

  chr <- gene_isoforms$seqnames[1]
  strand_val <- gene_isoforms$strand[1]

  # --------------------------------------------------------------------------
  # Expand list columns to per-exon rows (vectorized)
  # --------------------------------------------------------------------------
  n_exons_per_iso <- lengths(gene_isoforms$exon_starts)
  all_exon_starts <- unlist(gene_isoforms$exon_starts)
  all_exon_ends <- unlist(gene_isoforms$exon_ends)
  all_iso_ids <- rep(gene_isoforms$isoform_id, n_exons_per_iso)
  all_exon_nums <- unlist(lapply(n_exons_per_iso, seq_len))

  if (length(all_exon_starts) == 0) next

  # --------------------------------------------------------------------------
  # Atomic boundary splitting
  # --------------------------------------------------------------------------

  # Unique sorted boundaries
  all_boundaries <- sort(unique(c(all_exon_starts, all_exon_ends)))
  if (length(all_boundaries) < 2) next

  # Classify boundaries using set lookups (vectorized)
  starts_set <- unique(all_exon_starts)
  ends_set <- unique(all_exon_ends)
  is_start <- all_boundaries %in% starts_set
  is_end <- all_boundaries %in% ends_set

  # Build segments — pre-allocate at max possible size
  n_bounds <- length(all_boundaries)
  # Max segments: n_bounds - 1 (normal) + number of BOTH boundaries (extra split)
  n_both <- sum(is_start & is_end)
  max_segs <- n_bounds - 1 + n_both
  seg_s <- integer(max_segs)
  seg_e <- integer(max_segs)
  n_seg <- 0L

  current_start <- all_boundaries[1]

  if (n_bounds > 2) {
    for (k in 2:(n_bounds - 1)) {
      b <- all_boundaries[k]
      if (is_start[k] && is_end[k]) {
        # BOTH: three segments [..., b-1], [b, b], [b+1, ...]
        n_seg <- n_seg + 1L; seg_s[n_seg] <- current_start; seg_e[n_seg] <- b - 1L
        n_seg <- n_seg + 1L; seg_s[n_seg] <- b;             seg_e[n_seg] <- b
        current_start <- b + 1L
      } else if (is_start[k]) {
        # START: split before
        n_seg <- n_seg + 1L; seg_s[n_seg] <- current_start; seg_e[n_seg] <- b - 1L
        current_start <- b
      } else {
        # END: split after
        n_seg <- n_seg + 1L; seg_s[n_seg] <- current_start; seg_e[n_seg] <- b
        current_start <- b + 1L
      }
    }
  }

  # Final segment
  n_seg <- n_seg + 1L; seg_s[n_seg] <- current_start; seg_e[n_seg] <- all_boundaries[n_bounds]

  # Trim to actual size and filter phantoms (start > end)
  seg_s <- seg_s[1:n_seg]
  seg_e <- seg_e[1:n_seg]
  valid <- seg_s <= seg_e
  seg_s <- seg_s[valid]
  seg_e <- seg_e[valid]

  if (length(seg_s) == 0) next

  # --------------------------------------------------------------------------
  # Vectorized coverage check: is each segment covered by at least one exon?
  # Uses outer comparison — fine for per-gene sizes (typically small)
  # --------------------------------------------------------------------------
  n_segs <- length(seg_s)
  n_exons <- length(all_exon_starts)

  # For each segment, check if any exon fully contains it
  # seg_s[j] >= exon_start[k] AND seg_e[j] <= exon_end[k]
  covered <- logical(n_segs)
  for (j in seq_len(n_segs)) {
    covered[j] <- any(all_exon_starts <= seg_s[j] & all_exon_ends >= seg_e[j])
  }

  seg_s <- seg_s[covered]
  seg_e <- seg_e[covered]

  if (length(seg_s) == 0) next

  # --------------------------------------------------------------------------
  # Create union exon records
  # --------------------------------------------------------------------------
  n_ue <- length(seg_s)
  ue_ids <- paste0(gene, "_UE", seq_len(n_ue))
  ue_lengths <- seg_e - seg_s + 1L

  # Write union exons to TSV (tab-separated, no quotes)
  # Format: seqnames, start, end, union_exon_id, strand, gene_id, number, length
  ue_lines <- paste(chr, seg_s, seg_e, ue_ids, strand_val, gene,
                    seq_len(n_ue), ue_lengths, sep = "\t")
  writeLines(ue_lines, ue_con)

  total_ues <- total_ues + n_ue

  # --------------------------------------------------------------------------
  # Vectorized isoform-to-UE mapping via containment
  # For each isoform exon, find all UEs where ue_start >= exon_start & ue_end <= exon_end
  # --------------------------------------------------------------------------
  unique_isos <- unique(all_iso_ids)

  for (iso_id in unique_isos) {
    mask <- all_iso_ids == iso_id
    ex_starts <- all_exon_starts[mask]
    ex_ends <- all_exon_ends[mask]
    ex_nums <- all_exon_nums[mask]

    # For each exon in this isoform, find contained UEs
    map_lines <- character(0)
    for (e_idx in seq_along(ex_starts)) {
      contained_mask <- seg_s >= ex_starts[e_idx] & seg_e <= ex_ends[e_idx]
      if (any(contained_mask)) {
        which_contained <- which(contained_mask)
        # Format: ue_id, ue_number, ue_start, ue_end, ue_length, gene_id, chr, strand,
        #         isoform_id, exon_number, iso_exon_start, iso_exon_end
        new_lines <- paste(
          ue_ids[which_contained],
          which_contained,
          seg_s[which_contained],
          seg_e[which_contained],
          ue_lengths[which_contained],
          gene, chr, strand_val,
          iso_id, ex_nums[e_idx],
          ex_starts[e_idx], ex_ends[e_idx],
          sep = "\t"
        )
        map_lines <- c(map_lines, new_lines)
      }
    }

    if (length(map_lines) > 0) {
      writeLines(map_lines, map_con)
      total_mappings <- total_mappings + length(map_lines)
    }
  }

  genes_processed <- genes_processed + 1L

  if (i %% 1000 == 0) {
    setTxtProgressBar(pb, i)
  }
}

setTxtProgressBar(pb, length(genes))
close(pb)
close(ue_con)
close(map_con)

cat(sprintf("\n  Created %d atomic union exons across %d genes\n",
            total_ues, genes_processed))
cat(sprintf("  Created %d isoform-UE mappings\n", total_mappings))

# ==============================================================================
# 4. Create tabix-indexed union exon file
# ==============================================================================

cat("\nCreating tabix-indexed union exon file...\n")

# Read the TSV, sort by chr+start (required for tabix), write tabix format
ue_df <- read_tsv(ue_tmp_file, show_col_types = FALSE)

# Tabix file: chr, start, end, union_exon_id, strand, gene_id (matches development format)
tabix_file <- file.path(data_dir, "union_exons.tsv")
ue_sorted <- ue_df %>% arrange(seqnames, union_exon_start)

writeLines("#chr\tstart\tend\tunion_exon_id\tstrand\tgene_id", tabix_file)
tabix_lines <- paste(ue_sorted$seqnames, ue_sorted$union_exon_start,
                     ue_sorted$union_exon_end, ue_sorted$union_exon_id,
                     ue_sorted$strand, ue_sorted$gene_id, sep = "\t")
write(tabix_lines, tabix_file, append = TRUE)

# bgzip and tabix index
bgzip_status <- system(sprintf("bgzip -f %s", tabix_file))
if (bgzip_status != 0) {
  warning(sprintf("bgzip failed with status %d", bgzip_status))
} else {
  tabix_status <- system(sprintf("tabix -f -s 1 -b 2 -e 3 %s.gz", tabix_file))
  if (tabix_status != 0) {
    warning(sprintf("tabix failed with status %d", tabix_status))
  } else {
    cat(sprintf("  %s/union_exons.tsv.gz (tabix-indexed)\n", data_dir))
  }
}

# ==============================================================================
# 5. Save RDS outputs (for scripts that bulk-load)
# ==============================================================================

cat("\nSaving RDS outputs...\n")

# Union exons RDS (same columns as before for compatibility)
union_exons <- ue_df %>%
  select(union_exon_id, union_exon_number, union_exon_start, union_exon_end,
         union_exon_length, gene_id, seqnames, strand)

saveRDS(union_exons, file.path(data_dir, "union_exons.rds"))
cat(sprintf("  %s/union_exons.rds\n", data_dir))

# Isoform mapping RDS
map_df <- read_tsv(map_tmp_file, show_col_types = FALSE)

saveRDS(map_df, file.path(data_dir, "isoform_union_mapping.rds"))
cat(sprintf("  %s/isoform_union_mapping.rds\n", data_dir))

# ==============================================================================
# 6. Validation
# ==============================================================================

cat("\nValidating atomic union exons...\n")

# Check for overlapping union exons within genes (should be 0)
overlaps_check <- union_exons %>%
  group_by(gene_id) %>%
  arrange(union_exon_start) %>%
  mutate(
    next_start = lead(union_exon_start),
    overlaps_next = !is.na(next_start) & union_exon_end >= next_start
  ) %>%
  ungroup()

n_overlaps <- sum(overlaps_check$overlaps_next, na.rm = TRUE)

if (n_overlaps > 0) {
  cat(sprintf("  WARNING: %d overlapping union exons detected\n", n_overlaps))
} else {
  cat("  No overlapping union exons (correct)\n")
}

# Check mapping completeness
n_isoforms_with_mapping <- n_distinct(map_df$isoform_id)
n_isoforms_total <- nrow(isoform_structures)

cat(sprintf("  Mapping created for %d / %d isoforms (%.1f%%)\n",
            n_isoforms_with_mapping,
            n_isoforms_total,
            100 * n_isoforms_with_mapping / n_isoforms_total))

# Summary statistics
cat("\nAtomic union exon statistics:\n")
cat(sprintf("  Total union exons: %d\n", nrow(union_exons)))
ue_per_gene <- table(union_exons$gene_id)
cat(sprintf("  Union exons per gene - mean: %.1f, median: %.0f, range: %d-%d\n",
            mean(ue_per_gene), median(ue_per_gene),
            min(ue_per_gene), max(ue_per_gene)))
cat(sprintf("  Union exon length - mean: %.0f bp, median: %.0f bp\n",
            mean(union_exons$union_exon_length),
            median(union_exons$union_exon_length)))

cat("\nMapping statistics:\n")
cat(sprintf("  Total mappings: %d\n", nrow(map_df)))
cat(sprintf("  Union exons per isoform - mean: %.1f\n",
            nrow(map_df) / n_isoforms_with_mapping))

# Clean up temp files
file.remove(ue_tmp_file)
file.remove(map_tmp_file)

cat("\nStep 3 complete\n")
cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Genes processed: %d\n", genes_processed))
cat(sprintf("Atomic union exons created: %d\n", nrow(union_exons)))
cat(sprintf("Isoforms mapped: %d\n", n_isoforms_with_mapping))
cat(sprintf("Mapping entries: %d\n", nrow(map_df)))
cat(sprintf("Mean union exons per gene: %.1f\n", mean(ue_per_gene)))
cat(sprintf("Outputs: union_exons.rds, union_exons.tsv.gz (tabix), isoform_union_mapping.rds\n"))
cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")
