#!/usr/bin/env Rscript
# Gene Isoform Annotation Pipeline (Optimized)
# Extracts comprehensive isoform-level annotations for a target gene
# Integrates GENCODE, SQANTI, and expression data

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(readr)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Target gene
GENE_NAME <- "MTCL1"

# File paths
GENCODE_GTF <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.chrnamesedited.gtf"
SQANTI_CLASSIFICATION <- "/Users/petecastaldi/claude_projects/nmd/results/sqanti_runs/isoseq_sqanti3_filtered/sqanti3_classification.txt"
SQANTI_GTF <- "/Users/petecastaldi/claude_projects/nmd/results/sqanti_runs/isoseq_sqanti3_filtered/sqanti3_corrected.gtf"
SQANTI_GFF3 <- "/Users/petecastaldi/claude_projects/nmd/results/sqanti_runs/isoseq_sqanti3_filtered/sqanti3_corrected.cds.gff3"

# DGE files (all cell types)
DGE_DIR <- "/Users/petecastaldi/claude_projects/nmd/longread_dge"
DGE_FILES <- c(
  at2 = "nmd_dge_at2_2026.1.18.csv",
  dd = "nmd_dge_dd_2026.1.18.csv",
  ddali = "nmd_dge_ddali_2026.1.18.csv",
  doali = "nmd_dge_doali_2026.1.18.csv",
  fb = "nmd_dge_fb_2026.1.18.csv",
  mv = "nmd_dge_mv_2026.1.18.csv"
)

# Output
OUTPUT_FILE <- paste0("isoform_annotation_", GENE_NAME, ".tsv")

# =============================================================================
# FUNCTIONS
# =============================================================================

#' Extract exons for a transcript from GTF
get_exon_structure <- function(gtf, transcript_id) {
  tx_exons <- gtf[gtf$type == "exon" & gtf$transcript_id == transcript_id]

  if (length(tx_exons) == 0) {
    return(list(
      exons = NA,
      exon_types = NA,
      junctions = NA,
      n_exons = 0
    ))
  }

  # Sort by start position
  tx_exons <- tx_exons[order(start(tx_exons))]

  # Create union exon names: start_stop_strand
  strand_char <- as.character(strand(tx_exons)[1])

  # Determine exon type (CDS, UTR, etc) - default to "exon"
  # Note: GTF exon features don't typically have type annotations
  # We'll get CDS/UTR info from separate features
  exon_types <- rep("exon", length(tx_exons))
  exon_names <- paste0(start(tx_exons), "_", end(tx_exons), "_", strand_char)

  # Add type annotation to exon name
  exon_names_typed <- paste0(exon_names, ":", exon_types)

  # Create junction strings (between consecutive exons)
  if (length(tx_exons) > 1) {
    junctions <- paste0(
      end(tx_exons[-length(tx_exons)]),
      "_",
      start(tx_exons[-1]),
      "_",
      strand_char
    )
    junctions_str <- paste(junctions, collapse = ";")
  } else {
    junctions_str <- NA
  }

  return(list(
    exons = paste(exon_names_typed, collapse = ";"),
    exon_types = paste(exon_types, collapse = ";"),
    junctions = junctions_str,
    n_exons = length(tx_exons)
  ))
}

#' Get CDS information for a transcript with exon-level annotation
get_cds_info <- function(gtf, transcript_id) {
  tx_exons <- gtf[gtf$type == "exon" & gtf$transcript_id == transcript_id]
  cds <- gtf[gtf$type == "CDS" & gtf$transcript_id == transcript_id]
  utr5 <- gtf[gtf$type == "five_prime_utr" & gtf$transcript_id == transcript_id]
  utr3 <- gtf[gtf$type == "three_prime_utr" & gtf$transcript_id == transcript_id]

  if (length(cds) == 0) {
    return(list(
      has_cds = FALSE,
      cds_length = NA,
      protein_length = NA,
      exon_annotations = rep("exon", length(tx_exons))
    ))
  }

  # Annotate each exon as CDS, UTR5, UTR3, or mixed
  exon_annotations <- rep("exon", length(tx_exons))

  if (length(tx_exons) > 0) {
    for (i in seq_along(tx_exons)) {
      exon <- tx_exons[i]
      exon_range <- ranges(exon)

      # Check overlap with CDS
      cds_overlap <- any(overlapsAny(exon, cds))
      utr5_overlap <- any(overlapsAny(exon, utr5))
      utr3_overlap <- any(overlapsAny(exon, utr3))

      if (cds_overlap && (utr5_overlap || utr3_overlap)) {
        exon_annotations[i] <- "mixed"
      } else if (cds_overlap) {
        exon_annotations[i] <- "CDS"
      } else if (utr5_overlap) {
        exon_annotations[i] <- "5UTR"
      } else if (utr3_overlap) {
        exon_annotations[i] <- "3UTR"
      }
    }
  }

  cds_length <- sum(width(cds))
  protein_length <- floor(cds_length / 3)

  return(list(
    has_cds = TRUE,
    cds_length = cds_length,
    protein_length = protein_length,
    exon_annotations = exon_annotations
  ))
}

#' Get exon structure with type annotations
get_exon_structure_annotated <- function(gtf, transcript_id) {
  tx_exons <- gtf[gtf$type == "exon" & gtf$transcript_id == transcript_id]

  if (length(tx_exons) == 0) {
    return(list(
      exons = NA,
      junctions = NA,
      n_exons = 0
    ))
  }

  # Sort by start position
  tx_exons <- tx_exons[order(start(tx_exons))]

  # Get CDS info to annotate exons
  cds_info <- get_cds_info(gtf, transcript_id)

  # Create union exon names: start_stop_strand:type
  strand_char <- as.character(strand(tx_exons)[1])
  exon_names <- paste0(
    start(tx_exons), "_",
    end(tx_exons), "_",
    strand_char, ":",
    cds_info$exon_annotations
  )

  # Create junction strings
  if (length(tx_exons) > 1) {
    junctions <- paste0(
      end(tx_exons[-length(tx_exons)]),
      "_",
      start(tx_exons[-1]),
      "_",
      strand_char
    )
    junctions_str <- paste(junctions, collapse = ";")
  } else {
    junctions_str <- NA
  }

  return(list(
    exons = paste(exon_names, collapse = ";"),
    junctions = junctions_str,
    n_exons = length(tx_exons)
  ))
}

#' Get UTR lengths for a transcript
get_utr_lengths <- function(gtf, transcript_id) {
  utr5 <- gtf[gtf$type == "five_prime_utr" & gtf$transcript_id == transcript_id]
  utr3 <- gtf[gtf$type == "three_prime_utr" & gtf$transcript_id == transcript_id]

  return(list(
    utr5_length = ifelse(length(utr5) > 0, sum(width(utr5)), 0),
    utr3_length = ifelse(length(utr3) > 0, sum(width(utr3)), 0)
  ))
}

#' Get TSS and TES for a transcript
get_tss_tes <- function(gtf, transcript_id) {
  tx <- gtf[gtf$type == "transcript" & gtf$transcript_id == transcript_id]

  if (length(tx) == 0) {
    return(list(tss = NA, tes = NA))
  }

  if (as.character(strand(tx)) == "+") {
    tss <- start(tx)
    tes <- end(tx)
  } else {
    tss <- end(tx)
    tes <- start(tx)
  }

  return(list(tss = tss, tes = tes))
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

message(paste("Finding isoforms for gene:", GENE_NAME))

# =============================================================================
# Step 1: Find gene ID in GENCODE (using grep for speed)
# =============================================================================
message("Step 1: Finding gene ID in GENCODE...")

# Use grep to find gene lines matching MTCL1 (macOS compatible)
gene_lines <- system(
  paste0("grep 'gene_name \"", GENE_NAME, "\"' ", GENCODE_GTF, " | grep -w gene"),
  intern = TRUE
)

if (length(gene_lines) == 0) {
  stop(paste("Gene", GENE_NAME, "not found in GENCODE"))
}

# Extract gene_id from GTF attributes
gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', gene_lines[1])
message(paste("Found GENCODE gene ID:", gene_id))

# =============================================================================
# Step 2: Extract GENCODE isoforms (filtered)
# =============================================================================
message("Step 2: Extracting GENCODE isoforms...")

# Create temp filtered GTF for this gene only
temp_gtf <- tempfile(fileext = ".gtf")
system(paste0("grep '", gene_id, "' ", GENCODE_GTF, " > ", temp_gtf))

# Load only the filtered GTF
gencode <- import(temp_gtf)
unlink(temp_gtf)

message(paste("Loaded", length(gencode), "features from GENCODE"))

# Get all transcripts for this gene
gencode_transcripts <- gencode[gencode$type == "transcript" & gencode$gene_id == gene_id]
message(paste("Found", length(gencode_transcripts), "GENCODE transcripts"))

gencode_results <- data.frame()

for (i in seq_along(gencode_transcripts)) {
  tx <- gencode_transcripts[i]
  tx_id <- tx$transcript_id

  # Get exon structure with annotations
  exon_info <- get_exon_structure_annotated(gencode, tx_id)

  # Get CDS info
  cds_info <- get_cds_info(gencode, tx_id)

  # Get UTR lengths
  utr_info <- get_utr_lengths(gencode, tx_id)

  # Get TSS/TES
  tss_tes <- get_tss_tes(gencode, tx_id)

  # Determine protein coding status
  is_protein_coding <- FALSE
  coding_source <- NA

  if (!is.na(tx$transcript_type) && tx$transcript_type == "protein_coding") {
    is_protein_coding <- TRUE
    coding_source <- "GENCODE_biotype"
  } else if (cds_info$has_cds) {
    is_protein_coding <- TRUE
    coding_source <- "GENCODE_CDS"
  }

  gencode_results <- rbind(gencode_results, data.frame(
    isoform_id = tx_id,
    isoform_name = ifelse(!is.na(tx$transcript_name), tx$transcript_name, tx_id),
    source = "GENCODE",
    chromosome = as.character(seqnames(tx)),
    strand = as.character(strand(tx)),
    tss = tss_tes$tss,
    tes = tss_tes$tes,
    transcript_length = width(tx),
    n_exons = exon_info$n_exons,
    exons = exon_info$exons,
    junctions = exon_info$junctions,
    utr5_length = utr_info$utr5_length,
    utr3_length = utr_info$utr3_length,
    is_protein_coding = is_protein_coding,
    coding_source = coding_source,
    cds_length = cds_info$cds_length,
    protein_length = cds_info$protein_length,
    stringsAsFactors = FALSE
  ))
}

message(paste("Extracted", nrow(gencode_results), "GENCODE isoforms"))

# =============================================================================
# Step 3: Extract SQANTI/PacBio isoforms (filtered)
# =============================================================================
message("Step 3: Extracting SQANTI isoforms...")

# Pre-filter SQANTI classification for target gene using grep (much faster!)
message("  Pre-filtering SQANTI classification file...")
temp_sqanti_class <- tempfile(fileext = ".txt")
system(paste0("head -1 ", SQANTI_CLASSIFICATION, " > ", temp_sqanti_class))
system(paste0("grep '", gene_id, "' ", SQANTI_CLASSIFICATION, " >> ", temp_sqanti_class))

# Read filtered SQANTI classification
sqanti_class <- read_tsv(
  temp_sqanti_class,
  show_col_types = FALSE
)
unlink(temp_sqanti_class)

# Filter for target gene (should already be filtered, but just in case)
sqanti_target <- sqanti_class %>%
  filter(grepl(gene_id, associated_gene, ignore.case = TRUE))

message(paste("Found", nrow(sqanti_target), "SQANTI isoforms"))

if (nrow(sqanti_target) > 0) {
  # Extract PacBio isoform IDs
  pb_ids <- sqanti_target$isoform

  # Filter SQANTI GTF - write IDs to temp file for grep
  temp_ids <- tempfile()
  writeLines(pb_ids, temp_ids)

  temp_sqanti_gtf <- tempfile(fileext = ".gtf")
  system(paste0("grep -F -f ", temp_ids, " ", SQANTI_GTF, " > ", temp_sqanti_gtf))
  unlink(temp_ids)

  # Check if filtered file is not empty
  file_size <- file.info(temp_sqanti_gtf)$size

  if (is.na(file_size) || file_size == 0) {
    message("Warning: No SQANTI GTF features found for these isoforms")
    sqanti_gtf <- GRanges()
  } else {
    sqanti_gtf <- import(temp_sqanti_gtf)
    message(paste("Loaded", length(sqanti_gtf), "features from SQANTI GTF"))
  }
  unlink(temp_sqanti_gtf)

  # Try to load SQANTI GFF3 CDS info
  sqanti_gff3 <- tryCatch({
    temp_ids_gff3 <- tempfile()
    writeLines(pb_ids, temp_ids_gff3)

    temp_sqanti_gff3 <- tempfile(fileext = ".gff3")
    system(paste0("grep -F -f ", temp_ids_gff3, " ", SQANTI_GFF3, " > ", temp_sqanti_gff3))
    unlink(temp_ids_gff3)

    # Check if file has content
    file_size_gff3 <- file.info(temp_sqanti_gff3)$size

    if (is.na(file_size_gff3) || file_size_gff3 == 0) {
      message("Note: No SQANTI GFF3 features found")
      unlink(temp_sqanti_gff3)
      NULL
    } else {
      gff3_data <- import(temp_sqanti_gff3)
      unlink(temp_sqanti_gff3)
      message(paste("Loaded", length(gff3_data), "features from SQANTI GFF3"))
      gff3_data
    }
  }, error = function(e) {
    message(paste("Note: Could not load SQANTI GFF3:", e$message))
    NULL
  })

  sqanti_results <- data.frame()

  for (i in 1:nrow(sqanti_target)) {
    pb_id <- sqanti_target$isoform[i]

    # Get exon structure from SQANTI GTF (if available)
    if (length(sqanti_gtf) > 0) {
      exon_info <- get_exon_structure_annotated(sqanti_gtf, pb_id)
    } else {
      # If no GTF data, create basic exon info from SQANTI classification
      exon_info <- list(exons = NA, junctions = NA, n_exons = sqanti_target$exons[i])
    }

    # Get CDS info from SQANTI GFF3 if available
    cds_info <- list(has_cds = FALSE, cds_length = NA, protein_length = NA)
    if (!is.null(sqanti_gff3)) {
      cds_info <- get_cds_info(sqanti_gff3, pb_id)
    }

    # UTR info (usually not in SQANTI GTF)
    utr_info <- list(utr5_length = 0, utr3_length = 0)

    # Get TSS/TES (if GTF available)
    if (length(sqanti_gtf) > 0) {
      tss_tes <- get_tss_tes(sqanti_gtf, pb_id)
    } else {
      tss_tes <- list(tss = NA, tes = NA)
    }

    # Determine protein coding status
    is_protein_coding <- FALSE
    coding_source <- NA

    if (sqanti_target$coding[i] == "coding") {
      is_protein_coding <- TRUE
      coding_source <- "SQANTI"
    } else if (cds_info$has_cds) {
      is_protein_coding <- TRUE
      coding_source <- "SQANTI_GFF3_CDS"
    }

    # Use SQANTI classification info for protein length if available
    if (is_protein_coding && !is.na(sqanti_target$ORF_length[i])) {
      protein_length_sqanti <- sqanti_target$ORF_length[i]
    } else {
      protein_length_sqanti <- cds_info$protein_length
    }

    sqanti_results <- rbind(sqanti_results, data.frame(
      isoform_id = pb_id,
      isoform_name = pb_id,
      source = "SQANTI",
      chromosome = sqanti_target$chrom[i],
      strand = sqanti_target$strand[i],
      tss = tss_tes$tss,
      tes = tss_tes$tes,
      transcript_length = sqanti_target$length[i],
      n_exons = sqanti_target$exons[i],
      exons = exon_info$exons,
      junctions = exon_info$junctions,
      utr5_length = utr_info$utr5_length,
      utr3_length = utr_info$utr3_length,
      is_protein_coding = is_protein_coding,
      coding_source = coding_source,
      cds_length = ifelse(!is.na(sqanti_target$CDS_length[i]), sqanti_target$CDS_length[i], cds_info$cds_length),
      protein_length = protein_length_sqanti,
      stringsAsFactors = FALSE
    ))
  }
} else {
  sqanti_results <- data.frame()
}

# =============================================================================
# Step 4: Combine results
# =============================================================================
all_isoforms <- rbind(gencode_results, sqanti_results)
message(paste("Total isoforms:", nrow(all_isoforms)))

# =============================================================================
# Step 5: Add expression data
# =============================================================================
message("Step 4: Adding expression data from DMSO samples...")

for (cell_type in names(DGE_FILES)) {
  dge_file <- file.path(DGE_DIR, DGE_FILES[cell_type])

  if (!file.exists(dge_file)) {
    message(paste("Warning: DGE file not found:", dge_file))
    all_isoforms[[paste0("expr_DMSO_", cell_type)]] <- NA
    next
  }

  # Read only necessary columns
  dge <- read_csv(dge_file, show_col_types = FALSE, col_select = c(txid, AveExpr))

  # Match isoforms
  expr_col <- paste0("expr_DMSO_", cell_type)
  all_isoforms[[expr_col]] <- sapply(all_isoforms$isoform_id, function(tx_id) {
    match_row <- dge[dge$txid == tx_id, ]
    if (nrow(match_row) == 0) {
      return(NA)
    } else {
      return(match_row$AveExpr[1])
    }
  })

  message(paste("  ", cell_type, ":", sum(!is.na(all_isoforms[[expr_col]])), "isoforms with expression"))
}

# =============================================================================
# Step 6: Write output
# =============================================================================
message(paste("Writing output to:", OUTPUT_FILE))

write_tsv(all_isoforms, OUTPUT_FILE)

message("Done!")
message(paste("Total isoforms annotated:", nrow(all_isoforms)))
message(paste("  - GENCODE:", sum(all_isoforms$source == "GENCODE")))
message(paste("  - SQANTI:", sum(all_isoforms$source == "SQANTI")))
message(paste("  - Protein coding:", sum(all_isoforms$is_protein_coding, na.rm = TRUE)))
