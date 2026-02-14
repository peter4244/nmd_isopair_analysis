#!/usr/bin/env Rscript
# Gene Isoform Annotation Pipeline - Production Version
# Extracts comprehensive isoform annotations for genes of interest
#
# Features:
# - Combines GENCODE reference annotations with SQANTI PacBio isoforms
# - Detailed exon structures with CDS/UTR annotations
# - DMSO-specific expression per cell type (from DGEList)
# - Fast processing using grep pre-filtering and tabix indexing
# - Optional GTF and FASTA outputs
#
# Usage:
#   Rscript gene_isoform_annotation.R <GENE_NAME_OR_ID>
#   Example: Rscript gene_isoform_annotation.R MTCL1
#   Example: Rscript gene_isoform_annotation.R ENSG00000168502
#
# Prerequisites:
#   1. Run preprocess_sqanti_tabix.sh once to create tabix-indexed SQANTI GTF
#   2. Ensure all reference files and DGEList RDS are available

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(edgeR)
})

# =============================================================================
# COMMAND-LINE ARGUMENTS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("Usage: Rscript gene_isoform_annotation.R <GENE_NAME_OR_ID> [OPTIONS]\n")
  cat("\nArguments:\n")
  cat("  GENE_NAME_OR_ID    Gene name (e.g., MTCL1) or Ensembl ID (e.g., ENSG00000168502)\n")
  cat("\nOptions:\n")
  cat("  --gtf              Output GTF file with all isoforms\n")
  cat("  --fasta            Output FASTA file with isoform sequences\n")
  cat("  --no-expr          Skip expression calculation\n")
  cat("\nExamples:\n")
  cat("  Rscript gene_isoform_annotation.R MTCL1\n")
  cat("  Rscript gene_isoform_annotation.R MTCL1 --gtf --fasta\n")
  cat("  Rscript gene_isoform_annotation.R ENSG00000168502 --no-expr\n")
  stop("No gene name or ID provided")
}

GENE_INPUT <- args[1]

# Parse optional flags
OUTPUT_GTF <- "--gtf" %in% args
OUTPUT_FASTA <- "--fasta" %in% args
INCLUDE_EXPRESSION <- !"--no-expr" %in% args

# Determine if input is Ensembl ID or gene name
if (grepl("^ENSG[0-9]+", GENE_INPUT)) {
  GENE_ID <- sub("\\..*", "", GENE_INPUT)  # Remove version if present
  GENE_NAME <- NULL  # Will be determined from GTF
  SEARCH_BY <- "id"
} else {
  GENE_NAME <- GENE_INPUT
  GENE_ID <- NULL  # Will be determined from GTF
  SEARCH_BY <- "name"
}

# =============================================================================
# CONFIGURATION
# =============================================================================

# Note: OUTPUT_GTF, OUTPUT_FASTA, and INCLUDE_EXPRESSION are now set via command-line flags

# File paths (project standard paths)
# Note: Both GTF files must be tabix-indexed (run preprocess_sqanti_tabix.sh first)
GENCODE_GTF_UNINDEXED <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.chrnamesedited.gtf"
GENCODE_GTF_INDEXED <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.chrnamesedited.sorted.gtf.gz"
GENCODE_FASTA <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.transcripts.fa.gz"
SQANTI_CLASSIFICATION <- "/Users/petecastaldi/claude_projects/nmd/results/sqanti_runs/isoseq_sqanti3_filtered/sqanti3_classification.txt"
SQANTI_GTF_INDEXED <- "/Users/petecastaldi/claude_projects/nmd/results/sqanti_runs/isoseq_sqanti3_filtered/sqanti3_corrected.sorted.gtf.gz"
DGELIST_RDS <- "/Users/petecastaldi/claude_projects/nmd/results/rds/dge_isoform_2026.1.20.rds"
SQANTI_FASTA <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode49_merged_collapsed_2025.12.21.fa.gz"

# Output files (will be set after gene name is resolved)
OUTPUT_TSV <- NULL
OUTPUT_GTF_FILE <- NULL
OUTPUT_FASTA_FILE <- NULL

# Cell types for expression data
CELL_TYPES <- c("DD_ALI", "DD", "DO_ALI", "AT2", "FB", "MV")

# =============================================================================
# HELPER FUNCTIONS FOR GENCODE PARSING
# =============================================================================

get_exon_structure_annotated <- function(gtf, transcript_id) {
  # Extract exon positions with CDS/UTR annotations
  tx_exons <- gtf[gtf$type == "exon" & gtf$transcript_id == transcript_id]
  if (length(tx_exons) == 0) return(list(exons = NA, junctions = NA, n_exons = 0))

  tx_exons <- tx_exons[order(start(tx_exons))]
  strand_char <- as.character(strand(tx_exons)[1])

  # Get CDS/UTR annotations
  cds <- gtf[gtf$type == "CDS" & gtf$transcript_id == transcript_id]
  utr5 <- gtf[gtf$type == "five_prime_utr" & gtf$transcript_id == transcript_id]
  utr3 <- gtf[gtf$type == "three_prime_utr" & gtf$transcript_id == transcript_id]

  # Annotate each exon
  exon_types <- rep("exon", length(tx_exons))
  for (i in seq_along(tx_exons)) {
    exon <- tx_exons[i]
    cds_overlap <- any(overlapsAny(exon, cds))
    utr5_overlap <- any(overlapsAny(exon, utr5))
    utr3_overlap <- any(overlapsAny(exon, utr3))

    if (cds_overlap && (utr5_overlap || utr3_overlap)) {
      exon_types[i] <- "mixed"
    } else if (cds_overlap) {
      exon_types[i] <- "CDS"
    } else if (utr5_overlap) {
      exon_types[i] <- "5UTR"
    } else if (utr3_overlap) {
      exon_types[i] <- "3UTR"
    }
  }

  exon_names <- paste0(start(tx_exons), "_", end(tx_exons), "_", strand_char, ":", exon_types)

  # Junctions
  if (length(tx_exons) > 1) {
    junctions <- paste0(end(tx_exons[-length(tx_exons)]), "_", start(tx_exons[-1]), "_", strand_char)
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

get_exon_structure_simple <- function(gtf, transcript_id) {
  # Extract exon positions without CDS/UTR annotations (for SQANTI)
  tx_exons <- gtf[gtf$type == "exon" & gtf$transcript_id == transcript_id]
  if (length(tx_exons) == 0) return(list(exons = NA, junctions = NA, n_exons = 0))

  tx_exons <- tx_exons[order(start(tx_exons))]
  strand_char <- as.character(strand(tx_exons)[1])

  # Simple exon annotation (no CDS/UTR info for novel isoforms)
  exon_names <- paste0(start(tx_exons), "_", end(tx_exons), "_", strand_char, ":exon")

  # Junctions
  if (length(tx_exons) > 1) {
    junctions <- paste0(end(tx_exons[-length(tx_exons)]), "_", start(tx_exons[-1]), "_", strand_char)
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

get_cds_info <- function(gtf, transcript_id) {
  # Extract CDS length, protein length, and CDS genomic coordinates
  cds <- gtf[gtf$type == "CDS" & gtf$transcript_id == transcript_id]
  if (length(cds) == 0) {
    return(list(
      has_cds = FALSE,
      cds_length = NA,
      protein_length = NA,
      cds_start = NA,
      cds_end = NA
    ))
  }
  cds_length <- sum(width(cds))

  # CDS start and end are genomic coordinates (always start < end regardless of strand)
  cds_start <- min(start(cds))
  cds_end <- max(end(cds))

  return(list(
    has_cds = TRUE,
    cds_length = cds_length,
    protein_length = floor(cds_length / 3),
    cds_start = cds_start,
    cds_end = cds_end
  ))
}

calculate_utr_lengths <- function(tss, tes, cds_start, cds_end, strand) {
  # Calculate UTR lengths from TSS, TES, and CDS positions (strand-aware)
  # Returns NA if no CDS information available

  if (is.na(cds_start) || is.na(cds_end) || is.na(tss) || is.na(tes)) {
    return(list(utr5_length = NA, utr3_length = NA))
  }

  if (strand == "+") {
    # Plus strand: TSS < CDS_start < CDS_end < TES
    utr5_length <- max(0, cds_start - tss)
    utr3_length <- max(0, tes - cds_end)
  } else if (strand == "-") {
    # Minus strand: TES < CDS_start < CDS_end < TSS
    utr5_length <- max(0, tss - cds_end)
    utr3_length <- max(0, cds_start - tes)
  } else {
    return(list(utr5_length = NA, utr3_length = NA))
  }

  return(list(
    utr5_length = utr5_length,
    utr3_length = utr3_length
  ))
}

get_tss_tes <- function(gtf, transcript_id) {
  # Get transcription start and end sites
  tx <- gtf[gtf$type == "transcript" & gtf$transcript_id == transcript_id]
  if (length(tx) == 0) return(list(tss = NA, tes = NA))

  if (as.character(strand(tx)) == "+") {
    return(list(tss = start(tx), tes = end(tx)))
  } else {
    return(list(tss = end(tx), tes = start(tx)))
  }
}

# =============================================================================
# VALIDATION CHECKS
# =============================================================================

message("====================================================================")
message("Gene Isoform Annotation Pipeline")
message("====================================================================")
message(paste("Gene:", GENE_NAME))
message(paste("Date:", Sys.Date()))
message("")

# Check required files
message("Checking required files...")
files_to_check <- c(
  "GENCODE GTF (unindexed)" = GENCODE_GTF_UNINDEXED,
  "GENCODE GTF (indexed)" = GENCODE_GTF_INDEXED,
  "GENCODE GTF index" = paste0(GENCODE_GTF_INDEXED, ".tbi"),
  "SQANTI GTF (indexed)" = SQANTI_GTF_INDEXED,
  "SQANTI GTF index" = paste0(SQANTI_GTF_INDEXED, ".tbi"),
  "SQANTI Classification" = SQANTI_CLASSIFICATION
)

if (INCLUDE_EXPRESSION) {
  files_to_check["DGEList RDS"] <- DGELIST_RDS
}

if (OUTPUT_FASTA) {
  files_to_check["GENCODE FASTA"] <- GENCODE_FASTA
  files_to_check["SQANTI FASTA"] <- SQANTI_FASTA
}

missing_files <- c()
for (file_desc in names(files_to_check)) {
  file_path <- files_to_check[file_desc]
  if (file.exists(file_path)) {
    message(paste("  [OK]", file_desc))
  } else {
    message(paste("  [MISSING]", file_desc, "-", file_path))
    missing_files <- c(missing_files, file_path)
  }
}

if (length(missing_files) > 0) {
  stop(paste("ERROR: Missing required files. Please run preprocess_sqanti_tabix.sh first to create indexed GTF files."))
}

# Check for tabix command
if (Sys.which("tabix") == "") {
  stop("ERROR: tabix command not found. Install with: conda install -c bioconda htslib")
}

# Check for Biostrings if FASTA output requested
if (OUTPUT_FASTA) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("ERROR: Biostrings package required for FASTA output.\nInstall with: BiocManager::install('Biostrings')")
  }
}

message("")
message("All checks passed!")
message("")

# =============================================================================
# STEP 1: FIND GENE ID
# =============================================================================

message("====================================================================")
message("STEP 1: Finding gene in GENCODE")
message("====================================================================")

if (SEARCH_BY == "name") {
  message(paste("  Searching by gene name:", GENE_NAME))
  gene_lines <- system(
    paste0("grep 'gene_name \"", GENE_NAME, "\"' ", GENCODE_GTF_UNINDEXED, " | grep -w gene"),
    intern = TRUE
  )

  if (length(gene_lines) == 0) {
    stop(paste("ERROR: Gene name", GENE_NAME, "not found in GENCODE GTF"))
  }

  gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', gene_lines[1])
  gene_name_found <- sub('.*gene_name "([^"]+)".*', '\\1', gene_lines[1])

  message(paste("  Found gene ID:", gene_id))
  message(paste("  Gene name:", gene_name_found))

} else {
  message(paste("  Searching by Ensembl ID:", GENE_ID))
  gene_lines <- system(
    paste0("grep 'gene_id \"", GENE_ID, "' ", GENCODE_GTF_UNINDEXED, " | grep -w gene | head -1"),
    intern = TRUE
  )

  if (length(gene_lines) == 0) {
    stop(paste("ERROR: Gene ID", GENE_ID, "not found in GENCODE GTF"))
  }

  gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', gene_lines[1])
  gene_name_found <- sub('.*gene_name "([^"]+)".*', '\\1', gene_lines[1])

  GENE_NAME <- gene_name_found  # Set for output file naming

  message(paste("  Found gene ID:", gene_id))
  message(paste("  Gene name:", gene_name_found))
}

# Parse gene coordinates from the gene line for tabix query
gene_parts <- strsplit(gene_lines[1], "\t")[[1]]
gene_chr <- gene_parts[1]
gene_start <- as.numeric(gene_parts[4])
gene_end <- as.numeric(gene_parts[5])

message(paste("  Gene coordinates:", paste0(gene_chr, ":", gene_start, "-", gene_end)))

# Set output filenames now that GENE_NAME is determined
OUTPUT_TSV <- paste0("isoform_annotation_", GENE_NAME, ".tsv")
OUTPUT_GTF_FILE <- paste0("isoform_annotation_", GENE_NAME, ".gtf")
OUTPUT_FASTA_FILE <- paste0("isoform_annotation_", GENE_NAME, ".fasta")

message("")

# =============================================================================
# STEP 2: EXTRACT GENCODE ISOFORMS
# =============================================================================

message("====================================================================")
message("STEP 2: Extracting GENCODE isoforms")
message("====================================================================")

# Query tabix-indexed GENCODE GTF by coordinates (fast: <1 second)
message("  Querying tabix-indexed GENCODE GTF...")
tabix_cmd <- sprintf("tabix %s %s:%d-%d", GENCODE_GTF_INDEXED, gene_chr, gene_start, gene_end)
gencode_lines <- system(tabix_cmd, intern = TRUE)

if (length(gencode_lines) == 0) {
  stop(paste("ERROR: No features found for gene coordinates:", paste0(gene_chr, ":", gene_start, "-", gene_end)))
}

# Filter for this specific gene (tabix returns all features in region, may include overlapping genes)
gencode_lines <- gencode_lines[grepl(paste0('gene_id "', gene_id), gencode_lines)]

message(paste("  Retrieved", length(gencode_lines), "features from indexed GTF"))

# Parse into GRanges
temp_gtf <- tempfile(fileext = ".gtf")
writeLines(gencode_lines, temp_gtf)
gencode <- import(temp_gtf)
unlink(temp_gtf)

gencode_transcripts <- gencode[gencode$type == "transcript" & gencode$gene_id == gene_id]
message(paste("  Found", length(gencode_transcripts), "GENCODE transcripts"))

gencode_results <- data.frame()

message("  Extracting detailed annotations...")
for (i in seq_along(gencode_transcripts)) {
  tx <- gencode_transcripts[i]
  tx_id <- tx$transcript_id

  exon_info <- get_exon_structure_annotated(gencode, tx_id)
  cds_info <- get_cds_info(gencode, tx_id)
  tss_tes <- get_tss_tes(gencode, tx_id)

  # Calculate UTR lengths from TSS/TES/CDS positions (strand-aware)
  utr_info <- calculate_utr_lengths(
    tss_tes$tss,
    tss_tes$tes,
    cds_info$cds_start,
    cds_info$cds_end,
    as.character(strand(tx))
  )

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
    cds_start = cds_info$cds_start,
    cds_end = cds_info$cds_end,
    cds_length = cds_info$cds_length,
    protein_length = cds_info$protein_length,
    utr5_length = utr_info$utr5_length,
    utr3_length = utr_info$utr3_length,
    is_protein_coding = is_protein_coding,
    coding_source = coding_source,
    stringsAsFactors = FALSE
  ))
}

message(paste("  Extracted", nrow(gencode_results), "GENCODE isoforms with full annotations"))
message("")

# =============================================================================
# STEP 3: EXTRACT SQANTI ISOFORMS
# =============================================================================

message("====================================================================")
message("STEP 3: Extracting SQANTI isoforms")
message("====================================================================")

# Part A: Get isoform metadata from SQANTI classification
message("  Loading SQANTI classification...")
temp_sqanti <- tempfile(fileext = ".txt")
system(paste0("head -1 ", SQANTI_CLASSIFICATION, " > ", temp_sqanti))
system(paste0("grep '", gene_id, "' ", SQANTI_CLASSIFICATION, " >> ", temp_sqanti))

sqanti_class <- read_tsv(temp_sqanti, show_col_types = FALSE)
unlink(temp_sqanti)

sqanti_target <- sqanti_class %>%
  filter(grepl(paste0("\\b", gene_id, "\\b"), associated_gene, ignore.case = TRUE))

message(paste("  Found", nrow(sqanti_target), "SQANTI isoforms"))

sqanti_results <- data.frame()

if (nrow(sqanti_target) > 0) {
  # Part B: Get exon structures from tabix-indexed SQANTI GTF
  message("  Querying tabix-indexed SQANTI GTF for exon structures...")

  # Query tabix by gene coordinates
  tabix_cmd <- sprintf("tabix %s %s:%d-%d", SQANTI_GTF_INDEXED, gene_chr, gene_start, gene_end)
  sqanti_lines <- system(tabix_cmd, intern = TRUE)

  sqanti_features <- NULL
  if (length(sqanti_lines) > 0) {
    message(paste("    Retrieved", length(sqanti_lines), "GTF lines"))

    # Parse tabix output into GRanges
    temp_gtf <- tempfile(fileext = ".gtf")
    writeLines(sqanti_lines, temp_gtf)
    sqanti_features <- import(temp_gtf)
    unlink(temp_gtf)

    message("    Parsed SQANTI GTF features")
  } else {
    message("    No SQANTI features found in gene region")
  }

  # Process each SQANTI isoform
  message("  Extracting SQANTI isoform annotations...")
  for (i in 1:nrow(sqanti_target)) {
    row <- sqanti_target[i, ]
    pb_id <- row$isoform

    # Get exon structure and CDS info from GTF if available
    exon_info <- list(exons = NA, junctions = NA, n_exons = row$exons)
    tss_tes <- list(tss = NA, tes = NA)
    cds_info <- list(has_cds = FALSE, cds_start = NA, cds_end = NA, cds_length = NA, protein_length = NA)

    if (!is.null(sqanti_features)) {
      # Check if this isoform is in the GTF
      tx_features <- sqanti_features[sqanti_features$transcript_id == pb_id]

      if (length(tx_features) > 0) {
        exon_info <- get_exon_structure_simple(sqanti_features, pb_id)
        tss_tes <- get_tss_tes(sqanti_features, pb_id)
        cds_info <- get_cds_info(sqanti_features, pb_id)
      }
    }

    # Use SQANTI classification for CDS/protein if GTF doesn't have it
    if (!cds_info$has_cds) {
      cds_info$cds_length <- ifelse(!is.na(row$CDS_length), row$CDS_length, NA)
      cds_info$protein_length <- ifelse(!is.na(row$ORF_length), row$ORF_length, NA)
    }

    # Calculate UTR lengths from TSS/TES/CDS if available
    utr_info <- calculate_utr_lengths(
      tss_tes$tss,
      tss_tes$tes,
      cds_info$cds_start,
      cds_info$cds_end,
      row$strand
    )

    # Determine coding status
    is_protein_coding <- FALSE
    coding_source <- NA
    if (!is.na(row$coding) && row$coding == "coding") {
      is_protein_coding <- TRUE
      coding_source <- "SQANTI"
    }

    sqanti_results <- rbind(sqanti_results, data.frame(
      isoform_id = pb_id,
      isoform_name = pb_id,
      source = "SQANTI",
      chromosome = row$chrom,
      strand = row$strand,
      tss = tss_tes$tss,
      tes = tss_tes$tes,
      transcript_length = row$length,
      n_exons = exon_info$n_exons,
      exons = exon_info$exons,
      junctions = exon_info$junctions,
      cds_start = cds_info$cds_start,
      cds_end = cds_info$cds_end,
      cds_length = cds_info$cds_length,
      protein_length = cds_info$protein_length,
      utr5_length = utr_info$utr5_length,
      utr3_length = utr_info$utr3_length,
      is_protein_coding = is_protein_coding,
      coding_source = coding_source,
      stringsAsFactors = FALSE
    ))
  }

  message(paste("  Extracted", nrow(sqanti_results), "SQANTI isoforms"))
}

message("")

# =============================================================================
# STEP 4: CALCULATE DMSO-SPECIFIC EXPRESSION (OPTIONAL)
# =============================================================================

# Combine GENCODE and SQANTI results
all_isoforms <- rbind(gencode_results, sqanti_results)

if (INCLUDE_EXPRESSION) {
  message("====================================================================")
  message("STEP 4: Calculating DMSO-specific expression per cell type")
  message("====================================================================")
  message(paste("  Total isoforms:", nrow(all_isoforms)))

# Load DGEList
message("  Loading DGEList...")
dge_isoform <- readRDS(DGELIST_RDS)

# Get normalized CPM
message("  Calculating normalized CPM...")
cpm_normalized <- cpm(dge_isoform, normalized = TRUE, log = FALSE)

# Get transcript IDs from genes component
transcript_ids <- dge_isoform$genes$txid
rownames(cpm_normalized) <- transcript_ids

# Get sample metadata
sample_info <- dge_isoform$samples

# Validate required columns
if (!"bamid" %in% colnames(sample_info)) {
  stop("ERROR: 'bamid' column not found in DGEList$samples. Available columns: ",
       paste(colnames(sample_info), collapse = ", "))
}
if (!"treatment" %in% colnames(sample_info)) {
  stop("ERROR: 'treatment' column not found in DGEList$samples")
}
if (!"ct" %in% colnames(sample_info)) {
  stop("ERROR: 'ct' column not found in DGEList$samples")
}

# Map cell type names to DGEList codes
ct_map <- c(
  "DD_ALI" = "DD_ALI",
  "DD" = "DD",
  "DO_ALI" = "DO_ALI",
  "AT2" = "AT",
  "FB" = "FB",
  "MV" = "MV"
)

# Calculate mean DMSO expression per cell type
message("  Calculating DMSO means per cell type...")
for (ct_name in names(ct_map)) {
  ct_code <- ct_map[ct_name]

  # Filter for DMSO samples of this cell type
  dmso_samples <- sample_info %>%
    filter(treatment == "DMSO", ct == ct_code) %>%
    pull(bamid)

  expr_col_name <- paste0("expr_DMSO_", tolower(gsub("_", "", ct_name)))

  if (length(dmso_samples) > 0) {
    # Calculate mean across DMSO samples
    dmso_means <- rowMeans(cpm_normalized[, dmso_samples, drop = FALSE])

    # Match to isoforms by transcript ID
    matched <- all_isoforms$isoform_id %in% names(dmso_means)
    n_missing <- sum(!matched)
    all_isoforms[[expr_col_name]] <- dmso_means[all_isoforms$isoform_id]

    n_expressed <- sum(!is.na(all_isoforms[[expr_col_name]]) & all_isoforms[[expr_col_name]] > 0)
    message(paste("    ", ct_name, ":", n_expressed, "isoforms with DMSO expression > 0"))
    if (n_missing > 0) {
      message(paste("      Note:", n_missing, "isoforms not found in DGEList"))
    }
  } else {
    all_isoforms[[expr_col_name]] <- NA
    message(paste("    ", ct_name, ": No DMSO samples found"))
  }
}

  message("")
} else {
  message("====================================================================")
  message("STEP 4: Skipping expression calculation (INCLUDE_EXPRESSION = FALSE)")
  message("====================================================================")
  message(paste("  Total isoforms:", nrow(all_isoforms)))
  message("")
}

# =============================================================================
# STEP 5: WRITE OUTPUT FILES
# =============================================================================

message("====================================================================")
message("STEP 5: Writing output files")
message("====================================================================")

# Round expression values to 2 significant digits (if included)
if (INCLUDE_EXPRESSION) {
  expr_cols <- grep("^expr_DMSO_", colnames(all_isoforms), value = TRUE)
  for (col in expr_cols) {
    all_isoforms[[col]] <- signif(all_isoforms[[col]], 2)
  }
}

# Write TSV
message(paste("  Writing TSV:", OUTPUT_TSV))
write_tsv(all_isoforms, OUTPUT_TSV)
message(paste("    Wrote", nrow(all_isoforms), "isoforms"))

# Optional: Write GTF
if (OUTPUT_GTF) {
  message(paste("  Writing GTF:", OUTPUT_GTF_FILE))

  # Combine GENCODE and SQANTI GTF features
  output_gtf <- gencode

  if (!is.null(sqanti_features) && length(sqanti_features) > 0) {
    output_gtf <- c(output_gtf, sqanti_features)
  }

  export(output_gtf, OUTPUT_GTF_FILE, format = "gtf")
  message(paste("    Wrote", length(output_gtf), "GTF features"))
}

# Optional: Write FASTA
if (OUTPUT_FASTA) {
  message(paste("  Writing FASTA:", OUTPUT_FASTA_FILE))
  message("    (This may take 1-2 minutes for large FASTA files)")

  library(Biostrings)

  # Separate GENCODE and SQANTI isoforms
  gencode_ids <- all_isoforms$isoform_id[all_isoforms$source == "GENCODE"]
  sqanti_ids <- all_isoforms$isoform_id[all_isoforms$source == "SQANTI"]

  # Read GENCODE FASTA and extract GENCODE isoforms
  message("    Reading GENCODE transcripts FASTA...")
  gencode_fasta <- readDNAStringSet(GENCODE_FASTA)

  # GENCODE FASTA has pipe-delimited headers like: ENST00000832824.1|ENSG...|...
  # Extract transcript ID (first field before pipe)
  gencode_fasta_ids <- sub("\\|.*", "", names(gencode_fasta))

  # Match by transcript ID
  matches_idx <- gencode_fasta_ids %in% gencode_ids
  gencode_matches <- gencode_fasta[matches_idx]

  # Simplify names to just transcript ID
  names(gencode_matches) <- gencode_fasta_ids[matches_idx]

  message(paste("      Found", length(gencode_matches), "GENCODE sequences"))

  # Read SQANTI FASTA and extract SQANTI isoforms
  message("    Reading SQANTI transcripts FASTA...")
  sqanti_fasta <- readDNAStringSet(SQANTI_FASTA)
  sqanti_matches <- sqanti_fasta[names(sqanti_fasta) %in% sqanti_ids]
  message(paste("      Found", length(sqanti_matches), "SQANTI sequences"))

  # Combine sequences
  all_sequences <- c(gencode_matches, sqanti_matches)

  if (length(all_sequences) > 0) {
    writeXStringSet(all_sequences, OUTPUT_FASTA_FILE, format = "fasta")
    message(paste("    Wrote", length(all_sequences), "total sequences to", OUTPUT_FASTA_FILE))
  } else {
    message("    WARNING: No matching sequences found in FASTA files")
  }
}

message("")

# =============================================================================
# SUMMARY
# =============================================================================

message("====================================================================")
message("SUMMARY")
message("====================================================================")
message(paste("Gene:", GENE_NAME, "(", gene_id, ")"))
message("")
message("Isoform counts:")
message(paste("  GENCODE:", sum(all_isoforms$source == "GENCODE")))
message(paste("  SQANTI:", sum(all_isoforms$source == "SQANTI")))
message(paste("  Total:", nrow(all_isoforms)))
message("")
message(paste("  Protein coding:", sum(all_isoforms$is_protein_coding, na.rm = TRUE)))
message(paste("  Non-coding:", sum(!all_isoforms$is_protein_coding, na.rm = TRUE)))
message("")
message("Output files:")
message(paste("  ", OUTPUT_TSV))
if (OUTPUT_GTF) message(paste("  ", OUTPUT_GTF_FILE))
if (OUTPUT_FASTA) message(paste("  ", OUTPUT_FASTA_FILE))
message("")
message("====================================================================")
message("DONE!")
message("====================================================================")
