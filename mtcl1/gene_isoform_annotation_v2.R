#!/usr/bin/env Rscript
# Gene Isoform Annotation Pipeline v2 (Simplified & Fast)
# Uses GENCODE GTF for detailed parsing, SQANTI classification for PacBio isoforms

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(readr)

# =============================================================================
# CONFIGURATION
# =============================================================================

GENE_NAME <- "MTCL1"

# File paths
GENCODE_GTF <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.chrnamesedited.gtf"
SQANTI_CLASSIFICATION <- "/Users/petecastaldi/claude_projects/nmd/results/sqanti_runs/isoseq_sqanti3_filtered/sqanti3_classification.txt"

# DGE files
DGE_DIR <- "/Users/petecastaldi/claude_projects/nmd/longread_dge"
DGE_FILES <- c(
  at2 = "nmd_dge_at2_2026.1.18.csv",
  dd = "nmd_dge_dd_2026.1.18.csv",
  ddali = "nmd_dge_ddali_2026.1.18.csv",
  doali = "nmd_dge_doali_2026.1.18.csv",
  fb = "nmd_dge_fb_2026.1.18.csv",
  mv = "nmd_dge_mv_2026.1.18.csv"
)

OUTPUT_FILE <- paste0("isoform_annotation_", GENE_NAME, ".tsv")

# =============================================================================
# FUNCTIONS FOR GENCODE PARSING
# =============================================================================

get_exon_structure_annotated <- function(gtf, transcript_id) {
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

get_cds_info <- function(gtf, transcript_id) {
  cds <- gtf[gtf$type == "CDS" & gtf$transcript_id == transcript_id]
  if (length(cds) == 0) {
    return(list(has_cds = FALSE, cds_length = NA, protein_length = NA))
  }
  cds_length <- sum(width(cds))
  return(list(has_cds = TRUE, cds_length = cds_length, protein_length = floor(cds_length / 3)))
}

get_utr_lengths <- function(gtf, transcript_id) {
  utr5 <- gtf[gtf$type == "five_prime_utr" & gtf$transcript_id == transcript_id]
  utr3 <- gtf[gtf$type == "three_prime_utr" & gtf$transcript_id == transcript_id]
  return(list(
    utr5_length = ifelse(length(utr5) > 0, sum(width(utr5)), 0),
    utr3_length = ifelse(length(utr3) > 0, sum(width(utr3)), 0)
  ))
}

get_tss_tes <- function(gtf, transcript_id) {
  tx <- gtf[gtf$type == "transcript" & gtf$transcript_id == transcript_id]
  if (length(tx) == 0) return(list(tss = NA, tes = NA))

  if (as.character(strand(tx)) == "+") {
    return(list(tss = start(tx), tes = end(tx)))
  } else {
    return(list(tss = end(tx), tes = start(tx)))
  }
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

message(paste("=== Analyzing gene:", GENE_NAME, "==="))

# Find gene ID
message("Step 1: Finding gene ID in GENCODE...")
gene_lines <- system(
  paste0("grep 'gene_name \"", GENE_NAME, "\"' ", GENCODE_GTF, " | grep -w gene"),
  intern = TRUE
)
if (length(gene_lines) == 0) stop(paste("Gene", GENE_NAME, "not found"))

gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', gene_lines[1])
message(paste("  Found:", gene_id))

# Extract GENCODE isoforms
message("Step 2: Extracting GENCODE isoforms...")
temp_gtf <- tempfile(fileext = ".gtf")
system(paste0("grep '", gene_id, "' ", GENCODE_GTF, " > ", temp_gtf))
gencode <- import(temp_gtf)
unlink(temp_gtf)

gencode_transcripts <- gencode[gencode$type == "transcript" & gencode$gene_id == gene_id]
message(paste("  Found", length(gencode_transcripts), "GENCODE transcripts"))

gencode_results <- data.frame()

for (i in seq_along(gencode_transcripts)) {
  tx <- gencode_transcripts[i]
  tx_id <- tx$transcript_id

  exon_info <- get_exon_structure_annotated(gencode, tx_id)
  cds_info <- get_cds_info(gencode, tx_id)
  utr_info <- get_utr_lengths(gencode, tx_id)
  tss_tes <- get_tss_tes(gencode, tx_id)

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

message(paste("  Extracted", nrow(gencode_results), "GENCODE isoforms"))

# Extract SQANTI isoforms (using classification file only - much faster!)
message("Step 3: Extracting SQANTI isoforms...")
temp_sqanti <- tempfile(fileext = ".txt")
system(paste0("head -1 ", SQANTI_CLASSIFICATION, " > ", temp_sqanti))
system(paste0("grep '", gene_id, "' ", SQANTI_CLASSIFICATION, " >> ", temp_sqanti))

sqanti_class <- read_tsv(temp_sqanti, show_col_types = FALSE)
unlink(temp_sqanti)

sqanti_target <- sqanti_class %>%
  filter(grepl(gene_id, associated_gene, ignore.case = TRUE))

message(paste("  Found", nrow(sqanti_target), "SQANTI isoforms"))

if (nrow(sqanti_target) > 0) {
  sqanti_results <- sqanti_target %>%
    mutate(
      isoform_id = isoform,
      isoform_name = isoform,
      source = "SQANTI",
      chromosome = chrom,
      tss = NA,  # Not easily available without GTF parsing
      tes = NA,
      transcript_length = length,
      n_exons = exons,
      exons = NA,  # Would require GTF parsing (too slow)
      junctions = NA,
      utr5_length = NA,
      utr3_length = NA,
      is_protein_coding = (coding == "coding"),
      coding_source = ifelse(coding == "coding", "SQANTI", NA),
      cds_length = CDS_length,
      protein_length = ORF_length
    ) %>%
    select(isoform_id, isoform_name, source, chromosome, strand, tss, tes,
           transcript_length, n_exons, exons, junctions, utr5_length, utr3_length,
           is_protein_coding, coding_source, cds_length, protein_length)
} else {
  sqanti_results <- data.frame()
}

# Combine
all_isoforms <- rbind(gencode_results, sqanti_results)
message(paste("Total isoforms:", nrow(all_isoforms)))

# Add expression data
message("Step 4: Adding expression data...")

for (cell_type in names(DGE_FILES)) {
  dge_file <- file.path(DGE_DIR, DGE_FILES[cell_type])

  if (!file.exists(dge_file)) {
    all_isoforms[[paste0("expr_DMSO_", cell_type)]] <- NA
    next
  }

  dge <- read_csv(dge_file, show_col_types = FALSE, col_select = c(txid, AveExpr))

  expr_col <- paste0("expr_DMSO_", cell_type)
  all_isoforms[[expr_col]] <- sapply(all_isoforms$isoform_id, function(tx_id) {
    match_row <- dge[dge$txid == tx_id, ]
    ifelse(nrow(match_row) == 0, NA, match_row$AveExpr[1])
  })

  message(paste("  ", cell_type, ":", sum(!is.na(all_isoforms[[expr_col]])), "with expression"))
}

# Write output
message(paste("Step 5: Writing output to", OUTPUT_FILE))
write_tsv(all_isoforms, OUTPUT_FILE)

message("\n=== SUMMARY ===")
message(paste("Total isoforms:", nrow(all_isoforms)))
message(paste("  GENCODE:", sum(all_isoforms$source == "GENCODE")))
message(paste("  SQANTI:", sum(all_isoforms$source == "SQANTI")))
message(paste("  Protein coding:", sum(all_isoforms$is_protein_coding, na.rm = TRUE)))
message(paste("\nOutput:", OUTPUT_FILE))
message("Done!")
