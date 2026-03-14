#!/usr/bin/env Rscript
# 04_productive_frameshift_precompute.R
#
# Pre-computes external API queries for productive frameshift characterization.
# Saves results to data_mashr/analysis_cache/productive_frameshift_allsamples.rds
# so the report (03_nmd_analysis_mashr.Rmd) can load cached results without API calls.
#
# Sections:
#   0. Setup & data loading
#   1. Protein sequence extraction (from SQANTI protein FASTA)
#   2. Mass-spec evidence (EBI Proteins API)
#   3. Functional domain mapping (biomaRt + reference-based proxy for novel)
#   4. Conservation & annotation (APPRIS/TSL via biomaRt, UniProt REST API)
#
# Run from: results/isoform_transitions/Version_6.0/isopair_wrapper/
# Usage: Rscript 04_productive_frameshift_precompute.R [--test]
#   --test: run on first 5 pairs only (for validation)

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(biomaRt)
})

args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args

cat("\n")
cat("==================================================================\n")
cat("   04: Productive Frameshift Pre-computation\n")
if (test_mode) cat("   *** TEST MODE: 5 pairs only ***\n")
cat("==================================================================\n\n")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 0: Setup & data loading
# ═══════════════════════════════════════════════════════════════════════════════

cat("--- Section 0: Loading data ---\n\n")

data_dir <- "data_mashr"
cache_dir <- file.path(data_dir, "analysis_cache")

np <- readRDS(file.path(cache_dir, "novel_protein_c4_v2_allsamples.rds"))
fc <- readRDS(file.path(cache_dir, "fc_c4_allsamples.rds"))
fw <- readRDS(file.path(cache_dir, "fw_c4_allsamples.rds"))
cds <- readRDS(file.path(data_dir, "cds.rds"))
structures <- readRDS(file.path(data_dir, "structures.rds"))

cat(sprintf("  Novel protein pairs: %d\n", nrow(np)))

# Merge reference_isoform_id from fc$pair_summary
np <- merge(np, fc$pair_summary[, c("gene_id", "reference_isoform_id", "comparator_isoform_id")],
            by = c("gene_id", "comparator_isoform_id"))

# Gene symbol lookup from mashr DE files
de_dir <- "/Users/petecastaldi/claude_projects/nmd/isocall_dge/mashr"
gene_sym <- character(0)
for (f in list.files(de_dir, pattern = "^nmd_mashr_die_.*\\.csv$", full.names = TRUE)) {
  tmp <- read.csv(f, stringsAsFactors = FALSE)[, c("gene_id", "hgnc_symbol")]
  tmp <- tmp[!duplicated(tmp$gene_id) & !is.na(tmp$hgnc_symbol) & tmp$hgnc_symbol != "", ]
  new <- setdiff(tmp$gene_id, names(gene_sym))
  gene_sym[tmp$gene_id[tmp$gene_id %in% new]] <- tmp$hgnc_symbol[tmp$gene_id %in% new]
}
np$gene_symbol <- gene_sym[np$gene_id]
np$is_novel_isoform <- grepl("\\.novel[0-9]+$", np$comparator_isoform_id)

cat(sprintf("  With gene symbols: %d / %d\n", sum(!is.na(np$gene_symbol)), nrow(np)))
cat(sprintf("  Novel comparators: %d, ENST comparators: %d\n",
            sum(np$is_novel_isoform), sum(!np$is_novel_isoform)))

if (test_mode) {
  np <- np[1:min(5, nrow(np)), ]
  cat(sprintf("  TEST MODE: using %d pairs\n", nrow(np)))
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Protein sequence extraction
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 1: Extracting protein sequences ---\n\n")

protein_fasta <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.faa"
stopifnot(file.exists(protein_fasta))

# Collect all isoform IDs to extract
all_ids <- unique(c(np$reference_isoform_id, np$comparator_isoform_id))
cat(sprintf("  Unique isoforms to extract: %d\n", length(all_ids)))

# Write IDs to temp file, then use awk for single-pass extraction
id_tmpfile <- tempfile(fileext = ".txt")
writeLines(all_ids, id_tmpfile)

# awk: read IDs into array, then extract matching sequences
awk_cmd <- sprintf(
  "awk 'NR==FNR{ids[$1]; next} /^>/{id=$1; sub(/^>/,\"\",id); sub(/\\t.*/,\"\",id); p=(id in ids)} p' '%s' '%s'",
  id_tmpfile, protein_fasta
)
fasta_lines <- system(awk_cmd, intern = TRUE)
unlink(id_tmpfile)

# Parse FASTA lines into named sequences
seq_list <- list()
current_id <- NULL
current_seq <- character(0)
for (line in fasta_lines) {
  if (startsWith(line, ">")) {
    if (!is.null(current_id)) {
      seq_list[[current_id]] <- paste0(current_seq, collapse = "")
    }
    # Extract ID: first field after >, before tab
    header <- sub("^>", "", line)
    current_id <- sub("\t.*", "", header)
    current_seq <- character(0)
  } else {
    current_seq <- c(current_seq, line)
  }
}
if (!is.null(current_id)) {
  seq_list[[current_id]] <- paste0(current_seq, collapse = "")
}

cat(sprintf("  Extracted sequences: %d / %d\n", length(seq_list), length(all_ids)))
missing_ids <- setdiff(all_ids, names(seq_list))
if (length(missing_ids) > 0) {
  cat(sprintf("  WARNING: %d IDs not found in protein FASTA: %s\n",
              length(missing_ids), paste(head(missing_ids, 5), collapse = ", ")))
}

# Build characterization table
# NOTE: conserved_bp measures CDS nucleotides in the COMPARATOR upstream of the
# genomic frameshift boundary. conserved_bp/3 = position in the comparator's protein
# where the reading frame shifts. This is NOT the same as shared amino acids between
# ref and comp, because the isoforms can have different exonic content upstream of the
# frameshift (same reading frame, different protein sequence).
frameshift_position_aa <- np$conserved_bp / 3

characterization <- data.frame(
  gene_id = np$gene_id,
  gene_symbol = np$gene_symbol,
  reference_isoform_id = np$reference_isoform_id,
  comparator_isoform_id = np$comparator_isoform_id,
  is_novel_isoform = np$is_novel_isoform,
  total_aa = np$total_aa,
  frameshift_position_aa = round(frameshift_position_aa, 1),
  novel_aa = np$novel_aa,
  pct_novel = np$pct_novel,
  stringsAsFactors = FALSE
)

# Add protein sequences
characterization$ref_protein_seq <- seq_list[np$reference_isoform_id]
characterization$comp_protein_seq <- seq_list[np$comparator_isoform_id]

# Replace NULL with NA for missing sequences
characterization$ref_protein_seq[sapply(characterization$ref_protein_seq, is.null)] <- NA_character_
characterization$comp_protein_seq[sapply(characterization$comp_protein_seq, is.null)] <- NA_character_
characterization$ref_protein_seq <- unlist(characterization$ref_protein_seq)
characterization$comp_protein_seq <- unlist(characterization$comp_protein_seq)

# Compute actual_shared_aa: the real protein-level identity between ref and comp.
# This is the position where the protein sequences first diverge, determined by
# direct sequence comparison. Differs from frameshift_position_aa when isoforms
# have different exonic content upstream of the frameshift boundary.
characterization$actual_shared_aa <- NA_integer_
n_exact_match <- 0
n_early_diverge <- 0
for (i in seq_len(nrow(characterization))) {
  ref_seq <- characterization$ref_protein_seq[i]
  comp_seq <- characterization$comp_protein_seq[i]
  if (is.na(ref_seq) || is.na(comp_seq)) next

  ref_clean <- sub("\\*$", "", ref_seq)
  comp_clean <- sub("\\*$", "", comp_seq)
  ref_chars <- strsplit(ref_clean, "")[[1]]
  comp_chars <- strsplit(comp_clean, "")[[1]]
  min_len <- min(length(ref_chars), length(comp_chars))

  shared <- min_len  # default: identical up to shorter sequence
  for (j in seq_len(min_len)) {
    if (ref_chars[j] != comp_chars[j]) {
      shared <- j - 1L
      break
    }
  }
  characterization$actual_shared_aa[i] <- shared

  fs_pos <- floor(characterization$frameshift_position_aa[i])
  if (shared >= fs_pos) {
    n_exact_match <- n_exact_match + 1
  } else {
    n_early_diverge <- n_early_diverge + 1
  }
}

cat(sprintf("  Sequence comparison: %d pairs where proteins match up to frameshift position\n",
            n_exact_match))
cat(sprintf("  %d pairs where proteins diverge BEFORE frameshift (different upstream exons)\n",
            n_early_diverge))
cat(sprintf("  Median actual_shared_aa: %.0f (vs median frameshift_position_aa: %.0f)\n",
            median(characterization$actual_shared_aa, na.rm = TRUE),
            median(characterization$frameshift_position_aa, na.rm = TRUE)))

# Identify frameshift-causing events
fw_events <- fw$events
fs_events <- fw_events[fw_events$is_frameshift == TRUE, ]

# For each pair, get the event types that cause frameshifts
event_info <- lapply(seq_len(nrow(np)), function(i) {
  pair_fs <- fs_events[fs_events$gene_id == np$gene_id[i] &
                          fs_events$comparator_isoform_id == np$comparator_isoform_id[i], ]
  if (nrow(pair_fs) == 0) {
    return(data.frame(frameshift_event_types = NA_character_,
                      frameshift_event_directions = NA_character_,
                      n_frameshift_events = 0L, stringsAsFactors = FALSE))
  }
  data.frame(
    frameshift_event_types = paste(unique(pair_fs$event_type), collapse = ","),
    frameshift_event_directions = paste(unique(pair_fs$direction), collapse = ","),
    n_frameshift_events = nrow(pair_fs),
    stringsAsFactors = FALSE
  )
})
event_df <- do.call(rbind, event_info)
characterization <- cbind(characterization, event_df)

cat(sprintf("  Pairs with frameshift events: %d / %d\n",
            sum(characterization$n_frameshift_events > 0), nrow(characterization)))


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Mass-spec evidence (EBI Proteins API)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 2: Querying mass-spec evidence ---\n\n")

# In-silico tryptic digest: split at K/R not followed by P, keep 7-30 aa
tryptic_digest <- function(seq) {
  # Remove stop codon marker
  seq <- sub("\\*$", "", seq)
  # Split at K or R not followed by P
  # Insert separator after K/R when not followed by P
  marked <- gsub("([KR])([^P])", "\\1|\\2", seq)
  marked <- gsub("([KR])$", "\\1|", marked)  # handle terminal K/R
  peptides <- unlist(strsplit(marked, "\\|"))
  peptides <- peptides[nchar(peptides) >= 7 & nchar(peptides) <= 30]
  peptides
}

# Query EBI proteomics API for a peptide
query_ebi_proteomics <- function(peptide) {
  url <- sprintf("https://www.ebi.ac.uk/proteins/api/proteomics?peptide=%s", peptide)
  tryCatch({
    resp <- GET(url, accept("application/json"), timeout(30))
    if (status_code(resp) == 200) {
      body <- content(resp, as = "text", encoding = "UTF-8")
      if (nchar(body) > 2) {  # non-empty JSON array
        return(list(status = "success", data = fromJSON(body)))
      }
      return(list(status = "empty", data = NULL))
    }
    return(list(status = paste0("http_", status_code(resp)), data = NULL))
  }, error = function(e) {
    return(list(status = "error", data = NULL))
  })
}

massspec_rows <- list()
n_queries <- 0

for (i in seq_len(nrow(characterization))) {
  comp_seq <- characterization$comp_protein_seq[i]
  ref_seq <- characterization$ref_protein_seq[i]
  fs_pos <- floor(characterization$frameshift_position_aa[i])
  gene_id <- characterization$gene_id[i]
  comp_id <- characterization$comparator_isoform_id[i]
  gene_symbol <- characterization$gene_symbol[i]

  if (is.na(comp_seq) || is.na(ref_seq) || fs_pos < 1) next

  # Extract novel region (after frameshift position in comparator)
  novel_region <- substr(comp_seq, fs_pos + 1, nchar(sub("\\*$", "", comp_seq)))
  if (nchar(novel_region) < 7) next

  # Tryptic digest of novel region
  novel_peptides <- tryptic_digest(novel_region)
  if (length(novel_peptides) == 0) next

  # Filter to peptides unique to novel region (not in reference)
  ref_clean <- sub("\\*$", "", ref_seq)
  unique_peptides <- novel_peptides[!sapply(novel_peptides, function(p) grepl(p, ref_clean, fixed = TRUE))]

  if (length(unique_peptides) == 0) next

  cat(sprintf("  [%d/%d] %s (%s): %d unique peptides to query\n",
              i, nrow(characterization),
              ifelse(is.na(gene_symbol), gene_id, gene_symbol),
              comp_id, length(unique_peptides)))

  for (pep in unique_peptides) {
    result <- query_ebi_proteomics(pep)
    n_queries <- n_queries + 1

    ebi_hit <- FALSE
    ebi_gene_match <- FALSE
    ebi_accession <- NA_character_

    if (result$status == "success" && !is.null(result$data)) {
      ebi_hit <- TRUE
      # Check if any hit maps to the same gene
      if (is.data.frame(result$data) && "accession" %in% names(result$data)) {
        ebi_accession <- paste(unique(result$data$accession), collapse = ";")
        # Check gene name match (best effort)
        if ("gene" %in% names(result$data) && !is.na(gene_symbol)) {
          gene_names <- unlist(result$data$gene)
          if (any(toupper(gene_names) == toupper(gene_symbol))) {
            ebi_gene_match <- TRUE
          }
        }
      }
    }

    massspec_rows[[length(massspec_rows) + 1]] <- data.frame(
      gene_id = gene_id,
      comparator_isoform_id = comp_id,
      peptide = pep,
      peptide_length = nchar(pep),
      is_unique_to_novel = TRUE,
      ebi_hit = ebi_hit,
      ebi_gene_match = ebi_gene_match,
      ebi_accession = ebi_accession,
      query_status = result$status,
      stringsAsFactors = FALSE
    )

    Sys.sleep(0.5)  # Rate limiting
  }
}

massspec <- if (length(massspec_rows) > 0) do.call(rbind, massspec_rows) else {
  data.frame(gene_id = character(0), comparator_isoform_id = character(0),
             peptide = character(0), peptide_length = integer(0),
             is_unique_to_novel = logical(0), ebi_hit = logical(0),
             ebi_gene_match = logical(0), ebi_accession = character(0),
             query_status = character(0), stringsAsFactors = FALSE)
}

cat(sprintf("\n  Total API queries: %d\n", n_queries))
cat(sprintf("  Peptides with hits: %d / %d\n", sum(massspec$ebi_hit), nrow(massspec)))
cat(sprintf("  Gene-matched hits: %d\n", sum(massspec$ebi_gene_match)))


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Functional domain mapping
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 3: Querying domain annotations ---\n\n")

# 3a. ENST isoforms — bulk biomaRt query
all_enst <- unique(c(
  characterization$reference_isoform_id[grepl("^ENST", characterization$reference_isoform_id)],
  characterization$comparator_isoform_id[grepl("^ENST", characterization$comparator_isoform_id)]
))
cat(sprintf("  ENST isoforms for domain query: %d\n", length(all_enst)))

# Strip version for biomaRt (Ensembl uses unversioned IDs)
enst_base <- sub("\\.[0-9]+$", "", all_enst)

# Connect to Ensembl
cat("  Connecting to Ensembl biomaRt...\n")
mart <- tryCatch({
  useMart("ensembl", dataset = "hsapiens_gene_ensembl")
}, error = function(e) {
  cat(sprintf("  WARNING: biomaRt connection failed: %s\n", e$message))
  cat("  Trying archive...\n")
  tryCatch(
    useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast"),
    error = function(e2) { cat(sprintf("  Archive also failed: %s\n", e2$message)); NULL }
  )
})

domain_rows <- list()
domain_result <- NULL

if (!is.null(mart)) {
  cat("  Querying Pfam/InterPro domains...\n")
  domain_result <- tryCatch({
    getBM(
      attributes = c("ensembl_transcript_id", "pfam", "pfam_start", "pfam_end",
                      "interpro", "interpro_short_description",
                      "interpro_start", "interpro_end"),
      filters = "ensembl_transcript_id",
      values = enst_base,
      mart = mart
    )
  }, error = function(e) {
    cat(sprintf("  WARNING: domain query failed: %s\n", e$message))
    NULL
  })

  if (!is.null(domain_result) && nrow(domain_result) > 0) {
    cat(sprintf("  Domain annotations retrieved: %d rows\n", nrow(domain_result)))

    # Map back to versioned IDs and classify relative to frameshift boundary
    # For each isoform's OWN domains, classify using frameshift_position_aa
    # (the position in the comparator where the reading frame shifts).
    # For reference domains, this is an approximation — the frameshift position in
    # reference protein coordinates may differ if the isoforms have different upstream exons.
    for (idx in seq_len(nrow(characterization))) {
      ref_id <- characterization$reference_isoform_id[idx]
      comp_id <- characterization$comparator_isoform_id[idx]
      fs_pos <- floor(characterization$frameshift_position_aa[idx])

      for (role in c("reference", "comparator")) {
        iso_id <- if (role == "reference") ref_id else comp_id
        iso_base <- sub("\\.[0-9]+$", "", iso_id)

        if (!grepl("^ENST", iso_id)) next  # Skip novel isoforms

        iso_domains <- domain_result[domain_result$ensembl_transcript_id == iso_base, ]
        if (nrow(iso_domains) == 0) next

        # Process Pfam domains
        pfam_rows <- iso_domains[!is.na(iso_domains$pfam) & iso_domains$pfam != "", ]
        if (nrow(pfam_rows) > 0) {
          for (j in seq_len(nrow(pfam_rows))) {
            ds <- pfam_rows$pfam_start[j]
            de <- pfam_rows$pfam_end[j]
            pos_class <- if (de <= fs_pos) "conserved" else if (ds > fs_pos) "novel_region" else "disrupted"

            domain_rows[[length(domain_rows) + 1]] <- data.frame(
              gene_id = characterization$gene_id[idx],
              isoform_id = iso_id, isoform_role = role,
              domain_id = pfam_rows$pfam[j],
              domain_name = pfam_rows$interpro_short_description[j],
              domain_source = "Pfam",
              domain_start_aa = ds, domain_end_aa = de,
              position_class = pos_class,
              stringsAsFactors = FALSE
            )
          }
        }

        # Process InterPro domains (non-Pfam)
        ipr_rows <- iso_domains[!is.na(iso_domains$interpro) & iso_domains$interpro != "" &
                                  (is.na(iso_domains$pfam) | iso_domains$pfam == ""), ]
        if (nrow(ipr_rows) > 0) {
          for (j in seq_len(nrow(ipr_rows))) {
            ds <- ipr_rows$interpro_start[j]
            de <- ipr_rows$interpro_end[j]
            pos_class <- if (de <= fs_pos) "conserved" else if (ds > fs_pos) "novel_region" else "disrupted"

            domain_rows[[length(domain_rows) + 1]] <- data.frame(
              gene_id = characterization$gene_id[idx],
              isoform_id = iso_id, isoform_role = role,
              domain_id = ipr_rows$interpro[j],
              domain_name = ipr_rows$interpro_short_description[j],
              domain_source = "InterPro",
              domain_start_aa = ds, domain_end_aa = de,
              position_class = pos_class,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
}

# 3b. Novel isoform domains via reference proxy
# For novel comparators (and novel references), use the paired ENST isoform's
# domain annotations restricted to the region where BOTH proteins are actually
# identical (actual_shared_aa), not just the frameshift position.
n_proxy <- 0
n_proxy_skipped_early_diverge <- 0
for (idx in seq_len(nrow(characterization))) {
  ref_id <- characterization$reference_isoform_id[idx]
  comp_id <- characterization$comparator_isoform_id[idx]
  shared_aa <- characterization$actual_shared_aa[idx]
  if (is.na(shared_aa) || shared_aa < 1) next

  # Check if either isoform is novel
  ref_is_novel <- !grepl("^ENST", ref_id)
  comp_is_novel <- !grepl("^ENST", comp_id)

  if (!ref_is_novel && !comp_is_novel) next  # Both ENST, already handled
  if (ref_is_novel && comp_is_novel) next     # Both novel, no proxy available

  # Use the ENST partner's domains restricted to actually-shared region
  proxy_enst <- if (comp_is_novel) ref_id else comp_id
  proxy_base <- sub("\\.[0-9]+$", "", proxy_enst)
  novel_id <- if (comp_is_novel) comp_id else ref_id
  novel_role <- if (comp_is_novel) "comparator" else "reference"

  if (is.null(domain_result) || nrow(domain_result) == 0) next
  proxy_domains <- domain_result[domain_result$ensembl_transcript_id == proxy_base, ]
  if (nrow(proxy_domains) == 0) next

  # Only include domains entirely within the actually-shared protein region
  pfam_rows <- proxy_domains[!is.na(proxy_domains$pfam) & proxy_domains$pfam != "", ]
  if (nrow(pfam_rows) > 0) {
    for (j in seq_len(nrow(pfam_rows))) {
      ds <- pfam_rows$pfam_start[j]
      de <- pfam_rows$pfam_end[j]
      if (de > shared_aa) {
        n_proxy_skipped_early_diverge <- n_proxy_skipped_early_diverge + 1
        next
      }

      domain_rows[[length(domain_rows) + 1]] <- data.frame(
        gene_id = characterization$gene_id[idx],
        isoform_id = novel_id, isoform_role = novel_role,
        domain_id = pfam_rows$pfam[j],
        domain_name = pfam_rows$interpro_short_description[j],
        domain_source = "Pfam (ref proxy)",
        domain_start_aa = ds, domain_end_aa = de,
        position_class = "conserved",
        stringsAsFactors = FALSE
      )
      n_proxy <- n_proxy + 1
    }
  }

  ipr_rows <- proxy_domains[!is.na(proxy_domains$interpro) & proxy_domains$interpro != "" &
                              (is.na(proxy_domains$pfam) | proxy_domains$pfam == ""), ]
  if (nrow(ipr_rows) > 0) {
    for (j in seq_len(nrow(ipr_rows))) {
      ds <- ipr_rows$interpro_start[j]
      de <- ipr_rows$interpro_end[j]
      if (de > shared_aa) {
        n_proxy_skipped_early_diverge <- n_proxy_skipped_early_diverge + 1
        next
      }

      domain_rows[[length(domain_rows) + 1]] <- data.frame(
        gene_id = characterization$gene_id[idx],
        isoform_id = novel_id, isoform_role = novel_role,
        domain_id = ipr_rows$interpro[j],
        domain_name = ipr_rows$interpro_short_description[j],
        domain_source = "InterPro (ref proxy)",
        domain_start_aa = ds, domain_end_aa = de,
        position_class = "conserved",
        stringsAsFactors = FALSE
      )
      n_proxy <- n_proxy + 1
    }
  }
}
cat(sprintf("  Reference-proxy domain annotations for novel isoforms: %d\n", n_proxy))
cat(sprintf("  Proxy domains skipped (beyond actual_shared_aa boundary): %d\n", n_proxy_skipped_early_diverge))

domains <- if (length(domain_rows) > 0) do.call(rbind, domain_rows) else {
  data.frame(gene_id = character(0), isoform_id = character(0),
             isoform_role = character(0), domain_id = character(0),
             domain_name = character(0), domain_source = character(0),
             domain_start_aa = integer(0), domain_end_aa = integer(0),
             position_class = character(0), stringsAsFactors = FALSE)
}
# Deduplicate
domains <- unique(domains)
cat(sprintf("  Total domain annotations: %d\n", nrow(domains)))
cat(sprintf("  Position classes: %s\n",
            paste(names(table(domains$position_class)), table(domains$position_class),
                  sep = "=", collapse = ", ")))


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Conservation & annotation
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 4: Querying conservation & annotation ---\n\n")

# 4a. APPRIS and TSL via biomaRt
conservation <- data.frame(
  gene_id = characterization$gene_id,
  reference_isoform_id = characterization$reference_isoform_id,
  comparator_isoform_id = characterization$comparator_isoform_id,
  ref_appris = NA_character_, comp_appris = NA_character_,
  ref_tsl = NA_character_, comp_tsl = NA_character_,
  ref_is_canonical = NA, comp_is_canonical = NA,
  comp_uniprot_accession = NA_character_,
  comp_uniprot_reviewed = NA_character_,
  comp_protein_existence = NA_character_,
  stringsAsFactors = FALSE
)

if (!is.null(mart)) {
  cat("  Querying APPRIS/TSL annotations...\n")
  appris_result <- tryCatch({
    getBM(
      attributes = c("ensembl_transcript_id", "transcript_appris",
                      "transcript_tsl", "transcript_is_canonical"),
      filters = "ensembl_transcript_id",
      values = enst_base,
      mart = mart
    )
  }, error = function(e) {
    cat(sprintf("  WARNING: APPRIS query failed: %s\n", e$message))
    NULL
  })

  if (!is.null(appris_result) && nrow(appris_result) > 0) {
    cat(sprintf("  APPRIS/TSL annotations: %d rows\n", nrow(appris_result)))

    # Map to versioned IDs
    appris_lookup <- setNames(appris_result$transcript_appris,
                               appris_result$ensembl_transcript_id)
    tsl_lookup <- setNames(appris_result$transcript_tsl,
                            appris_result$ensembl_transcript_id)
    canon_lookup <- setNames(appris_result$transcript_is_canonical,
                              appris_result$ensembl_transcript_id)

    for (i in seq_len(nrow(conservation))) {
      ref_base <- sub("\\.[0-9]+$", "", conservation$reference_isoform_id[i])
      comp_base <- sub("\\.[0-9]+$", "", conservation$comparator_isoform_id[i])

      if (ref_base %in% names(appris_lookup)) {
        conservation$ref_appris[i] <- appris_lookup[ref_base]
        conservation$ref_tsl[i] <- tsl_lookup[ref_base]
        conservation$ref_is_canonical[i] <- canon_lookup[ref_base]
      }
      if (comp_base %in% names(appris_lookup)) {
        conservation$comp_appris[i] <- appris_lookup[comp_base]
        conservation$comp_tsl[i] <- tsl_lookup[comp_base]
        conservation$comp_is_canonical[i] <- canon_lookup[comp_base]
      }
    }
  }
}

# 4b. UniProt protein existence (REST API) — ENST comparators only
enst_comps <- characterization$comparator_isoform_id[!characterization$is_novel_isoform]
cat(sprintf("  Querying UniProt for %d ENST comparators...\n", length(enst_comps)))

for (ui in seq_along(enst_comps)) {
  comp_id <- enst_comps[ui]
  if (ui %% 10 == 1 || ui == length(enst_comps)) {
    cat(sprintf("  [%d/%d] UniProt query: %s\n", ui, length(enst_comps), comp_id))
  }
  enst_base_id <- sub("\\.[0-9]+$", "", comp_id)
  url <- sprintf(
    "https://rest.uniprot.org/uniprotkb/search?query=xref:ensembl-%s&fields=accession,protein_name,protein_existence,annotation_score,reviewed&format=json&size=5",
    enst_base_id
  )

  result <- tryCatch({
    resp <- GET(url, timeout(30))
    if (status_code(resp) == 200) {
      body <- content(resp, as = "text", encoding = "UTF-8")
      fromJSON(body)
    } else NULL
  }, error = function(e) NULL)

  idx <- which(conservation$comparator_isoform_id == comp_id)
  if (!is.null(result) && "results" %in% names(result) && length(result$results) > 0) {
    res <- result$results
    if (is.data.frame(res) && nrow(res) > 0) {
      conservation$comp_uniprot_accession[idx] <- res$primaryAccession[1]
      conservation$comp_uniprot_reviewed[idx] <- ifelse(
        isTRUE(res$entryType[1] == "UniProtKB reviewed (Swiss-Prot)"),
        "Swiss-Prot", "TrEMBL"
      )
      if ("proteinExistence" %in% names(res)) {
        conservation$comp_protein_existence[idx] <- res$proteinExistence[1]
      }
    }
  }

  Sys.sleep(1)  # Rate limiting
}

n_uniprot <- sum(!is.na(conservation$comp_uniprot_accession))
cat(sprintf("  UniProt matches: %d / %d ENST comparators\n", n_uniprot, length(enst_comps)))


# ═══════════════════════════════════════════════════════════════════════════════
# SAVE RESULTS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- Saving results ---\n\n")

api_failures <- list(
  massspec_errors = sum(massspec$query_status == "error"),
  domain_query_failed = is.null(mart),
  uniprot_missing = sum(is.na(conservation$comp_uniprot_accession[!characterization$is_novel_isoform]))
)

output <- list(
  characterization = characterization,
  massspec = massspec,
  domains = domains,
  conservation = conservation,
  metadata = list(
    run_date = Sys.time(),
    n_pairs = nrow(characterization),
    fasta_file = protein_fasta,
    api_failures = api_failures,
    test_mode = test_mode
  )
)

out_file <- file.path(cache_dir, "productive_frameshift_allsamples.rds")
saveRDS(output, out_file)
cat(sprintf("  Saved to: %s\n", out_file))
cat(sprintf("  Characterization: %d pairs\n", nrow(characterization)))
cat(sprintf("  Mass-spec queries: %d peptides, %d hits\n", nrow(massspec), sum(massspec$ebi_hit)))
cat(sprintf("  Domain annotations: %d entries\n", nrow(domains)))
cat(sprintf("  Conservation: %d pairs annotated\n", nrow(conservation)))
cat("\nDone.\n")
