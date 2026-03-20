#!/usr/bin/env Rscript
# ==============================================================================
# 05o_atg_attractiveness.R
#
# Train an ATG attractiveness model on MANE Select transcripts.
#
# For each MANE Select protein-coding transcript:
#   - Find all ATGs in the spliced transcript sequence
#   - Label the annotated CDS start ATG as positive (1), all others as negative (0)
#   - Compute per-ATG features: Kozak context, distance from 5' cap,
#     upstream ATG count, position as fraction of transcript length
#
# Train logistic regression to predict which ATG the ribosome selects.
# Then apply to all isocall isoforms to get per-ATG initiation probabilities.
#
# Inputs:
#   - GENCODE v49 GTF (MANE_Select CDS start positions + exon structure)
#   - GENCODE v49 transcript FASTA (spliced sequences)
#   - Isocall ORF landscape data (from 05m)
#   - Isocall corrected FASTA (spliced sequences)
#
# Output:
#   - data_mashr/analysis_cache/atg_attractiveness.rds
#
# Usage:
#   Rscript 05o_atg_attractiveness.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

GENCODE_GTF <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.gtf.gz"
GENCODE_FASTA <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.transcripts.fa.gz"
ISOCALL_FASTA <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta"
OUTPUT_PATH <- "data_mashr/analysis_cache/atg_attractiveness.rds"

cat("=== ATG Attractiveness Model ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Parse GENCODE GTF for MANE Select transcripts
# ==============================================================================
cat("Parsing GENCODE GTF for MANE Select transcripts...\n")

# Read GTF lines with start_codon and exon features for MANE_Select
gtf_con <- gzfile(GENCODE_GTF, "r")
gtf_lines <- readLines(gtf_con)
close(gtf_con)

# Filter to MANE_Select lines
mane_lines <- gtf_lines[grepl("MANE_Select", gtf_lines, fixed = TRUE)]
cat("  MANE_Select GTF lines:", length(mane_lines), "\n")

# Parse start_codon features to get CDS start positions
start_codon_lines <- mane_lines[grepl("\tstart_codon\t", mane_lines)]
cat("  start_codon features:", length(start_codon_lines), "\n")

parse_gtf_line <- function(line) {
  fields <- strsplit(line, "\t")[[1]]
  chr <- fields[1]
  feature <- fields[3]
  start <- as.integer(fields[4])
  end <- as.integer(fields[5])
  strand <- fields[7]
  attrs <- fields[9]

  # Extract transcript_id
  tx_match <- regmatches(attrs, regexpr('transcript_id "([^"]+)"', attrs))
  tx_id <- sub('transcript_id "', '', sub('"$', '', tx_match))

  # Extract gene_name
  gene_match <- regmatches(attrs, regexpr('gene_name "([^"]+)"', attrs))
  gene_name <- sub('gene_name "', '', sub('"$', '', gene_match))

  list(chr = chr, feature = feature, start = start, end = end,
       strand = strand, tx_id = tx_id, gene_name = gene_name)
}

# Parse start codons
start_codons <- do.call(rbind, lapply(start_codon_lines, function(l) {
  p <- parse_gtf_line(l)
  data.frame(tx_id = p$tx_id, chr = p$chr, strand = p$strand,
             start_codon_start = p$start, start_codon_end = p$end,
             gene_name = p$gene_name, stringsAsFactors = FALSE)
}))
# Keep first start_codon per transcript (some may have multiple exon parts)
start_codons <- start_codons[!duplicated(start_codons$tx_id), ]
cat("  Unique MANE transcripts with start_codon:", nrow(start_codons), "\n")

# Parse exon features for these transcripts
exon_lines <- mane_lines[grepl("\texon\t", mane_lines)]
mane_tx_ids <- start_codons$tx_id

exons_list <- list()
for (line in exon_lines) {
  p <- parse_gtf_line(line)
  if (p$tx_id %in% mane_tx_ids) {
    exons_list[[length(exons_list) + 1]] <- data.frame(
      tx_id = p$tx_id, strand = p$strand,
      exon_start = p$start, exon_end = p$end,
      stringsAsFactors = FALSE)
  }
}
exons_df <- do.call(rbind, exons_list)
cat("  Exon features parsed:", nrow(exons_df), "\n")

# Compute transcript-space CDS start position
# For + strand: CDS start = start_codon_start
# Map genomic CDS start to transcript position using exon structure

compute_tx_pos <- function(genomic_pos, exon_starts, exon_ends, strand) {
  # Exons in transcript order
  if (strand == "-") {
    # Reverse order for minus strand
    exon_starts <- rev(exon_starts)
    exon_ends <- rev(exon_ends)
  }

  cum_pos <- 0L
  for (i in seq_along(exon_starts)) {
    es <- exon_starts[i]
    ee <- exon_ends[i]
    if (strand == "+") {
      if (genomic_pos >= es && genomic_pos <= ee) {
        return(cum_pos + (genomic_pos - es) + 1L)
      }
    } else {
      # For minus strand, genomic position maps in reverse within exon
      if (genomic_pos >= es && genomic_pos <= ee) {
        return(cum_pos + (ee - genomic_pos) + 1L)
      }
    }
    cum_pos <- cum_pos + (ee - es + 1L)
  }
  return(NA_integer_)
}

# For each MANE transcript, compute CDS start in transcript space
cds_start_tx <- integer(nrow(start_codons))
for (i in seq_len(nrow(start_codons))) {
  tx <- start_codons$tx_id[i]
  strand <- start_codons$strand[i]
  ex <- exons_df[exons_df$tx_id == tx, ]
  ex <- ex[order(ex$exon_start), ]

  # CDS start genomic position
  if (strand == "+") {
    cds_genomic <- start_codons$start_codon_start[i]
  } else {
    # For minus strand, the start codon's last position is the 5' end
    cds_genomic <- start_codons$start_codon_end[i]
  }

  cds_start_tx[i] <- compute_tx_pos(cds_genomic, ex$exon_start, ex$exon_end, strand)
}
start_codons$cds_start_tx <- cds_start_tx
start_codons <- start_codons[!is.na(start_codons$cds_start_tx), ]
cat("  Transcripts with mapped CDS start:", nrow(start_codons), "\n")

rm(gtf_lines, mane_lines, exon_lines, exons_list)

# ==============================================================================
# 2. Extract MANE transcript sequences from GENCODE FASTA
# ==============================================================================
cat("\nExtracting MANE transcript sequences from GENCODE FASTA...\n")

# Read FASTA with gzip
fasta_con <- gzfile(GENCODE_FASTA, "r")
fasta_lines <- readLines(fasta_con)
close(fasta_con)

# Parse: headers have format >ENST00000832824.1|ENSG...|...
headers <- which(startsWith(fasta_lines, ">"))
fasta_ids <- sub("\\|.*", "", sub("^>", "", fasta_lines[headers]))

# Only keep MANE transcripts
mane_idx <- which(fasta_ids %in% start_codons$tx_id)
cat("  MANE transcripts found in FASTA:", length(mane_idx), "\n")

mane_seqs <- character(length(mane_idx))
names(mane_seqs) <- fasta_ids[mane_idx]

for (j in seq_along(mane_idx)) {
  h <- headers[mane_idx[j]]
  # Sequence runs from h+1 to next header-1 (or end of file)
  next_h <- if (mane_idx[j] < length(headers)) headers[mane_idx[j] + 1] - 1 else length(fasta_lines)
  mane_seqs[j] <- toupper(paste(fasta_lines[(h + 1):next_h], collapse = ""))
}

cat("  Sequences extracted:", length(mane_seqs), "\n")
rm(fasta_lines)

# Verify CDS start: check that ATG exists at the annotated position
n_verified <- 0L
for (i in seq_len(nrow(start_codons))) {
  tx <- start_codons$tx_id[i]
  if (!tx %in% names(mane_seqs)) next
  seq <- mane_seqs[tx]
  pos <- start_codons$cds_start_tx[i]
  if (pos + 2 <= nchar(seq)) {
    codon <- substr(seq, pos, pos + 2)
    if (codon == "ATG") n_verified <- n_verified + 1L
  }
}
cat("  CDS start verified as ATG:", n_verified, "of", nrow(start_codons),
    sprintf("(%.1f%%)\n", 100 * n_verified / nrow(start_codons)))

# Filter to verified
verified_tx <- character(0)
for (i in seq_len(nrow(start_codons))) {
  tx <- start_codons$tx_id[i]
  if (!tx %in% names(mane_seqs)) next
  seq <- mane_seqs[tx]
  pos <- start_codons$cds_start_tx[i]
  if (pos + 2 <= nchar(seq) && substr(seq, pos, pos + 2) == "ATG") {
    verified_tx <- c(verified_tx, tx)
  }
}
start_codons <- start_codons[start_codons$tx_id %in% verified_tx, ]
cat("  Final training transcripts:", nrow(start_codons), "\n")

# ==============================================================================
# 3. Build training data: all ATGs per transcript, label CDS start
# ==============================================================================
cat("\nBuilding ATG training data...\n")

# Kozak scoring (same as 05m)
score_kozak <- function(seq, atg_pos) {
  pos_m3 <- atg_pos - 3
  pos_p4 <- atg_pos + 3
  n <- nchar(seq)
  nt_m3 <- if (pos_m3 >= 1) substr(seq, pos_m3, pos_m3) else ""
  nt_p4 <- if (pos_p4 <= n) substr(seq, pos_p4, pos_p4) else ""
  has_m3 <- nt_m3 %in% c("A", "G")
  has_p4 <- nt_p4 == "G"
  if (has_m3 && has_p4) return(2L)
  if (has_m3 || has_p4) return(1L)
  return(0L)
}

stop_codons <- c("TAA", "TAG", "TGA")

training_rows <- list()

for (i in seq_len(nrow(start_codons))) {
  tx <- start_codons$tx_id[i]
  seq <- mane_seqs[tx]
  seq_len <- nchar(seq)
  cds_pos <- start_codons$cds_start_tx[i]

  # Find all ATGs
  atg_matches <- gregexpr("ATG", seq, fixed = TRUE)[[1]]
  if (atg_matches[1] == -1L) next
  atg_positions <- as.integer(atg_matches)

  for (j in seq_along(atg_positions)) {
    atg <- atg_positions[j]
    is_cds_start <- (atg == cds_pos)

    # Features
    kozak <- score_kozak(seq, atg)
    dist_from_cap <- atg  # 1-based distance from 5' end
    frac_position <- atg / seq_len
    n_upstream_atgs <- sum(atg_positions < atg)

    # ORF length from this ATG (first in-frame stop)
    orf_len <- NA_integer_
    pos <- atg + 3L
    while (pos + 2L <= seq_len) {
      codon <- substr(seq, pos, pos + 2L)
      if (codon %in% stop_codons) {
        orf_len <- pos - atg
        break
      }
      pos <- pos + 3L
    }

    training_rows[[length(training_rows) + 1]] <- data.frame(
      tx_id = tx,
      atg_position = atg,
      is_cds_start = as.integer(is_cds_start),
      kozak = kozak,
      log_dist_from_cap = log1p(dist_from_cap),
      frac_position = frac_position,
      n_upstream_atgs = n_upstream_atgs,
      log_orf_length = ifelse(is.na(orf_len), 0, log1p(orf_len)),
      has_orf = as.integer(!is.na(orf_len)),
      tx_length = seq_len,
      stringsAsFactors = FALSE
    )
  }

  if (i %% 5000 == 0) cat(sprintf("  Processed %d / %d transcripts\n", i, nrow(start_codons)))
}

train_df <- do.call(rbind, training_rows)
cat("  Total ATGs:", nrow(train_df), "\n")
cat("  CDS starts (positive):", sum(train_df$is_cds_start), "\n")
cat("  Ratio:", sprintf("1:%.0f\n", (nrow(train_df) - sum(train_df$is_cds_start)) / sum(train_df$is_cds_start)))

# ==============================================================================
# 4. Train logistic regression
# ==============================================================================
cat("\nTraining ATG attractiveness model...\n")

# Holdout: chr1, 3, 5, 7 (match other models)
tx_chr <- start_codons[, c("tx_id", "chr")]
train_df <- merge(train_df, tx_chr, by = "tx_id", all.x = TRUE)

holdout_chrs <- c("chr1", "chr3", "chr5", "chr7")
atg_train <- train_df[!train_df$chr %in% holdout_chrs, ]
atg_test <- train_df[train_df$chr %in% holdout_chrs, ]

cat("  Train ATGs:", nrow(atg_train), "(pos:", sum(atg_train$is_cds_start), ")\n")
cat("  Test ATGs:", nrow(atg_test), "(pos:", sum(atg_test$is_cds_start), ")\n")

features <- c("kozak", "log_dist_from_cap", "frac_position",
              "n_upstream_atgs", "log_orf_length", "has_orf")

f <- as.formula(paste("is_cds_start ~", paste(features, collapse = " + ")))
fit <- glm(f, data = atg_train, family = binomial())

cat("\n  Model coefficients:\n")
print(round(coef(fit), 4))

# Evaluate
pred_train <- predict(fit, atg_train, type = "response")
pred_test <- predict(fit, atg_test, type = "response")

roc_train <- pROC::roc(atg_train$is_cds_start, pred_train, quiet = TRUE)
roc_test <- pROC::roc(atg_test$is_cds_start, pred_test, quiet = TRUE)

cat("\n  Train AUC:", round(pROC::auc(roc_train), 4), "\n")
cat("  Test AUC:", round(pROC::auc(roc_test), 4), "\n")

# Within-transcript accuracy: for each transcript, does the highest-scoring
# ATG match the annotated CDS start?
top1_acc <- atg_test %>%
  group_by(tx_id) %>%
  mutate(pred = predict(fit, pick(everything()), type = "response")) %>%
  slice_max(pred, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  summarise(correct = mean(is_cds_start))

cat("  Top-1 accuracy (test):", round(top1_acc$correct, 4), "\n")

# ==============================================================================
# 5. Apply to isocall transcripts
# ==============================================================================
cat("\nApplying ATG attractiveness model to isocall isoforms...\n")

orf_df <- readRDS("data_mashr/analysis_cache/orf_landscape.rds")
structures <- readRDS("data_mashr/structures.rds")
nmd_class <- readRDS("data_mashr/nmd_classification.rds")

target_ids <- orf_df$isoform_id

# Extract isocall sequences
id_file <- tempfile(fileext = ".txt")
writeLines(target_ids, id_file)
awk_cmd <- sprintf(
  "awk 'BEGIN{while((getline id < \"%s\") > 0) want[id]=1} /^>/{if(seq && found) print hdr \"\\t\" seq; hdr=substr($1,2); found=(hdr in want); seq=\"\"; next} found{seq=seq$0} END{if(seq && found) print hdr \"\\t\" seq}' '%s'",
  id_file, ISOCALL_FASTA
)
raw_output <- system(awk_cmd, intern = TRUE)
unlink(id_file)

iso_seqs <- character(0)
if (length(raw_output) > 0) {
  split_lines <- strsplit(raw_output, "\t", fixed = TRUE)
  iso_seqs <- setNames(
    toupper(vapply(split_lines, `[`, character(1), 2)),
    vapply(split_lines, `[`, character(1), 1))
}
cat("  Isocall sequences:", length(iso_seqs), "\n")

# Compute junction positions (reuse logic from 05m)
target_structs <- structures[structures$isoform_id %in% target_ids, ]
junction_list <- setNames(vector("list", nrow(target_structs)), target_structs$isoform_id)
for (i in seq_len(nrow(target_structs))) {
  starts <- target_structs$exon_starts[[i]]
  ends <- target_structs$exon_ends[[i]]
  strand <- target_structs$strand[i]
  n_exons <- length(starts)
  if (n_exons <= 1) { junction_list[[i]] <- integer(0); next }
  exon_lengths <- ends - starts + 1L
  if (strand == "-") exon_lengths <- rev(exon_lengths)
  junction_list[[i]] <- cumsum(exon_lengths)[-n_exons]
}

# For each isocall isoform, score all ATGs and compute NMD risk
cat("  Scoring ATGs and computing NMD risk...\n")

iso_results <- vector("list", length(target_ids))
names(iso_results) <- target_ids

t_start <- proc.time()
for (idx in seq_along(target_ids)) {
  iso_id <- target_ids[idx]
  if (!iso_id %in% names(iso_seqs)) next

  seq <- iso_seqs[iso_id]
  seq_len_val <- nchar(seq)
  junctions <- junction_list[[iso_id]]
  if (is.null(junctions)) junctions <- integer(0)

  # Find all ATGs
  atg_matches <- gregexpr("ATG", seq, fixed = TRUE)[[1]]
  if (atg_matches[1] == -1L) {
    iso_results[[idx]] <- data.frame(
      isoform_id = iso_id, max_nmd_risk = 0, sum_nmd_risk = 0,
      n_atgs_scored = 0L, best_attractive_nmd = 0,
      prob_any_nmd = 0, stringsAsFactors = FALSE)
    next
  }

  atg_positions <- as.integer(atg_matches)

  # Score each ATG
  atg_data <- data.frame(
    kozak = integer(length(atg_positions)),
    log_dist_from_cap = numeric(length(atg_positions)),
    frac_position = numeric(length(atg_positions)),
    n_upstream_atgs = integer(length(atg_positions)),
    log_orf_length = numeric(length(atg_positions)),
    has_orf = integer(length(atg_positions))
  )

  downstream_ejcs <- integer(length(atg_positions))

  for (j in seq_along(atg_positions)) {
    atg <- atg_positions[j]
    atg_data$kozak[j] <- score_kozak(seq, atg)
    atg_data$log_dist_from_cap[j] <- log1p(atg)
    atg_data$frac_position[j] <- atg / seq_len_val
    atg_data$n_upstream_atgs[j] <- sum(atg_positions < atg)

    # Find first in-frame stop
    orf_len <- NA_integer_
    stop_pos <- NA_integer_
    pos <- atg + 3L
    while (pos + 2L <= seq_len_val) {
      codon <- substr(seq, pos, pos + 2L)
      if (codon %in% stop_codons) {
        orf_len <- pos - atg
        stop_pos <- pos
        break
      }
      pos <- pos + 3L
    }

    atg_data$log_orf_length[j] <- ifelse(is.na(orf_len), 0, log1p(orf_len))
    atg_data$has_orf[j] <- as.integer(!is.na(orf_len))

    # Downstream EJCs from stop
    if (!is.na(stop_pos)) {
      stop_end <- stop_pos + 2L
      downstream_ejcs[j] <- sum(junctions > stop_end)
    }
  }

  # ATG attractiveness scores
  attract <- predict(fit, atg_data, type = "response")

  # NMD risk per ATG: attractiveness × (stop has downstream EJCs)
  nmd_indicator <- as.integer(downstream_ejcs > 0)
  nmd_risk_per_atg <- attract * nmd_indicator

  # Transcript-level summaries
  max_nmd_risk <- max(nmd_risk_per_atg)
  sum_nmd_risk <- sum(nmd_risk_per_atg)
  best_attractive_nmd <- max(attract[nmd_indicator == 1], 0)

  # P(NMD across multiple rounds) — simplified: use max attractive NMD-competent ATG
  # P(ever select it in N rounds) ≈ 1 - (1 - p)^N
  # We don't know N, so report the per-round max probability
  prob_any_nmd <- max_nmd_risk

  iso_results[[idx]] <- data.frame(
    isoform_id = iso_id,
    max_nmd_risk = max_nmd_risk,
    sum_nmd_risk = sum_nmd_risk,
    n_atgs_scored = length(atg_positions),
    best_attractive_nmd = best_attractive_nmd,
    prob_any_nmd = prob_any_nmd,
    stringsAsFactors = FALSE
  )

  if (idx %% 10000 == 0) {
    cat(sprintf("    Scored %d / %d isoforms (%.0f sec)\n",
                idx, length(target_ids), (proc.time() - t_start)[3]))
  }
}

iso_scored <- do.call(rbind, iso_results[!sapply(iso_results, is.null)])
cat("  Scored isoforms:", nrow(iso_scored), "\n")

# Add NMD status
nmd_ids <- nmd_class$all_samples$nmd
iso_scored$is_nmd <- as.integer(iso_scored$isoform_id %in% nmd_ids)

# Evaluate as NMD predictor
roc_nmd <- pROC::roc(iso_scored$is_nmd, iso_scored$max_nmd_risk, quiet = TRUE)
roc_sum <- pROC::roc(iso_scored$is_nmd, iso_scored$sum_nmd_risk, quiet = TRUE)
roc_best <- pROC::roc(iso_scored$is_nmd, iso_scored$best_attractive_nmd, quiet = TRUE)

cat("\n=== ATG Attractiveness as NMD Predictor ===\n")
cat("  max_nmd_risk AUC:", round(pROC::auc(roc_nmd), 4), "\n")
cat("  sum_nmd_risk AUC:", round(pROC::auc(roc_sum), 4), "\n")
cat("  best_attractive_nmd AUC:", round(pROC::auc(roc_best), 4), "\n")

# ==============================================================================
# 6. Save
# ==============================================================================
results <- list(
  # ATG attractiveness model (trained on MANE)
  model = fit,
  features = features,
  training_summary = list(
    n_transcripts = nrow(start_codons),
    n_atgs_train = nrow(atg_train),
    n_atgs_test = nrow(atg_test),
    auc_train = as.numeric(pROC::auc(roc_train)),
    auc_test = as.numeric(pROC::auc(roc_test)),
    top1_accuracy = top1_acc$correct
  ),

  # Isocall application
  iso_scores = iso_scored,
  nmd_prediction = list(
    max_nmd_risk_auc = as.numeric(pROC::auc(roc_nmd)),
    sum_nmd_risk_auc = as.numeric(pROC::auc(roc_sum)),
    best_attractive_nmd_auc = as.numeric(pROC::auc(roc_best)),
    roc_max = roc_nmd,
    roc_sum = roc_sum,
    roc_best = roc_best
  ),

  metadata = list(
    holdout_chrs = holdout_chrs,
    run_timestamp = Sys.time()
  )
)

saveRDS(results, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
