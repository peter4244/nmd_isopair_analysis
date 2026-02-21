#!/usr/bin/env Rscript
################################################################################
# Test Script: End-to-end event detection + reconstruction on GENCODE test GTF
################################################################################
#
# Parses a GENCODE GTF, runs event detection on all isoform pairs (with
# randomly assigned dominant/comparator roles), and verifies reconstruction.
#
# Usage:
#   Rscript scripts/dev/test_08_on_gencode.R
#   Rscript scripts/dev/test_08_on_gencode.R --gtf path/to/file.gtf
#   Rscript scripts/dev/test_08_on_gencode.R --seed 123
#
# Default GTF: development/reconstruction/real_data/gencode_test.gtf
#
################################################################################

library(tidyverse)

cat("\n")
cat("========================================================================\n")
cat("  GENCODE TEST: End-to-end event detection + reconstruction\n")
cat("========================================================================\n")
cat("\n")

base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"
setwd(base_dir)

# Parse CLI args
args <- commandArgs(trailingOnly = TRUE)
gtf_file <- "development/reconstruction/real_data/gencode_test.gtf"
seed <- 42L
if ("--gtf" %in% args) {
  idx <- which(args == "--gtf")
  if (idx < length(args)) gtf_file <- args[idx + 1]
}
if ("--seed" %in% args) {
  idx <- which(args == "--seed")
  if (idx < length(args)) seed <- as.integer(args[idx + 1])
}

# ==============================================================================
# 0. Source core libraries
# ==============================================================================

source("scripts/core/event_detection_functions.R")
source("scripts/core/reconstruction_functions.R")
cat("  Loaded event_detection_functions.R + reconstruction_functions.R\n")

# ==============================================================================
# 1. Parse GTF -> isoform exon structures
# ==============================================================================

cat(sprintf("\nParsing GTF: %s\n", gtf_file))

lines <- readLines(gtf_file)
lines <- lines[!startsWith(lines, '#') & lines != '']

gtf_results <- list()
for (line in lines) {
  parts <- strsplit(line, '\t')[[1]]
  if (parts[3] != "exon") next

  attrs_str <- parts[9]
  gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', attrs_str)
  transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', attrs_str)

  gtf_results[[length(gtf_results) + 1]] <- tibble(
    seqnames = parts[1],
    exon_start = as.integer(parts[4]),
    exon_end = as.integer(parts[5]),
    strand = parts[7],
    gene_id = gene_id,
    isoform_id = transcript_id
  )
}

isoform_structures <- bind_rows(gtf_results)

# Assign exon numbers (transcriptional order)
isoform_structures <- isoform_structures %>%
  group_by(isoform_id, strand) %>%
  mutate(exon_number = if (first(strand) == "+") {
    rank(exon_start, ties.method = "first")
  } else {
    rank(desc(exon_start), ties.method = "first")
  }) %>%
  ungroup()

n_tx <- n_distinct(isoform_structures$isoform_id)
n_genes <- n_distinct(isoform_structures$gene_id)
cat(sprintf("  Parsed %d exons, %d transcripts, %d genes\n",
            nrow(isoform_structures), n_tx, n_genes))

# ==============================================================================
# 2. Form pairs with RANDOM dominant/comparator assignment
# ==============================================================================

cat(sprintf("\nForming pairs (random dominant assignment, seed=%d)...\n", seed))
set.seed(seed)

# For each gene with 2+ isoforms, randomly pick one as dominant
# Then pair it against all others
genes_multi <- isoform_structures %>%
  distinct(gene_id, isoform_id) %>%
  count(gene_id) %>%
  filter(n >= 2) %>%
  pull(gene_id)

all_pairs <- map_dfr(genes_multi, function(g) {
  isos <- isoform_structures %>%
    filter(gene_id == g) %>%
    distinct(isoform_id) %>%
    pull(isoform_id)
  dom <- sample(isos, 1)
  comps <- setdiff(isos, dom)
  tibble(
    gene_id = g,
    dominant_isoform_id = dom,
    comparator_isoform_id = comps
  )
})

cat(sprintf("  %d pairs across %d genes\n", nrow(all_pairs), length(genes_multi)))

# ==============================================================================
# 3. Build union exons (inline, same logic as core/03)
# ==============================================================================

cat("\nBuilding union exons...\n")

build_gene_union_exons <- function(gene_exons) {
  boundaries <- sort(unique(c(gene_exons$exon_start, gene_exons$exon_end)))
  if (length(boundaries) < 2) {
    return(tibble(
      start = boundaries[1], end = boundaries[1],
      union_exon_id = paste0(gene_exons$gene_id[1], "_UE_1")
    ))
  }

  starts_set <- unique(gene_exons$exon_start)
  ends_set <- unique(gene_exons$exon_end)

  segments <- list()
  current_start <- boundaries[1]

  for (b in boundaries[-c(1, length(boundaries))]) {
    is_start <- b %in% starts_set
    is_end <- b %in% ends_set
    is_both <- is_start && is_end

    if (is_both) {
      if (current_start <= b - 1) {
        segments[[length(segments) + 1]] <- c(current_start, b - 1)
      }
      segments[[length(segments) + 1]] <- c(b, b)
      current_start <- b + 1
    } else if (is_start) {
      if (current_start <= b - 1) {
        segments[[length(segments) + 1]] <- c(current_start, b - 1)
      }
      current_start <- b
    } else if (is_end) {
      segments[[length(segments) + 1]] <- c(current_start, b)
      current_start <- b + 1
    }
  }

  # Final segment
  if (current_start <= boundaries[length(boundaries)]) {
    segments[[length(segments) + 1]] <- c(current_start, boundaries[length(boundaries)])
  }

  # Filter to segments covered by at least one exon
  seg_df <- tibble(
    start = sapply(segments, `[`, 1),
    end = sapply(segments, `[`, 2)
  ) %>% filter(start <= end)

  covered <- map_lgl(seq_len(nrow(seg_df)), function(i) {
    any(gene_exons$exon_start <= seg_df$start[i] & gene_exons$exon_end >= seg_df$end[i])
  })
  seg_df <- seg_df[covered, ]

  gene <- gene_exons$gene_id[1]
  seg_df$union_exon_id <- paste0(gene, "_UE_", seq_len(nrow(seg_df)))
  seg_df
}

union_exons_list <- isoform_structures %>%
  group_by(gene_id) %>%
  group_split() %>%
  map(function(g) {
    ue <- build_gene_union_exons(g)
    ue$gene_id <- g$gene_id[1]
    ue$chr <- g$seqnames[1]
    ue$strand <- g$strand[1]
    ue
  })

union_exons <- bind_rows(union_exons_list)
cat(sprintf("  Built %d union exons for %d genes\n",
            nrow(union_exons), n_distinct(union_exons$gene_id)))

# ==============================================================================
# 4. Run event detection + reconstruction verification
# ==============================================================================

cat("\nRunning event detection + reconstruction...\n")

pass_count <- 0L
fail_count <- 0L
error_count <- 0L

results_list <- list()

for (i in seq_len(nrow(all_pairs))) {
  pair <- all_pairs[i, ]
  gene <- pair$gene_id
  dom_id <- pair$dominant_isoform_id
  comp_id <- pair$comparator_isoform_id

  result <- tryCatch({
    dom_exons <- isoform_structures %>%
      filter(isoform_id == dom_id) %>%
      arrange(exon_start)
    comp_exons <- isoform_structures %>%
      filter(isoform_id == comp_id) %>%
      arrange(exon_start)

    if (nrow(dom_exons) == 0 || nrow(comp_exons) == 0) {
      list(status = "ERROR", reason = "Missing exons", n_events = 0L)
    } else {
      gene_strand <- dom_exons$strand[1]

      # Detect events (uses exon structures directly, not UE comparison table)
      events <- detect_events_for_pair(
        dominant_exons = dom_exons,
        comparator_exons = comp_exons,
        gene_id = gene,
        dominant_tid = dom_id,
        comparator_tid = comp_id,
        strand = gene_strand
      )

      n_events <- if (!is.null(events)) nrow(events) else 0L

      # Reconstruct
      comp_for_recon <- comp_exons %>%
        rename(chr = seqnames, transcript_id = isoform_id) %>%
        select(chr, exon_start, exon_end, strand, gene_id, transcript_id)

      gene_ue_recon <- union_exons %>%
        filter(gene_id == gene)

      if (n_events == 0) {
        # No events: comparator should equal dominant
        vr <- verify_transcript(dom_exons, comp_for_recon, gene_strand)
        list(status = if (vr$pass) "PASS" else "FAIL",
             reason = vr$reason, n_events = 0L)
      } else {
        reconstructed <- reconstruct_dominant_v2(comp_for_recon, events, gene_ue_recon)
        if (nrow(reconstructed) == 0) {
          list(status = "ERROR", reason = "Reconstruction produced 0 exons",
               n_events = n_events)
        } else {
          vr <- verify_transcript(dom_exons, reconstructed, gene_strand)
          list(status = if (vr$pass) "PASS" else "FAIL",
               reason = vr$reason, n_events = n_events)
        }
      }
    }
  }, error = function(e) {
    list(status = "ERROR", reason = paste0("Exception: ", e$message), n_events = 0L)
  })

  if (result$status == "PASS") pass_count <- pass_count + 1L
  else if (result$status == "FAIL") fail_count <- fail_count + 1L
  else error_count <- error_count + 1L

  results_list[[i]] <- tibble(
    gene_id = gene,
    dominant_isoform_id = dom_id,
    comparator_isoform_id = comp_id,
    status = result$status,
    reason = result$reason,
    n_events = result$n_events
  )

  # Progress every 500
  if (i %% 500 == 0 || i == nrow(all_pairs)) {
    cat(sprintf("  %d/%d (%.1f%%) | PASS: %d  FAIL: %d  ERROR: %d\n",
                i, nrow(all_pairs),
                100 * i / nrow(all_pairs),
                pass_count, fail_count, error_count))
  }
}

results_df <- bind_rows(results_list)

# ==============================================================================
# 5. Summary
# ==============================================================================

total <- nrow(results_df)

cat("\n")
cat("========================================================================\n")
cat("  GENCODE TEST -- VERIFICATION SUMMARY\n")
cat("========================================================================\n\n")

cat(sprintf("  GTF:    %s\n", gtf_file))
cat(sprintf("  Seed:   %d (random dominant assignment)\n", seed))
cat(sprintf("  Genes:  %d (with 2+ isoforms)\n", length(genes_multi)))
cat(sprintf("  Pairs:  %d\n\n", total))

cat(sprintf("  PASS:   %d/%d (%.1f%%)\n", pass_count, total, 100 * pass_count / total))
cat(sprintf("  FAIL:   %d/%d (%.1f%%)\n", fail_count, total, 100 * fail_count / total))
cat(sprintf("  ERROR:  %d/%d (%.1f%%)\n", error_count, total, 100 * error_count / total))

if (fail_count > 0) {
  cat("\nFailed pairs:\n")
  failures <- results_df %>% filter(status == "FAIL")
  for (j in seq_len(min(nrow(failures), 30))) {
    f <- failures[j, ]
    cat(sprintf("  [FAIL] %s: %s vs %s -- %s (%d events)\n",
                f$gene_id, f$dominant_isoform_id, f$comparator_isoform_id,
                f$reason, f$n_events))
  }
  if (nrow(failures) > 30) cat(sprintf("  ... and %d more\n", nrow(failures) - 30))
}

if (error_count > 0) {
  cat("\nError summary:\n")
  error_reasons <- results_df %>%
    filter(status == "ERROR") %>%
    count(reason, sort = TRUE)
  for (j in seq_len(nrow(error_reasons))) {
    cat(sprintf("  %s: %d\n", error_reasons$reason[j], error_reasons$n[j]))
  }
}

cat("\n========================================================================\n")

# Save
saveRDS(results_df, "data/gencode_test_verification.rds")
cat(sprintf("Saved: data/gencode_test_verification.rds\n\n"))
