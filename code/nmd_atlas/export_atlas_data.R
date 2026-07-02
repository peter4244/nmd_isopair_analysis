#!/usr/bin/env Rscript
# =============================================================================
# export_atlas_data.R — Build the static JSON data files for the NMD lung atlas.
#
# Produces (under public/data/):
#   genes_index.json           — gene_id, hgnc_symbol, n_isoforms, has_nmd_iso
#   quantiles.json             — per-CT CPM quantile breakpoints (gene + isoform)
#   genes/<gene_id>.json       — per-gene shard with all isoforms' structures +
#                                 expression + log2FC + lfsr + GENCODE biotype
#
# Scope: intersection of structures.rds (long-read structures) and the per-CT
# mashr DIE tables (NMD-response statistics). Genes with at least one isoform
# in this intersection are exported.
#
# Cost target: the entire data payload is ~5 MB compressed; per-gene shards are
# typically <10 KB each. Static-host friendly.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

HERE <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grep("--file=", args)]
  if (length(fa) > 0) return(dirname(normalizePath(sub("^--file=", "", fa[1]))))
  normalizePath(getwd())
})()
REPO <- normalizePath(file.path(HERE, "..", ".."))

OUT_DIR       <- file.path(HERE, "public", "data")
GENE_SHARD_DIR <- file.path(OUT_DIR, "genes")
dir.create(GENE_SHARD_DIR, recursive = TRUE, showWarnings = FALSE)

CTS <- c("AT", "DD", "FB", "MV")
ct_lower <- tolower(CTS)
CACHE_RDS <- file.path(REPO, "nmd_orf_model_v5_4ct", "tmp",
                       "ref_orf_ptc_cache_with_nmd_2026.5.12.rds")
STR_RDS  <- NULL  # No longer used — SQANTI3 GTF is the primary structural source.
CDS_RDS  <- file.path(REPO, "results", "isoform_transitions", "Version_6.0",
                      "isopair_wrapper", "data_mashr", "cds.rds")
GTF      <- file.path(REPO, "reference_files",
                      "gencode.v49.primary_assembly.annotation.gtf.gz")
GTF_CACHE <- file.path(HERE, ".gencode_v49_cache.rds")
REF_ATG_RDS <- file.path(REPO, "results", "isoform_transitions", "Version_6.0",
                          "isopair_wrapper", "data_mashr", "analysis_cache",
                          "ref_atg_analysis.rds")
# SQANTI3 sources — the paper's canonical structural + DIE-tested universe
SQANTI_GTF   <- file.path(REPO, "sqanti", "nmd_lungcells", "results",
                          "nmd_lungcells_corrected.sorted.gtf.gz")
SQANTI_CLASS <- file.path(REPO, "sqanti", "nmd_lungcells", "results",
                          "nmd_lungcells_classification.txt")
SQANTI_CACHE <- file.path(HERE, ".sqanti3_cache.rds")
NF_TAG_RX <- "cds_start_NF|cds_end_NF|mRNA_start_NF|mRNA_end_NF"

opt_test <- "--test" %in% commandArgs(trailingOnly = TRUE)

# ── GENCODE v49 transcripts (cached; parsed once) ──
build_gencode_cache <- function(gtf_path, cache_path) {
  cat(sprintf("Parsing GENCODE v49 GTF (this takes a few minutes; cached to %s) ...\n",
              basename(cache_path)))
  suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
  })
  gr <- rtracklayer::import(gtf_path)
  # Transcript-level metadata
  tx_gr <- gr[gr$type == "transcript"]
  tx <- data.table(
    transcript_id = tx_gr$transcript_id,
    gene_id        = tx_gr$gene_id,
    gene_name      = tx_gr$gene_name,
    tx_type        = tx_gr$transcript_type,
    tx_name        = tx_gr$transcript_name,
    chr            = as.character(seqnames(tx_gr)),
    strand         = as.character(strand(tx_gr)),
    tx_start       = start(tx_gr),
    tx_end         = end(tx_gr),
    tags           = vapply(tx_gr$tag,
                             function(x) if (length(x) == 0) NA_character_ else paste(x, collapse = ","),
                             character(1))
  )
  # Per-transcript exon coordinate lists (in genomic order)
  exon_gr <- gr[gr$type == "exon"]
  exon_dt <- data.table(
    transcript_id = exon_gr$transcript_id,
    start = start(exon_gr),
    end   = end(exon_gr)
  )
  setorder(exon_dt, transcript_id, start)
  exon_agg <- exon_dt[, .(exon_starts = list(start),
                           exon_ends   = list(end),
                           n_exons     = .N), by = transcript_id]
  # CDS coordinates per transcript (may be absent)
  cds_gr <- gr[gr$type == "CDS"]
  cds_dt <- data.table(
    transcript_id = cds_gr$transcript_id,
    start = start(cds_gr),
    end   = end(cds_gr)
  )
  cds_agg <- cds_dt[, .(cds_start = min(start),
                         cds_stop  = max(end)),
                     by = transcript_id]
  tx <- merge(tx, exon_agg, by = "transcript_id", all.x = TRUE)
  tx <- merge(tx, cds_agg,  by = "transcript_id", all.x = TRUE)
  saveRDS(tx, cache_path)
  cat(sprintf("Cached %d GENCODE transcripts.\n", nrow(tx)))
  tx
}
load_gencode <- function() {
  if (file.exists(GTF_CACHE)) {
    cat(sprintf("Loading GENCODE cache from %s ...\n", basename(GTF_CACHE)))
    return(readRDS(GTF_CACHE))
  }
  build_gencode_cache(GTF, GTF_CACHE)
}
gencode <- load_gencode()

# ── SQANTI3 structural cache (broadened source for the DIE-tested universe) ──
# Every isoform tested for differential expression (162,800 total) has exon
# coordinates and a structural_category in the SQANTI3 output; we cache that.
build_sqanti_cache <- function(gtf_path, class_path, cache_path) {
  cat(sprintf("Parsing SQANTI3 GTF + classification (cached to %s) ...\n",
              basename(cache_path)))
  suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
  })
  gr <- rtracklayer::import(gtf_path)
  tx_gr <- gr[gr$type == "transcript"]
  tx <- data.table(
    transcript_id = tx_gr$transcript_id,
    gene_id       = tx_gr$gene_id,
    gene_name     = tx_gr$gene_name,
    chr           = as.character(seqnames(tx_gr)),
    strand        = as.character(strand(tx_gr)),
    tx_start      = start(tx_gr),
    tx_end        = end(tx_gr)
  )
  exon_gr <- gr[gr$type == "exon"]
  exon_dt <- data.table(transcript_id = exon_gr$transcript_id,
                         start = start(exon_gr), end = end(exon_gr))
  setorder(exon_dt, transcript_id, start)
  exon_agg <- exon_dt[, .(exon_starts = list(start),
                           exon_ends   = list(end),
                           n_exons     = .N), by = transcript_id]
  tx <- merge(tx, exon_agg, by = "transcript_id", all.x = TRUE)
  # Classification for structural_category + subcategory
  cls <- fread(class_path, sep = "\t", header = TRUE, select = c(
    "isoform", "structural_category", "subcategory"))
  setnames(cls, c("isoform"), c("transcript_id"))
  tx <- merge(tx, cls, by = "transcript_id", all.x = TRUE)
  saveRDS(tx, cache_path)
  cat(sprintf("Cached %d SQANTI3 transcripts.\n", nrow(tx)))
  tx
}
load_sqanti <- function() {
  if (file.exists(SQANTI_CACHE)) {
    cat(sprintf("Loading SQANTI3 cache from %s ...\n", basename(SQANTI_CACHE)))
    return(readRDS(SQANTI_CACHE))
  }
  build_sqanti_cache(SQANTI_GTF, SQANTI_CLASS, SQANTI_CACHE)
}
sqanti <- load_sqanti()

# ── Load core tables ──
cat("Loading CDS (TD2 called ORFs; used as fallback CDS source) ...\n")
cds <- as.data.table(readRDS(CDS_RDS))

cat("Loading per-CT mashr DIE ...\n")
die <- list()
for (i in seq_along(CTS)) {
  ct <- CTS[i]
  f  <- file.path(REPO, "isocall_dge", "mashr",
                   sprintf("nmd_mashr_die_%s_2026.3.10.csv", ct_lower[i]))
  d <- fread(f)
  setnames(d, c("logFC", "adj.P.Val", "nmd_responsive"),
              c(sprintf("logFC_%s", ct), sprintf("adjP_%s", ct),
                sprintf("nmd_resp_%s", ct)))
  die[[ct]] <- d[, .SD, .SDcols = c("txid", "gene_id", "hgnc_symbol",
                                     sprintf("logFC_%s", ct),
                                     sprintf("adjP_%s", ct),
                                     sprintf("nmd_resp_%s", ct))]
}
mashr <- Reduce(function(a, b) merge(a, b, by = c("txid", "gene_id", "hgnc_symbol"),
                                       all = TRUE), die)
setnames(mashr, "txid", "isoform_id")

cat("Loading lfsr (shared, 4-CT) ...\n")
lfsr <- fread(file.path(REPO, "isocall_dge", "mashr",
                         "mashr_isoform_lfsr_2026.3.10.csv"))
setnames(lfsr, "txid", "isoform_id")
for (ct in CTS) setnames(lfsr, sprintf("Smg1i_in_%s", ct), sprintf("lfsr_%s", ct))
lfsr <- lfsr[, .(isoform_id,
                  lfsr_AT, lfsr_DD, lfsr_FB, lfsr_MV)]

cat("Loading expression (CPM per CT) from ref-AUG cache ...\n")
cache <- readRDS(CACHE_RDS)
es <- as.data.table(cache$expr_score)
es <- es[, .(isoform_id, gene_id,
              cpm_AT, cpm_DD, cpm_FB, cpm_MV)]

# ── Reference-AUG projection lookup (paper's canonical CDS methodology) ──
# For each novel comparator paired to a reference isoform: project the
# reference's canonical AUG into the novel's transcript coordinates and walk
# to the first in-frame stop. `ref_atg_analysis.rds` (from the Isopair pipeline)
# gives us ref_atg_genomic, comp_stop_genomic, ref_atg_exonic_in_comp, category.
# We use projections only when they're TRUE-mapped and non-mapping-failed.
cat("Loading reference-AUG projection lookup ...\n")
ra <- readRDS(REF_ATG_RDS)
# Only c2 (NMD-substrate comparators) has the projected stop coord; c4
# (Control comparators) uses their canonical CDS from cds.rds, which is
# not subject to the TD2 PTC-avoidance bias since Controls are not
# NMD substrates.
ref_atg <- as.data.table(ra$c2)[, .(comparator_isoform_id, reference_isoform_id,
                                     ref_atg_genomic, comp_stop_genomic,
                                     ref_atg_exonic_in_comp, category)]
ref_atg <- ref_atg[
  ref_atg_exonic_in_comp == TRUE &
  !is.na(ref_atg_genomic) &
  !is.na(comp_stop_genomic) &
  category != "mapping_failed" &
  !is.na(category)
]
setkey(ref_atg, comparator_isoform_id)
ref_atg <- unique(ref_atg, by = "comparator_isoform_id")
setnames(ref_atg,
          c("ref_atg_genomic", "comp_stop_genomic", "reference_isoform_id"),
          c("refaug_cds_start_g", "refaug_cds_stop_g", "refaug_ref_iso"))
cat(sprintf("  Reference-AUG projections available: %d\n", nrow(ref_atg)))

# ── Stitch the master isoform table ──
cat("Building master isoform table ...\n")
m <- merge(mashr, lfsr, by = "isoform_id", all.x = TRUE)
m <- merge(m, es[, .(isoform_id, cpm_AT, cpm_DD, cpm_FB, cpm_MV)],
            by = "isoform_id", all.x = TRUE)
m <- merge(m, cds[, .(isoform_id, coding_status, cds_start, cds_stop, orf_length)],
            by = "isoform_id", all.x = TRUE)
m <- merge(m, sqanti[, .(isoform_id = transcript_id,
                          chr, strand, n_exons,
                          exon_starts, exon_ends, tx_start, tx_end,
                          sqanti_category = structural_category,
                          sqanti_subcategory = subcategory)],
            by = "isoform_id", all.x = TRUE)
m <- merge(m,
            ref_atg[, .(isoform_id = comparator_isoform_id,
                         refaug_cds_start_g, refaug_cds_stop_g,
                         refaug_ref_iso, refaug_category = category)],
            by = "isoform_id", all.x = TRUE)

# Scope: isoforms that have BOTH a structure AND any expression/mashr signal.
has_struct  <- lengths(m$exon_starts) > 0
has_expr    <- with(m, !(is.na(cpm_AT) & is.na(cpm_DD) & is.na(cpm_FB) & is.na(cpm_MV)))
has_mashr   <- with(m, !(is.na(logFC_AT) & is.na(logFC_DD) &
                          is.na(logFC_FB) & is.na(logFC_MV)))
# Exclude fusion transcripts — they span multiple genes and break the
# one-gene-per-shard mental model. Small share (~1.7% of DIE-tested).
not_fusion  <- is.na(m$sqanti_category) | m$sqanti_category != "fusion"
keep <- has_struct & (has_expr | has_mashr) & not_fusion
cat(sprintf("Isoforms kept: %d / %d (with structure AND expression-or-mashr)\n",
            sum(keep), nrow(m)))
m <- m[keep]

# Test scope: SRSF8 + 50 other random NMD-responsive genes for fast iteration
if (opt_test) {
  sample_genes <- c("ENSG00000263465.5",     # SRSF8 — featured NMD case
                     "ENSG00000128272.20",    # ATF4 — canonical uORF/ISR
                     "ENSG00000112081.19",    # SRSF3 — classic AS-NMD
                     "ENSG00000115875.20",    # SRSF7 — classic AS-NMD
                     "ENSG00000001626.19",    # CFTR
                     "ENSG00000157106.19",    # SMG1
                     "ENSG00000005007.15",    # UPF1
                     "ENSG00000215301.12",    # DDX3X
                     "ENSG00000116044.18",    # NFE2L2
                     "ENSG00000111057.12")    # KRT18
  m <- m[gene_id %in% sample_genes]
  gencode <- gencode[gene_id %in% sample_genes]
  cat(sprintf("TEST scope: %d isoforms across %d genes; %d GENCODE annotated transcripts\n",
              nrow(m), uniqueN(m$gene_id), nrow(gencode)))
}

# ── Gene index (small file for autocomplete + search) ──
# Union of expressed genes and GENCODE-annotated genes so users can also look up
# genes with no detected isoform.
cat("Writing gene index ...\n")
gi_expr <- m[, .(
  n_expr_isoforms = .N,
  hgnc_symbol_expr = data.table::first(na.omit(hgnc_symbol)),
  any_nmd_iso = any(nmd_resp_AT | nmd_resp_DD | nmd_resp_FB | nmd_resp_MV, na.rm = TRUE)
), by = gene_id]
gi_gc <- gencode[, .(
  n_gencode_isoforms = uniqueN(transcript_id),
  hgnc_symbol_gc = data.table::first(na.omit(gene_name))
), by = gene_id]
gi <- merge(gi_expr, gi_gc, by = "gene_id", all = TRUE)
gi[, hgnc_symbol := ifelse(!is.na(hgnc_symbol_expr) & nzchar(hgnc_symbol_expr),
                            hgnc_symbol_expr, hgnc_symbol_gc)]
gi[is.na(hgnc_symbol), hgnc_symbol := ""]
gi[is.na(n_expr_isoforms),    n_expr_isoforms    := 0L]
gi[is.na(n_gencode_isoforms), n_gencode_isoforms := 0L]
gi[is.na(any_nmd_iso), any_nmd_iso := FALSE]
gi[, n_isoforms := pmax(n_expr_isoforms, n_gencode_isoforms)]
gi[, detected   := n_expr_isoforms > 0L]
# Drop internal fields the UI doesn't need — keeps genes_index.json compact.
# `detected` is kept only through the sharding loop; stripped again before the
# genes_index.json write below.
gi <- gi[, .(gene_id, hgnc_symbol, n_isoforms, any_nmd_iso, detected)]
# Note: we don't write the index yet — we'll add the `chr` column after the
# per-chromosome bucketing pass and write the final index there.
cat(sprintf("  -> %d genes\n", nrow(gi)))

# ── Per-CT CPM quantile breakpoints (for "quantile of expression" display) ──
cat("Computing per-CT CPM quantile breakpoints ...\n")
qs <- list()
for (ct in CTS) {
  col <- sprintf("cpm_%s", ct)
  v <- m[[col]]
  v <- v[!is.na(v) & v > 0]
  qs[[ct]] <- as.list(round(quantile(v, probs = seq(0, 1, 0.1)), 3))
}
write_json(qs, file.path(OUT_DIR, "quantiles.json"), auto_unbox = TRUE,
            na = "null")

# ── Per-gene shard JSON ──
# Each shard contains all of the gene's isoforms with structure + expression +
# log2FC + lfsr + flags. Kept compact via short field names and rounded values.
cat("Writing per-gene shards ...\n")
round_or_null <- function(x, d = 3) {
  if (is.null(x) || length(x) == 0) return(NULL)
  if (all(is.na(x))) return(NULL)
  unname(round(x, d))
}
determine_cds_for_expr <- function(r, gc_meta) {
  # Priority for detected isoforms:
  #   1) GENCODE — if the isoform is in GENCODE and no *_NF tags
  #   2) Reference-AUG projection — paper's canonical CDS
  #   3) TD2 (cds.rds) — fallback with visible caveat
  #   4) None
  # GENCODE with any *_NF tag → no CDS block; user sees the biotype + tags
  if (!is.null(gc_meta)) {
    tags_str <- if (is.na(gc_meta$tags)) "" else gc_meta$tags
    has_nf <- nzchar(tags_str) && grepl(NF_TAG_RX, tags_str)
    if (has_nf) {
      return(list(cds_start = NA_integer_, cds_stop = NA_integer_,
                   cds_source = "gencode_nf_placeholder",
                   cds_source_detail = "GENCODE annotates a placeholder CDS but flags the start/end as \"not found\"; no reliable CDS is available."))
    }
    if (!is.na(gc_meta$cds_start) && !is.na(gc_meta$cds_stop)) {
      return(list(cds_start = as.integer(gc_meta$cds_start),
                   cds_stop  = as.integer(gc_meta$cds_stop),
                   cds_source = "gencode_v49",
                   cds_source_detail = "CDS coordinates from GENCODE v49 annotation."))
    }
  }
  if (!is.na(r$refaug_cds_start_g) && !is.na(r$refaug_cds_stop_g)) {
    return(list(cds_start = as.integer(min(r$refaug_cds_start_g, r$refaug_cds_stop_g)),
                 cds_stop  = as.integer(max(r$refaug_cds_start_g, r$refaug_cds_stop_g)),
                 cds_source = "ref_aug_projection",
                 cds_source_detail = sprintf(
                   "Reference-AUG projection: canonical start codon of %s traced into this novel isoform and walked to the first in-frame stop (Isopair; unbiased against PTC-containing ORFs).",
                   r$refaug_ref_iso)))
  }
  if (!is.na(r$cds_start) && !is.na(r$cds_stop) &&
      !is.na(r$coding_status) && r$coding_status == "coding") {
    return(list(cds_start = as.integer(r$cds_start),
                 cds_stop  = as.integer(r$cds_stop),
                 cds_source = "td2",
                 cds_source_detail = "TD2 / SQANTI3-called ORF. This caller is known to avoid PTC-containing ORFs, so the CDS shown for NMD-substrate isoforms may not be the biologically active one."))
  }
  list(cds_start = NA_integer_, cds_stop = NA_integer_,
        cds_source = "none",
        cds_source_detail = NA_character_)
}

iso_record_from_expr_row <- function(r, gc_meta = NULL) {
  cds <- determine_cds_for_expr(r, gc_meta)
  list(
    id      = r$isoform_id,
    n_exons = if (is.na(r$n_exons)) NULL else as.integer(r$n_exons),
    chr     = if (is.na(r$chr))    NULL else r$chr,
    strand  = if (is.na(r$strand)) NULL else r$strand,
    tx_start = if (is.na(r$tx_start)) NULL else as.integer(r$tx_start),
    tx_end   = if (is.na(r$tx_end))   NULL else as.integer(r$tx_end),
    exon_starts = if (is.null(r$exon_starts[[1]])) NULL else as.integer(r$exon_starts[[1]]),
    exon_ends   = if (is.null(r$exon_ends[[1]]))   NULL else as.integer(r$exon_ends[[1]]),
    cds_start = if (is.na(cds$cds_start)) NULL else cds$cds_start,
    cds_stop  = if (is.na(cds$cds_stop))  NULL else cds$cds_stop,
    cds_source        = cds$cds_source,
    cds_source_detail = if (is.na(cds$cds_source_detail)) NULL else cds$cds_source_detail,
    orf_length    = if (is.na(r$orf_length))    NULL else as.integer(r$orf_length),
    coding_status = if (is.na(r$coding_status)) NULL else r$coding_status,
    gencode_biotype = if (is.null(gc_meta)) NULL else gc_meta$tx_type,
    gencode_name    = if (is.null(gc_meta)) NULL else gc_meta$tx_name,
    tags            = if (is.null(gc_meta)) NULL else gc_meta$tags,
    sqanti_category = if (is.na(r$sqanti_category)) NULL else r$sqanti_category,
    sqanti_subcategory = if (is.na(r$sqanti_subcategory)) NULL else r$sqanti_subcategory,
    cpm = list(AT = round_or_null(r$cpm_AT), DD = round_or_null(r$cpm_DD),
                FB = round_or_null(r$cpm_FB), MV = round_or_null(r$cpm_MV)),
    logfc = list(AT = round_or_null(r$logFC_AT), DD = round_or_null(r$logFC_DD),
                  FB = round_or_null(r$logFC_FB), MV = round_or_null(r$logFC_MV)),
    adjP  = list(AT = round_or_null(r$adjP_AT, 6), DD = round_or_null(r$adjP_DD, 6),
                  FB = round_or_null(r$adjP_FB, 6), MV = round_or_null(r$adjP_MV, 6)),
    lfsr  = list(AT = round_or_null(r$lfsr_AT, 4), DD = round_or_null(r$lfsr_DD, 4),
                  FB = round_or_null(r$lfsr_FB, 4), MV = round_or_null(r$lfsr_MV, 4)),
    nmd_responsive = list(
      AT = isTRUE(r$nmd_resp_AT), DD = isTRUE(r$nmd_resp_DD),
      FB = isTRUE(r$nmd_resp_FB), MV = isTRUE(r$nmd_resp_MV))
  )
}
iso_record_from_gencode <- function(g) {
  tags_str <- if (is.na(g$tags)) "" else g$tags
  has_nf   <- nzchar(tags_str) && grepl(NF_TAG_RX, tags_str)
  if (has_nf) {
    cds_s <- NA_integer_; cds_e <- NA_integer_
    src <- "gencode_nf_placeholder"
    src_det <- "GENCODE annotates a placeholder CDS but flags the start/end as \"not found\"; no reliable CDS is available."
  } else if (!is.na(g$cds_start) && !is.na(g$cds_stop)) {
    cds_s <- as.integer(g$cds_start); cds_e <- as.integer(g$cds_stop)
    src <- "gencode_v49"
    src_det <- "CDS coordinates from GENCODE v49 annotation."
  } else {
    cds_s <- NA_integer_; cds_e <- NA_integer_
    src <- "none"; src_det <- NA_character_
  }
  list(
    id      = g$transcript_id,
    n_exons = if (is.na(g$n_exons)) NULL else as.integer(g$n_exons),
    chr     = g$chr, strand = g$strand,
    tx_start = as.integer(g$tx_start), tx_end = as.integer(g$tx_end),
    exon_starts = if (is.null(g$exon_starts[[1]])) NULL else as.integer(g$exon_starts[[1]]),
    exon_ends   = if (is.null(g$exon_ends[[1]]))   NULL else as.integer(g$exon_ends[[1]]),
    cds_start = if (is.na(cds_s)) NULL else cds_s,
    cds_stop  = if (is.na(cds_e)) NULL else cds_e,
    cds_source        = src,
    cds_source_detail = if (is.na(src_det)) NULL else src_det,
    orf_length    = NULL,
    coding_status = NULL,
    gencode_biotype = g$tx_type,
    gencode_name    = g$tx_name,
    tags            = g$tags,
    cpm   = list(AT = NULL, DD = NULL, FB = NULL, MV = NULL),
    logfc = list(AT = NULL, DD = NULL, FB = NULL, MV = NULL),
    adjP  = list(AT = NULL, DD = NULL, FB = NULL, MV = NULL),
    lfsr  = list(AT = NULL, DD = NULL, FB = NULL, MV = NULL),
    nmd_responsive = list(AT = FALSE, DD = FALSE, FB = FALSE, MV = FALSE),
    gencode_only = TRUE
  )
}
# Returns the per-gene shard list (does NOT write to disk — the caller
# aggregates all shards for a chromosome, then writes one file per chr).
# Also returns the chromosome so the aggregator can bucket correctly.
build_gene_shard <- function(gid, gene_rows, gencode_rows) {
  # Expressed / detected isoforms — annotate with GENCODE biotype when the id matches
  gc_by_id <- setNames(split(gencode_rows, gencode_rows$transcript_id),
                        gencode_rows$transcript_id)
  isos_expr <- lapply(seq_len(nrow(gene_rows)), function(i) {
    r <- gene_rows[i]
    gc <- gc_by_id[[r$isoform_id]]
    iso_record_from_expr_row(r, if (is.null(gc)) NULL else gc[1])
  })
  # GENCODE transcripts NOT in our expressed set → include with gencode_only = TRUE
  expr_ids <- gene_rows$isoform_id
  extra <- gencode_rows[!(transcript_id %in% expr_ids)]
  isos_extra <- lapply(seq_len(nrow(extra)), function(i) iso_record_from_gencode(extra[i]))
  hgnc <- data.table::first(na.omit(c(gene_rows$hgnc_symbol,
                                       if (nrow(gencode_rows) > 0) gencode_rows$gene_name)))
  # Pick a chromosome for this gene (nearly always uniform across its isoforms)
  chrs <- na.omit(c(gene_rows$chr, if (nrow(gencode_rows) > 0) gencode_rows$chr))
  gene_chr <- if (length(chrs) > 0) as.character(chrs[1]) else NA_character_
  out <- list(
    gene_id     = gid,
    hgnc_symbol = hgnc,
    chr         = gene_chr,
    isoforms    = c(isos_expr, isos_extra)
  )
  list(shard = out, chr = gene_chr)
}

# Set keys once — orders of magnitude faster than repeated equality filtering
setkey(m, gene_id)
setkey(gencode, gene_id)
all_gene_ids <- union(unique(m$gene_id), unique(gencode$gene_id))
n_built <- 0
n_skipped_undetected <- 0
gene_chr_map <- character(length(all_gene_ids)); names(gene_chr_map) <- all_gene_ids
detected_gids <- gi[detected == TRUE, gene_id]
detected_set  <- new.env(hash = TRUE, size = length(detected_gids))
for (g in detected_gids) detected_set[[g]] <- TRUE
# Write one JSON per gene. Each click is a single tiny fetch (~5-15 KB raw,
# ~1-3 KB gzipped) instead of pulling the whole chromosome. gene_chr_map is
# kept only to populate the index's `chr` field so the "not detected" state
# can name a location for GENCODE-only pseudogenes.
for (gid in all_gene_ids) {
  gene_rows    <- m[.(gid), nomatch = 0L]
  gencode_rows <- gencode[.(gid), nomatch = 0L]
  if (nrow(gene_rows) == 0 && nrow(gencode_rows) == 0) next
  if (is.null(detected_set[[gid]])) {
    n_skipped_undetected <- n_skipped_undetected + 1
    next
  }
  built <- build_gene_shard(gid, gene_rows, gencode_rows)
  chr_key <- if (is.na(built$chr) || is.null(built$chr)) "unknown" else built$chr
  write_json(built$shard,
             file.path(GENE_SHARD_DIR, sprintf("%s.json", gid)),
             auto_unbox = TRUE, na = "null", null = "null")
  gene_chr_map[gid] <- chr_key
  n_built <- n_built + 1
  if (n_built %% 5000 == 0) cat(sprintf("  ...%d gene shards written\n", n_built))
}
cat(sprintf("Skipped %d undetected genes (kept in index; not written to gene shards).\n",
            n_skipped_undetected))
cat(sprintf("Wrote %d per-gene shards to %s\n", n_built, GENE_SHARD_DIR))


# Emit gene_id → chr lookup into the gene index. Not used for fetching (each
# gene is its own file) — kept so the UI can display the location and the
# "not detected" placeholder can name a chromosome.
gene_chr_dt <- data.table(gene_id = names(gene_chr_map), chr = unname(gene_chr_map))
# Fill in chr for undetected GENCODE-only genes from the GENCODE cache.
gencode_chr <- gencode[!is.na(chr), .(gc_chr = chr[1]), by = gene_id]
gene_chr_dt <- merge(gene_chr_dt, gencode_chr, by = "gene_id", all.x = TRUE)
gene_chr_dt[is.na(chr) | chr == "" | chr == "unknown",
            chr := ifelse(is.na(gc_chr), "", gc_chr)]
gene_chr_dt[, gc_chr := NULL]
gi <- merge(gi, gene_chr_dt, by = "gene_id", all.x = TRUE)
# `detected` is emitted so the UI can skip the fetch entirely for GENCODE-only
# genes and render the "not detected" placeholder from the index alone.
write_json(gi, file.path(OUT_DIR, "genes_index.json"),
            auto_unbox = TRUE, dataframe = "rows", na = "null", null = "null")

# ── Manifest ──
manifest <- list(
  generated_at = Sys.time(),
  n_genes      = n_built,
  n_isoforms   = nrow(m),
  n_gene_shards = n_built,
  cell_types   = CTS,
  data_version = "2026.7.2"
)
write_json(manifest, file.path(OUT_DIR, "manifest.json"), auto_unbox = TRUE)
cat("Done.\n")
