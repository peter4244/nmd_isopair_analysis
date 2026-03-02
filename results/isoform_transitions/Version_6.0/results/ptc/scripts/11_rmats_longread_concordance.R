#!/usr/bin/env Rscript
# 11_rmats_longread_concordance.R
#
# rMATS–PacBio concordance analysis: systematic cross-platform comparison
# to understand why rMATS detects frame-disruption enrichment but PacBio
# shows no PTC enrichment in NMD-responsive isoforms.
#
# Tests two hypotheses:
#   H1 (Reconcilable): Frame-disrupting events at junction level don't
#       necessarily produce PTCs in full-length transcripts
#   H2 (Reconstruction error): PacBio isoform catalog doesn't capture
#       the relevant splice variants that rMATS detects
#
# Event types: SE, A3SS, A5SS
#
# Input:
#   ../../../results/rds/allrmats_jc_2025.8.20.rds  (maser objects)
#   data/isoform_structures.rds                      (exon coordinates)
#   data/isoform_cds_metadata.rds                    (CDS boundaries)
#   results/ptc/results/ptc_status.rds               (PTC status)
#   ../../../rds/dge_isoform_2026.2.7.rds            (DGEList, raw counts)
#   ../../../longread_dge/nmd_dge_{ct}_2026.1.18.csv (DE results)
#
# Output (results/ptc/results/):
#   rmats_junction_concordance.tsv
#   rmats_ptc_verification.tsv
#   rmats_psi_concordance.tsv
#   rmats_discrepancy_waterfall.tsv
#   rmats_concordance_summary.tsv
#
# Run from: results/isoform_transitions/Version_6.0/

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(maser)
  library(GenomicRanges)
})

# ─── Parse arguments ─────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
source_type <- "oarfish"  # default
if ("--source" %in% args) {
  idx <- which(args == "--source")
  if (idx < length(args)) source_type <- args[idx + 1]
}
stopifnot(source_type %in% c("oarfish", "isocall"))

source_label <- if (source_type == "isocall") "isocall" else "PacBio"

cat("\n")
cat("==================================================================\n")
cat(sprintf("   11: rMATS–%s Concordance Analysis\n", source_label))
cat("==================================================================\n\n")

base_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0"
setwd(base_dir)

out_dir <- if (source_type == "isocall") "results/ptc/results/isocall" else "results/ptc/results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Cell types to analyze: DD (primary), MV, AT (validation)
FOCUS_CTS <- c("DD", "MV", "AT")

# rMATS cell type name → pipeline cell type name
CT_MAP <- c(AT2 = "AT", DD = "DD", DD_ALI = "DD_ALI",
            DO_ALI = "DO", FB = "FB", MV = "MV")

# Pipeline cell type → DE file code
if (source_type == "isocall") {
  CT_DE_MAP <- c(AT = "at", DD = "dd", DD_ALI = "ddali",
                 DO = "do", FB = "fb", MV = "mv")
} else {
  CT_DE_MAP <- c(AT = "at2", DD = "dd", DD_ALI = "ddali",
                 DO = "doali", FB = "fb", MV = "mv")
}

# Significance thresholds (same as Script 06)
FDR_THRESH <- 0.05
DPSI_THRESH <- 0.05

# Event types for this analysis
EVENT_TYPES <- c("SE", "A5SS", "A3SS")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Load Data
# ═══════════════════════════════════════════════════════════════════════════════

cat("── Section 1: Loading data ────────────────────────────────────────\n\n")

# --- 1a: rMATS maser objects ---
rmats_jc <- readRDS("../../../results/rds/allrmats_jc_2025.8.20.rds")
cat(sprintf("  rMATS cell types: %s\n", paste(names(rmats_jc), collapse = ", ")))

obj0 <- rmats_jc[[1]]
cond_labels <- slot(obj0, "conditions")
cat(sprintf("  Conditions: cond1 = %s, cond2 = %s\n", cond_labels[1], cond_labels[2]))

# --- 1b: Isoform structures ---
iso_str_path <- if (source_type == "isocall") "data/isocall/isoform_structures.rds" else "data/isoform_structures.rds"
iso_str <- as.data.table(readRDS(iso_str_path))
cat(sprintf("  Isoform structures: %s isoforms\n",
            format(nrow(iso_str), big.mark = ",")))

# --- 1c: PTC status ---
ptc_path <- if (source_type == "isocall") "results/ptc/results/isocall/ptc_status.rds" else "results/ptc/results/ptc_status.rds"
ptc_status <- as.data.table(readRDS(ptc_path))
cat(sprintf("  PTC status: %s isoforms\n",
            format(nrow(ptc_status), big.mark = ",")))

# --- 1d: DGEList for raw counts ---
if (source_type == "isocall") {
  dge <- readRDS("data/isocall/dge_isocall_unfiltered.rds")
} else {
  dge <- readRDS("../../../rds/dge_isoform_2026.2.7.rds")
}
cat(sprintf("  DGEList: %s isoforms × %d samples\n",
            format(nrow(dge$counts), big.mark = ","), ncol(dge$counts)))

# Build iso → gene map from DGEList genes table
if (source_type == "isocall") {
  iso_gene_map <- as.data.table(dge$genes)[, .(isoform_id = transcript_id,
                                                gene_id_ensg = sub("\\.[0-9]+$", "", gene_id))]
} else {
  iso_gene_map <- as.data.table(dge$genes)[, .(isoform_id = txid,
                                                gene_id_ensg = gene_id_ens115_sqanti)]
}
cat(sprintf("  iso_gene_map: %s isoforms\n",
            format(nrow(iso_gene_map), big.mark = ",")))

# --- 1e: DE files ---
if (source_type == "isocall") {
  de_files <- list(
    AT = "../../../isocall_dge/nmd_isocall_dge_at_2026.3.1.csv",
    DD = "../../../isocall_dge/nmd_isocall_dge_dd_2026.3.1.csv",
    MV = "../../../isocall_dge/nmd_isocall_dge_mv_2026.3.1.csv"
  )
} else {
  de_files <- list(
    AT = "../../../longread_dge/nmd_dge_at2_2026.1.18.csv",
    DD = "../../../longread_dge/nmd_dge_dd_2026.1.18.csv",
    MV = "../../../longread_dge/nmd_dge_mv_2026.1.18.csv"
  )
}

de_results <- list()
for (ct in names(de_files)) {
  de_results[[ct]] <- fread(de_files[[ct]])
  # Normalize column name: isocall uses transcript_id, oarfish uses txid
  if (source_type == "isocall" && "transcript_id" %in% names(de_results[[ct]])) {
    setnames(de_results[[ct]], "transcript_id", "txid")
  }
  cat(sprintf("  DE %s: %s transcripts\n", ct,
              format(nrow(de_results[[ct]]), big.mark = ",")))
}

# --- 1f: CDS metadata for gene-level CDS envelope ---
cds_meta_path <- if (source_type == "isocall") "data/isocall/isoform_cds_metadata.rds" else "data/isoform_cds_metadata.rds"
cds_meta <- readRDS(cds_meta_path)
gencode_coding <- cds_meta %>%
  filter(coding_status == "coding", grepl("^ENST", isoform_id)) %>%
  inner_join(
    iso_str[, .(isoform_id, gene_id, seqnames, strand)] %>% as_tibble(),
    by = "isoform_id"
  )

# Strip gene_id version for isocall so it matches rMATS unversioned gene IDs
if (source_type == "isocall") {
  gencode_coding$gene_id <- sub("\\.[0-9]+$", "", gencode_coding$gene_id)
}

gene_cds <- gencode_coding %>%
  group_by(gene_id, seqnames) %>%
  summarise(
    cds_min = min(cds_start),
    cds_max = max(cds_stop),
    .groups = "drop"
  )
cat(sprintf("  Gene-level CDS envelopes: %s genes\n",
            format(nrow(gene_cds), big.mark = ",")))


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Extract Splice Junctions
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Section 2: Extract splice junctions ─────────────────────────────\n\n")

# --- 2a: Extract rMATS events with coordinates and stats ---

extract_rmats_events <- function(maser_obj, cell_type) {
  results <- list()
  for (etype in EVENT_TYPES) {
    events <- slot(maser_obj, paste0(etype, "_events")) %>% as_tibble()
    stats  <- slot(maser_obj, paste0(etype, "_stats"))  %>% as_tibble()
    gr     <- slot(maser_obj, paste0(etype, "_gr"))

    merged <- inner_join(events, stats, by = "ID")

    if (etype == "SE") {
      target <- gr[["exon_target"]]
      up_gr  <- gr[["exon_upstream"]]
      dn_gr  <- gr[["exon_downstream"]]
      merged$event_chr        <- as.character(seqnames(target))
      merged$event_start      <- start(target)
      merged$event_end        <- end(target)
      merged$event_width      <- width(target)
      merged$upstream_start   <- start(up_gr)
      merged$upstream_end     <- end(up_gr)
      merged$downstream_start <- start(dn_gr)
      merged$downstream_end   <- end(dn_gr)
    } else if (etype %in% c("A5SS", "A3SS")) {
      long_gr  <- gr[["exon_long"]]
      short_gr <- gr[["exon_short"]]
      flank_gr <- gr[["exon_flanking"]]
      merged$event_chr   <- as.character(seqnames(long_gr))
      merged$long_start  <- start(long_gr)
      merged$long_end    <- end(long_gr)
      merged$long_width  <- width(long_gr)
      merged$short_start <- start(short_gr)
      merged$short_end   <- end(short_gr)
      merged$short_width <- width(short_gr)
      merged$flank_start <- start(flank_gr)
      merged$flank_end   <- end(flank_gr)
      merged$event_start <- pmin(start(long_gr), start(short_gr))
      merged$event_end   <- pmax(end(long_gr), end(short_gr))
      merged$width_diff  <- abs(merged$long_width - merged$short_width)
    }

    merged$event_type <- etype
    merged$cell_type  <- CT_MAP[cell_type]
    merged$gene_id_unversioned <- sub("[.].*", "", merged$GeneID)

    results[[etype]] <- merged
  }
  bind_rows(results)
}

all_events <- map2_dfr(rmats_jc[names(CT_MAP) %in% names(rmats_jc)],
                        names(rmats_jc[names(CT_MAP) %in% names(rmats_jc)]),
                        extract_rmats_events)

# Add CDS overlap and frame disruption (same logic as Script 06)
all_events <- all_events %>%
  left_join(gene_cds, by = c("gene_id_unversioned" = "gene_id",
                              "event_chr" = "seqnames")) %>%
  mutate(
    overlaps_cds = !is.na(cds_min) &
                   event_start <= cds_max &
                   event_end >= cds_min,
    frame_disrupting = case_when(
      event_type == "SE" ~ (event_width %% 3) != 0,
      event_type %in% c("A5SS", "A3SS") ~ (width_diff %% 3) != 0,
      TRUE ~ NA
    ),
    frame_class = case_when(
      !overlaps_cds                   ~ "non_cds",
      overlaps_cds & frame_disrupting ~ "cds_frame_disrupting",
      overlaps_cds & !frame_disrupting ~ "cds_frame_preserving",
      TRUE                            ~ "unclassified"
    ),
    is_significant = !is.na(FDR) & FDR < FDR_THRESH &
                     !is.na(IncLevelDifference) & abs(IncLevelDifference) >= DPSI_THRESH
  )

cat(sprintf("  Total events extracted: %s\n",
            format(nrow(all_events), big.mark = ",")))

# Filter to focus cell types + significant frame-disrupting CDS events
sig_fd <- all_events %>%
  filter(cell_type %in% FOCUS_CTS,
         frame_class == "cds_frame_disrupting",
         is_significant)

cat(sprintf("  Significant frame-disrupting CDS events in focus cell types: %s\n",
            format(nrow(sig_fd), big.mark = ",")))
cat("  By cell type and event type:\n")
sig_fd %>% count(cell_type, event_type) %>%
  pivot_wider(names_from = event_type, values_from = n, values_fill = 0) %>%
  print(n = Inf)

# --- 2b: Extract rMATS junctions per event ---
# SE: 3 junctions (inclusion_1, inclusion_2, skipping)
# A3SS: 2 junctions (long_form, short_form)
# A5SS: 2 junctions (long_form, short_form)
#
# STRAND HANDLING: maser names exons by transcript order (upstream/downstream),
# but on - strand genes these are reversed in genomic space. We use pmin/pmax
# to always produce junctions as (left_coord, right_coord) in genomic order,
# matching PacBio junction extraction which is always genomic-order.

extract_rmats_junctions <- function(events_df) {
  jn_list <- list()

  se <- events_df %>% filter(event_type == "SE")
  if (nrow(se) > 0) {
    # SE: three exons in transcript order: upstream → target → downstream
    # On - strand, genomic order is reversed: downstream < target < upstream
    # Use pmin/pmax for strand-independent genomic-order junctions:
    #   left_flank_end  = pmin(upstream_end, downstream_end)
    #   right_flank_start = pmax(upstream_start, downstream_start)
    jn_list[["SE_incl1"]] <- se %>%
      transmute(event_key = paste(cell_type, event_type, ID, sep = ":"),
                cell_type, event_type, ID, gene_id_unversioned,
                junction_type = "inclusion_1",
                donor = pmin(upstream_end, downstream_end),
                acceptor = event_start)
    jn_list[["SE_incl2"]] <- se %>%
      transmute(event_key = paste(cell_type, event_type, ID, sep = ":"),
                cell_type, event_type, ID, gene_id_unversioned,
                junction_type = "inclusion_2",
                donor = event_end,
                acceptor = pmax(upstream_start, downstream_start))
    jn_list[["SE_skip"]] <- se %>%
      transmute(event_key = paste(cell_type, event_type, ID, sep = ":"),
                cell_type, event_type, ID, gene_id_unversioned,
                junction_type = "skipping",
                donor = pmin(upstream_end, downstream_end),
                acceptor = pmax(upstream_start, downstream_start))
  }

  # A3SS/A5SS: two adjacent exons (long/short and flanking) separated by intron.
  # On + strand flanking may be left; on - strand flanking may be right.
  # Junction = (end of left exon, start of right exon) in genomic order:
  #   donor   = pmin(flank_end, long_end)    = end of leftward exon
  #   acceptor = pmax(flank_start, long_start) = start of rightward exon
  a3 <- events_df %>% filter(event_type == "A3SS")
  if (nrow(a3) > 0) {
    jn_list[["A3SS_long"]] <- a3 %>%
      transmute(event_key = paste(cell_type, event_type, ID, sep = ":"),
                cell_type, event_type, ID, gene_id_unversioned,
                junction_type = "long_form",
                donor = pmin(flank_end, long_end),
                acceptor = pmax(flank_start, long_start))
    jn_list[["A3SS_short"]] <- a3 %>%
      transmute(event_key = paste(cell_type, event_type, ID, sep = ":"),
                cell_type, event_type, ID, gene_id_unversioned,
                junction_type = "short_form",
                donor = pmin(flank_end, short_end),
                acceptor = pmax(flank_start, short_start))
  }

  a5 <- events_df %>% filter(event_type == "A5SS")
  if (nrow(a5) > 0) {
    jn_list[["A5SS_long"]] <- a5 %>%
      transmute(event_key = paste(cell_type, event_type, ID, sep = ":"),
                cell_type, event_type, ID, gene_id_unversioned,
                junction_type = "long_form",
                donor = pmin(flank_end, long_end),
                acceptor = pmax(flank_start, long_start))
    jn_list[["A5SS_short"]] <- a5 %>%
      transmute(event_key = paste(cell_type, event_type, ID, sep = ":"),
                cell_type, event_type, ID, gene_id_unversioned,
                junction_type = "short_form",
                donor = pmin(flank_end, short_end),
                acceptor = pmax(flank_start, short_start))
  }

  bind_rows(jn_list)
}

rmats_junctions <- extract_rmats_junctions(sig_fd)
cat(sprintf("\n  rMATS junctions from sig FD events: %s\n",
            format(nrow(rmats_junctions), big.mark = ",")))

# --- 2c: Extract isoform junctions from isoform_structures ---
# For each isoform: junction = (exon_ends[i], exon_starts[i+1])
# Includes both PacBio and GENCODE isoforms at genes with rMATS events

# Get genes with rMATS events
rmats_genes <- unique(sig_fd$gene_id_unversioned)

# Map isoforms (PacBio + GENCODE) to ENSG gene IDs via DGEList
isos_at_genes <- iso_gene_map[gene_id_ensg %in% rmats_genes]
cat(sprintf("  %s isoforms at rMATS event genes: %s (via DGEList gene map)\n", source_label,
            format(nrow(isos_at_genes), big.mark = ",")))

# Pre-build named lookup: isoform_id → gene_id_ensg (already unversioned via iso_gene_map)
iso2gene <- setNames(isos_at_genes$gene_id_ensg, isos_at_genes$isoform_id)

# Extract junctions from isoform_structures for these isoforms
lr_isos <- iso_str[isoform_id %in% isos_at_genes$isoform_id]

lr_jn_list <- vector("list", nrow(lr_isos))
for (i in seq_len(nrow(lr_isos))) {
  iso <- lr_isos[i]
  starts <- iso$exon_starts[[1]]
  ends   <- iso$exon_ends[[1]]
  n_ex   <- length(starts)
  if (n_ex < 2) next

  lr_jn_list[[i]] <- data.table(
    isoform_id = iso$isoform_id,
    gene_id_ensg = iso2gene[iso$isoform_id],
    donor = ends[-n_ex],
    acceptor = starts[-1]
  )
}
pb_junctions <- rbindlist(lr_jn_list[!sapply(lr_jn_list, is.null)])

cat(sprintf("  %s junctions at event genes: %s\n", source_label,
            format(nrow(pb_junctions), big.mark = ",")))

# Key for fast lookup
setkey(pb_junctions, gene_id_ensg, donor, acceptor)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Junction Concordance (Tests H2)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Section 3: Junction concordance ─────────────────────────────────\n\n")

# For each rMATS event, check if PacBio has matching junctions.
# SE: need inclusion (both jn1+jn2 in same isoform) AND skipping
# A3SS/A5SS: need long_form AND short_form junctions

classify_event_junctions <- function(event_key, event_type, gene_id,
                                      junctions_df, pb_jn_dt) {
  ev_jn <- junctions_df[junctions_df$event_key == event_key, ]

  # Check if gene has any PacBio isoforms
  gene_pb <- pb_jn_dt[gene_id_ensg == gene_id]
  if (nrow(gene_pb) == 0) {
    return(list(category = "no_pb_gene",
                inclusion_isos = character(0),
                skipping_isos = character(0),
                long_isos = character(0),
                short_isos = character(0)))
  }

  if (event_type == "SE") {
    # Inclusion: isoform has BOTH jn1 (upstream→target) and jn2 (target→downstream)
    jn1 <- ev_jn[ev_jn$junction_type == "inclusion_1", ]
    jn2 <- ev_jn[ev_jn$junction_type == "inclusion_2", ]
    jn_skip <- ev_jn[ev_jn$junction_type == "skipping", ]

    incl1_isos <- gene_pb[donor == jn1$donor[1] & acceptor == jn1$acceptor[1], isoform_id]
    incl2_isos <- gene_pb[donor == jn2$donor[1] & acceptor == jn2$acceptor[1], isoform_id]
    inclusion_isos <- intersect(incl1_isos, incl2_isos)

    skip_isos <- gene_pb[donor == jn_skip$donor[1] & acceptor == jn_skip$acceptor[1], isoform_id]

    has_incl <- length(inclusion_isos) > 0
    has_skip <- length(skip_isos) > 0

    category <- case_when(
      has_incl & has_skip ~ "both",
      has_incl & !has_skip ~ "inclusion_only",
      !has_incl & has_skip ~ "skipping_only",
      TRUE ~ "neither"
    )

    return(list(category = category,
                inclusion_isos = inclusion_isos,
                skipping_isos = skip_isos,
                long_isos = character(0),
                short_isos = character(0)))

  } else {
    # A3SS or A5SS: check for long_form and short_form junctions
    jn_long  <- ev_jn[ev_jn$junction_type == "long_form", ]
    jn_short <- ev_jn[ev_jn$junction_type == "short_form", ]

    long_isos  <- gene_pb[donor == jn_long$donor[1] & acceptor == jn_long$acceptor[1], isoform_id]
    short_isos <- gene_pb[donor == jn_short$donor[1] & acceptor == jn_short$acceptor[1], isoform_id]

    has_long  <- length(long_isos) > 0
    has_short <- length(short_isos) > 0

    category <- case_when(
      has_long & has_short ~ "both",
      has_long & !has_short ~ "inclusion_only",
      !has_long & has_short ~ "skipping_only",
      TRUE ~ "neither"
    )

    return(list(category = category,
                inclusion_isos = character(0),
                skipping_isos = character(0),
                long_isos = long_isos,
                short_isos = short_isos))
  }
}

# Run junction classification for all significant frame-disrupting events
unique_events <- sig_fd %>%
  distinct(cell_type, event_type, ID, gene_id_unversioned,
           .keep_all = TRUE) %>%
  mutate(event_key = paste(cell_type, event_type, ID, sep = ":"))

cat(sprintf("  Classifying %d events...\n", nrow(unique_events)))

concordance_results <- vector("list", nrow(unique_events))
for (i in seq_len(nrow(unique_events))) {
  ev <- unique_events[i, ]
  res <- classify_event_junctions(
    ev$event_key, ev$event_type, ev$gene_id_unversioned,
    rmats_junctions, pb_junctions
  )

  concordance_results[[i]] <- tibble(
    event_key = ev$event_key,
    cell_type = ev$cell_type,
    event_type = ev$event_type,
    event_id = ev$ID,
    gene_id = ev$gene_id_unversioned,
    gene_symbol = ev$geneSymbol,
    dPSI = ev$IncLevelDifference,
    FDR = ev$FDR,
    junction_category = res$category,
    n_inclusion_isos = if (ev$event_type == "SE") length(res$inclusion_isos) else length(res$long_isos),
    n_skipping_isos = if (ev$event_type == "SE") length(res$skipping_isos) else length(res$short_isos),
    inclusion_isos = paste(if (ev$event_type == "SE") res$inclusion_isos else res$long_isos, collapse = ";"),
    skipping_isos = paste(if (ev$event_type == "SE") res$skipping_isos else res$short_isos, collapse = ";")
  )

  if (i %% 500 == 0) cat(sprintf("    %d / %d\n", i, nrow(unique_events)))
}

junction_concordance <- bind_rows(concordance_results)

cat("\n  Junction concordance summary:\n")
junction_concordance %>%
  count(cell_type, event_type, junction_category) %>%
  pivot_wider(names_from = junction_category, values_from = n, values_fill = 0) %>%
  print(n = Inf)

# Overall concordance rate
both_rate <- junction_concordance %>%
  summarise(
    total = n(),
    both = sum(junction_category == "both"),
    pct_both = round(100 * both / total, 1)
  )
cat(sprintf("\n  Overall: %d / %d (%.1f%%) have both forms in %s\n",
            both_rate$both, both_rate$total, both_rate$pct_both, source_label))

# --- Near-miss diagnostic ---
cat("\n  Near-miss diagnostic (±1-2bp)...\n")

check_near_miss <- function(event_key, event_type, gene_id,
                             junctions_df, pb_jn_dt, tolerance = 2) {
  ev_jn <- junctions_df[junctions_df$event_key == event_key, ]
  gene_pb <- pb_jn_dt[gene_id_ensg == gene_id]
  if (nrow(gene_pb) == 0) return(FALSE)

  if (event_type == "SE") {
    jn1 <- ev_jn[ev_jn$junction_type == "inclusion_1", ]
    jn2 <- ev_jn[ev_jn$junction_type == "inclusion_2", ]
    jn_skip <- ev_jn[ev_jn$junction_type == "skipping", ]

    near1 <- gene_pb[abs(donor - jn1$donor[1]) <= tolerance &
                      abs(acceptor - jn1$acceptor[1]) <= tolerance]
    near2 <- gene_pb[abs(donor - jn2$donor[1]) <= tolerance &
                      abs(acceptor - jn2$acceptor[1]) <= tolerance]
    near_skip <- gene_pb[abs(donor - jn_skip$donor[1]) <= tolerance &
                          abs(acceptor - jn_skip$acceptor[1]) <= tolerance]

    has_near_incl <- length(intersect(near1$isoform_id, near2$isoform_id)) > 0
    has_near_skip <- nrow(near_skip) > 0
    return(has_near_incl & has_near_skip)
  } else {
    jn_long  <- ev_jn[ev_jn$junction_type == "long_form", ]
    jn_short <- ev_jn[ev_jn$junction_type == "short_form", ]

    near_long  <- gene_pb[abs(donor - jn_long$donor[1]) <= tolerance &
                           abs(acceptor - jn_long$acceptor[1]) <= tolerance]
    near_short <- gene_pb[abs(donor - jn_short$donor[1]) <= tolerance &
                           abs(acceptor - jn_short$acceptor[1]) <= tolerance]
    return(nrow(near_long) > 0 & nrow(near_short) > 0)
  }
}

# Check near-misses only for events that failed exact match
missed_events <- unique_events %>%
  inner_join(
    junction_concordance %>% filter(junction_category != "both") %>%
      select(event_key),
    by = "event_key"
  )

if (nrow(missed_events) > 0) {
  near_miss_flags <- logical(nrow(missed_events))
  for (i in seq_len(nrow(missed_events))) {
    ev <- missed_events[i, ]
    near_miss_flags[i] <- check_near_miss(
      ev$event_key, ev$event_type, ev$gene_id_unversioned,
      rmats_junctions, pb_junctions, tolerance = 2
    )
  }
  n_near <- sum(near_miss_flags)
  cat(sprintf("  Events that fail exact match but would match at ±2bp: %d / %d (%.1f%%)\n",
              n_near, nrow(missed_events),
              100 * n_near / nrow(missed_events)))
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Frame-Disruption → PTC Verification (Tests H1)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Section 4: Frame-disruption → PTC verification ──────────────────\n\n")

# For events with both forms in PacBio, check PTC status of matched isoforms
both_events <- junction_concordance %>%
  filter(junction_category == "both")

cat(sprintf("  Events with both forms: %d\n", nrow(both_events)))

ptc_verification <- vector("list", nrow(both_events))
for (i in seq_len(nrow(both_events))) {
  ev <- both_events[i, ]

  # Get isoform IDs for both forms
  if (ev$event_type == "SE") {
    form1_isos <- strsplit(ev$inclusion_isos, ";")[[1]]  # inclusion = has the exon
    form2_isos <- strsplit(ev$skipping_isos, ";")[[1]]    # skipping = lacks the exon
  } else {
    form1_isos <- strsplit(ev$inclusion_isos, ";")[[1]]   # long form
    form2_isos <- strsplit(ev$skipping_isos, ";")[[1]]    # short form
  }

  # Look up PTC status
  form1_ptc <- ptc_status[isoform_id %in% form1_isos]
  form2_ptc <- ptc_status[isoform_id %in% form2_isos]

  # Determine if PTC status differs between forms
  form1_any_ptc <- any(form1_ptc$has_ptc == TRUE, na.rm = TRUE)
  form2_any_ptc <- any(form2_ptc$has_ptc == TRUE, na.rm = TRUE)
  form1_all_na  <- all(is.na(form1_ptc$has_ptc)) || nrow(form1_ptc) == 0
  form2_all_na  <- all(is.na(form2_ptc$has_ptc)) || nrow(form2_ptc) == 0

  ptc_differs <- if (form1_all_na || form2_all_na) {
    NA
  } else {
    form1_any_ptc != form2_any_ptc
  }

  ptc_verification[[i]] <- tibble(
    event_key = ev$event_key,
    cell_type = ev$cell_type,
    event_type = ev$event_type,
    gene_id = ev$gene_id,
    gene_symbol = ev$gene_symbol,
    n_form1_isos = length(form1_isos),
    n_form2_isos = length(form2_isos),
    form1_any_ptc = form1_any_ptc,
    form2_any_ptc = form2_any_ptc,
    form1_ptc_known = !form1_all_na,
    form2_ptc_known = !form2_all_na,
    ptc_differs = ptc_differs,
    ptc_in_inclusion = form1_any_ptc,
    ptc_in_skipping = form2_any_ptc
  )
}

ptc_verify_df <- bind_rows(ptc_verification)

cat("\n  PTC verification results (frame-disrupting events with both PB forms):\n")
ptc_verify_df %>%
  filter(!is.na(ptc_differs)) %>%
  count(cell_type, event_type, ptc_differs) %>%
  pivot_wider(names_from = ptc_differs, values_from = n, values_fill = 0,
              names_prefix = "ptc_differs_") %>%
  print(n = Inf)

ptc_diff_rate <- ptc_verify_df %>%
  filter(!is.na(ptc_differs)) %>%
  summarise(
    total = n(),
    ptc_differs = sum(ptc_differs),
    pct = round(100 * ptc_differs / total, 1)
  )
cat(sprintf("\n  Frame-disrupting events where PTC actually differs: %d / %d (%.1f%%)\n",
            ptc_diff_rate$ptc_differs, ptc_diff_rate$total, ptc_diff_rate$pct))

# Also check frame-preserving events as control
sig_fp <- all_events %>%
  filter(cell_type %in% FOCUS_CTS,
         frame_class == "cds_frame_preserving",
         is_significant)

if (nrow(sig_fp) > 0) {
  cat("\n  Running frame-preserving control...\n")

  fp_junctions <- extract_rmats_junctions(sig_fp)
  fp_unique <- sig_fp %>%
    distinct(cell_type, event_type, ID, gene_id_unversioned, .keep_all = TRUE) %>%
    mutate(event_key = paste(cell_type, event_type, ID, sep = ":"))

  fp_concordance <- vector("list", nrow(fp_unique))
  for (i in seq_len(nrow(fp_unique))) {
    ev <- fp_unique[i, ]
    res <- classify_event_junctions(
      ev$event_key, ev$event_type, ev$gene_id_unversioned,
      fp_junctions, pb_junctions
    )
    fp_concordance[[i]] <- tibble(
      event_key = ev$event_key,
      event_type = ev$event_type,
      gene_id = ev$gene_id_unversioned,
      junction_category = res$category,
      inclusion_isos = paste(if (ev$event_type == "SE") res$inclusion_isos else res$long_isos, collapse = ";"),
      skipping_isos = paste(if (ev$event_type == "SE") res$skipping_isos else res$short_isos, collapse = ";")
    )
  }
  fp_conc_df <- bind_rows(fp_concordance)

  fp_both <- fp_conc_df %>% filter(junction_category == "both")
  if (nrow(fp_both) > 0) {
    fp_ptc_results <- vector("list", nrow(fp_both))
    for (i in seq_len(nrow(fp_both))) {
      ev <- fp_both[i, ]
      form1_isos <- strsplit(ev$inclusion_isos, ";")[[1]]
      form2_isos <- strsplit(ev$skipping_isos, ";")[[1]]
      form1_ptc <- ptc_status[isoform_id %in% form1_isos]
      form2_ptc <- ptc_status[isoform_id %in% form2_isos]
      form1_any <- any(form1_ptc$has_ptc == TRUE, na.rm = TRUE)
      form2_any <- any(form2_ptc$has_ptc == TRUE, na.rm = TRUE)
      form1_na  <- all(is.na(form1_ptc$has_ptc)) || nrow(form1_ptc) == 0
      form2_na  <- all(is.na(form2_ptc$has_ptc)) || nrow(form2_ptc) == 0
      fp_ptc_results[[i]] <- tibble(
        ptc_differs = if (form1_na || form2_na) NA else form1_any != form2_any
      )
    }
    fp_ptc_df <- bind_rows(fp_ptc_results)
    fp_diff_rate <- fp_ptc_df %>%
      filter(!is.na(ptc_differs)) %>%
      summarise(total = n(), ptc_differs = sum(ptc_differs),
                pct = round(100 * ptc_differs / total, 1))
    cat(sprintf("  Frame-PRESERVING control: PTC differs in %d / %d (%.1f%%)\n",
                fp_diff_rate$ptc_differs, fp_diff_rate$total, fp_diff_rate$pct))
  }
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5: PSI Concordance (Tests Both)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Section 5: PSI concordance ──────────────────────────────────────\n\n")

# For events with both forms, compute PacBio PSI and compare to rMATS PSI

# Get sample metadata from DGEList (keep.rownames preserves sample IDs)
sample_meta <- as.data.table(dge$samples, keep.rownames = "sample_id")
# Ensure ct is character for matching (it may be a factor)
sample_meta[, ct := as.character(ct)]

# Get counts matrix (will subset columns/rows as needed)
counts_mat <- dge$counts
txids <- if (source_type == "isocall") dge$genes$transcript_id else dge$genes$txid

psi_results <- vector("list", nrow(both_events))
for (i in seq_len(nrow(both_events))) {
  ev <- both_events[i, ]

  # Get isoform IDs (inclusion_isos = SE inclusion or A3SS/A5SS long form)
  incl_isos <- strsplit(ev$inclusion_isos, ";")[[1]]
  skip_isos <- strsplit(ev$skipping_isos, ";")[[1]]

  # Find matching rows in counts matrix
  incl_rows <- which(txids %in% incl_isos)
  skip_rows <- which(txids %in% skip_isos)
  if (length(incl_rows) == 0 || length(skip_rows) == 0) next

  # Find matching samples for this cell type
  ct_samples <- sample_meta[ct == ev$cell_type]
  if (nrow(ct_samples) == 0) next

  sample_cols <- which(colnames(counts_mat) %in% ct_samples$sample_id)
  if (length(sample_cols) == 0) next

  # Compute per-sample counts
  incl_counts <- if (length(incl_rows) == 1) {
    counts_mat[incl_rows, sample_cols]
  } else {
    colSums(counts_mat[incl_rows, sample_cols, drop = FALSE])
  }

  total_counts <- incl_counts + if (length(skip_rows) == 1) {
    counts_mat[skip_rows, sample_cols]
  } else {
    colSums(counts_mat[skip_rows, sample_cols, drop = FALSE])
  }

  # PSI = inclusion / total (set NA where total == 0)
  pb_psi <- ifelse(total_counts > 0, incl_counts / total_counts, NA_real_)

  # Get treatment for each sample (match by sample_id)
  sample_treatment <- ct_samples[match(colnames(counts_mat)[sample_cols],
                                        ct_samples$sample_id), treatment]

  dmso_psi <- pb_psi[sample_treatment == "DMSO"]
  smg1i_psi <- pb_psi[sample_treatment == "Smg1i"]

  pb_dpsi <- mean(smg1i_psi, na.rm = TRUE) - mean(dmso_psi, na.rm = TRUE)

  psi_results[[i]] <- tibble(
    event_key = ev$event_key,
    cell_type = ev$cell_type,
    event_type = ev$event_type,
    gene_id = ev$gene_id,
    gene_symbol = ev$gene_symbol,
    rmats_dPSI = ev$dPSI,
    pb_dPSI = pb_dpsi,
    pb_mean_dmso_psi = mean(dmso_psi, na.rm = TRUE),
    pb_mean_smg1i_psi = mean(smg1i_psi, na.rm = TRUE),
    pb_total_counts_mean = mean(total_counts, na.rm = TRUE),
    n_incl_isos = length(incl_isos),
    n_skip_isos = length(skip_isos)
  )

  if (i %% 200 == 0) cat(sprintf("    %d / %d\n", i, nrow(both_events)))
}

psi_concordance <- bind_rows(psi_results[!sapply(psi_results, is.null)])

cat(sprintf("  PSI concordance computed for %d events\n", nrow(psi_concordance)))

if (nrow(psi_concordance) > 0) {
  # Filter to events with reasonable coverage
  psi_good <- psi_concordance %>%
    filter(!is.na(pb_dPSI), pb_total_counts_mean >= 5)

  if (nrow(psi_good) > 0) {
    # Spearman correlation of dPSI
    cor_test <- cor.test(psi_good$rmats_dPSI, psi_good$pb_dPSI,
                         method = "spearman")
    cat(sprintf("  dPSI Spearman correlation: rho = %.3f, p = %.2e (n = %d)\n",
                cor_test$estimate, cor_test$p.value, nrow(psi_good)))

    # Sign concordance
    sign_conc <- psi_good %>%
      mutate(same_sign = sign(rmats_dPSI) == sign(pb_dPSI)) %>%
      summarise(total = n(), concordant = sum(same_sign),
                pct = round(100 * concordant / total, 1))
    cat(sprintf("  Sign concordance: %d / %d (%.1f%%)\n",
                sign_conc$concordant, sign_conc$total, sign_conc$pct))

    # By cell type
    cat("\n  By cell type:\n")
    psi_good %>%
      group_by(cell_type) %>%
      summarise(
        n = n(),
        rho = cor(rmats_dPSI, pb_dPSI, method = "spearman"),
        sign_concordance = round(100 * mean(sign(rmats_dPSI) == sign(pb_dPSI)), 1),
        .groups = "drop"
      ) %>%
      print(n = Inf)
  }
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6: Discrepancy Waterfall
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Section 6: Discrepancy waterfall ────────────────────────────────\n\n")

# Sequential bucket assignment for ALL significant frame-disrupting CDS events

# Start from junction_concordance which has all sig FD events
waterfall <- junction_concordance %>%
  select(event_key, cell_type, event_type, gene_id, gene_symbol,
         dPSI, FDR, junction_category, inclusion_isos, skipping_isos)

# Bucket 1: No PacBio gene coverage
waterfall <- waterfall %>%
  mutate(bucket = case_when(
    junction_category == "no_pb_gene" ~ "1_no_pb_gene",
    TRUE ~ NA_character_
  ))

# Bucket 2: No junction match (gene present but neither/one form only)
waterfall <- waterfall %>%
  mutate(bucket = case_when(
    !is.na(bucket) ~ bucket,
    junction_category %in% c("neither", "inclusion_only", "skipping_only") ~ "2_no_junction_match",
    TRUE ~ bucket
  ))

# Bucket 3-5: For "both" events, check PTC and NMD sensitivity
# Join PTC verification
waterfall <- waterfall %>%
  left_join(
    ptc_verify_df %>% select(event_key, ptc_differs, ptc_in_inclusion, ptc_in_skipping),
    by = "event_key"
  )

# Bucket 3: Both forms exist but no PTC difference
waterfall <- waterfall %>%
  mutate(bucket = case_when(
    !is.na(bucket) ~ bucket,
    junction_category == "both" & (is.na(ptc_differs) | !ptc_differs) ~ "3_no_ptc_difference",
    TRUE ~ bucket
  ))

# For bucket 4-5, need to check if PTC+ isoform is NMD-sensitive in DE results
# Look up DE status for inclusion/skipping isoforms
check_nmd_sensitivity <- function(event_key, cell_type, inclusion_isos_str,
                                   skipping_isos_str, ptc_in_incl, ptc_in_skip,
                                   de_results_list) {
  de <- de_results_list[[cell_type]]
  if (is.null(de)) return(NA)

  # Get PTC+ isoforms
  incl_isos <- strsplit(inclusion_isos_str, ";")[[1]]
  skip_isos <- strsplit(skipping_isos_str, ";")[[1]]

  ptc_isos <- character(0)
  if (ptc_in_incl) ptc_isos <- c(ptc_isos, incl_isos)
  if (ptc_in_skip) ptc_isos <- c(ptc_isos, skip_isos)

  if (length(ptc_isos) == 0) return(NA)

  # Check DE status
  ptc_de <- de[txid %in% ptc_isos]
  if (nrow(ptc_de) == 0) return(FALSE)

  # NMD-sensitive: logFC > 0 AND adj.P.Val < 0.05
  any(ptc_de$logFC > 0 & ptc_de$adj.P.Val < 0.05, na.rm = TRUE)
}

# Apply to events with PTC difference
ptc_diff_events <- waterfall %>%
  filter(is.na(bucket), junction_category == "both", !is.na(ptc_differs), ptc_differs)

if (nrow(ptc_diff_events) > 0) {
  nmd_flags <- logical(nrow(ptc_diff_events))
  for (i in seq_len(nrow(ptc_diff_events))) {
    ev <- ptc_diff_events[i, ]
    nmd_flags[i] <- check_nmd_sensitivity(
      ev$event_key, ev$cell_type, ev$inclusion_isos, ev$skipping_isos,
      ev$ptc_in_inclusion, ev$ptc_in_skipping, de_results
    )
  }

  ptc_diff_events$is_nmd_sensitive <- nmd_flags

  # Bucket 4: PTC exists but isoform not NMD-sensitive
  # Bucket 5: Concordant — PTC+ isoform is NMD-sensitive
  ptc_diff_events <- ptc_diff_events %>%
    mutate(bucket = case_when(
      is.na(is_nmd_sensitive) | !is_nmd_sensitive ~ "4_ptc_not_nmd_sensitive",
      is_nmd_sensitive ~ "5_concordant",
      TRUE ~ NA_character_
    ))

  # Update waterfall
  waterfall <- waterfall %>%
    rows_update(
      ptc_diff_events %>% select(event_key, bucket),
      by = "event_key"
    )
}

# Any remaining unclassified
waterfall <- waterfall %>%
  mutate(bucket = if_else(is.na(bucket), "unclassified", bucket))

cat("  Waterfall results:\n")
waterfall %>%
  count(cell_type, event_type, bucket) %>%
  pivot_wider(names_from = bucket, values_from = n, values_fill = 0) %>%
  print(n = Inf)

# Overall waterfall percentages
cat("\n  Overall waterfall (all cell types combined):\n")
waterfall %>%
  count(bucket) %>%
  mutate(pct = round(100 * n / sum(n), 1),
         cumulative_pct = round(100 * cumsum(n) / sum(n), 1)) %>%
  print(n = Inf)

# By event type
cat("\n  By event type:\n")
waterfall %>%
  count(event_type, bucket) %>%
  group_by(event_type) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  print(n = Inf)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7: Summary and Save
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Section 7: Summary and save ─────────────────────────────────────\n\n")

# --- Save output files ---
write_tsv(junction_concordance, file.path(out_dir, "rmats_junction_concordance.tsv"))
cat("  Wrote rmats_junction_concordance.tsv\n")

write_tsv(ptc_verify_df, file.path(out_dir, "rmats_ptc_verification.tsv"))
cat("  Wrote rmats_ptc_verification.tsv\n")

if (nrow(psi_concordance) > 0) {
  write_tsv(psi_concordance, file.path(out_dir, "rmats_psi_concordance.tsv"))
  cat("  Wrote rmats_psi_concordance.tsv\n")
}

write_tsv(waterfall, file.path(out_dir, "rmats_discrepancy_waterfall.tsv"))
cat("  Wrote rmats_discrepancy_waterfall.tsv\n")

# --- Summary statistics ---
summary_stats <- tibble(
  metric = character(),
  value = character(),
  cell_type = character(),
  event_type = character()
)

# Junction concordance by cell type × event type
for (ct in FOCUS_CTS) {
  for (et in EVENT_TYPES) {
    subset <- junction_concordance %>%
      filter(cell_type == ct, event_type == et)
    if (nrow(subset) == 0) next

    n_total <- nrow(subset)
    n_both <- sum(subset$junction_category == "both")
    n_no_gene <- sum(subset$junction_category == "no_pb_gene")
    n_neither <- sum(subset$junction_category == "neither")
    n_incl_only <- sum(subset$junction_category == "inclusion_only")
    n_skip_only <- sum(subset$junction_category == "skipping_only")

    summary_stats <- bind_rows(summary_stats, tibble(
      metric = c("total_sig_fd_events", "both_forms_in_pb", "no_pb_gene",
                  "neither_junction", "inclusion_only", "skipping_only",
                  "pct_both_forms"),
      value = as.character(c(n_total, n_both, n_no_gene, n_neither,
                              n_incl_only, n_skip_only,
                              round(100 * n_both / n_total, 1))),
      cell_type = ct,
      event_type = et
    ))
  }
}

# PTC verification summary
ptc_known <- ptc_verify_df %>% filter(!is.na(ptc_differs))
if (nrow(ptc_known) > 0) {
  summary_stats <- bind_rows(summary_stats, tibble(
    metric = c("ptc_verification_total", "ptc_differs_count", "ptc_differs_pct"),
    value = as.character(c(nrow(ptc_known),
                            sum(ptc_known$ptc_differs),
                            round(100 * sum(ptc_known$ptc_differs) / nrow(ptc_known), 1))),
    cell_type = "all",
    event_type = "all"
  ))
}

# PSI concordance summary
if (exists("psi_good") && nrow(psi_good) > 0) {
  summary_stats <- bind_rows(summary_stats, tibble(
    metric = c("psi_concordance_n", "psi_spearman_rho", "psi_sign_concordance_pct"),
    value = as.character(c(nrow(psi_good),
                            round(cor_test$estimate, 3),
                            sign_conc$pct)),
    cell_type = "all",
    event_type = "all"
  ))
}

# Waterfall summary
wf_summary <- waterfall %>%
  count(bucket) %>%
  mutate(pct = round(100 * n / sum(n), 1))

summary_stats <- bind_rows(summary_stats,
  wf_summary %>%
    transmute(metric = paste0("waterfall_", bucket),
              value = sprintf("%d (%.1f%%)", n, pct),
              cell_type = "all", event_type = "all")
)

write_tsv(summary_stats, file.path(out_dir, "rmats_concordance_summary.tsv"))
cat("  Wrote rmats_concordance_summary.tsv\n")

# --- Console summary ---
cat("\n")
cat("==================================================================\n")
cat("   CONCORDANCE ANALYSIS SUMMARY\n")
cat("==================================================================\n\n")

cat("  H2 (Reconstruction error) — Junction concordance:\n")
junction_concordance %>%
  count(junction_category) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  mutate(label = sprintf("    %-20s %5d (%5.1f%%)", junction_category, n, pct)) %>%
  pull(label) %>%
  cat(sep = "\n")

cat("\n\n  H1 (Reconcilable) — Frame disruption → PTC:\n")
cat(sprintf("    Events with both %s forms and known PTC: %d\n", source_label, nrow(ptc_known)))
if (nrow(ptc_known) > 0) {
  cat(sprintf("    PTC actually differs between forms:      %d (%.1f%%)\n",
              sum(ptc_known$ptc_differs),
              100 * sum(ptc_known$ptc_differs) / nrow(ptc_known)))
  cat(sprintf("    PTC same (frame disruption ≠ PTC):       %d (%.1f%%)\n",
              sum(!ptc_known$ptc_differs),
              100 * sum(!ptc_known$ptc_differs) / nrow(ptc_known)))
}

cat("\n  Discrepancy waterfall (where rMATS signal is lost):\n")
wf_summary %>%
  mutate(label = sprintf("    %-30s %5d (%5.1f%%)", bucket, n, pct)) %>%
  pull(label) %>%
  cat(sep = "\n")

cat("\n\n  Done.\n")
