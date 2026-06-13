#!/usr/bin/env Rscript
# Data export for new Fig 3 (Isopair + PTC combined multipanel).
#
# 2026-06-13 PM rewrite — population structure based on Pete's clarified scope:
#
# Scope policy:
#   "TD2 CDS annotations are unreliable, so for analyses that depend on
#    identifying the stop codon, we limit to isoform pairs where reference
#    AUG tracing can be performed."
#
# Three populations used in this figure:
#
#   pop_BC (n=3,009 NMD / 3,009 Control):
#       Stage 2 gene-matched pairs. No CDS dependency. Used by Panels B
#       (sequence similarity) and C (splice event prevalence).
#
#   pop_traceable (n=2,289 NMD / 1,763 Control):
#       Pairs where ref-AUG tracing produced a valid ORF position
#       (categories: effectively_ptc + no_downstream_ejc + truncated_no_ejc).
#       Used by Panel D (stop-to-last-EJC distance density). Stop position
#       is the ref-AUG-traced comp_stop_tx_pos.
#
#   pop_ptc_plus (n=1,912 NMD / 288 Control):
#       Ref-AUG PTC+ pairs (effectively_ptc only). Used by Panels E and F.
#       Mechanism attribution combined from two sources:
#         - ref_atg_analysis$c2$attr_mechanism / attr_event for the
#           original_ptc=FALSE subset (900 NMD pairs — the ref-AUG-recovered)
#         - existing all_attr (current goal2-ptc-mechanisms output) for
#           original_ptc=TRUE pairs that are also effectively_ptc
#           (1,012 NMD pairs)

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Paths ----
script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) > 0) return(normalizePath(sub("^--file=", "", file_arg[1])))
  if (!is.null(sys.frame(1)$ofile)) return(normalizePath(sys.frame(1)$ofile))
  return(normalizePath(file.path(getwd(), "data_export.R")))
})()
HERE <- dirname(script_path)
ISOPAIR_DIR <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper"
DM_DIR <- file.path(ISOPAIR_DIR, "data_mashr")
CACHE_DIR <- file.path(DM_DIR, "analysis_cache")
OUT_DIR <- file.path(HERE, "data")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[data_export.R] Writing to: %s\n", OUT_DIR))

# ---- Load common data ----
cat("[setup] Loading profiles + cds + ptc + ref_atg + structures\n")
profiles_c2 <- as.data.table(readRDS(file.path(DM_DIR, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM_DIR, "profiles_c4_allsamples.rds")))
cds <- as.data.table(readRDS(file.path(DM_DIR, "cds.rds")))
ptc <- as.data.table(readRDS(file.path(DM_DIR, "ptc.rds")))
structures <- as.data.table(readRDS(file.path(DM_DIR, "structures.rds")))
ref_atg <- readRDS(file.path(CACHE_DIR, "ref_atg_analysis.rds"))

# ============================================================================
# Population 1 — pop_BC: Stage 2 gene-matched (B + C)
# ============================================================================
make_key <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")
shared_keys <- intersect(unique(make_key(profiles_c2)),
                         unique(make_key(profiles_c4)))
pop_BC_c2 <- profiles_c2[make_key(profiles_c2) %in% shared_keys]
pop_BC_c4 <- profiles_c4[make_key(profiles_c4) %in% shared_keys]
cat(sprintf("[pop_BC] Stage 2 gene-matched: NMD=%d, Control=%d\n",
            nrow(pop_BC_c2), nrow(pop_BC_c4)))

# ============================================================================
# Population 2 — pop_traceable: ref-AUG tracing succeeded
# ============================================================================
TRACEABLE_CATS <- c("effectively_ptc", "no_downstream_ejc", "truncated_no_ejc")
ra_c2 <- as.data.table(ref_atg$c2)
ra_c4 <- as.data.table(ref_atg$c4)
pop_trace_c2_ids <- ra_c2$comparator_isoform_id[ra_c2$category %in% TRACEABLE_CATS]
pop_trace_c4_ids <- ra_c4$comparator_isoform_id[ra_c4$category %in% TRACEABLE_CATS]
pop_trace_c2 <- pop_BC_c2[comparator_isoform_id %in% pop_trace_c2_ids]
pop_trace_c4 <- pop_BC_c4[comparator_isoform_id %in% pop_trace_c4_ids]
cat(sprintf("[pop_traceable] ref-AUG traceable: NMD=%d, Control=%d\n",
            nrow(pop_trace_c2), nrow(pop_trace_c4)))

# ============================================================================
# Population 3 — pop_ptc_plus: ref-AUG-defined PTC+ (effectively_ptc)
# ============================================================================
pop_ptc_c2_ids <- ra_c2$comparator_isoform_id[ra_c2$category == "effectively_ptc"]
pop_ptc_c4_ids <- ra_c4$comparator_isoform_id[ra_c4$category == "effectively_ptc"]
pop_ptc_c2 <- pop_BC_c2[comparator_isoform_id %in% pop_ptc_c2_ids]
pop_ptc_c4 <- pop_BC_c4[comparator_isoform_id %in% pop_ptc_c4_ids]
cat(sprintf("[pop_ptc_plus] ref-AUG PTC+: NMD=%d, Control=%d\n",
            nrow(pop_ptc_c2), nrow(pop_ptc_c4)))
cat(sprintf("  NMD PTC rate: %.1f%% | Control PTC rate: %.1f%%\n",
            100 * nrow(pop_ptc_c2) / nrow(pop_trace_c2),
            100 * nrow(pop_ptc_c4) / nrow(pop_trace_c4)))

# ============================================================================
# Panel B — sequence similarity NMD vs Control
# Population: pop_BC (3,009 each)
# ============================================================================
cat("\n[Panel B] Sequence similarity (pop_BC = 3,009 each)\n")
div_c2 <- readRDS(file.path(CACHE_DIR, "div_c2_allsamples.rds"))
div_c4 <- readRDS(file.path(CACHE_DIR, "div_c4_allsamples.rds"))
# div is row-aligned to pop_BC (Stage 2 gene-matched)
stopifnot(nrow(div_c2) == nrow(pop_BC_c2))
stopifnot(nrow(div_c4) == nrow(pop_BC_c4))
panelB <- rbind(
  data.table(comparison = "NMD",
             pct_shared = div_c2$pct_shared),
  data.table(comparison = "Control",
             pct_shared = div_c4$pct_shared)
)
panelB <- panelB[!is.na(pct_shared)]
fwrite(panelB, file.path(OUT_DIR, "panelB_sequence_similarity.tsv"), sep = "\t")
cat(sprintf("  -> NMD=%d, Control=%d non-NA pct_shared\n",
            sum(panelB$comparison == "NMD"), sum(panelB$comparison == "Control")))

# ============================================================================
# Panel C — splice event prevalence NMD vs Control
# Population: pop_BC (3,009 each)
# ============================================================================
cat("[Panel C] Event prevalence (pop_BC = 3,009 each)\n")
evt_types <- c("Alt_TSS", "Alt_TES", "A5SS", "A3SS", "SE",
               "Missing_Internal", "IR", "IR_diff",
               "Partial_IR_5", "Partial_IR_3")

compute_event_freq <- function(profiles) {
  n_pir5 <- sum(sapply(profiles$detailed_events, function(de)
    any(de$event_type == "Partial_IR_5")))
  n_pir3 <- sum(sapply(profiles$detailed_events, function(de)
    any(de$event_type == "Partial_IR_3")))
  counts <- c(
    sum(profiles$tss_changed), sum(profiles$tes_changed),
    sum(profiles$n_a5ss > 0), sum(profiles$n_a3ss > 0),
    sum(profiles$n_se > 0), sum(profiles$n_missing_internal > 0),
    sum(profiles$n_ir > 0), sum(profiles$n_ir_diff > 0),
    n_pir5, n_pir3
  )
  data.frame(event_type = evt_types, n_with_event = counts,
             pct = round(100 * counts / nrow(profiles), 1))
}

ef_c2 <- compute_event_freq(pop_BC_c2)
ef_c4 <- compute_event_freq(pop_BC_c4)
panelC <- merge(ef_c2, ef_c4, by = "event_type", suffixes = c("_NMD", "_Ctrl"))
N_BC <- nrow(pop_BC_c2)
panelC$fisher_OR <- NA_real_
panelC$fisher_p <- NA_real_
for (i in seq_len(nrow(panelC))) {
  a <- panelC$n_with_event_NMD[i]; b <- N_BC - a
  c <- panelC$n_with_event_Ctrl[i]; d <- N_BC - c
  ft <- tryCatch(fisher.test(matrix(c(a, b, c, d), nrow = 2)),
                 error = function(e) NULL)
  if (!is.null(ft)) {
    panelC$fisher_OR[i] <- round(ft$estimate, 2)
    panelC$fisher_p[i] <- ft$p.value
  }
}
panelC$direction <- ifelse(panelC$fisher_p > 0.05, "Non-significant",
  ifelse(panelC$fisher_OR > 1, "Enriched in NMD", "Enriched in Control"))
panelC <- panelC[order(panelC$fisher_p), ]
setnames(panelC,
         old = c("n_with_event_NMD", "pct_NMD", "n_with_event_Ctrl", "pct_Ctrl"),
         new = c("n_pairs_with_event_NMD", "pct_of_pairs_NMD",
                 "n_pairs_with_event_Ctrl", "pct_of_pairs_Ctrl"))
fwrite(panelC, file.path(OUT_DIR, "panelC_event_prevalence.tsv"), sep = "\t")
cat(sprintf("  -> %d event types (denom=%d each)\n", nrow(panelC), N_BC))

# ============================================================================
# Panel D — Stop-to-last-EJC distance NMD vs Control
# Population: pop_traceable (NMD=2,289 / Control=1,763)
# Stop: ref-AUG-traced comp_stop_tx_pos
# ============================================================================
cat("[Panel D] Stop-to-last-EJC distance (pop_traceable, ref-AUG stops)\n")

compute_last_ejc <- function(starts, ends, strand, n_exons) {
  if (n_exons < 2) return(NA_integer_)
  lens <- ends - starts + 1
  if (strand == "-") lens <- rev(lens)
  as.integer(sum(lens[1:(n_exons - 1)]))
}
structures[, last_ejc_tx_pos := mapply(compute_last_ejc,
                                       exon_starts, exon_ends, strand, n_exons)]
struct_lookup <- structures[, .(isoform_id, last_ejc_tx_pos)]

build_distance_df <- function(pop, ra, label) {
  d <- merge(pop[, .(comparator_isoform_id, gene_id, reference_isoform_id)],
             ra[, .(comparator_isoform_id, category, comp_stop_tx_pos)],
             by = "comparator_isoform_id")
  d <- merge(d, struct_lookup, by.x = "comparator_isoform_id",
             by.y = "isoform_id", all.x = TRUE)
  # ref-AUG distance: positive = stop upstream of last EJC = PTC direction
  d[, distance := last_ejc_tx_pos - comp_stop_tx_pos]
  d[, comparison := label]
  d[, .(comparator_isoform_id, comparison, distance, category)]
}
panelD <- rbind(
  build_distance_df(pop_trace_c2, ra_c2, "NMD"),
  build_distance_df(pop_trace_c4, ra_c4, "Control")
)
fwrite(panelD[!is.na(distance)],
       file.path(OUT_DIR, "panelD_stop_codon_distance.tsv"), sep = "\t")
cat(sprintf("  -> NMD=%d, Control=%d non-NA distances\n",
            sum(panelD$comparison == "NMD" & !is.na(panelD$distance)),
            sum(panelD$comparison == "Control" & !is.na(panelD$distance))))

# ============================================================================
# Panel E + F — PTC-causing event attribution + mechanism breakdown
# Population: pop_ptc_plus (NMD=1,912 effectively_ptc)
# Attribution sources:
#   - For original_ptc=FALSE (900 pairs): ref_atg_analysis$c2$attr_event / attr_mechanism
#   - For original_ptc=TRUE (1,012 pairs): goal2-ptc-mechanisms output (current all_attr)
# Control baseline: all events in pop_trace_c4 (1,763 pairs)
# ============================================================================
cat("[Panel E+F] PTC-causing attribution at pop_ptc_plus scope\n")
source(file.path(HERE, "panel_e_compute.R"))

cat("\n[data_export.R] done.\n")
cat(sprintf("Population summary:\n"))
cat(sprintf("  pop_BC:        NMD=%d Control=%d\n", nrow(pop_BC_c2), nrow(pop_BC_c4)))
cat(sprintf("  pop_traceable: NMD=%d Control=%d\n", nrow(pop_trace_c2), nrow(pop_trace_c4)))
cat(sprintf("  pop_ptc_plus:  NMD=%d Control=%d\n", nrow(pop_ptc_c2), nrow(pop_ptc_c4)))
