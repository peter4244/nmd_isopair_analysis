# Panel E + F recomputation — all-3-ENST GENCODE-CDS framework.
#
# Sourced by data_export.R after pop_ptc_c2 / pop_trace_c4 (all-3-ENST scope,
# PTC determined from each comparator's own GENCODE-annotated stop) are defined.
#
# Approach (all-3-ENST scope is fully GENCODE-anchored — no TD2 dependence
# anywhere; no ref-AUG projection needed):
#   For each PTC+ NMD comparator, attribute the splice event that introduced
#   the PTC using attribute_ptc_events() / attribute_3utr_splice() from
#   analysis_functions.R, keyed on the COMPARATOR'S OWN GENCODE-annotated stop
#   (genomic position from cds.rds; converted to genomic stop position from
#   cds_start/cds_stop based on strand).
#
# Control side: no per-pair mechanism attribution — only the all-events
# baseline from pop_trace_c4 (all-3-ENST) detailed_events.

source(file.path(ISOPAIR_DIR, "analysis_functions.R"))

# Get the comparator-own GENCODE stop genomic position per isoform.
# For + strand: cds_stop is the 3' end of CDS = stop codon position.
# For - strand: cds_start is the 3' end of CDS = stop codon position.
cds_dt <- as.data.table(cds)[coding_status == "coding"]
own_stop_genomic_lookup <- setNames(
  ifelse(cds_dt$strand == "+", cds_dt$cds_stop, cds_dt$cds_start),
  cds_dt$isoform_id)
own_atg_genomic_lookup <- setNames(
  ifelse(cds_dt$strand == "+", cds_dt$cds_start, cds_dt$cds_stop),
  cds_dt$isoform_id)
strand_lookup <- setNames(cds_dt$strand, cds_dt$isoform_id)

# Pop_ptc_c2 already has the PTC+ subset (own GENCODE stop > 50nt past last EJC).
pop_ptc_c2_dt <- as.data.table(pop_ptc_c2)
ep_ids <- pop_ptc_c2_dt$comparator_isoform_id

# is_frameshift flag from fc_c2_cache (frame comparison between comp & ref CDS)
fc_c2_cache <- readRDS(file.path(CACHE_DIR, "fc_c2_allsamples.rds"))
fw_c2 <- readRDS(file.path(CACHE_DIR, "fw_c2_allsamples.rds"))
fc_ps <- as.data.table(fc_c2_cache$pair_summary)
fc_ep <- fc_ps[comparator_isoform_id %in% ep_ids]
frameshift_cats <- c("same_start_frameshift",
                     "diff_start_diff_frame_with_frameshift")
fc_ep[, is_frameshift_cat := frame_category %in% frameshift_cats]
cat(sprintf("  [E/F-NMD] all-3-ENST PTC+: %d pairs\n", length(ep_ids)))

# Build attribute_ptc_events input — PTC stop = comparator's own GENCODE stop
ptc_genomic_pos <- own_stop_genomic_lookup[fc_ep$comparator_isoform_id]
atg_genomic_pos <- own_atg_genomic_lookup[fc_ep$comparator_isoform_id]
strand_vec <- strand_lookup[fc_ep$comparator_isoform_id]
names(ptc_genomic_pos) <- fc_ep$comparator_isoform_id
names(atg_genomic_pos) <- fc_ep$comparator_isoform_id
names(strand_vec) <- fc_ep$comparator_isoform_id
is_fs_lookup <- setNames(fc_ep$is_frameshift_cat, fc_ep$comparator_isoform_id)

attr_b <- attribute_ptc_events(
  pairs = fc_ep[, .(comparator_isoform_id, reference_isoform_id)],
  fw_events = fw_c2$events,
  profiles = as.data.frame(pop_ptc_c2),
  ptc_genomic_pos = ptc_genomic_pos,
  atg_genomic_pos = atg_genomic_pos,
  strand_vec = strand_vec,
  is_frameshift_vec = is_fs_lookup
)
attr_b <- as.data.table(attr_b)
cat(sprintf("  [E/F-NMD] attribute_ptc_events: %d direct + %d unresolved\n",
            sum(attr_b$attribution == "direct"),
            sum(attr_b$attribution == "unresolved")))

# Same-stop 3'UTR splice handling (pairs where comp's own stop == ref's own stop)
fc_ep[, strand_s := strand_lookup[comparator_isoform_id]]
fc_ep[, comp_stop_s := own_stop_genomic_lookup[comparator_isoform_id]]
fc_ep[, ref_stop_s := own_stop_genomic_lookup[reference_isoform_id]]
fc_ep[, same_stop_s := comp_stop_s == ref_stop_s]
same_stop_b <- fc_ep[same_stop_s == TRUE]
ss_events_b <- data.table()
if (nrow(same_stop_b) > 0) {
  ss_stop_vec <- setNames(same_stop_b$comp_stop_s, same_stop_b$comparator_isoform_id)
  ss_strand_vec <- setNames(same_stop_b$strand_s, same_stop_b$comparator_isoform_id)
  ss_events_b <- as.data.table(attribute_3utr_splice(
    pairs = same_stop_b[, .(comparator_isoform_id, reference_isoform_id)],
    profiles = as.data.frame(pop_ptc_c2),
    stop_genomic_pos = ss_stop_vec,
    strand_vec = ss_strand_vec
  ))[attribution == "direct"]
  cat(sprintf("  [E/F-NMD] same-stop 3'UTR splice attributions: %d\n",
              nrow(ss_events_b)))
}

# Combine: for diff-stop, use attribute_ptc_events direct; for same-stop,
# use attribute_3utr_splice (overwrites diff-stop attribution for these pairs).
same_stop_ids <- if (nrow(ss_events_b) > 0) ss_events_b$comparator_isoform_id else character(0)
diffstop_b <- attr_b[attribution == "direct" &
                     !comparator_isoform_id %in% same_stop_ids]
all_attr_new <- rbind(diffstop_b, ss_events_b, fill = TRUE)
all_attr_new <- all_attr_new[attribution == "direct"]
all_attr_new$event_short <- shorten_event_labels(all_attr_new$ptc_causing_event,
                                                  collapse_rare = 3)
cat(sprintf("  [E/F-NMD] all_attr_new (combined direct): %d pairs\n",
            nrow(all_attr_new)))

# ---- Panel E: PTC-causing event proportions vs Control baseline ----
n_ptc_attr <- nrow(all_attr_new)
ptc_evt_freq <- sort(table(all_attr_new$ptc_causing_event), decreasing = TRUE)
ptc_evt_df <- data.table(event_type = names(ptc_evt_freq),
                         n_ptc_events = as.integer(ptc_evt_freq),
                         pct_of_ptc = round(100 * as.integer(ptc_evt_freq) / n_ptc_attr, 1))

# Control baseline: all detailed_events in pop_trace_c4 (all-3-ENST scope).
# pop_trace_c4 here is the data.table from data_export.R augmentation step
# (own_stop annotated). Pull the matching pop_BC_c4_tri pair rows to get
# detailed_events.
ctrl_pop <- pop_BC_c4_tri[comparator_isoform_id %in% pop_trace_c4$comparator_isoform_id]
ctrl_all_events <- do.call(rbind, lapply(ctrl_pop$detailed_events, as.data.frame))
ctrl_event_freq <- table(ctrl_all_events$event_type)
ctrl_total <- sum(ctrl_event_freq)
cat(sprintf("  [E] Control baseline events: %d (across %d pop_trace_c4 pairs)\n",
            ctrl_total, nrow(ctrl_pop)))

ptc_evt_df[, n_ctrl_events := as.integer(ctrl_event_freq[event_type])]
ptc_evt_df[is.na(n_ctrl_events), n_ctrl_events := 0L]
ptc_evt_df[, pct_ctrl := round(100 * n_ctrl_events / ctrl_total, 1)]
ptc_evt_df[, enrichment := round(pct_of_ptc / pmax(pct_ctrl, 0.1), 1)]
ptc_evt_df[, fisher_p := NA_real_]
ptc_evt_df[, direction := NA_character_]
for (i in seq_len(nrow(ptc_evt_df))) {
  a <- ptc_evt_df$n_ptc_events[i]; b <- n_ptc_attr - a
  cv <- ptc_evt_df$n_ctrl_events[i]; dv <- ctrl_total - cv
  ft <- tryCatch(fisher.test(matrix(c(a, b, cv, dv), nrow = 2)),
                 error = function(e) NULL)
  if (!is.null(ft)) {
    ptc_evt_df$fisher_p[i] <- signif(ft$p.value, 3)
    ptc_evt_df$direction[i] <- if (ft$p.value > 0.05) "Non-significant"
      else if (ptc_evt_df$enrichment[i] > 1) "Enriched in PTC-causing"
      else "Depleted in PTC-causing"
  }
}
panelE <- ptc_evt_df[, .(event_type, n_ptc_events, pct_of_ptc,
                         n_ctrl_events, pct_ctrl, enrichment,
                         fisher_p, direction)]
fwrite(panelE, file.path(OUT_DIR, "panelE_ptc_event_attribution.tsv"), sep = "\t")
cat(sprintf("  [E] -> %d event types (n_ptc_attr=%d, ctrl_total=%d)\n",
            nrow(panelE), n_ptc_attr, ctrl_total))

# ---- Panel F: mechanism breakdown × event_type (long form for stacked bar) ----
all_attr_new[, mechanism := factor(mechanism,
                                   levels = c("Frameshift", "In-frame stop", "3'UTR splice"))]
sankey_agg <- all_attr_new[!is.na(mechanism),
                            .(n = .N), by = .(event_short, mechanism)]
evt_totals <- sankey_agg[, .(total = sum(n)), by = event_short][order(-total)]
sankey_agg[, event_short := factor(event_short, levels = evt_totals$event_short)]
panelF <- sankey_agg[order(event_short, mechanism),
                     .(event_type = as.character(event_short),
                       mechanism = as.character(mechanism), n)]
fwrite(panelF, file.path(OUT_DIR, "panelF_mechanism_breakdown.tsv"), sep = "\t")
cat(sprintf("  [F] -> %d (event_type x mechanism) rows; %d events; total attributed=%d\n",
            nrow(panelF), nrow(evt_totals), sum(panelF$n)))
