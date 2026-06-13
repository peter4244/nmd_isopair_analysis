# Panel E + F recomputation — under new ref-AUG-traceable scope.
#
# Sourced by data_export.R after pop_ptc_c2 / pop_trace_c4 are defined.
# Expects in scope: pop_ptc_c2 (1,912 effectively_ptc NMD pairs),
# pop_trace_c4 (1,763 ref-AUG-traceable Control pairs), ra_c2, ra_c4, ptc,
# CACHE_DIR, ISOPAIR_DIR, OUT_DIR.
#
# Attribution sources for NMD effectively_ptc pairs:
#   (a) original_ptc=FALSE (900 of 1,912): ref_atg_analysis$c2$attr_event /
#       attr_mechanism (set by 05r_ref_atg_analysis.R at ref-AUG runtime)
#   (b) original_ptc=TRUE (1,012 of 1,912): existing all_attr from
#       goal2-ptc-mechanisms chunk (TD2-stop-based attribution)
#
# Control side: no per-pair mechanism attribution needed — only the all-events
# baseline computed from pop_trace_c4's detailed_events.

source(file.path(ISOPAIR_DIR, "analysis_functions.R"))
cds_lookups <- build_cds_lookups(cds)
coding_ids <- cds_lookups$coding_ids
strand_lookup <- cds_lookups$strand_lookup
cds_start_lookup <- cds_lookups$cds_start_lookup
cds_stop_lookup  <- cds_lookups$cds_stop_lookup

# ---- (b): existing goal2-ptc-mechanisms attribution (TD2-based) ----
# Reproduce the relevant pieces of the Rmd chunk on pop_ptc_c2.

fc_c2_cache <- readRDS(file.path(CACHE_DIR, "fc_c2_allsamples.rds"))
fw_c2 <- readRDS(file.path(CACHE_DIR, "fw_c2_allsamples.rds"))

fc_ps <- as.data.table(fc_c2_cache$pair_summary)
# Restrict fc_ps to the effectively_ptc subset (canonical PTC+ under new scope)
ep_ids <- pop_ptc_c2$comparator_isoform_id
fc_ep <- fc_ps[comparator_isoform_id %in% ep_ids]
fc_ep[, has_ptc := comparator_isoform_id %in% ptc$isoform_id[ptc$has_ptc == TRUE]]
frameshift_cats <- c("same_start_frameshift",
                     "diff_start_diff_frame_with_frameshift")
fc_ep[, is_frameshift_cat := frame_category %in% frameshift_cats]

# Split into (b) TD2 PTC+ (original_ptc=TRUE) and (a) recovered (original_ptc=FALSE)
ra_c2_dt <- as.data.table(ra_c2)
orig_lookup <- setNames(ra_c2_dt$original_ptc, ra_c2_dt$comparator_isoform_id)
fc_ep[, original_ptc := orig_lookup[comparator_isoform_id]]
fc_ep_td2  <- fc_ep[original_ptc == TRUE]
fc_ep_recov <- fc_ep[original_ptc == FALSE]
cat(sprintf("  [E/F-NMD] effectively_ptc split: TD2+ %d | recovered %d\n",
            nrow(fc_ep_td2), nrow(fc_ep_recov)))

# Build attribution from (a): ref_atg_analysis attr_mechanism / attr_event
attr_a <- ra_c2_dt[category == "effectively_ptc" & original_ptc == FALSE,
                  .(comparator_isoform_id, reference_isoform_id,
                    attribution = ifelse(!is.na(attr_event) &
                                          attr_event != "Unresolved",
                                          "direct", "unresolved"),
                    ptc_causing_event = attr_event,
                    mechanism = attr_mechanism)]
cat(sprintf("  [E/F-NMD] (a) recovered attributed: %d directly + %d unresolved\n",
            sum(attr_a$attribution == "direct"),
            sum(attr_a$attribution == "unresolved")))

# Build attribution from (b): run attribute_ptc_events on fc_ep_td2
if (nrow(fc_ep_td2) > 0) {
  ptc_genomic_pos <- ifelse(strand_lookup[fc_ep_td2$comparator_isoform_id] == "-",
                            cds_start_lookup[fc_ep_td2$comparator_isoform_id],
                            cds_stop_lookup[fc_ep_td2$comparator_isoform_id])
  names(ptc_genomic_pos) <- fc_ep_td2$comparator_isoform_id
  is_fs_lookup <- setNames(fc_ep_td2$is_frameshift_cat, fc_ep_td2$comparator_isoform_id)
  attr_b <- attribute_ptc_events(
    pairs = fc_ep_td2[, .(comparator_isoform_id, reference_isoform_id)],
    fw_events = fw_c2$events,
    profiles = as.data.frame(pop_ptc_c2),
    ptc_genomic_pos = ptc_genomic_pos,
    atg_genomic_pos = NULL,
    strand_vec = strand_lookup,
    is_frameshift_vec = is_fs_lookup
  )
  attr_b <- as.data.table(attr_b)
  cat(sprintf("  [E/F-NMD] (b) TD2+ attributed: %d directly + %d unresolved\n",
              sum(attr_b$attribution == "direct"),
              sum(attr_b$attribution == "unresolved")))
} else {
  attr_b <- data.table()
}

# Same-stop 3'UTR splice handling on (b) subset
if (nrow(fc_ep_td2) > 0) {
  fc_ep_td2[, strand_s := strand_lookup[comparator_isoform_id]]
  fc_ep_td2[, comp_stop_s := ifelse(strand_s == "-",
                                    cds_start_lookup[comparator_isoform_id],
                                    cds_stop_lookup[comparator_isoform_id])]
  fc_ep_td2[, ref_stop_s := ifelse(strand_s == "-",
                                   cds_start_lookup[reference_isoform_id],
                                   cds_stop_lookup[reference_isoform_id])]
  fc_ep_td2[, same_stop_s := comp_stop_s == ref_stop_s]
  same_stop_b <- fc_ep_td2[same_stop_s == TRUE]
  if (nrow(same_stop_b) > 0) {
    ss_stop_vec <- setNames(same_stop_b$comp_stop_s, same_stop_b$comparator_isoform_id)
    ss_strand_vec <- setNames(same_stop_b$strand_s, same_stop_b$comparator_isoform_id)
    ss_events_b <- attribute_3utr_splice(
      pairs = same_stop_b[, .(comparator_isoform_id, reference_isoform_id)],
      profiles = as.data.frame(pop_ptc_c2),
      stop_genomic_pos = ss_stop_vec,
      strand_vec = ss_strand_vec
    )
    ss_events_b <- as.data.table(ss_events_b)[attribution == "direct"]
    cat(sprintf("  [E/F-NMD] (b) same-stop 3'UTR-splice attributions: %d\n",
                nrow(ss_events_b)))
  } else {
    ss_events_b <- data.table()
  }
} else {
  ss_events_b <- data.table()
}

# Combine (a) + (b) attribution. For (b), exclude same-stop pairs from diff-stop
# direct attribution (per Rmd convention) and add the same-stop 3'UTR-splice records.
if (nrow(attr_b) > 0) {
  same_stop_ids_b <- ss_events_b$comparator_isoform_id %||% character(0)
  diffstop_b <- attr_b[attribution == "direct" &
                       !comparator_isoform_id %in% same_stop_ids_b]
  attr_b_final <- rbind(diffstop_b, ss_events_b, fill = TRUE)
} else {
  attr_b_final <- data.table()
}

all_attr_new <- rbind(attr_a, attr_b_final, fill = TRUE)
all_attr_new <- all_attr_new[attribution == "direct"]
all_attr_new$event_short <- shorten_event_labels(all_attr_new$ptc_causing_event,
                                                  collapse_rare = 3)
cat(sprintf("  [E/F-NMD] all_attr_new (combined direct attribution): %d pairs\n",
            nrow(all_attr_new)))

# ---- Panel E: PTC-causing event proportions vs Control baseline ----
n_ptc_attr <- nrow(all_attr_new)
ptc_evt_freq <- sort(table(all_attr_new$ptc_causing_event), decreasing = TRUE)
ptc_evt_df <- data.table(event_type = names(ptc_evt_freq),
                         n_ptc_events = as.integer(ptc_evt_freq),
                         pct_of_ptc = round(100 * as.integer(ptc_evt_freq) / n_ptc_attr, 1))

# Control baseline: all detailed_events in pop_trace_c4
ctrl_all_events <- do.call(rbind, lapply(pop_trace_c4$detailed_events, as.data.frame))
ctrl_event_freq <- table(ctrl_all_events$event_type)
ctrl_total <- sum(ctrl_event_freq)
cat(sprintf("  [E] Control baseline events: %d (across %d pop_trace_c4 pairs)\n",
            ctrl_total, nrow(pop_trace_c4)))

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
