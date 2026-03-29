#!/usr/bin/env Rscript
# Paper Figure: PTC Reclassification (NMD+/PTC- Analysis)
# 6-panel figure (3×2 layout):
#   A) 5'UTR Length density
#   B) 5' Boundary Shifts (TSS vs CDS start)
#   C) Reclassification categories (most are effectively PTC)
#   D) TD2 selects longer ORF than truncated ref CDS
#   E) Kozak strength: reference ATG vs TD2
#   F) Splice event proportions: reclassified vs PTC+ vs Control

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
source("visualization_functions.R")
source("analysis_functions.R")
library(ggplot2)
library(patchwork)
library(dplyr)
library(Isopair)

# ---- Consistent theme matching FigX-PTCsinNMD ----
paper_theme <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        plot.margin = margin(5, 8, 5, 5))

fmt <- function(x) format(x, big.mark = ",")

# ---- Load data ----
profiles_c2 <- readRDS("data_mashr/profiles_c2_allsamples.rds")
profiles_c4 <- readRDS("data_mashr/profiles_c4_allsamples.rds")
cds <- readRDS("data_mashr/cds.rds")
ptc_data <- readRDS("data_mashr/ptc.rds")
structures <- readRDS("data_mashr/structures.rds")

cache_dir <- "data_mashr/analysis_cache"
div_c2_all <- readRDS(file.path(cache_dir, "div_c2_allsamples.rds"))
div_c4_all <- readRDS(file.path(cache_dir, "div_c4_allsamples.rds"))
ref_atg <- readRDS(file.path(cache_dir, "ref_atg_analysis.rds"))

# Gene-match
shared <- merge(
  unique(profiles_c2[, c("gene_id", "reference_isoform_id")]),
  unique(profiles_c4[, c("gene_id", "reference_isoform_id")]),
  by = c("gene_id", "reference_isoform_id"))
gm_c2 <- merge(profiles_c2, shared, by = c("gene_id", "reference_isoform_id"))
gm_c4 <- merge(profiles_c4, shared, by = c("gene_id", "reference_isoform_id"))

# Coding pairs
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
coding_c2 <- gm_c2[gm_c2$reference_isoform_id %in% coding_ids &
                      gm_c2$comparator_isoform_id %in% coding_ids, ]
coding_c4 <- gm_c4[gm_c4$reference_isoform_id %in% coding_ids &
                      gm_c4$comparator_isoform_id %in% coding_ids, ]

# PTC status
ptc_lookup <- setNames(ptc_data$has_ptc, ptc_data$isoform_id)
coding_c2$comp_has_ptc <- ptc_lookup[coding_c2$comparator_isoform_id]

ptc_neg_comp_ids <- coding_c2$comparator_isoform_id[!coding_c2$comp_has_ptc]
div_ptcneg <- div_c2_all[div_c2_all$comparator_isoform_id %in% ptc_neg_comp_ids, ]
n_ptc_neg <- length(ptc_neg_comp_ids)

# ---- Panel A: 5'UTR Length density ----
utr5_data <- rbind(
  data.frame(group = "NMD PTC\u2212 comp", utr5 = div_ptcneg$utr5_bp_comp),
  data.frame(group = "NMD PTC\u2212 ref", utr5 = div_ptcneg$utr5_bp_ref),
  data.frame(group = "Control comp", utr5 = div_c4_all$utr5_bp_comp))

pA <- ggplot(utr5_data[utr5_data$utr5 > 0, ], aes(x = utr5, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values = c("NMD PTC\u2212 comp" = "#ef8a62",
                                "NMD PTC\u2212 ref" = "#fddbc7",
                                "Control comp" = "#67a9cf")) +
  labs(title = "5\u2032UTR Length of NMD+/PTC\u2212 Isoforms versus Reference and Control Isoforms",
       subtitle = sprintf("PTC\u2212 comp median: %s bp | Ref: %s bp | Control: %s bp",
                          fmt(median(div_ptcneg$utr5_bp_comp)),
                          fmt(median(div_ptcneg$utr5_bp_ref)),
                          fmt(median(div_c4_all$utr5_bp_comp))),
       x = "5\u2032UTR length (bp, log scale)", y = "Density", fill = NULL) +
  paper_theme +
  theme(legend.position = c(0.80, 0.80),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.4, "cm"))

# ---- Panel B: 5' Boundary Shifts ----
ptc_neg_pairs <- coding_c2[!coding_c2$comp_has_ptc, ]

ref_structs <- structures[match(ptc_neg_pairs$reference_isoform_id, structures$isoform_id), ]
comp_structs <- structures[match(ptc_neg_pairs$comparator_isoform_id, structures$isoform_id), ]
strand_vec <- ref_structs$strand

# TSS: strand-aware 5' end
ref_tss <- ifelse(strand_vec == "+", ref_structs$tx_start, ref_structs$tx_end)
comp_tss <- ifelse(strand_vec == "+", comp_structs$tx_start, comp_structs$tx_end)

# CDS start: strand-aware
ref_cds5 <- ifelse(cds$strand[match(ptc_neg_pairs$reference_isoform_id, cds$isoform_id)] == "+",
                    cds$cds_start[match(ptc_neg_pairs$reference_isoform_id, cds$isoform_id)],
                    cds$cds_stop[match(ptc_neg_pairs$reference_isoform_id, cds$isoform_id)])
comp_cds5 <- ifelse(cds$strand[match(ptc_neg_pairs$comparator_isoform_id, cds$isoform_id)] == "+",
                     cds$cds_start[match(ptc_neg_pairs$comparator_isoform_id, cds$isoform_id)],
                     cds$cds_stop[match(ptc_neg_pairs$comparator_isoform_id, cds$isoform_id)])

# Genomic shifts (positive = extends 5'UTR)
tss_shift <- ifelse(strand_vec == "+",
  as.numeric(ref_tss) - as.numeric(comp_tss),
  as.numeric(comp_tss) - as.numeric(ref_tss))
cds_shift <- ifelse(strand_vec == "+",
  as.numeric(comp_cds5) - as.numeric(ref_cds5),
  as.numeric(ref_cds5) - as.numeric(comp_cds5))

# Control shifts
c4_ref_s <- structures[match(coding_c4$reference_isoform_id, structures$isoform_id), ]
c4_comp_s <- structures[match(coding_c4$comparator_isoform_id, structures$isoform_id), ]
c4_strand <- c4_ref_s$strand
c4_ref_tss <- ifelse(c4_strand == "+", c4_ref_s$tx_start, c4_ref_s$tx_end)
c4_comp_tss <- ifelse(c4_strand == "+", c4_comp_s$tx_start, c4_comp_s$tx_end)
c4_tss_shift <- ifelse(c4_strand == "+",
  as.numeric(c4_ref_tss) - as.numeric(c4_comp_tss),
  as.numeric(c4_comp_tss) - as.numeric(c4_ref_tss))

c4_ref_cds5 <- ifelse(c4_strand == "+",
  cds$cds_start[match(coding_c4$reference_isoform_id, cds$isoform_id)],
  cds$cds_stop[match(coding_c4$reference_isoform_id, cds$isoform_id)])
c4_comp_cds5 <- ifelse(c4_strand == "+",
  cds$cds_start[match(coding_c4$comparator_isoform_id, cds$isoform_id)],
  cds$cds_stop[match(coding_c4$comparator_isoform_id, cds$isoform_id)])
c4_cds_shift <- ifelse(c4_strand == "+",
  as.numeric(c4_comp_cds5) - as.numeric(c4_ref_cds5),
  as.numeric(c4_ref_cds5) - as.numeric(c4_comp_cds5))

genomic_data <- rbind(
  data.frame(group = "NMD PTC\u2212",
             panel = "TSS shift\n(upstream = extends)",
             shift = tss_shift),
  data.frame(group = "Control",
             panel = "TSS shift\n(upstream = extends)",
             shift = c4_tss_shift),
  data.frame(group = "NMD PTC\u2212",
             panel = "CDS start shift\n(downstream = extends)",
             shift = cds_shift),
  data.frame(group = "Control",
             panel = "CDS start shift\n(downstream = extends)",
             shift = c4_cds_shift))

genomic_data$panel <- factor(genomic_data$panel,
  levels = c("TSS shift\n(upstream = extends)",
             "CDS start shift\n(downstream = extends)"))

genomic_data$shift_trunc <- pmin(pmax(genomic_data$shift, -50000), 50000)
n_trunc_genomic <- sum(abs(genomic_data$shift) > 50000, na.rm = TRUE)

pB <- ggplot(genomic_data, aes(x = group, y = shift_trunc, fill = group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~panel, scales = "free_y") +
  scale_fill_manual(values = c("NMD PTC\u2212" = "#ef8a62", "Control" = "#67a9cf")) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "5\u2032 Boundary Shifts for NMD+/PTC\u2212 Isoforms",
       subtitle = sprintf("PTC\u2212 (n=%s) vs Control (n=%s) | Truncated at \u00b150 kb (%d clipped)",
                          fmt(nrow(ptc_neg_pairs)), fmt(nrow(coding_c4)), n_trunc_genomic),
       x = NULL, y = "Genomic shift (bp)", fill = NULL) +
  paper_theme +
  theme(legend.position = "none")

# ---- Panel C: ORF Truncation density (ref ATG traced through comparator) ----
c2_ref <- ref_atg$c2
ptc_neg_ref <- c2_ref[!c2_ref$original_ptc, ]

n_atg_avail <- sum(ptc_neg_ref$ref_atg_exonic_in_comp == TRUE, na.rm = TRUE)
atg_avail <- ptc_neg_ref[ptc_neg_ref$category == "effectively_ptc", ]
atg_avail$orf_ratio <- atg_avail$comp_orf_length_from_ref_atg / atg_avail$ref_orf_length

pC <- ggplot(atg_avail, aes(x = orf_ratio)) +
  geom_density(fill = "#ef8a62", alpha = 0.6) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = median(atg_avail$orf_ratio, na.rm = TRUE),
             color = "#d62728", linewidth = 1) +
  annotate("text", x = median(atg_avail$orf_ratio, na.rm = TRUE), y = Inf,
           label = sprintf("median = %.2f", median(atg_avail$orf_ratio, na.rm = TRUE)),
           vjust = 2, hjust = -0.1, color = "#d62728", size = 3.5) +
  scale_x_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.25)) +
  labs(title = "Most NMD+/PTC\u2212 Have a PTC-Truncated Reference CDS",
       subtitle = sprintf("Effectively-PTC cases (n = %s of %s with ref ATG available) | Ratio < 1 = truncated",
                          fmt(nrow(atg_avail)), fmt(n_atg_avail)),
       x = "ORF length ratio (comparator / reference)", y = "Density") +
  paper_theme

# ---- Panel D: TD2 selects longer ORF ----
eff_ptc <- ptc_neg_ref[ptc_neg_ref$category == "effectively_ptc", ]
eff_ptc_cds <- cds[match(eff_ptc$comparator_isoform_id, cds$isoform_id), ]
eff_ptc$td2_orf_nt <- eff_ptc_cds$orf_length * 3L
td2_longer <- sum(eff_ptc$td2_orf_nt > eff_ptc$comp_orf_length_from_ref_atg, na.rm = TRUE)
pct_td2_longer <- round(100 * td2_longer / nrow(eff_ptc), 1)

length_df <- data.frame(
  ORF = rep(c("TD2-selected CDS", "Reference ATG ORF"), each = nrow(eff_ptc)),
  length_nt = c(eff_ptc$td2_orf_nt, eff_ptc$comp_orf_length_from_ref_atg))

pD <- ggplot(length_df, aes(x = length_nt, fill = ORF)) +
  geom_density(alpha = 0.5) +
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values = c("TD2-selected CDS" = "#67a9cf",
                                "Reference ATG ORF" = "#ef8a62")) +
  labs(title = "TD2 Selects Longer ORF Than Truncated Reference CDS",
       subtitle = sprintf("%s reclassified pairs | TD2 CDS longer in %s%% of cases",
                          fmt(nrow(eff_ptc)), pct_td2_longer),
       x = "ORF length (nt, log scale)", y = "Density", fill = NULL) +
  paper_theme +
  theme(legend.position = c(0.75, 0.80),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.4, "cm"))

# ---- Panel E: Kozak strength comparison ----
# Need to compute Kozak PWM scores for reference ATG vs TD2 ATG
FASTA_PATH <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta"

# Structure lookup
structs_lookup <- setNames(
  lapply(seq_len(nrow(structures)), function(i)
    list(starts = structures$exon_starts[[i]], ends = structures$exon_ends[[i]])),
  structures$isoform_id)

# Genomic to transcript position converter
genomic_to_transcript_local <- function(genomic_pos, exon_starts, exon_ends, strand) {
  if (strand == "-") { exon_starts <- rev(exon_starts); exon_ends <- rev(exon_ends) }
  cum <- 0L
  for (j in seq_along(exon_starts)) {
    es <- exon_starts[j]; ee <- exon_ends[j]
    if (strand == "+") {
      if (genomic_pos >= es && genomic_pos <= ee) return(cum + (genomic_pos - es) + 1L)
    } else {
      if (genomic_pos >= es && genomic_pos <= ee) return(cum + (ee - genomic_pos) + 1L)
    }
    cum <- cum + (ee - es + 1L)
  }
  NA_integer_
}

# Kozak PWM scoring function
kozak_pwm_score <- function(seq_str, atg_tx_pos) {
  if (is.na(atg_tx_pos)) return(NA_real_)
  start <- atg_tx_pos - 6L; end <- atg_tx_pos + 4L
  if (start < 1 || end > nchar(seq_str)) return(NA_real_)
  w <- substr(seq_str, start, end)
  if (nchar(w) != 11 || substr(w, 7, 9) != "ATG") return(NA_real_)
  ch <- strsplit(w, "")[[1]]
  s <- 0
  if (ch[4] %in% c("A","G")) s <- s + 2    # -3
  if (ch[10] == "G") s <- s + 2             # +4
  if (ch[6] == "C") s <- s + 1              # -1
  if (ch[5] == "C") s <- s + 1              # -2
  if (ch[3] == "C") s <- s + 0.5            # -4
  if (ch[2] == "C") s <- s + 0.5            # -5
  if (ch[1] %in% c("A","G")) s <- s + 0.5   # -6
  s
}

# Identify pairs with different ATGs
ptc_neg_avail <- ptc_neg_ref[ptc_neg_ref$ref_atg_exonic_in_comp == TRUE &
                               !is.na(ptc_neg_ref$ref_atg_exonic_in_comp), ]
comp_cds_koz <- cds[match(ptc_neg_avail$comparator_isoform_id, cds$isoform_id), ]
td2_atg_koz <- ifelse(comp_cds_koz$strand == "+", comp_cds_koz$cds_start, comp_cds_koz$cds_stop)
diff_atg_mask <- ptc_neg_avail$ref_atg_genomic != td2_atg_koz & !is.na(td2_atg_koz)
diff_pairs_koz <- ptc_neg_avail[diff_atg_mask, ]
diff_td2_atg <- td2_atg_koz[diff_atg_mask]

# Load comparator sequences
comp_ids_koz <- unique(diff_pairs_koz$comparator_isoform_id)
id_file_koz <- tempfile(fileext = ".txt")
writeLines(comp_ids_koz, id_file_koz)
awk_koz <- sprintf(
  "awk 'BEGIN{while((getline id < \"%s\") > 0) want[id]=1} /^>/{if(seq && found) print hdr \"\\t\" seq; hdr=substr($1,2); found=(hdr in want); seq=\"\"; next} found{seq=seq$0} END{if(seq && found) print hdr \"\\t\" seq}' '%s'",
  id_file_koz, FASTA_PATH)
raw_koz <- system(awk_koz, intern = TRUE)
unlink(id_file_koz)
split_koz <- strsplit(raw_koz, "\t", fixed = TRUE)
seq_koz <- setNames(toupper(vapply(split_koz, `[`, character(1), 2)),
                      vapply(split_koz, `[`, character(1), 1))

# Score Kozak for both ATGs
ref_pwm <- numeric(nrow(diff_pairs_koz))
td2_pwm <- numeric(nrow(diff_pairs_koz))
for (k in seq_len(nrow(diff_pairs_koz))) {
  cid <- diff_pairs_koz$comparator_isoform_id[k]
  strand_k <- diff_pairs_koz$strand[k]
  if (!cid %in% names(seq_koz)) { ref_pwm[k] <- NA; td2_pwm[k] <- NA; next }
  cs <- structs_lookup[[cid]]
  if (is.null(cs)) { ref_pwm[k] <- NA; td2_pwm[k] <- NA; next }
  ref_tx <- genomic_to_transcript_local(diff_pairs_koz$ref_atg_genomic[k],
                                          cs$starts, cs$ends, strand_k)
  td2_tx <- genomic_to_transcript_local(diff_td2_atg[k], cs$starts, cs$ends, strand_k)
  ref_pwm[k] <- kozak_pwm_score(seq_koz[cid], ref_tx)
  td2_pwm[k] <- kozak_pwm_score(seq_koz[cid], td2_tx)
}

koz_valid <- !is.na(ref_pwm) & !is.na(td2_pwm)
n_koz_valid <- sum(koz_valid)
wt_koz <- wilcox.test(ref_pwm[koz_valid], td2_pwm[koz_valid], paired = TRUE)

koz_df <- rbind(
  data.frame(ATG_source = "Reference ATG", score = ref_pwm[koz_valid]),
  data.frame(ATG_source = "TD2-predicted ATG", score = td2_pwm[koz_valid]))
koz_df$ATG_source <- factor(koz_df$ATG_source,
  levels = c("Reference ATG", "TD2-predicted ATG"))

pE <- ggplot(koz_df, aes(x = ATG_source, y = score, fill = ATG_source)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, width = 0.5) +
  scale_fill_manual(values = c("Reference ATG" = "#4393c3",
                                "TD2-predicted ATG" = "#ef8a62")) +
  annotate("text", x = 1.5, y = max(koz_df$score, na.rm = TRUE) + 0.3,
           label = sprintf("paired Wilcoxon p = %s", signif(wt_koz$p.value, 3)),
           size = 3.5) +
  labs(title = "Reference CDS Start Site Kozak Strength Versus TD2 Selected CDS",
       subtitle = sprintf("%s NMD comparators with different ref and TD2 start codons",
                          fmt(n_koz_valid)),
       x = NULL, y = "Kozak PWM score") +
  paper_theme +
  theme(legend.position = "none")

# ---- Panel F: Splice event proportions ----
# Reclassified PTC events
eff_ptc_profiles <- profiles_c2[profiles_c2$comparator_isoform_id %in%
                                   eff_ptc$comparator_isoform_id, ]

nmd_reclass_events <- do.call(rbind, lapply(seq_len(nrow(eff_ptc_profiles)), function(i) {
  de <- eff_ptc_profiles$detailed_events[[i]]
  if (is.null(de) || nrow(de) == 0) return(NULL)
  de[, c("event_type", "direction", "bp_diff")]
}))

# Original PTC+ events
ptc_pos_comp_ids <- coding_c2$comparator_isoform_id[coding_c2$comp_has_ptc]
ptc_pos_profiles <- profiles_c2[profiles_c2$comparator_isoform_id %in% ptc_pos_comp_ids, ]

ptcpos_events <- do.call(rbind, lapply(seq_len(nrow(ptc_pos_profiles)), function(i) {
  de <- ptc_pos_profiles$detailed_events[[i]]
  if (is.null(de) || nrow(de) == 0) return(NULL)
  de[, c("event_type", "direction", "bp_diff")]
}))

# Control events
ctrl_events <- do.call(rbind, lapply(seq_len(nrow(coding_c4)), function(i) {
  de <- coding_c4$detailed_events[[i]]
  if (is.null(de) || nrow(de) == 0) return(NULL)
  de[, c("event_type", "direction", "bp_diff")]
}))

n_nmd_events <- nrow(nmd_reclass_events)
n_ptcpos_events <- nrow(ptcpos_events)
n_ctrl_events <- nrow(ctrl_events)

compute_evt_pct <- function(events, n_total, label) {
  freq <- as.data.frame(table(events$event_type), stringsAsFactors = FALSE)
  names(freq) <- c("event_type", "n")
  freq$pct <- round(100 * freq$n / n_total, 1)
  freq$group <- label
  freq
}

reclass_freq <- compute_evt_pct(nmd_reclass_events, n_nmd_events, "Reclassified PTC")
ptcpos_freq <- compute_evt_pct(ptcpos_events, n_ptcpos_events, "Original PTC+")
ctrl_freq <- compute_evt_pct(ctrl_events, n_ctrl_events, "Control")

all_evt <- rbind(reclass_freq, ptcpos_freq, ctrl_freq)
evt_order <- reclass_freq$event_type[order(-reclass_freq$pct)]
all_evt$event_type <- factor(all_evt$event_type, levels = evt_order)
all_evt$group <- factor(all_evt$group,
  levels = c("Original PTC+", "Reclassified PTC", "Control"))

pF <- ggplot(all_evt, aes(x = event_type, y = pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = c("Original PTC+" = "#fddbc7",
                                "Reclassified PTC" = "#ef8a62",
                                "Control" = "#67a9cf")) +
  labs(title = "Splice Event Proportions: NMD+ Reclassified PTC Similar to NMD+/PTC+",
       subtitle = sprintf("PTC+: %s pairs | Reclassified: %s | Control: %s",
                          fmt(nrow(ptc_pos_profiles)), fmt(nrow(eff_ptc)),
                          fmt(nrow(coding_c4))),
       x = NULL, y = "% of events", fill = NULL) +
  paper_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))

# ---- Compose 3×2 layout with full-width titles ----
# Move titles out of individual panels into full-width title strips
# so long titles can span the entire figure width.

# Store titles/subtitles, then strip them from panels
titles <- list(
  A = pA$labels$title, B = pB$labels$title,
  C = pC$labels$title, D = pD$labels$title,
  E = pE$labels$title, F = pF$labels$title)
subtitles <- list(
  A = pA$labels$subtitle, B = pB$labels$subtitle,
  C = pC$labels$subtitle, D = pD$labels$subtitle,
  E = pE$labels$subtitle, F = pF$labels$subtitle)

strip_title <- function(p) p + labs(title = NULL, subtitle = NULL)
pA <- strip_title(pA); pB <- strip_title(pB)
pC <- strip_title(pC); pD <- strip_title(pD)
pE <- strip_title(pE); pF <- strip_title(pF)

# Helper: full-width title strip for one row (two panels side by side)
make_row_title <- function(tag_l, title_l, sub_l, tag_r, title_r, sub_r) {
  wrap_elements(full = grid::grobTree(
    grid::textGrob(
      label = paste0(tag_l, "   ", gsub("\n", " ", title_l)),
      x = grid::unit(0.02, "npc"), y = grid::unit(0.62, "npc"),
      hjust = 0, gp = grid::gpar(fontsize = 14, fontface = "bold")),
    grid::textGrob(
      label = paste0("     ", gsub("\n", " ", sub_l)),
      x = grid::unit(0.02, "npc"), y = grid::unit(0.12, "npc"),
      hjust = 0, gp = grid::gpar(fontsize = 10, col = "grey30")),
    grid::textGrob(
      label = paste0(tag_r, "   ", gsub("\n", " ", title_r)),
      x = grid::unit(0.52, "npc"), y = grid::unit(0.62, "npc"),
      hjust = 0, gp = grid::gpar(fontsize = 14, fontface = "bold")),
    grid::textGrob(
      label = paste0("     ", gsub("\n", " ", sub_r)),
      x = grid::unit(0.52, "npc"), y = grid::unit(0.12, "npc"),
      hjust = 0, gp = grid::gpar(fontsize = 10, col = "grey30"))
  ))
}

t1 <- make_row_title("A", titles$A, subtitles$A, "B", titles$B, subtitles$B)
t2 <- make_row_title("C", titles$C, subtitles$C, "D", titles$D, subtitles$D)
t3 <- make_row_title("E", titles$E, subtitles$E, "F", titles$F, subtitles$F)

# Use patchwork design matrix: 9 plot slots (3 title strips + 6 panels)
# Each title strip spans full width; each panel row has two panels side by side
design <- c(
  area(1, 1, 1, 2),   # t1: row 1 title (full width)
  area(2, 1, 3, 1),   # pA: row 1 left
  area(2, 2, 3, 2),   # pB: row 1 right
  area(4, 1, 4, 2),   # t2: row 2 title (full width)
  area(5, 1, 6, 1),   # pC: row 2 left
  area(5, 2, 6, 2),   # pD: row 2 right
  area(7, 1, 7, 2),   # t3: row 3 title (full width)
  area(8, 1, 9, 1),   # pE: row 3 left
  area(8, 2, 9, 2)    # pF: row 3 right
)

combined <- t1 + pA + pB + t2 + pC + pD + t3 + pE + pF +
  plot_layout(design = design, heights = c(1, 3, 3, 1, 3, 3, 1, 3, 3))

ggsave("paper_figures/FigX-PTCReclassification_generated.png", combined,
       width = 16, height = 15, dpi = 300, bg = "white")
ggsave("paper_figures/FigX-PTCReclassification_generated.pdf", combined,
       width = 16, height = 15, bg = "white")

cat("Saved: paper_figures/FigX-PTCReclassification_generated.png and .pdf\n")
