#!/usr/bin/env Rscript
# Paper Figure: PTCs in NMD
# Panels: A) Stop codon distance  B) NMD response vs distance
#         C) PTC event proportions vs control  D) PTC-causing events stacked bar

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
source("visualization_functions.R")
source("analysis_functions.R")
library(ggplot2)
library(patchwork)
library(dplyr)
library(Isopair)

paper_theme <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        plot.margin = margin(5, 8, 5, 5))

pal <- c("NMD" = "#ef8a62", "Control" = "#67a9cf")

# ---- Load data ----
profiles_c2 <- readRDS("data_mashr/profiles_c2_allsamples.rds")
profiles_c4 <- readRDS("data_mashr/profiles_c4_allsamples.rds")
cds <- readRDS("data_mashr/cds.rds")
ptc_data <- readRDS("data_mashr/ptc.rds")
structures <- readRDS("data_mashr/structures.rds")
fw_data <- readRDS("data_mashr/analysis_cache/fw_c2_allsamples.rds")

# Gene-match
shared <- merge(
  unique(profiles_c2[, c("gene_id", "reference_isoform_id")]),
  unique(profiles_c4[, c("gene_id", "reference_isoform_id")]),
  by = c("gene_id", "reference_isoform_id"))
gm_c2 <- merge(profiles_c2, shared, by = c("gene_id", "reference_isoform_id"))
gm_c4 <- merge(profiles_c4, shared, by = c("gene_id", "reference_isoform_id"))

coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
coding_c2 <- gm_c2[gm_c2$reference_isoform_id %in% coding_ids &
                      gm_c2$comparator_isoform_id %in% coding_ids, ]
coding_c4 <- gm_c4[gm_c4$reference_isoform_id %in% coding_ids &
                      gm_c4$comparator_isoform_id %in% coding_ids, ]

# PTC status
ptc_lookup <- setNames(ptc_data$has_ptc, ptc_data$isoform_id)
ptc_dist_lookup <- setNames(ptc_data$ptc_distance, ptc_data$isoform_id)

# ---- Panel A: Stop codon distance density ----
dist_data <- rbind(
  data.frame(comparison = "NMD",
             distance = ptc_dist_lookup[coding_c2$comparator_isoform_id]),
  data.frame(comparison = "Control",
             distance = ptc_dist_lookup[coding_c4$comparator_isoform_id]))
dist_data <- dist_data[!is.na(dist_data$distance), ]

pA <- ggplot(dist_data, aes(x = distance, fill = comparison)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey40") +
  annotate("text", x = 60, y = Inf, vjust = 2, label = "PTC threshold (50 nt)",
           size = 3, color = "grey40", fontface = "italic") +
  scale_fill_manual(values = pal) +
  coord_cartesian(xlim = c(-1500, 2000)) +
  labs(title = "Stop Codon Distance to Last EJC",
       subtitle = "Dashed line = 50 nt PTC threshold",
       x = "Distance from stop codon to last EJC (nt)",
       y = "Density", fill = "Comparison") +
  paper_theme +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = NA))

# ---- Panel B: NMD response vs distance (scatter) ----
# Use mashr posterior means (same as report)
mashr_model <- readRDS("/Users/petecastaldi/claude_projects/nmd/isocall_dge/mashr/mashr_isoform_model_2026.3.10.rds")
pm <- ashr::get_pm(mashr_model)
mean_logfc_mashr <- data.frame(
  isoform_id = rownames(pm),
  mean_logFC = rowMeans(pm),
  stringsAsFactors = FALSE)

nmd_comp_ids <- coding_c2$comparator_isoform_id
scatter_df <- data.frame(
  isoform_id = nmd_comp_ids,
  distance = ptc_dist_lookup[nmd_comp_ids],
  has_ptc = ptc_lookup[nmd_comp_ids])
scatter_df <- merge(scatter_df, mean_logfc_mashr, by = "isoform_id")
scatter_df <- scatter_df[!is.na(scatter_df$distance), ]

# Trim to 1st-99th percentile for x-axis (matching report)
xlims <- quantile(scatter_df$distance, c(0.01, 0.99), na.rm = TRUE)

pB <- ggplot(scatter_df, aes(x = distance, y = mean_logFC)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_smooth(method = "loess", se = TRUE, color = "steelblue") +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_x_continuous(limits = xlims) +
  labs(title = "NMD Response vs Stop Codon Distance",
       subtitle = sprintf("n = %s NMD comparators | mashr logFC",
                          format(nrow(scatter_df), big.mark = ",")),
       x = "Distance from stop codon to last EJC (nt)",
       y = "Mean logFC (mashr)") +
  paper_theme

# ---- Panel C: PTC event proportions vs control baseline ----
# PTC-causing events
ptc_plus_c2 <- coding_c2[ptc_lookup[coding_c2$comparator_isoform_id] == TRUE, ]
cds_strand <- setNames(cds$strand, cds$isoform_id)
ptc_stop <- ifelse(cds_strand[ptc_plus_c2$comparator_isoform_id] == "-",
  setNames(cds$cds_start, cds$isoform_id)[ptc_plus_c2$comparator_isoform_id],
  setNames(cds$cds_stop, cds$isoform_id)[ptc_plus_c2$comparator_isoform_id])
names(ptc_stop) <- ptc_plus_c2$comparator_isoform_id
strand_vec <- cds_strand[ptc_plus_c2$comparator_isoform_id]
ref_stop_c <- ifelse(cds_strand[ptc_plus_c2$reference_isoform_id] == "-",
  setNames(cds$cds_start, cds$isoform_id)[ptc_plus_c2$reference_isoform_id],
  setNames(cds$cds_stop, cds$isoform_id)[ptc_plus_c2$reference_isoform_id])
same_stop <- ptc_stop == ref_stop_c

# Diff-stop attribution
cds_exons <- readRDS("data_mashr/cds_exons.rds")
fc <- compareIsoformFrames(gm_c2, cds_exons, cds)
fs_cats <- c("same_start_frameshift", "diff_start_diff_frame_with_frameshift")
is_fs <- setNames(fc$pair_summary$frame_category %in% fs_cats,
                    fc$pair_summary$comparator_isoform_id)

diffstop <- ptc_plus_c2[!same_stop & !is.na(same_stop), ]
diffstop_attr <- attributePtcEvents(
  pairs = diffstop[, c("comparator_isoform_id", "reference_isoform_id")],
  fw_events = fw_data$events, profiles = gm_c2,
  ptc_genomic_pos = ptc_stop[diffstop$comparator_isoform_id],
  strand_vec = strand_vec[diffstop$comparator_isoform_id],
  is_frameshift_vec = is_fs[diffstop$comparator_isoform_id])

# Same-stop attribution
samestop <- ptc_plus_c2[same_stop & !is.na(same_stop), ]
samestop_attr <- attribute3UtrSplice(
  pairs = samestop[, c("comparator_isoform_id", "reference_isoform_id")],
  profiles = gm_c2,
  stop_genomic_pos = ptc_stop[samestop$comparator_isoform_id],
  strand_vec = strand_vec[samestop$comparator_isoform_id])

# Combine attributed events
all_attr <- rbind(
  diffstop_attr[diffstop_attr$attribution == "direct", c("mechanism", "ptc_causing_event")],
  samestop_attr[samestop_attr$attribution == "direct", c("mechanism", "ptc_causing_event")])
all_attr$event_short <- shorten_event_labels(all_attr$ptc_causing_event, collapse_rare = 3)

# Control baseline: all events
ctrl_events <- do.call(rbind, lapply(seq_len(nrow(coding_c4)), function(i) {
  de <- coding_c4$detailed_events[[i]]
  if (is.null(de) || nrow(de) == 0) return(NULL)
  de[, "event_type", drop = FALSE]
}))
ctrl_events$event_short <- shorten_event_labels(ctrl_events$event_type)

# Compute proportions
n_ptc <- nrow(all_attr); n_ctrl <- nrow(ctrl_events)
ptc_freq <- as.data.frame(table(all_attr$event_short), stringsAsFactors = FALSE)
names(ptc_freq) <- c("event_type", "n")
ptc_freq$pct <- round(100 * ptc_freq$n / n_ptc, 1)
ptc_freq$group <- "PTC-causing (NMD)"

ctrl_freq <- as.data.frame(table(ctrl_events$event_short), stringsAsFactors = FALSE)
names(ctrl_freq) <- c("event_type", "n")
ctrl_freq$pct <- round(100 * ctrl_freq$n / n_ctrl, 1)
ctrl_freq$group <- "All events (Control)"

bar_ptc <- rbind(ptc_freq, ctrl_freq)
evt_order_ptc <- ptc_freq$event_type[order(-ptc_freq$pct)]
bar_ptc$event_type <- factor(bar_ptc$event_type, levels = evt_order_ptc)
bar_ptc$group <- factor(bar_ptc$group,
  levels = c("PTC-causing (NMD)", "All events (Control)"))

# Significance annotations (Fisher's exact test per event type)
ptc_fisher <- data.frame(event_type = ptc_freq$event_type, stringsAsFactors = FALSE)
ptc_fisher$fisher_p <- sapply(seq_len(nrow(ptc_fisher)), function(i) {
  et <- ptc_fisher$event_type[i]
  a <- ptc_freq$n[ptc_freq$event_type == et]
  b <- n_ptc - a
  c_val <- ctrl_freq$n[ctrl_freq$event_type == et]
  if (length(c_val) == 0) c_val <- 0L
  d_val <- n_ctrl - c_val
  ft <- tryCatch(fisher.test(matrix(c(a, b, c_val, d_val), nrow = 2)),
                  error = function(e) NULL)
  if (!is.null(ft)) ft$p.value else NA
})
ptc_fisher$sig_label <- ifelse(ptc_fisher$fisher_p < 0.001, "***",
  ifelse(ptc_fisher$fisher_p < 0.01, "**",
    ifelse(ptc_fisher$fisher_p < 0.05, "*", "ns")))
ptc_fisher$event_type <- factor(ptc_fisher$event_type, levels = evt_order_ptc)

# y position for significance labels
ptc_fisher$y_pos <- sapply(ptc_fisher$event_type, function(et) {
  vals <- bar_ptc$pct[bar_ptc$event_type == et]
  if (length(vals) == 0) return(2) else max(vals, na.rm = TRUE) + 2
})

pC <- ggplot(bar_ptc, aes(x = event_type, y = pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_text(data = ptc_fisher, aes(x = event_type, y = y_pos, label = sig_label),
            inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("PTC-causing (NMD)" = "#ef8a62",
                                "All events (Control)" = "#67a9cf")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "PTC-Causing Event Proportions: NMD versus Control Pairs",
       x = NULL, y = "% of events", fill = NULL) +
  paper_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom")

# ---- Panel D: PTC-causing events stacked bar ----
mech_colors <- c("Frameshift" = "#d6604d",
                  "In-frame stop" = "#4393c3",
                  "3'UTR splice" = "#66c2a5")

ptc_agg <- all_attr %>%
  group_by(event_short, mechanism) %>%
  summarise(freq = n(), .groups = "drop")

evt_totals <- ptc_agg %>%
  group_by(event_short) %>%
  summarise(total = sum(freq), .groups = "drop") %>%
  arrange(total)
ptc_agg$event_short <- factor(ptc_agg$event_short, levels = evt_totals$event_short)
ptc_agg$mechanism <- factor(ptc_agg$mechanism,
  levels = c("Frameshift", "In-frame stop", "3'UTR splice"))

pD <- ggplot(ptc_agg, aes(x = freq, y = event_short, fill = mechanism)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = ifelse(freq >= 10, freq, "")),
            position = position_stack(vjust = 0.5), size = 3.5,
            color = "white", fontface = "bold") +
  scale_fill_manual(values = mech_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "PTC-Causing Splice Events",
       subtitle = sprintf("%s attributed pairs",
                          format(sum(ptc_agg$freq), big.mark = ",")),
       x = "Number of pairs", y = NULL, fill = "Mechanism") +
  paper_theme +
  theme(legend.position = "bottom")

# ---- Compose ----
combined <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A",
    theme = theme(plot.margin = margin(10, 10, 10, 10))) &
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave("paper_figures/FigX-PTCsinNMD_generated.png", combined,
       width = 16, height = 9, dpi = 300, bg = "white")
ggsave("paper_figures/FigX-PTCsinNMD_generated.pdf", combined,
       width = 16, height = 9, bg = "white")

cat("Saved: paper_figures/FigX-PTCsinNMD_generated.png and .pdf\n")
