#!/usr/bin/env Rscript
# Conceptual figure: isoform pair framework for splice event detection
# Output: figures/concept_isoform_pairs.png / .pdf

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
library(ggplot2)

# ---- Gene model coordinates ----
# 6-exon gene, + strand, wide intron gaps for clear arrows
ref_starts  <- c(200, 600, 1000, 1400, 1800, 2200)
ref_ends    <- c(350, 800, 1200, 1550, 2000, 2500)

nmd_starts  <- c(200, 600, 1000, 1800, 2200)       # exon 4 skipped
nmd_ends    <- c(350, 800, 1200, 2000, 2500)

ctrl_starts <- c(200, 600, 940, 1400, 1800, 2200)   # A3SS on exon 3
ctrl_ends   <- c(350, 800, 1200, 1550, 2000, 2500)

# ---- Layout ----
ref_y  <- 3
nmd_y  <- 2
ctrl_y <- 1
exon_h <- 0.28

# ---- Exon rectangles ----
build_exons <- function(starts, ends, y, fill) {
  data.frame(xmin = starts, xmax = ends,
             ymin = y - exon_h/2, ymax = y + exon_h/2,
             fill = fill, stringsAsFactors = FALSE)
}

all_exons <- rbind(
  build_exons(ref_starts, ref_ends, ref_y, "#4393c3"),
  build_exons(nmd_starts, nmd_ends, nmd_y, "#ef8a62"),
  build_exons(ctrl_starts, ctrl_ends, ctrl_y, "#67a9cf"))

# ---- Intron lines ----
build_introns <- function(starts, ends, y) {
  n <- length(starts)
  if (n < 2) return(data.frame())
  data.frame(x = ends[-n], xend = starts[-1], y = y, yend = y)
}

all_introns <- rbind(
  build_introns(ref_starts, ref_ends, ref_y),
  build_introns(nmd_starts, nmd_ends, nmd_y),
  build_introns(ctrl_starts, ctrl_ends, ctrl_y))

# ---- Intron chevron arrows ----
build_chevrons <- function(starts, ends, y, arm = 0.07) {
  n <- length(starts)
  if (n < 2) return(data.frame())
  rows <- list()
  for (i in seq_len(n - 1)) {
    intron_s <- ends[i]; intron_e <- starts[i + 1]
    intron_len <- intron_e - intron_s
    n_chev <- max(1L, min(3L, floor(intron_len / 120)))
    step <- intron_len / (n_chev + 1)
    for (j in seq_len(n_chev)) {
      cx <- intron_s + step * j
      w <- min(25, intron_len * 0.12)
      rows[[length(rows) + 1]] <- data.frame(
        x = c(cx - w, cx - w), xend = c(cx, cx),
        y = c(y + arm, y - arm), yend = c(y, y))
    }
  }
  do.call(rbind, rows)
}

all_chevrons <- rbind(
  build_chevrons(ref_starts, ref_ends, ref_y),
  build_chevrons(nmd_starts, nmd_ends, nmd_y),
  build_chevrons(ctrl_starts, ctrl_ends, ctrl_y))

# ---- Pair brackets (ticks point LEFT, toward gene models) ----
bracket_nmd_x  <- 2600
bracket_ctrl_x <- 2780
tick_len <- 40

bracket_segs <- rbind(
  # NMD pair: vertical line spanning ref to NMD
  data.frame(x = bracket_nmd_x, xend = bracket_nmd_x,
             y = nmd_y, yend = ref_y),
  # NMD pair: top tick (pointing left)
  data.frame(x = bracket_nmd_x, xend = bracket_nmd_x - tick_len,
             y = ref_y, yend = ref_y),
  # NMD pair: bottom tick (pointing left)
  data.frame(x = bracket_nmd_x, xend = bracket_nmd_x - tick_len,
             y = nmd_y, yend = nmd_y),
  # Control pair: vertical line spanning ref to ctrl
  data.frame(x = bracket_ctrl_x, xend = bracket_ctrl_x,
             y = ctrl_y, yend = ref_y),
  # Control pair: top tick (pointing left)
  data.frame(x = bracket_ctrl_x, xend = bracket_ctrl_x - tick_len,
             y = ref_y, yend = ref_y),
  # Control pair: bottom tick (pointing left)
  data.frame(x = bracket_ctrl_x, xend = bracket_ctrl_x - tick_len,
             y = ctrl_y, yend = ctrl_y))

# ---- Build the figure ----
fig <- ggplot() +
  # Intron lines
  geom_segment(data = all_introns,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = "grey50", linewidth = 0.5) +
  # Chevron arrows
  geom_segment(data = all_chevrons,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = "grey50", linewidth = 0.5) +
  # Exon rectangles
  geom_rect(data = all_exons,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = all_exons$fill, color = "grey30", linewidth = 0.5) +
  # Brackets
  geom_segment(data = bracket_segs,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = "grey20", linewidth = 0.8) +
  # Track labels
  annotate("text", x = 140, y = ref_y,
           label = "Reference\n(dominant, non-NMD)",
           hjust = 1, size = 4.5, fontface = "bold", color = "#2166ac") +
  annotate("text", x = 140, y = nmd_y,
           label = "NMD comparator\n(SMG1i-responsive)",
           hjust = 1, size = 4.5, fontface = "bold", color = "#d6604d") +
  annotate("text", x = 140, y = ctrl_y,
           label = "Control comparator\n(non-NMD, next-best)",
           hjust = 1, size = 4.5, fontface = "bold", color = "#2b8cbe") +
  # Bracket labels (to the right of brackets)
  annotate("text", x = bracket_nmd_x + 20, y = (ref_y + nmd_y) / 2,
           label = "NMD\npair", size = 4.5, fontface = "bold",
           hjust = 0, color = "#d6604d") +
  annotate("text", x = bracket_ctrl_x + 20, y = (ref_y + ctrl_y) / 2,
           label = "Control\npair", size = 4.5, fontface = "bold",
           hjust = 0, color = "#2b8cbe") +
  # Skipped exon annotation (solid line below, matching report convention)
  annotate("segment", x = 1370, xend = 1580, y = nmd_y - 0.35, yend = nmd_y - 0.35,
           color = "#d6604d", linewidth = 1.2) +
  annotate("text", x = 1475, y = nmd_y - 0.50,
           label = "Skipped exon", size = 3.8, color = "#d6604d",
           fontface = "italic") +
  # A3SS annotation (solid line below, wider than the actual event region
  # so the line renders cleanly)
  annotate("segment", x = 910, xend = 1030, y = ctrl_y - 0.35, yend = ctrl_y - 0.35,
           color = "#2166ac", linewidth = 1.2) +
  annotate("text", x = 970, y = ctrl_y - 0.50,
           label = "A3SS", size = 3.8, color = "#2166ac",
           fontface = "italic") +
  # Title
  labs(title = "Isoform Pair Framework for Splice Event Detection",
       subtitle = "Each gene contributes two paired comparisons sharing the same reference isoform:\nNMD comparator vs Reference and Control comparator vs Reference. Splice events are detected by Isopair and validated through reference reconstruction.") +
  coord_cartesian(xlim = c(-250, 3000), ylim = c(0.2, 3.7)) +
  theme_void(base_size = 14) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(size = 16, face = "bold", hjust = 0,
                                   margin = margin(b = 5)),
        plot.subtitle = element_text(size = 11, hjust = 0, color = "grey30",
                                      margin = margin(b = 10)),
        plot.margin = margin(15, 20, 10, 10))

ggsave("figures/concept_isoform_pairs.png", fig,
       width = 16, height = 6, dpi = 300, bg = "white")
ggsave("figures/concept_isoform_pairs.pdf", fig,
       width = 16, height = 6, bg = "white")

cat("Saved: figures/concept_isoform_pairs.png and .pdf\n")
