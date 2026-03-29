#!/usr/bin/env Rscript
# Paper Figure: Splice Events in NMD
# Multipanel figure with consistent font sizes
# Panels: A) concept diagram  B) sequence shared  C) event prevalence  D) gain direction

setwd("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")
source("visualization_functions.R")
library(ggplot2)
library(patchwork)
library(dplyr)

# ---- Consistent theme for all panels ----
paper_theme <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        plot.margin = margin(5, 8, 5, 5))

pal <- c("NMD" = "#ef8a62", "Control" = "#67a9cf")

# ---- Panel A: Concept diagram ----
# Reuse the concept figure but with adjusted fonts for multipanel
ref_starts  <- c(200, 600, 1000, 1400, 1800, 2200)
ref_ends    <- c(350, 800, 1200, 1550, 2000, 2500)
nmd_starts  <- c(200, 600, 1000, 1800, 2200)
nmd_ends    <- c(350, 800, 1200, 2000, 2500)
ctrl_starts <- c(200, 600, 940, 1400, 1800, 2200)
ctrl_ends   <- c(350, 800, 1200, 1550, 2000, 2500)

ref_y <- 3; nmd_y <- 2; ctrl_y <- 1; exon_h <- 0.28

build_exons <- function(s, e, y, fill) {
  data.frame(xmin=s, xmax=e, ymin=y-exon_h/2, ymax=y+exon_h/2, fill=fill)
}
build_introns <- function(s, e, y) {
  n <- length(s); if(n<2) return(data.frame())
  data.frame(x=e[-n], xend=s[-1], y=y, yend=y)
}
build_chevrons <- function(s, e, y, arm=0.07) {
  n <- length(s); if(n<2) return(data.frame())
  rows <- list()
  for(i in seq_len(n-1)) {
    is <- e[i]; ie <- s[i+1]; il <- ie-is
    nc <- max(1L,min(3L,floor(il/120))); step <- il/(nc+1)
    for(j in seq_len(nc)) {
      cx <- is+step*j; w <- min(25,il*0.12)
      rows[[length(rows)+1]] <- data.frame(
        x=c(cx-w,cx-w),xend=c(cx,cx),y=c(y+arm,y-arm),yend=c(y,y))
    }
  }
  do.call(rbind, rows)
}

all_exons <- rbind(build_exons(ref_starts,ref_ends,ref_y,"#4393c3"),
                    build_exons(nmd_starts,nmd_ends,nmd_y,"#ef8a62"),
                    build_exons(ctrl_starts,ctrl_ends,ctrl_y,"#67a9cf"))
all_introns <- rbind(build_introns(ref_starts,ref_ends,ref_y),
                      build_introns(nmd_starts,nmd_ends,nmd_y),
                      build_introns(ctrl_starts,ctrl_ends,ctrl_y))
all_chevrons <- rbind(build_chevrons(ref_starts,ref_ends,ref_y),
                       build_chevrons(nmd_starts,nmd_ends,nmd_y),
                       build_chevrons(ctrl_starts,ctrl_ends,ctrl_y))

bracket_nmd_x <- 2600; bracket_ctrl_x <- 2800; tick_len <- 40
bracket_segs <- rbind(
  data.frame(x=bracket_nmd_x,xend=bracket_nmd_x,y=nmd_y,yend=ref_y),
  data.frame(x=bracket_nmd_x,xend=bracket_nmd_x-tick_len,y=ref_y,yend=ref_y),
  data.frame(x=bracket_nmd_x,xend=bracket_nmd_x-tick_len,y=nmd_y,yend=nmd_y),
  data.frame(x=bracket_ctrl_x,xend=bracket_ctrl_x,y=ctrl_y,yend=ref_y),
  data.frame(x=bracket_ctrl_x,xend=bracket_ctrl_x-tick_len,y=ref_y,yend=ref_y),
  data.frame(x=bracket_ctrl_x,xend=bracket_ctrl_x-tick_len,y=ctrl_y,yend=ctrl_y))

pA <- ggplot() +
  geom_segment(data=all_introns, aes(x=x,xend=xend,y=y,yend=yend),
               color="grey50", linewidth=0.4) +
  geom_segment(data=all_chevrons, aes(x=x,xend=xend,y=y,yend=yend),
               color="grey50", linewidth=0.4) +
  geom_rect(data=all_exons, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill=all_exons$fill, color="grey30", linewidth=0.4) +
  geom_segment(data=bracket_segs, aes(x=x,xend=xend,y=y,yend=yend),
               color="grey20", linewidth=0.6) +
  annotate("text",x=140,y=ref_y,label="Reference\n(dominant, non-NMD)",
           hjust=1,size=3.2,fontface="bold",color="#2166ac") +
  annotate("text",x=140,y=nmd_y,label="NMD comparator\n(SMG1i-responsive)",
           hjust=1,size=3.2,fontface="bold",color="#d6604d") +
  annotate("text",x=140,y=ctrl_y,label="Control comparator\n(non-NMD, next-best)",
           hjust=1,size=3.2,fontface="bold",color="#2b8cbe") +
  annotate("text",x=bracket_nmd_x+20,y=2.5,label="NMD\npair",
           size=3.2,fontface="bold",hjust=0,color="#d6604d") +
  annotate("text",x=bracket_ctrl_x+20,y=2,label="Control\npair",
           size=3.2,fontface="bold",hjust=0,color="#2b8cbe") +
  annotate("segment",x=1370,xend=1580,y=nmd_y-0.35,yend=nmd_y-0.35,
           color="#d6604d",linewidth=0.8) +
  annotate("text",x=1475,y=nmd_y-0.50,label="Skipped exon",
           size=2.8,color="#d6604d",fontface="italic") +
  annotate("segment",x=910,xend=1030,y=ctrl_y-0.35,yend=ctrl_y-0.35,
           color="#2166ac",linewidth=0.8) +
  annotate("text",x=970,y=ctrl_y-0.50,label="A3SS",
           size=2.8,color="#2166ac",fontface="italic") +
  coord_cartesian(xlim=c(-550,3100),ylim=c(0.3,3.5),clip="off") +
  theme_void() +
  theme(plot.background=element_rect(fill="white",color=NA),
        panel.background=element_rect(fill="white",color=NA),
        plot.margin=margin(5,5,5,40))

# ---- Load report data for panels B-D ----
profiles_c2 <- readRDS("data_mashr/profiles_c2_allsamples.rds")
profiles_c4 <- readRDS("data_mashr/profiles_c4_allsamples.rds")
structures <- readRDS("data_mashr/structures.rds")
nmd_class <- readRDS("data_mashr/nmd_classification.rds")

# Gene-match
shared <- merge(
  unique(profiles_c2[,c("gene_id","reference_isoform_id")]),
  unique(profiles_c4[,c("gene_id","reference_isoform_id")]),
  by=c("gene_id","reference_isoform_id"))
gm_c2 <- merge(profiles_c2, shared, by=c("gene_id","reference_isoform_id"))
gm_c4 <- merge(profiles_c4, shared, by=c("gene_id","reference_isoform_id"))

# Load cached divergence
set_c2 <- readRDS("data_mashr/analysis_cache/div_c2_allsamples.rds")
set_c4 <- readRDS("data_mashr/analysis_cache/div_c4_allsamples.rds")

# ---- Panel B: Sequence shared density ----
pct_shared <- rbind(
  data.frame(comparison="NMD", pct_shared=set_c2$pct_shared),
  data.frame(comparison="Control", pct_shared=set_c4$pct_shared))
pct_shared <- pct_shared[!is.na(pct_shared$pct_shared),]

pB <- ggplot(pct_shared, aes(x=pct_shared, fill=comparison)) +
  geom_density(alpha=0.5) +
  scale_fill_manual(values=pal) +
  labs(title="Exonic Sequence Shared with\nReference Isoform",
       x="Sequence Shared (%)", y="Density", fill="Comparison") +
  paper_theme +
  theme(legend.position=c(0.25,0.75),
        legend.background=element_rect(fill="white",color=NA))

# ---- Panel C: Event prevalence ----
source("analysis_functions.R")

compute_event_freq_local <- function(profiles) {
  evt_types <- c("Alt_TSS","Alt_TES","A5SS","A3SS","SE",
                  "Missing_Internal","IR","Partial_IR_5","Partial_IR_3")
  n <- nrow(profiles)
  freq <- data.frame(event_type=evt_types, n_pairs=0L, stringsAsFactors=FALSE)
  for (i in seq_along(evt_types)) {
    et <- evt_types[i]
    if (et == "Alt_TSS") freq$n_pairs[i] <- sum(profiles$tss_changed, na.rm=TRUE)
    else if (et == "Alt_TES") freq$n_pairs[i] <- sum(profiles$tes_changed, na.rm=TRUE)
    else {
      col <- paste0("n_", tolower(et))
      if (col %in% names(profiles)) {
        freq$n_pairs[i] <- sum(profiles[[col]] > 0, na.rm=TRUE)
      } else {
        # Check detailed_events
        freq$n_pairs[i] <- sum(sapply(seq_len(nrow(profiles)), function(j) {
          de <- profiles$detailed_events[[j]]
          if(is.null(de)||nrow(de)==0) return(FALSE)
          any(de$event_type == et)
        }))
      }
    }
  }
  freq$pct <- round(100 * freq$n_pairs / n, 1)
  freq
}

ef_c2 <- compute_event_freq_local(gm_c2)
ef_c4 <- compute_event_freq_local(gm_c4)
n_c2 <- nrow(gm_c2); n_c4 <- nrow(gm_c4)

enrichment <- merge(ef_c2, ef_c4, by="event_type", suffixes=c("_NMD","_Ctrl"))
enrichment$fisher_p <- NA_real_
for(i in seq_len(nrow(enrichment))) {
  a <- enrichment$n_pairs_NMD[i]; b <- n_c2-a
  c <- enrichment$n_pairs_Ctrl[i]; d <- n_c4-c
  ft <- tryCatch(fisher.test(matrix(c(a,b,c,d),nrow=2)),error=function(e)NULL)
  if(!is.null(ft)) enrichment$fisher_p[i] <- ft$p.value
}

bar_df <- rbind(
  data.frame(event_type=enrichment$event_type,comparison="NMD",
             pct=enrichment$pct_NMD),
  data.frame(event_type=enrichment$event_type,comparison="Control",
             pct=enrichment$pct_Ctrl))

evt_order <- enrichment$event_type[order(-enrichment$pct_NMD)]
bar_df$event_type <- factor(bar_df$event_type, levels=evt_order)

sig_df <- data.frame(
  event_type=factor(enrichment$event_type, levels=evt_order),
  y=pmax(enrichment$pct_NMD, enrichment$pct_Ctrl)+5,
  label=ifelse(enrichment$fisher_p<0.001,"***",
         ifelse(enrichment$fisher_p<0.01,"**",
           ifelse(enrichment$fisher_p<0.05,"*","ns"))))

pC <- ggplot(bar_df, aes(x=event_type,y=pct,fill=comparison)) +
  geom_col(position=position_dodge(width=0.7),width=0.6) +
  geom_text(data=sig_df,aes(x=event_type,y=y,label=label),
            inherit.aes=FALSE,size=3.5,fontface="bold") +
  scale_fill_manual(values=pal) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  labs(title="Splicing Event Prevalence",
       x=NULL, y="% of Pairs with Event", fill="Comparison") +
  paper_theme +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=9),
        legend.position="bottom")

# ---- Panel D: Gain direction ----
get_gain_pct <- function(profiles, label) {
  all_events <- do.call(rbind, lapply(seq_len(nrow(profiles)), function(i) {
    de <- profiles$detailed_events[[i]]
    if(is.null(de)||nrow(de)==0) return(NULL)
    de[,c("event_type","direction")]
  }))
  if(is.null(all_events)||nrow(all_events)==0) return(NULL)
  all_events %>%
    group_by(event_type) %>%
    summarise(n=n(), n_gain=sum(direction=="GAIN"),
              pct_gain=round(100*mean(direction=="GAIN"),1), .groups="drop") %>%
    mutate(comparison=label)
}

gl_c2 <- get_gain_pct(gm_c2, "NMD")
gl_c4 <- get_gain_pct(gm_c4, "Control")
gl_all <- rbind(gl_c2, gl_c4) %>% filter(n >= 20)

# Fisher test per event type
gl_wide <- merge(gl_c2[,c("event_type","n","n_gain")],
                  gl_c4[,c("event_type","n","n_gain")],
                  by="event_type", suffixes=c("_NMD","_Ctrl"))
gl_wide$fisher_p <- sapply(seq_len(nrow(gl_wide)), function(i) {
  ft <- tryCatch(fisher.test(matrix(c(gl_wide$n_gain_NMD[i],
    gl_wide$n_NMD[i]-gl_wide$n_gain_NMD[i],
    gl_wide$n_gain_Ctrl[i],
    gl_wide$n_Ctrl[i]-gl_wide$n_gain_Ctrl[i]),nrow=2)),error=function(e)NULL)
  if(!is.null(ft)) ft$p.value else NA
})

evt_order_d <- gl_all$event_type[gl_all$comparison=="NMD"][
  order(-gl_all$pct_gain[gl_all$comparison=="NMD"])]
gl_all$event_type <- factor(gl_all$event_type, levels=evt_order_d)

gl_wide_filt <- gl_wide[gl_wide$event_type %in% levels(gl_all$event_type),]
sig_d <- data.frame(
  event_type=factor(gl_wide_filt$event_type, levels=evt_order_d),
  y=sapply(gl_wide_filt$event_type, function(et) {
    vals <- gl_all$pct_gain[gl_all$event_type==et]
    if(length(vals)==0) return(5) else max(vals,na.rm=TRUE)+3
  }),
  label=ifelse(gl_wide_filt$fisher_p<0.001,"***",
         ifelse(gl_wide_filt$fisher_p<0.01,"**",
           ifelse(gl_wide_filt$fisher_p<0.05,"*","ns"))))

pD <- ggplot(gl_all, aes(x=event_type,y=pct_gain,fill=comparison)) +
  geom_col(position=position_dodge(width=0.7),width=0.6) +
  geom_text(data=sig_d,aes(x=event_type,y=y,label=label),
            inherit.aes=FALSE,size=3.5,fontface="bold") +
  scale_fill_manual(values=pal) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(title="GAIN Direction by Event Type",
       x=NULL, y="% Events with GAIN Direction", fill="Comparison") +
  paper_theme +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=9),
        legend.position="bottom")

# ---- Compose multipanel figure ----
combined <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A",
    theme=theme(plot.margin=margin(10,10,10,10))) &
  theme(plot.tag=element_text(size=18,face="bold"))

ggsave("paper_figures/FigX_SpliceEvents_generated.png", combined,
       width=16, height=9, dpi=300, bg="white")
ggsave("paper_figures/FigX_SpliceEvents_generated.pdf", combined,
       width=16, height=9, bg="white")

cat("Saved: paper_figures/FigX_SpliceEvents.png and .pdf\n")
