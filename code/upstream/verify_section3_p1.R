#!/usr/bin/env Rscript
# =============================================================================
# verify_section3_p1.R — Step-1 factual-accuracy check of Section 3, paragraph 1
# ("percent transcriptional output lost"), manuscript 2026.7.17.
#
# Independent recomputation from the long-read isoform DGEList. Does NOT read any
# cached transcriptional_output/ result. Algorithm ported verbatim from Yul's
# code/upstream/transcriptional_output.Rmd (chunks: load, matched-pairs,
# normalize, compute-frac-lost, global-metric, boxplot-isoform,
# gene-level-output, gene-level-global-metric) with BASE repointed to local paths.
#
# Manuscript claims under test (§3 ¶1):
#   C1 gene-level % output lost:    11.1% (AT2) -- 12.5% (MV)        [SF20]
#   C2 isoform-level % output lost: 14.4% (AT2) -- 18.2% (LAE)       [Fig 2A]
#   C3 per-isoform % lost:          median 60-67%, IQR 30-100%       [SF21]
#   C4 (context, old map 3.5) isoform:gene ratio 1.23--1.52
# =============================================================================
suppressPackageStartupMessages({
  library(edgeR); library(dplyr); library(tibble)
})
select <- dplyr::select; filter <- dplyr::filter; first <- dplyr::first

BASE              <- path.expand("~/claude_projects/nmd")
DGELIST_PATH      <- file.path(BASE, "nmd_fig_data/dge_isoform_longread_2026.3.3.rds")
DGELIST_FILT_PATH <- file.path(BASE, "nmd_fig_data/dge_isoform_longread_filtered_2026.3.3.rds")
CT_KEEP           <- c("LAE", "AT2", "FB", "MV")

cat("=== INPUTS ===\n")
cat("unfiltered DGEList:", DGELIST_PATH, "\n")
stopifnot(file.exists(DGELIST_PATH), file.exists(DGELIST_FILT_PATH))

# ---- load + relabel (Rmd chunk: load) ---------------------------------------
dge <- readRDS(DGELIST_PATH)
cat("DGEList dimensions:", dim(dge), "\n")

samp <- dge$samples %>%
  rownames_to_column("bam_id") %>%
  mutate(ct = as.character(ct), treatment = as.character(treatment), donor_id = id)
samp$ct[samp$ct == "AT"] <- "AT2"
samp$ct[samp$ct == "DD"] <- "LAE"
n_before <- nrow(samp)
samp <- samp %>% filter(ct %in% CT_KEEP)
cat(sprintf("samples: %d -> %d after CT_KEEP filter\n", n_before, nrow(samp)))
cat("cell types:", paste(sort(unique(samp$ct)), collapse=", "), "\n")

iso_info <- dge$genes %>% select(any_of(c("txid","gene_id","hgnc_symbol")))

# ---- matched donor pairs (Rmd chunk: matched-pairs) -------------------------
paired_donors <- samp %>%
  group_by(ct, donor_id) %>%
  summarise(has_smg1i = any(treatment=="Smg1i"), has_dmso = any(treatment=="DMSO"), .groups="drop") %>%
  filter(has_smg1i & has_dmso)
cat("\nmatched donor pairs per CT:\n")
print(paired_donors %>% count(ct, name="n_pairs"))

# ---- plain CPM, no TMM (Rmd chunk: normalize) -------------------------------
cpm_mat <- cpm(dge[, samp$bam_id], normalized.lib.sizes = FALSE, log = FALSE)
cat("\nCPM matrix:", nrow(cpm_mat), "isoforms x", ncol(cpm_mat), "samples\n")

# ---- per-isoform fraction lost (Rmd chunk: compute-frac-lost) ---------------
frac_one <- function(mat, ids) {
  out <- list()
  for (ct_i in sort(unique(paired_donors$ct))) {
    for (d in paired_donors$donor_id[paired_donors$ct == ct_i]) {
      s <- samp$bam_id[samp$ct==ct_i & samp$donor_id==d & samp$treatment=="Smg1i"]
      m <- samp$bam_id[samp$ct==ct_i & samp$donor_id==d & samp$treatment=="DMSO"]
      stopifnot(length(s)==1, length(m)==1)
      v_s <- mat[, s]; v_d <- mat[, m]
      delta <- v_s - v_d; delta_pos <- pmax(delta, 0)
      out[[paste(ct_i,d,sep="_")]] <- data.frame(
        id = ids, ct = ct_i, donor_id = d,
        smg1i_val = v_s, delta_pos = delta_pos,
        frac_lost = ifelse(v_s > 0, delta_pos / v_s, 0),
        stringsAsFactors = FALSE)
    }
  }
  bind_rows(out)
}

results_all <- frac_one(cpm_mat, rownames(cpm_mat))
cat("isoform x ct x donor rows:", nrow(results_all), "\n")

# ---- C2: isoform-level global metric (Rmd chunk: global-metric) -------------
iso_global <- results_all %>%
  group_by(ct) %>%
  summarise(pct_output_lost = 100*sum(delta_pos)/sum(smg1i_val), .groups="drop") %>%
  arrange(desc(pct_output_lost))

# ---- C1: gene-level (Rmd chunks: gene-level-output / -global-metric) --------
iso_to_gene <- iso_info %>% select(txid, gene_id, hgnc_symbol) %>%
  filter(!is.na(gene_id) & gene_id != "")
keep_iso  <- intersect(rownames(cpm_mat), iso_to_gene$txid)
cpm_gene  <- cpm_mat[keep_iso, samp$bam_id]
gene_ids  <- iso_to_gene$gene_id[match(keep_iso, iso_to_gene$txid)]
stopifnot(length(gene_ids) == nrow(cpm_gene), !any(is.na(gene_ids)))
cpm_gene_agg <- rowsum(cpm_gene, group = gene_ids)
cat("gene-level CPM matrix:", nrow(cpm_gene_agg), "genes x", ncol(cpm_gene_agg), "samples\n")

gene_results_all <- frac_one(cpm_gene_agg, rownames(cpm_gene_agg))
gene_global <- gene_results_all %>%
  group_by(ct) %>%
  summarise(pct_output_lost = 100*sum(delta_pos)/sum(smg1i_val), .groups="drop") %>%
  arrange(desc(pct_output_lost))

# ---- C3: per-isoform distribution (Rmd chunk: boxplot-isoform) --------------
dge_filt      <- readRDS(DGELIST_FILT_PATH)
filt_isoforms <- rownames(dge_filt$counts)
cat("filtered isoform set:", length(filt_isoforms), "\n")
plot_data <- results_all %>% filter(frac_lost > 0, id %in% filt_isoforms)
cat("isoform-donor rows after filter:", nrow(plot_data), "\n")
box_stats <- plot_data %>%
  group_by(ct) %>%
  summarise(n = n(),
            Q1     = round(100*quantile(frac_lost,0.25),1),
            median = round(100*median(frac_lost),1),
            Q3     = round(100*quantile(frac_lost,0.75),1),
            IQR    = round(100*IQR(frac_lost),1), .groups="drop")

# ---- report -----------------------------------------------------------------
cat("\n\n=========== RECOMPUTED RESULTS ===========\n")
cat("\n-- C2 ISOFORM-level % output lost (claim: 14.4 AT2 -- 18.2 LAE) --\n")
print(as.data.frame(iso_global %>% mutate(pct_output_lost=round(pct_output_lost,2))))
cat("\n-- C1 GENE-level % output lost (claim: 11.1 AT2 -- 12.5 MV) --\n")
print(as.data.frame(gene_global %>% mutate(pct_output_lost=round(pct_output_lost,2))))
cat("\n-- C4 isoform:gene ratio (old map 3.5 claim: 1.23--1.52) --\n")
cmp <- iso_global %>% select(ct, pct_isoform=pct_output_lost) %>%
  left_join(gene_global %>% select(ct, pct_gene=pct_output_lost), by="ct") %>%
  mutate(ratio = round(pct_isoform/pct_gene, 3),
         pct_isoform=round(pct_isoform,2), pct_gene=round(pct_gene,2))
print(as.data.frame(cmp))
cat("\n-- C3 per-isoform % lost, filtered set (claim: median 60-67, IQR 30-100) --\n")
print(as.data.frame(box_stats))

cat("\n=========== CLAIM CHECK ===========\n")
chk <- function(lbl, got, lo, hi) {
  cat(sprintf("%-46s observed %-14s claimed %s\n", lbl,
              sprintf("%.1f-%.1f", min(got), max(got)), sprintf("%.1f-%.1f", lo, hi)))
}
chk("C1 gene-level range (%)",    gene_global$pct_output_lost, 11.1, 12.5)
chk("C2 isoform-level range (%)", iso_global$pct_output_lost,  14.4, 18.2)
chk("C3 per-isoform median (%)",  box_stats$median,            60,   67)
cat("\nC1 argmin/argmax CT: ", gene_global$ct[which.min(gene_global$pct_output_lost)],
    " / ", gene_global$ct[which.max(gene_global$pct_output_lost)],
    "   (claim: AT2 low, MV high)\n", sep="")
cat("C2 argmin/argmax CT: ", iso_global$ct[which.min(iso_global$pct_output_lost)],
    " / ", iso_global$ct[which.max(iso_global$pct_output_lost)],
    "   (claim: AT2 low, LAE high)\n", sep="")
cat("\nsessionInfo R:", R.version.string, "| edgeR", as.character(packageVersion("edgeR")), "\n")
