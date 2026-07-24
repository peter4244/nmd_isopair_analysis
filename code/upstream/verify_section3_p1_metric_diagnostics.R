#!/usr/bin/env Rscript
# =============================================================================
# verify_section3_p1_metric_diagnostics.R — Step-2 (result correctness) probe of
# the "% transcriptional output lost" metric used in §3 ¶1 / Fig 2A / SF20-21
# and quoted in the Abstract ("NMD degraded 11-18% of transcriptional output").
#
# Companion to verify_section3_p1.R, which confirms the published numbers
# reproduce EXACTLY. This script asks a different question: does the statistic
# measure what its name claims? Three diagnostics:
#
#   D1 ATTRIBUTION  what share of the metric comes from isoforms that are
#                   statistically NMD susceptible (mashr lfsr<0.05 & PM>0)?
#   D2 SYMMETRY     recompute in the reverse direction (DMSO - Smg1i). Under a
#                   closed composition this is provably identical, so the metric
#                   cannot distinguish output "lost" from output "gained".
#   D3 NULL         same statistic between two DMSO samples from DIFFERENT donors
#                   (no treatment). NOTE: cross-donor contrasts also carry genuine
#                   donor biology, so this is an upper bound on a no-treatment
#                   baseline, not a matched technical null.
#
# Run:  Rscript code/upstream/verify_section3_p1_metric_diagnostics.R
# =============================================================================
suppressPackageStartupMessages({
  library(edgeR); library(dplyr); library(tibble); library(data.table)
})
select <- dplyr::select; filter <- dplyr::filter

BASE <- path.expand("~/claude_projects/nmd")
dge  <- readRDS(file.path(BASE, "nmd_fig_data/dge_isoform_longread_2026.3.3.rds"))
samp <- dge$samples %>% rownames_to_column("bam_id") %>%
  mutate(ct=as.character(ct), treatment=as.character(treatment), donor_id=id)
samp$ct[samp$ct=="AT"] <- "AT2"; samp$ct[samp$ct=="DD"] <- "LAE"
samp <- samp %>% filter(ct %in% c("LAE","AT2","FB","MV"))
cpm_mat <- cpm(dge[, samp$bam_id], normalized.lib.sizes=FALSE, log=FALSE)
pd <- samp %>% group_by(ct, donor_id) %>%
  summarise(a=any(treatment=="Smg1i"), b=any(treatment=="DMSO"), .groups="drop") %>% filter(a & b)

# ---- D0: closed-composition check (the mechanism behind D2) -----------------
cs <- colSums(cpm_mat)
cat("=== D0 closed composition ===\n")
cat("all CPM column sums == 1e6? ",
    isTRUE(all.equal(unname(cs), rep(1e6, length(cs)), tolerance=1e-9)), "\n")
cat("=> sum(delta)=0 for every pair => sum(max(d,0)) == sum(max(-d,0)) == 0.5*sum|d|\n")
cat("=> the metric IS the total variation distance between the two CPM vectors.\n\n")

# ---- D1 / D2 / D3 -----------------------------------------------------------
lf <- fread(file.path(BASE,"isocall_dge/mashr/mashr_isoform_lfsr_2026.3.10.csv"))
pm <- fread(file.path(BASE,"isocall_dge/mashr/mashr_isoform_posterior_means_2026.3.10.csv"))
idcol <- names(lf)[1]
ctmap <- c(AT2="Smg1i_in_AT", LAE="Smg1i_in_DD", FB="Smg1i_in_FB", MV="Smg1i_in_MV")

out <- list()
for (ct_i in sort(unique(pd$ct))) {
  col <- ctmap[[ct_i]]
  nmd_ids <- lf[[idcol]][ lf[[col]] < 0.05 & pm[[col]] > 0 ]
  inn <- rownames(cpm_mat) %in% nmd_ids
  ds <- pd$donor_id[pd$ct == ct_i]
  fn <- fd <- rn <- rd <- nn <- 0
  for (d in ds) {
    s <- samp$bam_id[samp$ct==ct_i & samp$donor_id==d & samp$treatment=="Smg1i"]
    m <- samp$bam_id[samp$ct==ct_i & samp$donor_id==d & samp$treatment=="DMSO"]
    delta <- cpm_mat[, s] - cpm_mat[, m]
    fn <- fn + sum(pmax(delta,0));  fd <- fd + sum(cpm_mat[, s])   # forward (published)
    rn <- rn + sum(pmax(-delta,0)); rd <- rd + sum(cpm_mat[, m])   # reverse
    nn <- nn + sum(pmax(delta,0)[inn])                             # NMD-susceptible only
  }
  xn <- xd <- 0; npair <- 0
  if (length(ds) >= 2) for (i in 1:(length(ds)-1)) {
    a <- samp$bam_id[samp$ct==ct_i & samp$donor_id==ds[i]   & samp$treatment=="DMSO"]
    b <- samp$bam_id[samp$ct==ct_i & samp$donor_id==ds[i+1] & samp$treatment=="DMSO"]
    xn <- xn + sum(pmax(cpm_mat[,b]-cpm_mat[,a],0)); xd <- xd + sum(cpm_mat[,b]); npair <- npair+1
  }
  out[[ct_i]] <- data.frame(
    ct = ct_i, n_nmd_iso = length(nmd_ids),
    D_forward_published        = round(100*fn/fd, 2),
    D2_reverse                 = round(100*rn/rd, 2),
    D1_from_nmd_susceptible    = round(100*nn/fd, 2),
    D1_share_of_signal_pct     = round(100*nn/fn, 1),
    D3_null_DMSOvsDMSO         = round(100*xn/xd, 2),
    n_null_pairs = npair)
}
res <- do.call(rbind, out)
cat("=== D1-D3 diagnostics (all values are % of Smg1i total CPM) ===\n")
print(res, row.names = FALSE)

cat("\n=== READING ===\n")
cat("D2 == D_forward exactly  -> metric is symmetric; 'lost' is not established by the statistic.\n")
cat("D1_share  ", paste0(range(res$D1_share_of_signal_pct), collapse="-"),
    "% of the signal traces to NMD-susceptible isoforms; the remainder is sub-threshold or noise.\n")
cat("D3 vs forward: null is LARGER in", sum(res$D3_null_DMSOvsDMSO > res$D_forward_published),
    "of", nrow(res), "cell types (caveat: cross-donor null includes donor biology).\n")
