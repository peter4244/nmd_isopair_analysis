#!/usr/bin/env Rscript
# =============================================================================
# verify_section1_claim113.R — identify the "additional expression filter" in
# §1 claim 1.13.
#
# Claim (manuscript 2026.7.17, §1): "Only 27.8% (LAE) to 51.7% (MV) of cell
# type-restricted isoforms passed the additional expression filter required for
# differential testing, but among those that did all reached significance."
#
# 1.12's restricted-isoform COUNTS are already verified exact (LAE 8,087 ·
# AT2 2,267 · MV 2,015 · FB 1,131). The open question is which filter turns those
# into 27.8% / 51.7%. The 2026-07-21 attempt used the DMSO-only filterByExpr set
# (105,938) and got LAE 19.5% / MV 31.4% — right ordering, both too low, i.e. the
# real filter is more PERMISSIVE.
#
# Hypothesis under test: the filter for differential testing used ALL samples
# (DMSO + Smg1i), not DMSO only. Restricted isoforms are disproportionately NMD
# targets, which are low in DMSO and rescued under Smg1i — so including Smg1i
# should pass substantially more of them.
#
# Targets: LAE 27.8% of 8,087 = ~2,248 · MV 51.7% of 2,015 = ~1,042
# =============================================================================
suppressPackageStartupMessages({ library(edgeR); library(dplyr); library(tibble) })
select <- dplyr::select; filter <- dplyr::filter

REPO <- path.expand("~/claude_projects/nmd")
dge  <- readRDS(file.path(REPO, "nmd_fig_data/dge_isoform_longread_filtered_2026.3.3.rds"))

samp <- dge$samples %>% rownames_to_column("bam_id") %>%
  mutate(ct = as.character(ct), treatment = as.character(treatment))
samp$ct[samp$ct == "AT"] <- "AT2"; samp$ct[samp$ct == "DD"] <- "LAE"
CT <- c("AT2","LAE","FB","MV")
samp <- samp %>% filter(ct %in% CT)
cnt  <- dge$counts[, samp$bam_id]
cat(sprintf("filtered DGEList: %d isoforms x %d samples (4 CTs)\n", nrow(cnt), ncol(cnt)))

# ---- 1.12: restricted sets (Methods definition; counts already verified) ----
expressed_dmso <- sapply(CT, function(c) {
  ids <- samp$bam_id[samp$ct == c & samp$treatment == "DMSO"]
  rowSums(cnt[, ids, drop = FALSE]) >= 1
})
n_ct_detected <- rowSums(expressed_dmso)
restricted <- lapply(setNames(CT, CT), function(c) rownames(cnt)[expressed_dmso[, c] & n_ct_detected == 1])
cat("\nrestricted isoform counts (1.12 — verified exact previously):\n")
print(sapply(restricted, length))

# ---- candidate filters ------------------------------------------------------
pass_rate <- function(keep_ids, label) {
  r <- sapply(CT, function(c) 100 * mean(restricted[[c]] %in% keep_ids))
  data.frame(filter = label, n_pass = length(keep_ids),
             AT2 = round(r["AT2"],1), LAE = round(r["LAE"],1),
             FB = round(r["FB"],1), MV = round(r["MV"],1),
             row.names = NULL)
}

res <- list()

# (A) DMSO-only filterByExpr — the 2026-07-21 attempt
idx <- samp$treatment == "DMSO"
y <- DGEList(cnt[, idx]); y$samples$group <- factor(samp$ct[idx])
k <- filterByExpr(y, design = model.matrix(~0 + factor(samp$ct[idx])), min.count = 5, min.total.count = 10)
res[[1]] <- pass_rate(rownames(cnt)[k], "A. DMSO-only filterByExpr (min.count=5)")

# (B) ALL samples, ~0+ct  <-- hypothesis
y2 <- DGEList(cnt); y2$samples$group <- factor(samp$ct)
k2 <- filterByExpr(y2, design = model.matrix(~0 + factor(samp$ct)), min.count = 5, min.total.count = 10)
res[[2]] <- pass_rate(rownames(cnt)[k2], "B. ALL samples filterByExpr (min.count=5)")

# (C) ALL samples, edgeR defaults
k3 <- filterByExpr(y2, design = model.matrix(~0 + factor(samp$ct)))
res[[3]] <- pass_rate(rownames(cnt)[k3], "C. ALL samples filterByExpr (defaults)")

# (D) ALL samples, group = ct x treatment
grp <- factor(paste(samp$ct, samp$treatment, sep = "_"))
k4 <- filterByExpr(DGEList(cnt, group = grp), group = grp)
res[[4]] <- pass_rate(rownames(cnt)[k4], "D. filterByExpr, group = ct x treatment")

# (E) per-CT filterByExpr within that CT's samples only
k5 <- unique(unlist(lapply(CT, function(c) {
  ids <- samp$bam_id[samp$ct == c]
  g <- factor(samp$treatment[samp$ct == c])
  rownames(cnt)[filterByExpr(DGEList(cnt[, ids], group = g), group = g)]
})))
res[[5]] <- pass_rate(k5, "E. per-CT filterByExpr (union)")

# (F) simple CPM thresholds
cpm_all <- cpm(DGEList(cnt))
for (thr in c(0.25, 1)) for (ns in c(2, 3)) {
  keep <- rownames(cnt)[rowSums(cpm_all >= thr) >= ns]
  res[[length(res)+1]] <- pass_rate(keep, sprintf("F. CPM>=%.2f in >=%d samples", thr, ns))
}

out <- bind_rows(res)
cat("\n=========== CANDIDATE FILTERS vs CLAIM ===========\n")
cat("claim: LAE 27.8%% (lowest) -> MV 51.7%% (highest)\n\n")
print(out, row.names = FALSE)

out$err <- abs(out$LAE - 27.8) + abs(out$MV - 51.7)
best <- out[which.min(out$err), ]
cat(sprintf("\nclosest match: %s  (LAE %.1f vs 27.8, MV %.1f vs 51.7, |err| %.1f)\n",
            best$filter, best$LAE, best$MV, best$err))
cat("\nordering check (claim: LAE lowest, MV highest):\n")
for (i in seq_len(nrow(out))) {
  v <- unlist(out[i, c("AT2","LAE","FB","MV")])
  cat(sprintf("  %-42s min=%s max=%s %s\n", out$filter[i], names(which.min(v)), names(which.max(v)),
              ifelse(names(which.min(v))=="LAE" && names(which.max(v))=="MV", "<-- ordering matches", "")))
}
