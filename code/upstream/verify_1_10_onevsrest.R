suppressMessages({library(edgeR); library(limma)})
setwd("/Users/petecastaldi/claude_projects/nmd")
CT_KEEP <- c("AT2","LAE","FB","MV")
y_iso <- readRDS("nmd_fig_data/dge_isoform_longread_filtered_2026.3.3.rds")
keep_s <- y_iso$samples$treatment=="DMSO" & y_iso$samples$ct %in% CT_KEEP
y <- y_iso[, keep_s]
ct  <- droplevels(factor(as.character(y$samples$ct), levels=CT_KEEP))
don <- droplevels(factor(as.character(y$samples$id)))
cat("samples:", ncol(y), " donors:", nlevels(don), "\n")
design <- model.matrix(~ 0 + ct); colnames(design) <- gsub("^ct","",colnames(design))
cm <- makeContrasts(
  AT2_vs_rest = AT2 - (LAE + FB + MV)/3,
  LAE_vs_rest = LAE - (AT2 + FB + MV)/3,
  FB_vs_rest  = FB  - (AT2 + LAE + MV)/3,
  MV_vs_rest  = MV  - (AT2 + LAE + FB)/3,
  levels = design)
keep <- filterByExpr(y, design=design, min.count=5, min.total.count=10)
y2 <- y[keep,, keep.lib.sizes=FALSE]; y2 <- calcNormFactors(y2, method="TMM")
cat("retained:", nrow(y2), "(claim 1.9: 105,938)\n")
v <- voom(y2, design, plot=FALSE)
cat("running duplicateCorrelation...\n")
cf <- duplicateCorrelation(v, design, block=don)
cat("consensus correlation:", round(cf$consensus.correlation,3), "\n")
fit <- lmFit(v, design, block=don, correlation=cf$consensus.correlation)
fit2 <- eBayes(contrasts.fit(fit, cm))
cat("\n=== 1.10: significant at 5% FDR per CT (claim: LAE 28,930 | AT2 22,131 | MV 21,429 | FB 18,422) ===\n")
for (cn in colnames(cm)) {
  tt <- topTable(fit2, coef=cn, number=Inf, sort.by="none")
  cat(sprintf("  %-12s %6d\n", cn, sum(tt$adj.P.Val < 0.05)))
}
