# Independent re-derivation of §1-§2 manuscript claims from the data bundle.
# Run from repo root. Requires nmd_fig_data/ (gitignored, local only).
suppressMessages({library(data.table)})
cts <- c(at="AT2", dd="LAE", fb="FB", mv="MV")

cat("== 1.2/1.4/1.5: expressed isoform set ==\n")
f  <- readRDS("nmd_fig_data/dge_isoform_longread_filtered_2026.3.3.rds")
ids <- rownames(f); if (is.null(ids)) ids <- f$genes$txid
cl <- fread("sqanti/nmd_lungcells/results/nmd_lungcells_classification.txt",
            select=c("isoform","associated_gene","structural_category"))
sub <- cl[isoform %in% ids]
cat("isoforms:", nrow(f), " genes(SQANTI associated_gene):", uniqueN(sub$associated_gene), "\n")
print(sub[, .N, by=structural_category][order(-N)])
n <- sub[, .N, by=associated_gene]$N
cat(sprintf("isoforms/gene: median %.0f mean %.1f IQR %.0f-%.0f max %d | >=2 %d (%.1f%%) | >=10 %d (%.1f%%)\n",
  median(n), mean(n), quantile(n,.25), quantile(n,.75), max(n),
  sum(n>=2), 100*mean(n>=2), sum(n>=10), 100*mean(n>=10)))

cat("\n== 2.1/2.2 SR gene, 2.3/2.4 LR isoform ==\n")
grab <- function(fmt, dir, idcands) {
  out <- list()
  for (k in names(cts)) {
    d <- fread(sprintf(fmt, dir, k))
    idc <- intersect(idcands, names(d))[1]
    cat(sprintf("%-4s tested %6d sig %6d nmd_resp %6d (%.2f%% of sig)\n", cts[k], nrow(d),
        sum(d$adj.P.Val<0.05), sum(d$nmd_responsive %in% TRUE),
        100*sum(d$nmd_responsive %in% TRUE)/sum(d$adj.P.Val<0.05)))
    out[[cts[k]]] <- unique(d[[idc]][d$nmd_responsive %in% TRUE])
  }
  out
}
sr <- grab("%s/nmd_mashr_dge_%s_2026.3.10.csv", "shortread_dge/mashr",
           c("ensembl_gene_id_version","gene_id","V1"))
lr <- grab("%s/nmd_mashr_die_%s_2026.3.10.csv", "isocall_dge/mashr", c("txid","V1"))

cat("\n== 2.5/2.6 union + sharing (LR isoform) ==\n")
tb <- table(unlist(lr))
cat(sprintf("union %d | all-4 %d (%.1f%%) | exactly-1 %d (%.1f%%)\n",
  length(tb), sum(tb==4), 100*sum(tb==4)/length(tb), sum(tb==1), 100*sum(tb==1)/length(tb)))
cat("partition by #CTs:", paste(sprintf("%dCT=%d", 1:4, tabulate(tb,4)), collapse=" "), "\n")

cat("\n== 2.12/2.13 CT-specificity ==\n")
spec <- function(l, label) { t <- table(unlist(l)); s <- names(t)[t==1]
  for (k in names(l)) cat(sprintf("  %s %-4s %5d (%.1f%%)\n", label, k,
    sum(s %in% l[[k]]), 100*sum(s %in% l[[k]])/length(l[[k]]))) }
spec(sr, "SR"); spec(lr, "LR")
