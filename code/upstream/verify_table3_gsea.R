suppressMessages({library(data.table); library(dplyr); library(fgsea); library(msigdbr); library(stringr)})
setwd("/Users/petecastaldi/claude_projects/nmd")
cat("msigdbr:", as.character(packageVersion("msigdbr")), "\n")
GSEA_MIN <- 15; GSEA_MAX <- 500
mk <- function(coll, sub=NULL) { a <- if(is.null(sub)) msigdbr(species="Homo sapiens", collection=coll) else msigdbr(species="Homo sapiens", collection=coll, subcollection=sub); a }
h  <- mk("H"); bp <- mk("C5","GO:BP"); kg <- mk("C2","CP:KEGG_LEGACY"); re <- mk("C2","CP:REACTOME")
allsets <- bind_rows(mutate(h,collection="Hallmark"), mutate(bp,collection="GO_BP"),
                     mutate(kg,collection="KEGG"), mutate(re,collection="Reactome"))
pathways <- split(allsets$gene_symbol, allsets$gs_name)
pf <- pathways[sapply(pathways,length)>=GSEA_MIN & sapply(pathways,length)<=GSEA_MAX]
cat("pathways after size filter:", length(pf), "\n")
# AT2 ranking
g <- fread("shortread_dge/mashr/nmd_mashr_dge_at_2026.3.10.csv")
g <- g[!is.na(hgnc_symbol) & hgnc_symbol != ""]
rk <- g[, .SD[which.max(abs(logFC))], by=hgnc_symbol][order(-logFC)]
rv <- setNames(rk$logFC, rk$hgnc_symbol)
set.seed(42)
res <- fgsea(pathways=pf, stats=rv, minSize=GSEA_MIN, maxSize=GSEA_MAX, nPermSimple=100000)
top <- res[NES>0 & padj<0.05][order(pval)][1:6, .(pathway, NES=round(NES,2), padj=signif(padj,2), size)]
cat("\n=== AT2 top gene-level pathways (compare to updated Table 3) ===\n")
print(top)
cat("\n-- Table 3 (updated) AT2 top rows: RNA splicing 1.48/0.007/462; RNA splicing via transesterif; capped intron; ribosome biogenesis; translation --\n")
