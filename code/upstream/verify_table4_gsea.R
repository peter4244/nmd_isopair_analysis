suppressMessages({library(data.table); library(dplyr); library(fgsea); library(msigdbr)})
setwd("/Users/petecastaldi/claude_projects/nmd")
GMIN<-15; GMAX<-500
mk <- function(coll,sub=NULL) if(is.null(sub)) msigdbr(species="Homo sapiens",collection=coll) else msigdbr(species="Homo sapiens",collection=coll,subcollection=sub)
allsets <- bind_rows(mutate(mk("H"),collection="Hallmark"), mutate(mk("C5","GO:BP"),collection="GO_BP"),
                     mutate(mk("C2","CP:KEGG_LEGACY"),collection="KEGG"), mutate(mk("C2","CP:REACTOME"),collection="Reactome"))
pw <- split(allsets$gene_symbol, allsets$gs_name)
pf <- pw[sapply(pw,length)>=GMIN & sapply(pw,length)<=GMAX]
# AT2 isoform-level ranking: dedupe to largest |logFC| per gene symbol, rank by signed logFC
g <- fread("isocall_dge/mashr/nmd_mashr_die_at_2026.3.10.csv")
setnames(g, "hgnc_symbol", "gene_symbol")
rk <- g %>% filter(!is.na(gene_symbol), gene_symbol!="", !is.na(logFC), !grepl("^ENSG", gene_symbol)) %>%
      arrange(desc(abs(logFC))) %>% distinct(gene_symbol, .keep_all=TRUE)
cat("genes ranked:", nrow(rk), "\n")
rv <- setNames(rk$logFC, rk$gene_symbol)
set.seed(42)
res <- fgsea(pathways=pf, stats=rv, minSize=GMIN, maxSize=GMAX)
top <- as.data.table(res)[NES>0 & padj<0.05][order(pval)][1:6, .(pathway, NES=round(NES,2), padj=signif(padj,2), size)]
cat("\n=== AT2 isoform-level top (compare Table 4: RNA splicing 1.33/9.0e-17/438; via transesterif 1.34/7.9e-14/316; capped intron 1.33/1.6e-11/291; ribosome biogenesis 1.28/1.7e-9/325; translation 1.26/3.3e-9/364) ===\n")
print(top)
