suppressMessages({library(SummarizedExperiment); library(data.table); library(edgeR)})
setwd("/Users/petecastaldi/claude_projects/nmd")
DATA<-"nmd_fig_data"; P<-function(f) file.path(DATA,f)
CT_KEEP<-c("AT","DD","FB","MV")
extract_ct<-function(s) sapply(strsplit(s,"_"),function(x) paste(x[3:length(x)],collapse="_"))
dge_gene<-readRDS(P("dge_gene_unfiltered_2026.1.2.rds"))
se<-readRDS(P("salmon.merged.gene_counts_length_scaled.FULL.rds")); cm<-assay(se,"counts")
pc<-function(cn){p<-strsplit(cn,"_")[[1]];f<-p[1]
 if(startsWith(f,"AT2")){ct<-"AT";sid<-substring(f,4)}else{ct<-substring(f,1,2);sid<-substring(f,3)}
 paste(sid,p[2],ct,sep="_")}
colnames(cm)<-sapply(colnames(cm),pc)
mct<-gsub("DO_ALI","DO",dge_gene$samples$ct); mct<-gsub("AT2","AT",mct)
cid<-paste(dge_gene$samples$sample_id,dge_gene$samples$treatment,mct,sep="_")
dd<-sapply(colnames(cm),function(nm){if(!grepl("_DD$",nm))return(FALSE); gsub("_DD$","_DD_ALI",nm)%in%cid})
if(sum(dd)>0) cm<-cm[,!dd]
ki<-extract_ct(colnames(cm))%in%CT_KEEP; if(sum(!ki)>0) cm<-cm[,ki]
rownames(cm)<-sub("\\..*","",rownames(cm)); counts_combined<-cm
ci<-fread(P("nmd_lungcells_filtered.count_matrix.txt"),data.table=FALSE)
di<-readRDS(P("dge_isoform_longread_2026.3.3.rds"))
t2g<-unique(di$genes[!is.na(di$genes$txid),c("txid","gene_id")])
rownames(ci)<-ci$id; ci$id<-NULL; mmat<-round(as.matrix(ci))
nms<-gsub("^Sample\\d+_","",colnames(mmat))
colnames(mmat)<-sapply(nms,function(x){p<-strsplit(x,"_")[[1]];tr<-p[length(p)];sid<-p[length(p)-1]
 ct<-paste(p[1:(length(p)-2)],collapse="_"); paste(sid,tr,ct,sep="_")})
kl<-extract_ct(colnames(mmat))%in%CT_KEEP; if(sum(!kl)>0) mmat<-mmat[,kl]
idx<-match(rownames(mmat),t2g$txid); gi<-t2g$gene_id[idx]; hg<-!is.na(gi)
cg<-mmat[hg,]; gi<-sub("\\..*","",gi[hg]); cgl<-rowsum(cg,group=gi)
common_genes<-intersect(rownames(counts_combined),rownames(cgl))
common_samples<-intersect(colnames(counts_combined),colnames(cgl))
cat(sprintf("common genes  : %d\n", length(common_genes)))
cat(sprintf("common samples: %d\n", length(common_samples)))
cs<-counts_combined[common_genes,common_samples]; cl<-cgl[common_genes,common_samples]
kg<-(rowSums(cs)>0)|(rowSums(cl)>0); cs<-cs[kg,]; cl<-cl[kg,]
cat(sprintf("after expressed-filter: %d genes\n", sum(kg)))
ds<-calcNormFactors(DGEList(counts=cs),method="TMM"); dl<-calcNormFactors(DGEList(counts=cl),method="TMM")
xs<-cpm(ds,normalized.lib.sizes=TRUE,log=TRUE); xl<-cpm(dl,normalized.lib.sizes=TRUE,log=TRUE)
pe<-sapply(seq_len(ncol(xs)),function(i) cor(xs[,i],xl[,i],method="pearson"))
sp<-sapply(seq_len(ncol(xs)),function(i) cor(xs[,i],xl[,i],method="spearman"))
cat(sprintf("\nmean Pearson  r = %.3f  (range %.3f - %.3f)\n", mean(pe),min(pe),max(pe)))
cat(sprintf("mean Spearman r = %.3f  (range %.3f - %.3f)\n", mean(sp),min(sp),max(sp)))
