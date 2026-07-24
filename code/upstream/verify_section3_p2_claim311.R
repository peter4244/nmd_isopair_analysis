#!/usr/bin/env Rscript
# =============================================================================
# verify_section3_p2_claim311.R — closes §3 ¶2 claim 3.11 (previously blocked on
# a missing topGO install).
#
# Claim (§3 ¶2): "Genes where productive output rose were enriched in the
# unfolded-protein response pathway (leading-edge genes including the ER
# chaperone HSPA5 (BiP) and the UPR transducers XBP1 and EIF2AK3) and for
# NF-kB/TNFa inflammatory signaling (ST4)."
#
# Two sub-claims, each replicated with Yul's own code:
#   (a) UPR  — chunk `leading-edge-check` (GO:0006986, response to unfolded
#              protein) over productive-UP NMD genes; plus a Fisher enrichment
#              test, since the prose asserts "enriched", which that chunk (a
#              membership listing) does not itself test.
#   (b) NFkB — chunk `step16-2-up-NFkB`, gsea_one() on
#              HALLMARK_TNFA_SIGNALING_VIA_NFKB via msigdbr.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(org.Hs.eg.db); library(AnnotationDbi); library(msigdbr)
})
REPO <- path.expand("~/claude_projects/nmd")
run  <- source(file.path(REPO, "code/upstream/_run_productive_response.R"))$value
env  <- run(REPO, cache = TRUE, quiet = FALSE)

prod_dir_mashr <- get("prod_dir_mashr", envir = env)
CELL_TYPES     <- get("CELL_TYPES",     envir = env)
svf            <- get("svf",            envir = env)

cat("\n=========== RECOMPUTED §3 ¶2 claim 3.11 ===========\n")

# ---- (a) UPR leading-edge genes (Yul's leading-edge-check, verbatim) --------
upr_ens <- unique(na.omit(
  AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0006986",
                        keytype = "GOALL", columns = "ENSEMBL")$ENSEMBL))
cat(sprintf("\nGO:0006986 (response to unfolded protein) ENSEMBL genes: %d\n", length(upr_ens)))

TARGETS <- c("HSPA5","XBP1","EIF2AK3")
cat("\n-- (a) UPR genes among productive-UP NMD genes, per CT --\n")
for (ct_i in CELL_TYPES) {
  up_g <- prod_dir_mashr %>% filter(ct==ct_i, gene_category=="NMD", dir=="up") %>% pull(gene_id)
  le   <- intersect(svf(up_g), upr_ens)
  syms <- sort(na.omit(AnnotationDbi::mapIds(org.Hs.eg.db, keys = le, column = "SYMBOL",
                                             keytype = "ENSEMBL", multiVals = "first")))
  hit <- TARGETS[TARGETS %in% syms]
  cat(sprintf("%-4s UPR-UP n=%-3d | named-in-paper present: %s\n", ct_i, length(le),
              ifelse(length(hit), paste(hit, collapse=", "), "NONE")))
  cat("      ", paste(syms, collapse=", "), "\n")
}

# ---- (a2) is UPR actually ENRICHED among the up genes? ---------------------
cat("\n-- (a2) Fisher test: UPR membership, productive-UP vs rest (NMD genes) --\n")
cat("   (the prose says \"enriched\"; leading-edge-check only lists membership)\n\n")
fish <- bind_rows(lapply(CELL_TYPES, function(ct_i) {
  d <- prod_dir_mashr %>% filter(ct==ct_i, gene_category=="NMD") %>%
    mutate(ens = svf(gene_id), is_upr = ens %in% upr_ens, is_up = dir == "up")
  tt <- table(factor(d$is_up, c(FALSE,TRUE)), factor(d$is_upr, c(FALSE,TRUE)))
  ft <- fisher.test(tt)
  data.frame(ct = ct_i, n_up = sum(d$is_up), upr_in_up = sum(d$is_up & d$is_upr),
             upr_total = sum(d$is_upr),
             OR = round(as.numeric(ft$estimate), 2), p = signif(ft$p.value, 3))
}))
print(as.data.frame(fish), row.names = FALSE)

# ---- (b) NF-kB / TNFa (Yul's gsea_one, verbatim) --------------------------
cat("\n-- (b) HALLMARK_TNFA_SIGNALING_VIA_NFKB, continuous test on productive pm --\n")
hm <- tryCatch(msigdbr(species="Homo sapiens", category="H"),
               error=function(e) msigdbr(species="Homo sapiens", collection="H"))
nfkb_set <- unique(hm$gene_symbol[hm$gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"])
cat(sprintf("   NF-kB/TNFa set size: %d\n\n", length(nfkb_set)))
fmt_p <- function(p) ifelse(is.na(p), NA, ifelse(p>=1e-3, formatC(p,format="f",digits=3), formatC(p,format="e",digits=1)))
gsea_one <- function(ct_i, class, set){
  d <- prod_dir_mashr %>% filter(ct==ct_i, gene_category==class, !is.na(pm))
  d <- d[!duplicated(d$hgnc_symbol),]
  inset <- d$hgnc_symbol %in% set; w <- wilcox.test(d$pm[inset], d$pm[!inset])
  data.frame(ct=ct_i, class=class, n_in=sum(inset),
             median_pm_in=round(median(d$pm[inset]),3),
             median_pm_out=round(median(d$pm[!inset]),3),
             dir=ifelse(median(d$pm[inset])>median(d$pm[!inset]),"UP","down"),
             p=fmt_p(w$p.value))
}
print(as.data.frame(bind_rows(lapply(CELL_TYPES, function(c)
  bind_rows(gsea_one(c,"non-NMD",nfkb_set), gsea_one(c,"NMD",nfkb_set))))), row.names = FALSE)
cat("\n")
