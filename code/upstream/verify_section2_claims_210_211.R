#!/usr/bin/env Rscript
# =============================================================================
# verify_section2_claims_210_211.R — unblocks §2 claims 2.10 and 2.11.
#
# These were BLOCKED in the 2026-07-21 pass because that attempt joined against
# `all_proportions_2026.1.20.rds` — a stale oarfish-era table covering <17% of
# NMD-susceptible isoforms. Yul has since identified the real producer, and it is
# now local: interpret_isoform_patterns_mashr_2026.3.10.Rmd, which builds
# proportions from the isocall DGEList (the correct universe).
#
# This runs HER Rmd (purled, paths repointed) through the `q3-dmso-vs-smg1i`
# chunk, then reads the claims off her objects.
#
# Claims (manuscript 2026.7.17, §2):
#   2.10 NMD-susceptible isoforms "concentrated at low isoform proportions"
#        (SF5) and "spanned more than four orders of magnitude of absolute
#        expression"
#   2.11 "Between DMSO and SMG1i, these isoforms showed a clear upward shift in
#        both proportion and absolute expression"
#
# DISCLOSED DEVIATIONS (lookup only):
#   - paths repointed to local copies of the same inputs
#   - CELLTYPE_MAP keys at2/lae -> at/dd: the repo's mashr CSVs predate the
#     AT->AT2 / DD->LAE rename. The Rmd deliberately DROPS rows with unmapped
#     codes (see its comment at :206), so without this patch every row would be
#     discarded and the report would come out empty.
# =============================================================================
suppressPackageStartupMessages({ library(knitr); library(dplyr) })

REPO <- path.expand("~/claude_projects/nmd")
RMD  <- file.path(REPO, "code/upstream/interpret_isoform_patterns_mashr_2026.3.10.Rmd")
WORK <- file.path(REPO, "tmp/verify_s2_210"); dir.create(WORK, recursive=TRUE, showWarnings=FALSE)

src <- file.path(WORK, "purled.R")
knitr::purl(RMD, output = src, quiet = TRUE, documentation = 0)
code <- readLines(src)
cat(sprintf("purled %d lines\n", length(code)))

sub1 <- function(pat, rep) { i <- grep(pat, code); stopifnot(length(i) >= 1); code[i] <<- rep; invisible(NULL) }
sub1('^INPUT_DGELIST +<-',          sprintf('INPUT_DGELIST <- "%s"',          file.path(REPO,"nmd_fig_data/dge_isoform_longread_2026.3.3.rds")))
sub1('^INPUT_DGELIST_FILTERED *<-', sprintf('INPUT_DGELIST_FILTERED <- "%s"', file.path(REPO,"nmd_fig_data/dge_isoform_longread_filtered_2026.3.3.rds")))
sub1('^INPUT_DGE_DIR *<-',          sprintf('INPUT_DGE_DIR <- "%s"',          file.path(REPO,"isocall_dge/mashr")))
sub1('^OUTPUT_DIR *<-',             sprintf('OUTPUT_DIR <- "%s"', WORK))
sub1('^OUT_DIR_FIG *<-',            sprintf('OUT_DIR_FIG <- "%s"', WORK))
sub1('^CELLTYPE_MAP *<-',           'CELLTYPE_MAP <- c("at"="AT2","dd"="LAE","fb"="FB","mv"="MV","at2"="AT2","lae"="LAE")')
cat("patched paths + CELLTYPE_MAP (at/dd aliases)\n")

# truncate after the DMSO-vs-Smg1i chunk (before Q6 GSEA)
stop_at <- grep('Q6|fgsea|GSEA_SEED *<- *42', code)
stop_at <- stop_at[stop_at > 200]
if (length(stop_at)) code <- code[1:(min(stop_at)-1)]
patched <- file.path(WORK, "run.R"); writeLines(code, patched)
cat(sprintf("executing %d lines\n\n", length(code)))

env <- new.env(); logf <- file.path(WORK, "log.txt")
sink(logf, split = FALSE)
ok <- tryCatch({ source(patched, local = env, echo = FALSE); TRUE },
               error = function(e) { sink(); message("ERROR: ", conditionMessage(e)); FALSE })
if (isTRUE(ok)) sink()
if (!ok) { cat("\n--- log tail ---\n"); cat(tail(readLines(logf), 25), sep="\n"); quit(status=1) }
cat("Yul's Rmd ran clean.\n")

ap <- get("all_proportions", envir = env)
cat(sprintf("\nall_proportions: %d rows | %d isoforms | CTs: %s | treatments: %s\n",
    nrow(ap), dplyr::n_distinct(ap$transcript_id),
    paste(sort(unique(ap$cell_type)), collapse=","), paste(sort(unique(ap$treatment)), collapse=",")))

# NMD-susceptible set — HER object from the q3-identify chunk.
# NOTE: q3-identify prints "(lfsr < 0.05)" but filters `adj.P.Val < 0.05 & logFC > 0`.
# Per CLAUDE.md those are different quantities (adj.P.Val = limma-voom FDR underlying
# mashr; the mashr lfsr rule is the `nmd_responsive` column). We use HER definition for
# the primary result, and cross-check against the mashr-column rule below.
all_dge <- get("all_dge", envir = env)
nmd <- get("nmd_responsive", envir = env) %>% select(cell_type, transcript_id) %>% distinct()
cat(sprintf("all_dge: %d rows | NMD-susceptible isoform-CT pairs (her rule): %d\n", nrow(all_dge), nrow(nmd)))
if ("nmd_responsive" %in% names(all_dge)) {
  alt <- all_dge %>% filter(nmd_responsive == TRUE) %>% select(cell_type, transcript_id) %>% distinct()
  cat(sprintf("cross-check, mashr-column rule (lfsr<0.05 & pm>0): %d pairs (overlap %d)\n",
              nrow(alt), nrow(dplyr::intersect(nmd, alt))))
}

stopifnot("expected normalized_cpm column" = "normalized_cpm" %in% names(ap))
cat("\n=========== RECOMPUTED §2 claims 2.10 / 2.11 ===========\n")
cat("(NB: the expression column is `normalized_cpm`; a bare `cpm` silently resolves to edgeR::cpm)\n")

d <- ap %>% inner_join(nmd, by = c("cell_type","transcript_id"))
cat(sprintf("\nNMD-susceptible rows in the proportion table: %d (coverage check — the blocked\n", nrow(d)))
cat(sprintf("2026-07-21 attempt covered <17%%; this join covers %.1f%% of NMD isoform-CT pairs)\n",
            100*dplyr::n_distinct(paste(d$cell_type,d$transcript_id))/nrow(nmd)))

# ---- 2.10 baseline (DMSO) proportions --------------------------------------
cat("\n-- 2.10a DMSO isoform proportion of NMD-susceptible isoforms --\n")
cat("   claim: \"concentrated at low isoform proportions\" (SF5)\n\n")
dm <- d %>% filter(treatment == "DMSO") %>%
  group_by(cell_type, transcript_id) %>%
  summarise(prop = mean(isoform_proportion), cpm = mean(normalized_cpm), .groups="drop")
print(as.data.frame(dm %>% group_by(cell_type) %>% summarise(
  n = n(),
  median_prop = round(100*median(prop),1),
  pct_under_50 = round(100*mean(prop < 0.50),1),
  pct_under_10 = round(100*mean(prop < 0.10),1), .groups="drop")), row.names = FALSE)

cat("\n-- 2.10b absolute expression span --\n")
cat("   claim: \"spanned more than four orders of magnitude of absolute expression\"\n\n")
print(as.data.frame(dm %>% filter(cpm > 0) %>% group_by(cell_type) %>% summarise(
  n = n(),
  min_cpm = signif(min(cpm),3), max_cpm = signif(max(cpm),3),
  log10_span = round(log10(max(cpm)/min(cpm)),2),
  span_p1_p99 = round(log10(quantile(cpm,0.99)/quantile(cpm,0.01)),2), .groups="drop")), row.names = FALSE)

# ---- 2.11 DMSO vs Smg1i shift ----------------------------------------------
cat("\n-- 2.11 DMSO vs Smg1i shift in proportion and absolute expression --\n")
cat("   claim: \"clear upward shift in both proportion and absolute expression\"\n\n")
pair <- d %>% group_by(cell_type, transcript_id, treatment) %>%
  summarise(prop = mean(isoform_proportion), cpm = mean(normalized_cpm), .groups="drop") %>%
  tidyr::pivot_wider(names_from = treatment, values_from = c(prop, cpm)) %>%
  filter(!is.na(prop_DMSO), !is.na(prop_Smg1i))
print(as.data.frame(pair %>% group_by(cell_type) %>% summarise(
  n = n(),
  med_prop_DMSO  = round(100*median(prop_DMSO),1),
  med_prop_SMG   = round(100*median(prop_Smg1i),1),
  pct_prop_up    = round(100*mean(prop_Smg1i > prop_DMSO),1),
  p_prop         = signif(wilcox.test(prop_Smg1i, prop_DMSO, paired=TRUE)$p.value, 3),
  med_cpm_DMSO   = round(median(cpm_DMSO),2),
  med_cpm_SMG    = round(median(cpm_Smg1i),2),
  pct_cpm_up     = round(100*mean(cpm_Smg1i > cpm_DMSO),1),
  p_cpm          = signif(wilcox.test(cpm_Smg1i, cpm_DMSO, paired=TRUE)$p.value, 3),
  .groups="drop")), row.names = FALSE)

saveRDS(list(dm=dm, pair=pair), file.path(WORK, "s2_210_211.rds"))
cat("\nsaved\n")
