#!/usr/bin/env Rscript
# =============================================================================
# verify_section3_p2.R — Step-1 factual-accuracy check of Section 3, paragraph 2
# ("custom NMD-anchored normalization + productive DE"), manuscript 2026.7.17.
#
# METHOD: rather than reimplement Yul's normalization (which would test my port,
# not her result), this purls productive_response.Rmd, repoints its /udd/reyle
# paths at local copies of the identical inputs, and executes her code verbatim
# through the `mashr` chunk. Claims are then read off the resulting objects.
#
# Claims under test (§3 ¶2, renumbered map IDs):
#   3.5  Spearman >= 0.99 between adjusted and uncorrected fold changes
#   3.6  non-NMD genes: Spearman >= 0.96, DMSO vs adjusted-Smg1i  (SF23)
#   3.8  median adjusted log2FC: AT2 +0.14, LAE -0.04, FB +0.10, MV +0.05
#   3.9  no significant change in ~83% of genes                    (Fig 2B; ST4)
#   3.10 down 88 (AT2) - 1337 (LAE); up 635 (MV) - 888 (LAE)       (ST4)
# =============================================================================
suppressPackageStartupMessages({ library(knitr); library(dplyr) })

REPO <- path.expand("~/claude_projects/nmd")
RMD  <- file.path(REPO, "code/upstream/productive_response.Rmd")
WORK <- file.path(REPO, "tmp/verify_s3p2");  dir.create(WORK, recursive=TRUE, showWarnings=FALSE)

# ---- 1. purl -----------------------------------------------------------------
src <- file.path(WORK, "productive_response_purled.R")
knitr::purl(RMD, output = src, quiet = TRUE, documentation = 0)
code <- readLines(src)
cat(sprintf("purled %d lines\n", length(code)))

# ---- 2. repoint paths (identical files, local copies) ------------------------
sub1 <- function(pat, rep) { i <- grep(pat, code); stopifnot(length(i) >= 1); code[i] <<- rep; invisible(length(i)) }
sub1('^proj *<-',          sprintf('proj <- "%s"', REPO))
sub1('^OUT_DIR *<-',       sprintf('OUT_DIR <- "%s"', WORK))
sub1('^DGE_ISO_PATH *<-',  sprintf('DGE_ISO_PATH <- "%s"',  file.path(REPO,"nmd_fig_data/dge_isoform_longread_filtered_2026.3.3.rds")))
sub1('^PM_ISO_PATH *<-',   sprintf('PM_ISO_PATH <- "%s"',   file.path(REPO,"isocall_dge/mashr/mashr_isoform_posterior_means_2026.3.10.csv")))
sub1('^LFSR_ISO_PATH *<-', sprintf('LFSR_ISO_PATH <- "%s"', file.path(REPO,"isocall_dge/mashr/mashr_isoform_lfsr_2026.3.10.csv")))
sub1('^SQANTI_PATH *<-',   sprintf('SQANTI_PATH <- "%s"',   file.path(REPO,"nmd_fig_data/nmd_lungcells_classification.txt")))
sub1('^SR_DIR *<-',        sprintf('SR_DIR <- "%s"',        file.path(REPO,"shortread_dge/mashr")))

# ---- 2b. DISCLOSED DEVIATION: cell-type label aliases in the mashr columns ---
# Yul's resolve_mash_col() expects display-name columns (AT2 / LAE). The mashr CSVs
# tracked in THIS repo predate the April-2026 rename and use the internal codes
# (Smg1i_in_AT / Smg1i_in_DD). Per CLAUDE.md the mapping is AT->AT2, DD->LAE. We
# widen the resolver to accept either spelling. This changes column LOOKUP only —
# no data, threshold, or computation is altered.
i <- grep('^resolve_mash_col *<- *function', code)
stopifnot("resolve_mash_col definition not found" = length(i) == 1)
j <- i; while (j < length(code) && !grepl('stop\\("no mash col', code[j])) j <- j + 1
code <- append(code[-(i:j)], values = c(
  'resolve_mash_col <- function(df, ct){',
  '  cn <- colnames(df)',
  '  alias <- c(AT2="AT", LAE="DD", AT="AT2", DD="LAE")',
  '  cands <- unique(c(ct, unname(alias[ct]))); cands <- cands[!is.na(cands)]',
  '  for (k in cands) {',
  '    if (k %in% cn) return(k)',
  '    h <- cn[grepl(paste0("_in_", k, "$"), cn)]; if (length(h)==1) return(h)',
  '    h <- cn[grepl(paste0("^", k, "[._:-]"), cn)]; if (length(h)==1) return(h)',
  '  }',
  '  stop("no mash col ", ct)',
  '}'), after = i - 1)
cat("applied alias patch to resolve_mash_col (AT2<->AT, LAE<->DD)\n")

# same alias problem in the short-read resolver (sr-load chunk)
i <- grep('^resolve_sr *<- *function', code)
stopifnot("resolve_sr definition not found" = length(i) == 1)
j <- i; while (j < length(code) && !grepl('stop\\("no SR col', code[j])) j <- j + 1
code <- append(code[-(i:j)], values = c(
  'resolve_sr <- function(df, ct){',
  '  cn <- colnames(df)',
  '  alias <- c(AT2="AT", LAE="DD", AT="AT2", DD="LAE")',
  '  cands <- unique(c(ct, unname(alias[ct]))); cands <- cands[!is.na(cands)]',
  '  for (k in cands) {',
  '    h <- cn[grepl(paste0("_in_", k, "$"), cn)]; if (length(h)==1) return(h)',
  '    if (k %in% cn) return(k)',
  '  }',
  '  stop("no SR col ", ct)',
  '}'), after = i - 1)
cat("applied alias patch to resolve_sr\n")

# ---- 3. truncate at the END of the `mashr` chunk ----------------------------
# purl(documentation=0) strips chunk headers, so anchor on the mashr chunk's last
# statement: saveRDS(m, .../mashr_custom_<date>.rds). Everything after is the
# exploratory noise battery, not needed for the paragraph-2 claims.
stop_at <- grep('saveRDS\\(m, *file\\.path\\(OUT_DIR', code)
stopifnot("mashr-chunk end anchor not found" = length(stop_at) >= 1)
code <- code[1:stop_at[1]]
patched <- file.path(WORK, "run_through_mashr.R")
writeLines(code, patched)
cat(sprintf("executing Yul's code through the mashr chunk (%d lines)\n\n", length(code)))

# ---- 4. execute --------------------------------------------------------------
env <- new.env()
sink(file.path(WORK, "run_log.txt"), split = FALSE)
ok <- tryCatch({ source(patched, local = env, echo = FALSE); TRUE },
               error = function(e) { sink(); message("ERROR: ", conditionMessage(e)); FALSE })
if (isTRUE(ok)) sink()
if (!ok) { cat("\n--- tail of run log ---\n"); cat(tail(readLines(file.path(WORK,"run_log.txt")), 25), sep="\n"); quit(status=1) }
cat("Yul's pipeline ran clean.\n")

base <- get("base", envir = env)
CT   <- get("CELL_TYPES", envir = env)

cat("\n=========== RECOMPUTED §3 ¶2 ===========\n")

# ---- 3.8 median adjusted productive log2FC (NMD genes) ----------------------
m38 <- sapply(CT, function(c) median(base$custom[base$ct==c & base$gene_category=="NMD"]))
cat("\n-- 3.8 median adjusted log2FC, NMD genes (claim AT2 +0.14, LAE -0.04, FB +0.10, MV +0.05) --\n")
print(round(m38, 3))

# ---- 3.5 adjusted vs uncorrected fold change --------------------------------
cat("\n-- 3.5 Spearman(custom, Kitagawa) = adjusted vs uncorrected logFC (claim >= 0.99) --\n")
c35 <- base %>% group_by(ct) %>%
  summarise(n = n(),
            spearman = round(cor(custom, Kitagawa, method="spearman"), 4),
            pearson  = round(cor(custom, Kitagawa), 4), .groups="drop")
print(as.data.frame(c35))
cat("   NMD genes only:\n")
print(as.data.frame(base %>% filter(gene_category=="NMD") %>% group_by(ct) %>%
  summarise(n=n(), spearman=round(cor(custom, Kitagawa, method="spearman"),4), .groups="drop")))

# ---- 3.6 control (non-NMD) genes: DMSO vs adjusted-Smg1i --------------------
cat("\n-- 3.6 non-NMD genes: Spearman(prodCPM_DMSO, prodCPM_SMG_adj) (claim >= 0.96) --\n")
c36 <- base %>% filter(gene_category=="non-NMD") %>% group_by(ct) %>%
  summarise(n = n(),
            spearman      = round(cor(prodCPM_DMSO, prodCPM_SMG_adj, method="spearman"), 4),
            pearson_log   = round(cor(log2(prodCPM_DMSO+1), log2(prodCPM_SMG_adj+1)), 4),
            median_log2fc = round(median(custom), 4), .groups="drop")
print(as.data.frame(c36))

# ---- 3.9 / 3.10 significant productive change (mashr on the custom fit) -----
pdm <- get("prod_dir_mashr", envir = env)
cat("\n-- 3.9 / 3.10 significant productive change (mashr lfsr) --\n")
cat("dir levels:", paste(sort(unique(as.character(pdm$dir))), collapse=" | "), "\n")
cat("gene_category levels:", paste(sort(unique(as.character(pdm$gene_category))), collapse=" | "), "\n")

tab_all <- pdm %>% count(ct, dir) %>% tidyr::pivot_wider(names_from=dir, values_from=n, values_fill=0)
cat("\n[all gene categories]\n"); print(as.data.frame(tab_all))

nmd_only <- pdm %>% filter(gene_category == "NMD")
tab_nmd <- nmd_only %>% count(ct, dir) %>% tidyr::pivot_wider(names_from=dir, values_from=n, values_fill=0)
cat("\n[NMD genes only — the manuscript's stated population]\n")
tab_nmd <- as.data.frame(tab_nmd)
tot <- rowSums(tab_nmd[, setdiff(names(tab_nmd), "ct"), drop=FALSE])
tab_nmd$total <- tot
nscol <- intersect(c("ns","none","flat","NS"), names(tab_nmd))
if (length(nscol)) tab_nmd$pct_no_change <- round(100*tab_nmd[[nscol[1]]]/tot, 1)
print(tab_nmd)

cat("\n   claim 3.9  : ~83% of genes show no significant change\n")
cat("   claim 3.10 : down 88 (AT2) - 1337 (LAE); up 635 (MV) - 888 (LAE)\n")

saveRDS(list(base=base, m38=m38, c35=c35, c36=c36, prod_dir_mashr=pdm, prod_dir=get("prod_dir",envir=env)),
        file.path(WORK,"p2_objects.rds"))
cat("\nsaved:", file.path(WORK,"p2_objects.rds"), "\n")
