#!/usr/bin/env Rscript
# =============================================================================
# _run_productive_response.R — shared harness for the §3 verification scripts.
#
# Executes Yul's productive_response.Rmd verbatim (purled) with its /udd/reyle
# paths repointed at local copies of the identical inputs, stopping at the end of
# the `mashr` chunk. Returns the resulting environment so each paragraph's
# verification script can read the real objects (base, prod_dir_mashr, sr_long,
# kit_nmd_table(), ...) rather than reimplementing the analysis.
#
# DISCLOSED DEVIATIONS (lookup only — no data, threshold, or computation altered):
#   The mashr CSVs tracked in this repo predate the April-2026 cell-type rename,
#   so they carry Smg1i_in_AT / Smg1i_in_DD where Yul's code expects AT2 / LAE.
#   resolve_mash_col() and resolve_sr() are widened to accept either spelling.
#
# Usage:  src <- source("code/upstream/_run_productive_response.R")$value
#         env <- run_productive_response()
# =============================================================================

run_productive_response <- function(repo = path.expand("~/claude_projects/nmd"),
                                    cache = TRUE, quiet = TRUE) {
  suppressPackageStartupMessages({ library(knitr) })
  RMD  <- file.path(repo, "code/upstream/productive_response.Rmd")
  WORK <- file.path(repo, "tmp/verify_s3p2")
  dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
  cache_f <- file.path(WORK, "env_cache.rds")

  if (cache && file.exists(cache_f) &&
      file.mtime(cache_f) > file.mtime(RMD)) {
    if (!quiet) cat("[runner] loading cached env\n")
    return(readRDS(cache_f))
  }

  src <- file.path(WORK, "productive_response_purled.R")
  knitr::purl(RMD, output = src, quiet = TRUE, documentation = 0)
  code <- readLines(src)

  sub1 <- function(pat, rep) {
    i <- grep(pat, code); stopifnot(length(i) >= 1); code[i] <<- rep; invisible(NULL)
  }
  sub1('^proj *<-',          sprintf('proj <- "%s"', repo))
  sub1('^OUT_DIR *<-',       sprintf('OUT_DIR <- "%s"', WORK))
  sub1('^DGE_ISO_PATH *<-',  sprintf('DGE_ISO_PATH <- "%s"',  file.path(repo,"nmd_fig_data/dge_isoform_longread_filtered_2026.3.3.rds")))
  sub1('^PM_ISO_PATH *<-',   sprintf('PM_ISO_PATH <- "%s"',   file.path(repo,"isocall_dge/mashr/mashr_isoform_posterior_means_2026.3.10.csv")))
  sub1('^LFSR_ISO_PATH *<-', sprintf('LFSR_ISO_PATH <- "%s"', file.path(repo,"isocall_dge/mashr/mashr_isoform_lfsr_2026.3.10.csv")))
  sub1('^SQANTI_PATH *<-',   sprintf('SQANTI_PATH <- "%s"',   file.path(repo,"nmd_fig_data/nmd_lungcells_classification.txt")))
  sub1('^SR_DIR *<-',        sprintf('SR_DIR <- "%s"',        file.path(repo,"shortread_dge/mashr")))

  patch_resolver <- function(code, fname, stop_pat, body) {
    i <- grep(paste0('^', fname, ' *<- *function'), code)
    stopifnot(length(i) == 1)
    j <- i; while (j < length(code) && !grepl(stop_pat, code[j])) j <- j + 1
    append(code[-(i:j)], values = body, after = i - 1)
  }
  alias_lines <- function(fname, loop_body) c(
    sprintf('%s <- function(df, ct){', fname),
    '  cn <- colnames(df)',
    '  alias <- c(AT2="AT", LAE="DD", AT="AT2", DD="LAE")',
    '  cands <- unique(c(ct, unname(alias[ct]))); cands <- cands[!is.na(cands)]',
    '  for (k in cands) {', loop_body, '  }',
    sprintf('  stop("no col ", ct)'), '}')

  code <- patch_resolver(code, 'resolve_mash_col', 'stop\\("no mash col',
    alias_lines('resolve_mash_col', paste0(
      '    if (k %in% cn) return(k)\n',
      '    h <- cn[grepl(paste0("_in_", k, "$"), cn)]; if (length(h)==1) return(h)\n',
      '    h <- cn[grepl(paste0("^", k, "[._:-]"), cn)]; if (length(h)==1) return(h)')))
  code <- patch_resolver(code, 'resolve_sr', 'stop\\("no SR col',
    alias_lines('resolve_sr', paste0(
      '    h <- cn[grepl(paste0("_in_", k, "$"), cn)]; if (length(h)==1) return(h)\n',
      '    if (k %in% cn) return(k)')))

  stop_at <- grep('saveRDS\\(m, *file\\.path\\(OUT_DIR', code)
  stopifnot("mashr-chunk end anchor not found" = length(stop_at) >= 1)
  code <- code[1:stop_at[1]]
  patched <- file.path(WORK, "run_through_mashr.R")
  writeLines(code, patched)

  env <- new.env()
  logf <- file.path(WORK, "run_log.txt")
  sink(logf, split = FALSE)
  ok <- tryCatch({ source(patched, local = env, echo = FALSE); TRUE },
                 error = function(e) { sink(); message("RUNNER ERROR: ", conditionMessage(e)); FALSE })
  if (isTRUE(ok)) sink()
  if (!ok) { cat(tail(readLines(logf), 20), sep = "\n"); stop("pipeline failed") }
  if (!quiet) cat("[runner] pipeline ran clean\n")

  if (cache) saveRDS(env, cache_f)
  env
}

run_productive_response
