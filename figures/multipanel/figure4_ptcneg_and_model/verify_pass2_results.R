#!/usr/bin/env Rscript
# Pass 2: Result correctness for Figure 4 (NEW 2x2 layout).
# For each panel, recompute the 3 pairwise Wilcoxon comparisons (HL shift, CI,
# Cliff δ, p) from long.tsv and verify they match panelX_pairwise.tsv.
# Then map p → significance marker and verify they would render correctly.

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/Users/petecastaldi/claude_projects/nmd/figures/multipanel/figure4_ptcneg_and_model"
fails <- list()
pass_ct <- 0L
fail_ct <- 0L

check <- function(label, expected, actual, tol = 1e-6) {
  ok <- isTRUE(all.equal(as.numeric(expected), as.numeric(actual), tolerance = tol))
  if (ok) {
    pass_ct <<- pass_ct + 1L
    cat(sprintf("PASS  %-72s exp=%s act=%s\n", label, format(expected, digits = 6),
                format(actual, digits = 6)))
  } else {
    fail_ct <<- fail_ct + 1L
    msg <- sprintf("FAIL  %-72s expected=%s actual=%s", label,
                   format(expected, digits = 8), format(actual, digits = 8))
    fails[[length(fails) + 1L]] <<- msg
    cat(msg, "\n")
  }
  invisible(ok)
}

sig_marker <- function(p) {
  if (is.na(p)) return("n.s.")
  if (p < 1e-4) return("***")
  if (p < 1e-3) return("**")
  if (p < 0.05) return("*")
  return("n.s.")
}

cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) == 0L || length(y) == 0L) return(NA_real_)
  (sum(outer(x, y, ">")) - sum(outer(x, y, "<"))) / (length(x) * length(y))
}

recompute_pairwise <- function(long_path, vcol) {
  d <- fread(long_path)
  d[, (vcol) := as.numeric(get(vcol))]
  groups <- unique(d$group)
  combos <- combn(groups, 2, simplify = FALSE)
  out <- do.call(rbind, lapply(combos, function(p) {
    x <- d[group == p[1], get(vcol)]
    y <- d[group == p[2], get(vcol)]
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    wt <- wilcox.test(x, y, conf.int = TRUE, exact = FALSE)
    data.table(group_x = p[1], group_y = p[2],
               n_x = length(x), n_y = length(y),
               median_x = median(x), median_y = median(y),
               hl_shift = unname(wt$estimate),
               hl_ci_lo = unname(wt$conf.int[1]),
               hl_ci_hi = unname(wt$conf.int[2]),
               cliffs_delta = cliffs_delta(x, y),
               wilcox_p = wt$p.value)
  }))
  out
}

panel_check2 <- function(letter, slug, vcol) {
  cat(sprintf("\n=== PANEL %s — %s ===\n", letter, vcol))
  long_path <- file.path(ROOT, sprintf("data/panel%s_%s_long.tsv", letter, slug))
  pair_path <- file.path(ROOT, sprintf("data/panel%s_%s_pairwise.tsv", letter, slug))
  stored <- fread(pair_path)
  recomp <- recompute_pairwise(long_path, vcol)
  cat("Stored:\n"); print(stored)
  cat("Recomputed:\n"); print(recomp)
  for (i in seq_len(nrow(stored))) {
    s <- stored[i]
    r <- recomp[(group_x == s$group_x & group_y == s$group_y) |
                (group_x == s$group_y & group_y == s$group_x)]
    if (nrow(r) == 0) {
      fail_ct <<- fail_ct + 1L
      fails[[length(fails) + 1L]] <<- sprintf("FAIL Panel %s pair %s vs %s missing in recompute",
                                              letter, s$group_x, s$group_y)
      next
    }
    flipped <- s$group_x != r$group_x
    sign <- if (flipped) -1 else 1
    tag <- sprintf("Panel %s %s vs %s", letter, s$group_x, s$group_y)
    check(paste(tag, "n_x"), s$n_x, if (flipped) r$n_y else r$n_x)
    check(paste(tag, "n_y"), s$n_y, if (flipped) r$n_x else r$n_y)
    check(paste(tag, "hl_shift"), s$hl_shift, sign * r$hl_shift, tol = 0.5)
    if (flipped) {
      check(paste(tag, "hl_ci_lo"), s$hl_ci_lo, -r$hl_ci_hi, tol = 0.5)
      check(paste(tag, "hl_ci_hi"), s$hl_ci_hi, -r$hl_ci_lo, tol = 0.5)
    } else {
      check(paste(tag, "hl_ci_lo"), s$hl_ci_lo, r$hl_ci_lo, tol = 0.5)
      check(paste(tag, "hl_ci_hi"), s$hl_ci_hi, r$hl_ci_hi, tol = 0.5)
    }
    check(paste(tag, "cliffs_delta"), s$cliffs_delta, sign * r$cliffs_delta, tol = 1e-3)
    rel_diff <- abs(s$wilcox_p - r$wilcox_p) / max(s$wilcox_p, r$wilcox_p, 1e-300)
    if (rel_diff < 1e-6) {
      pass_ct <<- pass_ct + 1L
      cat(sprintf("PASS  %s wilcox_p exp=%g act=%g\n", tag, s$wilcox_p, r$wilcox_p))
    } else {
      fail_ct <<- fail_ct + 1L
      msg <- sprintf("FAIL  %s wilcox_p exp=%g act=%g rel_diff=%g", tag, s$wilcox_p,
                     r$wilcox_p, rel_diff)
      fails[[length(fails) + 1L]] <<- msg
      cat(msg, "\n")
    }
    m_stored <- sig_marker(s$wilcox_p)
    m_recomp <- sig_marker(r$wilcox_p)
    if (m_stored == m_recomp) {
      pass_ct <<- pass_ct + 1L
      cat(sprintf("PASS  %s sig marker (%s)\n", tag, m_stored))
    } else {
      fail_ct <<- fail_ct + 1L
      msg <- sprintf("FAIL  %s sig marker stored=%s recomp=%s", tag, m_stored, m_recomp)
      fails[[length(fails) + 1L]] <<- msg
      cat(msg, "\n")
    }
  }
}

panel_check2("A", "5utr_length", "utr5_length_nt")
panel_check2("B", "longest_5utr_orf", "longest_5utr_orf_nt")
panel_check2("C", "5utr_length", "utr5_length_nt")
panel_check2("D", "longest_5utr_orf", "longest_5utr_orf_nt")

cat("\n=== Significance marker summary ===\n")
for (info in list(
  list(letter = "A", slug = "5utr_length"),
  list(letter = "B", slug = "longest_5utr_orf"),
  list(letter = "C", slug = "5utr_length"),
  list(letter = "D", slug = "longest_5utr_orf"))) {
  s <- fread(file.path(ROOT, sprintf("data/panel%s_%s_pairwise.tsv", info$letter, info$slug)))
  cat(sprintf("Panel %s pairwise markers:\n", info$letter))
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %s vs %s : p=%.3g → %s\n", s$group_x[i], s$group_y[i], s$wilcox_p[i],
                sig_marker(s$wilcox_p[i])))
  }
}

cat("\n=================================================\n")
cat(sprintf("PASS 2 SUMMARY: %d passes / %d failures\n", pass_ct, fail_ct))
cat("=================================================\n")
if (fail_ct > 0) {
  cat("\nFailures:\n")
  for (m in fails) cat(" -", m, "\n")
  quit(status = 1)
}
