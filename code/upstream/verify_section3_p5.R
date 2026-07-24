#!/usr/bin/env Rscript
# =============================================================================
# verify_section3_p5.R — Step-1 check of Section 3, paragraph 5
# (the three worked examples), manuscript 2026.7.17.
#
# Per-gene lookups from Yul's `base` + `sr_long` + `prod_dir_mashr` via the shared
# harness. All four cell types are printed for each gene so the numbers themselves
# identify which CT the manuscript is describing — this also resolves F-4, where
# the SRSF2 prose omits its cell type.
#
# Claims under test (map IDs):
#   3.21 SHMT2 (MV):  NMD susceptible; productive 86.5 -> 111.1 CPM;
#                     adjusted logFC +0.42; short-read logFC +0.56        (Fig 2C)
#   3.22 SRSF2 (CT unstated in prose; Fig 2D legend says MV):
#                     productive 84 -> 50 CPM; adjusted logFC -0.73;
#                     productive fraction 91% -> 26%                      (Fig 2D)
#   3.23 PCNA (LAE):  productive 251 -> 130 CPM; adjusted logFC -0.95;
#                     productive fraction 100% -> 93%; SR logFC -0.81     (Fig 2E)
# =============================================================================
suppressPackageStartupMessages({ library(dplyr) })
REPO <- path.expand("~/claude_projects/nmd")
run  <- source(file.path(REPO, "code/upstream/_run_productive_response.R"))$value
env  <- run(REPO, cache = TRUE, quiet = FALSE)

base           <- get("base",           envir = env)
sr_long        <- get("sr_long",        envir = env)
prod_dir_mashr <- get("prod_dir_mashr", envir = env)

look <- function(sym) {
  b <- base %>% filter(hgnc_symbol == sym) %>%
    select(ct, gene_id, gene_category, prodCPM_DMSO, prodCPM_SMG_adj, prodCPM_SMG_raw,
           custom, prodfrac_DMSO, prodfrac_SMG) %>%
    left_join(sr_long %>% select(ct, gene_id, sr_logFC), by = c("ct","gene_id")) %>%
    left_join(prod_dir_mashr %>% select(ct, gene_id, dir, pm, lfsr), by = c("ct","gene_id"))
  if (!nrow(b)) { cat("  (no rows for ", sym, ")\n", sep=""); return(invisible(NULL)) }
  b %>% mutate(
      prodCPM_DMSO    = round(prodCPM_DMSO, 1),
      prodCPM_SMG_adj = round(prodCPM_SMG_adj, 1),
      prodCPM_SMG_raw = round(prodCPM_SMG_raw, 1),
      custom          = round(custom, 3),
      prodfrac_DMSO   = round(100*prodfrac_DMSO, 1),
      prodfrac_SMG    = round(100*prodfrac_SMG, 1),
      sr_logFC        = round(sr_logFC, 3),
      pm              = round(pm, 3)) %>%
    select(-gene_id) %>% as.data.frame()
}

cat("\n=========== RECOMPUTED §3 ¶5 (worked examples) ===========\n")
cat("prodfrac columns are PERCENT. 'custom' = adjusted productive log2FC.\n")

cat("\n===== 3.21  SHMT2 =====\n")
cat("claim (MV): productive 86.5 -> 111.1 CPM | adjusted logFC +0.42 | SR logFC +0.56 | NMD susceptible\n\n")
print(look("SHMT2"), row.names = FALSE)

cat("\n===== 3.22  SRSF2 =====\n")
cat("claim (CT NOT stated in prose; Fig 2D legend = MV):\n")
cat("  productive 84 -> 50 CPM | adjusted logFC -0.73 | productive fraction 91% -> 26%\n\n")
print(look("SRSF2"), row.names = FALSE)

cat("\n===== 3.23  PCNA =====\n")
cat("claim (LAE): productive 251 -> 130 CPM | adjusted logFC -0.95 | fraction 100% -> 93% | SR logFC -0.81\n\n")
print(look("PCNA"), row.names = FALSE)

# ---- which CT best matches each claim? --------------------------------------
cat("\n=========== BEST-MATCHING CELL TYPE (by stated CPM pair) ===========\n")
match_ct <- function(sym, dmso, smg) {
  b <- base %>% filter(hgnc_symbol == sym)
  if (!nrow(b)) return(invisible(NULL))
  b <- b %>% mutate(err = abs(prodCPM_DMSO - dmso) + abs(prodCPM_SMG_adj - smg)) %>% arrange(err)
  cat(sprintf("%-6s claim %.1f -> %.1f  ==>  best = %s (%.1f -> %.1f, |err| %.2f)\n",
              sym, dmso, smg, b$ct[1], b$prodCPM_DMSO[1], b$prodCPM_SMG_adj[1], b$err[1]))
}
match_ct("SHMT2", 86.5, 111.1)
match_ct("SRSF2", 84, 50)
match_ct("PCNA", 251, 130)
cat("\n")
