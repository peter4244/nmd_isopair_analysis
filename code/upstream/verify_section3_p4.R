#!/usr/bin/env Rscript
# =============================================================================
# verify_section3_p4.R — Step-1/2 check of Section 3, paragraph 4
# (two mechanisms partitioning the productive-down genes), manuscript 2026.7.17.
#
# Replicates Yul's S6b chunk `down-arm-mechanism` using HER functions
# (classify_mech, theme_enrichment, THEME_SETS, test_set, sym_of) taken from the
# cached environment built by _run_productive_response.R.
#
# Claims under test (map IDs):
#   3.16 conclusion holds under the custom normalization + SR induction signal
#   3.17 composition-shift genes enriched for RNA splicing        (OR 2.5, p<1e-3)
#   3.18 ... and SR/hnRNP factors                                 (OR 6.8, p<1e-3)
#   3.19 of the 12 canonical SR proteins, EVERY significant productive change was
#        a decrease (8 LAE, 3 AT2, 2 MV)
#   3.20 transcriptional-repression genes enriched for DNA replication/repair
#                                                                 (OR 2.1, p<1e-3)
# Population per the code: dir=="down", FB excluded (provisional).
# =============================================================================
suppressPackageStartupMessages({ library(dplyr) })
REPO <- path.expand("~/claude_projects/nmd")
run  <- source(file.path(REPO, "code/upstream/_run_productive_response.R"))$value
env  <- run(REPO, cache = TRUE, quiet = FALSE)

prod_dir_mashr <- get("prod_dir_mashr", envir = env)
base           <- get("base",           envir = env)
sr_long        <- get("sr_long",        envir = env)
classify_mech  <- get("classify_mech",  envir = env)
theme_enrichment <- get("theme_enrichment", envir = env)
THEME_SETS     <- get("THEME_SETS",     envir = env)
CELL_TYPES     <- get("CELL_TYPES",     envir = env)

cat("\n=========== RECOMPUTED §3 ¶4 ===========\n")

# ---- S6b: down-arm mechanism split (Yul's code, verbatim) ------------------
pd_mech <- prod_dir_mashr %>% filter(dir == "down", ct != "FB") %>%
  left_join(base %>% select(ct, gene_id, prodfrac_DMSO, prodfrac_SMG), by = c("ct","gene_id")) %>%
  left_join(sr_long %>% select(ct, gene_id, sr_logFC), by = c("ct","gene_id")) %>%
  classify_mech("down")

cat("\n-- 3.16 down-arm mechanism split (custom normalization + SR induction) --\n")
split_tab <- pd_mech %>% count(mechanism) %>% mutate(pct = round(100*n/sum(n), 1))
print(as.data.frame(split_tab), row.names = FALSE)
cat("\n   per cell type:\n")
print(as.data.frame(pd_mech %>% count(ct, mechanism) %>% group_by(ct) %>%
        mutate(pct = round(100*n/sum(n),1)) %>% ungroup()), row.names = FALSE)

# ---- 3.17 / 3.18 / 3.20 theme enrichment ----------------------------------
trans_g    <- pd_mech %>% filter(mechanism == "transcriptional") %>% pull(gene_id) %>% unique()
comp_g     <- pd_mech %>% filter(mechanism == "composition")     %>% pull(gene_id) %>% unique()
universe_g <- unique(prod_dir_mashr$gene_id[prod_dir_mashr$ct != "FB"])
cat(sprintf("\ntranscriptional=%d | composition=%d | universe=%d\n",
            length(trans_g), length(comp_g), length(universe_g)))

cat("\n-- 3.17 / 3.18 / 3.20 theme enrichment --\n")
cat("   claim 3.17 composition x RNA splicing      : OR 2.5, p < 1e-3\n")
cat("   claim 3.18 composition x SR/hnRNP factors  : OR 6.8, p < 1e-3\n")
cat("   claim 3.20 transcriptional x DNA repl/repair: OR 2.1, p < 1e-3\n\n")
te <- theme_enrichment(trans_g, comp_g, universe_g)
print(as.data.frame(te), row.names = FALSE)

# ---- 3.19 canonical SR proteins ------------------------------------------
cat("\n-- 3.19 canonical SR proteins SRSF1-SRSF12 --\n")
cat("   claim: EVERY significant productive change was a decrease (8 LAE, 3 AT2, 2 MV)\n\n")
SRSF <- paste0("SRSF", 1:12)
sr_tab <- prod_dir_mashr %>% filter(hgnc_symbol %in% SRSF)
cat("all SRSF rows by ct x dir (all cell types incl. FB):\n")
print(as.data.frame(sr_tab %>% count(ct, dir) %>%
        tidyr::pivot_wider(names_from = dir, values_from = n, values_fill = 0)), row.names = FALSE)
cat("\nSIGNIFICANT SRSF changes only (dir != ns), per CT:\n")
sig_sr <- sr_tab %>% filter(dir != "ns")
print(as.data.frame(sig_sr %>% count(ct, dir)), row.names = FALSE)
cat("\ngene-level detail of every significant SRSF change:\n")
print(as.data.frame(sig_sr %>% select(ct, hgnc_symbol, dir, pm, lfsr, gene_category) %>%
        arrange(ct, hgnc_symbol)), row.names = FALSE)
cat(sprintf("\nany significant SRSF INCREASE? %s\n",
            ifelse(any(sig_sr$dir == "up"), "YES — contradicts the claim", "no — consistent with the claim")))
cat(sprintf("distinct SRSF genes detected: %d of 12 (%s)\n",
            length(unique(sr_tab$hgnc_symbol)), paste(sort(unique(sr_tab$hgnc_symbol)), collapse=", ")))

saveRDS(list(pd_mech = pd_mech, te = te, sr_tab = sr_tab),
        file.path(REPO, "tmp/verify_s3p2/p4_objects.rds"))
cat("\nsaved\n")
