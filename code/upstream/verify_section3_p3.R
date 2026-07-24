#!/usr/bin/env Rscript
# =============================================================================
# verify_section3_p3.R — Step-1/2 check of Section 3, paragraph 3
# (Kitagawa decomposition), manuscript 2026.7.17.
#
# Runs Yul's productive_response.Rmd via the shared harness, then calls her own
# kit_nmd_table() / kitagawa_decomp() to reproduce the paragraph's percentages.
#
# Claims under test (map IDs):
#   3.13 across significant NMD genes: level dominant 86.2%, composition 13.5%
#        (note 86.2 + 13.5 = 99.7 -> a 0.3% residual to account for)
#   3.14 productive-UP genes: 99.2% level, 0.3% composition (sums to 99.5)
#   3.15 productive-DOWN genes: 32.1% composition-dominated, 67.9% transcription
#   3.16 conclusion holds under the custom normalization w/ SR induction
#
# Population per the code: prod_dir_mashr dir in {up,down} & gene_category=="NMD",
# FB EXCLUDED (provisional) -> AT2 / LAE / MV.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr) })
REPO <- path.expand("~/claude_projects/nmd")
run  <- source(file.path(REPO, "code/upstream/_run_productive_response.R"))$value
env  <- run(REPO, cache = TRUE, quiet = FALSE)

kit_nmd_table <- get("kit_nmd_table", envir = env)
base <- get("base", envir = env); prod_dir_mashr <- get("prod_dir_mashr", envir = env)

cat("\n=========== RECOMPUTED §3 ¶3 (Kitagawa) ===========\n")

# ---- the manuscript population -------------------------------------------
kt <- kit_nmd_table()                    # FB excluded, significant NMD genes
cat(sprintf("\npopulation: significant NMD genes, FB excluded — n = %d\n", nrow(kt)))
print(kt %>% count(ct), row.names = FALSE)

# ---- decomposition identity (does level + composition == dP?) -------------
cat(sprintf("\nidentity check: max |dP - (level+composition)| = %.3e\n",
            max(abs(kt$dP - (kt$level + kt$composition)), na.rm = TRUE)))
cat(sprintf("NA driver rows: %d   NA level: %d   NA composition: %d\n",
            sum(is.na(kt$driver)), sum(is.na(kt$level)), sum(is.na(kt$composition))))

# ---- 3.13 overall driver split (Yul's S5 chunk verbatim) ------------------
cat("\n-- 3.13 overall (claim: level 86.2%, composition 13.5%) --\n")
s5 <- kt %>% count(driver) %>% mutate(pct = round(100 * n / sum(n), 1))
print(as.data.frame(s5), row.names = FALSE)

# ---- 3.14 / 3.15 split by direction --------------------------------------
cat("\n-- 3.14 productive-UP (claim: 99.2% level, 0.3% composition) --\n")
up <- kt %>% filter(dir == "up") %>% count(driver) %>% mutate(pct = round(100*n/sum(n), 1))
print(as.data.frame(up), row.names = FALSE)

cat("\n-- 3.15 productive-DOWN (claim: 32.1% composition, 67.9% transcription) --\n")
dn <- kt %>% filter(dir == "down") %>% count(driver) %>% mutate(pct = round(100*n/sum(n), 1))
print(as.data.frame(dn), row.names = FALSE)

# ---- does the residual come from ties? ------------------------------------
cat("\n-- residual probe: exact |composition| == |level| ties --\n")
cat(sprintf("ties: %d of %d (%.2f%%)\n", sum(abs(kt$composition) == abs(kt$level), na.rm=TRUE),
            nrow(kt), 100*mean(abs(kt$composition) == abs(kt$level), na.rm=TRUE)))

# ---- sensitivity: including FB, and using all NMD genes -------------------
cat("\n-- sensitivity --\n")
kt_fb <- kit_nmd_table(exclude_fb = FALSE)
cat("including FB:\n")
print(as.data.frame(kt_fb %>% count(driver) %>% mutate(pct=round(100*n/sum(n),1))), row.names=FALSE)

cat("\nper-CT driver split (FB excluded):\n")
print(as.data.frame(kt %>% count(ct, driver) %>% group_by(ct) %>%
  mutate(pct = round(100*n/sum(n),1)) %>% ungroup()), row.names = FALSE)

saveRDS(kt, file.path(REPO, "tmp/verify_s3p2/p3_kitagawa.rds"))
cat("\nsaved kitagawa table\n")
