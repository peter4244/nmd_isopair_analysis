#!/usr/bin/env Rscript
# =============================================================================
# build_4ct_pheno.R — sample metadata for the 26 manuscript samples.
#
# INTERNAL packaging step. Not shipped; see DEPOSIT_PROVENANCE_INTERNAL.md.
#
# WHY (2026-07-25)
# ---------------------------------------------------------------------------
# NMD_shortread_dge_fullmodel_2026.5.5.Rmd reads two phenotype sheets from a
# Dropbox path (.P$PHENO). They drive the design formula (~ ct + treatment +
# ct:treatment, blocked on donor) and the 38 -> 26 sample selection, and they are in
# neither the Zenodo deposit nor GEO. Without them §1 cannot be reproduced at all.
#
# The sheets do NOT distinguish submerged LAE from LAE-ALI. The published code infers
# it from a hard-coded donor list (001V, 027U, 029T), reclassifies those as DD_ALI and
# then drops them along with SAE. Depositing the raw sheets would ship that inference as
# something a reader must re-derive; depositing the resolved 26 samples matches how
# every other file in this record was built -- select at source, assert downstream.
#
# OUTPUT  pheno_4ct.csv — 26 rows, one per sample.
#   sample_name  <donor>_<treatment>_<cellType>, identical to the count-matrix columns
#   donor, treatment, cell_type   published names (AT2 / LAE / FB / MV)
#   bam          original alignment file, for traceability back to GEO
#
# Verified against the deposited short-read matrix: sample_name is exactly its column
# set, and the design is balanced 13 DMSO / 13 SMG1i.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

SRC  <- path.expand("~/Partners HealthCare Dropbox/Peter Castaldi/Research/NMD/full_dataset/pheno")
DEP  <- path.expand("~/claude_projects/nmd_deposit_2026/source_data")
OUT  <- file.path(DEP, "pheno_4ct.csv")
ALI_DONORS <- c("001V", "027U", "029T")   # submerged-vs-ALI split is not in the sheets

cat("=== reading the two phenotype sheets ===\n")
p1 <- fread(file.path(SRC, "FB_MV_pheno_2025.8.2.csv"), data.table = FALSE)
p2 <- fread(file.path(SRC, "DD_DO_AT_pheno_2025.8.15.csv"), data.table = FALSE)
cat(sprintf("  FB/MV %d rows | DD/DO/AT %d rows\n", nrow(p1), nrow(p2)))

p1$cell_type <- p1$ct
p2 <- p2[, setdiff(names(p2), "mapid")]
phe <- rbind(p1[, c("sample_id", "treatment", "cell_type", "bam")],
             p2[, c("sample_id", "treatment", "cell_type", "bam")])
cat(sprintf("  combined %d rows; cell types: %s\n", nrow(phe),
            paste(sort(unique(phe$cell_type)), collapse = ", ")))

# --- resolve the ALI cultures, then keep only the four manuscript cell types ---
is_ali <- phe$cell_type == "LAE" & phe$sample_id %in% ALI_DONORS
cat(sprintf("  LAE samples that are ALI cultures (donors %s): %d\n",
            paste(ALI_DONORS, collapse = ", "), sum(is_ali)))
phe <- phe[!is_ali & phe$cell_type %in% c("AT2", "LAE", "FB", "MV"), ]

phe$sample_name <- paste(phe$sample_id, phe$treatment, phe$cell_type, sep = "_")
phe <- phe[order(phe$cell_type, phe$sample_id, phe$treatment),
           c("sample_name", "sample_id", "treatment", "cell_type", "bam")]
names(phe)[names(phe) == "sample_id"] <- "donor"

cat("\n=== verification ===\n")
sr <- names(fread(file.path(DEP, "salmon_gene_counts_4ct.csv"), nrows = 0))[-1]
cat(sprintf("  rows                        : %d\n", nrow(phe)))
cat(sprintf("  cell types                  : %s\n", paste(sort(unique(phe$cell_type)), collapse = ", ")))
cat(sprintf("  treatment balance           : %s\n",
            paste(names(table(phe$treatment)), table(phe$treatment), sep = "=", collapse = " ")))
cat(sprintf("  matches short-read columns  : %s\n", setequal(phe$sample_name, sr)))
stopifnot("expected 26 manuscript samples"      = nrow(phe) == 26L,
          "sample_name must equal the deposited short-read column set" =
            setequal(phe$sample_name, sr),
          "design must be balanced 13/13"       = all(table(phe$treatment) == 13L),
          "unexpected cell type"                =
            setequal(unique(phe$cell_type), c("AT2", "LAE", "FB", "MV")))

fwrite(phe, OUT)
cat(sprintf("\nwrote %s (%.1f KB)\n", OUT, file.size(OUT) / 1024))
cat("verification passed.\n")
