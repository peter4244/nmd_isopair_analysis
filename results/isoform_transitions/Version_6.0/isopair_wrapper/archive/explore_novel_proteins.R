data_dir <- "data_mashr"
structures <- readRDS(file.path(data_dir, "structures.rds"))
cds <- readRDS(file.path(data_dir, "cds.rds"))
ptc <- readRDS(file.path(data_dir, "ptc.rds"))
fc_c4 <- readRDS(file.path(data_dir, "analysis_cache", "fc_c4_allsamples.rds"))
fc_c2 <- readRDS(file.path(data_dir, "analysis_cache", "fc_c2_allsamples.rds"))
np <- readRDS(file.path(data_dir, "analysis_cache", "novel_protein_c4_v2_allsamples.rds"))

# Gene symbol lookup
de_dir <- "/Users/petecastaldi/claude_projects/nmd/isocall_dge/mashr"
gene_sym <- character(0)
for (f in list.files(de_dir, pattern = "^nmd_mashr_die_.*\\.csv$", full.names = TRUE)) {
  tmp <- read.csv(f, stringsAsFactors = FALSE)[, c("gene_id", "hgnc_symbol")]
  tmp <- tmp[!duplicated(tmp$gene_id) & !is.na(tmp$hgnc_symbol) & tmp$hgnc_symbol != "", ]
  new <- setdiff(tmp$gene_id, names(gene_sym))
  gene_sym[tmp$gene_id[tmp$gene_id %in% new]] <- tmp$hgnc_symbol[tmp$gene_id %in% new]
}

cat("=== Novel protein candidates (C4 same_start_frameshift, no PTC) ===\n")
cat("Total pairs:", nrow(np), "\n")
cat("Unique genes:", length(unique(np$gene_id)), "\n\n")

np$gene_symbol <- gene_sym[np$gene_id]
np$is_novel_isoform <- grepl("\\.novel[0-9]+$", np$comparator_isoform_id)

cat("Novel vs ENST isoforms:\n")
print(table(novel = np$is_novel_isoform))
cat("\n")

cat("Summary of novel protein metrics:\n")
print(summary(np[, c("orf_length", "novel_aa", "total_aa", "pct_novel")]))
cat("\n")

# Top candidates: large proteins with substantial novel sequence
np_interesting <- np[np$pct_novel >= 20 & np$total_aa >= 100, ]
np_interesting <- np_interesting[order(-np_interesting$pct_novel), ]
cat("Candidates with >=20% novel sequence and >=100 aa total (n =",
    nrow(np_interesting), "):\n")
print(head(np_interesting[, c("gene_symbol", "gene_id", "comparator_isoform_id",
                               "total_aa", "novel_aa", "pct_novel", "is_novel_isoform")], 30))

cat("\n\n=== Also check diff_start_diff_frame without PTC (C2 NMD) ===\n")
ps2 <- fc_c2$pair_summary
ds2 <- ps2[ps2$frame_category %in% c("diff_start_diff_frame",
                                       "diff_start_diff_frame_with_frameshift"), ]
ds2_merged <- merge(ds2, ptc[, c("isoform_id", "has_ptc", "orf_length")],
                    by.x = "comparator_isoform_id", by.y = "isoform_id", all.x = TRUE)
ds2_noptc <- ds2_merged[!is.na(ds2_merged$has_ptc) & ds2_merged$has_ptc == FALSE, ]
cat("C2 diff_start_diff_frame without PTC:", nrow(ds2_noptc), "\n")
cat("Unique genes:", length(unique(ds2_noptc$gene_id)), "\n")
ds2_noptc$gene_symbol <- gene_sym[ds2_noptc$gene_id]
ds2_noptc <- ds2_noptc[order(-ds2_noptc$pct_shared_cds_frameshifted), ]
cat("\nTop diff-start novel protein genes:\n")
print(head(ds2_noptc[, c("gene_symbol", "gene_id", "comparator_isoform_id",
                          "frame_category", "pct_shared_cds_frameshifted",
                          "orf_length")], 20))

# What are the ENST IDs? We can look these up in UniProt
cat("\n\n=== ENST isoform IDs for UniProt lookup ===\n")
enst_novel <- np_interesting[!np_interesting$is_novel_isoform, ]
cat("ENST isoforms with novel protein (>=20%, >=100aa):", nrow(enst_novel), "\n")
# Strip version for UniProt lookup
enst_novel$enst_base <- sub("\\.[0-9]+$", "", enst_novel$comparator_isoform_id)
cat("\nTop ENST candidates:\n")
print(head(enst_novel[, c("gene_symbol", "enst_base", "total_aa", "novel_aa", "pct_novel")], 20))
