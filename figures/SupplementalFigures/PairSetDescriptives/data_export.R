#!/usr/bin/env Rscript
# =============================================================================
# Data export — pair-set descriptive numbers for the §4¶1 / §4¶2 supplements.
#
# CANONICAL SOURCE: this script mirrors the §1a / §1c chunks in the analysis
# report (`05_final_report_gencode_scope_2026-06-15.Rmd`) verbatim. The Rmd is
# the single source of truth; this script reproduces the same numbers so the
# figure can be rendered without rendering the full Rmd, but the computation
# logic is identical (gene-matched triplet construction; paired Wilcoxon for
# length; non-NMD-expressed-isoform filter for the count and reference-share).
#
# Outputs:
#   data/isoforms_per_gene.tsv            — per-gene non-NMD isoform count
#   data/ref_expression_fraction.tsv      — per-gene reference DMSO share
#   data/tx_length_by_role_long.tsv       — long form for the violin panel
#   data/descriptives_summary.tsv         — single-row summary
#   data/pairwise_tx_length.tsv           — paired Wilcoxon contrasts
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

HERE <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grep("--file=", args)]
  if (length(fa) > 0)
    return(dirname(normalizePath(sub("^--file=", "", fa[1]))))
  normalizePath(getwd())
})()
OUT_DIR <- file.path(HERE, "data")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DM <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr"

# ── Load (same objects as the Rmd setup chunk) ──
profiles_c2  <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
profiles_c4  <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))
structures   <- as.data.table(readRDS(file.path(DM, "structures.rds")))
expr_mat     <- readRDS(file.path(DM, "expression_data.rds"))
nmd_class    <- readRDS(file.path(DM, "nmd_classification.rds"))
dmso_samples <- readRDS(file.path(DM, "dmso_samples.rds"))

# ─── pop_BC construction — mirrors Rmd `sec1-pop-bc` chunk ───
pop_bc_shared <- merge(
  unique(profiles_c2[, .(gene_id, reference_isoform_id)]),
  unique(profiles_c4[, .(gene_id, reference_isoform_id)]),
  by = c("gene_id", "reference_isoform_id")
)
pop_bc_c2 <- merge(profiles_c2, pop_bc_shared,
                   by = c("gene_id", "reference_isoform_id"))
pop_bc_c4 <- merge(profiles_c4, pop_bc_shared,
                   by = c("gene_id", "reference_isoform_id"))

cat(sprintf("pop_BC: %d shared (gene, reference) pairs; %d NMD + %d Control comparator pairs\n",
            nrow(pop_bc_shared), nrow(pop_bc_c2), nrow(pop_bc_c4)))

# ─── §1a transcript-length — mirrors Rmd `sec1-tx-length` chunk ───
# (paired Wilcoxon on gene-matched triplets)
tx_lengths_all <- setNames(
  vapply(seq_len(nrow(structures)), function(i)
    sum(structures$exon_ends[[i]] - structures$exon_starts[[i]] + 1L),
    integer(1)),
  structures$isoform_id
)

pop_bc_triplets <- merge(
  pop_bc_c2[, .(gene_id, reference_isoform_id, comparator_isoform_id)],
  pop_bc_c4[, .(gene_id, reference_isoform_id, comparator_isoform_id)],
  by = c("gene_id", "reference_isoform_id"),
  suffixes = c("_nmd", "_ctrl")
)
pop_bc_triplets[, ref_len  := tx_lengths_all[reference_isoform_id]]
pop_bc_triplets[, nmd_len  := tx_lengths_all[comparator_isoform_id_nmd]]
pop_bc_triplets[, ctrl_len := tx_lengths_all[comparator_isoform_id_ctrl]]
pop_bc_triplets <- pop_bc_triplets[
  complete.cases(pop_bc_triplets[, .(ref_len, nmd_len, ctrl_len)])]

ref_lens  <- pop_bc_triplets$ref_len
nmd_lens  <- pop_bc_triplets$nmd_len
ctrl_lens <- pop_bc_triplets$ctrl_len

# Long form for the violin panel
tx_long <- rbindlist(list(
  data.table(role = "NMD comparator",     length_nt = nmd_lens),
  data.table(role = "Reference",          length_nt = ref_lens),
  data.table(role = "Control comparator", length_nt = ctrl_lens)
))
tx_long[, role := factor(role,
  levels = c("NMD comparator", "Reference", "Control comparator"))]
fwrite(tx_long, file.path(OUT_DIR, "tx_length_by_role_long.tsv"), sep = "\t")

# Paired Wilcoxon — three contrasts using the gene-matched triplet design
p_nmd_ref   <- wilcox.test(nmd_lens,  ref_lens,  paired = TRUE, exact = FALSE)$p.value
p_nmd_ctrl  <- wilcox.test(nmd_lens,  ctrl_lens, paired = TRUE, exact = FALSE)$p.value
p_ctrl_ref  <- wilcox.test(ctrl_lens, ref_lens,  paired = TRUE, exact = FALSE)$p.value
pw <- data.table(
  group_x  = c("NMD comparator", "NMD comparator", "Control comparator"),
  group_y  = c("Reference",      "Control comparator", "Reference"),
  n_pairs  = rep(nrow(pop_bc_triplets), 3),
  median_x = c(median(nmd_lens), median(nmd_lens), median(ctrl_lens)),
  median_y = c(median(ref_lens), median(ctrl_lens), median(ref_lens)),
  paired_wilcox_p = signif(c(p_nmd_ref, p_nmd_ctrl, p_ctrl_ref), 3)
)
fwrite(pw, file.path(OUT_DIR, "pairwise_tx_length.tsv"), sep = "\t")

# ─── §1c isoforms-per-gene — mirrors Rmd `sec1-iso-per-gene` chunk ───
nmd_iso_set <- nmd_class[["all_samples"]]$nmd
expr_iso    <- rownames(expr_mat)
sec1_iso    <- structures[
  gene_id %in% pop_bc_shared$gene_id &
    isoform_id %in% expr_iso &
    !isoform_id %in% nmd_iso_set,
  .(isoform_id, gene_id)]

iso_per_gene <- sec1_iso[, .(n_isoforms = .N), by = gene_id]
all_pop_genes <- data.table(gene_id = unique(pop_bc_shared$gene_id))
iso_per_gene <- merge(all_pop_genes, iso_per_gene, by = "gene_id", all.x = TRUE)
iso_per_gene[is.na(n_isoforms), n_isoforms := 0L]
fwrite(iso_per_gene, file.path(OUT_DIR, "isoforms_per_gene.tsv"), sep = "\t")

# ─── §1c reference-share — mirrors Rmd `sec1-ref-share` chunk ───
dmso_4ct_samples <- unique(c(dmso_samples$AT, dmso_samples$DD,
                              dmso_samples$FB, dmso_samples$MV))
stopifnot(all(dmso_4ct_samples %in% colnames(expr_mat)))

dmso_means <- rowMeans(expr_mat[, dmso_4ct_samples, drop = FALSE])
sec1_iso[, mean_dmso := dmso_means[isoform_id]]
sec1_iso[is.na(mean_dmso), mean_dmso := 0]
gene_total_dmso <- sec1_iso[, .(gene_total = sum(mean_dmso)), by = gene_id]

ref_share <- data.table(
  gene_id              = pop_bc_c2$gene_id,
  reference_isoform_id = pop_bc_c2$reference_isoform_id
)
ref_share <- merge(ref_share, gene_total_dmso, by = "gene_id", all.x = TRUE)
ref_share[, ref_dmso := dmso_means[reference_isoform_id]]
ref_share[is.na(ref_dmso), ref_dmso := 0]
ref_share[, ref_fraction_of_gene := ifelse(gene_total > 0,
                                            ref_dmso / gene_total, NA_real_)]
fwrite(ref_share, file.path(OUT_DIR, "ref_expression_fraction.tsv"), sep = "\t")

# ─── Descriptive summary ───
desc <- data.table(
  metric = c("n_genes_in_pop_BC",
             "median_isoforms_per_gene",
             "iqr_isoforms_per_gene",
             "median_ref_fraction_pct",
             "frac_refs_above_50pct",
             "median_nmd_comp_length_nt",
             "median_ref_length_nt",
             "median_control_comp_length_nt",
             "n_pop_bc_triplets",
             "paired_wilcox_p_nmd_vs_ref",
             "paired_wilcox_p_nmd_vs_ctrl",
             "paired_wilcox_p_ctrl_vs_ref"),
  value  = c(nrow(pop_bc_shared),
             median(iso_per_gene$n_isoforms),
             IQR(iso_per_gene$n_isoforms),
             round(100 * median(ref_share$ref_fraction_of_gene, na.rm = TRUE)),
             round(100 * mean(ref_share$ref_fraction_of_gene >= 0.50,
                              na.rm = TRUE)),
             median(nmd_lens),
             median(ref_lens),
             median(ctrl_lens),
             nrow(pop_bc_triplets),
             signif(p_nmd_ref,  3),
             signif(p_nmd_ctrl, 3),
             signif(p_ctrl_ref, 3))
)
fwrite(desc, file.path(OUT_DIR, "descriptives_summary.tsv"), sep = "\t")

cat("\n=== Descriptive summary (matches Rmd §1a/§1c) ===\n")
print(desc)
cat("\n=== Paired-Wilcoxon contrasts ===\n")
print(pw)
