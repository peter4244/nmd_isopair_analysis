#!/usr/bin/env Rscript
# =============================================================================
# Data export — pair-set descriptive numbers for Gap-§4-A + Gap-§4-B
#
# Three panels in the sibling SF (PairSetDescriptives):
#
#   A — Isoforms per gene (across the 3,009 genes in pop_BC).
#   B — Reference-isoform expression as a fraction of parent-gene expression
#       (DMSO baseline, 4-CT scope: AT, DD, FB, MV).
#   C — Transcript length, three roles within pop_BC:
#         NMD comparator (n = 3,009 from profiles_c2)
#         Reference      (n = 3,009 from either profiles_c2 or profiles_c4 — same set)
#         Control comparator (n = 3,009 from profiles_c4)
#
# Outputs:
#   data/isoforms_per_gene.tsv            — per-gene isoform count (3,009 rows)
#   data/ref_expression_fraction.tsv      — per-gene reference fraction (3,009 rows)
#   data/tx_length_by_role_long.tsv       — long form (9,027 rows) for panel C
#   data/descriptives_summary.tsv         — single-row table of medians cited in
#                                            Methods §1 and §4¶2
#   data/pairwise_tx_length.tsv           — pairwise Wilcoxon for panel C
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

# ── Load ──
profiles_c2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))

# Build the canonical pop_BC by intersecting the two arms on (gene_id,
# reference_isoform_id) — the same construction used in
# figures/multipanel/figure4_ptcneg_and_model/data_export.R.
make_key <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")
shared <- intersect(unique(make_key(profiles_c2)), unique(make_key(profiles_c4)))
profiles_c2 <- profiles_c2[make_key(profiles_c2) %in% shared]
profiles_c4 <- profiles_c4[make_key(profiles_c4) %in% shared]
# Order both arms on the shared key so paired sanity-checks below align.
setorder(profiles_c2, gene_id, reference_isoform_id)
setorder(profiles_c4, gene_id, reference_isoform_id)
expr   <- readRDS(file.path(DM, "expression_data.rds"))
gm     <- as.data.table(readRDS(file.path(DM, "gene_map.rds")))
ds     <- readRDS(file.path(DM, "dmso_samples.rds"))

# 4-CT DMSO sample union (AT + DD + FB + MV; excludes DD_ALI / DO_ALI per
# manuscript scope; see CLAUDE.md "Cell-type abbreviations").
dmso_4ct <- unique(c(ds$AT, ds$DD, ds$FB, ds$MV))
stopifnot(all(dmso_4ct %in% colnames(expr)))

cat(sprintf("pop_BC: profiles_c2 = %d rows; profiles_c4 = %d rows\n",
            nrow(profiles_c2), nrow(profiles_c4)))
cat(sprintf("DMSO 4-CT samples: %d (AT=%d DD=%d FB=%d MV=%d)\n",
            length(dmso_4ct), length(ds$AT), length(ds$DD),
            length(ds$FB), length(ds$MV)))

# ── Universe of genes (= pop_BC) ──
# Both arms must agree on the (gene_id, reference_isoform_id) key for a gene
# to be in pop_BC; profiles_c2 already gives us 3,009 such gene-reference
# entries (one row per gene because each gene has one reference isoform).
pop_genes <- unique(profiles_c2$gene_id)
cat(sprintf("Genes in pop_BC: %d\n", length(pop_genes)))
stopifnot(length(pop_genes) == nrow(profiles_c2))  # exactly 1 row per gene

# ── Panel A: isoforms per gene ──
# Matches the Isopair pipeline definition (05_final_report_mashr.Rmd §2):
# count isoforms per gene that are (a) in `structures` and `expression_data`,
# (b) excluded if they are NMD-sensitive (nmd_classification$all_samples$nmd).
# No additional relative-expression threshold — "in the expression matrix"
# is the filter.
structures <- as.data.table(readRDS(file.path(DM, "structures.rds")))
nmd_class  <- readRDS(file.path(DM, "nmd_classification.rds"))
nmd_iso    <- nmd_class$all_samples$nmd

exprs_iso <- rownames(expr)
sec1_iso  <- structures[gene_id %in% pop_genes &
                        isoform_id %in% exprs_iso &
                        !isoform_id %in% nmd_iso,
                        .(isoform_id, gene_id)]

iso_per_gene <- sec1_iso[, .(n_isoforms = .N), by = gene_id]
all_genes_dt <- data.table(gene_id = pop_genes)
iso_per_gene <- merge(all_genes_dt, iso_per_gene, by = "gene_id", all.x = TRUE)
iso_per_gene[is.na(n_isoforms), n_isoforms := 0L]

fwrite(iso_per_gene, file.path(OUT_DIR, "isoforms_per_gene.tsv"), sep = "\t")

# ── Panel B: reference-isoform fraction of parent gene ──
# Same denominator: sum of DMSO 4-CT mean expression across all non-NMD
# expressed isoforms of the gene (sec1_iso above). Numerator: the reference
# isoform's DMSO 4-CT mean.
dmso_mean <- rowMeans(expr[, dmso_4ct, drop = FALSE])
sec1_iso[, mean_dmso := dmso_mean[isoform_id]]
sec1_iso[is.na(mean_dmso), mean_dmso := 0]

gene_tot <- sec1_iso[, .(gene_total = sum(mean_dmso)), by = gene_id]
ref_fracs <- merge(
  data.table(gene_id              = profiles_c2$gene_id,
             reference_isoform_id = profiles_c2$reference_isoform_id),
  gene_tot, by = "gene_id", all.x = TRUE)
ref_fracs[, ref_dmso := dmso_mean[reference_isoform_id]]
ref_fracs[is.na(ref_dmso), ref_dmso := 0]
ref_fracs[, ref_fraction_of_gene := ifelse(gene_total > 0,
                                            ref_dmso / gene_total, NA_real_)]
fwrite(ref_fracs, file.path(OUT_DIR, "ref_expression_fraction.tsv"), sep = "\t")

# ── Panel C: transcript-length descriptives ──
# `length_ref` / `length_comp` in profiles_* are GENOMIC SPAN (tx_end − tx_start),
# NOT the spliced transcript length. We compute the spliced tx length from
# structures.rds = sum of exon lengths.
spliced_len <- vapply(seq_len(nrow(structures)), function(i) {
  s <- structures$exon_starts[[i]]; e <- structures$exon_ends[[i]]
  if (is.null(s) || length(s) == 0L) return(NA_integer_)
  as.integer(sum(e - s + 1L))
}, integer(1))
tx_len_lookup <- setNames(spliced_len, structures$isoform_id)

nmd_comp_len  <- tx_len_lookup[profiles_c2$comparator_isoform_id]
ref_len       <- tx_len_lookup[profiles_c2$reference_isoform_id]
ctrl_comp_len <- tx_len_lookup[profiles_c4$comparator_isoform_id]

# Sanity: same reference set on both arms after the pop_BC intersection
stopifnot(identical(profiles_c2$reference_isoform_id,
                    profiles_c4$reference_isoform_id))

tx_long <- rbindlist(list(
  data.table(role = "NMD comparator",     length_nt = nmd_comp_len),
  data.table(role = "Reference",          length_nt = ref_len),
  data.table(role = "Control comparator", length_nt = ctrl_comp_len)
))
tx_long[, role := factor(role,
  levels = c("NMD comparator", "Reference", "Control comparator"))]
fwrite(tx_long, file.path(OUT_DIR, "tx_length_by_role_long.tsv"), sep = "\t")

# Pairwise Wilcoxon
pairs <- list(c("NMD comparator", "Reference"),
              c("NMD comparator", "Control comparator"),
              c("Reference",      "Control comparator"))
pw <- rbindlist(lapply(pairs, function(p) {
  x <- tx_long[role == p[1], length_nt]
  y <- tx_long[role == p[2], length_nt]
  wt <- wilcox.test(x, y, exact = FALSE)
  data.table(group_x = p[1], group_y = p[2],
             n_x = length(x), n_y = length(y),
             median_x = median(x), median_y = median(y),
             wilcox_p = signif(wt$p.value, 3))
}))
fwrite(pw, file.path(OUT_DIR, "pairwise_tx_length.tsv"), sep = "\t")

# Kruskal–Wallis (test name that should be reported in the Methods if all
# three roles are compared simultaneously)
kw <- kruskal.test(length_nt ~ role, data = tx_long)

# ── Descriptive summary (for Methods + manuscript find/replace) ──
desc <- data.table(
  metric              = c("n_genes_in_pop_BC",
                          "median_isoforms_per_gene",
                          "iqr_isoforms_per_gene",
                          "median_ref_fraction_pct",
                          "frac_refs_above_50pct",
                          "median_nmd_comp_length_nt",
                          "median_ref_length_nt",
                          "median_control_comp_length_nt",
                          "kruskal_wallis_p"),
  value               = c(length(pop_genes),
                          median(iso_per_gene$n_isoforms),
                          IQR(iso_per_gene$n_isoforms),
                          round(100 * median(ref_fracs$ref_fraction_of_gene), 1),
                          round(mean(ref_fracs$ref_fraction_of_gene >= 0.50), 3),
                          median(nmd_comp_len),
                          median(ref_len),
                          median(ctrl_comp_len),
                          signif(kw$p.value, 3))
)
fwrite(desc, file.path(OUT_DIR, "descriptives_summary.tsv"), sep = "\t")

cat("\n=== Descriptive summary ===\n")
print(desc)

cat("\n=== Pairwise tx-length contrasts ===\n")
print(pw)

cat(sprintf("\nManuscript-cited medians (compare against §4¶2):\n"))
cat(sprintf("  NMD comparator   = %d nt  (manuscript: 3,049)\n",
            as.integer(median(nmd_comp_len))))
cat(sprintf("  Reference        = %d nt  (manuscript: 2,893)\n",
            as.integer(median(ref_len))))
cat(sprintf("  Control comp.    = %d nt  (manuscript: 2,762)\n",
            as.integer(median(ctrl_comp_len))))
cat(sprintf("  Kruskal–Wallis p = %s\n", format(kw$p.value, digits = 3)))

cat(sprintf("\nManuscript-cited descriptives (compare against §4¶1):\n"))
cat(sprintf("  median isoforms / gene             = %d  (manuscript: 7)\n",
            as.integer(median(iso_per_gene$n_isoforms))))
cat(sprintf("  median ref-isoform pct of gene     = %.1f%%  (manuscript: 70%%)\n",
            100 * median(ref_fracs$ref_fraction_of_gene)))
cat(sprintf("  fraction of refs above 50%% of gene = %.1f%%  (manuscript: 75%%)\n",
            100 * mean(ref_fracs$ref_fraction_of_gene >= 0.50)))
