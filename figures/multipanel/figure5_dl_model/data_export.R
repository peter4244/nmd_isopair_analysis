#!/usr/bin/env Rscript
# ==============================================================================
# Figure 5 — n=190 gencode_all3 isoform list (for Panel G: uORF attention)
#
# Output:
#   data/gencode_all3_n190_isoforms.tsv  (one row per comparator isoform)
#     columns: arm, group, gene_id, reference_isoform_id,
#              comparator_isoform_id, has_ptc
#     - arm ∈ {NMD, Control}
#     - group ∈ {NMD+/PTC+, NMD+/PTC-, Control}
#
# Joins to results_4ct/uorf_attention_metrics.tsv on comparator_isoform_id.
#
# Scope identical to Figure 4 Section A (gencode_all3 / Subset 1):
#   pop_BC ∩ all-3-ENST ∩ coding-CDS, gene-matched between NMD c2 and Control c4
# Yields n=130 NMD pairs + n=130 Control pairs (4-CT re-scope, 2026-07-11;
# the "n190" in the output filename is the stable subset-1 identifier).
#
# Run from this folder:  Rscript data_export.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Isopair)
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

cat("Loading profiles + cds + structures...\n")
cds        <- as.data.table(readRDS(file.path(DM, "cds.rds")))
structures <- as.data.table(readRDS(file.path(DM, "structures.rds")))
profiles_c2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))

is_enst  <- function(x) grepl("^ENST", x)
make_key <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")

# ── pop_BC (gene-matched) ────────────────────────────────────────────────
shared <- intersect(unique(make_key(profiles_c2)), unique(make_key(profiles_c4)))
pop_BC_c2 <- profiles_c2[make_key(profiles_c2) %in% shared]
pop_BC_c4 <- profiles_c4[make_key(profiles_c4) %in% shared]

# ── Section A (gencode_all3): all 3 ENST + coding-CDS, gene-matched ─────
pop_BC_c2_E <- pop_BC_c2[is_enst(reference_isoform_id) & is_enst(comparator_isoform_id)]
pop_BC_c4_E <- pop_BC_c4[is_enst(reference_isoform_id) & is_enst(comparator_isoform_id)]
joint1 <- intersect(make_key(pop_BC_c2_E), make_key(pop_BC_c4_E))
pop_BC_c2_E <- pop_BC_c2_E[make_key(pop_BC_c2_E) %in% joint1]
pop_BC_c4_E <- pop_BC_c4_E[make_key(pop_BC_c4_E) %in% joint1]

coding_ids <- cds[coding_status == "coding", isoform_id]
have_cds <- function(d) d$reference_isoform_id %in% coding_ids &
                        d$comparator_isoform_id %in% coding_ids
pop_BC_c2_E2 <- pop_BC_c2_E[have_cds(pop_BC_c2_E)]
pop_BC_c4_E2 <- pop_BC_c4_E[have_cds(pop_BC_c4_E)]
joint2 <- intersect(make_key(pop_BC_c2_E2), make_key(pop_BC_c4_E2))
gencode_all3_c2 <- pop_BC_c2_E2[make_key(pop_BC_c2_E2) %in% joint2]
gencode_all3_c4 <- pop_BC_c4_E2[make_key(pop_BC_c4_E2) %in% joint2]

# 4-CT re-scope + 25% ref-share floor (2026-07-11): n=190 -> n=130
cat(sprintf("[gencode_all3]  NMD=%d  Control=%d  (expect 130/130 4-CT)\n",
            nrow(gencode_all3_c2), nrow(gencode_all3_c4)))
stopifnot(nrow(gencode_all3_c2) == 130L, nrow(gencode_all3_c4) == 130L)

# ── PTC determination on NMD pairs (own GENCODE stop + 50-nt rule) ──────
coding_cds  <- cds[coding_status == "coding"]
struct_idx  <- match(coding_cds$isoform_id, structures$isoform_id)

own_stop_tx_lookup <- (function() {
  stop_tx <- rep(NA_integer_, nrow(coding_cds))
  for (i in seq_len(nrow(coding_cds))) {
    si <- struct_idx[i]
    if (is.na(si)) next
    strand <- coding_cds$strand[i]
    if (is.na(strand) || strand == "") next
    stop_g <- if (strand == "+") coding_cds$cds_stop[i] else coding_cds$cds_start[i]
    stop_tx[i] <- Isopair::genomicToTranscript(
      stop_g, structures$exon_starts[[si]], structures$exon_ends[[si]], strand)
  }
  data.table(isoform_id = coding_cds$isoform_id, own_stop_tx_pos = stop_tx)
})()

junctions_for <- function(iso_id) {
  si <- match(iso_id, structures$isoform_id)
  if (is.na(si)) return(integer(0))
  exon_starts <- structures$exon_starts[[si]]
  exon_ends   <- structures$exon_ends[[si]]
  strand      <- structures$strand[si]
  if (is.na(strand) || strand == "") return(integer(0))
  exon_lens <- exon_ends - exon_starts + 1L
  if (strand == "-") exon_lens <- rev(exon_lens)
  if (length(exon_lens) <= 1L) return(integer(0))
  cumsum(exon_lens[-length(exon_lens)])
}

augment_and_label <- function(pop) {
  d <- merge(pop[, .(gene_id, reference_isoform_id, comparator_isoform_id)],
             own_stop_tx_lookup,
             by.x = "comparator_isoform_id", by.y = "isoform_id",
             all.x = TRUE)
  d[, last_ejc_tx_pos := vapply(
        comparator_isoform_id,
        function(iso) { j <- junctions_for(iso); if (length(j) == 0L) NA_integer_ else max(j) },
        integer(1))]
  d[, distance := last_ejc_tx_pos - own_stop_tx_pos]
  d[, has_ptc := !is.na(distance) & distance > 50L]
  d
}

nmd_lab <- augment_and_label(gencode_all3_c2)
ctrl_lab <- augment_and_label(gencode_all3_c4)

nmd_lab[, arm := "NMD"]
nmd_lab[, group := ifelse(has_ptc, "NMD+/PTC+", "NMD+/PTC-")]

ctrl_lab[, arm := "Control"]
ctrl_lab[, group := "Control"]

out <- rbind(
  nmd_lab[,  .(arm, group, gene_id, reference_isoform_id, comparator_isoform_id, has_ptc)],
  ctrl_lab[, .(arm, group, gene_id, reference_isoform_id, comparator_isoform_id, has_ptc)]
)

cat("\nGroup counts:\n")
print(out[, .N, by = group])

# Cross-check breakdown (report; guards locked to the 4-CT re-scope run 2026-07-11).
have <- out[, .N, by = group][order(group)]
cat("\nBreakdown (4-CT):\n"); print(have)
stopifnot(have[group == "NMD+/PTC+", N] == 48L,
          have[group == "NMD+/PTC-", N] == 82L,
          have[group == "Control",   N] == 130L)

fwrite(out, file.path(OUT_DIR, "gencode_all3_n190_isoforms.tsv"), sep = "\t")
cat(sprintf("\nWrote: %s  (%d rows)\n",
            file.path(OUT_DIR, "gencode_all3_n190_isoforms.tsv"), nrow(out)))

# ──────────────────────────────────────────────────────────────────────────
# Panel G — uORF attention proportion (NMD+/PTC+ vs NMD+/PTC- vs Control)
# Joins our n=190 list to results_4ct/uorf_attention_metrics.tsv on
# comparator_isoform_id. Uses ALL splits (train + val + test) per project
# convention: model interpretability uses all data; test-only goes in supp.
# ──────────────────────────────────────────────────────────────────────────
UORF_METRICS <- "/Users/petecastaldi/claude_projects/NMD_orf_model_v5_4ct/results_4ct/uorf_attention_metrics.tsv"
metrics <- fread(UORF_METRICS,
                 select = c("isoform_id", "uorf_attention_frac",
                            "split", "subgroup", "label", "prob"))

panG <- merge(out, metrics,
              by.x = "comparator_isoform_id", by.y = "isoform_id",
              all.x = TRUE)
have <- panG[!is.na(uorf_attention_frac)]

cat(sprintf("\n[Panel G] join coverage: %d / %d (%.1f%%)\n",
            nrow(have), nrow(panG), 100 * nrow(have) / nrow(panG)))
cat("Group coverage:\n")
print(have[, .(n = .N, median_frac = round(median(uorf_attention_frac), 4),
               q25 = round(quantile(uorf_attention_frac, 0.25), 4),
               q75 = round(quantile(uorf_attention_frac, 0.75), 4)),
           by = group][order(factor(group, levels = c("NMD+/PTC+","NMD+/PTC-","Control")))])

# Long table — one row per isoform (joined only)
long <- have[, .(group, comparator_isoform_id, uorf_attention_frac, split)]
fwrite(long, file.path(OUT_DIR, "panelG_uorf_attention_long.tsv"), sep = "\t")

# Descriptives
desc <- have[, .(n = .N,
                 median = median(uorf_attention_frac),
                 q25 = quantile(uorf_attention_frac, 0.25),
                 q75 = quantile(uorf_attention_frac, 0.75),
                 mean = mean(uorf_attention_frac),
                 sd = sd(uorf_attention_frac)),
             by = group]
desc <- desc[order(factor(group, levels = c("NMD+/PTC+","NMD+/PTC-","Control")))]
fwrite(desc, file.path(OUT_DIR, "panelG_uorf_attention_descriptives.tsv"), sep = "\t")

# Pairwise Wilcoxon + Cliff's δ + Hodges–Lehmann shift
groups <- c("NMD+/PTC+", "NMD+/PTC-", "Control")
pairs <- combn(groups, 2, simplify = FALSE)
cliff_delta <- function(x, y) {
  gx <- sum(outer(x, y, `>`))
  lx <- sum(outer(x, y, `<`))
  (gx - lx) / (length(x) * length(y))
}
pw <- rbindlist(lapply(pairs, function(p) {
  x <- have[group == p[1], uorf_attention_frac]
  y <- have[group == p[2], uorf_attention_frac]
  wt <- wilcox.test(x, y, exact = FALSE, conf.int = TRUE)
  data.table(group_x = p[1], group_y = p[2],
             nx = length(x), ny = length(y),
             wilcox_p = wt$p.value,
             hl_shift = unname(wt$estimate),
             cliff_delta = cliff_delta(x, y))
}))
fwrite(pw, file.path(OUT_DIR, "panelG_uorf_attention_pairwise.tsv"), sep = "\t")

cat("\nPairwise tests:\n")
print(pw)

cat(sprintf("\nWrote:\n  %s\n  %s\n  %s\n",
            file.path(OUT_DIR, "panelG_uorf_attention_long.tsv"),
            file.path(OUT_DIR, "panelG_uorf_attention_descriptives.tsv"),
            file.path(OUT_DIR, "panelG_uorf_attention_pairwise.tsv")))
