#!/usr/bin/env Rscript
# ==============================================================================
# Figure 5 — n=1,166 Subset 2 (ref-AUG-traceable) isoform list + Path B-strict
# uORF attention computation.
#
# This is the broader pairwise scope from Figure 4 Section C — gene-matched
# 1:1 (NMD c2 vs Control c4), ENST-reference, ref-AUG-traceable (i.e.,
# mechanism_class ∈ {effectively_ptc, no_downstream_ejc, truncated_no_ejc}).
# Total: 1,166 NMD comparator pairs + 1,166 Control comparator pairs.
#
# Subgroup structure at this scope:
#   NMD+/PTC+         = category == "effectively_ptc"  (~1,050)
#   NMD+/PTC- retained = category ∈ {no_downstream_ejc, truncated_no_ejc}  (~116)
#   Control            = all 1,166 c4 comparators
#   (NMD+/PTC- lost not in this scope by construction — ref AUG is required
#   to project, and "lost" means it doesn't project)
#
# Joins to results_4ct/uorf_attention_metrics.tsv (v1) and applies
# Path B-strict structural reclassification locally using the priority-slot
# data (uorf_features_in_priority_slots.tsv) where available.
#
# Outputs:
#   data/subset2_n1166_isoforms.tsv       (one row per comparator with group)
#   data/panelG_uorf_attention_n1166.tsv  (per-isoform v1 + Path B-strict)
# ==============================================================================

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
CACHE <- file.path(DM, "analysis_cache")
LIB <- normalizePath(file.path(HERE, "..", "..", "lib"))
source(file.path(LIB, "mechanism_class.R"))

MODEL_RES <- "/Users/petecastaldi/claude_projects/NMD_orf_model_v5_4ct/results_4ct"

# ── Load ──
cat("Loading Isopair profiles + ref_atg cache...\n")
ref_atg     <- readRDS(file.path(CACHE, "ref_atg_analysis.rds"))
profiles_c2 <- as.data.table(readRDS(file.path(DM, "profiles_c2_allsamples.rds")))
profiles_c4 <- as.data.table(readRDS(file.path(DM, "profiles_c4_allsamples.rds")))

is_enst  <- function(x) grepl("^ENST", x)
make_key <- function(d) paste(d$gene_id, d$reference_isoform_id, sep = "::")

# ── Build pop_BC (gene-matched) ──
shared <- intersect(unique(make_key(profiles_c2)), unique(make_key(profiles_c4)))
pop_BC_c2 <- profiles_c2[make_key(profiles_c2) %in% shared]
pop_BC_c4 <- profiles_c4[make_key(profiles_c4) %in% shared]

# ── Section C scope: ENST-ref + ref-AUG-traceable, re-intersected 1:1 ──
sec_C_c2_base <- pop_BC_c2[is_enst(reference_isoform_id)]
sec_C_c4_base <- pop_BC_c4[is_enst(reference_isoform_id)]
c2_full <- as.data.table(ref_atg$c2)[comparator_isoform_id %in% sec_C_c2_base$comparator_isoform_id]
c4_full <- as.data.table(ref_atg$c4)[comparator_isoform_id %in% sec_C_c4_base$comparator_isoform_id]
TR <- c("effectively_ptc", "no_downstream_ejc", "truncated_no_ejc")
c2_tr <- c2_full[category %in% TR]
c4_tr <- c4_full[category %in% TR]
shared_keys <- merge(
  unique(c2_tr[, .(gene_id, reference_isoform_id)]),
  unique(c4_tr[, .(gene_id, reference_isoform_id)]),
  by = c("gene_id", "reference_isoform_id"))
c2_n1166 <- merge(c2_tr, shared_keys, by = c("gene_id", "reference_isoform_id"))
c4_n1166 <- merge(c4_tr, shared_keys, by = c("gene_id", "reference_isoform_id"))

# 4-CT re-scope + 25% ref-share floor (2026-07-11): n=1,166 -> n=819
cat(sprintf("[Subset 2] NMD c2 = %d, Control c4 = %d (expect 819/819 4-CT)\n",
            nrow(c2_n1166), nrow(c4_n1166)))
stopifnot(nrow(c2_n1166) == 819L, nrow(c4_n1166) == 819L)

# ── Assign group labels ──
c2_n1166[, group := fifelse(category == "effectively_ptc",
                             "NMD+/PTC+", "NMD+/PTC- retained")]
c2_n1166[, arm := "NMD"]
c4_n1166[, group := "Control"]
c4_n1166[, arm := "Control"]

out <- rbind(
  c2_n1166[, .(arm, group, gene_id, reference_isoform_id, comparator_isoform_id, category)],
  c4_n1166[, .(arm, group, gene_id, reference_isoform_id, comparator_isoform_id, category)]
)

cat("\nSubgroup counts:\n")
print(out[, .N, by = group])

fwrite(out, file.path(OUT_DIR, "subset2_n1166_isoforms.tsv"), sep = "\t")

# ── Join with model uORF metrics + priority-slot file ──
mtx <- fread(file.path(MODEL_RES, "uorf_attention_metrics.tsv"),
             select = c("isoform_id", "uorf_attention_frac", "top_is_uorf",
                        "main_orf_rank", "split", "label", "prob"))
slot <- fread(file.path(MODEL_RES, "uorf_features_in_priority_slots.tsv"))
tx   <- fread(file.path(MODEL_RES, "tx_summary.tsv"),
              select = c("isoform_id", "tx_length"))

# Path B-strict per slot
slot <- merge(slot, tx, by = "isoform_id")
slot[, orf_start_approx := frac_position * tx_length]
slot[, frac_stop := (orf_start_approx + orf_length) / tx_length]
slot[, is_strict := orf_length < 200 &
                     frac_position < 0.5 &
                     frac_stop < 0.6 &
                     kozak_score >= 1]

# Per isoform: sum attention over Path B-strict slots
# Isoforms NOT in slot file → strict attention = 0 (no flagged-uORF attention by v1)
pb <- slot[, .(strict_attn = sum(attention[is_strict])), by = isoform_id]

# Join: n=1,166 list + v1 metrics + strict attention
panG_n1166 <- merge(out, mtx, by.x = "comparator_isoform_id", by.y = "isoform_id",
                    all.x = TRUE)
panG_n1166 <- merge(panG_n1166, pb, by.x = "comparator_isoform_id", by.y = "isoform_id",
                    all.x = TRUE)
panG_n1166[is.na(strict_attn), strict_attn := 0]

cat(sprintf("\nModel coverage: %d / %d (%.1f%%) have v1 metrics from model\n",
            sum(!is.na(panG_n1166$uorf_attention_frac)),
            nrow(panG_n1166),
            100 * mean(!is.na(panG_n1166$uorf_attention_frac))))

# ── Subgroup table ──
gord <- c("NMD+/PTC+", "NMD+/PTC- retained", "Control")
have <- panG_n1166[!is.na(uorf_attention_frac)]

cat("\n=== Path B-strict vs v1: uORF attention by subgroup at n=1,166 ===\n\n")
tab <- have[, .(
  n              = .N,
  v1_mean        = round(mean(uorf_attention_frac), 3),
  v1_median      = round(median(uorf_attention_frac), 3),
  strict_mean    = round(mean(strict_attn), 3),
  strict_median  = round(median(strict_attn), 3),
  v1_top_is_uorf_pct     = round(100 * mean(top_is_uorf, na.rm = TRUE), 1),
  strict_any_attn_pct    = round(100 * mean(strict_attn > 0), 1),
  strict_attn_gt05_pct   = round(100 * mean(strict_attn > 0.05), 1),
  strict_attn_gt20_pct   = round(100 * mean(strict_attn > 0.20), 1)
), by = group][order(factor(group, levels = gord))]
print(tab)

cat("\n=== Pairwise Wilcoxon of Path B-strict attention ===\n")
pairs <- combn(gord, 2, simplify = FALSE)
pw <- rbindlist(lapply(pairs, function(p) {
  x <- have[group == p[1], strict_attn]
  y <- have[group == p[2], strict_attn]
  wt <- wilcox.test(x, y, exact = FALSE)
  data.table(group_x = p[1], group_y = p[2],
             nx = length(x), ny = length(y),
             wilcox_p = signif(wt$p.value, 3))
}))
print(pw)

# Descriptives — for the panel labels
desc <- have[, .(n = .N,
                 median = round(median(strict_attn), 4),
                 q25 = round(quantile(strict_attn, 0.25), 4),
                 q75 = round(quantile(strict_attn, 0.75), 4),
                 mean = round(mean(strict_attn), 4),
                 pct_any = round(100 * mean(strict_attn > 0), 1),
                 pct_gt05 = round(100 * mean(strict_attn > 0.05), 1),
                 pct_gt20 = round(100 * mean(strict_attn > 0.20), 1)),
             by = group][order(factor(group, levels = gord))]

fwrite(panG_n1166, file.path(OUT_DIR, "panelG_uorf_attention_n1166.tsv"), sep = "\t")
fwrite(pw,         file.path(OUT_DIR, "panelG_uorf_attention_n1166_pairwise.tsv"), sep = "\t")
fwrite(desc,       file.path(OUT_DIR, "panelG_uorf_attention_n1166_descriptives.tsv"), sep = "\t")

cat("\nWrote:\n  ", file.path(OUT_DIR, "subset2_n1166_isoforms.tsv"), "\n",
    "  ", file.path(OUT_DIR, "panelG_uorf_attention_n1166.tsv"), "\n",
    "  ", file.path(OUT_DIR, "panelG_uorf_attention_n1166_pairwise.tsv"), "\n",
    "  ", file.path(OUT_DIR, "panelG_uorf_attention_n1166_descriptives.tsv"), "\n",
    sep="")
