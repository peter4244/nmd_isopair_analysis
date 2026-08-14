suppressMessages(library(data.table))
# Resolve paths relative to THIS script, so it runs from any checkout. The first draft
# hardcoded a session scratchpad directory, which would have failed for anyone else.
args <- commandArgs(trailingOnly = FALSE)
HERE <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1])))
C <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr/analysis_cache"
x <- readRDS(file.path(C, "ref_atg_analysis.rds"))
ejc <- as.data.table(readRDS(file.path(HERE, "last_ejc.rds")))
bad <- c("mapping_failed", "ref_atg_lost")

arm <- function(d) as.data.table(d)[!category %in% bad & ref_source == "GENCODE"]
a2 <- arm(x$c2); a4 <- arm(x$c4)
k2 <- unique(a2[, .(gene_id, reference_isoform_id)]); k4 <- unique(a4[, .(gene_id, reference_isoform_id)])
sh <- merge(k2, k4, by = c("gene_id", "reference_isoform_id"))
i2 <- merge(a2, sh, by = c("gene_id", "reference_isoform_id"))[, comparison := "NMD"]
i4 <- merge(a4, sh, by = c("gene_id", "reference_isoform_id"))[, comparison := "Control"]

d <- rbind(i2, i4, fill = TRUE)
d <- merge(d, ejc, by.x = "comparator_isoform_id", by.y = "isoform_id", all.x = TRUE)
# SAME CONVENTION AS THE SHIPPED PANEL, verified on its 260 rows:
#   distance = last_ejc_tx_pos - stop_tx_pos ; positive = stop UPSTREAM of the last EJC
d[, distance := last_ejc_tx - comp_stop_tx_pos]
out <- d[!is.na(distance), .(comparator_isoform_id, gene_id, reference_isoform_id,
                             comparison, distance, category,
                             last_ejc_tx_pos = last_ejc_tx, own_stop_tx_pos = comp_stop_tx_pos)]
fwrite(out, file.path(HERE, "panelD_819.tsv"), sep = "\t")

cat(sprintf("rows written: %d\n", nrow(out)))
for (g in c("NMD", "Control")) {
  v <- out[comparison == g, distance]
  cat(sprintf("  %-8s n=%4d  median=%8.1f  >= +50nt: %5.1f%%   PTC category: %5.1f%%\n",
      g, length(v), median(v), 100*mean(v >= 50),
      100*mean(out[comparison==g, category] == "effectively_ptc")))
}
