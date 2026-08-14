suppressMessages(library(data.table))
# Resolve paths relative to THIS script, so it runs from any checkout. The first draft
# hardcoded a session scratchpad directory, which would have failed for anyone else.
args <- commandArgs(trailingOnly = FALSE)
HERE <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1])))
DM <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr"
s <- as.data.table(readRDS(file.path(DM, "structures.rds")))

# Last EJC in TRANSCRIPT coordinates = cumulative length of every exon but the last,
# taken in transcript order (genomic ascending on +, descending on -).
last_ejc <- function(st, en, strand) {
  if (length(st) < 2) return(NA_real_)
  o <- if (identical(strand, "-")) order(-st) else order(st)
  len <- (en - st + 1)[o]
  sum(len[seq_len(length(len) - 1L)])
}
s[, last_ejc_tx := mapply(last_ejc, exon_starts, exon_ends, strand)]

known <- fread("/Users/petecastaldi/claude_projects/nmd/figures/multipanel/figure3_isopair_and_ptc/data/panelD_stop_codon_distance.tsv")
chk <- merge(known[, .(comparator_isoform_id, shipped = last_ejc_tx_pos)],
             s[, .(comparator_isoform_id = isoform_id, mine = last_ejc_tx)], by = "comparator_isoform_id")
cat(sprintf("validation rows matched : %d of %d\n", nrow(chk), nrow(known)))
cat(sprintf("exact agreement         : %d (%.1f%%)\n", sum(chk$shipped == chk$mine),
            100*sum(chk$shipped == chk$mine)/nrow(chk)))
if (any(chk$shipped != chk$mine)) print(head(chk[shipped != mine], 4))
saveRDS(s[, .(isoform_id, last_ejc_tx)], file.path(HERE, "last_ejc.rds"))
