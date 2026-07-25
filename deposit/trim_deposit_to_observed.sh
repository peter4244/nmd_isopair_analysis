#!/usr/bin/env bash
# =============================================================================
# trim_deposit_to_observed.sh — restrict every deposit file to the features
# actually observed in this dataset.
#
# INTERNAL packaging step. Not shipped; see DEPOSIT_PROVENANCE_INTERNAL.md.
#
# Rationale: the deposit should describe exactly this dataset. A feature with
# zero counts in every sample is not part of it, and leaving such rows in the
# files invites questions the paper does not answer.
#
# Verified result-neutral before running (2026-07-25):
#   isoforms  645,272 -> 614,992   filterByExpr set identical; % output lost
#                                  identical to 4 dp
#   genes      78,899 ->  46,571   filterByExpr set identical; the published
#                                  filtered set contains zero all-zero genes
#
# Every file is keyed on the same universe, so all are trimmed to match:
#   counts (isoform, gene) · SQANTI classification · FASTA · GTF · CDS GFF3
# =============================================================================
set -euo pipefail
D="$HOME/claude_projects/nmd_deposit_2026/source_data"
W="$HOME/claude_projects/nmd_deposit_2026/.trim_work"
mkdir -p "$W"

echo "=== 1. keep-sets from the count matrices ==="
Rscript -e '
suppressPackageStartupMessages(library(data.table))
D <- path.expand("~/claude_projects/nmd_deposit_2026/source_data")
W <- path.expand("~/claude_projects/nmd_deposit_2026/.trim_work")
i <- fread(file.path(D,"nmd_lungcells_counts_4ct.csv"), showProgress=FALSE)
keep_i <- i[[1]][rowSums(as.matrix(i[,-1])) > 0]
fwrite(data.table(keep_i), file.path(W,"keep_isoforms.txt"), col.names=FALSE)
fwrite(i[i[[1]] %in% keep_i], file.path(D,"nmd_lungcells_counts_4ct.csv"))
cat(sprintf("  isoforms %d -> %d\n", nrow(i), length(keep_i)))
g <- fread(file.path(D,"salmon_gene_counts_4ct.csv"), showProgress=FALSE)
keep_g <- g[[1]][rowSums(as.matrix(g[,-1])) > 0]
fwrite(g[g[[1]] %in% keep_g], file.path(D,"salmon_gene_counts_4ct.csv"))
cat(sprintf("  genes    %d -> %d\n", nrow(g), length(keep_g)))
'

echo "=== 2. SQANTI classification ==="
awk -F'\t' 'NR==FNR{k[$1];next} FNR==1||($1 in k)' \
    "$W/keep_isoforms.txt" "$D/sqanti/nmd_lungcells_classification.txt" > "$W/cls.tmp"
mv "$W/cls.tmp" "$D/sqanti/nmd_lungcells_classification.txt"
echo "  rows now: $(($(wc -l < "$D/sqanti/nmd_lungcells_classification.txt")-1))"

echo "=== 3. FASTA ==="
awk 'NR==FNR{k[$1];next}
     /^>/{id=substr($1,2); keep=(id in k)}
     keep' "$W/keep_isoforms.txt" "$D/sqanti/nmd_lungcells_corrected.fasta" > "$W/fa.tmp"
mv "$W/fa.tmp" "$D/sqanti/nmd_lungcells_corrected.fasta"
echo "  records now: $(grep -c '^>' "$D/sqanti/nmd_lungcells_corrected.fasta")"

# GTF/GFF3: keep any line whose transcript_id is retained. Feature lines without a
# transcript_id (e.g. GTF "gene" rows) are kept when the gene still has >=1 retained
# transcript, so gene records are not orphaned.
trim_gxf() { # $1=path — TWO-PASS so original line order is preserved
  local f="$1"
  # pass 1: which genes still have at least one retained transcript?
  awk -v keepfile="$W/keep_isoforms.txt" '
    BEGIN{ while((getline l < keepfile)>0) k[l] }
    match($0,/transcript_id "[^"]+"/){
      tid=substr($0,RSTART+15,RLENGTH-16)
      if (tid in k && match($0,/gene_id "[^"]+"/))
        print substr($0,RSTART+9,RLENGTH-10)
    }' "$f" | sort -u > "$W/okgenes.txt"
  # pass 2: emit in original order
  awk -v keepfile="$W/keep_isoforms.txt" -v genefile="$W/okgenes.txt" '
    BEGIN{ while((getline l < keepfile)>0) k[l]
           while((getline g < genefile)>0) okgene[g] }
    {
      tid=""; gid=""
      if (match($0,/transcript_id "[^"]+"/)) tid=substr($0,RSTART+15,RLENGTH-16)
      if (match($0,/gene_id "[^"]+"/))       gid=substr($0,RSTART+9,RLENGTH-10)
      if (tid!="")      { if (tid in k) print }
      else if (gid!="") { if (gid in okgene) print }
      else print                      # headers / comment lines
    }' "$f" > "$W/gxf.tmp"
  mv "$W/gxf.tmp" "$f"
}

echo "=== 4. filtered GTF ==="
trim_gxf "$D/sqanti/nmd_lungcells_filtered.gtf"
echo "  transcripts now: $(grep -o 'transcript_id "[^\"]*"' "$D/sqanti/nmd_lungcells_filtered.gtf" | sort -u | wc -l)"

echo "=== 5. corrected CDS GFF3 ==="
trim_gxf "$D/sqanti/nmd_lungcells_corrected.cds.gff3"
echo "  transcripts now: $(grep -o 'transcript_id "[^\"]*"' "$D/sqanti/nmd_lungcells_corrected.cds.gff3" | sort -u | wc -l)"

echo "=== 6. refresh checksums ==="
cd "$D" && find . -type f -print0 | sort -z | xargs -0 shasum -a 256 > ../MANIFEST.sha256
echo "  $(wc -l < ../MANIFEST.sha256) files"
du -sh "$D"
echo "TRIM_DONE"
