#!/usr/bin/env bash
# =============================================================================
# trim_deposit_to_4ct_observed.sh — restrict the LONG-READ feature universe to the
# transcripts actually observed in the 26 four-cell-type samples.
#
# REPLACES trim_deposit_to_observed.sh, which is retired. Read this header before
# changing either: the difference between them is one published number.
#
# WHAT THE OLD SCRIPT GOT RIGHT, AND WHAT IT GOT WRONG.
#
#   RIGHT (isoform side): the SQANTI transcript universe was discovered across 38
#   samples spanning SIX cell types. Measured 2026-07-27: of the 30,280 transcripts
#   with zero counts in the 26 4-CT samples, 30,103 (99.4%) are nonzero ONLY in
#   DO/DD_ALI -- the two cell types this project excludes. They exist solely because
#   of samples the paper says are not in it. 614,992 IS the 4-CT-native universe.
#   Pete, 2026-07-27: the deposit should look "as if the other two don't exist".
#
#   WRONG (gene side): it applied the same predicate to the SHORT-READ gene matrix,
#   independently. That is a bug under ANY universe policy. The salmon gene matrix is
#   quantified against a REFERENCE transcriptome, so its 78,899-gene universe is not
#   discovered from these samples and is not cell-type dependent. Dropping genes that
#   are all-zero in short-read deleted 735 genes that ARE observed in the long-read
#   data, and correlation_analysis.Rmd:226-229 keeps a gene expressed in EITHER
#   platform. That is what drove the concordance gene count down to 27,916.
#
# SO: trim the isoform side, never the gene side.
#
#   isoform matrix + SQANTI annotation : 645,272 -> 614,992   (4-CT observed)
#   salmon gene matrix                 : 78,899  -> 78,899    (UNTOUCHED, reference)
#
# EXPECTED CONSEQUENCE, stated up front so it is not mistaken for a regression:
# correlation_analysis.Rmd reports 28,651 common genes, NOT the published 29,185.
# That is the point -- 29,185 depends on transcripts contributed by the excluded cell
# types, and the manuscript number changes to match the 4-CT-native deposit.
# =============================================================================
set -euo pipefail
D="$HOME/claude_projects/nmd_deposit_2026/source_data"
W="$HOME/claude_projects/nmd_deposit_2026/.trim4ct_work"
mkdir -p "$W"

echo "=== 0. pre-flight: refuse unless the deposit is currently FULL ==="
iso_rows=$(( $(wc -l < "$D/nmd_lungcells_counts_4ct.csv") - 1 ))
gene_rows=$(( $(wc -l < "$D/salmon_gene_counts_4ct.csv") - 1 ))
echo "  isoform matrix rows: $iso_rows | gene matrix rows: $gene_rows"
if [ "$iso_rows" != 645272 ]; then
  echo "  ABORT: expected the untrimmed 645,272 isoform matrix. Run revert_trim_to_full.sh first."
  exit 1
fi
if [ "$gene_rows" != 78899 ]; then
  echo "  ABORT: the gene matrix is not the full 78,899. It must NOT be trimmed;"
  echo "         run revert_trim_to_full.sh to restore it."
  exit 1
fi

echo "=== 1. keep-set: transcripts observed in the 26 4-CT samples ==="
Rscript -e '
suppressPackageStartupMessages(library(data.table))
D <- path.expand("~/claude_projects/nmd_deposit_2026/source_data")
W <- path.expand("~/claude_projects/nmd_deposit_2026/.trim4ct_work")
i <- fread(file.path(D,"nmd_lungcells_counts_4ct.csv"), showProgress=FALSE)
keep_i <- i[[1]][rowSums(as.matrix(i[,-1])) > 0]
fwrite(data.table(keep_i), file.path(W,"keep_isoforms.txt"), col.names=FALSE)
fwrite(i[i[[1]] %in% keep_i], file.path(D,"nmd_lungcells_counts_4ct.csv"))
cat(sprintf("  isoforms %d -> %d\n", nrow(i), length(keep_i)))
cat("  gene matrix: NOT TOUCHED (reference-derived universe; trimming it is the bug)\n")
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

trim_gxf() { # $1=path — TWO-PASS so original line order is preserved
  local f="$1"
  awk -v keepfile="$W/keep_isoforms.txt" '
    BEGIN{ while((getline l < keepfile)>0) k[l] }
    match($0,/transcript_id "[^"]+"/){
      tid=substr($0,RSTART+15,RLENGTH-16)
      if (tid in k && match($0,/gene_id "[^"]+"/))
        print substr($0,RSTART+9,RLENGTH-10)
    }' "$f" | sort -u > "$W/okgenes.txt"
  awk -v keepfile="$W/keep_isoforms.txt" -v genefile="$W/okgenes.txt" '
    BEGIN{ while((getline l < keepfile)>0) k[l]
           while((getline g < genefile)>0) okgene[g] }
    {
      tid=""; gid=""
      if (match($0,/transcript_id "[^"]+"/)) tid=substr($0,RSTART+15,RLENGTH-16)
      if (match($0,/gene_id "[^"]+"/))       gid=substr($0,RSTART+9,RLENGTH-10)
      if (tid!="")      { if (tid in k) print }
      else if (gid!="") { if (gid in okgene) print }
      else print
    }' "$f" > "$W/gxf.tmp"
  mv "$W/gxf.tmp" "$f"
}

echo "=== 4. filtered GTF ==="
trim_gxf "$D/sqanti/nmd_lungcells_filtered.gtf"
echo "  transcripts now: $(grep -o 'transcript_id "[^\"]*"' "$D/sqanti/nmd_lungcells_filtered.gtf" | sort -u | wc -l | tr -d ' ')"

echo "=== 5. corrected CDS GFF3 ==="
trim_gxf "$D/sqanti/nmd_lungcells_corrected.cds.gff3"
echo "  transcripts now: $(grep -o 'transcript_id \"[^\"]*\"' "$D/sqanti/nmd_lungcells_corrected.cds.gff3" | sort -u | wc -l | tr -d ' ')"

echo "=== 6. verify: isoform side trimmed, GENE SIDE UNTOUCHED ==="
ok=1
chk(){ printf '  %-34s %-9s (want %s)  %s\n' "$1" "$2" "$3" \
       "$([ "$2" = "$3" ] && echo OK || { ok=0; echo MISMATCH; })"; }
chk "isoform matrix rows" "$(( $(wc -l < "$D/nmd_lungcells_counts_4ct.csv") - 1 ))" 614992
chk "gene matrix rows"    "$(( $(wc -l < "$D/salmon_gene_counts_4ct.csv") - 1 ))"   78899
chk "classification rows" "$(( $(wc -l < "$D/sqanti/nmd_lungcells_classification.txt") - 1 ))" 614992
chk "fasta records"       "$(grep -c '^>' "$D/sqanti/nmd_lungcells_corrected.fasta")" 614992
[ "$ok" = 1 ] || { echo "  ABORT: unexpected dimensions."; exit 1; }

echo "=== 7. refresh checksums ==="
( cd "$D" && find . -type f -print0 | sort -z | xargs -0 shasum -a 256 ) > "$D/../MANIFEST.sha256"
echo "  $(wc -l < "$D/../MANIFEST.sha256") files"
du -sh "$D"
echo "TRIM_4CT_DONE — now regenerate the derived DGELists and rebuild zenodo_upload/"
