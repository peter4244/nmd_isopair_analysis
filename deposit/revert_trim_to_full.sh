#!/usr/bin/env bash
# =============================================================================
# revert_trim_to_full.sh — undo trim_deposit_to_observed.sh.
#
# INTERNAL packaging step, like the script it reverses. Not shipped.
#
# WHY. trim_deposit_to_observed.sh dropped every all-zero feature from each count
# matrix INDEPENDENTLY. correlation_analysis.Rmd:226-229 keeps a gene expressed in
# EITHER platform (keep_genes <- expressed_short | expressed_long) and drops only genes
# all-zero in BOTH -- so the trim destroyed exactly the category that rule depends on,
# and §1's concordance rebuilt as 27,916 genes instead of the published 29,185.
#
# The trim's "verified result-neutral" header was TRUE -- for the DE pipeline, where
# filterByExpr keeps an identical 25,955 genes / 162,800 isoforms either way. Neutrality
# was established on one path and generalised to all of them.
#
# Reverting reproduces the canonical render exactly: 29,447 common / 997 / 796 / 262
# all-zero / 29,185 final, Pearson 0.896 (0.834-0.906), Spearman 0.890 (0.849-0.901).
# Pete's decision, 2026-07-27: revert WHOLESALE. Every partial revert was measured to
# fail -- gene CSV alone 28,651, both CSVs alone still 28,651.
#
# The deposit is an UNPUBLISHED DRAFT Zenodo record, so files are replaced in place.
#
# WHAT THIS DOES NOT DO. The two DERIVED artifacts in the analysis repo are rebuilt from
# the deposit by their own Rmds and are NOT touched here:
#     tmp/out/dge_isoform_longread_2026.3.3.rds   (614,992 -> must become 645,272)
#     tmp/out/dge_shortread_gene_2026.3.2.rds     (46,571  -> must become 78,899)
# Regenerating them is a separate step, and skipping it makes this revert LOOK like it
# failed: correlation_analysis.Rmd reads the cached DGEList and still reports 28,651.
# =============================================================================
set -euo pipefail

DEP="$HOME/claude_projects/nmd_deposit_2026"
D="$DEP/source_data"
SRC="$HOME/claude_projects/nmd/sqanti/nmd_lungcells/results"
BUILD="$HOME/claude_projects/nmd/deposit"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="$DEP/.revert_${STAMP}"
mkdir -p "$LOG"

say() { printf '\n=== %s ===\n' "$1"; }

# ---------------------------------------------------------------------------
say "0. PRE-FLIGHT — refuse to start unless every source is present and full-sized"
# An exclusion is a claim: check each source explicitly rather than assuming the
# directory listing means they are the untrimmed vintage.
need_rows() { # $1=file $2=expected transcript/row count $3=label
  local n
  case "$3" in
    classification) n=$(( $(wc -l < "$1") - 1 )) ;;
    fasta)          n=$(grep -c '^>' "$1") ;;
    gxf)            n=$(grep -o 'transcript_id "[^"]*"' "$1" | sort -u | wc -l | tr -d ' ') ;;
  esac
  printf '  %-46s %s (want %s)\n' "$(basename "$1")" "$n" "$2"
  [ "$n" = "$2" ] || { echo "  ABORT: $1 is not the untrimmed vintage."; exit 1; }
}
for f in nmd_lungcells_classification.txt nmd_lungcells_corrected.fasta \
         nmd_lungcells_filtered.gtf nmd_lungcells_corrected.cds.gff3; do
  [ -f "$SRC/$f" ] || { echo "  ABORT: missing source $SRC/$f"; exit 1; }
done
need_rows "$SRC/nmd_lungcells_classification.txt" 645272 classification
need_rows "$SRC/nmd_lungcells_corrected.fasta"    645272 fasta

say "1. RECORD THE PRE-REVERT STATE (checksums, not copies)"
# The runbook's rule: freeze the substrate with a manifest of checksums. The trimmed
# state is fully reconstructible by re-running the trim, so copying gigabytes to back it
# up would buy nothing a hash does not.
( cd "$D" && find . -type f -print0 | sort -z | xargs -0 shasum -a 256 ) > "$LOG/pre_revert.sha256"
du -sh "$D" | tee "$LOG/pre_revert.size"
wc -l < "$LOG/pre_revert.sha256" | xargs printf '  %s files hashed\n'

say "2. RESTORE THE FOUR SQANTI FILES from the pre-trim SQANTI outputs"
for f in nmd_lungcells_classification.txt nmd_lungcells_corrected.fasta \
         nmd_lungcells_filtered.gtf nmd_lungcells_corrected.cds.gff3; do
  printf '  %s ... ' "$f"
  cp "$SRC/$f" "$D/sqanti/$f"
  printf 'done (%s)\n' "$(du -h "$D/sqanti/$f" | cut -f1)"
done

say "3. REGENERATE THE TWO COUNT MATRICES"
# Both builders already write the FULL matrix -- build_4ct_salmon_counts.R emits all
# 78,899 rows and self-verifies against the published DGEList. The trim ran AFTER them,
# so re-running them IS the revert. NOTE: both are hard-linked into zenodo_upload/, so
# rewriting the source_data copy updates what is staged for upload too.
Rscript "$BUILD/build_4ct_count_matrix.R"  2>&1 | tail -4
Rscript "$BUILD/build_4ct_salmon_counts.R" 2>&1 | tail -4

say "4. VERIFY — every file back to the full universe"
ok=1
chk() { # $1=label $2=actual $3=expected
  printf '  %-42s %-10s (want %s)  %s\n' "$1" "$2" "$3" \
    "$([ "$2" = "$3" ] && echo OK || { ok=0; echo MISMATCH; })"
}
chk "isoform counts rows" "$(( $(wc -l < "$D/nmd_lungcells_counts_4ct.csv") - 1 ))" 645272
chk "gene counts rows"    "$(( $(wc -l < "$D/salmon_gene_counts_4ct.csv") - 1 ))"   78899
chk "classification rows" "$(( $(wc -l < "$D/sqanti/nmd_lungcells_classification.txt") - 1 ))" 645272
chk "fasta records"       "$(grep -c '^>' "$D/sqanti/nmd_lungcells_corrected.fasta")" 645272
[ "$ok" = 1 ] || { echo "  ABORT: the revert did not restore the full universe."; exit 1; }

say "5. REFRESH CHECKSUMS"
( cd "$D" && find . -type f -print0 | sort -z | xargs -0 shasum -a 256 ) > "$DEP/MANIFEST.sha256"
cp "$DEP/MANIFEST.sha256" "$LOG/post_revert.sha256"
printf '  %s files\n' "$(wc -l < "$DEP/MANIFEST.sha256")"
du -sh "$D" | tee "$LOG/post_revert.size"

say "6. WHAT CHANGED"
diff <(awk '{print $2}' "$LOG/pre_revert.sha256") <(awk '{print $2}' "$LOG/post_revert.sha256") \
  > /dev/null && echo "  file LIST unchanged (contents differ where expected)" \
              || echo "  NOTE: the file list itself changed — inspect $LOG"
join -j 2 <(sort -k2 "$LOG/pre_revert.sha256") <(sort -k2 "$LOG/post_revert.sha256") \
  | awk '$2 != $3 {print "  CHANGED: " $1}'

cat <<NEXT

=== REVERT_DONE ===
Log: $LOG

STILL REQUIRED — this revert alone does NOT make correlation_analysis.Rmd print 29,185:
  1. rebuild zenodo_upload/  ->  bash $DEP/build_zenodo_upload.sh
     (sqanti.zip is a BUILT artifact, not a hard link, so it still holds trimmed files)
  2. regenerate the two derived DGELists in the analysis repo:
       tmp/out/dge_isoform_longread_2026.3.3.rds   (must become 645,272 x 26)
       tmp/out/dge_shortread_gene_2026.3.2.rds     (must become 78,899 x 26)
  3. GATE: correlation_analysis.Rmd must print
       29,447 common / 997 / 796 / 262 all-zero / 29,185 final
       Pearson 0.896 (0.834-0.906), Spearman 0.890 (0.849-0.901)
  4. Do NOT let trim_deposit_to_observed.sh run again on this deposit.
NEXT
