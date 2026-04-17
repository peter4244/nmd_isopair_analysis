#!/bin/bash
#
# stage_shortread.sh
#
# Stages the short-read (Illumina) SubSeries files for GEO submission.
#   1. Reads the two short-read pheno CSVs (FB_MV + DD_DO_AT) and matches
#      each sample to its paired-end FASTQs under the two raw directories.
#   2. Symlinks FASTQs into the staging dir with GEO-friendly names:
#        {sample_alias}_R1.fastq.gz / {sample_alias}_R2.fastq.gz
#   3. Copies the single combined Salmon gene-count matrix from nf-core/rnaseq
#      (already merged over both short-read batches) to staging.
#   4. Computes MD5 checksums for every staged file.
#   5. Emits samples_shortread.tsv, raw_files_shortread.tsv,
#      processed_files_shortread.tsv, and paired_end_shortread.tsv.
#
# Run on the cluster. No arguments.
#
# Dependencies: bash, coreutils, awk, md5sum (or md5 -q on macOS).
#
set -euo pipefail

# ============================================================
# CONFIGURATION
# ============================================================
PROJECT_REPO="/udd/repjc/claude_projects/nmd"                      # TODO: confirm repo checkout path on cluster
STAGING_ROOT="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/geo_submission"
STAGING_DIR="${STAGING_ROOT}/shortread"

# Pheno files (short-read)
PHENO_FB_MV="/udd/repjc/RESEARCH/NMD/NMD_Cell_Lines/pheno/FB_MV_pheno_2025.8.2.csv"
PHENO_DD_DO_AT="/udd/repjc/RESEARCH/NMD/NMD_Cell_Lines/pheno/DD_DO_AT_pheno_2025.8.15.csv"

# FASTQ raw directories (paired-end). The 20250718 batch stores each sample
# inside a per-sample subdirectory (PS_STC1505_*_ds.<hash>/…), so find_fastq()
# must recurse. The later batches may be flat or nested; the recursive search
# below handles either.
FASTQ_DIRS=(
  "/proj/regeps/regep00/studies/ExternalCellLines/data/rnaseq/Randell_Lung_Cells_2025/data/raw/20250718"
  "/proj/regeps/regep00/studies/ExternalCellLines/data/rnaseq/Randell_Lung_Cells_2025/data/raw/20250811"
  "/proj/regeps/regep00/studies/ExternalCellLines/data/rnaseq/Randell_Lung_Cells_2025/data/raw/20260129"
)

# Processed data (nf-core/rnaseq Salmon gene-level counts)
# Single combined nf-core run spanning both short-read batches (20250811 + 20260129).
COUNTS_MATRIX="/proj/regeps/regep00/studies/ExternalCellLines/data/rnaseq/Randell_Lung_Cells_2025/results/nfcore-rnaseq/20250811_plus_20260129/star_salmon/salmon.merged.gene_counts_length_scaled.tsv"

# TODO: confirm the short-read sequencing instrument (e.g., Illumina NovaSeq X, NovaSeq 6000)
INSTRUMENT_MODEL="Illumina NovaSeq 6000"     # TODO_CONFIRM
READ_LENGTH="150"                             # TODO_CONFIRM

# ============================================================
# SETUP
# ============================================================
mkdir -p "${STAGING_DIR}/raw" "${STAGING_DIR}/processed" "${STAGING_DIR}/metadata"

RAW_MD5="${STAGING_DIR}/md5sums_raw.txt"
PROCESSED_MD5="${STAGING_DIR}/md5sums_processed.txt"
SAMPLES_TSV="${STAGING_DIR}/metadata/samples_shortread.tsv"
RAW_TSV="${STAGING_DIR}/metadata/raw_files_shortread.tsv"
PROCESSED_TSV="${STAGING_DIR}/metadata/processed_files_shortread.tsv"
PAIRED_TSV="${STAGING_DIR}/metadata/paired_end_shortread.tsv"
SAMPLE_META_CSV="${STAGING_DIR}/processed/sample_metadata_shortread.csv"

: > "${RAW_MD5}"
: > "${PROCESSED_MD5}"

md5_of() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | awk '{print $1}'
  else
    md5 -q "$1"
  fi
}

# ============================================================
# STEP 1 — merge the two pheno files into a unified table and
#          derive GEO sample aliases.
# ============================================================
# Expected pheno columns (from the short-read DGE Rmd):
#   lib.size, norm.factors, id, sample, treatment, ct, sample_num, bam
# and a 'name' column derived in the Rmd as basename(bam) w/o '.markdup'.
#
# We'll work in a scratch TSV: alias<TAB>ct<TAB>donor<TAB>treatment<TAB>name
# where:
#   alias    = {ct}_{donor}_{treatment}_SR
#   name     = stem used to match FASTQ files (basename of BAM w/o .markdup...)
echo "[shortread] Building unified sample table from pheno files..."

SAMPLES_SCRATCH="${STAGING_DIR}/metadata/.samples_scratch.tsv"
printf "sample_alias\tct\tdonor\ttreatment\tname\tlib_size\tnorm_factors\n" > "${SAMPLES_SCRATCH}"

parse_pheno() {
  local pheno="$1"
  [[ ! -f "${pheno}" ]] && { echo "  WARN: pheno missing: ${pheno}" >&2; return 0; }
  # Use python if available for robust CSV parsing; fall back to awk.
  if command -v python3 >/dev/null 2>&1; then
    python3 - "${pheno}" <<'PYEOF' >> "${SAMPLES_SCRATCH}"
import csv, os, sys, re
p = sys.argv[1]
with open(p, newline='') as fh:
    reader = csv.DictReader(fh)
    # Try common column names
    for row in reader:
        ct   = (row.get('ct')   or '').strip()
        donor= (row.get('sample_id') or row.get('id') or row.get('donor') or '').strip()
        trt  = (row.get('treatment') or '').strip()
        bam  = (row.get('bam')  or '').strip()
        lib  = (row.get('lib.size') or '').strip()
        nf   = (row.get('norm.factors') or '').strip()
        if not (ct and donor and trt and bam):
            continue
        name = os.path.basename(bam)
        # Strip common suffixes used by nf-core STAR output: .markdup.sorted.bam etc.
        name = re.sub(r'\.markdup.*$', '', name)
        name = re.sub(r'\.aligned.*$',  '', name)
        name = re.sub(r'\.bam$',        '', name)
        # NOTE: do NOT rewrite dashes to underscores — FASTQ filenames reuse the
        # BAM stem verbatim (e.g., IL21522-001_N4UD-F10_L002), and sanitizing
        # dashes would break find_fastq() matching.
        alias = f"{ct}_{donor}_{trt}_SR"
        print(f"{alias}\t{ct}\t{donor}\t{trt}\t{name}\t{lib}\t{nf}")
PYEOF
  else
    # Minimal awk fallback (assumes no embedded commas in fields)
    awk -F, 'NR==1 {
               for (i=1;i<=NF;i++) h[$i]=i; next
             }
             {
               ct=$h["ct"]; d=$h["id"]; t=$h["treatment"]; b=$h["bam"];
               ls=$h["lib.size"]; nf=$h["norm.factors"];
               if (ct=="" || d=="" || t=="" || b=="") next
               n=b; sub(/.*\//,"",n); sub(/\.markdup.*/,"",n);
                    sub(/\.aligned.*/,"",n); sub(/\.bam$/,"",n); gsub(/-/,"_",n);
               printf "%s_%s_%s_SR\t%s\t%s\t%s\t%s\t%s\t%s\n", ct,d,t, ct,d,t,n,ls,nf
             }' "${pheno}" >> "${SAMPLES_SCRATCH}"
  fi
}

parse_pheno "${PHENO_FB_MV}"
parse_pheno "${PHENO_DD_DO_AT}"

NSAMP=$(( $(wc -l < "${SAMPLES_SCRATCH}") - 1 ))
echo "  ${NSAMP} samples parsed from pheno files"

# ============================================================
# STEP 2 — locate R1/R2 FASTQs for each sample and symlink
# ============================================================
echo ""
echo "[shortread] Locating FASTQs in raw directories..."

missing=0

# Helper: find a FASTQ for a given {name} and {read} in any of the raw dirs.
# Searches each directory *recursively* (handles per-sample nested subdirs
# like PS_STC1505_*_ds.<hash>/) and tries common Illumina/nf-core filename
# patterns, strict match first, then a looser substring match.
find_fastq() {
  local name="$1"
  local read="$2"    # 1 or 2
  local d f
  for d in "${FASTQ_DIRS[@]}"; do
    [[ -d "${d}" ]] || continue
    # Strict exact-name match (fastest, preferred).
    f=$(find "${d}" -type f \( \
          -name "${name}_R${read}_001.fastq.gz" -o \
          -name "${name}_R${read}.fastq.gz"     -o \
          -name "${name}_${read}.fastq.gz"      -o \
          -name "${name}_${read}.fq.gz" \
        \) -print 2>/dev/null | head -n 1)
    [[ -n "${f}" && -f "${f}" ]] && { echo "${f}"; return 0; }
    # Loose substring fallback.
    f=$(find "${d}" -type f \( \
          -name "*${name}*_R${read}*.fastq.gz" -o \
          -name "*${name}*_R${read}*.fq.gz" \
        \) -print 2>/dev/null | head -n 1)
    [[ -n "${f}" && -f "${f}" ]] && { echo "${f}"; return 0; }
  done
  return 1
}

# Rewrite the scratch table with R1/R2 paths resolved.
SAMPLES_RESOLVED="${STAGING_DIR}/metadata/.samples_resolved.tsv"
printf "sample_alias\tct\tdonor\ttreatment\tname\tlib_size\tnorm_factors\tr1\tr2\n" > "${SAMPLES_RESOLVED}"

tail -n +2 "${SAMPLES_SCRATCH}" | while IFS=$'\t' read -r alias ct donor trt name lib nf; do
  r1=$(find_fastq "${name}" 1 || true)
  r2=$(find_fastq "${name}" 2 || true)
  if [[ -z "${r1}" || -z "${r2}" ]]; then
    echo "  MISSING  ${alias}  (name=${name})  r1=${r1:-?}  r2=${r2:-?}" >&2
    missing=$((missing + 1))
    continue
  fi
  ln -snf "${r1}" "${STAGING_DIR}/raw/${alias}_R1.fastq.gz"
  ln -snf "${r2}" "${STAGING_DIR}/raw/${alias}_R2.fastq.gz"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${alias}" "${ct}" "${donor}" "${trt}" "${name}" "${lib}" "${nf}" "${r1}" "${r2}" >> "${SAMPLES_RESOLVED}"
  echo "  OK       ${alias}"
done

if (( missing > 0 )); then
  echo ""
  echo "ERROR: ${missing} sample(s) had missing FASTQs." >&2
  echo "Check the patterns in find_fastq() against the actual filenames in:" >&2
  for d in "${FASTQ_DIRS[@]}"; do echo "  ${d}" >&2; done
  exit 1
fi

# ============================================================
# STEP 3 — stage the (already-combined) Salmon count matrix
# ============================================================
echo ""
echo "[shortread] Staging Salmon gene-count matrix..."

merged_counts="${STAGING_DIR}/processed/salmon_gene_counts_length_scaled.tsv"
if [[ ! -f "${COUNTS_MATRIX}" ]]; then
  echo "  ERROR: counts matrix missing: ${COUNTS_MATRIX}" >&2
  exit 1
fi
cp -f "${COUNTS_MATRIX}" "${merged_counts}"
echo "  OK  $(basename "${merged_counts}")  ← ${COUNTS_MATRIX}"

# Stage the per-sample metadata CSV too
{
  printf "sample_alias,ct,donor,treatment,lib_size,norm_factors\n"
  tail -n +2 "${SAMPLES_RESOLVED}" | awk -F'\t' -v OFS=',' '{print $1,$2,$3,$4,$6,$7}'
} > "${SAMPLE_META_CSV}"
echo "  wrote ${SAMPLE_META_CSV}"

# ============================================================
# STEP 4 — MD5 checksums
# ============================================================
echo ""
echo "[shortread] Computing MD5 checksums..."

for f in "${STAGING_DIR}/raw/"*.fastq.gz; do
  [[ -e "${f}" ]] || continue
  real=$(readlink -f "${f}")
  printf "%s  %s\n" "$(md5_of "${real}")" "$(basename "${f}")" >> "${RAW_MD5}"
done
for f in "${STAGING_DIR}/processed/"*; do
  [[ -f "${f}" ]] || continue
  real=$(readlink -f "${f}")
  printf "%s  %s\n" "$(md5_of "${real}")" "$(basename "${f}")" >> "${PROCESSED_MD5}"
done

# ============================================================
# STEP 5 — GEO metadata TSVs (samples, raw files, processed files, paired-end)
# ============================================================
echo ""
echo "[shortread] Writing metadata TSVs..."

# Build a lookup: alias → md5 from RAW_MD5 for _R1 / _R2
declare -A MD5
while read -r m f; do MD5["${f}"]="${m}"; done < "${RAW_MD5}"

# --- samples_shortread.tsv ---
printf "Sample name\ttitle\tsource name\torganism\tcharacteristics: cell type\tcharacteristics: donor\tcharacteristics: treatment\tmolecule\tdescription\tprocessed data file\traw file read 1\traw file read 2\n" > "${SAMPLES_TSV}"

tail -n +2 "${SAMPLES_RESOLVED}" | while IFS=$'\t' read -r alias ct donor trt name lib nf r1 r2; do
  # Human-readable source name per cell type
  case "${ct}" in
    AT|AT2)  src="Alveolar Type 2 cells (AT2)" ;;
    DD)      src="Differentiated airway epithelial cells (submerged)" ;;
    DD_ALI)  src="Differentiated airway epithelial cells (air-liquid interface)" ;;
    DO)      src="Undifferentiated donor airway epithelial cells" ;;
    FB)      src="Lung fibroblasts" ;;
    MV)      src="Lung microvascular endothelial cells" ;;
    *)       src="${ct}" ;;
  esac
  title="${src} donor ${donor} ${trt} (short-read)"
  desc="Primary human ${src} treated with ${trt}; Illumina paired-end RNA-seq"
  printf "%s\t%s\t%s\tHomo sapiens\t%s\t%s\t%s\tpolyA RNA\t%s\tsalmon_gene_counts_length_scaled.tsv\t%s_R1.fastq.gz\t%s_R2.fastq.gz\n" \
    "${alias}" "${title}" "${src}" "${ct}" "${donor}" "${trt}" "${desc}" "${alias}" "${alias}" >> "${SAMPLES_TSV}"
done

# --- raw_files_shortread.tsv ---
printf "raw file\tfile type\tfile checksum\tinstrument model\tread length\tsingle or paired-end\n" > "${RAW_TSV}"
while read -r m f; do
  printf "%s\tfastq\t%s\t%s\t%s\tpaired-end\n" "${f}" "${m}" "${INSTRUMENT_MODEL}" "${READ_LENGTH}" >> "${RAW_TSV}"
done < "${RAW_MD5}"

# --- processed_files_shortread.tsv ---
printf "file name\tfile type\tfile checksum\n" > "${PROCESSED_TSV}"
while read -r m f; do
  case "${f}" in
    *.tsv) ftype="tab-delimited text" ;;
    *.csv) ftype="comma-delimited text" ;;
    *)     ftype="text" ;;
  esac
  printf "%s\t%s\t%s\n" "${f}" "${ftype}" "${m}" >> "${PROCESSED_TSV}"
done < "${PROCESSED_MD5}"

# --- paired_end_shortread.tsv ---
printf "file name 1\tfile name 2\n" > "${PAIRED_TSV}"
tail -n +2 "${SAMPLES_RESOLVED}" | awk -F'\t' -v OFS='\t' '{printf "%s_R1.fastq.gz\t%s_R2.fastq.gz\n",$1,$1}' >> "${PAIRED_TSV}"

echo "  samples:   ${SAMPLES_TSV}         ($(( $(wc -l < "${SAMPLES_TSV}") - 1 )) rows)"
echo "  raw:       ${RAW_TSV}             ($(( $(wc -l < "${RAW_TSV}") - 1 )) rows)"
echo "  processed: ${PROCESSED_TSV}       ($(( $(wc -l < "${PROCESSED_TSV}") - 1 )) rows)"
echo "  paired:    ${PAIRED_TSV}          ($(( $(wc -l < "${PAIRED_TSV}") - 1 )) rows)"

echo ""
echo "[shortread] Done."
echo "  Staging dir: ${STAGING_DIR}"
echo "  Paste each metadata TSV into the matching section of the short-read seq_template.xlsx."
