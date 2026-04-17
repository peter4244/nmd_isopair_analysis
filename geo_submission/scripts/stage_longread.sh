#!/bin/bash
#
# stage_longread.sh
#
# Stages the long-read (PacBio IsoSeq) SubSeries files for GEO submission.
#   1. Symlinks per-sample uBAMs from bams/ALL/ into the staging dir,
#      renamed to {sample_alias}.bam for GEO clarity.
#   2. Copies the long-read processed data file(s) (isocall count matrix
#      and isoform GTF) into the staging dir.
#   3. Computes MD5 checksums for every staged file.
#   4. Emits raw_files_longread.tsv and processed_files_longread.tsv
#      with the columns GEO expects (filename + md5sum + instrument, etc.).
#
# Run on the cluster. No arguments.
#
# Dependencies: bash, coreutils, md5sum (or on macOS, `md5 -q`).
#
set -euo pipefail

# ============================================================
# CONFIGURATION — edit paths here if they change
# ============================================================
PROJECT_REPO="/udd/repjc/claude_projects/nmd"                      # TODO: confirm repo checkout path on cluster
UBAM_SOURCE_DIR="/proj/regeps/regep00/studies/ExternalCellLines/data/longread/mrna/Randell_Lung_Cells_2025/data/raw/bams/ALL"
STAGING_ROOT="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/geo_submission"
STAGING_DIR="${STAGING_ROOT}/longread"

SAMPLE_META="${PROJECT_REPO}/geo_submission/sample_metadata_longread.csv"

# Processed data file(s) produced by the isocall pipeline.
# These should be the final count matrix + isoform GTF written to disk.
# TODO: confirm these paths point at the final isocall outputs you want to upload.
ISOCALL_RESULTS_DIR="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/isocall/nmd_lungcells"
ISOCALL_COUNT_MATRIX="${ISOCALL_RESULTS_DIR}/call/nmd_isocall.count_matrix.txt"
ISOCALL_ISOFORMS_GTF="${ISOCALL_RESULTS_DIR}/call/nmd_isocall.isoforms.gtf.gz"

# ============================================================
# SETUP
# ============================================================
mkdir -p "${STAGING_DIR}/raw" "${STAGING_DIR}/processed" "${STAGING_DIR}/metadata"

RAW_MD5="${STAGING_DIR}/md5sums_raw.txt"
PROCESSED_MD5="${STAGING_DIR}/md5sums_processed.txt"
RAW_TSV="${STAGING_DIR}/metadata/raw_files_longread.tsv"
PROCESSED_TSV="${STAGING_DIR}/metadata/processed_files_longread.tsv"

: > "${RAW_MD5}"
: > "${PROCESSED_MD5}"

# md5 helper: portable between Linux (md5sum) and macOS (md5 -q)
md5_of() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | awk '{print $1}'
  else
    md5 -q "$1"
  fi
}

# ============================================================
# STEP 1 — symlink uBAMs into staging/raw as {alias}.bam
# ============================================================
echo "[longread] Staging uBAMs → ${STAGING_DIR}/raw/"
missing=0

# Skip header, read CSV columns: 1=sample_alias, 14=source_bam_basename
# (column numbers are 1-based; see sample_metadata_longread.csv)
while IFS=, read -r alias title ct ctdesc donor trt platform instr strat src sel layout mol src_bam; do
  [[ "${alias}" == "sample_alias" ]] && continue
  src="${UBAM_SOURCE_DIR}/${src_bam}"
  dst="${STAGING_DIR}/raw/${alias}.bam"
  if [[ ! -f "${src}" ]]; then
    echo "  MISSING  ${alias}  (expected ${src})"
    missing=$((missing + 1))
    continue
  fi
  ln -snf "${src}" "${dst}"
  echo "  OK       ${alias}  ← $(basename "${src}")"
done < "${SAMPLE_META}"

if (( missing > 0 )); then
  echo ""
  echo "ERROR: ${missing} uBAM(s) missing from ${UBAM_SOURCE_DIR}." >&2
  echo "Verify the 'source_bam_basename' column of sample_metadata_longread.csv" >&2
  echo "and the contents of the bams/ALL directory, then rerun." >&2
  exit 1
fi

# ============================================================
# STEP 2 — copy processed files into staging/processed
# ============================================================
echo ""
echo "[longread] Staging processed data → ${STAGING_DIR}/processed/"

stage_processed_file() {
  local src="$1"
  local staged_name="$2"
  if [[ ! -f "${src}" ]]; then
    echo "  WARNING: processed source missing: ${src}" >&2
    return 1
  fi
  cp -f "${src}" "${STAGING_DIR}/processed/${staged_name}"
  echo "  OK  ${staged_name}  ← ${src}"
}

stage_processed_file "${ISOCALL_COUNT_MATRIX}" "isocall_transcript_counts.tsv"
stage_processed_file "${ISOCALL_ISOFORMS_GTF}"  "isocall_isoforms.gtf.gz"

# Also copy the sample metadata CSV so it travels with the submission.
cp -f "${SAMPLE_META}" "${STAGING_DIR}/processed/sample_metadata_longread.csv"
echo "  OK  sample_metadata_longread.csv  ← ${SAMPLE_META}"

# ============================================================
# STEP 3 — md5 checksums
# ============================================================
echo ""
echo "[longread] Computing MD5 checksums..."

for f in "${STAGING_DIR}/raw/"*.bam; do
  [[ -e "${f}" ]] || continue
  # Use -L / readlink so md5 reads the underlying file behind the symlink.
  real=$(readlink -f "${f}")
  printf "%s  %s\n" "$(md5_of "${real}")" "$(basename "${f}")" >> "${RAW_MD5}"
done

for f in "${STAGING_DIR}/processed/"*; do
  [[ -f "${f}" ]] || continue
  real=$(readlink -f "${f}")
  printf "%s  %s\n" "$(md5_of "${real}")" "$(basename "${f}")" >> "${PROCESSED_MD5}"
done

echo "  raw md5 → ${RAW_MD5}"
echo "  processed md5 → ${PROCESSED_MD5}"

# ============================================================
# STEP 4 — build GEO-format raw files TSV
# ============================================================
echo ""
echo "[longread] Writing ${RAW_TSV}"

printf "raw file\tfile type\tfile checksum\tinstrument model\tread length\tsingle or paired-end\n" > "${RAW_TSV}"

while read -r md5 fname; do
  printf "%s\tBAM\t%s\tPacBio Revio\tN/A (HiFi)\tsingle\n" "${fname}" "${md5}" >> "${RAW_TSV}"
done < "${RAW_MD5}"

echo "  rows written: $(( $(wc -l < "${RAW_TSV}") - 1 ))"

# ============================================================
# STEP 5 — build GEO-format processed files TSV
# ============================================================
echo ""
echo "[longread] Writing ${PROCESSED_TSV}"

printf "file name\tfile type\tfile checksum\n" > "${PROCESSED_TSV}"

while read -r md5 fname; do
  case "${fname}" in
    *.tsv)       ftype="tab-delimited text" ;;
    *.gtf.gz)    ftype="gzipped GTF annotation" ;;
    *.csv)       ftype="comma-delimited text" ;;
    *)           ftype="text" ;;
  esac
  printf "%s\t%s\t%s\n" "${fname}" "${ftype}" "${md5}" >> "${PROCESSED_TSV}"
done < "${PROCESSED_MD5}"

echo "  rows written: $(( $(wc -l < "${PROCESSED_TSV}") - 1 ))"

echo ""
echo "[longread] Done."
echo "  Staging dir: ${STAGING_DIR}"
echo "  Paste the contents of metadata/raw_files_longread.tsv into the RAW FILES section"
echo "  and metadata/processed_files_longread.tsv into the PROCESSED DATA FILES section"
echo "  of the long-read seq_template.xlsx."
