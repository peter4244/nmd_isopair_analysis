#!/bin/bash
#SBATCH --job-name=mm2_align
#SBATCH --mem=80G
#SBATCH --cpus-per-task=5
#SBATCH --time=7-00:00:00
#SBATCH --output=/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/code/logs/minimap2_genome/mm2_align_%A_%a.log
#SBATCH --partition=linux01
#SBATCH --array=1,2,3,4,5,6,10,11,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38

source $(conda info --base)/etc/profile.d/conda.sh
conda activate oarfish_12.2025

# ---- Configuration ----
BASEDIR="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025"
INDEX="${BASEDIR}/results/reference_files/minimap_genome_index/GRCh38_splice.mmi"
JUNCTIONS="${BASEDIR}/results/reference_files/minimap_genome_index/known_junctions.bed"
OUTPUT_DIR="${BASEDIR}/data/bams/genome_hg38"
SAMPLE_MAP="${BASEDIR}/data/fastq_sample_map3.tsv"
THREADS=5

mkdir -p "${OUTPUT_DIR}"

# ---- Get input file and sample name for this array task ----
# fastq_sample_map3.tsv: two-column TSV (fastq_path<TAB>sample_name)
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_MAP}")
INPUT_FASTQ=$(echo "${LINE}" | cut -f1)
SAMPLE=$(echo "${LINE}" | cut -f2)

if [ -z "${INPUT_FASTQ}" ] || [ ! -f "${INPUT_FASTQ}" ]; then
    echo "ERROR: Input file not found for task ${SLURM_ARRAY_TASK_ID}: ${INPUT_FASTQ}"
    exit 1
fi

if [ -z "${SAMPLE}" ]; then
    echo "ERROR: Sample name not found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# ---- Skip if already aligned ----
if [ -f "${OUTPUT_DIR}/${SAMPLE}.aligned.bam.bai" ]; then
    echo "Already completed: ${SAMPLE} — skipping."
    exit 0
fi

echo "============================================"
echo "Sample: ${SAMPLE}"
echo "Input:  ${INPUT_FASTQ}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Start time: $(date)"
echo "============================================"

# ---- Align with minimap2 ----
# -x splice:hq    — splice-aware, high-fidelity (PacBio HiFi/Kinnex)
# --junc-bed       — known splice junctions
# -a               — output SAM
# --MD             — include MD tag for downstream use
# --cs=long        - include cs long tag
# -uf              — find forward-strand GTF junctions for cDNA

minimap2 \
    -ax splice:hq \
    --junc-bed "${JUNCTIONS}" \
    --MD \
    --cs=long \
    -uf \
    -t $((THREADS - 1)) \
    "${INDEX}" \
    "${INPUT_FASTQ}" \
| samtools sort -@ 0 -m 2G -o "${OUTPUT_DIR}/${SAMPLE}.aligned.bam"

# ---- Index the BAM ----
samtools index -@ ${THREADS} "${OUTPUT_DIR}/${SAMPLE}.aligned.bam"

echo "Finished: ${SAMPLE}"
echo "Output: ${OUTPUT_DIR}/${SAMPLE}.aligned.bam"
echo "End time: $(date)"
