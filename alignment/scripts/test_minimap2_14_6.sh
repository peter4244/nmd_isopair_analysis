#!/bin/bash
#SBATCH --job-name=test_mm2_14_6
#SBATCH --output=/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merd_collapse2025.12.21/test_14_6.slurm.out
#SBATCH --error=/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merd_collapse2025.12.21/test_14_6.slurm.err
#SBATCH --mem=80G
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00
#SBATCH --partition=linux12hr

# ============================================================================
# MINIMAP2 THREAD OPTIMIZATION TEST - Configuration 14/6 (CURRENT DEFAULT)
# ============================================================================
# Thread allocation: 14 minimap2 threads + 6 samtools threads
# Tests: Shared filesystem temp (avoiding SLURM temp management)
# ============================================================================

# Enable strict error handling
set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================

# Thread allocation for this test
MINIMAP2_THREADS=14
SAMTOOLS_THREADS=6

# File paths
REFERENCE="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merged_collapsed_2025.12.21/index/minimap_gencode49_merged_collapsed_2025.12.21.mmi"
INPUT_FASTQ="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merd_collapse2025.12.21/sample1_first_100k_reads.fastq.gz"
OUTPUT_BAM="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merd_collapse2025.12.21/test_alignment_14_6.bam"
LOG_FILE="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merd_collapse2025.12.21/test_alignment_14_6.log"

# Temp directory on shared filesystem (avoids SLURM temp management)
TEMP_DIR="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merd_collapse2025.12.21/tmp_test_${SLURM_JOB_ID:-$$}"

# ============================================================================
# STARTUP REPORT
# ============================================================================

echo "=========================================="
echo "MINIMAP2 ALIGNMENT TEST"
echo "Thread Configuration: ${MINIMAP2_THREADS} minimap2 / ${SAMTOOLS_THREADS} samtools"
echo "=========================================="
echo "Testing with shared filesystem temp directory"
echo "Avoiding SLURM temp management (node-local /tmp)"
echo "=========================================="
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-not_submitted_via_slurm}"
echo ""

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

echo "Setting up environment..."

# Initialize conda - try RHEL 8 first (for r8/stratus nodes), then RHEL 7
if [ -f /app/conda-24.3.0@i86-rhel8.0/etc/profile.d/conda.sh ]; then
    source /app/conda-24.3.0@i86-rhel8.0/etc/profile.d/conda.sh
    echo "Using RHEL 8 conda"
elif [ -f /app/conda-24.3.0@i86-rhel7.0/etc/profile.d/conda.sh ]; then
    source /app/conda-24.3.0@i86-rhel7.0/etc/profile.d/conda.sh
    echo "Using RHEL 7 conda"
else
    echo "ERROR: Could not find conda initialization script"
    exit 1
fi

# Activate conda environment
conda activate oarfish_12.2025
echo "Activated conda environment: oarfish_12.2025"

# Verify samtools is available
if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found in conda environment"
    echo "Please install samtools: conda install -n oarfish_12.2025 samtools"
    exit 1
fi
echo "Samtools version: $(samtools --version | head -1)"
echo ""

# ============================================================================
# TEMP DIRECTORY SETUP
# ============================================================================

echo "Setting up temp directory..."
mkdir -p ${TEMP_DIR}
echo "Temp directory: ${TEMP_DIR}"
echo "This uses shared filesystem, avoiding SLURM temp management"
echo ""

# ============================================================================
# CONFIGURATION REPORT
# ============================================================================

echo "=========================================="
echo "Configuration:"
echo "  Minimap2 threads: ${MINIMAP2_THREADS}"
echo "  Samtools threads: ${SAMTOOLS_THREADS}"
echo "  Samtools memory: 1500M per thread"
echo "  Temp directory: ${TEMP_DIR} (shared filesystem)"
echo "  Reference: ${REFERENCE}"
echo "  Input: ${INPUT_FASTQ}"
echo "  Output: ${OUTPUT_BAM}"
echo "=========================================="
echo ""

# ============================================================================
# ALIGNMENT WITH DETAILED TIMING
# ============================================================================

echo "Starting alignment..."
TOTAL_START=$(date +%s)

# ----------------------------------------------------------------------------
# Minimap2 Phase
# ----------------------------------------------------------------------------

echo "[$(date +%H:%M:%S)] Running minimap2 alignment..."
MINIMAP2_START=$(date +%s)

minimap2 -t ${MINIMAP2_THREADS} \
    -ax map-hifi \
    --split-prefix ${TEMP_DIR}/minimap_split \
    --eqx \
    -N 100 \
    ${REFERENCE} \
    ${INPUT_FASTQ} \
    > ${TEMP_DIR}/minimap2_output.sam \
    2> ${LOG_FILE}

MINIMAP2_END=$(date +%s)
MINIMAP2_TIME=$((MINIMAP2_END - MINIMAP2_START))
echo "[$(date +%H:%M:%S)] Minimap2 completed in ${MINIMAP2_TIME} seconds"

# ----------------------------------------------------------------------------
# Samtools Sort Phase
# ----------------------------------------------------------------------------

echo "[$(date +%H:%M:%S)] Running samtools sort..."
SAMTOOLS_START=$(date +%s)

samtools sort -n \
    -@ ${SAMTOOLS_THREADS} \
    -m 1500M \
    -T ${TEMP_DIR}/samtools_sort \
    -o ${OUTPUT_BAM} \
    ${TEMP_DIR}/minimap2_output.sam

SAMTOOLS_END=$(date +%s)
SAMTOOLS_TIME=$((SAMTOOLS_END - SAMTOOLS_START))
echo "[$(date +%H:%M:%S)] Samtools sort completed in ${SAMTOOLS_TIME} seconds"

# Cleanup intermediate SAM file
echo "[$(date +%H:%M:%S)] Cleaning up intermediate SAM file..."
rm -f ${TEMP_DIR}/minimap2_output.sam

# ============================================================================
# VALIDATION
# ============================================================================

echo ""
echo "=========================================="
echo "Validating output..."
echo "=========================================="

# Check BAM file exists
if [ ! -f ${OUTPUT_BAM} ]; then
    echo "ERROR: Output BAM file was not created: ${OUTPUT_BAM}"
    echo "Check for samtools sort failures or disk space issues"
    exit 1
fi
echo "✓ BAM file created"

# Validate BAM integrity
echo "Running samtools quickcheck..."
if ! samtools quickcheck ${OUTPUT_BAM}; then
    echo "ERROR: Output BAM file is corrupted: ${OUTPUT_BAM}"
    exit 1
fi
echo "✓ BAM file is valid"

# Report file size
BAM_SIZE=$(du -h ${OUTPUT_BAM} | cut -f1)
echo "✓ Output BAM size: ${BAM_SIZE}"

# Index the BAM file
VALIDATION_START=$(date +%s)
echo "Indexing BAM file..."
samtools index -@ ${SAMTOOLS_THREADS} ${OUTPUT_BAM}
echo "✓ BAM file indexed"

VALIDATION_END=$(date +%s)

# ============================================================================
# TEMP DIRECTORY CLEANUP
# ============================================================================

echo ""
echo "=========================================="
echo "Cleanup"
echo "=========================================="

# Report temp directory usage before cleanup
echo "Temp directory usage:"
du -sh ${TEMP_DIR}

# Cleanup temp directory
echo "Cleaning up temp directory..."
rm -rf ${TEMP_DIR}
echo "✓ Temp directory removed"

# ============================================================================
# TIMING SUMMARY
# ============================================================================

TOTAL_END=$(date +%s)
TOTAL_TIME=$((TOTAL_END - TOTAL_START))
VALIDATION_TIME=$((VALIDATION_END - VALIDATION_START))

echo ""
echo "=========================================="
echo "TIMING SUMMARY (${MINIMAP2_THREADS} minimap2 / ${SAMTOOLS_THREADS} samtools)"
echo "=========================================="
echo "Minimap2 alignment:   ${MINIMAP2_TIME} seconds"
echo "Samtools sort:        ${SAMTOOLS_TIME} seconds"
echo "Validation/indexing:  ${VALIDATION_TIME} seconds"
echo "Total:                ${TOTAL_TIME} seconds"
echo "=========================================="

# ============================================================================
# SUCCESS REPORT
# ============================================================================

echo ""
echo "=========================================="
echo "TEST COMPLETED SUCCESSFULLY"
echo "=========================================="
echo "Finished: $(date)"
echo "Output: ${OUTPUT_BAM}"
echo "Index: ${OUTPUT_BAM}.bai"
echo "Log: ${LOG_FILE}"
echo ""
echo "Note: No SLURM temp management issues!"
echo "Temp files used shared filesystem successfully"
echo "=========================================="

exit 0
