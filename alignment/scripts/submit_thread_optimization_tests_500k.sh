#!/bin/bash

# ============================================================================
# MASTER SCRIPT: Submit Minimap2 Thread Optimization Tests
# ============================================================================
# This script submits 5 different thread configurations to test optimal
# CPU allocation between minimap2 and samtools.
# ============================================================================

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================

# Base configuration paths
REFERENCE="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merged_collapsed_2025.12.21/index/minimap_gencode49_merged_collapsed_2025.12.21.mmi"
INPUT_FASTQ="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merged_collapsed_2025.12.21/sample1_first_500k_reads.fastq.gz"
OUTPUT_DIR="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merged_collapsed_2025.12.21"

# Thread configurations to test: "minimap2_threads:samtools_threads"
CONFIGS=(
    "18:2"
    "16:4"
    "14:6"
    "12:8"
    "10:10"
)

echo "=========================================="
echo "MINIMAP2 THREAD OPTIMIZATION TESTS"
echo "=========================================="
echo "Submitting 5 thread configuration tests"
echo "Partition: r8"
echo "Test file: $(basename ${INPUT_FASTQ})"
echo ""

# ============================================================================
# FUNCTION: Create and submit test job
# ============================================================================

submit_test_job() {
    local MM2_THREADS=$1
    local SAM_THREADS=$2
    local CONFIG_NAME="${MM2_THREADS}_${SAM_THREADS}"

    # Create temporary job script
    local JOB_SCRIPT=$(mktemp)

    cat > ${JOB_SCRIPT} <<EOF
#!/bin/bash
#SBATCH --job-name=test_mm2_${CONFIG_NAME}
#SBATCH --output=${OUTPUT_DIR}/test_${CONFIG_NAME}.slurm.out
#SBATCH --error=${OUTPUT_DIR}/test_${CONFIG_NAME}.slurm.err
#SBATCH --mem=80G
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00
#SBATCH --partition=r8

# ============================================================================
# MINIMAP2 THREAD OPTIMIZATION TEST - Configuration ${CONFIG_NAME}
# ============================================================================
# Thread allocation: ${MM2_THREADS} minimap2 threads + ${SAM_THREADS} samtools threads
# Tests: Shared filesystem temp (avoiding SLURM temp management)
# ============================================================================

# Enable strict error handling
set -euo pipefail

# Configuration
MINIMAP2_THREADS=${MM2_THREADS}
SAMTOOLS_THREADS=${SAM_THREADS}
REFERENCE="${REFERENCE}"
INPUT_FASTQ="${INPUT_FASTQ}"
OUTPUT_BAM="${OUTPUT_DIR}/test_alignment_${CONFIG_NAME}.bam"
LOG_FILE="${OUTPUT_DIR}/test_alignment_${CONFIG_NAME}.log"
TEMP_DIR="${OUTPUT_DIR}/tmp_test_\${SLURM_JOB_ID}"

# Trap to print diagnostic info on error
trap 'echo ""; echo "========================================"; echo "ERROR at line \${LINENO}"; echo "Command: \${BASH_COMMAND}"; echo "========================================"; echo ""' ERR

# ============================================================================
# STARTUP REPORT
# ============================================================================

echo "=========================================="
echo "MINIMAP2 ALIGNMENT TEST"
echo "Thread Configuration: \${MINIMAP2_THREADS} minimap2 / \${SAMTOOLS_THREADS} samtools"
echo "=========================================="
echo "Testing with shared filesystem temp directory"
echo "Avoiding SLURM temp management (node-local /tmp)"
echo "=========================================="
echo "Started: \$(date)"
echo "Node: \$(hostname)"
echo "Job ID: \${SLURM_JOB_ID}"
echo ""

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

echo "Setting up environment..."

# Initialize conda - try RHEL 8 first, then RHEL 7
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

conda activate oarfish_12.2025
echo "Activated conda environment: oarfish_12.2025"

if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found in conda environment"
    exit 1
fi
echo "Samtools version: \$(samtools --version | head -1)"
echo ""

# ============================================================================
# TEMP DIRECTORY SETUP
# ============================================================================

echo "Setting up temp directory..."
mkdir -p \${TEMP_DIR}
echo "Temp directory: \${TEMP_DIR}"
echo "This uses shared filesystem, avoiding SLURM temp management"
echo ""

# ============================================================================
# CONFIGURATION REPORT
# ============================================================================

echo "=========================================="
echo "Configuration:"
echo "  Minimap2 threads: \${MINIMAP2_THREADS}"
echo "  Samtools threads: \${SAMTOOLS_THREADS}"
echo "  Samtools memory: 1500M per thread"
echo "  Temp directory: \${TEMP_DIR} (shared filesystem)"
echo "  Reference: \${REFERENCE}"
echo "  Input: \${INPUT_FASTQ}"
echo "  Output: \${OUTPUT_BAM}"
echo "=========================================="
echo ""

# ============================================================================
# ALIGNMENT WITH DETAILED TIMING
# ============================================================================

echo "Starting alignment..."
TOTAL_START=\$(date +%s)

# Minimap2 Phase
echo "[\$(date +%H:%M:%S)] Running minimap2 alignment..."
MINIMAP2_START=\$(date +%s)

minimap2 -t \${MINIMAP2_THREADS} \\
    -ax map-hifi \\
    --split-prefix \${TEMP_DIR}/minimap_split \\
    --eqx \\
    -N 100 \\
    \${REFERENCE} \\
    \${INPUT_FASTQ} \\
    > \${TEMP_DIR}/minimap2_output.sam \\
    2> \${LOG_FILE}

MINIMAP2_END=\$(date +%s)
MINIMAP2_TIME=\$((MINIMAP2_END - MINIMAP2_START))
echo "[\$(date +%H:%M:%S)] Minimap2 completed in \${MINIMAP2_TIME} seconds"

# Samtools Sort Phase
echo "[\$(date +%H:%M:%S)] Running samtools sort..."
SAMTOOLS_START=\$(date +%s)

samtools sort -n \\
    -@ \${SAMTOOLS_THREADS} \\
    -m 1500M \\
    -T \${TEMP_DIR}/samtools_sort \\
    -o \${OUTPUT_BAM} \\
    \${TEMP_DIR}/minimap2_output.sam

SAMTOOLS_END=\$(date +%s)
SAMTOOLS_TIME=\$((SAMTOOLS_END - SAMTOOLS_START))
echo "[\$(date +%H:%M:%S)] Samtools sort completed in \${SAMTOOLS_TIME} seconds"

# Cleanup intermediate SAM file
echo "[\$(date +%H:%M:%S)] Cleaning up intermediate SAM file..."
rm -f \${TEMP_DIR}/minimap2_output.sam

# ============================================================================
# VALIDATION
# ============================================================================

echo ""
echo "=========================================="
echo "Validating output..."
echo "=========================================="

if [ ! -f \${OUTPUT_BAM} ]; then
    echo "ERROR: Output BAM file was not created: \${OUTPUT_BAM}"
    exit 1
fi
echo "✓ BAM file created"

echo "Running samtools quickcheck..."
if ! samtools quickcheck \${OUTPUT_BAM}; then
    echo "ERROR: Output BAM file is corrupted: \${OUTPUT_BAM}"
    exit 1
fi
echo "✓ BAM file is valid"

BAM_SIZE=\$(du -h \${OUTPUT_BAM} | cut -f1)
echo "✓ Output BAM size: \${BAM_SIZE}"

VALIDATION_START=\$(date +%s)
echo "Indexing BAM file..."
if samtools index -@ \${SAMTOOLS_THREADS} \${OUTPUT_BAM} 2>&1; then
    echo "✓ BAM file indexed"
else
    echo "ERROR: BAM indexing failed (exit code: \$?)"
    echo "Continuing anyway to complete timing analysis..."
fi

VALIDATION_END=\$(date +%s)

# ============================================================================
# TEMP DIRECTORY CLEANUP
# ============================================================================

echo ""
echo "=========================================="
echo "Cleanup"
echo "=========================================="

echo "Temp directory usage:"
du -sh \${TEMP_DIR}

echo "Cleaning up temp directory..."
if rm -rf \${TEMP_DIR} 2>&1; then
    echo "✓ Temp directory removed"
else
    echo "WARNING: Failed to remove temp directory (may require manual cleanup)"
fi

# ============================================================================
# TIMING SUMMARY
# ============================================================================

TOTAL_END=\$(date +%s)
TOTAL_TIME=\$((TOTAL_END - TOTAL_START))
VALIDATION_TIME=\$((VALIDATION_END - VALIDATION_START))

echo ""
echo "=========================================="
echo "TIMING SUMMARY (\${MINIMAP2_THREADS} minimap2 / \${SAMTOOLS_THREADS} samtools)"
echo "=========================================="
echo "Minimap2 alignment:   \${MINIMAP2_TIME} seconds"
echo "Samtools sort:        \${SAMTOOLS_TIME} seconds"
echo "Validation/indexing:  \${VALIDATION_TIME} seconds"
echo "Total:                \${TOTAL_TIME} seconds"
echo "=========================================="

# ============================================================================
# SUCCESS REPORT
# ============================================================================

echo ""
echo "=========================================="
echo "TEST COMPLETED SUCCESSFULLY"
echo "=========================================="
echo "Finished: \$(date)"
echo "Output: \${OUTPUT_BAM}"
echo "Index: \${OUTPUT_BAM}.bai"
echo "Log: \${LOG_FILE}"
echo ""
echo "Note: No SLURM temp management issues!"
echo "Temp files used shared filesystem successfully"
echo "=========================================="

exit 0
EOF

    # Submit the job
    local JOBID=$(sbatch --parsable ${JOB_SCRIPT})

    # Clean up temp script
    rm ${JOB_SCRIPT}

    # Return job ID
    echo "${JOBID}"
}

# ============================================================================
# SUBMIT ALL TESTS
# ============================================================================

declare -a JOB_IDS

echo "Submitting tests..."
echo ""

for CONFIG in "${CONFIGS[@]}"; do
    MM2_THREADS=${CONFIG%:*}
    SAM_THREADS=${CONFIG#*:}

    echo -n "  Submitting ${MM2_THREADS} minimap2 / ${SAM_THREADS} samtools... "
    JOBID=$(submit_test_job ${MM2_THREADS} ${SAM_THREADS})
    JOB_IDS+=("${JOBID}")
    echo "Job ${JOBID}"
done

echo ""
echo "=========================================="
echo "All tests submitted successfully!"
echo "=========================================="
echo ""
echo "Job IDs: ${JOB_IDS[@]}"
echo ""
echo "Monitor jobs with:"
echo "  squeue -u \$USER"
echo ""
echo "Check status with:"
echo "  sacct -j $(IFS=,; echo "${JOB_IDS[*]}") --format=JobID,JobName,State,Elapsed,MaxRSS"
echo ""
echo "Output files will be in:"
echo "  ${OUTPUT_DIR}/test_*.slurm.out"
echo ""
echo "Compare results when done:"
echo "  grep -A 5 'TIMING SUMMARY' ${OUTPUT_DIR}/test_*.slurm.out"
echo ""
echo "=========================================="
