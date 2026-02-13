#!/bin/bash

# ============================================================================
# DYNAMIC SLURM JOB SUBMISSION FOR OARFISH QUANTIFICATION
# ============================================================================
#
# Description:
#   Monitors queue and intelligently submits jobs to maximize CPU usage
#   Uses oarfish for direct FASTQ->quantification without intermediate BAMs
#
# Features:
#   - Dynamic job submission with priority order (r8 and linux01 partitions)
#   - Job state monitoring and automatic resubmission of failed jobs
#   - Comprehensive status tracking across multiple files
#   - Atomic status file updates to prevent race conditions
#
# ============================================================================

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input/Output paths
BAM_LIST="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merged_collapsed_2025.12.21/scripts/bam_sample_map2.fofn"
OARFISH_INDEX="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merged_collapsed_2025.12.21/index/oarfish_gencode49_merged_collapsed_2025.12.21.mmi"
OUTPUT_DIR="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merged_collapsed_2025.12.21/oarfish_aligned_counts/oarfish_gencode49_merged_collapsed_2025.12.21"
LOG_DIR="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merged_collapsed_2025.12.21/logs/oarfish_logs/oarfish_gencode49_merged_collapsed_2025.12.21"

# Create output directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOG_DIR}

# Job submission parameters
MAX_RUNNING_JOBS=20                # Target number of concurrent running jobs
CPUS_PER_JOB=20                    # CPUs for oarfish
MEMORY=80                          # Memory in GB
TIME_LIMIT="7-00:00:00"            # 7 days

# Partitions to try (in order of preference)
PARTITIONS=("r8" "linux01")

# ============================================================================
# STATUS FILE INITIALIZATION
# ============================================================================

# Job tracking file - records all submission attempts
JOB_TRACKING_FILE="${LOG_DIR}/job_tracking.txt"
touch ${JOB_TRACKING_FILE}

# Status files
# 1. Submission status - job assignment and current state
STATUS_SUBMISSION="${LOG_DIR}/status_submission.txt"
# Format: SAMPLE_ID|STATUS|JOBID|CPUS|PARTITION|NODE|TIMESTAMP
if [ ! -f ${STATUS_SUBMISSION} ]; then
    if ! echo "# SAMPLE_ID|STATUS|JOBID|CPUS|PARTITION|NODE|TIMESTAMP" > ${STATUS_SUBMISSION} 2>/dev/null; then
        echo "ERROR: Cannot create status file ${STATUS_SUBMISSION}"
        echo "This may indicate a disk quota issue. Check: quota -s"
        exit 1
    fi
fi

# 2. Commands used - actual oarfish commands
STATUS_COMMANDS="${LOG_DIR}/status_commands.txt"
# Format: SAMPLE_ID|OARFISH_CMD
if [ ! -f ${STATUS_COMMANDS} ]; then
    echo "# SAMPLE_ID|OARFISH_CMD" > ${STATUS_COMMANDS}
fi

# 3. Monitoring log - timestamped status checks
STATUS_MONITORING="${LOG_DIR}/status_monitoring.txt"
# Format: TIMESTAMP|SAMPLE_ID|JOBID|STATUS|ELAPSED_TIME|NODE
if [ ! -f ${STATUS_MONITORING} ]; then
    echo "# TIMESTAMP|SAMPLE_ID|JOBID|STATUS|ELAPSED_TIME|NODE" > ${STATUS_MONITORING}
fi

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

# Atomically update status files with new entry
update_status_file() {
    local FILE=$1
    local HEADER=$2
    local SAMPLE_ID=$3
    local NEW_ENTRY=$4

    # Create temp file with header
    if ! echo "$HEADER" > ${FILE}.tmp 2>/dev/null; then
        echo "ERROR: Cannot write to ${FILE}.tmp - possible disk quota issue!" >&2
        echo "Attempting to write directly to file..." >&2
        # Try direct append as fallback (less safe but better than nothing)
        echo "$NEW_ENTRY" >> ${FILE} 2>/dev/null || {
            echo "CRITICAL: Cannot update status file at all! Check disk quota." >&2
            return 1
        }
        return 0
    fi

    # Copy all lines except the one for this sample (skip header)
    tail -n +2 ${FILE} 2>/dev/null | grep -v "^${SAMPLE_ID}|" >> ${FILE}.tmp || true

    # Add new entry
    echo "$NEW_ENTRY" >> ${FILE}.tmp

    # Replace original file atomically
    mv ${FILE}.tmp ${FILE} 2>/dev/null || {
        echo "ERROR: Cannot move temp file - possible disk quota issue!" >&2
        rm -f ${FILE}.tmp 2>/dev/null
        return 1
    }
}

# Record a job submission to tracking file
record_job() {
    local SAMPLE_ID=$1
    local JOBID=$2
    echo "${SAMPLE_ID}|${JOBID}|$(date +%s)" >> ${JOB_TRACKING_FILE}
}

# Update sample status with full job information
update_sample_status() {
    local SAMPLE_ID=$1
    local STATUS=$2
    local JOBID=$3
    local CPUS=$4
    local PARTITION=$5
    local NODE=${6:-"none"}
    local FASTQ_PATH=${7:-""}
    local OUTPUT_FILE=${8:-""}
    local LOG_FILE=${9:-""}

    # Update submission status file
    local ENTRY="${SAMPLE_ID}|${STATUS}|${JOBID}|${CPUS}|${PARTITION}|${NODE}|$(date '+%Y-%m-%d %H:%M:%S')"
    update_status_file "${STATUS_SUBMISSION}" "# SAMPLE_ID|STATUS|JOBID|CPUS|PARTITION|NODE|TIMESTAMP" "${SAMPLE_ID}" "${ENTRY}"

    # Update commands file if we have command info
    if [ -n "$FASTQ_PATH" ]; then
        # Build command that EXACTLY matches what runs in the SLURM job
        local OARFISH_CMD="oarfish -j ${CPUS} --seq-tech pac-bio-hifi --reads ${FASTQ_PATH} --index ${OARFISH_INDEX} --output ${OUTPUT_FILE} --filter-group no-filters --model-coverage 2> ${LOG_FILE}"

        # Update commands file
        local CMD_ENTRY="${SAMPLE_ID}|${OARFISH_CMD}"
        update_status_file "${STATUS_COMMANDS}" "# SAMPLE_ID|OARFISH_CMD" "${SAMPLE_ID}" "${CMD_ENTRY}"
    fi
}

# Mark a job as completed (preserving existing job info)
mark_completed() {
    local SAMPLE_ID=$1
    local SAMPLE_STATUS=$2

    local LAST_JOBID=$(echo "$SAMPLE_STATUS" | cut -d'|' -f3)
    local CPUS=$(echo "$SAMPLE_STATUS" | cut -d'|' -f4)
    local PARTITION=$(echo "$SAMPLE_STATUS" | cut -d'|' -f5)
    local NODE=$(echo "$SAMPLE_STATUS" | cut -d'|' -f6)

    local ENTRY="${SAMPLE_ID}|COMPLETED|${LAST_JOBID}|${CPUS}|${PARTITION}|${NODE}|$(date '+%Y-%m-%d %H:%M:%S')"
    update_status_file "${STATUS_SUBMISSION}" "# SAMPLE_ID|STATUS|JOBID|CPUS|PARTITION|NODE|TIMESTAMP" "${SAMPLE_ID}" "${ENTRY}"
}

# ============================================================================
# STATUS QUERY FUNCTIONS
# ============================================================================

# Get current status of a sample
get_sample_status() {
    local SAMPLE_ID=$1
    grep "^${SAMPLE_ID}|" ${STATUS_SUBMISSION} 2>/dev/null | tail -1
}

# Check if sample is already completed
is_completed() {
    local SAMPLE_ID=$1
    local STATUS=$(get_sample_status "$SAMPLE_ID")
    [[ "$STATUS" =~ \|COMPLETED\| ]]
}

# Get count of currently running jobs (only R state)
get_running_jobs() {
    local COUNT=0

    # Get all RUNNING samples from status file
    local RUNNING_SAMPLES=$(grep "|RUNNING|" ${STATUS_SUBMISSION} 2>/dev/null)

    while IFS='|' read -r SAMPLE_ID STATUS JOBID REST; do
        [[ -z "$JOBID" ]] && continue

        # Check if job is actually running (R state only)
        local STATE=$(squeue -j $JOBID -h -o "%t" 2>/dev/null)
        if [ "$STATE" = "R" ]; then
            ((COUNT++))
        fi
    done <<< "$RUNNING_SAMPLES"

    echo $COUNT
}

# Check if a sample needs submission
needs_submission() {
    local SAMPLE_ID=$1

    # Check if completed
    if is_completed "$SAMPLE_ID"; then
        echo "[$(date +%H:%M:%S)] ℹ ${SAMPLE_ID} already completed"
        return 1
    fi

    # Check if currently running (and actually in R state)
    local STATUS=$(get_sample_status "$SAMPLE_ID")
    if [[ "$STATUS" =~ \|RUNNING\| ]]; then
        local JOBID=$(echo "$STATUS" | cut -d'|' -f3)
        local STATE=$(squeue -j $JOBID -h -o "%t" 2>/dev/null)
        if [ "$STATE" = "R" ]; then
            echo "[$(date +%H:%M:%S)] ℹ ${SAMPLE_ID} already running (Job ${JOBID})"
            return 1
        else
            echo "[$(date +%H:%M:%S)] ⚠ ${SAMPLE_ID} marked RUNNING but job not in queue, will resubmit"
        fi
    fi

    # Check if pending
    if [[ "$STATUS" =~ \|PENDING\| ]]; then
        local JOBID=$(echo "$STATUS" | cut -d'|' -f3)
        local STATE=$(squeue -j $JOBID -h -o "%t" 2>/dev/null)
        if [ "$STATE" = "PD" ] || [ "$STATE" = "R" ]; then
            echo "[$(date +%H:%M:%S)] ℹ ${SAMPLE_ID} already pending (Job ${JOBID})"
            return 1
        else
            echo "[$(date +%H:%M:%S)] ⚠ ${SAMPLE_ID} marked PENDING but job not in queue, will resubmit"
        fi
    fi

    return 0
}

# ============================================================================
# JOB SUBMISSION FUNCTION
# ============================================================================

submit_job() {
    local SAMPLE_ID=$1
    local BASENAME=$2
    local FASTQ_PATH=$3
    local OUTPUT_FILE=$4
    local LOG_FILE=$5
    local CPUS=$6
    local PARTITION=$7

    # Create SLURM job script
    local JOB_SCRIPT=$(mktemp /tmp/oarfish_job.XXXXXX.sh)

    cat > ${JOB_SCRIPT} <<'SLURM_SCRIPT'
#!/bin/bash
#SBATCH --job-name=oarfish_SAMPLE_ID_PLACEHOLDER
#SBATCH --cpus-per-task=CPUS_PLACEHOLDER
#SBATCH --mem=MEMORY_PLACEHOLDER
#SBATCH --time=TIME_LIMIT_PLACEHOLDER
#SBATCH --partition=PARTITION_PLACEHOLDER
#SBATCH --output=LOG_DIR_PLACEHOLDER/SAMPLE_ID_PLACEHOLDER.slurm.out
#SBATCH --error=LOG_DIR_PLACEHOLDER/SAMPLE_ID_PLACEHOLDER.slurm.err

# Enable strict error handling
set -euo pipefail

echo "=========================================="
echo "OARFISH QUANTIFICATION"
echo "Sample: SAMPLE_ID_PLACEHOLDER"
echo "=========================================="
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo ""

# Initialize conda
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

if ! command -v oarfish &> /dev/null; then
    echo "ERROR: oarfish not found in conda environment"
    exit 1
fi
echo "Oarfish version: $(oarfish --version 2>&1 | head -1)"
echo ""

echo "=========================================="
echo "Configuration:"
echo "  Threads: CPUS_PLACEHOLDER"
echo "  Input FASTQ: FASTQ_PATH_PLACEHOLDER"
echo "  Index: INDEX_PLACEHOLDER"
echo "  Output: OUTPUT_FILE_PLACEHOLDER"
echo "  Log: LOG_FILE_PLACEHOLDER"
echo "=========================================="
echo ""

# Run oarfish
echo "Starting oarfish quantification..."
START_TIME=$(date +%s)

oarfish -j CPUS_PLACEHOLDER \
    --seq-tech pac-bio-hifi \
    --reads FASTQ_PATH_PLACEHOLDER \
    --index INDEX_PLACEHOLDER \
    --output OUTPUT_FILE_PLACEHOLDER \
    --filter-group no-filters \
    --model-coverage \
    2> LOG_FILE_PLACEHOLDER

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "=========================================="
echo "Oarfish completed in ${ELAPSED} seconds"
echo "=========================================="
echo ""

# Validate output
if [ ! -d "OUTPUT_FILE_PLACEHOLDER" ]; then
    echo "ERROR: Output directory was not created"
    exit 1
fi
echo "✓ Output directory created"

# Check for key output files
if [ -f "OUTPUT_FILE_PLACEHOLDER/quant.json" ]; then
    echo "✓ quant.json found"
else
    echo "WARNING: quant.json not found"
fi

echo ""
echo "=========================================="
echo "QUANTIFICATION COMPLETED SUCCESSFULLY"
echo "=========================================="
echo "Finished: $(date)"
echo "Output: OUTPUT_FILE_PLACEHOLDER"
echo "Log: LOG_FILE_PLACEHOLDER"
echo "=========================================="

exit 0
SLURM_SCRIPT

    # Replace placeholders
    sed -i "s|SAMPLE_ID_PLACEHOLDER|${SAMPLE_ID}|g" ${JOB_SCRIPT}
    sed -i "s|CPUS_PLACEHOLDER|${CPUS}|g" ${JOB_SCRIPT}
    sed -i "s|MEMORY_PLACEHOLDER|${MEMORY}G|g" ${JOB_SCRIPT}
    sed -i "s|TIME_LIMIT_PLACEHOLDER|${TIME_LIMIT}|g" ${JOB_SCRIPT}
    sed -i "s|PARTITION_PLACEHOLDER|${PARTITION}|g" ${JOB_SCRIPT}
    sed -i "s|LOG_DIR_PLACEHOLDER|${LOG_DIR}|g" ${JOB_SCRIPT}
    sed -i "s|FASTQ_PATH_PLACEHOLDER|${FASTQ_PATH}|g" ${JOB_SCRIPT}
    sed -i "s|INDEX_PLACEHOLDER|${OARFISH_INDEX}|g" ${JOB_SCRIPT}
    sed -i "s|OUTPUT_FILE_PLACEHOLDER|${OUTPUT_FILE}|g" ${JOB_SCRIPT}
    sed -i "s|LOG_FILE_PLACEHOLDER|${LOG_FILE}|g" ${JOB_SCRIPT}

    # Submit job
    local JOBID=$(sbatch --parsable ${JOB_SCRIPT} 2>&1)

    # Clean up temporary script
    rm -f ${JOB_SCRIPT}

    # Record submission
    record_job "$SAMPLE_ID" "$JOBID"

    echo "$JOBID"
}

# ============================================================================
# SAMPLE SUBMISSION LOGIC
# ============================================================================

try_submit_sample() {
    local SAMPLE_NUM=$1
    local BAM_PATH=$2

    # Convert .bam to .fastq.gz to get the input file
    FASTQ_PATH="${BAM_PATH%.bam}.fastq.gz"

    # Extract the base filename for the output
    BASENAME=$(basename "$BAM_PATH" .bam)

    # Create sample-numbered output names
    SAMPLE_ID="sample${SAMPLE_NUM}"
    OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE_ID}_${BASENAME}"
    LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_${BASENAME}.log"

    # Check if this sample needs submission
    if ! needs_submission "${SAMPLE_ID}"; then
        return 1  # Already handled
    fi

    # Try each partition in order
    for PARTITION in "${PARTITIONS[@]}"; do
        # Check current running jobs
        RUNNING=$(get_running_jobs)

        if [ $RUNNING -ge $MAX_RUNNING_JOBS ]; then
            echo "[$(date +%H:%M:%S)] At max running jobs (${RUNNING}/${MAX_RUNNING_JOBS}), waiting..."
            return 2  # Signal to wait
        fi

        # Submit the job
        echo "[$(date +%H:%M:%S)] Submitting ${SAMPLE_ID} (${BASENAME}) - ${CPUS_PER_JOB} CPUs on ${PARTITION}..."
        JOBID=$(submit_job "$SAMPLE_ID" "$BASENAME" "$FASTQ_PATH" "$OUTPUT_FILE" "$LOG_FILE" "$CPUS_PER_JOB" "$PARTITION")

        # Wait 10 seconds to see if job gets scheduled
        sleep 10

        # Check if job is running or pending
        JOB_STATE=$(squeue -j $JOBID -h -o "%t %R" 2>/dev/null)
        STATE=$(echo "$JOB_STATE" | awk '{print $1}')
        NODE=$(squeue -j $JOBID -h -o "%N" 2>/dev/null)

        # Handle case where job might have already finished or failed
        if [ -z "$STATE" ]; then
            echo "[$(date +%H:%M:%S)] Job ${JOBID} not found in queue (may have failed quickly), trying next partition..."
            sleep 2
            continue
        fi

        # If running, mark as RUNNING and return success
        if [ "$STATE" = "R" ]; then
            echo "[$(date +%H:%M:%S)] ✓ Job ${JOBID} is RUNNING on ${NODE}!"
            update_sample_status "$SAMPLE_ID" "RUNNING" "$JOBID" "$CPUS_PER_JOB" "$PARTITION" "$NODE" "$FASTQ_PATH" "$OUTPUT_FILE" "$LOG_FILE"

            # Log to monitoring file
            echo "$(date '+%Y-%m-%d %H:%M:%S')|${SAMPLE_ID}|${JOBID}|STARTED|00:00:00|${NODE}" >> ${STATUS_MONITORING}

            # Re-check if we've hit the limit after this job started
            RUNNING=$(get_running_jobs)
            if [ $RUNNING -ge $MAX_RUNNING_JOBS ]; then
                echo "[$(date +%H:%M:%S)] Reached max running jobs (${RUNNING}/${MAX_RUNNING_JOBS})"
                return 2  # Signal to wait before processing more samples
            fi
            return 0
        else
            # Pending or other state - mark as PENDING and continue
            echo "[$(date +%H:%M:%S)] Job ${JOBID} in state '${STATE}', leaving pending"
            update_sample_status "$SAMPLE_ID" "PENDING" "$JOBID" "$CPUS_PER_JOB" "$PARTITION" "none" "$FASTQ_PATH" "$OUTPUT_FILE" "$LOG_FILE"
            return 0
        fi
    done

    echo "[$(date +%H:%M:%S)] All partitions tried for ${SAMPLE_ID}"
    return 0
}

# ============================================================================
# MAIN EXECUTION LOOP
# ============================================================================

echo "=========================================="
echo "Dynamic SLURM Job Submission - OARFISH"
echo "Target: ${MAX_RUNNING_JOBS} running jobs"
echo "Resources: ${CPUS_PER_JOB} CPUs, ${MEMORY}G memory, ${TIME_LIMIT}"
echo "Partitions: ${PARTITIONS[@]}"
echo "=========================================="
echo ""

SAMPLE_NUM=1

while read -r BAM_PATH; do
    # Skip empty lines
    [[ -z "$BAM_PATH" ]] && continue

    # Try to submit this sample
    try_submit_sample $SAMPLE_NUM "$BAM_PATH"
    RESULT=$?

    # If at max jobs, wait before continuing
    WAIT_ITERATIONS=0
    while [ $RESULT -eq 2 ]; do
        RUNNING=$(get_running_jobs)
        # Only print waiting message every 5 minutes (10 iterations of 30s)
        if [ $((WAIT_ITERATIONS % 10)) -eq 0 ]; then
            echo "[$(date +%H:%M:%S)] Waiting for jobs to finish... (${RUNNING}/${MAX_RUNNING_JOBS} running)"
        fi
        sleep 30
        ((WAIT_ITERATIONS++))

        RUNNING=$(get_running_jobs)
        if [ $RUNNING -lt $MAX_RUNNING_JOBS ]; then
            try_submit_sample $SAMPLE_NUM "$BAM_PATH"
            RESULT=$?
        fi
    done

    # Increment sample number
    ((SAMPLE_NUM++))

    # Small delay between submissions
    sleep 2

done < ${BAM_LIST}

# ============================================================================
# SUBMISSION SUMMARY
# ============================================================================

echo ""
echo "=========================================="
echo "Submission Summary"
echo "=========================================="
echo ""

# Count samples by status
RUNNING=$(grep "|RUNNING|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)
COMPLETED=$(grep "|COMPLETED|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)
PENDING=$(grep "|PENDING|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)
FAILED=$(grep "|FAILED|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)

echo "Status breakdown:"
echo "  RUNNING: ${RUNNING}"
echo "  COMPLETED: ${COMPLETED}"
echo "  PENDING: ${PENDING}"
echo "  FAILED: ${FAILED}"
echo ""
echo "Sample status saved to: ${STATUS_SUBMISSION}"
echo "Commands saved to: ${STATUS_COMMANDS}"
echo "Monitoring log: ${STATUS_MONITORING}"
echo ""
echo "=========================================="
echo "Monitor jobs with:"
echo "  squeue -u \$USER"
echo ""
echo "Check job status:"
echo "  ./oarfish_complete_pipeline_dynamic_monitor.sh"
echo "=========================================="
