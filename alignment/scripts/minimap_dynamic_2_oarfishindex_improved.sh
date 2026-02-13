#!/bin/bash

# ============================================================================
# DYNAMIC SLURM JOB SUBMISSION FOR MINIMAP2 ALIGNMENT
# ============================================================================
#
# Description:
#   Monitors queue and intelligently submits jobs to maximize CPU usage
#   IMPROVED VERSION: Enhanced samtools sort robustness with error checking
#   and validation. Commands recorded in status files match actual execution.
#   Uses node-local /tmp with reduced memory settings (300M) for limited temp space.
#
# Features:
#   - Dynamic job submission with priority order (20CPU/10CPU, linux01/r8)
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
REFERENCE="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merged_collapsed_2025.12.21/index/minimap_gencode49_merged_collapsed_2025.12.21.mmi"
OUTPUT_DIR="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merged_collapsed_2025.12.21/bams/oarfish_gencode49_merged_collapsed_2025.12.21"
LOG_DIR="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merged_collapsed_2025.12.21/logs/oarfish_gencode49_merged_collapsed_2025.12.21"

# Create output directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOG_DIR}

# Job submission parameters
MAX_RUNNING_JOBS=20                # Target number of concurrent running jobs
BASE_MEMORY=40                     # Base memory for minimap2 (GB)
MEMORY_PER_SAMTOOLS_THREAD=2       # Memory per samtools thread (GB)

# IMPORTANT: Nodes have limited temp storage (~5GB in /tmp)
# Each job uses ~4-5GB temp space, so jobs must spread across different nodes
# The SLURM scheduler should naturally distribute jobs, but monitor carefully

# ============================================================================
# STATUS FILE INITIALIZATION
# ============================================================================

# Job tracking file - records all submission attempts
JOB_TRACKING_FILE="${LOG_DIR}/job_tracking.txt"
touch ${JOB_TRACKING_FILE}

# Status files - reorganized for clarity
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

# 2. Commands used - actual minimap2 and samtools commands
STATUS_COMMANDS="${LOG_DIR}/status_commands.txt"
# Format: SAMPLE_ID|MINIMAP2_CMD|SAMTOOLS_CMD
if [ ! -f ${STATUS_COMMANDS} ]; then
    echo "# SAMPLE_ID|MINIMAP2_CMD|SAMTOOLS_CMD" > ${STATUS_COMMANDS}
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

# Calculate thread allocation between minimap2 and samtools
# Returns: "MINIMAP2_THREADS SAMTOOLS_THREADS"
calculate_threads() {
    local CPUS=$1
    local SAMTOOLS_THREADS=6
    local MINIMAP2_THREADS=$((CPUS - 6))

    # For low-CPU jobs, adjust
    if [ $CPUS -lt 8 ]; then
        SAMTOOLS_THREADS=4
        MINIMAP2_THREADS=$((CPUS - 4))
    fi

    # Ensure minimum threads
    if [ $MINIMAP2_THREADS -lt 2 ]; then
        MINIMAP2_THREADS=2
        SAMTOOLS_THREADS=$((CPUS - 2))
    fi

    echo "${MINIMAP2_THREADS} ${SAMTOOLS_THREADS}"
}

# ============================================================================
# STATUS FILE MANAGEMENT FUNCTIONS
# ============================================================================

# Atomically update status files with new entry
# Consolidates the repeated pattern of updating status files
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
    local OUTPUT_BAM=${8:-""}
    local LOG_FILE=${9:-""}

    # Update submission status file
    local ENTRY="${SAMPLE_ID}|${STATUS}|${JOBID}|${CPUS}|${PARTITION}|${NODE}|$(date '+%Y-%m-%d %H:%M:%S')"
    update_status_file "${STATUS_SUBMISSION}" "# SAMPLE_ID|STATUS|JOBID|CPUS|PARTITION|NODE|TIMESTAMP" "${SAMPLE_ID}" "${ENTRY}"

    # Update commands file if we have command info
    if [ -n "$FASTQ_PATH" ]; then
        # Calculate thread allocation
        read MINIMAP2_THREADS SAMTOOLS_THREADS <<< $(calculate_threads $CPUS)

        # Build commands that EXACTLY match what runs in the SLURM job
        # Note: SLURM_TMPDIR is unique per job and automatically cleaned up
        local MINIMAP2_CMD="minimap2 -t ${MINIMAP2_THREADS} -ax map-hifi --split-prefix \${SLURM_TMPDIR}/tmp --eqx -N 100 ${REFERENCE} ${FASTQ_PATH} 2> ${LOG_FILE}"
        local SAMTOOLS_CMD="samtools sort -n -@ ${SAMTOOLS_THREADS} -m 300M -T \${SLURM_TMPDIR}/temp_sort -o ${OUTPUT_BAM}"

        # Update commands file
        local CMD_ENTRY="${SAMPLE_ID}|${MINIMAP2_CMD}|${SAMTOOLS_CMD}"
        update_status_file "${STATUS_COMMANDS}" "# SAMPLE_ID|MINIMAP2_CMD|SAMTOOLS_CMD" "${SAMPLE_ID}" "${CMD_ENTRY}"
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
    grep "^${SAMPLE_ID}|" ${STATUS_SUBMISSION} 2>/dev/null | tail -1 || echo ""
}

# Count current running jobs for this user
get_running_jobs() {
    squeue -u $USER -h -t R | wc -l
}

# Check if sample needs submission based on current status
# Returns: 0 if needs submission, 1 if already handled
needs_submission() {
    local SAMPLE_ID=$1
    local SAMPLE_STATUS=$(get_sample_status "${SAMPLE_ID}")

    # No status record means needs submission
    [ -z "$SAMPLE_STATUS" ] && return 0

    # Parse status fields
    local STATUS=$(echo "$SAMPLE_STATUS" | cut -d'|' -f2)
    local LAST_JOBID=$(echo "$SAMPLE_STATUS" | cut -d'|' -f3)

    case "$STATUS" in
        COMPLETED)
            echo "[$(date +%H:%M:%S)] ${SAMPLE_ID}: Already completed (Job ${LAST_JOBID}), skipping"
            return 1
            ;;
        RUNNING)
            # Verify job is actually still running
            local CURRENT_STATE=$(squeue -j ${LAST_JOBID} -h -o "%t" 2>/dev/null)
            if [ "$CURRENT_STATE" == "R" ]; then
                echo "[$(date +%H:%M:%S)] ${SAMPLE_ID}: Job ${LAST_JOBID} still running, skipping"
                return 1
            fi

            # Job not running - check final state
            local JOB_STATE=$(sacct -j ${LAST_JOBID} -n -o State -X 2>/dev/null | head -1 | tr -d ' ')
            if [ "$JOB_STATE" == "COMPLETED" ]; then
                mark_completed "${SAMPLE_ID}" "${SAMPLE_STATUS}"
                echo "[$(date +%H:%M:%S)] ${SAMPLE_ID}: Job ${LAST_JOBID} completed, skipping"
                return 1
            else
                echo "[$(date +%H:%M:%S)] ${SAMPLE_ID}: Job ${LAST_JOBID} was running but now ${JOB_STATE}, needs resubmission"
                return 0
            fi
            ;;
        PENDING|NOT_SUBMITTED|FAILED)
            echo "[$(date +%H:%M:%S)] ${SAMPLE_ID}: Status is ${STATUS}, attempting submission"
            return 0
            ;;
        *)
            echo "[$(date +%H:%M:%S)] ${SAMPLE_ID}: Unknown status ${STATUS}, skipping"
            return 1
            ;;
    esac
}

# ============================================================================
# JOB SUBMISSION FUNCTIONS
# ============================================================================

# Function to submit a job with given parameters
submit_job() {
    local SAMPLE_ID=$1
    local BASENAME=$2
    local FASTQ_PATH=$3
    local OUTPUT_BAM=$4
    local LOG_FILE=$5
    local CPUS=$6
    local PARTITION=$7
    
    # Create the SLURM job script
    JOB_SCRIPT=$(mktemp)
    cat > ${JOB_SCRIPT} <<'EOF'
#!/bin/bash
#SBATCH --job-name=minimap2_SAMPLE_ID_PLACEHOLDER
#SBATCH --output=LOG_DIR_PLACEHOLDER/SAMPLE_ID_PLACEHOLDER_BASENAME_PLACEHOLDER.slurm.out
#SBATCH --error=LOG_DIR_PLACEHOLDER/SAMPLE_ID_PLACEHOLDER_BASENAME_PLACEHOLDER.slurm.err
#SBATCH --mem=80G
#SBATCH --cpus-per-task=CPUS_PLACEHOLDER
#SBATCH --time=7-4:00:00
#SBATCH --partition=PARTITION_PLACEHOLDER

# Enable strict error handling
set -euo pipefail

echo "Starting minimap2 alignment for SAMPLE_ID_PLACEHOLDER"
echo "Input: FASTQ_PATH_PLACEHOLDER"
echo "Output: OUTPUT_BAM_PLACEHOLDER"
echo "Started at: $(date)"

# Initialize conda - try RHEL 8 first (for r8/stratus nodes), then RHEL 7
if [ -f /app/conda-24.3.0@i86-rhel8.0/etc/profile.d/conda.sh ]; then
    source /app/conda-24.3.0@i86-rhel8.0/etc/profile.d/conda.sh
elif [ -f /app/conda-24.3.0@i86-rhel7.0/etc/profile.d/conda.sh ]; then
    source /app/conda-24.3.0@i86-rhel7.0/etc/profile.d/conda.sh
else
    echo "ERROR: Could not find conda initialization script"
    exit 1
fi

# Activate conda environment
conda activate oarfish_12.2025

# Verify samtools is available
if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found in conda environment or modules"
    echo "Please install samtools: conda install -n oarfish_12.2025 samtools"
    exit 1
fi

# Set up temp directory - use SLURM_TMPDIR if available, otherwise use node-local /tmp
if [ -z "${SLURM_TMPDIR:-}" ]; then
    SLURM_TMPDIR="/tmp/slurm_job_${SLURM_JOB_ID}_${USER}"
    mkdir -p ${SLURM_TMPDIR}
    echo "WARNING: SLURM_TMPDIR not set by SLURM, using node-local /tmp: ${SLURM_TMPDIR}"
    echo "Consider adding '#SBATCH --tmp=50G' to job script for better temp space management"
else
    echo "Using SLURM temp directory: ${SLURM_TMPDIR}"
fi

echo "Using ${SLURM_CPUS_PER_TASK} CPUs for alignment"

# IMPROVED: Fixed samtools threads to control memory usage
# NOTE: Using 300M per thread due to limited node temp storage (~5GB /tmp)
SAMTOOLS_THREADS=6
MINIMAP2_THREADS=$((SLURM_CPUS_PER_TASK - 6))

# For lower CPU jobs, scale down
if [ $SLURM_CPUS_PER_TASK -lt 8 ]; then
    SAMTOOLS_THREADS=4
    MINIMAP2_THREADS=$((SLURM_CPUS_PER_TASK - 4))
fi

# Ensure minimum threads
if [ $MINIMAP2_THREADS -lt 2 ]; then
    MINIMAP2_THREADS=2
    SAMTOOLS_THREADS=$((SLURM_CPUS_PER_TASK - 2))
fi

echo "Thread allocation: minimap2=${MINIMAP2_THREADS}, samtools=${SAMTOOLS_THREADS} (fixed for memory control)"
echo "Using SLURM temp directory: ${SLURM_TMPDIR}"

# Run minimap2 and pipe to samtools with robust configuration
# Using SLURM_TMPDIR ensures unique temp files per job and automatic cleanup
# NOTE: -m reduced to 300M due to limited node-local temp storage (only ~5GB available)
minimap2 -t ${MINIMAP2_THREADS} -ax map-hifi --split-prefix ${SLURM_TMPDIR}/tmp --eqx -N 100 \
    REFERENCE_PLACEHOLDER \
    FASTQ_PATH_PLACEHOLDER \
    2> LOG_FILE_PLACEHOLDER | \
    samtools sort -n \
    -@ ${SAMTOOLS_THREADS} \
    -m 300M \
    -T ${SLURM_TMPDIR}/temp_sort \
    -o OUTPUT_BAM_PLACEHOLDER

# Verify the output BAM was created successfully
if [ ! -f OUTPUT_BAM_PLACEHOLDER ]; then
    echo "ERROR: Output BAM file was not created: OUTPUT_BAM_PLACEHOLDER"
    echo "Check for samtools sort failures or disk space issues"
    exit 1
fi

# Check if BAM file is valid
echo "Validating BAM file..."
if ! samtools quickcheck OUTPUT_BAM_PLACEHOLDER; then
    echo "ERROR: Output BAM file is corrupted: OUTPUT_BAM_PLACEHOLDER"
    exit 1
fi

# Get file size for verification
BAM_SIZE=$(du -h OUTPUT_BAM_PLACEHOLDER | cut -f1)
echo "Output BAM size: ${BAM_SIZE}"

# Index the BAM file
echo "Indexing BAM file..."
samtools index -@ ${SAMTOOLS_THREADS} OUTPUT_BAM_PLACEHOLDER

# Clean up custom temp directory if we created it (not managed by SLURM)
if [[ ${SLURM_TMPDIR} == /tmp/slurm_job_* ]]; then
    echo "Cleaning up temp directory: ${SLURM_TMPDIR}"
    rm -rf ${SLURM_TMPDIR}
fi

echo "Successfully completed alignment for SAMPLE_ID_PLACEHOLDER"
echo "Output BAM: OUTPUT_BAM_PLACEHOLDER"
echo "Finished at: $(date)"
echo "Exit status: $?"
EOF

    # Replace placeholders
    sed -i "s|SAMPLE_ID_PLACEHOLDER|${SAMPLE_ID}|g" ${JOB_SCRIPT}
    sed -i "s|BASENAME_PLACEHOLDER|${BASENAME}|g" ${JOB_SCRIPT}
    sed -i "s|FASTQ_PATH_PLACEHOLDER|${FASTQ_PATH}|g" ${JOB_SCRIPT}
    sed -i "s|OUTPUT_BAM_PLACEHOLDER|${OUTPUT_BAM}|g" ${JOB_SCRIPT}
    sed -i "s|LOG_FILE_PLACEHOLDER|${LOG_FILE}|g" ${JOB_SCRIPT}
    sed -i "s|REFERENCE_PLACEHOLDER|${REFERENCE}|g" ${JOB_SCRIPT}
    sed -i "s|LOG_DIR_PLACEHOLDER|${LOG_DIR}|g" ${JOB_SCRIPT}
    sed -i "s|CPUS_PLACEHOLDER|${CPUS}|g" ${JOB_SCRIPT}
    sed -i "s|PARTITION_PLACEHOLDER|${PARTITION}|g" ${JOB_SCRIPT}
    
    # Submit the job
    JOBID=$(sbatch --parsable ${JOB_SCRIPT})
    
    # Record this submission
    record_job "${SAMPLE_ID}" "${JOBID}"
    
    # Clean up temporary job script
    rm ${JOB_SCRIPT}
    
    # Echo to stderr so it doesn't interfere with return value
    echo "[$(date +%H:%M:%S)] Submitted ${SAMPLE_ID}: ${BASENAME} -> ${PARTITION} with ${CPUS} CPUs (Job ${JOBID})" >&2
    
    # Return only the job ID
    echo "$JOBID"
}

# Function to try submitting a job with priority order
try_submit_sample() {
    local SAMPLE_NUM=$1
    local BAM_PATH=$2
    
    # Convert .bam to .fastq.gz to get the input file
    FASTQ_PATH="${BAM_PATH%.bam}.fastq.gz"
    
    # Extract the base filename for the output
    BASENAME=$(basename "$BAM_PATH" .bam)
    
    # Create sample-numbered output names
    SAMPLE_ID="sample${SAMPLE_NUM}"
    OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE_ID}_${BASENAME}.bam"
    LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_${BASENAME}.log"
    
    # Check if this sample needs submission (based on job tracking)
    if ! needs_submission "${SAMPLE_ID}"; then
        return 1  # Already handled
    fi
    
    # Priority order: CPU first (20 > 10), then partition (linux01 > r8)
    # Try in this order:
    # 1. 20 CPUs on linux01
    # 2. 20 CPUs on r8
    # 3. 10 CPUs on linux01
    # 4. 10 CPUs on r8
    
    local STRATEGIES=("20:linux01" "20:r8" "10:linux01" "10:r8")
    
    for STRATEGY in "${STRATEGIES[@]}"; do
        CPUS="${STRATEGY%:*}"
        PARTITION="${STRATEGY#*:}"
        
        # Check current running jobs
        RUNNING=$(get_running_jobs)
        
        if [ $RUNNING -ge $MAX_RUNNING_JOBS ]; then
            echo "[$(date +%H:%M:%S)] At max running jobs (${RUNNING}/${MAX_RUNNING_JOBS}), waiting..."
            return 2  # Signal to wait
        fi
        
        # Submit the job
        echo "[$(date +%H:%M:%S)] Trying ${CPUS} CPUs on ${PARTITION}..."
        JOBID=$(submit_job "$SAMPLE_ID" "$BASENAME" "$FASTQ_PATH" "$OUTPUT_BAM" "$LOG_FILE" "$CPUS" "$PARTITION")
        
        # Wait 10 seconds to see if job gets scheduled
        sleep 10
        
        # Check if job is running (R state only)
        JOB_STATE=$(squeue -j $JOBID -h -o "%t %R" 2>/dev/null)
        STATE=$(echo "$JOB_STATE" | awk '{print $1}')
        NODE=$(squeue -j $JOBID -h -o "%N" 2>/dev/null)
        
        # Handle case where job might have already finished or failed
        if [ -z "$STATE" ]; then
            echo "[$(date +%H:%M:%S)] Job ${JOBID} not found in queue (may have failed quickly), trying next strategy..."
            sleep 2
            continue
        fi
        
        if [ "$STATE" == "R" ]; then
            echo "[$(date +%H:%M:%S)] ✓ Job ${JOBID} is RUNNING on ${NODE}!"
            update_sample_status "$SAMPLE_ID" "RUNNING" "$JOBID" "$CPUS" "$PARTITION" "$NODE" "$FASTQ_PATH" "$OUTPUT_BAM" "$LOG_FILE"
            
            # Log to monitoring file immediately
            echo "$(date '+%Y-%m-%d %H:%M:%S')|${SAMPLE_ID}|${JOBID}|STARTED|00:00:00|${NODE}" >> ${STATUS_MONITORING}
            
            # Re-check if we've hit the limit after this job started
            RUNNING=$(get_running_jobs)
            if [ $RUNNING -ge $MAX_RUNNING_JOBS ]; then
                echo "[$(date +%H:%M:%S)] Reached max running jobs (${RUNNING}/${MAX_RUNNING_JOBS})"
                return 2  # Signal to wait before processing more samples
            fi
            return 0
        else
            # Not running (still PD or other state), cancel and try next strategy
            echo "[$(date +%H:%M:%S)] Job ${JOBID} in state '${STATE}' (not running), cancelling and trying next strategy..."
            scancel $JOBID 2>/dev/null
            sleep 2
        fi
    done
    
    # If all strategies failed, resubmit with lowest CPU request and leave it pending
    # Better to eventually run with 10 CPUs than wait forever for 20
    echo "[$(date +%H:%M:%S)] All strategies tried without immediate start, submitting with 10 CPUs on linux01 and leaving pending"
    JOBID=$(submit_job "$SAMPLE_ID" "$BASENAME" "$FASTQ_PATH" "$OUTPUT_BAM" "$LOG_FILE" "10" "linux01")
    update_sample_status "$SAMPLE_ID" "PENDING" "$JOBID" "10" "linux01" "none" "$FASTQ_PATH" "$OUTPUT_BAM" "$LOG_FILE"
    return 0
}

# ============================================================================
# MAIN EXECUTION LOOP
# ============================================================================

echo "=========================================="
echo "Dynamic SLURM Job Submission (IMPROVED)"
echo "Target: ${MAX_RUNNING_JOBS} running jobs"
echo "Enhanced samtools sort robustness"
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
echo "Job tracking saved to: ${JOB_TRACKING_FILE}"
echo "Monitoring log saved to: ${STATUS_MONITORING}"
echo ""

# ============================================================================
# JOB MONITORING
# ============================================================================

# Monitor running jobs until completion
if [ $RUNNING -gt 0 ]; then
    echo "=========================================="
    echo "Monitoring Running Jobs"
    echo "=========================================="
    echo ""
    echo "Checking job status every 2 minutes..."
    echo "Press Ctrl+C to exit monitoring (jobs will continue running)"
    echo ""
    
    # Track last counts to only show when they change
    LAST_RUNNING=$RUNNING
    LAST_COMPLETED=$COMPLETED
    LAST_FAILED=$FAILED

    # Track failures in memory in case status file updates fail
    declare -A FAILED_JOBS
    declare -A COMPLETED_JOBS

    while true; do
        # Get all RUNNING samples
        RUNNING_SAMPLES=$(tail -n +2 ${STATUS_SUBMISSION} | grep "|RUNNING|" | cut -d'|' -f1)
        
        if [ -z "$RUNNING_SAMPLES" ]; then
            echo "[$(date +%H:%M:%S)] All jobs completed or failed!"
            break
        fi
        
        # Check each running sample
        CHANGES_DETECTED=false
        for SAMPLE_ID in $RUNNING_SAMPLES; do
            SAMPLE_STATUS=$(get_sample_status "${SAMPLE_ID}")
            JOBID=$(echo "$SAMPLE_STATUS" | cut -d'|' -f3)
            CPUS=$(echo "$SAMPLE_STATUS" | cut -d'|' -f4)
            PARTITION=$(echo "$SAMPLE_STATUS" | cut -d'|' -f5)
            NODE=$(echo "$SAMPLE_STATUS" | cut -d'|' -f6)
            
            # Check if still running and get elapsed time
            CURRENT_STATE=$(squeue -j ${JOBID} -h -o "%t" 2>/dev/null)
            ELAPSED_TIME=$(squeue -j ${JOBID} -h -o "%M" 2>/dev/null)
            
            # Log the status check
            echo "$(date '+%Y-%m-%d %H:%M:%S')|${SAMPLE_ID}|${JOBID}|${CURRENT_STATE:-FINISHED}|${ELAPSED_TIME:-N/A}|${NODE}" >> ${STATUS_MONITORING}
            
            if [ "$CURRENT_STATE" != "R" ]; then
                # Job finished, check final state
                FINAL_STATE=$(sacct -j ${JOBID} -n -o State -X 2>/dev/null | head -1 | tr -d ' ')
                
                case "$FINAL_STATE" in
                    COMPLETED)
                        echo "[$(date +%H:%M:%S)] ✓ ${SAMPLE_ID} (Job ${JOBID}) COMPLETED successfully"
                        CHANGES_DETECTED=true
                        COMPLETED_JOBS["${SAMPLE_ID}"]="${JOBID}"
                        mark_completed "${SAMPLE_ID}" "${SAMPLE_STATUS}"
                        ;;
                    FAILED|CANCELLED|TIMEOUT|NODE_FAIL|OUT_OF_MEMORY)
                        echo "[$(date +%H:%M:%S)] ✗ ${SAMPLE_ID} (Job ${JOBID}) FAILED with state: ${FINAL_STATE}"
                        CHANGES_DETECTED=true
                        FAILED_JOBS["${SAMPLE_ID}"]="${FINAL_STATE}"
                        ENTRY="${SAMPLE_ID}|FAILED|${JOBID}|${CPUS}|${PARTITION}|${NODE}|$(date '+%Y-%m-%d %H:%M:%S')"
                        if ! update_status_file "${STATUS_SUBMISSION}" "# SAMPLE_ID|STATUS|JOBID|CPUS|PARTITION|NODE|TIMESTAMP" "${SAMPLE_ID}" "${ENTRY}"; then
                            echo "WARNING: Status file update failed for ${SAMPLE_ID}, but tracking in memory" >&2
                        fi
                        ;;
                    *)
                        if [ -n "$FINAL_STATE" ]; then
                            echo "[$(date +%H:%M:%S)] ? ${SAMPLE_ID} (Job ${JOBID}) in unknown state: ${FINAL_STATE}"
                            CHANGES_DETECTED=true
                        fi
                        ;;
                esac
            fi
        done
        
        # Get current counts from status file
        FILE_RUNNING=$(grep "|RUNNING|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)
        FILE_COMPLETED=$(grep "|COMPLETED|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)
        FILE_FAILED=$(grep "|FAILED|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)

        # Add in-memory tracked jobs (in case file updates failed)
        MEMORY_COMPLETED=${#COMPLETED_JOBS[@]}
        MEMORY_FAILED=${#FAILED_JOBS[@]}

        # Calculate total counts
        COMPLETED=$((FILE_COMPLETED + MEMORY_COMPLETED))
        FAILED=$((FILE_FAILED + MEMORY_FAILED))
        RUNNING=$((FILE_RUNNING - MEMORY_COMPLETED - MEMORY_FAILED))

        # Also get actual running count from SLURM as ground truth
        SLURM_RUNNING=$(get_running_jobs)

        # Only show status if something changed
        if [ "$CHANGES_DETECTED" = true ] || [ $RUNNING -ne $LAST_RUNNING ] || [ $COMPLETED -ne $LAST_COMPLETED ] || [ $FAILED -ne $LAST_FAILED ]; then
            echo "[$(date +%H:%M:%S)] Status - Running: ${RUNNING} (SLURM: ${SLURM_RUNNING}), Completed: ${COMPLETED}, Failed: ${FAILED}"
            if [ $MEMORY_FAILED -gt 0 ] || [ $MEMORY_COMPLETED -gt 0 ]; then
                echo "  Note: ${MEMORY_FAILED} failures and ${MEMORY_COMPLETED} completions tracked in memory only (status file update issues)"
            fi
            LAST_RUNNING=$RUNNING
            LAST_COMPLETED=$COMPLETED
            LAST_FAILED=$FAILED
        fi
        
        # Wait 2 minutes before next check
        sleep 120
    done
    
    echo ""
    echo "=========================================="
    echo "Monitoring Complete"
    echo "=========================================="

    # Report any jobs tracked only in memory
    if [ ${#FAILED_JOBS[@]} -gt 0 ]; then
        echo ""
        echo "Jobs that FAILED (tracked in memory due to disk quota issues):"
        for sample in "${!FAILED_JOBS[@]}"; do
            echo "  - ${sample}: ${FAILED_JOBS[$sample]}"
        done
    fi

    if [ ${#COMPLETED_JOBS[@]} -gt 0 ]; then
        echo ""
        echo "Jobs that COMPLETED (tracked in memory due to disk quota issues):"
        for sample in "${!COMPLETED_JOBS[@]}"; do
            echo "  - ${sample}: Job ${COMPLETED_JOBS[$sample]}"
        done
    fi
fi

# Final summary
RUNNING=$(grep "|RUNNING|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)
COMPLETED=$(grep "|COMPLETED|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)
PENDING=$(grep "|PENDING|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)
FAILED=$(grep "|FAILED|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)

echo ""
echo "Final Status:"
echo "  RUNNING: ${RUNNING}"
echo "  COMPLETED: ${COMPLETED}"
echo "  PENDING: ${PENDING}"
echo "  FAILED: ${FAILED}"
echo ""
if [ $PENDING -gt 0 ] || [ $FAILED -gt 0 ]; then
    echo "To retry PENDING/FAILED samples, run this script again."
fi
echo "=========================================="
