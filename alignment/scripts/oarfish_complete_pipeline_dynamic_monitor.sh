#!/bin/bash

# ============================================================================
# OARFISH JOB MONITORING SCRIPT
# ============================================================================
#
# Monitors running oarfish jobs and updates status files
#
# ============================================================================

# Paths (must match submission script)
LOG_DIR="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merged_collapsed_2025.12.21/logs/oarfish_logs/oarfish_gencode49_merged_collapsed_2025.12.21"
OUTPUT_DIR="/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/oarfish/gencode49_merged_collapsed_2025.12.21/oarfish_aligned_counts/oarfish_gencode49_merged_collapsed_2025.12.21"

STATUS_SUBMISSION="${LOG_DIR}/status_submission.txt"
STATUS_MONITORING="${LOG_DIR}/status_monitoring.txt"

# ============================================================================
# FUNCTIONS
# ============================================================================

update_status_file() {
    local FILE=$1
    local HEADER=$2
    local SAMPLE_ID=$3
    local NEW_ENTRY=$4

    echo "$HEADER" > ${FILE}.tmp
    tail -n +2 ${FILE} 2>/dev/null | grep -v "^${SAMPLE_ID}|" >> ${FILE}.tmp || true
    echo "$NEW_ENTRY" >> ${FILE}.tmp
    mv ${FILE}.tmp ${FILE}
}

# ============================================================================
# MAIN MONITORING LOOP
# ============================================================================

echo "=========================================="
echo "OARFISH JOB MONITOR"
echo "=========================================="
echo "Time: $(date)"
echo ""

# Get all samples that are marked RUNNING or PENDING
ACTIVE_JOBS=$(grep -E "\|(RUNNING|PENDING)\|" ${STATUS_SUBMISSION} 2>/dev/null)

if [ -z "$ACTIVE_JOBS" ]; then
    echo "No active jobs found in status file"
    echo ""
fi

RUNNING_COUNT=0
PENDING_COUNT=0
COMPLETED_COUNT=0
FAILED_COUNT=0

while IFS='|' read -r SAMPLE_ID STATUS JOBID CPUS PARTITION NODE TIMESTAMP; do
    [[ -z "$JOBID" ]] && continue

    # Get current job state from SLURM
    JOB_INFO=$(sacct -j ${JOBID} --format=JobID,State,Elapsed,NodeList --noheader --parsable2 2>/dev/null | head -1)

    if [ -z "$JOB_INFO" ]; then
        # Job not found in sacct - might be very recent or very old
        QUEUE_STATE=$(squeue -j ${JOBID} -h -o "%t" 2>/dev/null)
        if [ -z "$QUEUE_STATE" ]; then
            echo "[${SAMPLE_ID}] Job ${JOBID}: NOT FOUND (may have completed or been deleted)"
            # Check if output exists
            OUTPUT_PATH="${OUTPUT_DIR}/${SAMPLE_ID}_*"
            if ls ${OUTPUT_PATH}/quant.json 1> /dev/null 2>&1; then
                echo "  → Output found, marking as COMPLETED"
                ENTRY="${SAMPLE_ID}|COMPLETED|${JOBID}|${CPUS}|${PARTITION}|${NODE}|$(date '+%Y-%m-%d %H:%M:%S')"
                update_status_file "${STATUS_SUBMISSION}" "# SAMPLE_ID|STATUS|JOBID|CPUS|PARTITION|NODE|TIMESTAMP" "${SAMPLE_ID}" "${ENTRY}"
                echo "$(date '+%Y-%m-%d %H:%M:%S')|${SAMPLE_ID}|${JOBID}|COMPLETED|unknown|${NODE}" >> ${STATUS_MONITORING}
                ((COMPLETED_COUNT++))
            else
                echo "  → No output found, marking as FAILED"
                ENTRY="${SAMPLE_ID}|FAILED|${JOBID}|${CPUS}|${PARTITION}|${NODE}|$(date '+%Y-%m-%d %H:%M:%S')"
                update_status_file "${STATUS_SUBMISSION}" "# SAMPLE_ID|STATUS|JOBID|CPUS|PARTITION|NODE|TIMESTAMP" "${SAMPLE_ID}" "${ENTRY}"
                echo "$(date '+%Y-%m-%d %H:%M:%S')|${SAMPLE_ID}|${JOBID}|FAILED|unknown|${NODE}" >> ${STATUS_MONITORING}
                ((FAILED_COUNT++))
            fi
        elif [ "$QUEUE_STATE" = "R" ]; then
            echo "[${SAMPLE_ID}] Job ${JOBID}: RUNNING"
            ((RUNNING_COUNT++))
        elif [ "$QUEUE_STATE" = "PD" ]; then
            echo "[${SAMPLE_ID}] Job ${JOBID}: PENDING"
            ((PENDING_COUNT++))
        fi
    else
        # Parse sacct output
        JOB_STATE=$(echo "$JOB_INFO" | cut -d'|' -f2)
        ELAPSED=$(echo "$JOB_INFO" | cut -d'|' -f3)
        CURRENT_NODE=$(echo "$JOB_INFO" | cut -d'|' -f4)

        case "$JOB_STATE" in
            COMPLETED)
                echo "[${SAMPLE_ID}] Job ${JOBID}: COMPLETED (${ELAPSED})"
                ENTRY="${SAMPLE_ID}|COMPLETED|${JOBID}|${CPUS}|${PARTITION}|${CURRENT_NODE}|$(date '+%Y-%m-%d %H:%M:%S')"
                update_status_file "${STATUS_SUBMISSION}" "# SAMPLE_ID|STATUS|JOBID|CPUS|PARTITION|NODE|TIMESTAMP" "${SAMPLE_ID}" "${ENTRY}"
                echo "$(date '+%Y-%m-%d %H:%M:%S')|${SAMPLE_ID}|${JOBID}|COMPLETED|${ELAPSED}|${CURRENT_NODE}" >> ${STATUS_MONITORING}
                ((COMPLETED_COUNT++))
                ;;
            RUNNING)
                echo "[${SAMPLE_ID}] Job ${JOBID}: RUNNING (${ELAPSED} on ${CURRENT_NODE})"
                # Update node if it changed
                if [ "$CURRENT_NODE" != "$NODE" ] && [ -n "$CURRENT_NODE" ]; then
                    ENTRY="${SAMPLE_ID}|RUNNING|${JOBID}|${CPUS}|${PARTITION}|${CURRENT_NODE}|$(date '+%Y-%m-%d %H:%M:%S')"
                    update_status_file "${STATUS_SUBMISSION}" "# SAMPLE_ID|STATUS|JOBID|CPUS|PARTITION|NODE|TIMESTAMP" "${SAMPLE_ID}" "${ENTRY}"
                fi
                echo "$(date '+%Y-%m-%d %H:%M:%S')|${SAMPLE_ID}|${JOBID}|RUNNING|${ELAPSED}|${CURRENT_NODE}" >> ${STATUS_MONITORING}
                ((RUNNING_COUNT++))
                ;;
            PENDING)
                echo "[${SAMPLE_ID}] Job ${JOBID}: PENDING"
                echo "$(date '+%Y-%m-%d %H:%M:%S')|${SAMPLE_ID}|${JOBID}|PENDING|00:00:00|none" >> ${STATUS_MONITORING}
                ((PENDING_COUNT++))
                ;;
            FAILED|NODE_FAIL|TIMEOUT|CANCELLED|OUT_OF_MEMORY)
                echo "[${SAMPLE_ID}] Job ${JOBID}: FAILED (${JOB_STATE})"
                ENTRY="${SAMPLE_ID}|FAILED|${JOBID}|${CPUS}|${PARTITION}|${CURRENT_NODE}|$(date '+%Y-%m-%d %H:%M:%S')"
                update_status_file "${STATUS_SUBMISSION}" "# SAMPLE_ID|STATUS|JOBID|CPUS|PARTITION|NODE|TIMESTAMP" "${SAMPLE_ID}" "${ENTRY}"
                echo "$(date '+%Y-%m-%d %H:%M:%S')|${SAMPLE_ID}|${JOBID}|FAILED-${JOB_STATE}|${ELAPSED}|${CURRENT_NODE}" >> ${STATUS_MONITORING}
                ((FAILED_COUNT++))
                ;;
            *)
                echo "[${SAMPLE_ID}] Job ${JOBID}: ${JOB_STATE} (${ELAPSED})"
                ;;
        esac
    fi
done <<< "$ACTIVE_JOBS"

echo ""
echo "=========================================="
echo "Summary:"
echo "  Running:   ${RUNNING_COUNT}"
echo "  Pending:   ${PENDING_COUNT}"
echo "  Completed: ${COMPLETED_COUNT} (this check)"
echo "  Failed:    ${FAILED_COUNT} (this check)"
echo "=========================================="
echo ""

# Overall status from file
TOTAL_RUNNING=$(grep "|RUNNING|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)
TOTAL_COMPLETED=$(grep "|COMPLETED|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)
TOTAL_PENDING=$(grep "|PENDING|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)
TOTAL_FAILED=$(grep "|FAILED|" ${STATUS_SUBMISSION} 2>/dev/null | wc -l)

echo "Overall Status (from status file):"
echo "  Running:   ${TOTAL_RUNNING}"
echo "  Completed: ${TOTAL_COMPLETED}"
echo "  Pending:   ${TOTAL_PENDING}"
echo "  Failed:    ${TOTAL_FAILED}"
echo ""
echo "Status file: ${STATUS_SUBMISSION}"
echo "Monitor log: ${STATUS_MONITORING}"
echo "=========================================="
