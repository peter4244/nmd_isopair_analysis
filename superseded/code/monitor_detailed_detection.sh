#!/bin/bash
# Monitor detailed event detection progress

TASK_OUTPUT="/private/tmp/claude-502/-Users-petecastaldi-claude-projects/tasks/be4af30.output"
RESULTS_DIR="/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions"
CHECK_INTERVAL=300  # Check every 5 minutes

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║      MONITORING DETAILED EVENT DETECTION PROGRESS              ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Started monitoring at: $(date)"
echo "Checking every ${CHECK_INTERVAL} seconds (5 minutes)"
echo ""

while true; do
    # Check if output file exists
    if [ ! -f "$TASK_OUTPUT" ]; then
        echo "[$(date +%H:%M:%S)] Task output file not found. Waiting..."
        sleep $CHECK_INTERVAL
        continue
    fi
    
    # Get latest progress line
    LATEST_PROGRESS=$(grep -E "\[[[:space:]]*[0-9]+/[[:space:]]*[0-9]+\]" "$TASK_OUTPUT" | tail -1)
    
    if [ -n "$LATEST_PROGRESS" ]; then
        echo "[$(date +%H:%M:%S)] $LATEST_PROGRESS"
    fi
    
    # Check if complete (look for "Complete!" or "COMPLETE" in output)
    if grep -q "COMPLETE" "$TASK_OUTPUT" 2>/dev/null; then
        echo ""
        echo "╔════════════════════════════════════════════════════════════════╗"
        echo "║           DETAILED EVENT DETECTION COMPLETE!                   ║"
        echo "╚════════════════════════════════════════════════════════════════╝"
        echo ""
        echo "Completed at: $(date)"
        
        # Check for output file
        if [ -f "$RESULTS_DIR/detailed_event_vectors.rds" ]; then
            echo "✅ Output file created: detailed_event_vectors.rds"
            ls -lh "$RESULTS_DIR/detailed_event_vectors.rds"
        else
            echo "⚠️  Output file not found yet. Checking logs..."
        fi
        
        echo ""
        echo "Next steps:"
        echo "  1. Run Q1-Q3 analysis:"
        echo "     Rscript code/detailed_event_analysis_q1_q3.R"
        echo ""
        echo "  2. Generate visualizations:"
        echo "     Rscript code/detailed_event_visualizations.R"
        echo ""
        echo "  3. Or run full pipeline:"
        echo "     bash code/run_detailed_analysis_pipeline.sh"
        echo ""
        
        break
    fi
    
    # Check if process failed
    if grep -q "Error\|ERROR" "$TASK_OUTPUT" 2>/dev/null; then
        echo ""
        echo "⚠️  Possible error detected. Last 20 lines of output:"
        tail -20 "$TASK_OUTPUT"
        echo ""
        break
    fi
    
    sleep $CHECK_INTERVAL
done
