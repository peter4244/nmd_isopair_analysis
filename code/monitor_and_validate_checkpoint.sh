#!/bin/bash
# Monitor event detection and run validation at checkpoint

OUTPUT_FILE="/private/tmp/claude-502/-Users-petecastaldi-claude-projects/tasks/b0c54cb.output"
CHECKPOINT_FILE="results/isoform_transitions/v4.0_reference_based/event_vectors_full.rds"

echo "Monitoring for checkpoint at 1000 genes..."
echo "Checking every 30 seconds..."
echo ""

while true; do
    # Check if checkpoint message appears
    if grep -q "CHECKPOINT.*Saving at gene 1000" "$OUTPUT_FILE" 2>/dev/null; then
        echo "[$(date)] Checkpoint 1000 detected!"
        break
    fi

    # Show current progress
    CURRENT=$(tail -1 "$OUTPUT_FILE" 2>/dev/null | grep -oE '\[[0-9]+/[0-9]+\]' | head -1)
    if [ ! -z "$CURRENT" ]; then
        echo "[$(date)] Progress: $CURRENT"
    fi

    sleep 30
done

echo ""
echo "Waiting 10 seconds for file to be written..."
sleep 10

echo "Running validation..."
echo ""
Rscript code/validate_checkpoint_events.R

echo ""
echo "[$(date)] Validation complete!"
