#!/bin/bash
# Monitor event detection to see if minus strand Alt_TSS/TES are being detected

echo "Monitoring for minus strand Alt_TSS/TES detection..."
echo "Checking every 30 seconds..."
echo

sleep 120  # Wait 2 minutes for some genes to be processed

while true; do
  timestamp=$(date)

  # Check if results file exists
  if [ -f "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v4.0_reference_based/event_vectors_full.rds" ]; then

    # Run R code to check strand distribution
    result=$(Rscript -e "
      suppressPackageStartupMessages(library(tidyverse))
      results <- readRDS('/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v4.0_reference_based/event_vectors_full.rds')

      # Check if gene_strand column exists
      if ('gene_strand' %in% names(results)) {
        # Get Alt_TSS/TES by strand
        summary <- results %>%
          group_by(gene_strand) %>%
          summarise(
            n_genes = n_distinct(gene_id),
            n_transitions = n(),
            n_alt_tss = sum(n_Alt_TSS > 0),
            n_alt_tes = sum(n_Alt_TES > 0),
            .groups = 'drop'
          )

        cat('\\n[$timestamp]\\n')
        print(summary)

        # Check if we have minus strand Alt_TSS/TES
        minus_data <- filter(summary, gene_strand == '-')
        if (nrow(minus_data) > 0 && (minus_data\$n_alt_tss > 0 || minus_data\$n_alt_tes > 0)) {
          cat('\\n✅ MINUS STRAND Alt_TSS/TES DETECTED!\\n')
          cat('Stopping monitor...\\n')
          quit(status = 99)  # Special exit code to stop monitoring
        }
      } else {
        cat('[$timestamp] gene_strand column not found yet\\n')
      }
    " 2>&1)

    echo "$result"

    # Check exit code
    if [ $? -eq 99 ]; then
      echo
      echo "SUCCESS: Minus strand Alt_TSS/TES detection confirmed!"
      exit 0
    fi

  else
    echo "[$timestamp] Waiting for results file to be created..."
  fi

  sleep 30
done
