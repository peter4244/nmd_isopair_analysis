#!/usr/bin/env Rscript
# Show plus and minus strand examples for terminal boundary processing

library(tidyverse)

# Load data
exon_structures <- readRDS('results/isoform_transitions/v4.0_reference_based/exon_structures_by_isoform_full.rds')

cat('\n')
cat('═══════════════════════════════════════════════════════════════\n')
cat('PLUS STRAND EXAMPLE: WDR6\n')
cat('═══════════════════════════════════════════════════════════════\n\n')

# Get a representative isoform (WDR6 should be plus strand)
suox <- exon_structures %>%
  filter(grepl('WDR6', gene_id, fixed = TRUE)) %>%
  slice(1)

if (nrow(suox) > 0) {
  cat('Isoform:', suox$isoform_id, '\n')
  cat('Strand:', ifelse(suox$strand == 1, '+', '-'), '\n')
  cat('Number of exons:', suox$n_exons, '\n\n')

  starts <- suox$exon_starts[[1]]
  ends <- suox$exon_ends[[1]]
  n <- length(starts)

  # First exon
  cat('FIRST EXON (Exon 1):\n')
  cat(sprintf('  Coordinates: %d - %d\n', starts[1], ends[1]))
  cat('  ┌─────────────────────────────────────────┐\n')
  cat('  │ START: TSS (Transcription Start Site)  │\n')
  cat('  │        → SKIP for A5SS/A3SS detection  │\n')
  cat('  │        → Used for Alt_TSS detection    │\n')
  cat('  │                                         │\n')
  cat('  │ END:   5\'SS (Donor site)               │\n')
  cat('  │        → PROCESS for A5SS detection    │\n')
  cat('  │        → A5SS CAN occur here!          │\n')
  cat('  └─────────────────────────────────────────┘\n\n')

  # Last exon
  cat('LAST EXON (Exon', n, '):\n')
  cat(sprintf('  Coordinates: %d - %d\n', starts[n], ends[n]))
  cat('  ┌─────────────────────────────────────────┐\n')
  cat('  │ START: 3\'SS (Acceptor site)            │\n')
  cat('  │        → PROCESS for A3SS detection    │\n')
  cat('  │        → A3SS CAN occur here!          │\n')
  cat('  │                                         │\n')
  cat('  │ END:   TES (Transcript End Site)       │\n')
  cat('  │        → SKIP for A5SS/A3SS detection  │\n')
  cat('  │        → Used for Alt_TES detection    │\n')
  cat('  └─────────────────────────────────────────┘\n\n')

  cat('Plus strand logic:\n')
  cat('  process_start <- !is_first_exon  // First exon: FALSE (skip TSS)\n')
  cat('  process_end   <- !is_last_exon   // Last exon:  FALSE (skip TES)\n\n')
}

cat('\n')
cat('═══════════════════════════════════════════════════════════════\n')
cat('MINUS STRAND EXAMPLE: SUOX\n')
cat('═══════════════════════════════════════════════════════════════\n\n')

# Get SUOX (confirmed minus strand)
ggta1 <- exon_structures %>%
  filter(gene_id == 'ENSG00000139531.14', isoform_id == 'ENST00000930304.1')

if (nrow(ggta1) > 0) {
  cat('Isoform:', ggta1$isoform_id, '\n')
  cat('Gene:', ggta1$gene_id, '\n')
  cat('Strand:', ifelse(ggta1$strand == 1, '+', '-'), '\n')
  cat('Number of exons:', ggta1$n_exons, '\n\n')

  starts <- ggta1$exon_starts[[1]]
  ends <- ggta1$exon_ends[[1]]
  n <- length(starts)

  cat('NOTE: On minus strand, transcription goes 3\' → 5\' in genomic coordinates\n')
  cat('      (higher coordinates → lower coordinates)\n\n')

  # First exon (transcript 5' end)
  cat('FIRST EXON (Exon 1, transcript 5\' end):\n')
  cat(sprintf('  Coordinates: %d - %d\n', starts[1], ends[1]))
  cat('  ┌─────────────────────────────────────────┐\n')
  cat('  │ START: 5\'SS (Donor site)               │\n')
  cat('  │        → PROCESS for A5SS detection    │\n')
  cat('  │        → A5SS CAN occur here!          │\n')
  cat('  │                                         │\n')
  cat('  │ END:   TSS (Transcription Start Site)  │\n')
  cat('  │        → SKIP for A5SS/A3SS detection  │\n')
  cat('  │        → Used for Alt_TSS detection    │\n')
  cat('  │                                         │\n')
  cat('  │ ⚠️  COORDINATES REVERSED vs biology!    │\n')
  cat('  └─────────────────────────────────────────┘\n\n')

  # Last exon (transcript 3' end)
  cat('LAST EXON (Exon', n, ', transcript 3\' end):\n')
  cat(sprintf('  Coordinates: %d - %d\n', starts[n], ends[n]))
  cat('  ┌─────────────────────────────────────────┐\n')
  cat('  │ START: TES (Transcript End Site)       │\n')
  cat('  │        → SKIP for A5SS/A3SS detection  │\n')
  cat('  │        → Used for Alt_TES detection    │\n')
  cat('  │                                         │\n')
  cat('  │ END:   3\'SS (Acceptor site)            │\n')
  cat('  │        → PROCESS for A3SS detection    │\n')
  cat('  │        → A3SS CAN occur here!          │\n')
  cat('  │                                         │\n')
  cat('  │ ⚠️  COORDINATES REVERSED vs biology!    │\n')
  cat('  └─────────────────────────────────────────┘\n\n')

  cat('Minus strand logic:\n')
  cat('  process_start <- !is_last_exon   // Last exon:  FALSE (skip TES at START)\n')
  cat('  process_end   <- !is_first_exon  // First exon: FALSE (skip TSS at END)\n\n')
}

cat('\n')
cat('═══════════════════════════════════════════════════════════════\n')
cat('SUMMARY\n')
cat('═══════════════════════════════════════════════════════════════\n\n')

cat('Key insight: A5SS and A3SS CAN occur at first/last exons!\n')
cat('We only skip the TSS/TES boundaries themselves.\n\n')

cat('Plus strand:\n')
cat('  - First exon has TSS + donor  → skip TSS, process donor (A5SS possible)\n')
cat('  - Last exon has acceptor + TES → process acceptor (A3SS possible), skip TES\n\n')

cat('Minus strand:\n')
cat('  - First exon has donor + TSS  → process donor (A5SS possible), skip TSS\n')
cat('  - Last exon has TES + acceptor → skip TES, process acceptor (A3SS possible)\n\n')

cat('The strand-aware logic correctly handles this by swapping which\n')
cat('coordinate (START vs END) corresponds to TSS/TES.\n')
