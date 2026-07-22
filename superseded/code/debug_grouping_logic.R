#!/usr/bin/env Rscript
# Debug the grouping logic with detailed iteration tracing

library(tidyverse)

TSS_TES_TOLERANCE <- 20

# Manual implementation with debugging
group_first_exons_debug <- function(first_exons, dominant_iso) {
  if (nrow(first_exons) == 0) return(list())

  dominant_first <- first_exons %>% filter(isoform_id == dominant_iso) %>% slice(1)
  if (nrow(dominant_first) == 0) dominant_first <- first_exons %>% slice(1)

  cat("Dominant first exon:\n")
  print(dominant_first %>% select(isoform_id, start, end))
  cat("\n")

  first_groups <- list()
  processed <- rep(FALSE, nrow(first_exons))
  group_id <- 1

  for (i in seq_len(nrow(first_exons))) {
    cat(sprintf("Iteration %d:\n", i))
    cat(sprintf("  Current: %s (start=%d, end=%d)\n",
                first_exons$isoform_id[i],
                first_exons$start[i],
                first_exons$end[i]))

    if (processed[i]) {
      cat("  → SKIPPED (already processed)\n\n")
      next
    }

    current <- first_exons[i, ]

    # Show the grouping condition
    cat(sprintf("  Grouping condition: end==%d AND abs(start-%d)<=20\n",
                current$end, dominant_first$start))

    same_group <- which(first_exons$end == current$end &
                       abs(first_exons$start - dominant_first$start) <= TSS_TES_TOLERANCE)

    cat(sprintf("  Matched rows: %s\n", paste(same_group, collapse=", ")))
    cat("  Matched isoforms:\n")
    print(first_exons[same_group, ] %>% select(isoform_id, start, end))

    first_groups[[group_id]] <- first_exons[same_group, ]
    processed[same_group] <- TRUE

    cat(sprintf("  → Created Group %d\n", group_id))
    cat(sprintf("  → Marked rows %s as processed\n\n", paste(same_group, collapse=", "))
    )

    group_id <- group_id + 1
  }

  return(first_groups)
}

# Load SUOX first exons from previous debug
output_dir <- "results/isoform_transitions/v3.0_reference_based"
exon_structures <- readRDS(file.path(output_dir, "exon_structures_by_isoform_full.rds"))

# Get SUOX first exons manually
suox_first_exons <- tribble(
  ~isoform_id,           ~start,   ~end,
  "ENST00000886413.1",  55997022, 55997723,
  "ENST00000886414.1",  55997259, 55997339,
  "ENST00000886426.1",  55997301, 55997339,
  "ENST00000886434.1",  55997330, 55997537,
  "ENST00000886436.1",  55997333, 55997537,
  "ENST00000930303.1",  55997330, 55997723,
  "ENST00000930304.1",  55997320, 55997339,
  "ENST00000930305.1",  55997320, 55997339,
  "ENST00000954615.1",  55997286, 55997339
)

cat("═══════════════════════════════════════════════════════════════\n")
cat("DEBUGGING FIRST EXON GROUPING FOR SUOX\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

dominant <- "ENST00000930304.1"
groups <- group_first_exons_debug(suox_first_exons, dominant)

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("FINAL RESULT:\n")
cat(sprintf("Number of groups: %d\n\n", length(groups)))

for (i in seq_along(groups)) {
  cat(sprintf("Group %d:\n", i))
  print(groups[[i]])
  cat("\n")
}
