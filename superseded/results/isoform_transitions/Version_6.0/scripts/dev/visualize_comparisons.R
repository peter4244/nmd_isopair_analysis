#!/usr/bin/env Rscript
#
# Visualize random sample of isoform pairs from comparison pipeline
#
# Usage:
#   Rscript scripts/dev/visualize_comparisons.R                                    # defaults
#   Rscript scripts/dev/visualize_comparisons.R --profiles path.rds --n 30         # custom
#   Rscript scripts/dev/visualize_comparisons.R --failures                         # only FAIL/ERROR
#   Rscript scripts/dev/visualize_comparisons.R --seed 123 --output out.pdf        # reproducible
#

library(tidyverse)
library(ggplot2)
library(patchwork)

# Source core libraries
source("scripts/core/event_detection_functions.R")
source("scripts/core/reconstruction_functions.R")
source("scripts/core/visualization_functions.R")

# Parse CLI arguments
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default) {
  if (flag %in% args) {
    idx <- which(args == flag)
    if (idx < length(args)) return(args[idx + 1])
  }
  default
}

profiles_rds  <- get_arg("--profiles", "comparisons/deduplicated/all_splicing_profiles.rds")
structures_rds <- get_arg("--structures", "data/isoform_structures_filtered.rds")
ue_rds        <- get_arg("--ue", "data/union_exons_filtered.rds")
output_pdf    <- get_arg("--output", "comparisons/deduplicated/visual_validation.pdf")
n_sample      <- as.integer(get_arg("--n", "20"))
seed          <- as.integer(get_arg("--seed", "42"))
failures_only <- "--failures" %in% args

visualize_profiles(
  profiles_rds = profiles_rds,
  structures_rds = structures_rds,
  union_exons_rds = ue_rds,
  output_pdf = output_pdf,
  n_sample = n_sample,
  seed = seed,
  failures_only = failures_only
)
