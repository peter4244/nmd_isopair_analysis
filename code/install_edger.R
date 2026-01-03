#!/usr/bin/env Rscript
# Install edgeR for DGEList handling

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cat("Installing edgeR...\n")
BiocManager::install("edgeR", update = FALSE, ask = FALSE)

cat("\nTesting edgeR installation:\n")
library(edgeR)
cat("edgeR version:", as.character(packageVersion("edgeR")), "\n")
