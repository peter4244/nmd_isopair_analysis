# Install required packages for isoform proportion analysis

packages_needed <- c("tidyverse", "data.table", "DT")

for (pkg in packages_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat(pkg, "is already installed\n")
  }
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  cat("Installing rtracklayer...\n")
  BiocManager::install("rtracklayer", update = FALSE, ask = FALSE)
} else {
  cat("rtracklayer is already installed\n")
}

cat("\nPackage installation complete!\n")
