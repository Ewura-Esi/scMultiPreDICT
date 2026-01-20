#!/usr/bin/env Rscript
# ============================================================================
# scMultiPreDICT - R Package Installation Script
# ============================================================================
#
# Description:
#   Automated installation of all R package dependencies required for the
#   scMultiPreDICT analysis pipeline.
#
# Usage:
#   source("install_packages.R")  # From R console
#   Rscript install_packages.R    # From command line
#
# ============================================================================

cat("============================================================\n")
cat("      scMultiPreDICT R Package Installation                 \n")
cat("============================================================\n\n")

# ============================================================================
# CRAN Packages
# ============================================================================

cran_packages <- c(
  # Data manipulation
  "dplyr",
  "tidyr",
  "readr",
  "purrr",
  "stringr",
  "tibble",
  "forcats",
  
  # Visualization
  "ggplot2",
  "patchwork",
  "viridis",
  "scales",
  "ggrepel",
  "Cairo",
  "RColorBrewer",
  "cowplot",
  "reshape2",
  
  # Single-cell analysis
  "Seurat",
  "Matrix",
  
  # Machine learning
  "glmnet",
  "ranger",
  "caret",
  
  # Utilities
  "RANN",
  "irlba",
  "reticulate",
  "parallel",
  "doParallel",
  "foreach",
  "matrixStats"
)

cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    cat(sprintf("  Already installed: %s\n", pkg))
  }
}

# ============================================================================
# Bioconductor Packages
# ============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "Signac",
  "GenomicRanges",
  "GenomeInfoDb",
  "rtracklayer",
  "IRanges",
  "S4Vectors"
)

cat("\nInstalling Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing: %s\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(sprintf("  Already installed: %s\n", pkg))
  }
}

# ============================================================================
# Species-Specific Annotation Packages
# ============================================================================

cat("\n============================================================\n")
cat("Species-Specific Annotation Packages\n")
cat("============================================================\n")
cat("Install ONE of the following based on your organism:\n\n")

cat("For Mus musculus (mouse):\n")
cat("  BiocManager::install('EnsDb.Mmusculus.v79')\n")
cat("  BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')\n\n")

cat("For Homo sapiens (human):\n")
cat("  BiocManager::install('EnsDb.Hsapiens.v86')\n")
cat("  BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')\n\n")

# Interactive species selection
species_choice <- readline(prompt = "Install annotation packages? (mouse/human/skip): ")

if (tolower(species_choice) == "mouse") {
  cat("Installing mouse annotation packages...\n")
  BiocManager::install("EnsDb.Mmusculus.v79", update = FALSE, ask = FALSE)
  BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", update = FALSE, ask = FALSE)
} else if (tolower(species_choice) == "human") {
  cat("Installing human annotation packages...\n")
  BiocManager::install("EnsDb.Hsapiens.v86", update = FALSE, ask = FALSE)
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update = FALSE, ask = FALSE)
} else {
  cat("Skipping species-specific packages. Install manually when needed.\n")
}

# ============================================================================
# Installation Verification
# ============================================================================

cat("\n============================================================\n")
cat("Installation Verification\n")
cat("============================================================\n")

all_packages <- c(cran_packages, bioc_packages)
missing <- c()

for (pkg in all_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing <- c(missing, pkg)
  }
}

if (length(missing) > 0) {
  cat("\nWARNING: The following packages failed to install:\n")
  cat(paste("  -", missing, collapse = "\n"), "\n")
  cat("\nPlease install these packages manually.\n")
} else {
  cat("\nAll core packages installed successfully.\n")
}

# ============================================================================
# Session Information
# ============================================================================

cat("\n============================================================\n")
cat("Session Information\n")
cat("============================================================\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Platform: %s\n", R.version$platform))
if (requireNamespace("BiocManager", quietly = TRUE)) {
  cat(sprintf("BiocManager version: %s\n", as.character(packageVersion("BiocManager"))))
}
if (requireNamespace("Seurat", quietly = TRUE)) {
  cat(sprintf("Seurat version: %s\n", as.character(packageVersion("Seurat"))))
}
if (requireNamespace("Signac", quietly = TRUE)) {
  cat(sprintf("Signac version: %s\n", as.character(packageVersion("Signac"))))
}

cat("\n============================================================\n")
cat("Installation Complete\n")
cat("============================================================\n")
