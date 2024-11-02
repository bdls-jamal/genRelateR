#' Install and load required packages for genetic analysis
#'
#' This script installs and loads all necessary packages for the genetic analysis functions.
#' It checks for existing installations and only installs missing packages.

# Function to safely install and load packages
setupGeneticPackages <- function() {
  # List of required CRAN packages
  cranPackages <- c(
    "data.table",
    "ggplot2",
    "reshape2",
    "maps",
    "viridis"
  )

  # List of required Bioconductor packages
  biocPackages <- c(
    "VariantAnnotation",
    "GenomicRanges",
    "SNPRelate"
  )

  # Install BiocManager if not present
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  # Install missing CRAN packages
  for (pkg in cranPackages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }

  # Install missing Bioconductor packages
  for (pkg in biocPackages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, update = FALSE)
    }
  }

  # Load all required packages
  library(VariantAnnotation)
  library(GenomicRanges)
  library(SNPRelate)
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(maps)
  library(viridis)

  # Return success message
  cat("All required packages have been installed and loaded successfully.\n")
}

# Run the setup
setupGeneticPackages()
