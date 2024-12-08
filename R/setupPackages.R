#' Install and load required packages for genetic analysis
#'
#' This script installs and loads all necessary packages for the genetic analysis functions.
#' It checks for existing installations and only installs missing packages.
#' @export
setupGeneticPackages <- function() {
  # List of required CRAN packages
  cranPackages <- c(
    "data.table",
    "ggplot2",
    "maps",
    "viridis",
    "tidyverse",
    "readr",
    "dplyr",
    "tidyr",
    "stringr",
    "ggrepel",
    "igraph",
    "RColorBrewer",
    "plotly"
  )
  # List of required Bioconductor packages
  biocPackages <- c(
    "VariantAnnotation",
    "GenomicRanges"
  )
  # Install BiocManager if not present
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  # Install missing CRAN packages
  for (pkg in cranPackages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Installing %s from CRAN...", pkg))
      install.packages(pkg)
    }
  }
  # Install missing Bioconductor packages
  for (pkg in biocPackages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Installing %s from Bioconductor...", pkg))
      BiocManager::install(pkg, update = FALSE)
    }
  }
  # Load all required packages with error handling
  packages_to_load <- c(
    # Bioconductor packages
    "VariantAnnotation",
    "GenomicRanges",
    # CRAN packages
    "data.table",
    "ggplot2",
    "maps",
    "viridis",
    "tidyverse",
    "readr",
    "dplyr",
    "tidyr",
    "stringr",
    "ggrepel",
    "igraph",
    "RColorBrewer",
    "plotly"
  )
  for (pkg in packages_to_load) {
    tryCatch({
      library(pkg, character.only = TRUE)
      message(sprintf("Successfully loaded %s", pkg))
    }, error = function(e) {
      warning(sprintf("Failed to load %s: %s", pkg, e$message))
    })
  }
  # Check if all packages were successfully loaded
  loaded_packages <- (.packages())
  missing_packages <- setdiff(packages_to_load, loaded_packages)
  if (length(missing_packages) > 0) {
    warning(sprintf(
      "The following packages failed to load: %s",
      paste(missing_packages, collapse = ", ")
    ))
  } else {
    message("All required packages have been installed and loaded successfully.")
  }
  # Return invisibly whether all packages were loaded successfully
  invisible(length(missing_packages) == 0)
}
#' Function to verify all required packages are available
#' @return Logical indicating if all required packages are available
#' @export
checkGeneticPackages <- function() {
  required_packages <- c(
    # Bioconductor
    "VariantAnnotation",
    "GenomicRanges",
    # CRAN
    "data.table",
    "ggplot2",
    "maps",
    "viridis",
    "tidyverse",
    "readr",
    "dplyr",
    "tidyr",
    "stringr",
    "ggrepel",
    "igraph",
    "RColorBrewer",
    "plotly"
  )
  missing_packages <- required_packages[
    !sapply(required_packages, requireNamespace, quietly = TRUE)
  ]
  if (length(missing_packages) > 0) {
    warning(sprintf(
      "Missing required packages: %s\nPlease run setupGeneticPackages()",
      paste(missing_packages, collapse = ", ")
    ))
    return(FALSE)
  }
  return(TRUE)
}
# If this file is being sourced directly, run the setup
if (sys.nframe() == 0) {
  setupGeneticPackages()
  source("R/loadAndCleanData.R")
  source("R/computationAndAnalysis.R")
  source("R/visualizationScripts.R")
}
