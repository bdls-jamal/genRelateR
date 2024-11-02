#' Load and Format Genetic Data
#'
#' This function loads genetic data from the specified file path and formats it for analysis.
#'
#' @param file_path Path to VCF file
#' @param genome Reference genome version
#' @param param ScanVcfParam object for filtering
#' @param yield_size Number of records to read at once
#' @return Data frame of genetic data
#'
#' @examples
#' # Example 1:
#' Example usage with a yieldSize and specific region of interest
#' vcf_file <- "../../data/vcf/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
#' param <- ScanVcfParam(which = GRanges("1", IRanges(1e6, 1.5e6)))  # Specify a range
#' genetic_data <- loadGeneticData(vcf_file, genome = "b37", param = param, yield_size = 1000)
#'
#' @references
#' Example:
#'Akaike, H. (1973). Information theory and an extension of the maximum
#'likelihood principle. In \emph{Second International Symposium on Information
#'Theory}, New York, NY, USA, pp. 267–281. Springer Verlag. \href{https://link.springer.com/chapter/10.1007/978-1-4612-1694-0_15}{Link}
#' @export
#' @import VariantAnnotation GenomicRanges data.table


# Install the package if you haven't already
source("R/setupPackages.R")

# loadGeneticData <- function(vcf_path, regions = NULL, samples = NULL) {
#   # Validate input file
#   if (!file.exists(vcf_path)) {
#     stop("VCF file does not exist: ", vcf_path)
#   }
#
#   # Set up basic scanning parameters
#   scan_params <- ScanVcfParam()
#
#   # Add region filtering if specified
#   if (!is.null(regions)) {
#     scan_params <- ScanVcfParam(which=regions)
#   }
#
#   # Add sample filtering if specified
#   if (!is.null(samples)) {
#     scan_params <- ScanVcfParam(samples=samples)
#   }
#
#   # Load the VCF file
#   message("Loading VCF file...")
#   vcf_data <- readVcf(vcf_path, genome="hg19", param=scan_params)
#
#   return(vcf_data)
# }

loadGeneticData <- function(vcf_path, regions = NULL, samples = NULL) {
  # Validate input file
  if (!file.exists(vcf_path)) {
    stop("VCF file does not exist: ", vcf_path)
  }

  # Check for tabix index
  tabix_path <- paste0(vcf_path, ".tbi")
  if (!file.exists(tabix_path)) {
    stop("Tabix index (.tbi) not found. Please ensure the VCF file is properly indexed.")
  }

  # Create TabixFile object without yieldSize for compatibility with range queries
  tbx <- TabixFile(vcf_path)

  # Set up scanning parameters
  if (!is.null(regions) && !is.null(samples)) {
    scan_params <- ScanVcfParam(which=regions, samples=samples)
  } else if (!is.null(regions)) {
    scan_params <- ScanVcfParam(which=regions)
  } else if (!is.null(samples)) {
    scan_params <- ScanVcfParam(samples=samples)
  } else {
    scan_params <- ScanVcfParam()
  }

  # Load the VCF file
  message("Loading VCF file using tabix index...")
  vcf_data <- readVcf(tbx, genome="hg19", param=scan_params)

  return(vcf_data)
}


#' Filter Population Data
#'
#' This function filters the genetic data based on population, geographic region, or genetic criteria.
#'
#' @param genetic_data A data frame containing the genetic data to be filtered.
#' @param population Optional; A string specifying the population to filter by.
#' @param region Optional; A string specifying the geographic region to filter by.
#' @param criteria Optional; A list of named vectors representing additional genetic criteria (e.g., SNP values).
#' @return A filtered data frame.
#' @export
#' @import dplyr

filterPopulation <- function(genetic_data, population = NULL, region = NULL, criteria = NULL) {
  # Load necessary library
  library(dplyr)  # for data manipulation

  # Start filtering the data
  filtered_data <- genetic_data

  # Filter by population
  if (!is.null(population)) {
    filtered_data <- filtered_data %>% filter(Population == population)
  }

  # Filter by geographic region
  if (!is.null(region)) {
    filtered_data <- filtered_data %>% filter(GeographicRegion == region)
  }

  # Filter by additional criteria
  if (!is.null(criteria)) {
    for (criterion in names(criteria)) {
      filtered_data <- filtered_data %>% filter(!!sym(criterion) %in% criteria[[criterion]])
    }
  }

  # Return the filtered data
  return(filtered_data)
}
