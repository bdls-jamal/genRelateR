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
# Advise user to run setupPackages.R first

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
#' This function filters VCF data based on population metadata from the 1000 Genomes Project.
#' Includes optional filtering of related individuals.
#'
#' @param vcf_data A VCF object loaded using loadGeneticData
#' @param pop_file Path to population metadata file (tab-delimited with columns: Sample, Population, SuperPop)
#' @param rel_file Path to relationship metadata file (tab-delimited with columns: Sample, Population, Gender, Reason_for_exclusion)
#' @param population Vector of population codes to include (e.g., c("CEU", "YRI"))
#' @param super_pop Vector of super-population codes (e.g., c("EUR", "AFR"))
#' @param remove_related Logical indicating whether to remove related individuals (default: FALSE)
#' @param prioritize_gender Character indicating which gender to prioritize when removing related pairs ("male" or "female", default: NULL)
#' @return Filtered VCF object
#' @export
#' @import VariantAnnotation dplyr readr stringr
filterPopulation <- function(vcf_data, pop_file, rel_file, population = NULL,
                             super_pop = NULL, prioritize_gender = NULL) {
  # Load population metadata
  pop_metadata <- read_tsv(pop_file, col_types = cols())
  pop_metadata<- pop_metadata[, c(1, 2, 3, 4)]

  # Validate input
  if (!is(vcf_data, "VCF")) {
    stop("Input must be a VCF object")
  }

  # Get sample IDs from VCF
  vcf_samples <- colnames(vcf_data)

  # Filter samples based on population criteria
  selected_samples <- pop_metadata %>%
    filter(sample %in% vcf_samples)

  if (!is.null(population)) {
    selected_samples <- selected_samples %>%
      filter(pop %in% population)
  }

  if (!is.null(super_pop)) {
    selected_samples <- selected_samples %>%
      filter(super_pop %in% super_pop)
  }

  # Handle relationship-based filtering if requested
  if (!is.null(rel_file)) {
    rel_metadata <- read_tsv(rel_file, col_types = cols())

    # Extract related pairs
    related_pairs <- rel_metadata %>%
      mutate(
        related_sample = str_extract("Reason for exclusion", "(?<=:)[^:]+$")
      ) %>%
      select(Sample, related_sample, Gender, Population)

    # Determine which samples to remove based on relationships
    samples_to_remove <- character(0)

    if (!is.null(prioritize_gender)) {
      # Prioritize keeping specific gender when possible
      for (i in 1:nrow(related_pairs)) {
        pair <- related_pairs[i,]
        if (pair$Gender == prioritize_gender) {
          samples_to_remove <- c(samples_to_remove, pair$related_sample)
        } else {
          samples_to_remove <- c(samples_to_remove, pair$Sample)
        }
      }
    } else {
      # If no gender preference, remove the first sample of each pair
      samples_to_remove <- related_pairs$Sample
    }

    # Remove related samples from selection
    selected_samples <- selected_samples %>%
      filter(!(sample %in% samples_to_remove))
  }

  # Subset VCF to selected samples
  filtered_vcf <- vcf_data[, selected_samples$sample]

  # Add attributes about filtering
  attr(filtered_vcf, "n_samples_removed") <- length(vcf_samples) - ncol(filtered_vcf)
  attr(filtered_vcf, "populations_included") <- unique(selected_samples$pop)
  if (exists("samples_to_remove")) {
    attr(filtered_vcf, "related_samples_removed") <- samples_to_remove
  }

  return(list(vcf_data = filtered_vcf, pop_metadata = pop_metadata, rel_data = rel_metadata))
}
