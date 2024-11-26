#' Load and Format Genetic Data
#'
#' This function loads genetic data from the specified file path and formats it for analysis.
#' By default, it loads a manageable subset of the chromosome to prevent memory issues.
#' If the default region is still too large, please specify a smaller region using the
#' regions parameter (e.g., GRanges("1", IRanges(1e6, 2e6)) for chromosome 1, 1-2Mb).
#'
#' Warning: If you are receiving a tbi file is older than the vcf file error,
#' Run the following commands to recreate the tbi file:
#' library(VariantAnnotation)
#' indexTabix("path/to/vcf.gz/file", format = "vcf")
#'
#' @param vcf_path Path to VCF file
#' @param regions GRanges object specifying regions to load. If NULL, automatically selects
#'               a default 1Mb region based on the chromosome in the filename.
#' @param samples Vector of sample IDs to include. If NULL, includes all samples.
#' @return A VCF object containing the genetic data
#' @examples
#' # Example with default region:
#' vcf_file <- "../inst/ext/data/vcf/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
#' genetic_data <- genRelateR::loadGeneticData(vcf_file)
#'
#' # Example with custom region:
#' regions <- GenomicRanges::GRanges("1", IRanges(1e6, 2e6))
#' genetic_data <- genRelateR::loadGeneticData(vcf_file, regions = regions)
#'
#' @export
#' @importFrom VariantAnnotation readVcf
#' @importFrom GenomicRanges GRanges
#' @importFrom data.table fread
loadGeneticData <- function(vcf_path, regions = NULL, samples = NULL) {
  # Validate input file path
  vcf_path <- normalizePath(vcf_path, mustWork = TRUE)

  # Check for tabix index
  tabix_path <- paste0(vcf_path, ".tbi")

  # Debug print statements (you can remove these in production)
  message("VCF Path: ", vcf_path)
  message("Tabix Path: ", tabix_path)

  # Enhanced file existence checks
  if (!file.exists(vcf_path)) {
    stop("VCF file does not exist: ", vcf_path)
  }

  if (!file.exists(tabix_path)) {
    # Additional check for alternative tabix index naming
    alternative_tabix_path <- file.path(dirname(vcf_path),
                                        paste0(basename(vcf_path), ".tbi"))

    if (file.exists(alternative_tabix_path)) {
      tabix_path <- alternative_tabix_path
    } else {
      stop("Tabix index (.tbi) not found.
            Checked paths: \n",
           tabix_path, "\n",
           alternative_tabix_path, "\n",
           "Please ensure the VCF file is properly indexed.")
    }
  }

  # If no region specified, try to determine chromosome from filename and set default region
  if (is.null(regions)) {
    # Extract chromosome number from filename
    chr_match <- regexpr("chr([0-9XY]+)", basename(vcf_path), perl=TRUE)
    if (chr_match > 0) {
      chr <- regmatches(basename(vcf_path), chr_match)[[1]]
      chr_num <- sub("chr", "", chr)

      # Set default region to first 1Mb of the chromosome
      message(sprintf("No region specified. Loading first 1Mb of %s (position 1-1,000,000)...", chr))
      regions <- GenomicRanges::GRanges(chr_num, IRanges(1, 1e6))
    } else {
      stop("Could not determine chromosome from filename. Please specify a region.")
    }
  }

  # Create TabixFile object
  tbx <- TabixFile(vcf_path)

  # Set up scanning parameters
  if (!is.null(regions) && !is.null(samples)) {
    scan_params <- VariantAnnotation::ScanVcfParam(which=regions, samples=samples)
  } else if (!is.null(regions)) {
    scan_params <- VariantAnnotation::ScanVcfParam(which=regions)
  } else if (!is.null(samples)) {
    scan_params <- VariantAnnotation::ScanVcfParam(samples=samples)
  } else {
    scan_params <- VariantAnnotation::ScanVcfParam()
  }

  # Load the VCF file
  message("Loading VCF file using tabix index...")
  vcf_data <- VariantAnnotation::readVcf(tbx, genome="hg19", param=scan_params)

  return(vcf_data)
}


#' Filter Population Data
#'
#' This function filters VCF data based on population metadata from the 1000 Genomes Project.
#'
#' @param vcf_data A VCF object loaded using loadGeneticData
#' @param pop_file Path to population metadata file (tab-delimited with columns: Sample, Population, SuperPop)
#' @param population Vector of population codes to include (e.g., c("CEU", "YRI"))
#' @param super_pop Vector of super-population codes (e.g., c("EUR", "AFR"))
#' @return Filtered VCF object
#' @examples
#' # Load genetic data and metadata
#' genRelateR::loadGeneticDada(vcf_file)
#' pop_metadata <- "data/population_metadata.txt"
#'
#' # Filter populations
#' populations <- c(  "CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "GBR",
#' "FIN", "IBS", "YRI", "LWK", "GWD", "MSL", "ESN","ASW", "ACB", "MXL", "PUR",
#' "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU")
#'
#' filtered_data <- genRelateR::filterPopulation(genetic_data, pop_metadata, populations)
#'
#' @export
#' @importFrom VariantAnnotation readVcf
#' @importFrom dplyr filter select
#' @importFrom readr read_csv
#' @importFrom stringr str_detect
filterPopulation <- function(vcf_data, pop_file, population = NULL,
                             super_pop = NULL) {
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
    dplyr::filter(sample %in% vcf_samples)

  if (!is.null(population)) {
    selected_samples <- selected_samples %>%
      dplyr::filter(pop %in% population)
  }

  if (!is.null(super_pop)) {
    selected_samples <- selected_samples %>%
      dplyr::filter(super_pop %in% super_pop)
  }

  # Append pop code and name mapping to metadata
  pop_names_df <- read.table(
    system.file("extdata", "population_codes.txt", package = "genRelateR"),
    header = FALSE,
    col.names = c("Code", "Name")
  )
  pop_metadata$population <- pop_names_df$Name[match(pop_metadata$pop, pop_names_df$Code)]

  # Access population longitude and latitude
  pop_long_lat <- read.table(
    system.file("extdata", "population_long_lat.txt", package = "genRelateR"),
    header = TRUE,
    sep = "\t"
  )
  pop_metadata$Longitude <- pop_long_lat$Longitude[match(pop_metadata$pop, pop_long_lat$Name)]
  pop_metadata$Latitude <- pop_long_lat$Latitude[match(pop_metadata$pop, pop_long_lat$Name)]

  # Subset VCF to selected samples
  filtered_vcf <- vcf_data[, selected_samples$sample]

  # Add attributes about filtering
  attr(filtered_vcf, "n_samples_removed") <- length(vcf_samples) - ncol(filtered_vcf)
  attr(filtered_vcf, "populations_included") <- unique(selected_samples$pop)

  return(list(vcf_data = filtered_vcf, pop_metadata = pop_metadata))
}
