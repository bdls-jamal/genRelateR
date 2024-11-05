#' Load and Format Genetic Data
#'
#' This function loads genetic data from the specified file path and formats it for analysis.
#' By default, it loads a manageable subset of the chromosome to prevent memory issues.
#' If the default region is still too large, please specify a smaller region using the
#' regions parameter (e.g., GRanges("1", IRanges(1e6, 2e6)) for chromosome 1, 1-2Mb).
#'
#' @param vcf_path Path to VCF file
#' @param regions GRanges object specifying regions to load. If NULL, automatically selects
#'               a default 1Mb region based on the chromosome in the filename.
#' @param samples Vector of sample IDs to include. If NULL, includes all samples.
#' @return A VCF object containing the genetic data
#' @examples
#' # Example with default region:
#' vcf_file <- "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
#' genetic_data <- loadGeneticData(vcf_file)
#'
#' # Example with custom region:
#' regions <- GRanges("1", IRanges(1e6, 2e6))
#' genetic_data <- loadGeneticData(vcf_file, regions = regions)
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

  # If no region specified, try to determine chromosome from filename and set default region
  if (is.null(regions)) {
    # Extract chromosome number from filename
    chr_match <- regexpr("chr([0-9XY]+)", basename(vcf_path), perl=TRUE)
    if (chr_match > 0) {
      chr <- regmatches(basename(vcf_path), chr_match)[[1]]
      chr_num <- sub("chr", "", chr)

      # Set default region to first 1Mb of the chromosome
      message(sprintf("No region specified. Loading first 1Mb of %s (position 1-1,000,000)...", chr))
      regions <- GRanges(chr_num, IRanges(1, 1e6))
    } else {
      stop("Could not determine chromosome from filename. Please specify a region.")
    }
  }

  # Create TabixFile object
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
#'
#' @param vcf_data A VCF object loaded using loadGeneticData
#' @param pop_file Path to population metadata file (tab-delimited with columns: Sample, Population, SuperPop)
#' @param population Vector of population codes to include (e.g., c("CEU", "YRI"))
#' @param super_pop Vector of super-population codes (e.g., c("EUR", "AFR"))
#' @return Filtered VCF object
#' @export
#' @import VariantAnnotation dplyr readr stringr
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
    filter(sample %in% vcf_samples)

  if (!is.null(population)) {
    selected_samples <- selected_samples %>%
      filter(pop %in% population)
  }

  if (!is.null(super_pop)) {
    selected_samples <- selected_samples %>%
      filter(super_pop %in% super_pop)
  }

  # Append pop code and name mapping to metadata
  pop_names_df <- read.table("data/population_codes.txt", header = FALSE, col.names = c("Code", "Name"))
  pop_metadata$population <- pop_names_df$Name[match(pop_metadata$pop, pop_names_df$Code)]

  # Append long and lat of each pop code to metadata
  pop_long_lat <- read.table("data/population_long_lat.txt", header = TRUE, sep = "\t")
  pop_metadata$Longitude <- pop_long_lat$Longitude[match(pop_metadata$pop, pop_long_lat$Name)]
  pop_metadata$Latitude <- pop_long_lat$Latitude[match(pop_metadata$pop, pop_long_lat$Name)]

  # Subset VCF to selected samples
  filtered_vcf <- vcf_data[, selected_samples$sample]

  # Add attributes about filtering
  attr(filtered_vcf, "n_samples_removed") <- length(vcf_samples) - ncol(filtered_vcf)
  attr(filtered_vcf, "populations_included") <- unique(selected_samples$pop)

  return(list(vcf_data = filtered_vcf, pop_metadata = pop_metadata))
}
