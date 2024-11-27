# Helper functions needed for the main functions to work

#' Helper function to validate data frame columns
#' @param df Data frame to validate
#' @param required_cols Expected columns in the frame
#' @param df_name Name of the data frame for error message
validate_columns <- function(df, required_cols, df_name) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in %s: %s",
                 df_name, paste(missing_cols, collapse = ", ")))
  }
}

#' Filter variants by variance
#' @param geno_mat Numeric genotype matrix
#' @param min_var Minimum variance threshold
#'
#' @return Filtered genotype matrix
filterVariantsByVariance <- function(geno_mat, min_var = 1e-10) {
  # Calculate variance for each variant (row)
  vars <- apply(geno_mat, 1, var, na.rm = TRUE)

  # Keep variants with variance above threshold
  keep_variants <- vars > min_var & !is.na(vars)

  if (sum(keep_variants) == 0) {
    stop("No variants remain after variance filtering")
  }

  return(geno_mat[keep_variants, , drop = FALSE])
}

#' Function to extract genotype data from CollapsedVCF
#' @param collapsed_vcf Vcf data in the collapsed form
#'
#' @return Genotype data from a collapsedVCF
extractGenotypeFromCollapsed <- function(collapsed_vcf) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Please install SummarizedExperiment: BiocManager::install('SummarizedExperiment')")
  }

  tryCatch({
    # Try getting GT from assays directly
    geno_data <- SummarizedExperiment::assay(collapsed_vcf, "GT")
    return(geno_data)
  }, error = function(e) {
    tryCatch({
      geno_data <- assays(collapsed_vcf)$GT
      return(geno_data)
    }, error = function(e) {
      if (!is.null(collapsed_vcf@assays$GT)) {
        return(collapsed_vcf@assays$GT)
      } else {
        stop("Could not find GT data in the CollapsedVCF object")
      }
    })
  })
}

#' Function to validate and process genotype data
#' @param gt_data Genotype data results from extractGenotypeFromCollapsed
#'
#' @return Converted genotype data in a matrix
convertGTtoNumeric <- function(gt_data) {
  # Convert to character matrix once
  gt_chars <- as.character(gt_data)

  # Create output matrix
  geno_mat <- matrix(NA_real_, nrow = nrow(gt_data), ncol = ncol(gt_data))
  dimnames(geno_mat) <- dimnames(gt_data)

  # Vectorized conversion using regex
  geno_mat[grep("0[|/]0", gt_chars)] <- 0
  geno_mat[grep("1[|/]1", gt_chars)] <- 2
  geno_mat[grep("(0[|/]1|1[|/]0)", gt_chars)] <- 1

  attr(geno_mat, "encoding") <- c("0" = "Reference homozygous",
                                  "1" = "Heterozygous",
                                  "2" = "Alternate homozygous")

  return(geno_mat)
}

#' Analyze population structure
#'
#' A function that runs Principal Component Analysis or PCA(can me modified for other analyses).
#' Output is ready for plotting various plot types.
#' Warning: When sample size is too small, you may receive the following Warning messages:
#' 1: In MASS::cov.trob(data[, vars]) : Probable convergence failure, or similar
#'
#' @param vcf_data A CollapsedVCF object
#' @param pop_metadata Population metadata data.frame
#' @param method One of "pca"
#' @param n_components A positive integer representing the Number of components that will be used in PCA
#' @param min_var Minimum variance threshold for filtering variants
#'
#' @return List containing analysis results for use in plotting
#'
#' @examples
#' # Grab filtered_data from previous load/filter functions
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
#' pca_results <- genRelateR::analyzePopulationStructure(
#'  filtered_data$vcf_data,
#'  filtered_data$pop_metadata,
#'  method = "pca"
#' )
#'
#' @export
analyzePopulationStructure <- function(vcf_data, pop_metadata,
                                       method = "pca", n_components = 2,
                                       min_var = 0.01) {
  # Extract and convert genotype data
  gt_data <- extractGenotypeFromCollapsed(vcf_data)
  geno_mat <- convertGTtoNumeric(gt_data)

  # Get sample information
  samples <- colnames(vcf_data)
  if (is.null(samples)) {
    stop("No sample names found in VCF data")
  }

  # Filter data
  filtered_geno <- geno_mat[, samples, drop = FALSE]

  # Filter variants by variance
  filtered_geno <- filterVariantsByVariance(filtered_geno, 0.01)
  filtered_geno <- filtered_geno[!rowSums(is.na(filtered_geno)),]

  filtered_pop <- pop_metadata$population[match(samples, pop_metadata$sample)]

  if (method == "pca") {
    # Handle potential scaling issues
    tryCatch({
      # Perform PCA with scaling
      pca_result <- prcomp(t(filtered_geno), scale = TRUE)

      # Create plot data
      plot_data <- data.frame(
        PC1 = pca_result$x[,1],
        PC2 = pca_result$x[,2],
        Population = filtered_pop
      )

      # Calculate variance
      percent_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

    }, error = function(e) {
      # If scaling fails, try without scaling
      message("Scaling failed, performing PCA without scaling")
      pca_result <- prcomp(t(filtered_geno), scale = FALSE)

      plot_data <- data.frame(
        PC1 = pca_result$x[,1],
        PC2 = pca_result$x[,2],
        Population = filtered_pop
      )

      percent_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
    })

    return(list(
      plot_data = plot_data,
      percent_var = percent_var,
      samples = samples
    ))
  }
  stop("Unsupported population structure analysis method.
       Feel free to suggest new methods to include!")
}
