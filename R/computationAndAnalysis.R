#' Helper function to validate data frame columns
validate_columns <- function(df, required_cols, df_name) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in %s: %s",
                 df_name, paste(missing_cols, collapse = ", ")))
  }
}

# Helper functions needed for the main functions to work

#' Filter variants by variance
#' @param geno_mat Numeric genotype matrix
#' @param min_var Minimum variance threshold
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

# Function to extract genotype data from CollapsedVCF
extractGenotypeFromCollapsed <- function(collapsed_vcf) {
  # First try accessing through assays
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Please install SummarizedExperiment: BiocManager::install('SummarizedExperiment')")
  }

  tryCatch({
    # Try getting GT from assays directly
    geno_data <- SummarizedExperiment::assay(collapsed_vcf, "GT")
    return(geno_data)
  }, error = function(e) {
    # Alternative method using assays slot
    tryCatch({
      geno_data <- assays(collapsed_vcf)$GT
      return(geno_data)
    }, error = function(e) {
      # Last resort - try accessing the internal structure
      if (!is.null(collapsed_vcf@assays$GT)) {
        return(collapsed_vcf@assays$GT)
      } else {
        stop("Could not find GT data in the CollapsedVCF object")
      }
    })
  })
}

# Function to validate and process genotype data
convertGTtoNumeric <- function(gt_data) {
  # Create output matrix
  geno_mat <- matrix(NA, nrow = nrow(gt_data), ncol = ncol(gt_data))

  # Copy dimension names if they exist
  if (!is.null(dimnames(gt_data))) {
    dimnames(geno_mat) <- dimnames(gt_data)
  }

  # Conversion with input validation
  for(i in 1:nrow(gt_data)) {
    for(j in 1:ncol(gt_data)) {
      gt <- as.character(gt_data[i,j])

      # Handle different possible formats
      if (grepl("[0-9]\\|[0-9]", gt) || grepl("[0-9]/[0-9]", gt)) {
        # Split the genotype into alleles
        alleles <- as.numeric(strsplit(gt, "[|/]")[[1]])

        if (length(alleles) == 2) {
          if (all(alleles == 0)) geno_mat[i,j] <- 0      # 0|0 -> 0
          else if (all(alleles == 1)) geno_mat[i,j] <- 2 # 1|1 -> 2
          else if (sum(alleles) == 1) geno_mat[i,j] <- 1 # 0|1 or 1|0 -> 1
          else geno_mat[i,j] <- NA                       # Invalid format
        }
      } else {
        geno_mat[i,j] <- NA  # Invalid format
      }
    }
  }

  # Add informative attributes
  attr(geno_mat, "encoding") <- c("0" = "Reference homozygous",
                                  "1" = "Heterozygous",
                                  "2" = "Alternate homozygous")

  return(geno_mat)
}


#' Calculate IBS Matrix
#' @param geno_mat Numeric genotype matrix
#' @return IBS similarity matrix
calculateIBSMatrix <- function(geno_mat) {
  if (!is.matrix(geno_mat)) {
    stop("Input must be a matrix")
  }

  if (!is.numeric(geno_mat)) {
    stop("Input matrix must be numeric")
  }

  n_samples <- ncol(geno_mat)
  ibs_mat <- matrix(0, n_samples, n_samples)
  colnames(ibs_mat) <- colnames(geno_mat)
  rownames(ibs_mat) <- colnames(geno_mat)

  for(i in 1:n_samples) {
    for(j in i:n_samples) {
      shared <- sum(!is.na(geno_mat[,i]) & !is.na(geno_mat[,j]))
      if(shared > 0) {
        ibs <- sum(geno_mat[,i] == geno_mat[,j], na.rm=TRUE) / shared
        ibs_mat[i,j] <- ibs_mat[j,i] <- ibs
      }
    }
  }
  return(ibs_mat)
}

calculateFstMatrix <- function(geno_mat, pop_codes) {
  # Simple implementation of Weir and Cockerham's Fst
  unique_pops <- unique(pop_codes)
  n_pops <- length(unique_pops)
  fst_mat <- matrix(0, n_pops, n_pops)
  rownames(fst_mat) <- colnames(fst_mat) <- unique_pops

  # Convert genotypes to allele frequencies per population
  for(i in 1:n_pops) {
    for(j in i:n_pops) {
      pop1_samples <- pop_codes == unique_pops[i]
      pop2_samples <- pop_codes == unique_pops[j]

      # Calculate allele frequencies
      freq1 <- rowMeans(geno_mat[,pop1_samples] == "1|1", na.rm=TRUE)
      freq2 <- rowMeans(geno_mat[,pop2_samples] == "1|1", na.rm=TRUE)

      # Simple Fst calculation
      fst <- mean(abs(freq1 - freq2), na.rm=TRUE)
      fst_mat[i,j] <- fst_mat[j,i] <- fst
    }
  }
  return(fst_mat)
}

#' Compute genetic relatedness
#' @param vcf_data A CollapsedVCF object
#' @param pop_metadata Population metadata data.frame
#' @param method One of "ibs" or "fst"
#' @return List containing relatedness matrix and sample information
computeRelatedness <- function(vcf_data, pop_metadata, method = "ibs") {
  # Validate inputs
  if (!is(vcf_data, "CollapsedVCF")) {
    stop("vcf_data must be a CollapsedVCF object")
  }

  # Extract and convert genotype data
  gt_data <- extractGenotypeFromCollapsed(vcf_data)
  geno_mat <- convertGTtoNumeric(gt_data)

  # Get sample information
  samples <- colnames(vcf_data)
  if (is.null(samples)) {
    stop("No sample names found in VCF data")
  }

  # Match population data
  pop_codes <- pop_metadata$pop[match(samples, pop_metadata$sample)]

  # Calculate relatedness matrix
  relatedness_matrix <- if (method == "ibs") {
    calculateIBSMatrix(geno_mat)
  } else if (method == "fst") {
    calculateFstMatrix(geno_mat, pop_codes)
  } else {
    stop("Unsupported relatedness method")
  }

  return(list(
    relatedness_matrix = relatedness_matrix,
    samples = samples
  ))
}

#' Analyze population structure
#' @param vcf_data A CollapsedVCF object
#' @param pop_metadata Population metadata data.frame
#' @param method One of "pca" or "admixture"
#' @param n_components Number of components
#' @param min_var Minimum variance threshold for filtering variants
#' @return List containing analysis results
analyzePopulationStructure <- function(vcf_data, pop_metadata,
                                       method = "pca", n_components = 2,
                                       min_var = 1e-10) {
  # Extract and convert genotype data
  geno_mat <- extractGenotypeMatrix(vcf_data)

  # Get sample information
  samples <- colnames(vcf_data)
  if (is.null(samples)) {
    stop("No sample names found in VCF data")
  }

  # Filter data
  filtered_geno <- geno_mat[, samples, drop = FALSE]

  # Filter variants by variance
  filtered_geno <- filterVariantsByVariance(filtered_geno, min_var)

  filtered_pop <- pop_metadata$pop[match(samples, pop_metadata$sample)]

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

      # Calculate variance explained
      percent_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

    }, error = function(e) {
      # If scaling fails, try without scaling
      message("Warning: Scaling failed, performing PCA without scaling")
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

  } else if (method == "admixture") {
    # Perform MDS
    dist_mat <- dist(t(filtered_geno))
    mds_result <- cmdscale(dist_mat, k = n_components)

    plot_data <- data.frame(
      MDS1 = mds_result[,1],
      MDS2 = mds_result[,2],
      Population = filtered_pop
    )

    return(list(
      plot_data = plot_data,
      samples = samples
    ))
  }

  stop("Unsupported population structure analysis method")
}
