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

#' Extract genotype matrix from VCF data
#' @param vcf_data A CollapsedVCF object
#' @return A numeric matrix of genotypes
extractGenotypeMatrix <- function(vcf_data) {
  # Extract the GT (genotype) field from the VCF
  geno <- geno(vcf_data)$GT

  if (is.null(geno)) {
    stop("No GT (genotype) field found in VCF data")
  }

  # Convert genotypes to numeric matrix
  geno_num <- matrix(NA, nrow=nrow(geno), ncol=ncol(geno))
  colnames(geno_num) <- colnames(geno)
  rownames(geno_num) <- rownames(geno)

  for(i in 1:nrow(geno)) {
    for(j in 1:ncol(geno)) {
      if (!is.na(geno[i,j])) {
        alleles <- strsplit(geno[i,j], "/|\\|")[[1]]
        geno_num[i,j] <- sum(as.numeric(alleles))
      }
    }
  }
  return(geno_num)
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
#' @param rel_data Relatedness information data.frame
#' @param method One of "ibs" or "fst"
#' @return List containing relatedness matrix and sample information
computeRelatedness <- function(vcf_data, pop_metadata, rel_data, method = "ibs") {
  # Validate inputs
  if (!is(vcf_data, "CollapsedVCF")) {
    stop("vcf_data must be a CollapsedVCF object")
  }

  # Extract and convert genotype data
  geno_mat <- extractGenotypeMatrix(vcf_data)

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

  # Get excluded samples
  excluded_samples <- rel_data$Sample
  filtered_samples <- setdiff(samples, excluded_samples)

  return(list(
    relatedness_matrix = relatedness_matrix,
    filtered_samples = filtered_samples,
    excluded_samples = excluded_samples
  ))
}

#' Analyze population structure
#' @param vcf_data A CollapsedVCF object
#' @param pop_metadata Population metadata data.frame
#' @param rel_data Relatedness information data.frame
#' @param method One of "pca" or "admixture"
#' @param n_components Number of components
#' @param min_var Minimum variance threshold for filtering variants
#' @return List containing analysis results
analyzePopulationStructure <- function(vcf_data, pop_metadata, rel_data,
                                       method = "pca", n_components = 2,
                                       min_var = 1e-10) {
  # Extract and convert genotype data
  geno_mat <- extractGenotypeMatrix(vcf_data)

  # Get sample information
  samples <- colnames(vcf_data)
  if (is.null(samples)) {
    stop("No sample names found in VCF data")
  }

  # Get excluded samples
  excluded_samples <- rel_data$Sample
  filtered_samples <- setdiff(samples, excluded_samples)

  # Filter data
  filtered_geno <- geno_mat[, filtered_samples, drop = FALSE]

  # Filter variants by variance
  filtered_geno <- filterVariantsByVariance(filtered_geno, min_var)

  filtered_pop <- pop_metadata$pop[match(filtered_samples, pop_metadata$sample)]

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
      filtered_samples = filtered_samples,
      excluded_samples = excluded_samples
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
      filtered_samples = filtered_samples,
      excluded_samples = excluded_samples
    ))
  }

  stop("Unsupported population structure analysis method")
}
