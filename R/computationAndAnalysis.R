#' Helper function to validate data frame columns
validate_columns <- function(df, required_cols, df_name) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in %s: %s",
                 df_name, paste(missing_cols, collapse = ", ")))
  }
}

# Helper functions needed for the main functions to work
calculateIBSMatrix <- function(geno_mat) {
  # Convert genotypes to numeric (0,1,2)
  geno_num <- matrix(NA, nrow=nrow(geno_mat), ncol=ncol(geno_mat))
  for(i in 1:nrow(geno_mat)) {
    for(j in 1:ncol(geno_mat)) {
      gt <- geno_mat[i,j]
      geno_num[i,j] <- sum(as.numeric(strsplit(gt, "/|\\|")[[1]]))
    }
  }

  # Calculate IBS matrix
  n_samples <- ncol(geno_num)
  ibs_mat <- matrix(0, n_samples, n_samples)
  for(i in 1:n_samples) {
    for(j in i:n_samples) {
      shared <- sum(!is.na(geno_num[,i]) & !is.na(geno_num[,j]))
      if(shared > 0) {
        ibs <- sum(geno_num[,i] == geno_num[,j], na.rm=TRUE) / shared
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

genotypeToNumeric <- function(geno_mat) {
  # Convert genotypes to numeric values
  geno_num <- matrix(NA, nrow=nrow(geno_mat), ncol=ncol(geno_mat))
  for(i in 1:nrow(geno_mat)) {
    for(j in 1:ncol(geno_mat)) {
      gt <- geno_mat[i,j]
      geno_num[i,j] <- sum(as.numeric(strsplit(gt, "/|\\|")[[1]]))
    }
  }
  return(geno_num)
}

#' Compute Genetic Relatedness
#'
#' @param vcf_data A CollapsedVCF object already loaded and filtered
#' @param pop_metadata Population metadata data.frame
#' @param rel_data Relatedness information data.frame
#' @param method One of "ibs" or "fst"
#' @param n_threads Number of threads for parallel processing
#' @return List containing relatedness matrix and filtered samples
#' @import VariantAnnotation
computeRelatedness <- function(vcf_data, pop_metadata, rel_data, method = "ibs") {
  # Validate inputs
  validate_columns(rel_data, c("Sample", "Population", "Gender", "Reason for exclusion"), "relationship data")
  validate_columns(pop_metadata, c("sample", "pop", "super_pop", "gender"), "population metadata")

  # Get genotype matrix
  geno_mat <- geno(vcf_data)

  # Get sample names and populations
  samples <- colnames(vcf_data)
  pop_codes <- pop_metadata$pop[match(samples, pop_metadata$sample)]

  # Calculate relatedness matrix based on method
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

  list(
    relatedness_matrix = relatedness_matrix,
    filtered_samples = filtered_samples,
    excluded_samples = excluded_samples
  )
}

#' Analyze Population Structure
#'
#' @param vcf_data A CollapsedVCF object already loaded and filtered
#' @param pop_metadata Population metadata data.frame
#' @param rel_data Relatedness information data.frame
#' @param method One of "pca" or "admixture"
#' @param n_components Number of components for PCA or K for admixture
#' @return List containing analysis results and plotting data
#' @import VariantAnnotation
analyzePopulationStructure <- function(vcf_data, pop_metadata, rel_data,
                                       method = "pca", n_components = 2) {
  # Validate inputs
  validate_columns(rel_data, c("Sample", "Population", "Gender", "Reason for exclusion"), "relationship data")
  validate_columns(pop_metadata, c("sample", "pop", "super_pop", "gender"), "population metadata")

  # Convert genotypes to numeric matrix
  geno_num <- genotypeToNumeric(geno(vcf_data))

  # Get excluded samples
  excluded_samples <- rel_data$Sample
  filtered_samples <- setdiff(colnames(vcf_data), excluded_samples)

  # Filter data
  filtered_geno <- geno_num[, filtered_samples]
  filtered_pop <- pop_metadata$pop[match(filtered_samples, pop_metadata$sample)]

  if (method == "pca") {
    # Perform PCA
    pca_result <- prcomp(t(filtered_geno), scale = TRUE)

    # Create plot data
    plot_data <- data.frame(
      PC1 = pca_result$x[,1],
      PC2 = pca_result$x[,2],
      Population = filtered_pop
    )

    # Calculate percent variance
    percent_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

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
