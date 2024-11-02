#' Compute Genetic Relatedness
#'
#' @param genetic_data Data frame of genetic data
#' @param method Relatedness calculation method
#' @return Matrix of relatedness coefficients
#' @import SNPRelate
computeRelatedness <- function(genetic_data, method = "IBS") {
  # Convert genotype data to numeric format (0,1,2 coding)
  geno_matrix <- apply(genetic_data[,-(1:4)], 2, function(x) {
    ifelse(x == "0/0", 0,
           ifelse(x == "0/1" | x == "1/0", 1,
                  ifelse(x == "1/1", 2, NA)))
  })

  # Compute relatedness matrix
  if(method == "IBS") {
    n_samples <- ncol(geno_matrix)
    rel_matrix <- matrix(0, n_samples, n_samples)

    for(i in 1:n_samples) {
      for(j in i:n_samples) {
        # Calculate IBS (Identity By State)
        shared_alleles <- sum(geno_matrix[,i] == geno_matrix[,j], na.rm = TRUE)
        total_sites <- sum(!is.na(geno_matrix[,i]) & !is.na(geno_matrix[,j]))
        rel_matrix[i,j] <- rel_matrix[j,i] <- shared_alleles/total_sites
      }
    }

    colnames(rel_matrix) <- rownames(rel_matrix) <- colnames(genetic_data[,-(1:4)])
    return(rel_matrix)
  }
}

#' Analyze Population Structure
#'
#' @param genetic_data Data frame of genetic data
#' @param n_components Number of principal components
#' @return List containing PCA results and variance explained
#' @import stats
analyzePopulationStructure <- function(genetic_data, n_components = 3) {
  # Convert genotypes to numeric
  geno_matrix <- apply(genetic_data[,-(1:4)], 2, function(x) {
    ifelse(x == "0/0", 0,
           ifelse(x == "0/1" | x == "1/0", 1,
                  ifelse(x == "1/1", 2, NA)))
  })

  # Handle missing values
  geno_matrix[is.na(geno_matrix)] <- colMeans(geno_matrix, na.rm = TRUE)

  # Perform PCA
  pca_result <- prcomp(t(geno_matrix), scale. = TRUE)

  # Calculate variance explained
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

  return(list(
    pca = pca_result,
    variance_explained = var_explained,
    scores = pca_result$x[,1:n_components]
  ))
}
