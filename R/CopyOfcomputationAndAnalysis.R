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
#' @param collapsed_vcf vcf data in the collapsed form
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

# Function to validate and process genotype data
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

#' Calculate IBS Matrix
#' @param geno_mat Numeric genotype matrix
#' @return IBS similarity matrix
calculateIBSMatrixOptimized <- function(geno_mat) {
  if (!is.matrix(geno_mat) || !is.numeric(geno_mat)) {
    stop("Input must be a numeric matrix")
  }

  n_samples <- ncol(geno_mat)

  # Convert to sparse matrix for efficient operations
  sparse_geno <- Matrix::Matrix(geno_mat, sparse = TRUE)

  # Initialize result matrix
  ibs_mat <- Matrix::Matrix(0, n_samples, n_samples)
  dimnames(ibs_mat) <- list(colnames(geno_mat), colnames(geno_mat))

  # Calculate in parallel if possible
  if (requireNamespace("parallel", quietly = TRUE)) {
    # Determine number of cores (leave one for system)
    n_cores <- max(1, parallel::detectCores() - 1)

    # Create cluster
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))

    # Export necessary data
    parallel::clusterExport(cl, c("sparse_geno"), envir = environment())

    # Split work into chunks
    chunks <- split(1:n_samples, cut(1:n_samples, n_cores))

    # Process chunks in parallel
    results <- parallel::parLapply(cl, chunks, function(idx_range) {
      result_chunk <- Matrix::Matrix(0, n_samples, length(idx_range))

      for (i in seq_along(idx_range)) {
        idx <- idx_range[i]
        col_i <- sparse_geno[, idx]
        valid_i <- !is.na(col_i)

        for (j in idx:n_samples) {
          col_j <- sparse_geno[, j]
          valid_j <- !is.na(col_j)

          # Calculate shared valid positions
          valid_both <- valid_i & valid_j
          shared <- sum(valid_both)

          if (shared > 0) {
            # Calculate matches only for valid positions
            matches <- sum(col_i[valid_both] == col_j[valid_both])
            result_chunk[j, i] <- matches / shared
          }
        }
      }
      return(result_chunk)
    })

    # Combine results
    for (i in seq_along(chunks)) {
      idx_range <- chunks[[i]]
      chunk_result <- results[[i]]
      ibs_mat[, idx_range] <- chunk_result
    }
  } else {
    # Non-parallel fallback with optimizations
    for (i in 1:n_samples) {
      col_i <- sparse_geno[, i]
      valid_i <- !is.na(col_i)

      for (j in i:n_samples) {
        col_j <- sparse_geno[, j]
        valid_j <- !is.na(col_j)

        valid_both <- valid_i & valid_j
        shared <- sum(valid_both)

        if (shared > 0) {
          matches <- sum(col_i[valid_both] == col_j[valid_both])
          ibs_mat[i, j] <- ibs_mat[j, i] <- matches / shared
        }
      }
    }
  }

  # Make symmetric
  ibs_mat <- as.matrix(ibs_mat + t(ibs_mat) - Diagonal(n_samples) * diag(ibs_mat))
  return(ibs_mat)
}

# Calculate FST Matrix
#' @param geno_mat Numeric genotype matrix
#' @param pop_codes table of 3 letter population codes mapped to names
#' @return FST similarity matrix
calculateFSTMatrixOptimized <- function(geno_mat, pop_codes) {
  if (!is.matrix(geno_mat) || !is.numeric(geno_mat)) {
    stop("Input must be a numeric matrix")
  }

  populations <- unique(pop_codes)
  n_pops <- length(populations)

  # Pre-calculate population indices and frequencies
  pop_indices <- lapply(populations, function(pop) which(pop_codes == pop))

  # Pre-calculate frequencies in parallel if possible
  if (requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
    on.exit(parallel::stopCluster(cl))

    parallel::clusterExport(cl, c("geno_mat", "pop_indices"), envir = environment())

    pop_freqs <- parallel::parLapply(cl, pop_indices, function(idx) {
      rowMeans(geno_mat[, idx, drop = FALSE], na.rm = TRUE) / 2
    })
  } else {
    pop_freqs <- lapply(pop_indices, function(idx) {
      rowMeans(geno_mat[, idx, drop = FALSE], na.rm = TRUE) / 2
    })
  }

  # Initialize matrix
  fst_mat <- matrix(0, n_pops, n_pops)
  dimnames(fst_mat) <- list(populations, populations)

  # Calculate FST values
  for (i in 1:n_pops) {
    freq_i <- pop_freqs[[i]]
    n_i <- length(pop_indices[[i]])

    for (j in i:n_pops) {
      freq_j <- pop_freqs[[j]]
      n_j <- length(pop_indices[[j]])

      # Vectorized calculations
      freq_total <- (freq_i * n_i + freq_j * n_j) / (n_i + n_j)

      # Vectorized heterozygosity calculations
      Hs_i <- mean(2 * freq_i * (1 - freq_i), na.rm = TRUE)
      Hs_j <- mean(2 * freq_j * (1 - freq_j), na.rm = TRUE)
      Hs <- (Hs_i + Hs_j) / 2

      Ht <- mean(2 * freq_total * (1 - freq_total), na.rm = TRUE)

      # Calculate FST
      fst <- if (Ht > 0) max((Ht - Hs) / Ht, 0) else 0
      fst_mat[i, j] <- fst_mat[j, i] <- fst
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
  if (!is(vcf_data, "CollapsedVCF")) {
    stop("vcf_data must be a CollapsedVCF object")
  }

  # Extract and convert genotype data efficiently
  gt_data <- extractGenotypeFromCollapsed(vcf_data)
  geno_mat <- convertGTtoNumeric(gt_data)

  # Get sample information
  samples <- colnames(vcf_data)
  if (is.null(samples)) {
    stop("No sample names found in VCF data")
  }

  # Match population data efficiently
  pop_codes <- pop_metadata$pop[match(samples, pop_metadata$sample)]

  # Calculate relatedness matrix using optimized functions
  relatedness_matrix <- if (method == "ibs") {
    calculateIBSMatrixOptimized(geno_mat)
  } else if (method == "fst") {
    calculateFSTMatrixOptimized(geno_mat, pop_codes)
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
#' @param method One of "pca"
#' @param n_components Number of components
#' @param min_var Minimum variance threshold for filtering variants
#' @return List containing analysis results
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
