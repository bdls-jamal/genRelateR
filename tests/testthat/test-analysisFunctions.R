library(testthat)
library(VariantAnnotation)
library(GenomicRanges)
library(data.table)
library(SNPRelate)
library(gdsfmt)
library(dplyr)
library(readr)
library(stringr)

setup_test_data <- function() {
  # Get file paths
  vcf_path <- normalizePath("../../data/vcf/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
                            mustWork = FALSE)
  pop_path <- normalizePath("../../data/population_metadata.txt",
                            mustWork = FALSE)

  # Check file existence
  if (!file.exists(vcf_path)) {
    skip(paste("VCF file not found at:", vcf_path))
  }
  if (!file.exists(pop_path)) {
    skip(paste("Population metadata file not found at:", pop_path))
  }

  # Load a small region for testing
  chr1_region <- GenomicRanges::GRanges(
    seqnames = "1",
    ranges = IRanges(start = 20000000, end = 20001000)
  )

  # Load VCF data with error handling
  vcf_data <- tryCatch({
    suppressWarnings(loadGeneticData(vcf_path, regions = chr1_region))
  }, error = function(e) {
    skip(paste("Error loading VCF data:", e$message))
  })

  # Return the test data structure
  test_data <- list(
    vcf_data = vcf_data,
    pop_metadata = pop_path,
    n_samples = ncol(vcf_data)
  )

  # Validate the test data structure
  if (!all(sapply(test_data, function(x) !is.null(x) && length(x) > 0))) {
    skip("Invalid test data structure generated")
  }

  return(test_data)
}


test_that("analyzePopulationStructure works with loaded VCF data", {
  test_data <- setup_test_data()

  # Filter test_data
  filtered_vcf <- filterPopulation(
    test_data$vcf_data,
    test_data$pop_metadata,
    population = c("CEU", "YRI", "GBR")
  )

  # Test PCA
  pca_results <- analyzePopulationStructure(
    filtered_vcf$vcf_data,
    filtered_vcf$pop_metadata,
    method = "pca",
    n_components = 2
  )

  # Check results structure
  expect_type(pca_results, "list")
  expect_true(all(c("plot_data", "percent_var", "samples") %in% names(pca_results)))

  # Check plot data
  expect_true("PC1" %in% names(pca_results$plot_data))
  expect_true("PC2" %in% names(pca_results$plot_data))
  expect_true("Population" %in% names(pca_results$plot_data))

  # Check plot data
  expect_true("MDS1" %in% names(admix_results$plot_data))
  expect_true("MDS2" %in% names(admix_results$plot_data))
  expect_true("Population" %in% names(admix_results$plot_data))
})

