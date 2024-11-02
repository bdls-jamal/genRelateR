# Load testthat library
library(testthat)

# Load the function to be tested
source("../../R/loadAndCleanData.R")

# Test case: Ensure loadGeneticData loads the VCF file correctly
test_that("loadGeneticData loads VCF file correctly", {
  # Sample VCF file for testing (adjust the path as necessary)
  vcf_file <- "../../data/vcf/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

  # Define a ScanVcfParam for a specific range in chromosome 1
  param <- ScanVcfParam(which = GRanges("1", IRanges(1e6, 1.5e6)))  # Adjust range as needed

  # Call the loadGeneticData function
  genetic_data <- loadGeneticData(file_path = vcf_file, genome = "b37", param = param)

  # Check if the returned object is a data frame
  expect_s3_class(genetic_data, "data.frame")
})


