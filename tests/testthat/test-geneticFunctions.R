#' Test Suite for Genetic Analysis Functions
#'
#' This script contains unit tests for all genetic analysis functions
#' using real VCF data from the 1000 Genomes Project.

library(testthat)
library(VariantAnnotation)
library(GenomicRanges)
library(data.table)

# # Create a simple test VCF file
# create_test_vcf <- function() {
#   vcf_content <- '##fileformat=VCFv4.2
# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
# 1\t100\trs1\tA\tT\t100\tPASS\t.\tGT\t0/0\t0/1
# 1\t200\trs2\tC\tG\t100\tPASS\t.\tGT\t0/1\t1/1'
#
#   # Create temporary VCF file
#   tmp_file <- tempfile(fileext = ".vcf")
#   writeLines(vcf_content, tmp_file)
#   return(tmp_file)
# }
#
# test_that("loadGeneticData loads VCF data correctly", {
#   # Create test file
#   vcf_file <- create_test_vcf()
#
#   # Test basic loading
#   result <- loadGeneticData(vcf_file)
#
#   # Check that we get a VCF object back
#   expect_s4_class(result, "VCF")
#
#   # Check dimensions
#   expect_equal(nrow(result), 2)  # 2 variants
#   expect_equal(ncol(result), 2)  # 2 samples
#
#   # Clean up
#   unlink(vcf_file)
# })


test_that("loadGeneticData handles errors appropriately", {
  expect_error(loadGeneticData("nonexistent.vcf.gz"))
})

test_that("loadGeneticData works with 1000 Genomes Project data", {
  vcf_path <- "../../data/vcf/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

  skip_if_not(
    file.exists(vcf_path),
    "Skip test: 1000 Genomes Project VCF file not available"
  )

  # Test loading specific region
  chr22_region <- GRanges(
    seqnames = "22",
    ranges = IRanges(start = 20000000, end = 20001000)
  )

  # Suppress warning about index age as it's expected
  suppressWarnings({
    result <- loadGeneticData(
      vcf_path = vcf_path,
      regions = chr22_region
    )
  })

  expect_s4_class(result, "VCF")
  expect_true(nrow(result) > 0)

  # Test with sample filtering
  suppressWarnings({
    result_samples <- loadGeneticData(
      vcf_path = vcf_path,
      regions = chr22_region,
      samples = c("HG00096", "HG00097")  # Example sample IDs from 1000G
    )
  })

  expect_s4_class(result_samples, "VCF")
  expect_equal(ncol(result_samples), 2)
})
