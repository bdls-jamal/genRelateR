#' Test Suite for Genetic Analysis Functions
#'
#' This script contains unit tests for all genetic analysis functions
#' using real VCF data from the 1000 Genomes Project.

library(testthat)
library(VariantAnnotation)
library(GenomicRanges)
library(data.table)
library(SNPRelate)
library(gdsfmt)
library(dplyr)
library(readr)
library(stringr)

test_that("loadGeneticData handles errors appropriately", {
  expect_error(loadGeneticData("nonexistent.vcf.gz"))
})

test_that("loadGeneticData works with 1000 Genomes Project data", {
  vcf_path <- "../../data/vcf/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

  skip_if_not(
    file.exists(vcf_path),
    "Skip test: 1000 Genomes Project VCF file not available"
  )

  # Test loading specific region
  chr1_region <- GenomicRanges::GRanges(
    seqnames = "1",
    ranges = IRanges(start = 20000000, end = 20001000)
  )

  # Suppress warning about index age as it's expected
  suppressWarnings({
    result <- loadGeneticData(
      vcf_path = vcf_path,
      regions = chr1_region
    )
  })

  expect_s4_class(result, "VCF")
  expect_true(nrow(result) > 0)

  # Test with sample filtering
  suppressWarnings({
    result_samples <- loadGeneticData(
      vcf_path = vcf_path,
      regions = chr1_region,
      samples = c("HG00096", "HG00097")  # Example sample IDs from 1000G
    )
  })

  expect_s4_class(result_samples, "VCF")
  expect_equal(ncol(result_samples), 2)
})


# Test data setup helper function
setup_test_data <- function() {
  vcf_path <- "../../data/vcf/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  pop_file <- "../../data/population_metadata.txt"
  rel_file <- "../../data/20140625_related_individuals.txt"

  skip_if_not(
    file.exists(vcf_path) && file.exists(pop_file) && file.exists(rel_file),
    "Skip test: Required test files not available"
  )

  # Load a small region for testing
  chr1_region <- GenomicRanges::GRanges(
    seqnames = "1",
    ranges = IRanges(start = 20000000, end = 20001000)
  )

  vcf_data <- suppressWarnings(
    loadGeneticData(vcf_path, regions = chr1_region)
  )

  return(list(vcf_data = vcf_data, pop_file = pop_file, rel_file = rel_file))
}

test_that("filterPopulation works with 1000 Genomes Project data", {
  test_data <- setup_test_data()

  # Test basic filtering
  filtered_vcf <- filterPopulation(
    test_data$vcf_data,
    test_data$pop_file,
    test_data$rel_file,
    population = c("CEU", "YRI")
  )

  expect_s4_class(filtered_vcf$vcf_data, "VCF")
  expect_true(ncol(filtered_vcf$vcf_data) < ncol(test_data$vcf_data))
})



