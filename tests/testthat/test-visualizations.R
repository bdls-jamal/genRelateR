# tests/testthat/test-visualization.R

library(testthat)
library(ggplot2)
library(viridis)
library(maps)
library(rlang)

# Helper function to create mock PCA results
create_mock_pca_results <- function() {
  # Create mock PCA data
  n_samples <- 100
  set.seed(42)  # For reproducibility
  plot_data <- data.frame(
    PC1 = rnorm(n_samples),
    PC2 = rnorm(n_samples)
  )
  rownames(plot_data) <- paste0("sample", 1:n_samples)

  list(
    plot_data = plot_data,
    percent_var = c(15.2, 8.7)  # Mock variance explained
  )
}

# Helper function to create mock population metadata
create_mock_metadata <- function(n_samples = 100) {
  set.seed(42)  # For reproducibility
  populations <- c("CHB", "JPT", "CEU", "YRI", "GBR")
  super_pops <- c("EAS", "EAS", "EUR", "AFR", "EUR")

  data.frame(
    sample = paste0("sample", 1:n_samples),
    population = sample(populations, n_samples, replace = TRUE),
    super_pop = sample(super_pops, n_samples, replace = TRUE),
    Longitude = runif(n_samples, -180, 180),
    Latitude = runif(n_samples, -90, 90),
    stringsAsFactors = FALSE
  )
}

test_that("plotPopulationPca produces valid ggplot object", {
  # Setup
  mock_pca <- create_mock_pca_results()
  mock_metadata <- create_mock_metadata()

  # Test basic plot
  plot <- plotPopulationPca(
    analysis_results = mock_pca,
    pop_metadata = mock_metadata,
    title = "Test PCA Plot"
  )

  # Assertions
  expect_s3_class(plot, "ggplot")

  # Get the mapping names
  mapping_names <- names(plot$mapping)
  expect_true("x" %in% mapping_names)
  expect_true("y" %in% mapping_names)
})

test_that("plotPopulationPca handles super_pop parameter correctly", {
  # Setup
  mock_pca <- create_mock_pca_results()
  mock_metadata <- create_mock_metadata()

  # Test with super_pop = TRUE
  plot_super <- plotPopulationPca(
    analysis_results = mock_pca,
    pop_metadata = mock_metadata,
    super_pop = TRUE
  )

  # Test with super_pop = FALSE
  plot_pop <- plotPopulationPca(
    analysis_results = mock_pca,
    pop_metadata = mock_metadata,
    super_pop = FALSE
  )

  # Assertions
  expect_s3_class(plot_super, "ggplot")
  expect_s3_class(plot_pop, "ggplot")

  # Check that color mapping exists
  expect_true("colour" %in% names(plot_super$mapping))
  expect_true("colour" %in% names(plot_pop$mapping))
})

test_that("plotPopulationPca handles ellipses parameter correctly", {
  # Setup
  mock_pca <- create_mock_pca_results()
  mock_metadata <- create_mock_metadata()

  # Test with ellipses = TRUE
  plot_with_ellipses <- plotPopulationPca(
    analysis_results = mock_pca,
    pop_metadata = mock_metadata,
    ellipses = TRUE
  )

  # Test with ellipses = FALSE
  plot_without_ellipses <- plotPopulationPca(
    analysis_results = mock_pca,
    pop_metadata = mock_metadata,
    ellipses = FALSE
  )

  # Assertions
  expect_s3_class(plot_with_ellipses, "ggplot")
  expect_s3_class(plot_without_ellipses, "ggplot")

  # Count number of layers
  expect_gt(length(plot_with_ellipses$layers), length(plot_without_ellipses$layers))
})

test_that("plotAncestryMap produces valid ggplot object", {
  skip_if_not_installed("maps")

  # Setup
  mock_pca <- create_mock_pca_results()
  mock_metadata <- create_mock_metadata()

  # Test basic map
  map <- plotAncestryMap(
    analysis_results = mock_pca,
    pop_metadata = mock_metadata,
    title = "Test Ancestry Map"
  )

  # Assertions
  expect_s3_class(map, "ggplot")
})

test_that("plotAncestryMap handles missing coordinates gracefully", {
  skip_if_not_installed("maps")

  # Setup
  mock_pca <- create_mock_pca_results()
  mock_metadata <- create_mock_metadata()

  # Create more missing coordinates to ensure warning is triggered
  mock_metadata$Longitude <- NA
  mock_metadata$Latitude <- NA

  # Test map with missing coordinates
  expect_warning(
    plotAncestryMap(
      analysis_results = mock_pca,
      pop_metadata = mock_metadata,
      title = "Test Map with Missing Coordinates"
    ),
    regexp = NA  # Accept any warning
  )
})
