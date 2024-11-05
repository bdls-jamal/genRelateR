# Load required packages
library(ggplot2)
library(viridis)
library(reshape2)
library(maps)

source("R/setupPackages.R")
source("R/loadAndCleanData.R")

# Example data setup
# Load genetic data
vcf_file <- "data/vcf/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
genetic_data <- loadGeneticData(vcf_file)

# Create sample population metadata
pop_metadata <- "data/population_metadata.txt"

# Filter populations
filtered_data <- filterPopulation(genetic_data, pop_metadata, population = c("CEU", "YRI", "GBR"))

# Compute relatedness
relatedness_results <- computeRelatedness(
  filtered_data$vcf_data,
  filtered_data$pop_metadata
)

# Perform population structure analysis
pca_results <- analyzePopulationStructure(
  filtered_data$vcf_data,
  filtered_data$pop_metadata,
  method = "pca"
)

# Create visualizations
# PCA plot
pca_plot <- plotPopulationPca(
  analysis_results = pca_results,
  title = "Population Structure PCA",
  ellipses = TRUE
)
print(pca_plot)

# Relatedness heatmap
heatmap_plot <- plotRelatednessHeatmap(
  relatedness_results,
  title = "Population Relatedness"
)
print(heatmap_plot)

# Ancestry map
ancestry_map <- plotAncestryMap(
  pca_results,
  filtered_data$pop_metadata,
  title = "Global Population Distribution"
)
print(ancestry_map)

# Migration paths
migration_plot <- plotMigrationPaths(
  relatedness_results,
  filtered_data$pop_metadata,
  threshold = 0.2,
  title = "Inferred Migration Paths"
)
print(migration_plot)
