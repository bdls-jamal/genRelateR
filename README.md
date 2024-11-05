
<!-- README.md is generated from README.Rmd. Please edit that file -->

# genRelateR

## Description

genRelateR is a tool for analyzing global human genomic relatedness
across populations. It enhances workflows in bioinformatics by providing
exploratory analysis tools for genomic data. The `genRelateR` package
was developed using `R version 4.4.1 (2024-06-14 ucrt)`,
`Platform: x86_64-w64-mingw32/x64` and
`Running under: Windows 10 x64 (build 19045)`.

## Installation

To install the latest version of the package:

``` r
install.packages("devtools")
library(devtools)
devtools::install_github("bdls-jamal/genRelateR", build_vignettes = TRUE)
library(genRelateR) 
```

To run the Shiny app(not yet implemented):

``` r
rungenRelateR() # not yet implemented
```

## Overview

``` r
ls("package:genRelateR")
data(package = "genRelateR")
browseVignettes("genRelateR")
```

`genRelateR` contains 6 functions.

1.  ***loadGeneticData()*** loads and formats genetic vcf data from
    various sources
2.  ***filterPopulation()*** filters data based on population,
    geographic region, or genetic criteria
3.  ***analyzePopulationStructure()*** performs PCA analysis on
    population data
4.  ***plotPopulationPca()*** plots population structure using PCA
5.  ***plotAncestryMap()*** maps the distribution of genetic ancestry
    components
6.  ***createRelatednessDashboard()*** generates an interactive shiny
    dashboard for exploring relatedness data

Add overview image later: Refer to package vignettes for more details.
An overview of the package is illustrated below.

## Contributions

Provide a paragraph clearly indicating the name of the author of the
package and contributions from the author. Outline contributions from
other packages/sources for each function. Outline contributions from
generative AI tool(s) for each function. Include how the tools were used
and how the results from AI tools were incorporated. Remember your
individual contributions to the package are important. E.g.,

The author of the package is Kobi Jamal Schmalenberg.

## References

## Acknowledgements

This package was developed as part of an assessment for 2024 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
Canada. `genRelateR` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/bdls-jamal/genRelateR/issues).

## Example Usage

``` r
# Load genetic data from compressed vcf file(i.e. data from 1000 genomes project)
# Ensure this data has an associated tbi file in the same folder, or generate it using Tabix
vcf_file <- "data/vcf/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

# Define a region(range) for analysis(required for very large data files)
chr22_region <- GRanges(
  seqnames = "22",
  ranges = IRanges(start = 19999000, end = 20001000)
)

# Call loadGeneticData function with specified region and vcf file
genetic_data <- loadGeneticData(vcf_file, regions = chr22_region)

# Create population metadata object from metadata file
# File is provided for 1000 genome project use, for other data you must use with your own data
pop_metadata <- "data/population_metadata.txt"

# Choose populations to filter by(all are selected in this example)
populations <- c(
  "CHB", "JPT", "CHS", "CDX", "KHV",
  "CEU", "TSI", "GBR", "FIN", "IBS",
  "YRI", "LWK", "GWD", "MSL", "ESN",
  "ASW", "ACB", "MXL", "PUR", "CLM",
  "PEL", "GIH", "PJL", "BEB", "STU", "ITU"
)

# Call filterPopulation with loaded data, metadata, and specified populations
filtered_data <- filterPopulation(genetic_data, pop_metadata, populations)

# Perform population structure analysis by calling analyzePopulationStructure
pca_results <- analyzePopulationStructure(
  filtered_data$vcf_data,
  filtered_data$pop_metadata,
  method = "pca"
)

# Create PCA data cluster visualization
# Call plotPopulationPca then print to view
pca_plot <- plotPopulationPca(
  analysis_results = pca_results,
  filtered_data$pop_metadata,
  title = "Population Structure PCA",
  ellipses = TRUE,
  super_pop = TRUE
)
print(pca_plot)

# Create PCA data map visualization
# Call plotAncestryMap then print to view
ancestry_map <- plotAncestryMap(
  pca_results,
  filtered_data$pop_metadata,
  title = "Global Population Distribution"
)
print(ancestry_map)
```
