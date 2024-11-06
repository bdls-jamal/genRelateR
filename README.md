
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

This package also contains sample 1000 genome project data found in the
data folder. You may use these for initial analyses and then download
other files from the project later.

Refer to package vignettes for more details. An overview of the package
is illustrated below.

## Contributions

The author of this package is Kobi Jamal Schmalenberg. The author wrote
all of the functions in this package. Generative AI tools were used in
all functions for proofing, optimization suggestions, and to help easily
replicate styling conventions.

## References

- Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S
  Language. Wadsworth & Brooks/Cole.
- Mardia, K. V., J. T. Kent, and J. M. Bibby (1979) Multivariate
  Analysis, London: Academic Press.
- Susan Fairley, Ernesto Lowy-Gallego, Emily Perry, Paul Flicek, The
  International Genome Sample Resource (IGSR) collection of open human
  genomic variation resources, Nucleic Acids Research, Volume 48, Issue
  D1, 08 January 2020, Pages D941–D947,
  <https://doi.org/10.1093/nar/gkz836>
- The 1000 Genomes Project Consortium. A global reference for human
  genetic variation. Nature 526, 68–74 (2015).
  <https://doi.org/10.1038/nature15393>
- Venables, W. N. and B. D. Ripley (2002) Modern Applied Statistics with
  S, Springer-Verlag.

### Packages

- Barrett T, Dowle M, Srinivasan A, Gorecki J, Chirico M, Hocking T,
  Schwendinger B (2024). *data.table: Extension of `data.frame`*. R
  package version 1.16.0,
  <https://CRAN.R-project.org/package=data.table>.
- Becker OScbRA, Minka ARWRvbRBEbTP, Deckmyn. A (2023). *maps: Draw
  Geographical Maps*. R package version 3.4.2,
  <https://CRAN.R-project.org/package=maps>.
- Csardi G, Nepusz T (2006). “The igraph software package for complex
  network research.” *InterJournal*, *Complex Systems*, 1695.
  <https://igraph.org>.
- C. Sievert. Interactive Web-Based Data Visualization with R, plotly,
  and shiny. Chapman and Hall/CRC Florida, 2020.
- Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013)
  Software for Computing and Annotating Genomic Ranges. PLoS Comput Biol
  9(8): e1003118. <doi:10.1371/journal.pcbi.1003118>
- Neuwirth E (2022). *RColorBrewer: ColorBrewer Palettes*. R package
  version 1.1-3, <https://CRAN.R-project.org/package=RColorBrewer>.
- Obenchain V, Lawrence M, Carey V, Gogarten S, Shannon P, Morgan M
  (2014). “VariantAnnotation: a Bioconductor package for exploration and
  annotation of genetic variants.” *Bioinformatics*, *30*(14),
  2076-2078. <doi:10.1093/bioinformatics/btu168>
  <https://doi.org/10.1093/bioinformatics/btu168>.
- Simon Garnier, Noam Ross, Robert Rudis, Antônio P. Camargo, Marco
  Sciaini, and Cédric Scherer (2024). viridis(Lite) -
  Colorblind-Friendly Color Maps for R. viridis package version 0.6.5.
- Slowikowski K (2024). *ggrepel: Automatically Position Non-Overlapping
  Text Labels with ‘ggplot2’*. R package version 0.9.6,
  <https://CRAN.R-project.org/package=ggrepel>.
- Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R,
  Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller
  E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V,
  Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to
  the tidyverse.” *Journal of Open Source Software*, *4*(43), 1686.
  <doi:10.21105/joss.01686> <https://doi.org/10.21105/joss.01686>.
- Wickham H, François R, Henry L, Müller K, Vaughan D (2023). *dplyr: A
  Grammar of Data Manipulation*. R package version 1.1.4,
  <https://CRAN.R-project.org/package=dplyr>.
- Wickham H. ggplot2: Elegant Graphics for Data Analysis.
  Springer-Verlag New York, 2016.
- Wickham H, Hester J, Bryan J (2024). *readr: Read Rectangular Text
  Data*. R package version 2.1.5,
  <https://CRAN.R-project.org/package=readr>.
- Wickham H, Vaughan D, Girlich M (2024). *tidyr: Tidy Messy Data*. R
  package version 1.3.1, <https://CRAN.R-project.org/package=tidyr>.
- Wickham H (2023). *stringr: Simple, Consistent Wrappers for Common
  String Operations*. R package version 1.5.1,
  <https://CRAN.R-project.org/package=stringr>.

## Acknowledgements

This package was developed as part of an assessment for 2024 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
Canada. `genRelateR` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/bdls-jamal/genRelateR/issues).

## Example Usage

``` r
Navigate to R/setupPackages.R and run this entire file to install and load necessary packages
You may download the necessary files(vcf.gz/tbi files and metadata.txt file) directly from the github repository or from the 1000 genomes project.

# Load genetic data from compressed vcf file(i.e. data from 1000 genomes project)
# Ensure this data has an associated tbi file in the same folder, or generate it using Tabix
# Example: vcf_file <- "inst/extdata/vcf/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
vcf_file <- "path/to/your/vcf.gz/data"

# Define a region(range) for analysis(required for very large data files)
chr1_region <- GRanges(
  seqnames = "1",
  ranges = IRanges(start = 20000000, end = 20001000)
)

# Call loadGeneticData function with specified region and vcf file
genetic_data <- genRelateR::loadGeneticData(vcf_file, regions = chr1_region)

# Create population metadata object from metadata file
# File is provided for 1000 genome project use, for other data you must use with your own data
# Example: pop_metadata <- "inst/extdata/population_metadata.txt"
pop_metadata <- "path/to/your/metadata.txt/file"

# Choose populations to filter by(all are selected in this example)
populations <- c(
  "CHB", "JPT", "CHS", "CDX", "KHV",
  "CEU", "TSI", "GBR", "FIN", "IBS",
  "YRI", "LWK", "GWD", "MSL", "ESN",
  "ASW", "ACB", "MXL", "PUR", "CLM",
  "PEL", "GIH", "PJL", "BEB", "STU", "ITU"
)

# Call filterPopulation with loaded data, metadata, and specified populations
filtered_data <- genRelateR::filterPopulation(genetic_data, pop_metadata, populations)

# Perform population structure analysis by calling analyzePopulationStructure
pca_results <- genRelateR::analyzePopulationStructure(
  filtered_data$vcf_data,
  filtered_data$pop_metadata,
  method = "pca"
)

# Create PCA data cluster visualization
# Call plotPopulationPca then print to view
pca_plot <- genRelateR::plotPopulationPca(
  analysis_results = pca_results,
  filtered_data$pop_metadata,
  title = "Population Structure PCA",
  ellipses = TRUE,
  super_pop = TRUE
)
print(pca_plot)

# Create PCA data map visualization
# Call plotAncestryMap then print to view
ancestry_map <- genRelateR::plotAncestryMap(
  pca_results,
  filtered_data$pop_metadata,
  title = "Global Population Distribution"
)
print(ancestry_map)
```
