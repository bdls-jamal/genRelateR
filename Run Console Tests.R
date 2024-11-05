# Example data setup
# Load genetic data
vcf_file <- "data/vcf/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
chr22_region <- GRanges(
  seqnames = "22",
  ranges = IRanges(start = 19999000, end = 20001000)
)
genetic_data <- loadGeneticData(vcf_file, regions = chr22_region)

# Create sample population metadata
pop_metadata <- "data/population_metadata.txt"

# Filter populations
populations <- c(
  "CHB", "JPT", "CHS", "CDX", "KHV",
  "CEU", "TSI", "GBR", "FIN", "IBS",
  "YRI", "LWK", "GWD", "MSL", "ESN",
  "ASW", "ACB", "MXL", "PUR", "CLM",
  "PEL", "GIH", "PJL", "BEB", "STU", "ITU"
)

filtered_data <- filterPopulation(genetic_data, pop_metadata, populations)

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
  filtered_data$pop_metadata,
  title = "Population Structure PCA",
  ellipses = TRUE,
  super_pop = TRUE
)
print(pca_plot)

# Ancestry map
ancestry_map <- plotAncestryMap(
  pca_results,
  filtered_data$pop_metadata,
  title = "Global Population Distribution"
)
print(ancestry_map)

