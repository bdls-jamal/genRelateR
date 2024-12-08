% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadAndCleanData.R
\name{loadGeneticData}
\alias{loadGeneticData}
\title{Load and Format Genetic Data}
\usage{
loadGeneticData(vcf_path, regions = NULL, samples = NULL)
}
\arguments{
\item{vcf_path}{Path to VCF file}

\item{regions}{GRanges object specifying regions to load. If NULL, automatically selects
a default 1Mb region based on the chromosome in the filename.}

\item{samples}{Vector of sample IDs to include. If NULL, includes all samples.}
}
\value{
A VCF object containing the genetic data
}
\description{
This function loads genetic data from the specified file path and formats it for analysis.
By default, it loads a manageable subset of the chromosome to prevent memory issues.
If the default region is still too large, please specify a smaller region using the
regions parameter (e.g., GRanges("1", IRanges(1e6, 2e6)) for chromosome 1, 1-2Mb).
}
\details{
Warning: If you are receiving a tbi file is older than the vcf file error,
Run the following commands to recreate the tbi file:
library(VariantAnnotation)
indexTabix("path/to/vcf.gz/file", format = "vcf")
}
\examples{
# Example with default region:
vcf_file <- "../inst/ext/data/vcf/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
genetic_data <- genRelateR::loadGeneticData(vcf_file)

# Example with custom region:
regions <- GenomicRanges::GRanges("1", IRanges(1e6, 2e6))
genetic_data <- genRelateR::loadGeneticData(vcf_file, regions = regions)

}
