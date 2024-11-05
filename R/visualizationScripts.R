#' Plot Population Structure PCA
#'
#' Creates a scatter plot of the first two principal components from population structure analysis
#' with populations colored distinctly and confidence ellipses.
#'
#' @param analysis_results Results from analyzePopulationStructure() with method="pca"
#' @param title Plot title (optional)
#' @param ellipses Logical indicating whether to draw confidence ellipses (default: TRUE)
#' @param labels Logical indicating whether to show sample labels (default: FALSE)
#' @return ggplot object
#' @examples
#' # Grab pca_results from previous analysis functions
#' pca_results <- genRelateR::analyzePopulationStructure(
#'  filtered_data$vcf_data,
#'  filtered_data$pop_metadata,
#'  method = "pca"
#' )
#'
#' pca_plot <- genRelateR::plotPopulationPca(
#'   analysis_results = pca_results,
#'   filtered_data$pop_metadata,
#'   title = "Population Structure PCA",
#'   ellipses = TRUE,
#'   super_pop = TRUE
#'   )
#' print(pca_plot)
#'
#' @export
#' @import ggplot2 ggrepel
plotPopulationPca <- function(analysis_results, pop_metadata, title = NULL, ellipses = TRUE, super_pop = FALSE) {
  # Convert row names to a column in analysis_results$plot_data for merging
  analysis_results$plot_data <- data.frame(
    sample = row.names(analysis_results$plot_data),
    analysis_results$plot_data
  )

  # Merge analysis results with population metadata
  plot_data <- merge(
    analysis_results$plot_data,
    pop_metadata[, c("sample", "Longitude", "Latitude", "population", "super_pop")],
    by = "sample",
    all.x = TRUE
  )

  # Determine whether to color by Population or Super Population
  color_var <- if (super_pop) "super_pop" else "population"

  # Get number of unique populations/super populations
  unique_groups <- unique(plot_data[[color_var]])
  n_groups <- length(unique_groups)

  # Generate distinct colors using a simple rainbow scheme
  colors <- rainbow(n_groups)
  names(colors) <- unique_groups

  # Create the plot with improved aesthetics
  p <- ggplot(plot_data,
              aes(x = PC1, y = PC2, color = .data[[color_var]])) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(
      values = colors,
      name = if (super_pop) "Super Population" else "Population"
    ) +
    theme_bw(base_size = 12) +
    labs(
      x = sprintf("PC1 (%.1f%%)", analysis_results$percent_var[1]),
      y = sprintf("PC2 (%.1f%%)", analysis_results$percent_var[2]),
      title = title
    ) +
    theme(
      legend.position = "right",
      legend.text = element_text(linewidth = 8),
      legend.title = element_text(linewidth = 10, face = "bold"),
      legend.key.size = unit(1, "lines"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      plot.title = element_text(face = "bold", linewidth = 14),
      axis.title = element_text(face = "bold")
    )

  # Add confidence ellipses if requested
  if (ellipses) {
    p <- p + stat_ellipse(aes(group = .data[[color_var]]),
                          level = 0.95,
                          alpha = 0.3,
                          linewidth = 0.7)
  }

  return(p)
}

#' Plot Ancestry Map
#'
#' Creates a map visualization of genetic ancestry components across geographic regions
#'
#' @param analysis_results Analysis results from population structure analysis
#' @param pop_metadata Population metadata with geographic information
#' @param map_data Optional pre-loaded map data
#' @param title Plot title (optional)
#' @param individual Boolean to plot individual populations vs super populations
#' @return Either a ggplot object or plotly object depending on interactive parameter
#' @examples
#' # Grab pca_results from previous analysis functions
#' pca_results <- genRelateR::analyzePopulationStructure(
#'  filtered_data$vcf_data,
#'  filtered_data$pop_metadata,
#'  method = "pca"
#' )
#'
#' ancestry_map <- genRelateR::plotAncestryMap(
#'  pca_results,
#'  filtered_data$pop_metadata,
#'  title = "Global Population Distribution"
#' )
#' print(ancestry_map)
#' @export
#' @import ggplot2 maps
plotAncestryMap <- function(analysis_results, pop_metadata, map_data = NULL,
                            title = NULL, individual = FALSE) {
  # Validate inputs
  if (!all(c("Longitude", "Latitude", "super_pop") %in% colnames(pop_metadata))) {
    stop("Population metadata must include Longitude, Latitude, and super_pop columns")
  }

  # Convert row names to a column in analysis_results$plot_data for merging
  analysis_results$plot_data <- data.frame(
    sample = row.names(analysis_results$plot_data),
    analysis_results$plot_data
  )

  # Merge analysis results with geographic data
  plot_data <- merge(
    analysis_results$plot_data,
    pop_metadata[, c("sample", "Longitude", "Latitude", "population", "super_pop")],
    by.x = "sample",
    by.y = "sample",
    all.x = TRUE
  )

  # Create base map
  world <- if (is.null(map_data)) {
    map_data("world")
  } else {
    map_data
  }

  # Create base plot
  p <- ggplot() +
    geom_map(data = world, map = world,
             aes(map_id = region),
             color = "gray70", fill = "gray90", linewidth = 0.2) +
    xlim(-180, 180) +
    ylim(-90, 90)

  if (individual) {
    # Plot by individual populations
    p <- p +
      geom_point(data = plot_data,
                 aes(x = Longitude, y = Latitude, color = population),
                 alpha = 0.7, size = 3) +
      scale_color_viridis_d(name = "Population") +
      guides(color = guide_legend(override.aes = list(linewidth = 4, alpha = 1),
                                  title.position = "top",
                                  ncol = 2))  # Adjust number of columns as needed
  } else {
    # Plot by super-populations
    p <- p +
      geom_point(data = plot_data,
                 aes(x = Longitude, y = Latitude, color = super_pop),
                 alpha = 0.7, size = 3) +
      scale_color_viridis_d(name = "Super\nPopulation") +  # Added line break in title
      guides(color = guide_legend(override.aes = list(linewidth = 4, alpha = 1),
                                  title.position = "top",
                                  ncol = 2))  # Adjust number of columns as needed
  }

  # Add common theme elements with modified legend
  p <- p +
    theme_minimal() +
    labs(title = title) +
    coord_fixed(1.3)  +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      aspect.ratio = 0.5
    )

  return(p)
}
