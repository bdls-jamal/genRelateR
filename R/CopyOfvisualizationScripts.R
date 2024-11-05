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
#' @import ggplot2 ggrepel
#' @export
plotPopulationPca <- function(analysis_results, title = NULL, ellipses = TRUE) {
  # Get number of populations
  n_pops <- length(unique(analysis_results$plot_data$Population))

  # Generate distinct colors using a simple rainbow scheme
  colors <- rainbow(n_pops)

  # Create the plot with improved aesthetics
  p <- ggplot(analysis_results$plot_data,
              aes(x = PC1, y = PC2, color = Population)) +
    geom_point(alpha = 0.7, size = 2) +  # Increased point size
    scale_color_manual(values = colors) +
    theme_bw(base_size = 12) +  # Increased base font size
    labs(
      x = sprintf("PC1 (%.1f%%)", analysis_results$percent_var[1]),
      y = sprintf("PC2 (%.1f%%)", analysis_results$percent_var[2]),
      title = title
    ) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, face = "bold"),
      legend.key.size = unit(1, "lines"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold")
    )

  # Add confidence ellipses if requested
  if (ellipses) {
    p <- p + stat_ellipse(level = 0.95, alpha = 0.3, linewidth = 0.7)  # Reduced alpha for better visibility
  }

  return(p)
}

#' Plot Relatedness Heatmap
#'
#' Creates a heatmap visualization of genetic relatedness between populations
#'
#' @param relatedness_results Results from computeRelatedness()
#' @param pop_metadata Population metadata dataframe
#' @param cluster Logical indicating whether to cluster populations (default: TRUE)
#' @param title Plot title (optional)
#' @return ggplot object
#' @import ggplot2 reshape2
#' @export
plotRelatednessHeatmap <- function(relatedness_results, pop_metadata,
                                   cluster = TRUE, title = NULL) {
  if (!is.matrix(relatedness_results$relatedness_matrix)) {
    stop("Input must contain a relatedness matrix")
  }

  # Melt the matrix for ggplot
  rel_melted <- reshape2::melt(relatedness_results$relatedness_matrix)
  names(rel_melted) <- c("Pop1", "Pop2", "Relatedness")

  # Create heatmap
  p <- ggplot(rel_melted, aes(x = Pop1, y = Pop2, fill = Relatedness)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(title = title)

  if (cluster) {
    # Perform hierarchical clustering
    hc <- hclust(dist(relatedness_results$relatedness_matrix))
    pop_order <- relatedness_results$relatedness_matrix[hc$order, hc$order]
    p <- p + scale_x_discrete(limits = rownames(pop_order)) +
      scale_y_discrete(limits = rownames(pop_order))
  }

  return(p)
}

#' Plot Ancestry Map
#'
#' Creates a map visualization of genetic ancestry components across geographic regions
#' @param analysis_results Analysis results from population structure analysis
#' @param pop_metadata Population metadata with geographic information
#' @param map_data Optional pre-loaded map data
#' @param title Plot title
#' @param individual Boolean to plot individual populations vs super populations
#' @return Either a ggplot object or plotly object depending on interactive parameter
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
             color = "gray70", fill = "gray90", size = 0.2) +
    xlim(-180, 180) +
    ylim(-90, 90)

  if (individual) {
    # Plot by individual populations
    p <- p +
      geom_point(data = plot_data,
                 aes(x = Longitude, y = Latitude, color = population),
                 alpha = 0.7, size = 3) +
      scale_color_viridis_d(name = "Population") +
      guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                                  title.position = "top",
                                  ncol = 2))  # Adjust number of columns as needed
  } else {
    # Plot by super-populations
    p <- p +
      geom_point(data = plot_data,
                 aes(x = Longitude, y = Latitude, color = super_pop),
                 alpha = 0.7, size = 3) +
      scale_color_viridis_d(name = "Super\nPopulation") +  # Added line break in title
      guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
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
#' Plot Migration Paths
#'
#' Visualizes inferred migration paths between populations based on genetic similarity
#'
#' @param relatedness_results Results from computeRelatedness()
#' @param pop_metadata Population metadata including geographic coordinates
#' @param threshold Minimum relatedness threshold for drawing paths (default: 0.1)
#' @param title Plot title (optional)
#' @return ggplot object
#' @import ggplot2 sf igraph
#' @export
plotMigrationPaths <- function(relatedness_results, pop_metadata,
                               threshold = 0.1, title = NULL) {
  # Validate inputs
  if (!all(c("Longitude", "Latitude") %in% colnames(pop_metadata))) {
    stop("Population metadata must include Longitude and Latitude columns")
  }

  # Create edges based on relatedness threshold
  rel_mat <- relatedness_results$relatedness_matrix
  edges <- which(rel_mat > threshold & upper.tri(rel_mat), arr.ind = TRUE)

  if (nrow(edges) == 0) {
    stop("No connections found above the threshold")
  }

  # Create edge data frame
  edge_data <- data.frame(
    Pop1 = rownames(rel_mat)[edges[,1]],
    Pop2 = colnames(rel_mat)[edges[,2]],
    weight = rel_mat[edges]
  )

  # Get geographic coordinates
  coords <- merge(edge_data, pop_metadata[, c("Population", "Longitude", "Latitude")],
                  by.x = "Pop1", by.y = "Population")
  coords <- merge(coords, pop_metadata[, c("Population", "Longitude", "Latitude")],
                  by.x = "Pop2", by.y = "Population", suffixes = c(".1", ".2"))

  # Create base world map
  world <- map_data("world")

  # Create migration path plot
  p <- ggplot() +
    geom_map(data = world, map = world,
             aes(long, lat, map_id = region),
             color = "gray70", fill = "gray90", size = 0.2) +
    geom_curve(data = coords,
               aes(x = Longitude.1, y = Latitude.1,
                   xend = Longitude.2, yend = Latitude.2,
                   alpha = weight),
               curvature = 0.2, color = "red") +
    geom_point(data = pop_metadata,
               aes(x = Longitude, y = Latitude, color = Population),
               size = 3) +
    scale_color_viridis_d() +
    scale_alpha_continuous(range = c(0.2, 0.8)) +
    theme_minimal() +
    labs(title = title) +
    coord_fixed(1.3)

  return(p)
}
