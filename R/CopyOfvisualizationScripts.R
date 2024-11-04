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
plotPopulationPca <- function(analysis_results, title = NULL,
                              ellipses = TRUE, labels = FALSE) {
  if (!all(c("plot_data", "percent_var") %in% names(analysis_results))) {
    stop("Input must be results from analyzePopulationStructure() with method='pca'")
  }

  # Create base plot
  p <- ggplot(analysis_results$plot_data,
              aes(x = PC1, y = PC2, color = Population)) +
    geom_point(alpha = 0.7) +
    scale_color_viridis_d() +
    theme_bw() +
    labs(
      x = sprintf("PC1 (%.1f%%)", analysis_results$percent_var[1]),
      y = sprintf("PC2 (%.1f%%)", analysis_results$percent_var[2]),
      title = title
    )

  # Add confidence ellipses if requested
  if (ellipses) {
    p <- p + stat_ellipse(level = 0.95, alpha = 0.2)
  }

  # Add sample labels if requested
  if (labels) {
    p <- p + geom_text_repel(aes(label = rownames(analysis_results$plot_data)),
                             size = 3, max.overlaps = 20)
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
#'
#' @param analysis_results Results from analyzePopulationStructure() with method="admixture"
#' @param pop_metadata Population metadata including geographic coordinates
#' @param map_data Optional map data for custom boundaries
#' @param title Plot title (optional)
#' @return ggplot object
#' @import ggplot2 sf
#' @export
plotAncestryMap <- function(analysis_results, pop_metadata, map_data = NULL,
                            title = NULL) {
  # Validate inputs
  if (!all(c("Longitude", "Latitude") %in% colnames(pop_metadata))) {
    stop("Population metadata must include Longitude and Latitude columns")
  }

  # Merge analysis results with geographic data
  plot_data <- merge(
    analysis_results$plot_data,
    pop_metadata[, c("sample", "Longitude", "Latitude", "Population")],
    by.x = "rownames(analysis_results$plot_data)",
    by.y = "sample"
  )

  # Create base map
  world <- if (is.null(map_data)) {
    map_data("world")
  } else {
    map_data
  }

  # Create map plot
  p <- ggplot() +
    geom_map(data = world, map = world,
             aes(long, lat, map_id = region),
             color = "gray70", fill = "gray90", size = 0.2) +
    geom_point(data = plot_data,
               aes(x = Longitude, y = Latitude, color = Population),
               alpha = 0.7, size = 3) +
    scale_color_viridis_d() +
    theme_minimal() +
    labs(title = title) +
    coord_fixed(1.3)

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
