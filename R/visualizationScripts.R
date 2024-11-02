#' Plot Population PCA
#'
#' @param pca_results Output from analyzePopulationStructure
#' @param population_labels Vector of population labels
#' @return ggplot object
#' @import ggplot2
plotPopulationPca <- function(pca_results, population_labels) {
  # Create data frame for plotting
  plot_data <- data.frame(
    PC1 = pca_results$scores[,1],
    PC2 = pca_results$scores[,2],
    Population = population_labels
  )

  # Create PCA plot
  ggplot(plot_data, aes(x = PC1, y = PC2, color = Population)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_bw() +
    xlab(paste0("PC1 (", round(pca_results$variance_explained[1] * 100, 1), "%)")) +
    ylab(paste0("PC2 (", round(pca_results$variance_explained[2] * 100, 1), "%)")) +
    theme(legend.position = "right")
}

#' Plot Relatedness Heatmap
#'
#' @param relatedness_matrix Output from computeRelatedness
#' @param population_labels Vector of population labels
#' @return ggplot object
#' @import ggplot2 reshape2
plotRelatednessHeatmap <- function(relatedness_matrix, population_labels) {
  # Prepare data for plotting
  melted_matrix <- reshape2::melt(relatedness_matrix)
  melted_matrix$Population1 <- population_labels[melted_matrix$Var1]
  melted_matrix$Population2 <- population_labels[melted_matrix$Var2]

  # Create heatmap
  ggplot(melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0.5, limits = c(0, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "", fill = "Relatedness")
}

#' Plot Ancestry Map
#'
#' @param genetic_data Data frame with ancestry components
#' @param geographic_data Data frame with latitude/longitude
#' @return ggplot object
#' @import ggplot2 maps
plotAncestryMap <- function(genetic_data, geographic_data) {
  # Create world map
  world <- map_data("world")

  # Create base map
  ggplot() +
    geom_map(data = world, map = world,
             aes(x = long, y = lat, map_id = region),
             fill = "white", color = "grey50", size = 0.1) +
    geom_point(data = geographic_data,
               aes(x = longitude, y = latitude, color = ancestry_component),
               size = 3, alpha = 0.7) +
    theme_minimal() +
    coord_fixed(1.3) +
    scale_color_viridis_c() +
    labs(color = "Ancestry Component")
}

#' Plot Migration Paths
#'
#' @param migration_data Data frame with source/destination coordinates
#' @param flow_strength Vector of migration strengths
#' @return ggplot object
#' @import ggplot2 maps geosphere
plotMigrationPaths <- function(migration_data, flow_strength) {
  # Create world map
  world <- map_data("world")

  # Calculate great circle paths
  paths <- lapply(1:nrow(migration_data), function(i) {
    gcIntermediate(
      c(migration_data$source_lon[i], migration_data$source_lat[i]),
      c(migration_data$dest_lon[i], migration_data$dest_lat[i]),
      n = 100, addStartEnd = TRUE
    )
  })

  # Create base map with paths
  ggplot() +
    geom_map(data = world, map = world,
             aes(x = long, y = lat, map_id = region),
             fill = "white", color = "grey50", size = 0.1) +
    geom_path(data = do.call(rbind, paths),
              aes(x = lon, y = lat, group = group,
                  alpha = flow_strength[group],
                  size = flow_strength[group]),
              color = "red") +
    theme_minimal() +
    coord_fixed(1.3) +
    scale_alpha_continuous(range = c(0.2, 0.8)) +
    scale_size_continuous(range = c(0.5, 2)) +
    labs(alpha = "Migration Strength", size = "Migration Strength")
}
