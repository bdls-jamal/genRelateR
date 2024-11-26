#' Launch genRelateR Shiny Application
#'
#' @description
#' Launches the interactive Shiny web application for genomic relatedness analysis
#'
#' @return Launches the Shiny application
#'
#' @export
#' @importFrom shiny runApp
#'
#' @examples
#' \dontrun{
#' rungenRelateR()
#' }
rungenRelateR <- function() {
  # Find the Shiny app directory
  appDir <- system.file("shiny-scripts", package = "genRelateR")

  # Error handling if app directory not found
  if (appDir == "") {
    stop("Could not find Shiny app directory. Please reinstall genRelateR.", call. = FALSE)
  }

  # Launch the Shiny app
  shiny::runApp(appDir, display.mode = "normal")
}

