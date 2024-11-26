# Load required libraries
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(DT)
library(plotly)
library(genRelateR)

# UI Definition
ui <- dashboardPage(
  # Dashboard Header
  dashboardHeader(title = "genRelateR: Genomic Relatedness Explorer"),

  # Dashboard Sidebar
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("Population Selection", tabName = "population", icon = icon("users")),
      menuItem("Analysis", tabName = "analysis", icon = icon("chart-bar")),
      menuItem("Visualization", tabName = "visualization", icon = icon("image"))
    )
  ),

  # Dashboard Body
  dashboardBody(
    tabItems(
      # Data Upload Tab
      tabItem(tabName = "upload",
              fluidRow(
                box(
                  title = "Upload Genetic Data", status = "primary", solidHeader = TRUE,
                  fileInput("vcf_file", "Upload VCF File",
                            accept = c(".vcf", ".vcf.gz")),
                  helpText("Supported file types: .vcf, .vcf.gz")
                ),
                box(
                  title = "Upload Metadata", status = "primary", solidHeader = TRUE,
                  fileInput("metadata_file", "Upload Population Metadata",
                            accept = c(".txt")),
                  helpText("Population metadata file")
                )
              ),
              fluidRow(
                box(
                  title = "Demo Data", status = "info", solidHeader = TRUE,
                  actionButton("use_demo_data", "Use Demo Data")
                )
              )
      ),

      # Population Selection Tab
      tabItem(tabName = "population",
              fluidRow(
                box(
                  title = "Population Filter", status = "primary", solidHeader = TRUE, width = 12,
                  checkboxGroupInput("selected_populations",
                                     "Select Populations",
                                     choices = c(
                                       "East Asian" = "EAS",
                                       "European" = "EUR",
                                       "African" = "AFR",
                                       "American" = "AMR",
                                       "South Asian" = "SAS"
                                     ),
                                     selected = NULL,
                                     inline = TRUE)
                )
              )
      ),

      # Analysis Tab
      tabItem(tabName = "analysis",
              fluidRow(
                box(
                  title = "Analysis Options", status = "primary", solidHeader = TRUE, width = 12,
                  selectInput("analysis_type", "Select Analysis",
                              choices = c(
                                "Population Structure PCA" = "pca",
                                "Ancestry Distribution Map" = "map",
                                "Relatedness Dashboard" = "dashboard"
                              ))
                )
              )
      ),

      # Visualization Tab
      tabItem(tabName = "visualization",
              fluidRow(
                box(
                  title = "Analysis Results", status = "primary",
                  plotlyOutput("analysis_plot", height = "600px"),
                  width = 12
                )
              )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  # Reactive values to store data
  rv <- reactiveValues(
    vcf_data = NULL,
    metadata = NULL,
    filtered_data = NULL,
    analysis_results = NULL
  )

  # Demo Data Handler
  observeEvent(input$use_demo_data, {
    # Use package's demo data paths
    demo_vcf <- system.file("extdata", "demo_vcf", "demo_data.vcf.gz", package = "genRelateR")
    demo_metadata <- system.file("extdata", "demo_metadata", "demo_metadata.txt", package = "genRelateR")

    # Load demo data
    rv$vcf_data <- loadGeneticData(demo_vcf)
    rv$metadata <- read.table(demo_metadata, header = TRUE)

    # Notify user
    showNotification("Demo data loaded successfully!", type = "message")
  })

  # VCF File Upload
  observeEvent(input$vcf_file, {
    req(input$vcf_file)
    rv$vcf_data <- loadGeneticData(input$vcf_file$datapath)
    showNotification("VCF file loaded successfully!", type = "message")
  })

  # Metadata File Upload
  observeEvent(input$metadata_file, {
    req(input$metadata_file)
    rv$metadata <- read.table(input$metadata_file$datapath, header = TRUE)
    showNotification("Metadata file loaded successfully!", type = "message")
  })

  # Population Filtering
  observe({
    req(rv$vcf_data, rv$metadata, input$selected_populations)
    rv$filtered_data <- filterPopulation(
      rv$vcf_data,
      rv$metadata,
      populations = input$selected_populations
    )
  })

  # Analysis Execution
  observe({
    req(rv$filtered_data, input$analysis_type)

    rv$analysis_results <- switch(input$analysis_type,
                                  "pca" = analyzePopulationStructure(
                                    rv$filtered_data$vcf_data,
                                    rv$filtered_data$pop_metadata,
                                    method = "pca"
                                  ),
                                  "map" = plotAncestryMap(
                                    rv$vcf_data,
                                    rv$metadata
                                  ),
                                  "dashboard" = createRelatednessDashboard(
                                    rv$filtered_data$vcf_data,
                                    rv$filtered_data$pop_metadata
                                  )
    )
  })

  # Visualization Output
  output$analysis_plot <- renderPlotly({
    req(rv$analysis_results, input$analysis_type)

    plot <- switch(input$analysis_type,
                   "pca" = plotPopulationPca(
                     rv$analysis_results,
                     rv$filtered_data$pop_metadata,
                     title = "Population Structure PCA"
                   ),
                   "map" = rv$analysis_results,  # Already a plot from plotAncestryMap
                   "dashboard" = rv$analysis_results  # Assuming this returns a plot/dashboard
    )

    # Convert ggplot to plotly for interactivity
    ggplotly(plot)
  })
}

# Run the Shiny App
shinyApp(ui, server)
