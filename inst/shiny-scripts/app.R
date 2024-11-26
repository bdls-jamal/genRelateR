library(shiny)
library(shinydashboard)
library(shinyjs)
library(genRelateR)
library(VariantAnnotation)
library(GenomicRanges)
library(Rsamtools)

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Genomic Data Analyzer"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("Population Selection", tabName = "populations", icon = icon("users")),
      menuItem("Analysis", tabName = "analysis", icon = icon("chart-bar"))
    )
  ),

  dashboardBody(
    shinyjs::useShinyjs(),
    tags$head(
      tags$style(HTML("
        .shiny-input-container {
          max-width: 100%;
        }
        .large-file-warning {
          color: orange;
          font-weight: bold;
        }
      "))
    ),
    tabItems(
      # Data Upload Tab
      tabItem(tabName = "upload",
              fluidRow(
                box(
                  title = "VCF File Upload", status = "primary", solidHeader = TRUE, width = 12,
                  fileInput("vcf_file", "Upload VCF File",
                            accept = c(".vcf", ".vcf.gz"),
                            multiple = FALSE),
                  tags$div(
                    class = "large-file-warning",
                    "Note: For very large files, use region-based processing."
                  )
                ),
                box(
                  title = "Metadata Upload", status = "primary", solidHeader = TRUE, width = 12,
                  fileInput("metadata_file", "Upload Population Metadata",
                            accept = c(".txt", ".csv"))
                )
              )
      ),

      # Population Selection Tab
      tabItem(tabName = "populations",
              fluidRow(
                box(
                  title = "Chromosome and Region", status = "primary", solidHeader = TRUE,
                  selectInput("chromosome", "Select Chromosome",
                              choices = c(1:22, "X", "Y", "MT"),
                              selected = "1"),
                  sliderInput("region_range", "Select Genomic Region",
                              min = 1, max = 250000000,
                              value = c(20000000, 20001000), step = 1000000)
                ),
                box(
                  title = "Population Selection", status = "primary", solidHeader = TRUE,
                  checkboxGroupInput("populations", "Select Populations",
                                     choices = c(
                                       "East Asian" = "CHB,JPT,CHS,CDX,KHV",
                                       "European" = "CEU,TSI,GBR,FIN,IBS",
                                       "African" = "YRI,LWK,GWD,MSL,ESN,ASW,ACB",
                                       "American" = "MXL,PUR,CLM,PEL",
                                       "South Asian" = "GIH,PJL,BEB,STU,ITU"
                                     ),
                                     selected = c("CHB,JPT,CHS,CDX,KHV",
                                                  "CEU,TSI,GBR,FIN,IBS"))
                )
              )
      ),

      # Analysis Tab
      tabItem(tabName = "analysis",
              fluidRow(
                box(
                  title = "Analysis Options", status = "primary", solidHeader = TRUE,
                  selectInput("analysis_type", "Select Analysis",
                              choices = c(
                                "Population Structure PCA" = "pca",
                                "Ancestry Map" = "ancestry_map"
                              )),
                  actionButton("run_analysis", "Run Analysis", icon = icon("play"))
                ),
                box(
                  title = "Analysis Results", status = "primary", solidHeader = TRUE,
                  plotOutput("analysis_plot")
                )
              )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  # Increase file upload limit
  options(shiny.maxRequestSize = 10 * 1024^3)  # 10 GB
  # Reactive Values
  rv <- reactiveValues(
    vcf_file = NULL,
    metadata_file = NULL,
    genetic_data = NULL,
    filtered_data = NULL,
    pca_results = NULL
  )

  # Setup Genetic Packages
  observe({
    genRelateR::setupGeneticPackages()
  })

  # VCF File Upload
  observeEvent(input$vcf_file, {
    req(input$vcf_file)
    rv$vcf_file <- input$vcf_file$datapath

    # Check and potentially create tabix index
    tbi_file <- paste0(rv$vcf_file, ".tbi")

    tryCatch({
      # Attempt to create tabix index if it doesn't exist
      if (!file.exists(tbi_file)) {
        message("Creating tabix index...")

        # Use Rsamtools to create tabix index
        indexTabix(rv$vcf_file, format = "vcf")
      }

      # Load the genetic data with the updated function
      rv$genetic_data <- genRelateR::loadGeneticData(
        rv$vcf_file,
        regions = GenomicRanges::GRanges(seqnames = "1", ranges = IRanges(1, 1e6))
      )

    }, error = function(e) {
      showNotification(
        paste("Error processing VCF file:", e$message),
        type = "error",
        duration = 10
      )
      rv$genetic_data <- NULL
    })
  })

  # Metadata File Upload
  observeEvent(input$metadata_file, {
    req(input$metadata_file)
    rv$metadata_file <- input$metadata_file$datapath
  })

  # Load Genetic Data
  observe({
    req(rv$vcf_file, input$chromosome, input$region_range)

    # Create GRanges object for region
    region <- GRanges(
      seqnames = as.character(input$chromosome),
      ranges = IRanges(
        start = input$region_range[1],
        end = input$region_range[2]
      )
    )

    # Load Genetic Data
    rv$genetic_data <- tryCatch({
      genRelateR::loadGeneticData(rv$vcf_file, regions = region)
    }, error = function(e) {
      showNotification(
        paste("Error loading genetic data:", e$message),
        type = "error"
      )
      NULL
    })
  })

  # Filter Populations
  observe({
    req(rv$genetic_data, rv$metadata_file, input$populations)

    # Flatten and split populations
    populations <- unlist(strsplit(input$populations, ","))

    # Filter Populations
    rv$filtered_data <- tryCatch({
      genRelateR::filterPopulation(
        rv$genetic_data,
        rv$metadata_file,
        populations
      )
    }, error = function(e) {
      showNotification(
        paste("Error filtering populations:", e$message),
        type = "error"
      )
      NULL
    })
  })

  # Run Analysis
  observeEvent(input$run_analysis, {
    req(rv$filtered_data)

    # Perform Analysis
    rv$pca_results <- tryCatch({
      genRelateR::analyzePopulationStructure(
        rv$filtered_data$vcf_data,
        rv$filtered_data$pop_metadata,
        method = "pca"
      )
    }, error = function(e) {
      showNotification(
        paste("Error in population structure analysis:", e$message),
        type = "error"
      )
      NULL
    })
  })

  # Plot Results
  output$analysis_plot <- renderPlot({
    req(rv$pca_results, rv$filtered_data)

    if (input$analysis_type == "pca") {
      genRelateR::plotPopulationPca(
        analysis_results = rv$pca_results,
        rv$filtered_data$pop_metadata,
        title = "Population Structure PCA",
        ellipses = TRUE,
        super_pop = TRUE
      )
    } else if (input$analysis_type == "ancestry_map") {
      genRelateR::plotAncestryMap(
        rv$pca_results,
        rv$filtered_data$pop_metadata,
        title = "Global Population Distribution"
      )
    }
  })
}

# Run the application
shinyApp(ui, server)
