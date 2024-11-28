library(shiny)
library(shinydashboard)
library(shinyjs)
library(genRelateR)
library(GenomicRanges)
library(shinyWidgets)

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
        .population-group {
          font-weight: bold;
          margin-top: 10px;
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
                            accept = c(".vcf.gz"),
                            multiple = FALSE),
                  fileInput("tbi_file", "Upload Tabix Index File (.tbi)",
                            accept = c(".tbi"),
                            multiple = FALSE)
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
                  # Add select/deselect all buttons for entire populations
                  fluidRow(
                    column(6, actionButton("selectAllPopulations", "Select All Populations", width = "100%")),
                    column(6, actionButton("deselectAllPopulations", "Deselect All Populations", width = "100%"))
                  ),
                  # Updated population selection
                  div(class = "population-group", "East Asian"),
                  checkboxGroupInput("populations_east_asian", "East Asian Populations",
                                     choices = c("CHB", "JPT", "CHS", "CDX", "KHV")),

                  div(class = "population-group", "European"),
                  checkboxGroupInput("populations_european", "European Populations",
                                     choices = c("CEU", "TSI", "GBR", "FIN", "IBS")),

                  div(class = "population-group", "African"),
                  checkboxGroupInput("populations_african", "African Populations",
                                     choices = c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB")),

                  div(class = "population-group", "American"),
                  checkboxGroupInput("populations_american", "American Populations",
                                     choices = c("MXL", "PUR", "CLM", "PEL")),

                  div(class = "population-group", "South Asian"),
                  checkboxGroupInput("populations_south_asian", "South Asian Populations",
                                     choices = c("GIH", "PJL", "BEB", "STU", "ITU"))
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
                  textInput("custom_title", "Custom Plot Title (optional)", value = ""),
                  conditionalPanel(
                    condition = "input.analysis_type == 'pca'",
                    checkboxInput("show_ellipses", "Show Confidence Ellipses", value = TRUE),
                    checkboxInput("show_super_pop", "Show Superpopulation", value = TRUE)
                  ),
                  conditionalPanel(
                    condition = "input.analysis_type == 'ancestry_map'",
                    checkboxInput("show_individuals", "Show Individual Populations", value = FALSE)
                  ),
                  actionButton("run_analysis", "Run Analysis", icon = icon("play"))
                ),
                box(
                  title = "Analysis Progress", status = "primary", solidHeader = TRUE,
                  div(id = "analysis_progress",
                      tags$b("Analysis Status:"),
                      verbatimTextOutput("analysis_status")
                  ),
                  progressBar(
                    id = "analysis_progress_bar",
                    value = 0,
                    total = 100,
                    display_pct = TRUE
                  )
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
  # Reactive values to store uploaded files
  rv <- reactiveValues(
    vcf_file = NULL,
    tbi_file = NULL,
    metadata_file = NULL,
    genetic_data = NULL,
    filtered_data = NULL,
    pca_results = NULL
  )

  # Status and progress tracking
  analysis_status <- reactiveVal("Waiting to start analysis")
  analysis_progress <- reactiveVal(0)

  # Output status and progress
  output$analysis_status <- renderText({
    analysis_status()
  })

  # Update progress bar
  observe({
    updateProgressBar(
      session,
      id = "analysis_progress_bar",
      value = analysis_progress()
    )
  })

  # VCF File Upload
  observeEvent(input$vcf_file, {
    req(input$vcf_file)
    # Store the full path, ensuring it's the actual file path
    rv$vcf_file <- input$vcf_file$datapath

    # Extract chromosome from filename
    filename <- input$vcf_file$name
    chr_match <- regexpr("chr[0-9XYM]+", filename, ignore.case = TRUE)
    if (chr_match != -1) {
      chromosome_str <- substr(filename, chr_match, chr_match + attr(chr_match, "match.length") - 1)
      # Remove 'chr' prefix and convert to uppercase
      chromosome_str <- toupper(sub("^chr", "", chromosome_str))

      # Force the chromosome selection to match the extracted chromosome
      updateSelectInput(session, "chromosome", selected = chromosome_str)
    }
  })

  # Tabix Index File Upload
  observeEvent(input$tbi_file, {
    req(input$tbi_file)
    # Store the full path, ensuring it's the actual file path
    rv$tbi_file <- input$tbi_file$datapath
  })


  # Metadata File Upload
  observeEvent(input$metadata_file, {
    req(input$metadata_file)
    rv$metadata_file <- input$metadata_file$datapath
  })

  # Select/Deselect All Populations
  observeEvent(input$selectAllPopulations, {
    updateCheckboxGroupInput(session, "populations_east_asian", selected = c("CHB", "JPT", "CHS", "CDX", "KHV"))
    updateCheckboxGroupInput(session, "populations_european", selected = c("CEU", "TSI", "GBR", "FIN", "IBS"))
    updateCheckboxGroupInput(session, "populations_african", selected = c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"))
    updateCheckboxGroupInput(session, "populations_american", selected = c("MXL", "PUR", "CLM", "PEL"))
    updateCheckboxGroupInput(session, "populations_south_asian", selected = c("GIH", "PJL", "BEB", "STU", "ITU"))
  })

  observeEvent(input$deselectAllPopulations, {
    updateCheckboxGroupInput(session, "populations_east_asian", selected = character(0))
    updateCheckboxGroupInput(session, "populations_european", selected = character(0))
    updateCheckboxGroupInput(session, "populations_african", selected = character(0))
    updateCheckboxGroupInput(session, "populations_american", selected = character(0))
    updateCheckboxGroupInput(session, "populations_south_asian", selected = character(0))
  })

  # Reactive to combine selected populations
  selected_populations <- reactive({
    c(
      input$populations_east_asian,
      input$populations_european,
      input$populations_african,
      input$populations_american,
      input$populations_south_asian
    )
  })

  # Load Genetic Data
  observe({
    req(rv$vcf_file, rv$tbi_file, input$chromosome, input$region_range)

    # Update status
    analysis_status("Preparing to load genetic data...")
    analysis_progress(10)

    # Create GRanges object for region
    region <- GRanges(
      seqnames = as.character(input$chromosome),
      ranges = IRanges(
        start = input$region_range[1],
        end = input$region_range[2]
      )
    )

    # Load Genetic Data using genRelateR's function
    rv$genetic_data <- tryCatch({
      # Verify .tbi file exists in the same directory as VCF
      tbi_expected_path <- paste0(rv$vcf_file, ".tbi")

      # If uploaded .tbi is not in the expected location, copy it
      if (!file.exists(tbi_expected_path)) {
        file.copy(rv$tbi_file, tbi_expected_path, overwrite = TRUE)
      }

      # Update status
      analysis_status("Loading genetic data...")
      analysis_progress(30)

      # Load genetic data
      genRelateR::loadGeneticData(rv$vcf_file, regions = region)
    }, error = function(e) {
      # Update status
      analysis_status(paste("Error loading genetic data:", e$message))
      analysis_progress(0)
      NULL
    })
  })

  # Filter Populations
  observe({
    req(rv$genetic_data, rv$metadata_file)

    # Get selected populations
    populations <- selected_populations()

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

    # Update status and progress
    analysis_status("Starting population structure analysis...")
    analysis_progress(50)

    # Perform Analysis
    rv$pca_results <- tryCatch({
      # Update status
      analysis_status("Analyzing population structure...")
      analysis_progress(75)

      genRelateR::analyzePopulationStructure(
        rv$filtered_data$vcf_data,
        rv$filtered_data$pop_metadata,
        method = "pca"
      )
    }, error = function(e) {
      # Update status
      analysis_status(paste("Error in analysis:", e$message))
      analysis_progress(0)

      showNotification(
        paste("Error in population structure analysis:", e$message),
        type = "error"
      )
      NULL
    })

    # Update final status
    if (!is.null(rv$pca_results)) {
      analysis_status("Analysis complete")
      analysis_progress(100)
    }
  })

  # Plot Results
  output$analysis_plot <- renderPlot({
    req(rv$pca_results, rv$filtered_data)

    # Determine the title
    default_pca_title <- "Population Structure PCA"
    default_ancestry_title <- "Global Population Distribution"

    # Use custom title if provided, otherwise use default
    plot_title <- if(input$custom_title != "") input$custom_title else (
      if(input$analysis_type == "pca") default_pca_title else default_ancestry_title
    )

    if (input$analysis_type == "pca") {
      genRelateR::plotPopulationPca(
        analysis_results = rv$pca_results,
        rv$filtered_data$pop_metadata,
        title = plot_title,
        ellipses = input$show_ellipses,
        super_pop = input$show_super_pop
      )
    } else if (input$analysis_type == "ancestry_map") {
      genRelateR::plotAncestryMap(
        rv$pca_results,
        rv$filtered_data$pop_metadata,
        title = plot_title,
        individual = input$show_individuals
      )
    }
  })
}

# Run the application
shinyApp(ui, server)
