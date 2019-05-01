# Biome-shiny 0.4

# For Permanova
  # 1. Make it so that it takes a distance matrix (output from beta diversity?) instead of abundances
  # Which implies making it so Beta Diversity outputs a distance matrix on user input.

library(shiny)
library(shinydashboard)
library(microbiome)
library(phyloseq)
library(DT)
library(ggplot2)
library(plotly)
library(knitr)
library(dplyr)
library(ggpubr)
library(hrbrthemes)
library(reshape2)
library(DirichletMultinomial)
library(vegan)
#library(limma)


# Load sample datasets #
data("dietswap")
data("atlas1006")
data("peerj32")
peerj32 <- peerj32$phyloseq

# UI
ui <- dashboardPage(
  dashboardHeader(title = "biome-shiny v0.4"),
  dashboardSidebar(
    sidebarMenu(
      menuItem(
        "Introduction",
        tabName = "intro",
        icon = icon("dashboard")
      ),
      br(),
      paste0("Microbiome analysis"),
      menuItem("Core microbiota", tabName = "coremicrobiota"),
      menuItem("Community composition", tabName = "communitycomposition"),
      menuItem("Alpha diversity", tabName = "alphadiversity"),
      menuItem("Beta diversity", tabName = "betadiversity"),
      menuItem("Community landscape", tabName = "landscaping"),
      menuItem("DMM Clustering", tabName = "dirichlet"),
      br(),
      paste0("Statistical analysis"),
      menuItem("PERMANOVA", tabName = "permanova"),
      menuItem("ANOSIM", tabName = "anosim"),
      br(),
      paste0("Outputs and Results"),
      menuItem("Results", tabName = "results")
    ),
    
    br(),
    h5(
      "Made with ",
      img(src = "shiny.png", height = "35"),
      "by ",
      img(src = "biodata.png", height = "25")
    )
  ),
  
  dashboardBody(tabItems(
    #Introduction tab#
    tabItem(
      tabName = "intro",
      h1("Introduction to biome-shiny"),
      br(),
      box(
        width = "2",
        radioButtons(
          "datasetChoice",
          "Dataset",
          c("Upload dataset", "Use sample dataset"),
          selected = "Use sample dataset"
        ),
        conditionalPanel(
          condition = "input.datasetChoice == 'Upload dataset'",
          fileInput(
            "dataset",
            "Please upload a dataset:",
            multiple = TRUE,
            accept = c("text/csv",
                       "text/comma-separated-values,text/plain",
                       ".csv")
          )
        ),
        conditionalPanel(
          condition = "input.datasetChoice == 'Use sample dataset'",
          selectInput(
            "datasetSample",
            "Or choose a sample dataset:",
            choices = c("dietswap", "atlas1006", "peerj32")
          )
        )
      ),
      box(
        paste0(
          "biome-shiny is a microbiome analysis pipeline developed with the Shiny library for R, and based, primarily, on the \"microbiome\" and \"phyloseq\" libraries for analysis.\n\n\n\nBiome-shiny is being developed by BioData.pt for ELIXIR."
        )
      )#, #Test variables please ignore
      # box(title = "Test Box",
      #   paste0("Phyloseq summary"),
      #   verbatimTextOutput("testOutput"),
      #   paste0("RadioButtons tests"),
      #   verbatimTextOutput("testOutput2"),
      #   paste0("Uploaded file stuff"),
      #   verbatimTextOutput("testOutput3")
      # )
    ),
    
    # Core microbiota #
    tabItem(
      tabName = "coremicrobiota",
      box(
        title = "Variables",
        width = "2",
        collapsible = TRUE,
        collapsed = TRUE,
        
        numericInput("detectionPrevalence", "For prevalences: Choose detection value", min = 0.00, max = 100, value = 0.01, step = 0.01),
        numericInput("prevalencePrevalence","Input a prevalence value", min = 0, max = 1, value = 0.5, step = 0.05),
        selectInput("prevalenceSelection", multiple = TRUE, selectize = TRUE, choices = c(seq(0,1,by=0.01)), label = "Choose prevalence values for lineplot and heatmap"),
        textInput("detectionForLineplot", label = "Enter detection values, separated by comma, for the lineplot"),
        textInput("detectionMin", label = "Enter minimum limit for detection"),
        textInput("detectionMax", label = "Enter maximum limit for detection"),
        textInput("maxLength", label = "Enter length (number of values)")
      ),
      
      tabsetPanel(
        tabPanel("Prevalence (relative)",
                 # Output the prevalence in relatives and absolutes (Counts)
                 dataTableOutput("prevalenceRelativeOutput")
        ),
        tabPanel("Prevalence (absolute)",
                 dataTableOutput("prevalenceAbsoluteOutput")
        ),
        tabPanel("Core Taxa Summary",
                 verbatimTextOutput("corePhyloSummary"),
                 verbatimTextOutput("coreTaxa")
        ),
        tabPanel("Core Taxa Visualization",
                 plotlyOutput("coreLineplot"),
                 plotlyOutput("coreHeatmap")
        )
      )
    ),
    
    
    # Community composition #
    tabItem(
      tabName = "communitycomposition",
      box(
        title = "Variables",
        width = "2",
        collapsible = TRUE,
        collapsed = TRUE,
        selectInput(
          "z1",
          "Choose a metadata column:",
          # For subsetting data, #1
          choices = colnames("datasetMetadata"),
          selected = "bmi_group"
        ),
        selectInput(
          "v1",
          "Choose a metadata value:",
          # For subsetting data, metadata value
          choices = sapply("datasetMetadata", levels),
          selected = "lean"
        ),
        selectInput(
          "z2",
          "Choose a seocond metadata column:",
          # For subsetting data, #2
          choices = colnames("datasetMetadata"),
          selected = "nationality"
        ),
        selectInput(
          "v2",
          "Choose a second metadata value:",
          # For subsetting data, metadata #2 value
          choices = sapply("datasetMetadata", levels),
          selected = "AAM"
        ),
        selectInput(
          "z3",
          "Choose the intended timepoint:",
          # For subsetting data, timepoint data
          choices = colnames("datasetMetadata"),
          selected = "timepoint.within.group"
        ),
        selectInput(
          "v3",
          "Choose a timepoint value:",
          # For subsetting data, timepoint data value
          choices = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
          selected = "1"
        ),
        selectInput(
          "v4",
          "Choose a taxonomy rank:",
          # Tax rank to analyze
          choices = c("Phylum", "Class", "Order", "Family", "Genus"),
          selected = "Phylum"
        )
      ),
      
      tabsetPanel(
        tabPanel(title = "Abundance in samples by taxa",
                 plotOutput("communityPlot")),
        tabPanel(title = "Relative abundance",
                 plotOutput("communityPlotGenus")),
        tabPanel(title = "Community Plot #3",
                 plotOutput("communityBarplot")),
        tabPanel(title = "Taxonomy Prevalence Plot",
                 plotlyOutput("communityPrevalence"))
      )
    ),
    
    # Alpha Diversity #
    tabItem(
      tabName = "alphadiversity",
      box(
        width = "2",
        collapsible = TRUE,
        title = "Variables",
        collapsed = TRUE,
        # #X (the metadata) and Y (the diversity measure)
        selectInput(
          "x",
          "Choose a metadata column test:",
          choices = colnames("datasetMetadata")
        ),
        selectInput("y", "Choose a diversity measure:",
                    choices = colnames("datasetMetadata"))
      ),
      br(),
      tabsetPanel(
        type = "tabs",
        id = "tabsetpanel",
        
        # Tab 1: Phyloseq summary
        tabPanel(title = "Phyloseq Summary",
                 verbatimTextOutput("summary")),
        
        
        # Tab 2: Alpha diversity measures with metadata
        tabPanel(title = "Table",
                 DT::dataTableOutput("view")),
        
        # Tab 3: The exceptionally lewd Violin Plots
        tabPanel(title = "Violin Plot",
                 plotlyOutput("violinPlot")),
        
        # Tab 4: A phyloseq richness plot
        tabPanel(title = "Richness Plot",
                 plotlyOutput("richnessPlot"))
      )
    ),
    
    #Beta diversity#
    tabItem(
      tabName = "betadiversity",
      box(
        title = "Variables",
        width = "2",
        collapsible = TRUE,
        collapsed = TRUE,
        selectInput(
          "xb", "Choose a metadata column:", choices = colnames("datasetMetadata"), selected = "bmi_group"
        ),
        
        selectInput(
          "yb",
          "For metadata/metadata split plots, choose a metadata column:",
          choices = colnames("datasetMetadata"),
          selected = "bmi_group"
        ),
        
        selectInput(
          "zb",
          "For metadata/tax rank split plot, choose a taxonomy rank:",
          choices = c("Phylum", "Class", "Order", "Family", "Genus")
        ),
        
        selectInput(
          "ordinate.method",
          "Choose an ordination method:",
          choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
          selected = "CCA"
        ),
        
        selectInput(
          "ordinate.distance",
          "Choose a distance method:",
          choices = c("bray", "jaccard", "unifrac"),
          selected = "unifrac"
        ),
        #There's 44 of these in total
        
        sliderInput(
          "geom.size",
          "Plot geometry point size:",
          min = 1,
          max = 10,
          step = 0.5,
          value = "3"
        )
      ),
      
      tabsetPanel(
        tabPanel(title = "Ordination Plot",
                 plotlyOutput("ordinatePlot"),
                 textOutput("ordinatePrint")
                ),
        
        tabPanel(title = "Split Ordination Plot (Metadata/Metadata)",
                 plotlyOutput("splitOrd")),
        
        tabPanel(title = "Taxa Plot",
                 plotlyOutput("taxaOrd"))
      )
      
    ),
    
    tabItem(
       tabName = "landscaping",
       box(title = "Variables", width = "2", collapsible = TRUE, collapsed = TRUE,
           selectInput("ordinateMethodLandscape", "Choose an ordination method:", 
                       choices=c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"), selected = "CCA"),
           selectInput("ordinateDistanceLandscape", "Choose a distance method:", 
                       choices=c("bray","jaccard","unifrac"), selected = "unifrac"), #There's 44 of these in total
           numericInput("detectionLandscape", label="Input detection threshold for landscape analysis", min = 0, max = 100, step = 0.1, value = "0.1"),
           numericInput("prevalenceLandscape", label="Input prevalence percentage of samples for landscape analysis", min = 0, max = 100, step = 0.1, value = "50"),
           selectInput("metadataLandscape1", "Choose first metadata column to subset data:", choices=colnames("datasetMetadata")),
           selectInput("metadataLandscapeValue1", "Choose first metadata value:", choices=sapply("datasetMetadata", levels)),
           selectInput("metadataLandscape2", "Choose second metadata column to subset data:", choices=colnames("datasetMetadata")),
           selectInput("metadataLandscapeValue2", "Choose second metadata value:", choices=sapply("datasetMetadata", levels))
       ),
       tabsetPanel(
         tabPanel("PCA",
                  plotlyOutput("landscapePCA")
         ),
         tabPanel("PCoA/MDS/NMDS",
                  plotlyOutput("landscapeOrdination"),
                  plotlyOutput("landscapeOrdinationSamplenames")
         ),
         tabPanel("t-SNE",
                  plotlyOutput("landscapeTSne")
         ),
         tabPanel("Abundance histogram (one-dimensional landscape)",
                  plotlyOutput("landscapeAbundanceHistAbs"),
                  plotlyOutput("landscapeAbundanceHistRel")
         )
       )
       
    ),
    tabItem(
      tabName = "dirichlet",
      box( title = "Variables", width = "2", collapsible = TRUE, collapsed = TRUE,
           selectInput("datasetDMM","Choose the dataset to analyze.",
                       choices = c("dietswap","atlas1006","peerj32")),
           numericInput("detectionDMM", label="Input detection threshold for DMM", min = 0, max = 100, step = 0.1, value = "0.1"),
           numericInput("prevalenceDMM", label="Input prevalence percentage of samples for DMM0", min = 0, max = 100, step = 0.1, value = "50"),
           numericInput("maxModelsDMM", label="Input maximum number of community types for the DMM model, or leave empty for infinite types at the cost of performance.", value = "3")
      ),
      tabsetPanel(
        tabPanel("Model Verification Lineplot",
                 plotOutput("dmmModelCheck")
        ),
        
        tabPanel("Alpha and Theta Parameters",
                 dataTableOutput("dmmParameters"),
                 dataTableOutput("sampleAssignments")
        ),
        
        tabPanel("Taxa Contribution Per Community Type",
                 plotlyOutput("taxaContributionPerComponent")
        )
      )
    ),
    tabItem(
      tabName = "permanova",
      box( title = "Variables", width = "2", collapsible = TRUE, collapsed = TRUE,
           selectInput("permanovaDistanceMethod","Select distance method", choices = c("bray","jacard","unifrac"), selected = "unifrac"),
           selectInput("permanovaMethod","Select ordination method",
                       choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                       selected = "CCA"),
           selectInput("permanovaColumn","Select metadata for density plot", choices = colnames("datasetMetadata")),
           sliderInput("permanovaPlotSize", "Plot point size", min = 0.5, max = 10, step = 0.5, value = "3"),
           numericInput("permanovaPermutations", "Number of permutations", min = 1, step = 1, value = 99)
      ),
      tabsetPanel(
        tabPanel( title = "Population density plot",
           plotlyOutput("densityPlot")
        ),
        tabPanel( title = "P-Value",
                  dataTableOutput("pValue"),
                  dataTableOutput("homogeniety")
        ),
        tabPanel ( title = "Top Factors",
                   plotOutput("topFactorPlot") 
        ),
        tabPanel ( title = "Network Plot",
                   plotlyOutput("netPlot")
        ),
        tabPanel ( title = "Heatmap",
                   plotlyOutput("permaHeatmap")
        )
      )
      ),
      
      tabItem(tabName = "anosim",
              box( title = "Variables", width = "2", collapsible = TRUE, collapsed = TRUE,
                   selectInput("anosimDistanceMethod","Select distance method", choices = c("bray","jacard","unifrac"), selected = "unifrac"),
                   selectInput("anosimMethod","Select ordination method",
                               choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                               selected = "CCA"),
                   selectInput("anosimColumn","Select metadata for density plot", choices = colnames("datasetMetadata")),
                   sliderInput("anosimPlotSize", "Plot point size", min = 0.5, max = 10, step = 0.5, value = "3"),
                   numericInput("anosimPermutations", "Number of permutations", min = 1, step = 1, value = 99)
              ),
              textOutput("pValueAnosim"),
              textOutput("ordinateDataPrint")
      ),
    # tabItem(tabName = "limma",
    #     box( title = "Variables", width = "2", collapsible = TRUE, collapsed = TRUE,
    #          selectInput("limmaColumn","Select metadata for density plot", choices = colnames("datasetMetadata"))
    #     ),
    #     dataTableOutput("limmaTable")
    #   ),
    
    tabItem(tabName = "results"
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  datasetInput <- reactive({
    if (input$datasetChoice == "Use sample dataset") {
      switch(
        input$datasetSample,
        "dietswap" = dietswap,
        "atlas1006" = atlas1006,
        "peerj32" = peerj32
      )
    } else {
      req(input$dataset)
      datapath <- input$dataset$datapath
      biomfile <- import_biom(datapath)
      return(biomfile)
    }
  })
  
  ## Core Microbiota ##
  
  
  prevalenceAbsolute <- reactive({
    as.data.frame(prevalence(compositionalInput(), detection = input$detectionPrevalence/100, sort = TRUE, count = TRUE))
  })
  prevalenceRelative <- reactive({
    as.data.frame(prevalence(compositionalInput(), detection = input$detectionPrevalence/100, sort = TRUE))
  })
  
  # Prevalence output - DT is absolutely fantastic
  output$prevalenceAbsoluteOutput <- renderDT({
    datatable(prevalenceAbsolute())
  })
  output$prevalenceRelativeOutput <- renderDT({
    datatable(prevalenceRelative())
  })
  
  # Core analysis
  corePhylo <- reactive({ # Generates a phyloseq file from input
    core(compositionalInput(), detection = input$detectionPrevalence, prevalence = input$prevalencePrevalence )
  })
  output$corePhyloSummary <- renderPrint({ # Summary of corePhylo file
    summarize_phyloseq(corePhylo())  
  })
  output$coreTaxa <- renderPrint({ # Reports the taxa in corePhylo
    taxa(corePhylo())
  })
  
  # Visualization (lineplots and heatmaps)
  output$coreLineplot <- renderPlotly({
    prevalences <- as.numeric(input$prevalenceSelection)
    detections <- as.numeric(unlist(strsplit(input$detectionForLineplot, split = ",")))/100
    print(detections)
    lineplot <- plot_core(compositionalInput(), prevalences = prevalences, detections = detections, plot.type = "lineplot") + xlab("Relative Abundance (%)")
    plotly_build(lineplot)
  })
  
  output$coreHeatmap <- renderPlotly({
    # Core with compositionals:
    prevalences <- as.numeric(input$prevalenceSelection)
    detections <- 10^seq(log10(as.numeric(input$detectionMin)), log10(as.numeric(input$detectionMax)), length = as.numeric(input$maxLength))
    gray <- gray(seq(0,1,length=5))
    coreplot <- plot_core(compositionalInput(), plot.type = "heatmap", colours = gray, prevalences = prevalences, detections = detections, min.prevalence = .5) + xlab("Detection Threshold (Relative Abundance (%))")
    plotly_build(coreplot)
  })
  
  
  
  ## Community Composition ##
  
  # Updating SelectInputs when database changes #
  observe({
    updateSelectInput(session, "z1",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "z2",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "z3",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "v4",
                      choices = colnames(tax_table(datasetInput())))
  })
  
  #Update metadata value selectInputs for CC Analysis
  observeEvent(input$z1, {
    updateSelectInput(session, "v1",
                      choices = (sample_data(datasetInput())[[input$z1]]))
  })
  observeEvent(input$z2, {
    updateSelectInput(session, "v2",
                      choices = (sample_data(datasetInput())[[input$z2]]))
  })
  observeEvent(input$z3, {
    updateSelectInput(session, "v3",
                      choices = (sample_data(datasetInput())[[input$z3]]))
  })
  
  # Subset the database #
  datasetSubsetInput <- reactive({
    subset1 <-
      prune_samples(sample_data(datasetInput())[[input$z1]] == input$v1, datasetInput())
    subset1 <-
      prune_samples(sample_data(subset1)[[input$z2]] == input$v2, subset1)
    subset1 <-
      prune_samples(sample_data(subset1)[[input$z3]] == input$v3, subset1)
    subset2 <-
      prune_samples(sample_data(subset1)[[input$z1]] == input$v1, subset1) %>% aggregate_taxa(level = input$v4)
    microbiome::transform(subset2, "compositional")
  })
  
  
  # Make plots
  output$communityPlot <- renderPlot ({
    theme_set(theme_bw(21))
    
    communityplot <-
      datasetSubsetInput() %>% plot_composition(sample.sort = "Bacteroidetes", otu.sort = "abundance") +
      scale_fill_manual(values = default_colors("Phylum")[taxa(datasetSubsetInput())])
    print(communityplot)
    #plotly_build(communityplot)
  }) #otu.sort and and sample.sort need to be selectable
  
  # Make plot
  output$communityPlotGenus <- renderPlot({
    compositionplot <- plot_composition(
      datasetSubsetInput(),
      taxonomic.level = "Genus",
      sample.sort = "nationality",
      x.label = "nationality"
    ) +
      guides(fill = guide_legend(ncol = 1)) +
      scale_y_percent() +
      labs(
        x = "Samples",
        y = "Relative abundance (%)",
        title = "Relative abundance data",
        subtitle = "Subtitle",
        caption = "Caption text."
      ) +
      theme_ipsum(grid = "Y")
    print(compositionplot)
    #plotly_build(compositionplot)
  })
  
  # Barplot, averaged by z1 group - replace with user pick later
  output$communityBarplot <- renderPlot({
    compplot <- plot_composition(datasetSubsetInput(), average_by = input$z1, plot.type = "barplot")
    print(compplot)
    #plotly_build(compplot)
  })
  
  # And top it off with a taxa prevalence plot
  output$communityPrevalence <- renderPlotly({
    prevplot <- plot_taxa_prevalence(datasetSubsetInput(), input$v4) #Can be changed to whatever taxonomic rank the input files have
    plotly_build(prevplot)
  })
  
  
  
  ## Alpha Diversity ##
  
  # Updating SelectInputs when database changes #
  observe({
    updateSelectInput(session, "x",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "y",
                      choices = colnames(alpha(datasetInput())))
  })
  
  # Phyloseq Summary #
  output$summary <- renderPrint({
    req(datasetInput())
    summarize_phyloseq(datasetInput())
  })
  
  # Merged table - generate and output #
  mergedTable <- reactive({
    merge(meta(datasetInput()), alpha(datasetInput()), all.y = TRUE)
  })
  output$view <-
    DT::renderDataTable({
      #Issue with search -> search by specific term (like "male") doesn't work -> causes problem when searching between female (which gives all female samples) and male (which give both male and female samples). Still, a slightly dysfunctional search function is better than none at all
      datatable(mergedTable())
    })
  
  # Violin plot #
  output$violinPlot <- renderPlotly({
    violin <- ggviolin(
      mergedTable(),
      x = input$x,
      y = input$y,
      add = "boxplot",
      fill = input$x,
      palette = c("#a6cee3", "#b2df8a", "#fdbf6f")
    )
    plotly_build(violin)
  })
  
  # Richness Plot #
  output$richnessPlot <- renderPlotly({
    richnessplot <- plot_richness(
      datasetInput(),
      x = input$x,
      measures = c("Shannon", "Simpson"),
      color = input$x
    )
    plotly_build(richnessplot)
  })
  
  ## Report - Alpha Diversity ##
  
  # For PDF output, change this to "report.pdf"
  output$report <- downloadHandler(
    filename = "report.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(n = input$slider)
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  ## Beta Diversity ##
  
  # Updating SelectInputs #
  observe({
    updateSelectInput(session, "xb",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "yb",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "zb",
                      choices = colnames(tax_table(datasetInput())))
  })
  
  compositionalInput <- reactive({
    microbiome::transform(datasetInput(), "compositional")
  })
  
  ordinateData <- reactive({
    ordinate(
      compositionalInput(),
      method = input$ordinate.method,
      distance = input$ordinate.distance
    ) #note to self, make method choosable from dropdown
  })
  
  output$ordinatePrint <- renderPrint({
    ordinateData()
  })
  
  output$ordinatePlot <- renderPlotly({
    p <- phyloseq::plot_ordination(datasetInput(), ordinateData(), color = input$xb ) + geom_point(size = input$geom.size)
    plotly_build(p)
  })
  
  output$splitOrd <- renderPlotly({
    splitOrdplot <-
      plot_ordination(
        datasetInput(),
        ordinateData(),
        type = "split",
        shape = input$xb,
        color = input$yb,
        label = input$zb
      )
    plotly_build(splitOrdplot)
  })
  
  output$taxaOrd <- renderPlotly({
    taxaOrdplot <-
      plot_ordination(
        datasetInput(),
        ordinateData(),
        type = "taxa",
        color = input$zb,
        label = input$xb
      )
    plotly_build(taxaOrdplot)
  })
  
  
  # LANDSCAPE ANALYSIS #
  #Update metadata selectInputs when swapping datasets - Landscape analysis
  observe({
    updateSelectInput(session, "metadataLandscape1",
                      choices=colnames(meta(datasetInput())))
  })
  observe({
    updateSelectInput(session, "metadataLandscape2",
                      choices=colnames(meta(datasetInput())))
  })
  
  #Update metadata values for landscape analysis 
  observeEvent(input$metadataLandscape1,{
    updateSelectInput(session, "metadataLandscapeValue1",
                      choices=(sample_data(datasetInput())[[input$metadataLandscape1]]))
  })
  observeEvent(input$metadataLandscape2,{
    updateSelectInput(session, "metadataLandscapeValue2",
                      choices=(sample_data(datasetInput())[[input$metadataLandscape2]]))
  })
  
  # 1 - Select the column 
  datasetSubsetLandscape <- reactive({
    subset <- prune_samples(sample_data(datasetInput())[[input$metadataLandscape1]] == input$metadataLandscapeValue1, datasetInput())
    subset <- prune_samples(sample_data(subset)[[input$metadataLandscape2]] == input$metadataLandscapeValue2, subset)
    microbiome::transform(subset, "compositional")
  })
  
  #PCA plot
  
  output$landscapePCA <- renderPlotly({ #Change transformation and method type to user input
    p <- plot_landscape(datasetInput(), method = "PCA", transformation = "clr") +
      labs(title = paste("PCA / CLR"))
    plotly_build(p)
  })
  
  #NMDS/PCoa/MDS
  output$landscapeOrdination <- renderPlotly({
    x <- datasetInput()
    quiet(x.ord <- ordinate(x, input$ordinateMethodLandscape, input$ordinateDistanceLandscape))
    proj <- phyloseq::plot_ordination(x, x.ord, justDF=TRUE)
    names(proj)[1:2] <- paste("Comp", 1:2, sep=".")
    p <- plot_landscape(proj[, 1:2], col = proj$nationality, legend = T)
    plotly_build(p)
  })
  output$landscapeOrdinationSamplenames <- renderPlotly({
    x <- datasetInput()
    quiet(x.ord <- ordinate(x, input$ordinateMethodLandscape, input$ordinateDistanceLandscape))
    proj <- phyloseq::plot_ordination(x, x.ord, justDF=TRUE)
    names(proj)[1:2] <- paste("Comp", 1:2, sep=".")
    ax1 <- names(proj)[[1]]
    ax2 <- names(proj)[[2]]
    p <- ggplot(aes_string(x = ax1, y = ax2, label = "sample"), data = proj) +
      geom_text(size = 2)
    plotly_build(p)
  })
  
  #
  output$landscapeTSne <- renderPlotly({
    p <- plot_landscape(datasetInput(), "t-SNE",
                   distance = "euclidean", transformation = "hellinger") +
      labs(title = paste("t-SNE / Hellinger / Euclidean"))
    plotly_build(p)
  })
  
  output$landscapeAbundanceHistAbs <- renderPlotly({
    # Visualize population densities for specific taxa
    p <- plot_density(datasetInput(), "Dialister") + ggtitle("Absolute abundance")
    plotly_build(p)
  })
  
  output$landscapeAbundanceHistRel <- renderPlotly({
    # Visualize population densities for specific taxa
    x <- microbiome::transform(datasetInput(), "compositional")
    tax <- "Dialister"
    p <- plot_density(x, tax, log10 = TRUE) +
      ggtitle("Relative abundance") +
      xlab("Relative abundance (%)")
    plotly_build(p)
  })
  
  
  
  #DMM COMMUNITY TYPING#
  
  # Pick OTU count matrix, convert into samplex x taxa format and fit the DMM model
  dmmModelFit <- reactive({
    taxa <- core_members(compositionalInput(), detection = input$detectionDMM/100, prevalence = input$prevalenceDMM/100) #need to make this related to core comp, along with the coremembers in beta diversity
    pseq <- prune_taxa(taxa, datasetInput())
    dat <- abundances(pseq)
    count <- as.matrix(t(dat))
    mclapply(1:input$maxModelsDMM, dmn, count = count, verbose=TRUE)
  })
  
  output$dmmModelCheck <- renderPlot({
    lplc <- sapply(dmmModelFit(), laplace) # AIC / BIC / Laplace
    aic  <- sapply(dmmModelFit(), AIC) # AIC / BIC / Laplace
    bic  <- sapply(dmmModelFit(), BIC) # AIC / BIC / Laplace
    plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
    lines(aic, type="b", lty = 2)
    lines(bic, type="b", lty = 3)
  })
  
  dmmModelBestFit <- reactive({
    lplc <- sapply(dmmModelFit(), laplace) # AIC / BIC / Laplace
    aic  <- sapply(dmmModelFit(), AIC) # AIC / BIC / Laplace
    bic  <- sapply(dmmModelFit(), BIC) # AIC / BIC / Laplace
    dmmModelFit()[[which.min(lplc)]]
  })
  
  output$dmmParameters <- renderDT({
    mixturewt(dmmModelBestFit())
  })
  output$sampleAssignments <- renderDT({
    apply(mixture(dmmModelBestFit()), 1, which.max)
    
  })
  
  output$taxaContributionPerComponent <- renderPlotly({
    for (k in seq(ncol(fitted(dmmModelBestFit())))) {
      d <- melt(fitted(dmmModelBestFit()))
      colnames(d) <- c("OTU", "cluster", "value")
      d <- subset(d, cluster == k) %>%
        # Arrange OTUs by assignment strength
        arrange(value) %>%
        mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
        # Only show the most important drivers
        filter(abs(value) > quantile(abs(value), 0.8))     
      
      p <- ggplot(d, aes(x = OTU, y = value)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title = paste("Top drivers: community type", k))
    }
    ggplotly(p)
  })
  
  
  ###########################################
  # output$testOutput2 <- renderPrint({
  #     print(typeof(input$datasetChoice))
  #     print(input$datasetChoice)
  # })
  # output$testOutput2 <- renderPrint({
  #   print(typeof(input$datasetChoice))
  #   print(input$datasetChoice)
  # })
  # output$testOutput3 <- renderPrint({
  #   print(typeof(input$dataset))
  #   print(input$dataset)
  #   print(datasetInput())
  # })
  ###########################################
  
  ###########################
  ## Statistical analysis ###
  ###########################
  
  ## PERMANOVA ##
  
  #Update metadata column#
  observe({
    updateSelectInput(session, "permanovaColumn",
                      choices = colnames(meta(datasetInput())))
  })
  
  permanova <- reactive({
    otu <- abundances(compositionalInput())
    meta <- meta(compositionalInput())
    permnumber <- input$permanovaPermutations
    metadata <- input$permanovaColumn
    adonis(t(otu) ~ meta[[metadata]],
                        data = meta, permutations = permnumber, method = "bray", parallel = getOption("mc.cores")
    )
  })
  
  output$permaprint <- renderText({
      print(abundances(compositionalInput()))
      print(meta(datasetInput()))
      print(typeof(input$permanovaColumn))
      print(typeof(meta(datasetInput())$sex))
      print(input$permanovaPermutations)
  })
  
  output$densityPlot <- renderPlotly({
    p <- plot_landscape(compositionalInput(), method = input$permanovaMethod, distance = input$permanovaDistanceMethod, col = input$permanovaColumn, size = input$permanovaPlotSize)
    plotly_build(p)
  })
  
  # output$pValue <- renderPrint({
  #   as.data.frame(permanova()$aov.tab)[paste0('"',input$permanovaColumn,'"'), "Pr(>F)"]
  # })
  
  output$pValue <- renderDataTable({
    as.data.frame(permanova()$aov.tab)
  })
    
  output$homogeniety <- renderDataTable({
    otu <- abundances(compositionalInput())
    meta <- meta(compositionalInput())
    dist <- vegdist(t(otu))
    metadata <- input$permanovaColumn
    anova(betadisper(dist, meta[[metadata]]))
  })
  
  output$topFactorPlot <- renderPlot({
    otu <- abundances(compositionalInput())
    meta <- meta(compositionalInput())
    permnumber <- input$permanovaPermutations
    metadata <- input$permanovaColumn
    column <- meta[[metadata]]
    permanova <- adonis(t(otu) ~ column,
           data = meta, permutations = permnumber, method = "bray"
    )
    coef <- coefficients(permanova)["column1",]
    top.coef <- coef[rev(order(abs(coef)))[1:20]] #top 20 coefficients
    par(mar = c(3, 14, 2, 1))
    p <- barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
    print(p)
  })
  
  output$netPlot <- renderPlotly({
    n <- make_network(compositionalInput(), type = "otu", distance = "bray")
    p <- plot_network(n)
    plotly_build(p)  
  })
  
  output$permaHeatmap <- renderPlotly({
    p <- plot_heatmap(compositionalInput(), distance = ordinate(compositionalInput(), distance = input$ordinate.distance), method = input$ordinate.method)
    plotly_build(p, color = meta[[metadata]])
  })
  
  # ANOSIM #
  #Update metadata column#
  observe({
    updateSelectInput(session, "anosimColumn",
                      choices = colnames(meta(datasetInput())))
  })
  
  anosim <- reactive({
    otu <- abundances(compositionalInput())
    meta <- meta(compositionalInput())
    permnumber <- input$anosimPermutations
    metadata <- input$anosimColumn
  })
  
  output$pValueAnosim <- renderPrint({
    otu <- abundances(compositionalInput())
    meta <- meta(compositionalInput())
    metadata <- input$anosimColumn
    vegan::anosim(otu, meta[[metadata]], permutations = 99, distance = "bray", parallel = getOption("mc.cores"))
  })
    
  
  # LIMMA #
  
  # observe({
  #   updateSelectInput(session, "limmaColunm",
  #                     choices = colnames(meta(datasetInput())))
  # })
  
  # fitModel <- reactive({
  #   pseq <- datasetInput() # Rename the example data
  #   otu <- abundances(microbiome::transform(pseq, "log10"))
  #   meta <- meta(pseq)
  #   design <- cbind(intercept = 1, Grp2vs1 = meta[[input$limmaColumn]])
  #   rownames(design) <- rownames(meta)
  #   design <- design[colnames(otu), ]
  #   coef.index <- 2
  #   fit <- lmFit(otu, design)
  #   fit <- eBayes(fit)
  #   pvalues.limma = fit$p.value[, 2]
  #   efs.limma <-  fit$coefficients[, "Grp2vs1"]
  #   kable(topTable(fit, coef = coef.index, p.value=0.1), digits = 2)
  # })
  # 
  # output$limmaTable <- renderDT({
  #   fitModel()
  # })
  
}

# Run the application
shinyApp(ui = ui, server = server, options = list(port = 5000))
