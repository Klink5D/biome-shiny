# biome-shiny: A Shiny R app for metagenomics analysis, built around the "microbiome" package
# current version: v0.0.0.3

#Considerations: modal dialog to fullscreen plots on click/download plots
#                Separate plots into different tabs, along with the sidebar choices you need to generate said plots
#                File input that takes in a BIOM/Mothur/QIIME file and loads it as the active dataset
#                Log that dumps the code running in the background to a text file
#                Also a report function that dumps all plots, summaries and so on onto an HTML page or PDF/word document file.


# Load the required libraries
library("shiny")
library("shinythemes")
library("microbiome")
library("knitr")
library("dplyr")
library("ggpubr")
library("ggplot2")
library("hrbrthemes")
library("reshape2")
library("DT")
library("DirichletMultinomial")

#load the example datasets
data("dietswap")
data("atlas1006")
data("peerj32")
peerj32 <- peerj32$phyloseq

#UI#
ui <- navbarPage("biome-shiny v0.0.0.4", fluid = TRUE,
                 
                 theme = shinytheme("united"),
                 
                 tabPanel("Introduction",
                          titlePanel("Welcome to biome-shiny", windowTitle = "biome-shiny v0.0.2"),
                          
                          sidebarPanel(
                            fileInput("datase", "Please upload a dataset (BIOM, QIIME or MOTHUR file).",
                                      multiple = FALSE,
                                      accept = c("text/csv",
                                                 "text/comma-separated-values,text/plain",
                                                 ".csv")),
                            
                            
                            h5("Made with ",
                               img(src = "shiny.png", height = "50"), "by ", img(src = "biodata.png", height = "30"))
                          ),
                          
                          mainPanel(
                            verbatimTextOutput("introText")
                          )
                   ),
                 
                 tabPanel( "Alpha Diversity",
                           titlePanel("Alpha Diversity", windowTitle = "biome-shiny v0.0.0.4"), 

                           #The sidebar
                           sidebarPanel(
                             selectInput("dataset","Choose the dataset to analyze.",
                                         choices = c("dietswap","atlas1006","peerj32")),
                            # #X (the metadata) and Y (the diversity measure)
                             selectInput("x", "Choose a metadata column test:", 
                                         choices=colnames("testies")),
                             selectInput("y", "Choose a diversity measure:", 
                                         choices=colnames("testies")),
                             #Linebreaks
                             br(),
                             
                             h5("Made with ",
                                img(src = "shiny.png", height = "50"), "by ", img(src = "biodata.png", height = "30"))
                           ),
                           
                           mainPanel(
                             tabsetPanel( type = "tabs",
                                          id = "tabsetpanel",
                                          
                                          # Tab 1: Phyloseq summary
                                          tabPanel( title = "Phyloseq Summary",
                                                    verbatimTextOutput("summary")
                                          ),
                                          
                                          # Tab 2: Alpha diversity measures with metadata
                                          tabPanel( title = "Table",
                                                    DT::dataTableOutput("view")
                                          ),
                                          
                                          # Tab 3: The exceptionally lewd Violin Plots
                                          tabPanel( title = "Violin Plot",
                                                    plotOutput("violinPlot")
                                          ),
                                          
                                          # Tab 4: A phyloseq richness plot
                                          tabPanel( title = "Richness Plot",
                                                    plotOutput("richnessPlot")
                                          )
                             )
                           )
                 ),
                 tabPanel("Beta Diversity",
                          fluidPage(
                            titlePanel("Beta Diversity", windowTitle = "biome-shiny v0.0.0.4"),
                            
                            #The sidebar
                            sidebarPanel(
                              selectInput("datasetBeta","Choose the dataset to analyze.", #Bug where selections shared with Alpha Div. rely on Alpha's setting to function
                                          choices = c("dietswap","atlas1006","peerj32")),
                              
                              numericInput("beta.sneed","Set seed for outcome replication:",
                                           value = 1, min = 0),
                              
                              selectInput("xb", "Choose a metadata column:", 
                                          choices=colnames("testies"), selected = "bmi_group"),
                              
                              selectInput("xc", "For metadata/metadata split plots, choose a metadata column:", 
                                          choices=colnames("testies"), selected = "bmi_group"),
                              
                              selectInput("xd", "For metadata/tax rank split plot, choose a taxonomy rank:", 
                                          choices=c("Phylum", "Class", "Order", "Family","Genus")),
                              
                              selectInput("ordinate.method", "Choose an ordination method:", 
                                          choices=c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"), selected = "CCA"),
                              
                              selectInput("ordinate.distance", "Choose a distance method:", 
                                          choices=c("bray","jaccard","unifrac"), selected = "unifrac"), #There's 44 of these in total
                              
                              sliderInput("geom.size", "Plot geometry point size:", 
                                          min = 1, max = 10, step = 0.5, value = "3"),
                              
                              # Linebreaks
                              br(),
                              
                              h5("Made with ",
                                 img(src = "shiny.png", height = "50"), "by ", img(src = "biodata.png", height = "30"))
                            ),
                            
                            mainPanel(
                              
                              textOutput("betaseed"),
                              
                              tabsetPanel(
                              tabPanel( title = "Ordination Plot", 
                                        plotOutput("ordinatePlot", hover = "plot_hover")),
                              
                              tabPanel( title = "Split Ordination Plot (Metadata/Metadata)", 
                                        plotOutput("metaMetaSplitOrd", hover = "plot_hover_split")),
                              
                              tabPanel( title = "Taxa Plot",
                                        plotOutput("metaTaxSplitOrd", hover = "plot_hover_split_tax"))
                            )
                          )
                        )
                 ),
                 tabPanel("Community Composition",
                          
                          fluidPage(
                            titlePanel("Community Composition", windowTitle = "biome-shiny v0.0.0.4"),
                            
                            #The sidebar
                            sidebarPanel(
                              selectInput("datasetComp","Choose the dataset to analyze.",
                                          choices = c("dietswap","atlas1006","peerj32")),
                              selectInput("z1", "Choose a metadata column:", # For subsetting data, #1
                                          choices=colnames("testies"), selected = "bmi_group"),
                              selectInput("v1", "Choose a metadata value:", # For subsetting data, metadata value
                                          choices=sapply("testies", levels), selected = "lean"),
                              selectInput("z2", "Choose a seocond metadata column:", # For subsetting data, #2
                                          choices=colnames("testies"), selected = "nationality"),
                              selectInput("v2", "Choose a second metadata value:", # For subsetting data, metadata #2 value
                                          choices=sapply("testies", levels), selected = "AAM"),
                              selectInput("z3", "Choose the intended timepoint:", # For subsetting data, timepoint data
                                          choices=colnames("testies"), selected = "timepoint.within.group"),
                              selectInput("v3", "Choose a timepoint value:", # For subsetting data, timepoint data value
                                          choices=c("0","1","2","3","4","5","6","7","8","9","10"), selected = "1"),
                              selectInput("v4", "Choose a taxonomy rank:", # Tax rank to analyze
                                          choices=c("Phylum", "Class", "Order", "Family","Genus"), selected = "Phylum"),
                              
                              
                              br(),
                              
                              h5("Made with ",
                                 img(src = "shiny.png", height = "50"), "by ", img(src = "biodata.png", height = "30"))
                            ),
                            
                            mainPanel(
                            tabsetPanel(
                              tabPanel( title = "Abundance in samples by taxa",
                              plotOutput("communityPlot")
                              ),
                              tabPanel(title = "Relative abundance",
                              plotOutput("communityPlotGenus") 
                              ),
                              tabPanel(title = "Community Plot #3",
                              plotOutput("CommunityHeatmap2") 
                              ),
                              tabPanel(title = "Taxonomy Prevalence Plot",
                              plotOutput("communityPrevalence")) 
                            )
                            )
                          )
                 ),
            tabPanel("Core microbiota analysis", #Something I want to to is use this to subset the data
                fluidPage(
                    titlePanel("Core microbiota analysis", windowTitle = "biome-shiny v0.0.0.4"),
                    sidebarPanel(
                      
                      selectInput("datasetCore","Choose the dataset to analyze.",
                                  choices = c("dietswap","atlas1006","peerj32")),
                      numericInput("detectionPrevalence", "For prevalences: Choose detection value", min = 0.00, max = 100, value = 0.01, step = 0.01),
                      numericInput("prevalencePrevalence","Input a prevalence value", min = 0, max = 1, value = 0.5, step = 0.05),
                      selectInput("prevalenceSelection", multiple = TRUE, selectize = TRUE, choices = c(seq(0,1,by=0.01)), label = "Choose prevalence values for lineplot and heatmap"),
                      textInput("detectionForLineplot", label = "Enter detection values, separated by comma, for the lineplot"),
                      textInput("detectionMin", label = "Enter minimum limit for detection"),
                      textInput("detectionMax", label = "Enter maximum limit for detection"),
                      textInput("maxLength", label = "Enter length (number of values)"),
                      br(),
                      h5("Made with ",
                         img(src = "shiny.png", height = "50"), "by ", img(src = "biodata.png", height = "30"))
                    ),
                      mainPanel(
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
                                  plotOutput("coreLineplot"),
                                  plotOutput("coreHeatmap")
                          )
                        )
                      )
                    )
            ),
          tabPanel("Dirichlet Multinomial Mixtures",
                   titlePanel("Community typing with DMM", windowTitle = "biome-shiny v0.0.0.4"),
                   sidebarPanel(
                     selectInput("datasetDMM","Choose the dataset to analyze.",
                                 choices = c("dietswap","atlas1006","peerj32")),
                     numericInput("detectionDMM", label="Input detection threshold for DMM", min = 0, max = 100, step = 0.1, value = "0.1"),
                     numericInput("prevalenceDMM", label="Input prevalence percentage of samples for DMM0", min = 0, max = 100, step = 0.1, value = "50"),
                     numericInput("maxModelsDMM", label="Input maximum number of community types for the DMM model, or leave empty for infinite types at the cost of performance.", value = "3"),
                     h5("Made with ",
                        img(src = "shiny.png", height = "50"), "by ", img(src = "biodata.png", height = "30"))
                   ),
                   mainPanel(
                      tabsetPanel(
                        tabPanel("Model Verification Lineplot",
                                 plotOutput("dmmModelCheck")
                        ),
                        
                        tabPanel("Alpha and Theta Parameters",
                                 dataTableOutput("dmmParameters"),
                                 dataTableOutput("sampleAssignments")
                        ),
                        
                        tabPanel("Taxa Contribution Per Community Type",
                                 plotOutput("taxaContributionPerComponent")
                        )
                      )
                   )
          ),
          tabPanel("Landscape analysis",
                titlePanel("Landscape analysis", windowTitle = "biome-shiny v0.0.0.4"),
                sidebarPanel(
                  selectInput("datasetLandscape","Choose the dataset to analyze.",
                              choices = c("dietswap","atlas1006","peerj32")),
                  
                  selectInput("ordinateMethodLandscape", "Choose an ordination method:", 
                              choices=c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"), selected = "CCA"),
                  
                  selectInput("ordinateDistanceLandscape", "Choose a distance method:", 
                              choices=c("bray","jaccard","unifrac"), selected = "unifrac"), #There's 44 of these in total
                  numericInput("detectionLandscape", label="Input detection threshold for landscape analysis", min = 0, max = 100, step = 0.1, value = "0.1"),
                  numericInput("prevalenceLandscape", label="Input prevalence percentage of samples for landscape analysis", min = 0, max = 100, step = 0.1, value = "50"),
                  selectInput("metadataLandscape1", "Choose first metadata column to subset data:", choices=colnames("testies")),
                  selectInput("metadataLandscapeValue1", "Choose first metadata value:", choices=sapply("testies", levels)),
                  selectInput("metadataLandscape2", "Choose second metadata column to subset data:", choices=colnames("testies")),
                  selectInput("metadataLandscapeValue2", "Choose second metadata value:", choices=sapply("testies", levels)),
                  
                  h5("Made with ",
                     img(src = "shiny.png", height = "50"), "by ", img(src = "biodata.png", height = "30"))
                ),
                mainPanel(
                  tabsetPanel(
                    tabPanel("PCA",
                      plotOutput("landscapePCA")
                    ),
                    tabPanel("PCoA/MDS/NMDS",
                       plotOutput("landscapeOrdination"),
                       plotOutput("landscapeOrdinationSamplenames")
                    ),
                    tabPanel("t-SNE",
                       plotOutput("landscapeTSne")
                    ),
                    tabPanel("Abundance histogram (one-dimensional landscape)",
                       plotOutput("landscapeAbundanceHistAbs"),
                       plotOutput("landscapeAbundanceHistRel")
                    )
                  )
                )
              )
)

#SERVER#
server <- shinyServer(function(input, output, session){
  
  #Introduction text
  output$introText <- renderText({
    paste0("biome-shiny is a microbiome analysis pipeline developed with the Shiny library for R, and based, primarily, on the \"microbiome\" and \"phyloseq\" libraries for analysis.\n\nThe app is in its earliest stages, and right now can only perform an alpha diversity, beta diversity and community composition analysis. In the future, it will hopefully include more of microbiome and phyloseq's functions, and more types of visualizations.\n\nbiome-shiny is being developed by BioData.pt for ELIXIR.")
  })
  
  #Load dataset from file
  datasetInputFromFile <- reactive({
    req(input$datase)
    import_biom(input$datase)
  })
  
  # Pick the sample dataset
  datasetInput <- reactive({ #For alpha diversity
    switch(input$dataset,
           "dietswap" = dietswap,
           "atlas1006" = atlas1006,
           "peerj32" = peerj32
          )
  })
  
  datasetInputBeta <- reactive({ #For beta diversity
    switch(input$datasetBeta,
           "dietswap" = dietswap,
           "atlas1006" = atlas1006,
           "peerj32" = peerj32
    )
  })
  
  datasetInputComposition <- reactive({ #For CC
    switch(input$datasetComp,
           "dietswap" = dietswap,
           "atlas1006" = atlas1006,
           "peerj32" = peerj32
    )
  })
  
  datasetInputCore <- reactive({ #For Core Microbiota
    switch(input$datasetCore,
           "dietswap" = dietswap,
           "atlas1006" = atlas1006,
           "peerj32" = peerj32
    )
  })

  datasetInputDMM <- reactive({ #For DMM community typing
    switch(input$datasetDMM,
           "dietswap" = dietswap,
           "atlas1006" = atlas1006,
           "peerj32" = peerj32
    )
  })
  
  datasetInputLandscape <- reactive({ #For landscape analysis
    switch(input$datasetLandscape,
           "dietswap" = dietswap,
           "atlas1006" = atlas1006,
           "peerj32" = peerj32
    )
  })
  
  
    # Metadata
  output$testies <- renderTable({
    meta(datasetInput())
  })
  
  
  #Updating SelectInputs when swapping datasets - Alpha Diversity
  observe({
      updateSelectInput(session, "x",
                        choices=colnames(meta(datasetInput())))
      updateSelectInput(session, "y",
                      choices=colnames(alpha(datasetInput())))
  })
  
  #Updating SelectInputs when swapping datasets - Beta Diversity
  observe({
    updateSelectInput(session, "xb",
                      choices=colnames(meta(datasetInputBeta())))
    updateSelectInput(session, "xc",
                      choices=colnames(meta(datasetInputBeta())))
    updateSelectInput(session, "xd",
                      choices=colnames(tax_table(datasetInputBeta())))
  })
  
  #Updating selectInputs when swapping datasets - Community Composition Analysis
  observe({
    updateSelectInput(session, "z1",
                      choices=colnames(meta(datasetInputComposition())))
    updateSelectInput(session, "z2",
                      choices=colnames(meta(datasetInputComposition())))
    updateSelectInput(session, "z3",
                      choices=colnames(meta(datasetInputComposition())))
    updateSelectInput(session, "v4",
                      choices=colnames(tax_table(datasetInputComposition())))
    })
  
  #Update metadata value selectInputs for CC Analysis
  observeEvent(input$z1,{
    updateSelectInput(session, "v1",
                      choices=(sample_data(datasetInputComposition())[[input$z1]]))
  })
  observeEvent(input$z2,{
    updateSelectInput(session, "v2",
                      choices=(sample_data(datasetInputComposition())[[input$z2]]))
  })
  observeEvent(input$z3,{
    updateSelectInput(session, "v3",
                      choices=(sample_data(datasetInputComposition())[[input$z3]]))
  })
  
  #Update metadata selectInputs when swapping datasets - Landscape analysis
  observe({
    updateSelectInput(session, "metadataLandscape1",
                      choices=colnames(meta(datasetInputLandscape())))
  })
  observe({
    updateSelectInput(session, "metadataLandscape2",
                      choices=colnames(meta(datasetInputLandscape())))
  })
  
  #Update metadata values for landscape analysis 
  observeEvent(input$metadataLandscape1,{
    updateSelectInput(session, "metadataLandscapeValue1",
                      choices=(sample_data(datasetInputComposition())[[input$metadataLandscape1]]))
  })
  observeEvent(input$metadataLandscape2,{
    updateSelectInput(session, "metadataLandscapeValue2",
                      choices=(sample_data(datasetInputComposition())[[input$metadataLandscape2]]))
  })
  
  output$summary <- renderPrint({
    summarize_phyloseq(datasetInput())
  })
  
  mergedTable <- reactive({
    merge(meta(datasetInput()),alpha(datasetInput()), all.y=TRUE )
  })

    # Renders the merged table, number of rows defined by the user  
   output$view <- DT::renderDataTable({ #Issue with search -> search by specific term (like "male") doesn't work -> causes problem when searching between female (which gives all female samples) and male (which give both male and female samples). Still, a slightly dysfunctional search function is better than none at all
   datatable(mergedTable())
  })
   
   output$violinPlot <- renderPlot({
     # Basic violin plot
     ggviolin(mergedTable(), x = input$x, y = input$y, add = "boxplot", fill = input$x, palette = c("#a6cee3", "#b2df8a", "#fdbf6f")) 
     })
    # Plot richness
    output$richnessPlot <- renderPlot({
      plot_richness(datasetInput(), x=input$x, measures=c("Shannon","Simpson"), color=input$x)
    })

   
   #BETA DIV#
   
   # -1 - Extract taxonomic ranks from the selected database
   extractTaxa <- reactive({
     colnames(tax_table(dietswap))
   })
   
   # And make an output out of it.
   output$taxRanks <- renderText({
     c(extractTaxa())
   })
   
   # 0 - Set seed for outcome replication
   sneed <- reactive({
     req(input$beta.sneed)
     set.seed(input$beta.sneed)
     paste("Seed: ", input$beta.sneed)
   })
   
   output$betaseed <- renderText({
     sneed()
   })
   
   # 1 - Convert to relative values if necessary (It says if it is in the summary)
   compositionalInput <- reactive({
     microbiome::transform(datasetInput(),"compositional")
   })
   
   # 2.1 - Pick core taxa (optional, but recommended)
   compositionalCore <- reactive({
     core(compositionalInput(), detection = .1/100, prevalence = 90/100)
   })
   
   # 2.2 - Convert core taxa to relatives
   coreConvert <- reactive({
     microbiome::transform(compositionalInput(),"compositional")
   })
   
   # 3 - Ordinate the data
   ordinateData <- reactive({
     ordinate(coreConvert(), method = input$ordinate.method, distance = input$ordinate.distance) #note to self, make method choosable from dropdown
   })
   
   # 4 - Display the plot, colored by metadata selected
   output$ordinatePlot <- renderPlot({
     plot_ordination(datasetInput(), ordinateData(), color = input$xb ) + geom_point(size = input$geom.size)
   })

   # 5 - Make split plots as well, available as metadata/metadata and metadata/taxonomy level
   # Meta/Meta
   output$metaMetaSplitOrd <- renderPlot({
      plot_ordination(datasetInput(), ordinateData(), type = "split", shape = input$xb, color = input$xc, label = input$xb)
   })
  # Meta/Tax Level
   output$metaTaxSplitOrd <- renderPlot({
     plot_ordination(datasetInput(), ordinateData(), type = "taxa", color = input$xd, label = input$xb)
   })
   
   
   # Microbial community composition #
   # 1 - Select the column 
   datasetSubsetInput <- reactive({
     subset1 <- prune_samples(sample_data(datasetInputComposition())[[input$z1]] == input$v1, datasetInputComposition())
     subset1 <- prune_samples(sample_data(subset1)[[input$z2]] == input$v2, subset1)
     subset1 <- prune_samples(sample_data(subset1)[[input$z3]] == input$v3, subset1)
     subset2 <- prune_samples(sample_data(subset1)[[input$z1]] == input$v1, subset1) %>% aggregate_taxa(level = input$v4)
     microbiome::transform(subset2, "compositional")
  })
   
   
   # 3.1 - Make plots
  output$communityPlot <- renderPlot ({
    theme_set(theme_bw(21))
    datasetSubsetInput() %>% plot_composition(sample.sort = "Bacteroidetes", otu.sort = "abundance") +
    scale_fill_manual(values = default_colors("Phylum")[taxa(datasetSubsetInput())]) 
  }) #otu.sort and and sample.sort need to be selectable
  
   # 3.2 - Make plot
  output$communityPlotGenus <- renderPlot({
  plot_composition(datasetSubsetInput(),
                        taxonomic.level = "Genus",
                        sample.sort = "nationality",
                        x.label = "nationality") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_y_percent() +
    labs(x = "Samples", y = "Relative abundance (%)",
         title = "Relative abundance data",
         subtitle = "Subtitle",
         caption = "Caption text.") + 
   theme_ipsum(grid="Y")
   })
  
  # 3.3 - Make a heatmap too
  
  #Composition barplot
  output$communityHeatmap <- renderPlot({ plot_composition(datasetSubsetInput(),
                        plot.type = "heatmap",
                        sample.sort = "neatmap",
                        otu.sort = "neatmap")
  })
  
  # 3.4 - Also a barplot, but averaged by z1 group - replace with user pick later
  output$CommunityHeatmap2 <- renderPlot({ 
    plot_composition(datasetSubsetInput(), average_by = input$z1 )
  })
  
  # 7 - And top it off with a taxa prevalence plot
  output$communityPrevalence <- renderPlot({  
    plot_taxa_prevalence(datasetSubsetInput(), input$v4 ) #Can be changed to whatever taxonomic rank the input files have
  })
  
  ## CORE MICROBIOTA ##
  #Convert to compositional
  compositionalInput2 <- reactive({
    microbiome::transform(datasetInputCore(),"compositional")
  })
  
  #Prevalence calculation
  prevalenceAbsolute <- reactive({
    as.data.frame(prevalence(compositionalInput2(), detection = input$detectionPrevalence/100, sort = TRUE, count = TRUE))
  })
  prevalenceRelative <- reactive({
    as.data.frame(prevalence(compositionalInput2(), detection = input$detectionPrevalence/100, sort = TRUE))
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
    core(compositionalInput2(), detection = input$detectionPrevalence, prevalence = input$prevalencePrevalence )
  })
  output$corePhyloSummary <- renderPrint({ # Summary of corePhylo file
    summarize_phyloseq(corePhylo())  
  })
  output$coreTaxa <- renderPrint({ # Reports the taxa in corePhylo
    taxa(corePhylo())
  })
  
  # Visualization (lineplots and heatmaps)
  output$coreLineplot <- renderPlot({
    prevalences <- as.numeric(input$prevalenceSelection)
    detections <- as.numeric(unlist(strsplit(input$detectionForLineplot, split = ",")))/100
    print(detections)
    plot_core(compositionalInput2(), prevalences = prevalences, detections = detections, plot.type = "lineplot") + xlab("Relative Abundance (%)")
  })
  
  output$coreHeatmap <- renderPlot({
    # Core with compositionals:
    prevalences <- as.numeric(input$prevalenceSelection)
    detections <- 10^seq(log10(as.numeric(input$detectionMin)), log10(as.numeric(input$detectionMax)), length = as.numeric(input$maxLength))
    gray <- gray(seq(0,1,length=5))
    plot_core(compositionalInput2(), plot.type = "heatmap", colours = gray, prevalences = prevalences, detections = detections, min.prevalence = .5) + xlab("Detection Threshold (Relative Abundance (%))")
  })
  
  #DMM COMMUNITY TYPING#
  
  compositionalInputDMM <- reactive({
    microbiome::transform(datasetInputDMM(),"compositional")
  })
  
  # Pick OTU count matrix, convert into samplex x taxa format and fit the DMM model
  dmmModelFit <- reactive({
    taxa <- core_members(compositionalInputDMM(), detection = input$detectionDMM/100, prevalence = input$prevalenceDMM/100) #need to make this related to core comp, along with the coremembers in beta diversity
    pseq <- prune_taxa(taxa, datasetInputDMM())
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
    fit[[which.min(lplc)]]
  })
  
  output$dmmParameters <- renderDT({
    mixturewt(dmmModelBestFit())
  })
  
  output$taxaContributionPerComponent <- renderPlot({
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
      print(p)
    }
  })

  #Bug: Only the final plot is rendering
  
  # LANDSCAPE ANALYSIS #
  
  # 1 - Select the column 
  datasetSubsetLandscape <- reactive({
    subset <- prune_samples(sample_data(datasetInputLandscape())[[input$metadataLandscape1]] == input$metadataLandscapeValue1, datasetInputLandscape())
    subset <- prune_samples(sample_data(subset)[[input$metadataLandscape2]] == input$metadataLandscapeValue2, subset)
    microbiome::transform(subset, "compositional")
  })
  
  #PCA plot
  
  output$landscapePCA <- renderPlot({ #Change transformation and method type to user input
    p <- plot_landscape(datasetInputLandscape(), method = "PCA", transformation = "clr") +
      labs(title = paste("PCA / CLR"))
    print(p)
  })
  
  #NMDS/PCoa/MDS
  output$landscapeOrdination <- renderPlot({
    x <- datasetInputLandscape()
    quiet(x.ord <- ordinate(x, input$ordinateMethodLandscape, input$ordinateDistanceLandscape))
    proj <- phyloseq::plot_ordination(x, x.ord, justDF=TRUE)
    names(proj)[1:2] <- paste("Comp", 1:2, sep=".")
    plot_landscape(proj[, 1:2], col = proj$nationality, legend = T)
  })
  output$landscapeOrdinationSamplenames <- renderPlot({
    x <- datasetInputLandscape()
    quiet(x.ord <- ordinate(x, input$ordinateMethodLandscape, input$ordinateDistanceLandscape))
    proj <- phyloseq::plot_ordination(x, x.ord, justDF=TRUE)
    names(proj)[1:2] <- paste("Comp", 1:2, sep=".")
    ax1 <- names(proj)[[1]]
    ax2 <- names(proj)[[2]]
    ggplot(aes_string(x = ax1, y = ax2, label = "sample"), data = proj) +
      geom_text(size = 2)
  })
  
  #
  output$landscapeTSne <- renderPlot({
    plot_landscape(datasetInputLandscape(), "t-SNE",
    distance = "euclidean", transformation = "hellinger") +
    labs(title = paste("t-SNE / Hellinger / Euclidean"))       
  })
  
  output$landscapeAbundanceHistAbs <- renderPlot({
    # Visualize population densities for specific taxa
    plot_density(datasetInputLandscape(), "Dialister") + ggtitle("Absolute abundance")
  })
  
  output$landscapeAbundanceHistRel <- renderPlot({
    # Visualize population densities for specific taxa
    x <- microbiome::transform(datasetInputLandscape(), "compositional")
    tax <- "Dialister"
    plot_density(x, tax, log10 = TRUE) +
      ggtitle("Relative abundance") +
      xlab("Relative abundance (%)")
  })
  
})
# Run the application 

shinyApp(ui = ui, server = server)