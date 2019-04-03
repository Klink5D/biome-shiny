# biome-shiny: A Shiny R app for metagenomics analysis, built around the "microbiome" package
# current version: v0.0.0.2

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
library("randomcoloR")

#load the example datasets
data("dietswap")
data("atlas1006")
data("peerj32")
peerj32 <- peerj32$phyloseq

#UI#
ui <- navbarPage("biome-shiny v0.0.0.2", fluid = TRUE,
                 
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
                           titlePanel("Alpha Diversity", windowTitle = "biome-shiny v0.0.0.2"), 

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
                            titlePanel("Beta Diversity", windowTitle = "biome-shiny v0.0.0.2"),
                            
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
                            titlePanel("Community Composition", windowTitle = "biome-shiny v0.0.0.2"),
                            
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
                    titlePanel("Core microbiota analysis", windowTitle = "biome-shiny v0.0.0.2"),
                    sidebarPanel(
                      #selectInput("datasetCore", "Choose the dataset to analyze"),
                      sliderInput("detectionPrevalence", "For prevalences: Choose detection value", min = 0.00, max = 1, value = 0.01, step = 0.01),
                      numericInput("prevalencePrevalence","Input a prevalence value", min = 0, max = 1, value = 0.5, step = 0.05),
                      br(),
                      h5("Made with ",
                         img(src = "shiny.png", height = "50"), "by ", img(src = "biodata.png", height = "30"))
                    ),
                      mainPanel(
                        tabsetPanel(
                          tabPanel("Prevalence",
                            # Output the prevalence in relatives
                            verbatimTextOutput("prevalenceRelative"),
                            # And in absolutes
                            verbatimTextOutput("prevalenceAbsolute")
                          ),
                          tabPanel("Core taxa")
                        )
                      )
                    )
            )
              
)

#SERVER#
server <- shinyServer(function(input, output, session){
  
  #Introduction text
  output$introText <- renderText({
    paste0("biome-shiny is a metagenomics pipeline developed with the Shiny library for R, and based, primarily, on the \"microbiome\" and \"phyloseq\" libraries for analysis.\n\nThe app is in its earliest stages, and right now can only perform an alpha diversity, beta diversity and community composition analysis. In the future, it will hopefully include more of microbiome and phyloseq's functions, and more types of visualizations.\n\nbiome-shiny is being developed for BioData.pt and ELIXIR.")
  })
  
  #Load dataset from file
  datasetInputFromFile <- reactive({
    req(input$datase)
    input$datase
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
   # 1 - Sebolect the column 
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
    microbiome::transform(datasetInput(),"compositional")
  })
  
  # Prevalences
  prevalenceAbsolute <- renderText({
    prevalence(compositionalInput2(), detection = input$detectionPrevalence, sort = TRUE, count = TRUE)
  })
  
  prevalenceRelative <- renderText({
    relativeData <- microbiome::transform(datasetInput(), "compositional")
    prevalence(compositionalInput2(), detection = input$detectionPrevalence, sort = TRUE)
  })
  
  # Core analysis
  
  coreMembers <- renderText({
    core_members(compositionalInput2(), detection = input$detectionPrevalence, prevalence = input$prevalencePrevalence )
  })
  
  corePhylo <- reactive({
    core(compositionalInput2(), detection = input$detectionPrevalence, prevalence = input$prevalencePrevalence )
  })
  
  coreTaxa <- renderText({
    taxa(corePhylo())
  })
  
  
  
  # Core abundance per sample
  
  # Visualization (lineplots and heatmaps)
  
})
# Run the application 

shinyApp(ui = ui, server = server)