#####################################################################################################
# This Shiny web app is designed to work as a metagenomics pipeline                                 #
# It is based, primarily, around the microbiome package                                             #
# It was written while listening to Death Grips, the Sunless Skies OST, and Clair de Lune on a loop #
# It would behoove you to go listen to those after you're done fiddling about with the app          #
#####################################################################################################

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

#load the example datasets
data("dietswap")
data("atlas1006")
data("peerj32")
peerj32 <- peerj32$phyloseq

#UI#
ui <- navbarPage("meta-shiny v.0.0.0.1", fluid = TRUE,
                 
                 theme = shinytheme("superhero"),
                 
                 tabPanel( "Alpha Diversity",
                           #Page theme. Could add a selector to the sidebar instead, but for now, you get a dark theme and you'll like it.
                           #Take a guess
                           titlePanel("Alpha Diversity", windowTitle = "meta-shiny v.0.0.0.0«1"),
                           
                           #The sidebar
                           sidebarPanel(
                             selectInput("dataset","Choose the dataset to analyze.",
                                         choices = c("dietswap","atlas1006","peerj32")),
                             numericInput("obs", paste("Number of samples to view"), 10, 1, 222),
                             
                             #X (the metadata) and Y (the diversity measure)
                             selectInput("x", "Choose a metadata column:", 
                                         choices=colnames(testies)),
                             selectInput("y", "Choose a diversity measure:", 
                                         choices=colnames(testies)),
                             ##And the fill (also metadata)
                             # selectInput("fill", "Metadata for the color fill:", 
                             #            choices=colnames(testies)),
                             
                             #Needs a sort by panel - with subset_samples
                             
                             #Linebreaks
                             br(),
                             
                             h5("Made by ya boi Bobo", #This'll get replaced with a proper logo.
                                img(src = "bobo.png", height = "100"))
                           ),
                           
                           mainPanel(
                             tabsetPanel( type = "tabs",
                                          id = "tabsetpanel",
                                          
                                          # Tab 1: Phyloseq summary
                                          tabPanel( title = "Phyloseq Summary",
                                                    textOutput("summary")
                                          ),
                                          
                                          # Tab 2: Alpha diversity measures with metadata
                                          tabPanel( title = "Table",
                                                    tableOutput("view")
                                          ),
                                          
                                          # Tab 3: The exceptionally lewd Violin Plots
                                          tabPanel( title = "Plot",
                                                    plotOutput("alpha_graph")
                                          )
                             )
                           )
                 ),
                 tabPanel("Beta Diversity",
                          fluidPage(
                            titlePanel("Beta Diversity", windowTitle = "meta-shiny v.0.0.0.1"),
                            
                            #The sidebar
                            sidebarPanel(
                              selectInput("dataset","Choose the dataset to analyze.", #Bug where selections shared with Alpha Div. rely on Alpha's setting to function
                                          choices = c("dietswap","atlas1006","peerj32")),
                              
                              numericInput("beta.sneed","Set seed for outcome replication:",
                                           value = 1, min = 0),
                              
                              selectInput("xb", "Choose a metadata column:", 
                                          choices=colnames(testies), selected = "bmi_group"),
                              
                              selectInput("xc", "For metadata/metadata split plots, choose a metadata column:", 
                                          choices=colnames(testies), selected = "bmi_group"),
                              
                              selectInput("xd", "For metadata/tax rank split plot, choose a taxonomy rank:", 
                                          choices=c("Phylum", "Class", "Order", "Family","Genus")), # Needs to be able to populate the choices with
                              # a list obtained from database's colnames(tax_table())
                              
                              selectInput("ordinate.method", "Choose an ordination method:", 
                                          choices=c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"), selected = "CCA"),
                              
                              selectInput("ordinate.distance", "Choose a distance method:", 
                                          choices=c("bray","jaccard","unifrac"), selected = "unifrac"), #There's like 40 of these things. I could add them all, but...
                              
                              sliderInput("geom.size", "Plot geometry point size:", 
                                          min = 1, max = 10, step = 0.5, value = "3"),
                              
                              # Linebreaks
                              br(),
                              
                              h5("Made by ya boi Bobo", #This'll get replaced with a proper logo.
                                 img(src = "bobo.png", height = "100"))
                            ),
                            
                            mainPanel(
                              
                              textOutput("betaseed"),
                              
                              plotOutput("ordinatePlot", hover = "plot_hover"),
                              
                              br(),
                              
                              plotOutput("metaMetaSplitOrd", hover = "plot_hover_split"),
                              
                              br(),
                              
                              plotOutput("metaTaxSplitOrd", hover = "plot_hover_split_tax")
                            )
                          )
                 ),
                 tabPanel("Community Composition",
                          
                          fluidPage(
                            titlePanel("Community Composition", windowTitle = "meta-shiny v.0.0.0.1"),
                            
                            #The sidebar
                            sidebarPanel(
                              selectInput("datasetComp","Choose the dataset to analyze.", #Bug where selections shared with Alpha Div. rely on Alpha's setting to function
                                          choices = c("dietswap","atlas1006","peerj32")),
                              selectInput("z1", "Choose a metadata column:", # For subsetting data, #1
                                          choices=as.vector(colnames(testies)), selected = "bmi_group"),
                              selectInput("v1", "Choose a metadata value:", # For subsetting data, metadata value
                                          choices=sapply(testies, levels), selected = "lean"),
                              selectInput("z2", "Choose a seocond metadata column:", # For subsetting data, #2
                                          choices=as.character(colnames(testies)), selected = "nationality"),
                              selectInput("v2", "Choose a second metadata value:", # For subsetting data, metadata #2 value
                                          choices=sapply(testies, levels), selected = "AAM"),
                              selectInput("z3", "Choose the intended timepoint:", # For subsetting data, timepoint data
                                          choices=as.character(colnames(testies)), selected = "timepoint.within.group"),
                              selectInput("v3", "Choose a timepoint value:", # For subsetting data, timepoint data value
                                          choices=c("0","1","2","3","4","5","6","7","8","9","10"), selected = "1"),
                              selectInput("v4", "Choose a taxonomy rank:", # Tax rank to analyze
                                          choices=c("Phylum", "Class", "Order", "Family","Genus"), selected = "Phylum"),
                              
                              
                              br(),
                              
                              h5("Made by ya boi Bobo", #This'll get replaced with a proper logo.
                                 img(src = "bobo.png", height = "100"))
                            ),
                            
                            mainPanel(
                              plotOutput("communityPlot"), #Working
                              br(),
                              plotOutput("communityPlotGenus"), #Working
                              br(),
                              plotOutput("CommunityHeatmap2"), #Working
                              br(),
                              plotOutput("communityPrevalence"), #Kind of working
                              verbatimTextOutput("datasetAggregateTable")
                            )
                          )
                 )
)

#SERVER#
server <- shinyServer(function(input, output, session){
  
  
  # Pick the dataset
  datasetInput <- reactive({
    switch(input$dataset,
           "dietswap" = dietswap,
           "atlas1006" = atlas1006,
           "peerj32" = peerj32
          )
  })
  
  datasetInputComposition <- reactive({
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
  
  #Merges the alpha diversity table and the metadata of the dataset into one big table.
  output$summary <- renderPrint({
    summarize_phyloseq(datasetInput())
  })
  
  mergedTable <- reactive({
    merge(meta(datasetInput()),alpha(datasetInput()), all.y=TRUE )
  })

    # Renders the merged table, number of rows defined by the user  
   output$view <- renderTable({
   head(mergedTable(), n = input$obs)
  })
   
   output$alpha_graph <- renderPlot({
     # Basic violin plot
     ggviolin(testies, x = input$x, y = input$y, add = "boxplot", fill = input$x, palette = c("#a6cee3", "#b2df8a", "#fdbf6f")) 
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

   # 5 - Make split plots as well, available as metadata/metadata and metadata/taxonomy level (very fun part)
   # Meta/Meta
   output$metaMetaSplitOrd <- renderPlot({
      plot_ordination(datasetInput(), ordinateData(), type = "split", shape = input$xb, color = input$xc, label = input$xb)
   })
  # Meta/Tax Level
   output$metaTaxSplitOrd <- renderPlot({
     plot_ordination(datasetInput(), ordinateData(), type = "split", shape = input$xb, color = input$xd, label = input$xb)
   })
   
   # Microbial community composition #
   
   
#   # 1 - Select the column -at least the Vs work
   datasetSubsetInput <- reactive({
     subset1 <- prune_samples(sample_data(datasetInputComposition())$bmi_group == input$v1, datasetInputComposition())
     subset1 <- prune_samples(sample_data(subset1)$nationality == input$v2, subset1)
     subset1 <- prune_samples(sample_data(subset1)$timepoint.within.group == input$v3, subset1)
     subset2 <- prune_samples(sample_data(subset1)$bmi_group == input$v1, subset1) %>% aggregate_taxa(level = input$v4)
     microbiome::transform(subset2, "compositional")
  })
   
  # # 1 - Select the column
  # datasetSubsetInput <- reactive({
  #   subset1 <- subset_samples(datasetInputComposition(), bmi_group == "obese" & nationality == "AAM" & timepoint.within.group == "1") %>% microbiome::transform(transform = "compositional")
  #   subset2 <- subset_samples(subset1, bmi_group == "obese") %>% aggregate_taxa(level = input$v4)
  #   microbiome::transform(subset2, "compositional")
  #   #Very important bug: inputs aren't recognized because of subset_samples.
  # }) 
   
   output$datasetSubsetTable <- renderText({
      summarize_phyloseq(datasetSubsetInput())
      summarize_phyloseq(dietswap)
   })
   
   # 3 - Make plot
  output$communityPlot <- renderPlot ({
    theme_set(theme_bw(21))
    datasetSubsetInput() %>% plot_composition(sample.sort = "Firmicutes", otu.sort = "abundance") +
    scale_fill_manual(values = default_colors(input$v4)[taxa(datasetSubsetInput())]) 
  }) #otu.sort and and sample.sort need to be selectable
  
   # 4 - Make more plot
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
  
  # 6 - Make a heatmap too
  
  #Composition heatmap
  output$communityHeatmap <- renderPlot({ plot_composition(datasetSubsetInput(),
                        plot.type = "heatmap",
                        sample.sort = "neatmap",
                        otu.sort = "neatmap")
  })
  
  # 6.1 - Also a heatmap, but averaged by group
  output$CommunityHeatmap2 <- renderPlot({ 
    plot_composition(datasetSubsetInput(), average_by = "bmi_group")
  })
  
  # 7 - And top it off with a taxa prevalence plot
  output$communityPrevalence <- renderPlot({  
    plot_taxa_prevalence(datasetSubsetInput(), "Phylum") #Can be changed to whatever taxonomic rank the input files have
  })
})
# Run the application 

shinyApp(ui = ui, server = server)


#Microbiome Composition#
#Sample subset
#micro_composition <- subset_samples(alpha0, group == "DI" & nationality == "AFR" & timepoint.within.group == 1)

#data(dietswap)
#pseq3 <- dietswap %>% subset_samples(nationality == "AFR") %>% aggregate_taxa(level = "Phylum") %>%  microbiome::transform(transform = "compositional")

#Compositional barplot of analyzed data
#theme_set(theme_bw(21))
#p <- pseq3 %>% plot_composition(sample.sort = "Firmicutes", otu.sort = "abundance") +
  
#  # Set custom colors
#  scale_fill_manual(values = default_colors("Phylum")[taxa(pseq3)]) 

#Take a guess
#print(p)

# Limit the analysis on core taxa and specific sample group
#q <- plot_composition(alpha0,
#                      taxonomic.level = "Genus",
#                      sample.sort = "nationality",
#                      x.label = "nationality") +
#  guides(fill = guide_legend(ncol = 1)) +
#  scale_y_percent() +
#  labs(x = "Samples", y = "Relative abundance (%)",
#       title = "Relative abundance data",
#       subtitle = "Subtitle",
#       caption = "Caption text.") + 
 # theme_ipsum(grid="Y")
#print(q)  


# Averaged by group
#s <- plot_composition(alpha0,
#                      average_by = "bmi_group", transform = "compositional")
#s <- plot_composition(microbiome::transform(pseq, "compositional"),
#                      plot.type = "heatmap",
#                      sample.sort = "neatmap",
#                      otu.sort = "neatmap")
#print(s)

#Composition heatmap
#t <- plot_composition(microbiome::transform(pseq3, "compositional"),
#                      plot.type = "heatmap",
#                      sample.sort = "neatmap",
#                      otu.sort = "neatmap")
#print(t)