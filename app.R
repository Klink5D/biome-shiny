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

#load the example datasets
data("dietswap")
data("atlas1006")
data("peerj32")
peerj32 <- peerj32$phyloseq


#SERVER#
server <- shinyServer(function(input, output) {
  
  
  # Pick the dataset
  datasetInput <- reactive({
    switch(input$dataset,
           "dietswap" = dietswap,
           "atlas1006" = atlas1006,
           "peerj32" = peerj32
          )
  })
  
  # Metadata
  output$testies <- renderTable({
    meta(dataTable())
  })
  
  observeEvent(input$controller, {
    x <- input$controller
    
    up
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
k})

#UI#
ui <- navbarPage("meta-shiny v.0.0.0.0", fluid = TRUE,
  
  theme = shinytheme("superhero"),
                 
  tabPanel( "Alpha Diversity",
  #Page theme. Could add a selector to the sidebar instead, but for now, you get a dark theme and you'll like it.
  #Take a guess
  titlePanel("Alpha Diversity", windowTitle = "meta-shiny v.0.0.0.0"),
  
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
       img(src = "bobo.jpg", height = "100"))
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
      titlePanel("Beta Diversity", windowTitle = "meta-shiny v.0.0.0.0"),
      
      #The sidebar
      sidebarPanel(
        selectInput("dataset","Choose the dataset to analyze.", #Bug where selections shared with Alpha Div. rely on Alpha's setting to function
                    choices = c("dietswap","atlas1006","peerj32")),
        numericInput("beta.sneed","Set seed for outcome replication:",
                    value = 1, min = 0),
        selectInput("xb", "Choose a metadata column:", 
                    choices=colnames(testies), selected = "bmi_group"),
        selectInput("ordinate.method", "Choose an ordination method:", 
                    choices=c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"), selected = "CCA"),
        selectInput("ordinate.distance", "Choose a distance method:", 
                    choices=c("bray","jaccard","unifrac"), selected = "unifrac"), #There's like 40 of these things. I could add them all, but...
        sliderInput("geom.size", "Plot geometry point size:", 
                  min = 1, max = 10, step = 0.5, value = "3"),
        #Linebreaks
        br(),
        
        h5("Made by ya boi Bobo", #This'll get replaced with a proper logo.
           img(src = "bobo.jpg", height = "100"))
      ),
      
      mainPanel(
        plotOutput("ordinatePlot"),
        textOutput("betaseed")
      )
  )
)
)
# Run the application 

shinyApp(ui = ui, server = server)

