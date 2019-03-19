#####################################################################################################
# This Shiny web app is designed to work as a metagenomics pipeline                                 #
# It is based, primarily, around the microbiome package                                             #
# It was written while listening to Death Grips, the Sunless Skies OST, and Clair de Lune on a loop #
# It would behoove you to go listen to those after you're done fiddling about with the app          #
#####################################################################################################

# Load the required libraries
library("shiny")
library("shinythemes")
library("Cairo")
library("microbiome")
library("microbiomeutilities")
library("knitr")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("vegan")
library("Rtsne")
library("hrbrthemes")
library("tidyverse")
library("reshape2")

#load the datasets
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
  
  # Diversity measures
  output$diversity <- renderTable({
    alpha(dataTable())    
  })
  
  
  #Merges the alpha diversity table and the metadata of the dataset into one big table.
  mergedTable <- reactive({
    merge(meta(datasetInput()),alpha(datasetInput()), all.y=TRUE )
  })

    # Renders the merged table, number of rows defined by the user  
   output$view <- renderTable({
   head(mergedTable(), n = input$obs)
  })
   
   output$alpha_graph <- renderPlot({
     # Basic box plot
     ggviolin(testies, x = input$x, y = input$y, add = "boxplot", fill = input$fill, palette = c("#a6cee3", "#b2df8a", "#fdbf6f")) 
   })
   
})

#UI#
ui <- fluidPage(
  
  #Page theme. Could add a selector to the sidebar instead, but for now, you get a dark theme and you'll like it.
  theme = shinytheme("superhero"),
  
  #Take a guess
  titlePanel("Alpha Diversity", windowTitle = "Microbiome pipeline v.0.0.0.0"),
  
  #The sidebar
  sidebarPanel(
    selectInput("dataset","Choose the dataset to analyze",
                choices = c("dietswap","atlas1006","peerj32")),
    numericInput("obs", "Number of samples to view:", 10),
    
  #X (the metadata) and Y (the diversity measure)
    selectInput("x", "X:", 
                choices=colnames(testies)),
    selectInput("y", "Y:", 
                choices=colnames(testies)),
  #And the fill (also metadata)
    selectInput("fill", "Fill:", 
                choices=colnames(testies)),
    
    #Needs a sort by panel - with subset_samples
    
    #Linebreaks
    br(), br(),
    
    h5("Made by ya boi Bobo", #This'll get replaced with a proper logo.
       img(src = "bobo.jpg", height = "100"))
  ),
  
  mainPanel(
    tabsetPanel( type = "tabs",
                 id = "tabsetpanel",
                 
                 #Tab 1: Alpha diversity measures with metadata
                 tabPanel( title = "Table",
                           tableOutput("view")
                 ),
                 
                 #Tab 2: The exceptionally lewd Violin Plots
                 tabPanel( title = "Plot",
                           plotOutput("alpha_graph")
                 )
                 
    )
  )
)

# Run the application 

shinyApp(ui = ui, server = server)

