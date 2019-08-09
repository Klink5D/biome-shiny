# Biome-shiny 0.7

library(shiny)
library(shinydashboard)
library(microbiome)
library(phyloseq)
library(rmarkdown)
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
library(biomformat)
library(ggplotify)
#library(limma)

#Plot_ordered_bar function | Function created by pjames1 @ https://github.com/pjames1 | All credit to him
plot_ordered_bar<-function (physeq, x = "Sample", 
                            y = "Abundance", 
                            fill = NULL, 
                            leg_size = 0.5,
                            title = NULL) {
  require(ggplot2)
  require(phyloseq)
  require(plyr)
  require(grid)
  bb <- psmelt(physeq)
  
  
  samp_names <- aggregate(bb$Abundance, by=list(bb$Sample), FUN=sum)[,1]
  .e <- environment()
  bb[,fill]<- factor(bb[,fill], rev(sort(unique(bb[,fill])))) #fill to genus
  
  
  bb<- bb[order(bb[,fill]),] # genus to fill
  p = ggplot(bb, aes_string(x = x, y = y, 
                            fill = fill), 
             environment = .e, ordered = FALSE)
  
  
  p = p +geom_bar(stat = "identity", 
                  position = "stack", 
                  color = "black") 
  
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  
  p = p + guides(fill = guide_legend(override.aes = list(colour = NULL), reverse=TRUE)) + 
    theme(legend.key = element_rect(colour = "black")) 
  
  p = p + theme(legend.key.size = unit(leg_size, "cm"))
  
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

# Functions to dynamically generate chunks for the final report 
tidy_function_body <- function(fun) {
  paste(tidy_source(text = as.character(body(fun))[-1])$text.tidy, collapse="\n")
}

make_chunk_from_function_body <- function(fun, chunk.name="", chunk.options=list()) {
  opts <- paste(paste(names(chunk.options), chunk.options, sep="="), collapse=", ")
  header <- paste0("```{r ", chunk.name, " ", chunk.options, "}")
  paste(header, tidy_function_body(fun), "```", sep="\n")
}

report.source <- reactive({
  req(sessionData$import.params(),
      sessionData$filter.params())
  
  report <- readLines("sc_report_base.Rmd")
  
  insert.function <- function(report, tag, fun, chunk.name = "", chunk.options = list()) {
    w <- which(report == tag)
    report[w] <- make_chunk_from_function_body(fun, chunk.name = chunk.name, chunk.options = chunk.options)
    
    return(report)
  }
  
  # Import
  report <- insert.function(report, "<!-- import.fun -->", sessionData$import.fun(), chunk.name = "import")
  
  # Filter
  report <- insert.function(report, "<!-- filter.fun -->", sessionData$filter.fun(), chunk.name = "filter")
  
  
  return(report)
})


# Load sample datasets #
data("dietswap")
data("atlas1006")
data("peerj32")
peerj32 <- peerj32$phyloseq

# UI
ui <- dashboardPage(
  dashboardHeader(title = "biome-shiny v0.7"),
  dashboardSidebar(
    sidebarMenu(
      menuItem(
        "Introduction | Data Upload",
        tabName = "intro"
      ),
      menuItem("Phyloseq Summary", tabName="phyloseqsummary"),
      br(),
      paste0("Microbiome analysis"),
      menuItem("Core microbiota", tabName = "coremicrobiota"),
      menuItem("Community composition", tabName = "communitycomposition"),
      menuItem("Alpha diversity", tabName = "alphadiversity"),
      menuItem("Beta diversity", tabName = "betadiversity"),
      #menuItem("Density analysis", tabName = "landscaping"),
      #menuItem("DMM Clustering", tabName = "dirichlet"),
      br(),
      paste0("Statistical analysis"),
      menuItem("PERMANOVA", tabName = "permanova"),
      # menuItem("ANOSIM", tabName = "anosim"),
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
      h1("Introduction and Data Upload"),
      br(),
      box(
        width = "3",
        radioButtons(
          "datasetChoice",
          "Dataset Upload",
          c("Upload dataset", "Use sample dataset"),
          selected = "Use sample dataset"
        ),
        conditionalPanel( condition = "input.datasetChoice == 'Upload dataset'",
                          radioButtons("datasetType", "Select dataset characteristics:", c(".biom file including sample variables",".biom file with .csv mapping file",".biom file without .csv mapping file"))                  
        ),
        conditionalPanel(
          condition = "input.datasetType == '.biom file including sample variables'",
          fileInput(
            "dataset",
            "Dataset:",
            multiple = FALSE,
            accept = c(".biom"), placeholder="Phyloseq .biom files"
          )
        ),
        conditionalPanel(
          condition = "input.datasetType == '.biom file without .csv mapping file'",
          fileInput(
            "dataset",
            "Dataset:",
            multiple = FALSE,
            accept = c(".biom"), placeholder="Phyloseq .biom files"
          ),
          checkboxInput("samplesAreColumns","OTU Table: Samples are columns", value = FALSE)
        ),
        conditionalPanel(
          condition = "input.datasetType == '.biom file with .csv mapping file'",
          fileInput(
            "dataset",
            "Dataset:",
            multiple = FALSE,
            accept = c(".biom"), placeholder="Phyloseq .biom files"
          ),
          fileInput("datasetMetadata", ".csv mapping file (sample variables):",
                    multiple = FALSE,
                    accept = c(".csv"), placeholder=".csv files"
          )
        ),
        conditionalPanel(
          condition = "input.datasetChoice == 'Use sample dataset'",
          selectInput(
            "datasetSample",
            "Or choose a sample dataset:",
            choices = c("dietswap", "atlas1006", "peerj32")
          )
        ),
        actionButton("datasetUpdate", "Update Dataset")
      ),
      box(
        paste0(
          "biome-shiny is a microbiome analysis pipeline developed with the Shiny library for R, and based, primarily, on the \"microbiome\" and \"phyloseq\" libraries for analysis.\n\n\n\nThe application takes a .biom file, generated by programs such as QIIME, as an input. If necessary, it is possible to upload a .csv file containing the dataset's sample data. Finally, if the user does not happen to have any sample data, the application can generate sample data out of the sample headers. For more information on the .biom file format please visit the following link: http://biom-format.org/ "
        )
      )
    ),
    
    tabItem(tabName = "phyloseqsummary",
            h1("Phyloseq Summary"),
            verbatimTextOutput("summary")
    ),
    
    
    # Core microbiota #
    tabItem(
      tabName = "coremicrobiota",
      tabsetPanel(
        tabPanel("Core Taxa Summary",
                 tabsetPanel(
                   tabPanel("Variables",
                            box( title = "Variables", collapsible = TRUE,
                              numericInput("detectionPrevalence2", "Detection (Relative Abundance %):", min = 0.00, max = 100, value = 0.01, step = 0.01),
                              numericInput("prevalencePrevalence","Prevalence:", min = 0, max = 1, value = 0.5, step = 0.05),
                              checkboxInput("coreFilterDataset", "Set as active dataset", value = FALSE)
                            )
                   ),
                   tabPanel("Absolute prevalence", dataTableOutput("prevalenceAbsoluteOutput")),
                   tabPanel("Relative prevalence", dataTableOutput("prevalenceRelativeOutput")),
                   tabPanel("Summary", verbatimTextOutput("corePhyloSummary")),
                   tabPanel("Taxa", verbatimTextOutput("coreTaxa"))
                 )
        ),
        tabPanel("Core Taxa Visualization",
                   fixedRow(
                            box( width = "2", collapsible = TRUE,
                                textInput("detectionMin", label = "Minimum detection threshold (Relative Abundance(%))", value = "0.0000001"),
                                checkboxInput("transparentCoreHeatmap", "Transparent background", value = TRUE)
                              ),
                            box(plotlyOutput("coreHeatmap", width = "1000px", height = "500px"))
                   )
        )
      )
    ),
    
    # Community composition #
    tabItem(
      tabName = "communitycomposition",
      tabsetPanel(
        tabPanel(title = "Abundance in samples by taxa", #Absolute abundance/counts
                 tabsetPanel(       
                   tabPanel("Variables",
                            box(
                              width = "2",
                              collapsible = TRUE,
                              selectInput(
                                "z1",
                                "Sample variable:",
                                choices = colnames("datasetMetadata"),
                                selected = "sample"
                              ),
                              checkboxInput("transparentCommunityPlot", "Transparent background", value = TRUE),
                              checkboxInput("communityPlotFacetWrap", "Group samples by metadata variable (facet_wrap)", value = FALSE),
                              conditionalPanel(condition = "input.communityPlotFacetWrap == 1",
                                               selectInput(
                                                 "z2",
                                                 "Metadata:",
                                                 choices = colnames("datasetMetadata"),
                                                 selected = "nationality"
                                               ) 
                              )
                            ),
                            box("Taxonomy", width = "2", collapsible = TRUE,
                                selectInput(
                                  "v4",
                                  "Choose a taxonomy rank:",
                                  # Tax rank to analyze
                                  choices = c("Phylum", "Class", "Order", "Family", "Genus"),
                                  selected = "Phylum"
                                )
                            )
                   ),
                   tabPanel("Absolute Abundance Plot", plotlyOutput("communityPlot")),
                   tabPanel("Relative Abundance Plot", plotlyOutput("communityPlotGenus"))
                 )
        ),
        tabPanel(title = "Taxonomy Prevalence Plot",
                 plotlyOutput("communityPrevalence")
        )
      )
    ),
    
    # Alpha Diversity #
    tabItem(
      tabName = "alphadiversity",
      br(),
      tabsetPanel(
        type = "tabs",
        id = "tabsetpanel",
        
        tabPanel(title = "Evenness Table",
                 dataTableOutput("evennessTable")),
        tabPanel(title = "Abundance Table",
                 tabsetPanel(
                   tabPanel( title = "Abundance (Counts)",
                             dataTableOutput("absoluteAbundanceTable")  
                   ),
                   tabPanel( title = "Abundance (%)",
                             dataTableOutput("relativeAbundanceTable")
                   )
                 )
        ),        
        
        # Alpha diversity measures with metadata
        tabPanel(title = "Metadata Table with diversity measures",
                 DT::dataTableOutput("view")),
        
        # A phyloseq richness plot
        tabPanel(title = "Richness Plot",
                 tabsetPanel(
                   tabPanel( title = "Variables",
                             box(
                               width = "2",
                               title = "Variables",
                               # #X (the metadata) and Y (the diversity measure)
                               selectInput(
                                 "x2",
                                 "Samples:",
                                 choices = colnames("datasetMetadata")
                               ),
                               selectInput("x3", "Point color metadata:", choices = colnames("datasetMetadata"), selected = "subject"),
                               checkboxInput("transparentRichness", "Transparent background", value = TRUE),
                               checkboxInput("richnessPlotGridWrap", "Sort samples by metadata variable", value = FALSE),
                               conditionalPanel(condition = "input.richnessPlotGridWrap == 1",
                                                selectInput("x", "Sample sorting metadata:", choices = colnames("datasetMetadata"), selected = "nationality") 
                               )
                             ),
                             box(radioButtons("richnessChoices", "Choose diversity measures" ,choices = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), selected = "Shannon"))),
                   tabPanel( title = "Plot",
                             plotlyOutput("richnessPlot"))))
      )
    ),
    
    #Beta diversity#
    tabItem(
      tabName = "betadiversity",
      tabsetPanel(
        tabPanel(title = "Ordination Plot",
                 tabsetPanel(
                   tabPanel(title = "Variables",
                            box(
                              title = "Variables",
                              width = "2",
                              collapsible = TRUE,
                              collapsed = FALSE,
                              selectInput(
                                "xb", "Sample variable:", choices = colnames("datasetMetadata"), selected = "bmi_group"
                              ),
                              selectInput(
                                "ordinate.method",
                                "Ordination method:",
                                choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                selected = "CCA"
                              ),
                              selectInput(
                                "ordinate.distance",
                                "Distance:",
                                choices = c("bray", "jaccard", "unifrac"),
                                selected = "jacard"
                              ),
                              sliderInput(
                                "geom.size",
                                "Plot geometry point size:",
                                min = 1,
                                max = 10,
                                step = 0.5,
                                value = "3"
                              ),
                              checkboxInput("transparentOrdinatePlot", "Transparent background", value = TRUE)
                            )),
                   tabPanel(title = "Plot",
                            plotlyOutput("ordinatePlot"),
                            textOutput("ordinatePrint"))
                 )),
        
        tabPanel(title = "Split Ordination Plot",
                 tabsetPanel(
                   tabPanel(title = "Variables",
                            box(
                              title = "Variables",
                              width = "2",
                              collapsible = TRUE,
                              collapsed = FALSE,
                              selectInput(
                                "xb2", "Sample variable:", choices = colnames("datasetMetadata"), selected = "bmi_group"
                              ),
                              
                              selectInput(
                                "zbsplit",
                                "Taxonomy rank:",
                                choices = c("Phylum", "Class", "Order", "Family", "Genus")
                              ),
                              
                              selectInput(
                                "ordinate.method2",
                                "Ordination method:",
                                choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                selected = "CCA"
                              ),
                              
                              selectInput(
                                "ordinate.distance2",
                                "Distance:",
                                choices = c("bray", "jaccard", "unifrac"),
                                selected = "unifrac"
                              ),
                              
                              sliderInput(
                                "geom.size2",
                                "Point size:",
                                min = 1,
                                max = 10,
                                step = 0.5,
                                value = "3"
                              ),
                              checkboxInput("transparentSplitOrd", "Transparent background", value = TRUE)
                            )),
                   tabPanel( title = "Plot", plotlyOutput("splitOrd")))),
        
        tabPanel(title = "Taxa Plot",
                 tabsetPanel(
                   tabPanel( title = "Variables",
                             box(
                               title = "Variables",
                               width = "2",
                               collapsible = TRUE,
                               collapsed = FALSE,
                               # selectInput(
                               #   "xb3", "Choose a metadata column:", choices = colnames("datasetMetadata"), selected = "bmi_group"
                               # ),
                               selectInput(
                                 "zb",
                                 "Taxonomy rank:",
                                 choices = c("Phylum", "Class", "Order", "Family", "Genus")
                               ),
                               selectInput(
                                 "ordinate.method3",
                                 "Ordination method:",
                                 choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                 selected = "CCA"
                               ),
                               selectInput(
                                 "ordinate.distance3",
                                 "Distance:",
                                 choices = c("bray", "jaccard", "unifrac"),
                                 selected = "unifrac"
                               ),
                               sliderInput(
                                 "geom.size3",
                                 "Plot geometry point size:",
                                 min = 1,
                                 max = 10,
                                 step = 0.5,
                                 value = "3"
                               ),
                               checkboxInput("transparentTaxaOrd", "Transparent background", value = TRUE)
                             )),
                   tabPanel(title = "Plot", plotlyOutput("taxaOrd"))
                 )
        )
      )),
    
    tabItem(
      tabName = "landscaping",
      tabsetPanel(
        tabPanel("PCA",
                 plotlyOutput("landscapePCA")
        ),
        tabPanel("PCoA/MDS/NMDS",
                 tabsetPanel(
                   tabPanel("Variables",
                            box("Data Subsetting", width = "2", collapsible = TRUE,
                                numericInput("detectionLandscape", label="Input detection threshold for landscape analysis", min = 0, max = 100, step = 0.1, value = "0.1"),
                                numericInput("prevalenceLandscape", label="Input prevalence percentage of samples for landscape analysis", min = 0, max = 100, step = 0.1, value = "50")
                            ),
                            box("Data Subsetting", width = "2", collapsible = TRUE,
                                selectInput("metadataLandscape1", "Choose first metadata column to subset data:", choices=colnames("datasetMetadata")),
                                selectInput("metadataLandscapeValue1", "Choose first metadata value:", choices=sapply("datasetMetadata", levels)),
                                selectInput("metadataLandscape2", "Choose second metadata column to subset data:", choices=colnames("datasetMetadata")),
                                selectInput("metadataLandscapeValue2", "Choose second metadata value:", choices=sapply("datasetMetadata", levels))                    
                            ),
                            box(title = "Variables", width = "2", collapsible = TRUE,
                                selectInput("ordinateMethodLandscape", "Choose an ordination method:", 
                                            choices=c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"), selected = "CCA"),
                                selectInput("ordinateDistanceLandscape", "Choose a distance method:", 
                                            choices=c("bray","jaccard","unifrac"), selected = "unifrac")
                            )
                   ),
                   tabPanel("Plot", plotlyOutput("landscapeOrdination")
                            #,plotlyOutput("landscapeOrdinationSamplenames")
                   )
                 )
        ),
        tabPanel("t-SNE",
                 plotlyOutput("landscapeTSne")
        ),
        tabPanel("Abundance histogram",
                 tabsetPanel(
                   tabPanel("Taxa",
                            box(
                              selectInput("landscapeOTU", "Choose OTU:", choices = "")
                            )
                   ),
                   tabPanel("Absolute Abundance",
                            plotlyOutput("landscapeAbundanceHistAbs")
                   ),
                   tabPanel("Relative Abundance",
                            plotlyOutput("landscapeAbundanceHistRel")
                   )
                 )
        )
      )
      
    ),
    tabItem(
      tabName = "dirichlet",
      tabsetPanel(
        tabPanel("Variables",
                 box( title = "Variables", width = "2", collapsible = TRUE,
                      numericInput("detectionDMM", label="Input detection threshold for DMM", min = 0, max = 100, step = 0.1, value = "0.1"),
                      numericInput("prevalenceDMM", label="Input prevalence percentage of samples for DMM0", min = 0, max = 100, step = 0.1, value = "50"),
                      numericInput("maxModelsDMM", label="Input maximum number of community types for the DMM model.", value = "3")
                 )
        ),
        
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
      tabsetPanel(
        tabPanel( title = "P-Value",
                  tabsetPanel(
                    tabPanel(title = "Variables",
                             box( title = "Variables", width= "2", collapsible = TRUE,
                                  selectInput("permanovaDistanceMethodP","Dissimilarity index:", choices = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis"), selected = "bray"),
                                  selectInput("permanovaColumnP","Sample variable:", choices = colnames("datasetMetadata")),
                                  numericInput("permanovaPermutationsP", "Number of permutations:", min = 1, step = 1, value = 99)
                             )              
                    ),
                    tabPanel(title = "Data Tables",
                             h2("P-Value table"),
                             dataTableOutput("pValue"),
                             h2("Homogeniety table"),
                             dataTableOutput("homogeniety")
                    )
                  )
        ),
        tabPanel ( title = "Top Factors",
                   tabsetPanel(
                     tabPanel(title = "Variables",
                              box( title = "Variables", width= "2", collapsible = TRUE,
                                   selectInput("permanovaDistanceMethodFac","Distance method:", choices = c("bray","jacard","unifrac"), selected = "bray"),
                                   selectInput("permanovaColumnFac","Sample variable:", choices = colnames("datasetMetadata")),
                                   numericInput("permanovaPermutationsFac", "Number of permutations:", min = 1, step = 1, value = 99)
                              )              
                     ),
                     tabPanel(title = "Plot",
                              plotlyOutput("topFactorPlot")
                     )
                   )
        ),
        tabPanel ( title = "Network Plot",
                   tabsetPanel(
                     tabPanel(title = "Variables",
                              box( title = "Variables", width= "2", collapsible = TRUE,
                                   selectInput("permanovaPlotTypeNet", "Network Plot Type:", c("Sample Network Plot" = "samples", "Taxa/OTU Network Plot" = "taxa"), selected = "samples"),
                                   selectInput("permanovaDistanceMethodNet","Distance method (not required by all ordination methods):", choices = c("bray","jacard","unifrac"), selected = "bray"),
                                   selectInput("permanovaMethodNet","Ordination method:",
                                               choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                               selected = "CCA"),
                                   # sliderInput("permanovaPlotSizeNet", "Plot point size:", min = 0.5, max = 10, step = 0.5, value = "3"),
                                   numericInput("permanovaPermutationsNet", "Number of permutations:", min = 1, step = 1, value = 99),
                                   conditionalPanel(condition = "input.permanovaPlotTypeNet = 'samples'",
                                      selectInput("permanovaMetadataNet", "Sample variable to cluster data samples:", c("Update")),
                                      selectInput("permanovaMetaShapeNet", "Sample variable to set different point shapes:", c("Update"))
                                    )         
                                  )
                              )
                     ),
                     tabPanel(title ="Plot",
                              plotlyOutput("netPlot")
                     )
                   )
        )
      ),
    tabItem( tabName = "results",
             tabsetPanel(
               tabPanel("Core Microbiota",
                        box(
                          radioButtons("renderPrevTables", "Render Prevalence Tables", choices = c("Yes", "No")),
                          radioButtons("renderCoreTaxaSummary", "Render Core Taxa Summary", choices = c("Yes", "No")),
                          radioButtons("renderCoreHeatmap", "Render Core Taxa Heatmap", choices = c("Yes", "No"))
                        )
               ),            
               tabPanel("Community Composition",
                        box(
                          radioButtons("renderSampleAbundancePlot", "Render Abundance in Samples Plot", choices = c("Yes", "No")),
                          radioButtons("renderRelativeAbundancePlot", "Render Relative Abundance in Samples Plot", choices = c("Yes", "No")),
                          radioButtons("renderPrevalencePlot", "Render Taxa Prevalence Plot", choices = c("Yes", "No"))
                        )
               ),
               tabPanel("Alpha Diversity",
                        box(
                          radioButtons("renderRichness", "Render Richness Plot", choices = c("Yes", "No")),
                          radioButtons("printTableAlpha", "Print Table Head", choices = c("Yes", "No"))
                        )
               ),
               tabPanel("Beta Diversity",
                        box(
                          radioButtons("renderOrdplot", "Render Ordination Plot", choices = c("Yes", "No")), 
                          radioButtons("renderSplitOrdplot", "Render Split Ordination Plot", choices = c("Yes", "No")),
                          radioButtons("renderTaxaplot", "Render Taxa Plot", choices = c("Yes", "No"))
                        )
               ),
               
               tabPanel("PERMANOVA",
                        box(
                          radioButtons("renderDensityPlot", "Render density plot", choices = c("Yes","No")),
                          radioButtons("renderPValueTables", "Render P-Value tables", choices = c("Yes","No")),
                          radioButtons("renderFactorPlot", "Render top factors barplot", choices = c("Yes","No")),
                          radioButtons("renderNetworkMap", "Render network map", choices = c("Yes","No"))
                        )
               )
             ),
             fluidRow( 
               box( radioButtons('format', 'Document format (HTML only for now)', c('HTML'), inline = TRUE, selected = 'HTML'),
                    downloadButton('downloadReportAlpha', label = "Download report")
               )
             )  
          )   
    )
  )
)



# Server
server <- function(input, output, session) {
  datasetChoice <- reactive({
    if (input$datasetChoice == "Use sample dataset") {
      switch(
        input$datasetSample,
        "dietswap" = dietswap,
        "atlas1006" = atlas1006,
        "peerj32" = peerj32
      )
    } else {
      req(input$dataset)
      if(input$datasetType == ".biom file including sample variables") { #Simple .biom upload with sample_data() already set
        datapath <- input$dataset$datapath
        biomfile <- import_biom(datapath)
        return(biomfile)
      }
      if(input$datasetType == ".biom file with .csv mapping file"){ #Loads a .csv along with the .biom
        datapath <- input$dataset$datapath
        datapathMetadata <- input$datasetMetadata$datapath
        a <- import_biom(datapath)
        b <- sample_data(as.data.frame(read.csv(datapathMetadata, skipNul = TRUE)))
        biomfile <- merge_phyloseq(a,b)
        return(biomfile)
      }
      if(input$datasetType == ".biom file without .csv mapping file"){ #Loads a .biom file and generates sample metadata
          datapath <- input$dataset$datapath
          a <- import_biom(datapath)
          if(input$samplesAreColumns == TRUE){ 
            samples.out <- colnames(otu_table(a))
          }
          if(input$samplesAreColumns == FALSE){
            samples.out <- rownames(otu_table(a))
          }
          subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
          samdf <- data.frame(Subject=subject)
          rownames(samdf) <- samples.out
          b <- sample_data(samdf)
          biomfile <- merge_phyloseq(a, b)
          return(biomfile)
        }
      }
    }  
  )
  
  # New DatasetInput function works as an intermediary that checks if the dataset has been altered
  datasetInput <- reactive({
    if(input$coreFilterDataset == TRUE ){ # Filters the dataset
      dataset <- corePhylo()
    }
    if(input$coreFilterDataset == FALSE ) { # Standard dataset input without filtering applied
      dataset <- datasetChoice()
    }
    return(dataset)
  })
  
  
  ## Core Microbiota ##
  #Filtering the dataset
  prevalenceAbsolute <- reactive({
    as.data.frame(prevalence(compositionalInput(), detection = input$detectionPrevalence2/100, sort = TRUE, count = TRUE))
  })
  prevalenceRelative <- reactive({
    as.data.frame(prevalence(compositionalInput(), detection = input$detectionPrevalence2/100, sort = TRUE))
  })
  
  # Produce phyloseq file with core OTUs only
  corePhylo <- reactive({
    core(datasetChoice(), detection = input$detectionPrevalence2, prevalence = input$prevalencePrevalence )
  })
  
  output$corePhyloSummary <- renderPrint({ # Summary of corePhylo file
    summarize_phyloseq(corePhylo())  
  })
  output$coreTaxa <- renderPrint({ # Reports the taxa in corePhylo
    taxa(corePhylo())
  })
  
  output$prevalenceAbsoluteOutput <- renderDT({
    datatable(prevalenceAbsolute())
  })
  output$prevalenceRelativeOutput <- renderDT({
    datatable(prevalenceRelative())
  })
  
  # Visualization (lineplots and heatmaps)
  coreLineplotParams <- reactive({
    prevalences <- as.numeric(input$prevalenceSelection)
    detections <- as.numeric(unlist(strsplit(input$detectionForLineplot, split = ",")))/100
    print(detections)
    lineplot <- plot_core(compositionalInput(), prevalences = prevalences, detections = detections, plot.type = "lineplot") + xlab("Relative Abundance (%)")
    plotly_build(lineplot)
  })
  output$coreLineplot <- renderPlotly({
    coreLineplotParams()
  })
  coreHeatmapParams <- reactive({
    # Core with compositionals:
    #prevalences <- as.numeric(input$prevalenceSelectionHeat)
    detections <- 10^seq(log10(as.numeric(input$detectionMin)), log10(1), length = 10)
    gray <- heat.colors(3)
    coreplot <- plot_core(compositionalInput(), plot.type = "heatmap", colours = gray, prevalences = 0, detections = detections) + xlab("Detection Threshold (Relative Abundance (%))")
    if(input$transparentCoreHeatmap == TRUE){
      coreplot <- coreplot +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    ggplotly(coreplot, height = 500, width = 900)
  })
  output$coreHeatmap <- renderPlotly({
    coreHeatmapParams()
  })
  
  ## Community Composition ##
  
  # Updating SelectInputs when database changes #
  observeEvent(input$datasetUpdate, {
    updateSelectInput(session, "z1",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "z2",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "z3",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "v4",
                      choices = colnames(tax_table(datasetInput())))
    updateSelectInput(session, "z1Average",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "v4Plot",
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
  
  # Abundance of taxa in sample variable by taxa
  communityPlotParams <- reactive ({
    if(input$communityPlotFacetWrap == FALSE){
      compositionplot <- plot_ordered_bar(datasetInput(), x=input$z1, y="Abundance", fill=input$v4, title=paste0("Abundance by ", input$v4, " in ", input$z1))  + geom_bar(stat="identity") + theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90)
    } else {
      compositionplot <- plot_ordered_bar(datasetInput(), x=input$z1, y="Abundance", fill=input$v4, title=paste0("Abundance by ", input$v4, " in ", input$z1))  + geom_bar(stat="identity") + facet_grid(paste('~',input$z2), scales = "free", space = "free") + theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90)
    }
    if(input$transparentCommunityPlot == TRUE){
      compositionplot <- compositionplot + 
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    ggplotly(compositionplot, height = 500, width = 1060)
  })
  output$communityPlot <- renderPlotly({
    communityPlotParams()  
  })
  
  communityPlotGenusParams <- reactive({
    if(input$communityPlotFacetWrap == FALSE){
      compositionplot <- plot_ordered_bar(compositionalInput(), x=input$z1, fill=input$v4, title=paste0("Relative abundance by ", input$v4, " in ", input$z1))  + geom_bar(stat="identity") +
        guides(fill = guide_legend(ncol = 1)) +
        scale_y_percent() +
        theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90)
    } else {
      compositionplot <- plot_ordered_bar(compositionalInput(), x="Sample",  fill=input$v4, title=paste0("Relative abundance by ", input$v4, " in ", input$z1))  + geom_bar(stat="identity") +
        guides(fill = guide_legend(ncol = 1)) +
        scale_y_percent() +
        theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90) + facet_grid(paste('~',input$z2),scales = "free", space = "free")
    }
    if(input$transparentCommunityPlot == TRUE){
      compositionplot <- compositionplot + 
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    ggplotly(compositionplot, height = 500, width = 1050)
  })
  output$communityPlotGenus <- renderPlotly({
    communityPlotGenusParams()
  })
  
  # Barplot, averaged by z1Average group
  communityBarplotParams <- reactive({
    compplot <- plot_composition(datasetSubsetInput(), average_by = input$z1Average, plot.type = input$plotTypeZ1Average )
    print(compplot)
    #plotly_build(compplot)
  })
  output$communityBarplot <- renderPlotly({
    communityBarplotParams()
  })
  
  # Taxa prevalence plot
  communityPrevalenceParams <- reactive({
    prevplot <- plot_taxa_prevalence(compositionalInput(), input$v4) + theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90) #If OTUs > 25 it fails
    if(input$transparentCommunityPlot == TRUE){
      prevplot <- prevplot + 
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    ggplotly(prevplot, height = 500, width = 1000)
  })
  
  output$communityPrevalence <- renderPlotly({
    communityPrevalenceParams()
  })
  
  
  # Phyloseq Summary #
  summaryParams <- reactive({
    req(datasetInput())
    summarize_phyloseq(datasetInput())
  })
  
  output$summary <- renderPrint({
    summaryParams()
  })
  
  ## Alpha Diversity ##
  
  #Abundance and Evenness tables#
  
  evennessParams <- reactive({
    datatable(evenness(datasetInput()),options = list(scrollX = TRUE))
  })
  output$evennessTable <- renderDataTable({
    evennessParams()
  })
  
  absoluteAbundanceParams <- reactive({
    datatable(abundances(datasetInput()), options = list(scrollX = TRUE))
  })
  output$absoluteAbundanceTable <- renderDataTable({
    absoluteAbundanceParams()
  })
  
  relativeAbundanceParams <- reactive({
    datatable(abundances(datasetInput(), transform = "compositional"), options = list(scrollX = TRUE))
  })
  output$relativeAbundanceTable <- renderDataTable({
    relativeAbundanceParams()
  })
  
  # Updating SelectInputs when database changes #
  observeEvent(input$datasetUpdate, {

    
    updateSelectInput(session, "x",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "x2", choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "x3", choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "y",
                      choices = colnames(alpha(datasetInput())))
  })
  
  
  # Merged table - generate and output #
  mergedTable <- reactive({
    merge(meta(datasetInput()), alpha(datasetInput()), all.y = TRUE)
  })
  viewParams <- reactive({
    datatable(mergedTable(),options = list(scrollX = TRUE))
  })

  output$view <- DT::renderDataTable({
    viewParams()
  })
  
  # Alpha Diversity Richness Plot #
  richnessPlotParams <- reactive({
    if(input$richnessPlotGridWrap == FALSE){
      richnessplot <- plot_richness(
        datasetInput(),
        x = input$x2,
        measures = input$richnessChoices,
        color = input$x3
      ) + theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90) + ylab(paste("Alpha Diversity Measure (", input$richnessChoices , ")"))
    } else {
      richnessplot <- plot_richness(
        datasetInput(),
        x = input$x2,
        measures = input$richnessChoices,
        color = input$x3
      ) + facet_grid(paste('~',input$x),scales = "free", space = "free") + theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90) + ylab(paste("Alpha Diversity Measure (", input$richnessChoices , ")"))
    }
    if(input$transparentRichness == TRUE){
      richnessplot <- richnessplot +  
      theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    ggplotly(richnessplot, height = 500, width = 1050)
  })
  
  output$richnessPlot <- renderPlotly({
    richnessPlotParams()
  })
  
  ## Beta Diversity ##
  
  # Updating SelectInputs when dataset changes#
  observeEvent(input$datasetUpdate, {

    
    updateSelectInput(session, "xb",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "xb2",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "xb3",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "yb",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "zb",
                      choices = colnames(tax_table(datasetInput())))
    updateSelectInput(session, "zbsplit",
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
    )
  })
  
  ordinateDataSplit <- reactive({
    ordinate(
      compositionalInput(),
      method = input$ordinate.method2,
      distance = input$ordinate.distance2
    )
  })
  
  ordinateDataTaxa <- reactive({
    ordinate(
      compositionalInput(),
      method = input$ordinate.method3,
      distance = input$ordinate.distance3
    )
  })
  
  ordinatePlotParams <- reactive({
    p <- phyloseq::plot_ordination(datasetInput(), ordinateData(), color = input$xb, label = input$yb ) + geom_point(size = input$geom.size) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
    if(input$transparentOrdinatePlot){
      p <- p +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    ggplotly(p, height = 500, width = 1050)
  })
  
  output$ordinatePlot <- renderPlotly({
    ordinatePlotParams()
  })
  
  splitOrdParams <- reactive({
    splitOrdplot <-
      plot_ordination(
        datasetInput(),
        ordinateDataSplit(),
        type = "split",
        shape = input$xb,
        #color = input$yb,
        color = input$zbsplit
      ) + geom_point(size = input$geom.size2) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
    if(input$transparentSplitOrd){
      splitOrdplot <- splitOrdplot +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    ggplotly(splitOrdplot, height = 500, width = 1050)
  })
  
  output$splitOrd <- renderPlotly({
    splitOrdParams()
  })
  
  taxaOrdParams <- reactive({
    taxaOrdplot <-
      plot_ordination(
        datasetInput(),
        ordinateDataTaxa(),
        type = "taxa",
        color = input$zb,
        label = input$xb
      ) + geom_point(size = input$geom.size3) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
    if(input$transparentTaxaOrd){
      taxaOrdplot <- taxaOrdplot +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    ggplotly(taxaOrdplot, height = 500, width = 1050)
  })
  
  output$taxaOrd <- renderPlotly({
    taxaOrdParams()
  })
  
  
  # # LANDSCAPE ANALYSIS # - REMOVED FOR NOW
  # # #Update metadata selectInputs when swapping datasets - Landscape analysis
  # observeEvent(input$datasetUpdate, {
  # 
  # 
  #   updateSelectInput(session, "metadataLandscape1",
  #                     choices=colnames(meta(datasetInput())))
  # })
  # observeEvent(input$datasetUpdate, {
  # 
  # 
  #   updateSelectInput(session, "metadataLandscape2",
  #                     choices=colnames(meta(datasetInput())))
  # })
  # 
  # #Update metadata values for landscape analysis
  # observeEvent(input$metadataLandscape1,{
  #   updateSelectInput(session, "metadataLandscapeValue1",
  #                     choices=(sample_data(datasetInput())[[input$metadataLandscape1]]))
  # })
  # observeEvent(input$metadataLandscape2,{
  #   updateSelectInput(session, "metadataLandscapeValue2",
  #                     choices=(sample_data(datasetInput())[[input$metadataLandscape2]]))
  # })
  # observeEvent(input$datasetUpdate, {
  # 
  # 
  #   updateSelectInput(session, "landscapeOTU",
  #                     choices=(rownames(otu_table(datasetInput()))))
  # })
  # 
  # # Data subsetting
  # datasetSubsetLandscape <- reactive({
  #   subset1 <- prune_samples(sample_data(datasetInput())[[input$metadataLandscape1]] == input$metadataLandscapeValue1, datasetInput())
  #   subset <- prune_samples(sample_data(subset1)[[input$metadataLandscape2]] == input$metadataLandscapeValue2, subset1)
  #   microbiome::transform(subset, "compositional")
  # })
  # 
  # #PCA plot
  # landscapePCAParams <- reactive({
  #   p <- plot_landscape(datasetInput(), method = "PCA", transformation = "clr", col = "sample") +
  #     labs(title = paste("PCA / CLR")) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
  #   ggplotly(p, height = 500, width = 1050)
  # })
  # 
  # output$landscapePCA <- renderPlotly({
  #   landscapePCAParams()
  # })
  # 
  # #PCoA/MDS
  # landscapePCoAMDSParams <- reactive({
  #   p <- plot_landscape(datasetSubsetLandscape, method = "PCoA", distance = "bray") +
  #     labs(title = paste("PCoA / Compositional / Bray-Curtis")) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
  #   ggplotly(p, height = 500, width = 1050)
  # })
  # 
  # output$landscapePCoAMDS <- renderPlotly({
  #   landscapePCoAMDSParams()
  # })
  # 
  # #NMDS
  # landscapeOrdinationParams <- reactive({
  #   x <- datasetInput()
  #   quiet(x.ord <- ordinate(x, input$ordinateMethodLandscape, input$ordinateDistanceLandscape))
  #   proj <- phyloseq::plot_ordination(x, x.ord, justDF=TRUE)
  #   names(proj)[1:2] <- paste("Comp", 1:2, sep=".")
  #   p <- plot_landscape(proj[, 1:2], col = proj[[input$metadataLandscape1]], add.points = TRUE, legend = TRUE) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
  #   ggplotly(p, height = 500, width = 1050)
  # })
  # 
  # output$landscapeOrdination <- renderPlotly({
  #   landscapeOrdinationParams()
  # })
  # 
  # 
  # #t-SNE
  # landscapeTSneParams <- reactive({
  #   p <- plot_landscape(datasetInput(), "t-SNE",
  #                       distance = "euclidean", transformation = "hellinger") +
  #     labs(title = paste("t-SNE / Hellinger / Euclidean")) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
  #   ggplotly(p, height = 500, width = 1050)
  # })
  # output$landscapeTSne <- renderPlotly({
  #   landscapeTSneParams()
  # })
  # 
  # #Abundance histogram (Absolute/Relative)
  # landscapeAbundanceHistAbsParams <- reactive({
  #   # Visualize population densities for specific taxa
  #   p <- plot_density(datasetInput(), input$landscapeOTU) + ggtitle("Abundance density")
  #   ggplotly(p, height = 500, width = 1050)
  # })
  # output$landscapeAbundanceHistAbs <- renderPlotly({
  #   landscapeAbundanceHistAbsParams()
  # })
  # 
  # landscapeAbundanceHistRelParams <- reactive({
  #   # Visualize population densities for specific taxa
  #   p <- plot_density(compositionalInput(), input$landscapeOTU, log10 = TRUE, fill = "gray") +
  #     ggtitle("Relative abundance") +
  #     xlab("Relative abundance (%)")
  #   ggplotly(p, height = 500, width = 1050)
  # })
  # output$landscapeAbundanceHistRel <- renderPlotly({
  #   landscapeAbundanceHistRelParams()
  # })
  # 
  # #Output code - Landscape analysis
  # 
  # landscapePCACodeParams <- reactive({
  #   paste0("p <- plot_landscape(",input$datasetSample,", method = 'PCA', transformation = 'clr') + labs(title = paste('PCA / CLR')) plotly_build(p)")
  # })
  # 
  # landscapeOrdinationCodeParams <- reactive({
  #   paste0("x <- ",input$datasetSample,"quiet(x.ord <- ordinate(x, ",input$ordinateMethodLandscape,", ",input$ordinateDistanceLandscape,")) proj <- phyloseq::plot_ordination(x, x.ord, justDF=TRUE) names(proj)[1:2] <- paste('Comp0, 1:2, sep='.') p <- plot_landscape(proj[, 1:2], col = proj[[",input$metadataLandscape1,"]], add.points = TRUE, legend = TRUE) plotly_build(p)")
  # })
  # 
  # landscapeTSneCodeParams <- reactive({
  #   paste0("p <- plot_landscape(",input$datasetSample,", 't-SNE', distance = 'euclidean', transformation = 'hellinger') + labs(title = paste('t-SNE / Hellinger / Euclidean')) plotly_build(p)")
  # })
  # 
  # landscapeAbundanceHistAbsCodeParams <- reactive({
  #   paste0("p <- plot_density(",input$datasetSample,", 'Dialister') + ggtitle('Absolute abundance') plotly_build(p)")
  # })
  # landscapeAbundanceHistRelCodeParams <- reactive({
  #   paste0("p <- plot_density(microbiome::transform(",input$datasetSample,", 'compositional'), 'Dialister', log10 = TRUE) + ggtitle('Relative abundance') + xlab('Relative abundance (%)') plotly_build(p)")
  # })
  # 
  # #DMM COMMUNITY TYPING#
  # 
  # # Pick OTU count matrix, convert into samplex x taxa format and fit the DMM model
  # dmmModelFit <- reactive({
  #   taxa <- core_members(compositionalInput(), detection = input$detectionDMM/100, prevalence = input$prevalenceDMM/100) #need to make this related to core comp, along with the coremembers in beta diversity
  #   pseq <- prune_taxa(taxa, datasetInput())
  #   dat <- abundances(pseq)
  #   count <- as.matrix(t(dat))
  #   mclapply(1:input$maxModelsDMM, dmn, count = count, verbose=TRUE)
  # })
  # 
  # dmmModelCheckParams <- reactive({
  #   lplc <- sapply(dmmModelFit(), laplace) # AIC / BIC / Laplace
  #   aic  <- sapply(dmmModelFit(), AIC) # AIC / BIC / Laplace
  #   bic  <- sapply(dmmModelFit(), BIC) # AIC / BIC / Laplace
  #   plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
  #   lines(aic, type="b", lty = 2)
  #   lines(bic, type="b", lty = 3)
  # })
  # 
  # output$dmmModelCheck <- renderPlot({
  #   dmmModelCheckParams()
  # })
  # 
  # dmmModelBestFit <- reactive({
  #   lplc <- sapply(dmmModelFit(), laplace) # AIC / BIC / Laplace
  #   aic  <- sapply(dmmModelFit(), AIC) # AIC / BIC / Laplace
  #   bic  <- sapply(dmmModelFit(), BIC) # AIC / BIC / Laplace
  #   dmmModelFit()[[which.min(lplc)]]
  # })
  # 
  # output$dmmParameters <- renderDT({
  #   mixturewt(dmmModelBestFit())
  # })
  # output$sampleAssignments <- renderDT({
  #   apply(mixture(dmmModelBestFit()), 1, which.max)
  # })
  # 
  # taxaContributionPerComponentParams <- reactive({
  #   for (k in seq(ncol(fitted(dmmModelBestFit())))) {
  #     d <- melt(fitted(dmmModelBestFit()))
  #     colnames(d) <- c("OTU", "cluster", "value")
  #     d <- subset(d, cluster == k) %>%
  #       # Arrange OTUs by assignment strength
  #       arrange(value) %>%
  #       mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  #       # Only show the most important drivers
  #       filter(abs(value) > quantile(abs(value), 0.8))
  # 
  #     p <- ggplot(d, aes(x = OTU, y = value)) +
  #       geom_bar(stat = "identity") +
  #       coord_flip() +
  #       labs(title = paste("Top drivers: community type", k))
  #   }
  #   ggplotly(p)
  # })
  # 
  # output$taxaContributionPerComponent <- renderPlotly({
  #   taxaContributionPerComponentParams()
  # })
  # 
  # #Output code - DMM community typing
  # dmmModelFitCode1 <- reactive({
  #   paste0("taxa <- core_members(microbiome::transform(",input$datasetSample,",'compositional'), detection = ",input$detectionDMM,"/100, prevalence = ",input$prevalenceDMM,"/100)")
  # })
  # dmmModelFitCode2 <- reactive({
  #   paste0("pseq <- prune_taxa(taxa, ",input$datasetSample,")")
  # })
  # dmmModelFitCode3 <- reactive({
  #   paste0("dat <- abundances(pseq)")
  # })
  # dmmModelFitCode4 <- reactive({
  #   paste0("count <- as.matrix(t(dat))")
  # })
  # dmmModelFitCode5 <- reactive({
  #   paste0("fit <- mclapply(1:",input$maxModelsDMM,", dmn, count = count, verbose=TRUE)")
  # })
  # 
  # dmmModelCheckCodeParams1 <- reactive({
  #   paste0("lplc <- sapply(fit, laplace)") # AIC / BIC / Laplace
  # })
  # 
  # dmmModelCheckCodeParams2 <- reactive({
  #   paste0("aic  <- sapply(fit, AIC)") # AIC / BIC / Laplace
  # })
  # dmmModelCheckCodeParams3 <- reactive({
  #   paste0("bic  <- sapply(fit, BIC)") # AIC / BIC / Laplace
  # })
  # dmmModelCheckCodeParams4 <- reactive({
  #   paste0("plot(lplc, type='b', xlab='Number of Dirichlet Components', ylab='Model Fit')")
  # })
  # dmmModelCheckCodeParams5 <- reactive({
  #   paste0("lines(aic, type='b', lty = 2)")
  # })
  # dmmModelCheckCodeParams6 <- reactive({
  #   paste0("lines(bic, type='b', lty = 3)")
  # })
  # 
  # dmmModelBestFitCode1 <- reactive({
  #   paste0("lplc <- sapply(fit, laplace)") # AIC / BIC / Laplace
  # })
  # 
  # dmmModelBestFitCode2 <- reactive({
  #   paste0("aic  <- sapply(fit, AIC)") # AIC / BIC / Laplace
  # })
  # dmmModelBestFitCode3 <- reactive({
  #   paste0("bic  <- sapply(fit, BIC)") # AIC / BIC / Laplace
  # 
  # })
  # dmmModelBestFitCode4 <- reactive({
  #   paste0("bestFit <- fit[[which.min(lplc)]]")
  # 
  # })
  # 
  # dmmParametersCode <- reactive({
  #   paste0("datatable(mixturewt(bestFit))")
  # })
  # sampleAssignmentsCode <- reactive({
  #   paste0("datatable(apply(mixture(bestFit), 1, which.max))")
  # })
  # 
  # taxaContributionPerComponentCodeParams1 <- reactive({
  #   paste0( "for (k in seq(ncol(fitted(bestFit)))) {")
  # })
  # 
  # taxaContributionPerComponentCodeParams2 <- reactive({
  #   paste0( "d <- melt(fitted(bestFit))")
  # })
  # taxaContributionPerComponentCodeParams3 <- reactive({
  #   paste0( "colnames(d) <- c('OTU', 'cluster', 'value')")
  # })
  # taxaContributionPerComponentCodeParams4 <- reactive({
  #   paste0("d <- subset(d, cluster == k) %>% arrange(value) %>% mutate(OTU = factor(OTU, levels = unique(OTU))) %>% filter(abs(value) > quantile(abs(value), 0.8))")
  # })
  # taxaContributionPerComponentCodeParams5 <- reactive({
  #   paste0("p <- ggplot(d, aes(x = OTU, y = value)) + geom_bar(stat = 'identity') + coord_flip() + labs(title = paste('Top drivers: community type', k))")
  # })
  # taxaContributionPerComponentCodeParams6 <- reactive({
  #   paste0("print(p)}")
  # })

  ###########################
  ## Statistical analysis ###
  ###########################

  ## PERMANOVA ##
  #Update metadata column when dataset changes
  observeEvent(input$datasetUpdate, {
    updateSelectInput(session, "permanovaColumn",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "permanovaColumnP",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "permanovaColumnFac",
                      choices = colnames(meta(datasetInput())))
  })

  permanova <- reactive({
    otu <- abundances(compositionalInput())
    meta <- meta(compositionalInput())
    permnumber <- input$permanovaPermutationsP
    metadata <- input$permanovaColumnP
    adonis(t(otu) ~ meta[[metadata]],
           data = meta, permutations = permnumber, method = input$permanovaDistanceMethodP, parallel = getOption("mc.cores")
    )
  })
  
  output$pValue <- renderDataTable({
    as.data.frame(permanova()$aov.tab)
  })

  homogenietyParams <- reactive({
    otu <- abundances(compositionalInput())
    meta <- meta(compositionalInput())
    dist <- vegdist(t(otu))
    metadata <- input$permanovaColumnP
    anova(betadisper(dist, meta[[metadata]]))
  })

  output$homogeniety <- renderDataTable({
    homogenietyParams()
  })

  topFactorPlotParams <- reactive({
    otu <- abundances(compositionalInput())
    meta <- meta(compositionalInput())
    permnumber <- input$permanovaPermutationsFac
    metadata <- input$permanovaColumnFac
    column <- meta[[metadata]]
    permanova <- adonis(t(otu) ~ column,
                        data = meta, permutations = permnumber, method = "bray"
    )
    coef <- coefficients(permanova)["column1",]
    top.coef <- coef[rev(order(abs(coef)))[1:20]] #top 20 coefficients
    par(mar = c(3, 14, 2, 1))
     p <- barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
    if(input$transparentPermanova == TRUE){
      p <- p +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    print(p)
  })

  output$topFactorPlot <- renderPlot({
    topFactorPlotParams()
  })

  netPlotParams <- reactive({
    n <- make_network(compositionalInput(), type = input$permanovaPlotTypeNet, distance = input$permanovaDistanceMethodNet)
    #p <- plot_network(n, compositionalInput(), type = "taxa", color= ntaxa(otu_table(compositionalInput())))
    p <- plot_network(n, compositionalInput(), shape = "nationality", color= "nationality")
    
        if(input$transparentPermanova == TRUE){
      p <- p +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    ggplotly(p, height = 500, width = 1050)
  })

  output$netPlot <- renderPlotly({
    netPlotParams()
  })

  output$permaHeatmap <- renderPlotly({
    xp <- plot_heatmap(compositionalInput(), distance = ordinate(compositionalInput(), distance = input$permanovaDistanceMethodHeat), method = input$permanovaMethodHeat)
    ggplotly(p, height = 500, width = 1050)
  })
  
  # output$permaHeatmap <- renderPlotly({
  #   permaHeatmapParams()
  # })
  
  # ANOSIM #
  #Update metadata column#
  # observe({
  #   updateSelectInput(session, "anosimColumn",
  #                     choices = colnames(meta(datasetInput())))
  # })
  # 
  # anosim <- reactive({
  #   otu <- abundances(compositionalInput())
  #   meta <- meta(compositionalInput())
  #   permnumber <- input$anosimPermutations
  #   metadata <- input$anosimColumn
  # })
  # 
  # output$pValueAnosim <- renderPrint({
  #   otu <- abundances(compositionalInput())
  #   meta <- meta(compositionalInput())
  #   metadata <- input$anosimColumn
  #   vegan::anosim(otu, meta[[metadata]], permutations = 99, distance = "bray", parallel = getOption("mc.cores"))
  # })
  
  
  
  ### RESULTS ###
  
  # Core Micrboiota  
  output$downloadReportCoreMicro <- downloadHandler(
    filename = function() {
      paste('core-microbiota-report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html'
      ))
    },
    content = function(file) {
      src <- normalizePath('core_microbiota_report.Rmd')
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'core_microbiota_report.Rmd', overwrite = TRUE)
      
      out <- rmarkdown::render('core_microbiota_report.Rmd',
                               switch(input$format,
                                      PDF = pdf_document(), 
                                      HTML = html_document() 
                               ))
      file.rename(out, file)
    }
  )
  
  # Community Composition #  
  
  output$downloadReportComposition <- downloadHandler(
    filename = function() {
      paste('community-composition-report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html'
      ))
    },
    content = function(file) {
      src <- normalizePath('community_comp_microbiota_report.Rmd')
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'community_comp_microbiota_report.Rmd', overwrite = TRUE)
      
      out <- rmarkdown::render('community_comp_microbiota_report.Rmd',
                               switch(input$format,
                                      PDF = pdf_document(), 
                                      HTML = html_document() 
                               ))
      file.rename(out, file)
    }
  )  
  
  
  output$downloadReportAlpha <- downloadHandler(
    filename = function() {
      paste('report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html'
      ))
    },
    content = function(file) {
      src <- normalizePath('alpha_diversity_report.Rmd')
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'alpha_diversity_report.Rmd', overwrite = TRUE)
      
      out <- rmarkdown::render('alpha_diversity_report.Rmd',
                               switch(input$format,
                                      PDF = pdf_document(), 
                                      HTML = html_document() 
                               ))
      file.rename(out, file)
    }
  )
  
  # Beta Diversity
  
  output$downloadReportBeta <- downloadHandler(
    filename = function() {
      paste('beta-report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html'
      ))
    },
    content = function(file) {
      src <- normalizePath('beta_diversity_report.Rmd')
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'beta_diversity_report.Rmd', overwrite = TRUE)
      
      out <- rmarkdown::render('beta_diversity_report.Rmd',
                               switch(input$format,
                                      PDF = pdf_document(), 
                                      HTML = html_document() 
                               ))
      file.rename(out, file)
    }
  )
  
  # #Community Landscape#
  # output$downloadReportLandscape <- downloadHandler(
  #   filename = function() {
  #     paste('landscape-report', sep = '.', switch(
  #       input$format, PDF = 'pdf', HTML = 'html'
  #     ))
  #   },
  #   content = function(file) {
  #     src <- normalizePath('landscape_report.Rmd')
  #     
  #     # temporarily switch to the temp dir, in case you do not have write
  #     # permission to the current working directory
  #     owd <- setwd(tempdir())
  #     on.exit(setwd(owd))
  #     file.copy(src, 'landscape_report.Rmd', overwrite = TRUE)
  #     
  #     out <- rmarkdown::render('landscape_report.Rmd',
  #                              switch(input$format,
  #                                     PDF = pdf_document(), 
  #                                     HTML = html_document() 
  #                              ))
  #     file.rename(out, file)
  #   }
  # )
  # 
  # #DMM Community Typing#
  # output$downloadReportDMM <- downloadHandler(
  #   filename = function() {
  #     paste('dmm-report', sep = '.', switch(
  #       input$format, PDF = 'pdf', HTML = 'html'
  #     ))
  #   },
  #   content = function(file) {
  #     src <- normalizePath('dmm_report.Rmd')
  #     
  #     # temporarily switch to the temp dir, in case you do not have write
  #     # permission to the current working directory
  #     owd <- setwd(tempdir())
  #     on.exit(setwd(owd))
  #     file.copy(src, 'dmm_report.Rmd', overwrite = TRUE)
  #     
  #     out <- rmarkdown::render('dmm_report.Rmd',
  #                              switch(input$format,
  #                                     PDF = pdf_document(), 
  #                                     HTML = html_document() 
  #                              ))
  #     file.rename(out, file)
  #   }
  # )
  # 
  #PERMANOVA#
  output$downloadReportPermanova <- downloadHandler(
    filename = function() {
      paste('permanova-report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html'
      ))
    },
    content = function(file) {
      src <- normalizePath('permanova_report.Rmd')
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'permanova_report.Rmd', overwrite = TRUE)
      
      out <- rmarkdown::render('permanova_report.Rmd',
                               switch(input$format,
                                      PDF = pdf_document(), 
                                      HTML = html_document() 
                               ))
      file.rename(out, file)
    }
  )
  
}

# Run the application
shinyApp(ui = ui, server = server, options = list(port = 5000))
