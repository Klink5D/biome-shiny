# Biome-shiny 0.6

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
#library(limma)


# Load sample datasets #
data("dietswap")
data("atlas1006")
data("peerj32")
peerj32 <- peerj32$phyloseq

# UI
ui <- dashboardPage(
  dashboardHeader(title = "biome-shiny v0.6"),
  dashboardSidebar(
    sidebarMenu(
      menuItem(
        "Introduction",
        tabName = "intro",
        icon = icon("dashboard")
      ),
      menuItem("Phyloseq Summary", tabName="phyloseqsummary"),
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
      h1("Introduction to biome-shiny"),
      br(),
      box(
        width = "3",
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
            accept = c("text/biom"), placeholder="Phyloseq .biom files"
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
      )
    ),
    
    tabItem(tabName = "phyloseqsummary",
            h1("Phyloseq Summary"),
            verbatimTextOutput("summary")
    ),
    
    
    # Core microbiota #
    tabItem(
      tabName = "coremicrobiota",
      # box(
      #   title = "Variables",
      #   width = "2",
      #   collapsible = TRUE,
      #   collapsed = TRUE,
      #   
        #numericInput("detectionPrevalence", "For prevalences: Choose detection value", min = 0.00, max = 100, value = 0.01, step = 0.01),
        # numericInput("prevalencePrevalence","Input a prevalence value", min = 0, max = 1, value = 0.5, step = 0.05),
        # selectInput("prevalenceSelection", multiple = TRUE, selectize = TRUE, choices = c(seq(0,1,by=0.01)), label = "Choose prevalence values for lineplot and heatmap"),
        # textInput("detectionForLineplot", label = "Enter detection values, separated by comma, for the lineplot"),
        # textInput("detectionMin", label = "Enter minimum limit for detection"),
        # textInput("detectionMax", label = "Enter maximum limit for detection"),
        # textInput("maxLength", label = "Enter length (number of values)")
      # ),
      
      tabsetPanel(
        tabPanel("Prevalence (absolute/relative)",
                 # Output the prevalence in relatives and absolutes (Counts)
          tabsetPanel(
            tabPanel("Variables",
              box( width ="2", collapsible = TRUE,
                numericInput("detectionPrevalence", "For prevalences: Choose detection value", min = 0.00, max = 100, value = 0.01, step = 0.01)              
              )
            ),
            tabPanel("Absolute prevalence", dataTableOutput("prevalenceAbsoluteOutput")),
            tabPanel("Relative prevalence", dataTableOutput("prevalenceRelativeOutput"))
          )
        ),
        tabPanel("Core Taxa Summary",
          tabsetPanel(
            tabPanel("Variables",
                     numericInput("detectionPrevalence2", "For prevalences: Choose detection value", min = 0.00, max = 100, value = 0.01, step = 0.01),
                     numericInput("prevalencePrevalence","Input a prevalence value", min = 0, max = 1, value = 0.5, step = 0.05)
            ),
            tabPanel("Summary", verbatimTextOutput("corePhyloSummary")),
            tabPanel("Taxa", verbatimTextOutput("coreTaxa"))
          )
        ),
        tabPanel("Core Taxa Visualization",
            tabsetPanel(
              tabPanel("Variables (Lineplot)",
                 box(width = "2", collapsible = TRUE,
                   selectInput("prevalenceSelection", multiple = TRUE, selectize = TRUE, choices = c(seq(0,1,by=0.01)), label = "Choose prevalence values for lineplot"),
                   textInput("detectionForLineplot", label = "Enter detection values, separated by comma, for the lineplot")
                 )
              ),
              tabPanel("Variables (Heatmap)",
                 box(width = "2", collapsible = TRUE,
                     selectInput("prevalenceSelectionHeat", multiple = TRUE, selectize = TRUE, choices = c(seq(0,1,by=0.01)), label = "Choose prevalence values for heatmap"),
                     textInput("detectionMin", label = "Enter minimum limit for detection"),
                     textInput("detectionMax", label = "Enter maximum limit for detection"),
                     numericInput("maxLength", label = "Enter length (number of values)", min = "1", value= "10")
                )
              ),
              tabPanel("Lineplot", plotlyOutput("coreLineplot")),
              tabPanel("Heatmap", plotlyOutput("coreHeatmap"))
            )
        )
      )
    ),
    
    
    # Community composition #
    tabItem(
      tabName = "communitycomposition",
      tabsetPanel(
        tabPanel(title = "Abundance in samples by taxa",
          tabsetPanel(       
                 tabPanel("Variables",
                          box(
                            title = "Data Subsetting",
                            width = "2",
                            collapsible = TRUE,
                            selectInput(
                              "z1",
                              "Choose a metadata column:",
                              # For subsetting data, #1
                              choices = colnames("datasetMetadata"),
                              selected = "bmi_group"
                            ),
                            selectInput(
                              "z2",
                              "Choose a seocond metadata column:",
                              # For subsetting data, #2
                              choices = colnames("datasetMetadata"),
                              selected = "nationality"
                            ),
                            selectInput(
                              "z3",
                              "Choose the intended timepoint (if applicable):",
                              # For subsetting data, timepoint data
                              choices = colnames("datasetMetadata"),
                              selected = "timepoint.within.group"
                            )
                          ),
                          box("Metadata Values", width = "2", collapsible = TRUE,
                              selectInput(
                                "v1",
                                "Choose a metadata value:",
                                # For subsetting data, metadata value
                                choices = sapply("datasetMetadata", levels),
                                selected = "lean"
                              ),
                              selectInput(
                                "v2",
                                "Choose a second metadata value:",
                                # For subsetting data, metadata #2 value
                                choices = sapply("datasetMetadata", levels),
                                selected = "AAM"
                              ),
                              selectInput(
                                "v3",
                                "Choose a timepoint value (if applicable):",
                                # For subsetting data, timepoint data value
                                choices = c("","0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
                                selected = "1"
                              ),
                              selectInput(
                                "v4",
                                "Choose a taxonomy rank:",
                                # Tax rank to analyze
                                choices = c("Phylum", "Class", "Order", "Family", "Genus"),
                                selected = "Phylum"
                              )
                          )
                 ),
                 tabPanel("Plot", plotOutput("communityPlot"))
                 )
          ),
        tabPanel(title = "Relative abundance",
                 plotOutput("communityPlotGenus")),
        tabPanel(title = "Relative abundance, averaged by metadata",
          tabsetPanel(
            tabPanel("Variables",
                  box(width = "2", collapsible = TRUE,
                     selectInput(
                       "z1Average",
                       "Choose a metadata column for averaging:",
                       # For subsetting data, #1
                       choices = colnames("datasetMetadata"),
                       selected = "bmi_group"
                     ),
                     selectInput(
                       "plotTypeZ1Average",
                       "Plot type:",
                       choices = c("barplot","heatmap"),
                       selected = "barplot"
                     )
                  )
            ),
             tabPanel("Plot", plotOutput("communityBarplot"))
          )
        ),
        tabPanel(title = "Taxonomy Prevalence Plot",
           tabsetPanel(
              tabPanel("Variables",
                box(title = "Taxonomy rank", collapsible = TRUE, width = "2",
                  selectInput("v4Plot", label = "Choose a taxonomy rank", choices = c("Phylum", "Class", "Order", "Family", "Genus"), selected = "Phylum")
                )
              ),
              tabPanel("Plot", plotlyOutput("communityPrevalence"))
          )
        )
      )
    ),
    
    # Alpha Diversity #
    tabItem(
      tabName = "alphadiversity",
      # box(
      #   width = "2",
      #   collapsible = TRUE,
      #   title = "Variables",
      #   collapsed = TRUE,
      #   # #X (the metadata) and Y (the diversity measure)
      #   selectInput(
      #     "x",
      #     "Choose a metadata column test:",
      #     choices = colnames("datasetMetadata")
      #   ),
      #   selectInput("y", "Choose a diversity measure:",
      #               choices = colnames("datasetMetadata"))
      # ),
      br(),
      tabsetPanel(
        type = "tabs",
        id = "tabsetpanel",
      
        tabPanel(title = "Evenness Table",
                 dataTableOutput("evennessTable")),
        tabPanel(title = "Abundance Table",
                 dataTableOutput("abundanceTable")),        
        
        # Alpha diversity measures with metadata
        tabPanel(title = "Metadata Table with diversity measures",
                 DT::dataTableOutput("view")),
        
        # # Violin Plots
        # tabPanel(
        #       title = "Violin Plot",
        #     tabsetPanel(
        #       tabPanel(title = "Variables",
        #       box(
        #         width = "2",
        #         collapsible = TRUE,
        #         title = "Variables",
        #         collapsed = TRUE,
        #         # #X (the metadata) and Y (the diversity measure)
        #         selectInput(
        #           "x",
        #           "Choose a metadata column test:",
        #           choices = colnames("datasetMetadata")
        #         ),
        #         selectInput("y", "Choose a diversity measure:",
        #                     choices = colnames("datasetMetadata"))
        #       )),
        #       tabPanel(title = "Plot", plotlyOutput("violinPlot")))),
        
        # A phyloseq richness plot
        tabPanel(title = "Richness Plot",
                 tabsetPanel(
                  tabPanel( title = "Variables",
                     box(
                       width = "2",
                       collapsible = TRUE,
                       title = "Variables",
                       collapsed = TRUE,
                       # #X (the metadata) and Y (the diversity measure)
                       selectInput(
                         "x2",
                         "Choose a metadata column test:",
                         choices = colnames("datasetMetadata")
                       )
                     ),
                     box(checkboxGroupInput("richnessChoices", "Choose diversity measures" ,choices = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), selected = c("Shannon", "Simpson")))),
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
                              collapsed = TRUE,
                              selectInput(
                                "xb", "Choose a metadata column:", choices = colnames("datasetMetadata"), selected = "bmi_group"
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
                              sliderInput(
                                "geom.size",
                                "Plot geometry point size:",
                                min = 1,
                                max = 10,
                                step = 0.5,
                                value = "3"
                              )                            
                   )),
                   tabPanel(title = "Plot",
                            plotlyOutput("ordinatePlot"),
                            textOutput("ordinatePrint"))
                   )),
        
        tabPanel(title = "Split Ordination Plot (Metadata/Metadata)",
                 tabsetPanel(
                   tabPanel(title = "Variables",
                            box(
                              title = "Variables",
                              width = "2",
                              collapsible = TRUE,
                              collapsed = TRUE,
                              selectInput(
                                "xb2", "Choose a metadata column:", choices = colnames("datasetMetadata"), selected = "bmi_group"
                              ),
                              
                              selectInput(
                                "yb",
                                "For metadata/metadata split plots, choose a metadata column:",
                                choices = colnames("datasetMetadata"),
                                selected = "bmi_group"
                              ),
                              
                              selectInput(
                                "ordinate.method2",
                                "Choose an ordination method:",
                                choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                selected = "CCA"
                              ),
                              
                              selectInput(
                                "ordinate.distance2",
                                "Choose a distance method:",
                                choices = c("bray", "jaccard", "unifrac"),
                                selected = "unifrac"
                              ),
                              #There's 44 of these in total
                              
                              sliderInput(
                                "geom.size2",
                                "Plot geometry point size:",
                                min = 1,
                                max = 10,
                                step = 0.5,
                                value = "3"
                              )          
                  )),
                 tabPanel( title = "Plot", plotlyOutput("splitOrd")))),
        
        tabPanel(title = "Taxa Plot",
                 tabsetPanel(
                   tabPanel( title = "Variables",
                             box(
                               title = "Variables",
                               width = "2",
                               collapsible = TRUE,
                               collapsed = TRUE,
                               selectInput(
                                 "xb3", "Choose a metadata column:", choices = colnames("datasetMetadata"), selected = "bmi_group"
                               ),
                               selectInput(
                                 "zb",
                                 "For tax rank plot, choose a taxonomy rank:",
                                 choices = c("Phylum", "Class", "Order", "Family", "Genus")
                               ),
                               selectInput(
                                 "ordinate.method3",
                                 "Choose an ordination method:",
                                 choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                 selected = "CCA"
                               ),
                               selectInput(
                                 "ordinate.distance3",
                                 "Choose a distance method:",
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
                               )
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
         tabPanel("Abundance histogram (one-dimensional landscape)",
                  plotlyOutput("landscapeAbundanceHistAbs"),
                  plotlyOutput("landscapeAbundanceHistRel")
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
        tabPanel( title = "Population density plot",
           tabsetPanel(
              tabPanel(title = "Variables",
                box( title = "Variables", width= "2", collapsible = TRUE,
                     selectInput("permanovaDistanceMethod","Select distance method", choices = c("bray","jacard","unifrac"), selected = "unifrac"),
                     selectInput("permanovaMethod","Select ordination method",
                                 choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                 selected = "CCA"),
                     selectInput("permanovaColumn","Select metadata for density plot", choices = colnames("datasetMetadata")),
                     sliderInput("permanovaPlotSize", "Plot point size", min = 0.5, max = 10, step = 0.5, value = "3")
                )
              ),
              tabPanel(title = "Plot", plotlyOutput("densityPlot"))
           )
        ),
        tabPanel( title = "P-Value",
           tabsetPanel(
             tabPanel(title = "Variables",
                      box( title = "Variables", width= "2", collapsible = TRUE,
                           selectInput("permanovaDistanceMethodP","Select distance method", choices = c("bray","jacard","unifrac"), selected = "unifrac"),
                           selectInput("permanovaMethodP","Select ordination method",
                                       choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                       selected = "CCA"),
                           selectInput("permanovaColumnP","Select metadata for density plot", choices = colnames("datasetMetadata")),
                           numericInput("permanovaPermutationsP", "Number of permutations", min = 1, step = 1, value = 99)
                      )              
            ),
             tabPanel(title = "Data Tables",
                  dataTableOutput("pValue"),
                  dataTableOutput("homogeniety")
             )
           )
        ),
        tabPanel ( title = "Top Factors",
          tabsetPanel(
            tabPanel(title = "Variables",
                     box( title = "Variables", width= "2", collapsible = TRUE,
                          selectInput("permanovaDistanceMethodFac","Select distance method", choices = c("bray","jacard","unifrac"), selected = "unifrac"),
                          selectInput("permanovaColumnFac","Select metadata for density plot", choices = colnames("datasetMetadata")),
                          numericInput("permanovaPermutationsFac", "Number of permutations", min = 1, step = 1, value = 99)
                      )              
            ),
            tabPanel(title = "Plot",
                   plotOutput("topFactorPlot")
            )
          )
        ),
        tabPanel ( title = "Network Plot",
          tabsetPanel(
            tabPanel(title = "Variables",
                     box( title = "Variables", width= "2", collapsible = TRUE,
                          selectInput("permanovaDistanceMethodNet","Select distance method", choices = c("bray","jacard","unifrac"), selected = "unifrac"),
                          selectInput("permanovaMethodNet","Select ordination method",
                                      choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                      selected = "CCA"),
                          sliderInput("permanovaPlotSizeNet", "Plot point size", min = 0.5, max = 10, step = 0.5, value = "3"),
                          numericInput("permanovaPermutationsNet", "Number of permutations", min = 1, step = 1, value = 99)
                      )
            ),
            tabPanel(title ="Plot",
                     plotlyOutput("netPlot")
            )
          )
        ),
        tabPanel ( title = "Heatmap",
            tabsetPanel(
              tabPanel(title = "Variables",
                       box( title = "Variables", width= "2", collapsible = TRUE,
                            selectInput("permanovaDistanceMethodHeat","Select distance method", choices = c("bray","jacard","unifrac"), selected = "unifrac"),
                            selectInput("permanovaMethodHeat","Select ordination method",
                                        choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                        selected = "CCA"),
                            numericInput("permanovaPermutationsHeat", "Number of permutations", min = 1, step = 1, value = 99)
                        )              
                       
              ),
              tabPanel(title = "Plot",
                   plotlyOutput("permaHeatmap")
              ) 
          )
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
    tabItem(tabName = "results",
       tabsetPanel(
         tabPanel("Category Options",
           tabsetPanel(
             tabPanel("Core Microbiota",
                    box(
                      radioButtons("renderPrevTables", "Render Prevalence Tables", choices = c("Yes", "No")),
                      radioButtons("renderCoreTaxaSummary", "Render Core Taxa Summary", choices = c("Yes", "No")),
                      radioButtons("renderCoreLineplot", "Render Core Taxa Lineplot", choices = c("Yes", "No")),
                      radioButtons("renderCoreHeatmap", "Render Core Taxa Heatmap", choices = c("Yes", "No"))
                    ),
                    box(
                      paste0("Format, button to download as format"),
                      radioButtons('format', 'Document format', c('HTML'),
                                   inline = TRUE, selected = 'HTML'),
                      downloadButton('downloadReportCoreMicro')
                    )
             ),            
             tabPanel("Community Composition",
                    box(
                      radioButtons("renderSampleAbundancePlot", "Render Abundance in Samples Plot", choices = c("Yes", "No")),
                      radioButtons("renderRelativeAbundancePlot", "Render Relative Abundance Plot", choices = c("Yes", "No")),
                      radioButtons("renderAveragedRelativeAbundancePlot", "Render Averaged Relative Abundance Plot", choices = c("Yes", "No")),
                      radioButtons("renderPrevalencePlot", "Render Taxa Prevalence Plot", choices = c("Yes", "No"))
                    ),
                    box(title = "Download Options",
                        paste0("Format, button to download as format"),
                        radioButtons('format', 'Document format', c('HTML'),
                                     inline = TRUE, selected = 'HTML'),
                        downloadButton('downloadReportComposition')
                    )
             ),
             tabPanel("Alpha Diversity",
                box(
                  #radioButtons("renderViolin", "Render Violin Plot", choices = c("Yes", "No")),
                  radioButtons("renderRichness", "Render Richness Plot", choices = c("Yes", "No")),
                  radioButtons("printTableAlpha", "Print Table Head", choices = c("Yes", "No")),
                  radioButtons("printSummary", "Print Phyloseq Summary", choices = c("Yes", "No"))
                ),
                box(title = "Download Options",
                    paste0("Format, button to download as format"),
                    radioButtons('format', 'Document format', c('HTML'),
                                 inline = TRUE, selected = 'HTML'),
                    downloadButton('downloadReportAlpha')
                )
                
              ),
             tabPanel("Beta Diversity",
                box(
                  radioButtons("renderOrdplot", "Render Ordination Plot", choices = c("Yes", "No")), 
                  radioButtons("renderSplitOrdplot", "Render Split Ordination Plot", choices = c("Yes", "No")),
                  radioButtons("renderTaxaplot", "Render Taxa Plot", choices = c("Yes", "No"))
                ),
                box(
                  paste0("Format, button to download as format"),
                  radioButtons('format', 'Document format', c('HTML'),
                               inline = TRUE, selected = 'HTML'),
                  downloadButton('downloadReportBeta')
                )
              ),
             tabPanel("Community Landscape",
                box(
                  radioButtons("renderPCA", "Render PCA Plot", choices = c("Yes","No")),
                  radioButtons("renderPCOAMDS", "Render PCoA/MDS/NMDS Plot", choices = c("Yes","No")),
                  radioButtons("renderTSNE", "Render t-SNE Plot", choices = c("Yes","No")),
                  radioButtons("renderAbundanceHistograms", "Render Abundance Histograms", choices = c("Yes","No"))
                ),
                box(
                  paste0("Format, button to download as format"),
                  radioButtons('format', 'Document format', c('HTML'),
                               inline = TRUE, selected = 'HTML'),
                  downloadButton('downloadReportLandscape')
                )
              ),
             tabPanel("DMM Clustering",
                box(
                  radioButtons("renderModelVerification", "Render Model Verification Plot", choices = c("Yes", "No")),
                  radioButtons("renderParameters", "Render Alpha and Theta Parameters", choices = c("Yes", "No")),
                  radioButtons("renderTaxaContributionPlot", "Render Taxa Contribution Plots", choices = c("Yes", "No"))
                ),
                box(
                  paste0("Format, button to download as format"),
                  radioButtons('format', 'Document format', c('HTML'),
                               inline = TRUE, selected = 'HTML'),
                  downloadButton('downloadReportDMM')
                )
              ),
             tabPanel("PERMANOVA",
              box(
                  radioButtons("renderDensityPlot", "Render density plot", choices = c("Yes","No")),
                  radioButtons("renderPValueTables", "Render P-Value tables", choices = c("Yes","No")),
                  radioButtons("renderFactorPlot", "Render top factors barplot", choices = c("Yes","No")),
                  radioButtons("renderNetworkMap", "Render network map", choices = c("Yes","No"))
                ),
              box(
                paste0("Format, button to download as format"),
                radioButtons('format', 'Document format', c('HTML'),
                             inline = TRUE, selected = 'HTML'),
                downloadButton('downloadReportPermanova')
              )
            )
          ) 
         )
       )
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
  
  output$testThing <- renderPrint({
    if (input$datasetChoice == "Use sample dataset") {
      paste0(input$datasetSample)
    } else {
        deparse(substitute(biomfile))
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
    core(compositionalInput(), detection = input$detectionPrevalence2, prevalence = input$prevalencePrevalence )
  })
  output$corePhyloSummary <- renderPrint({ # Summary of corePhylo file
    summarize_phyloseq(corePhylo())  
  })
  output$coreTaxa <- renderPrint({ # Reports the taxa in corePhylo
    taxa(corePhylo())
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
    prevalences <- as.numeric(input$prevalenceSelectionHeat)
    detections <- 10^seq(log10(as.numeric(input$detectionMin)), log10(as.numeric(input$detectionMax)), length = as.numeric(input$maxLength))
    gray <- gray(seq(0,1,length=5))
    coreplot <- plot_core(compositionalInput(), plot.type = "heatmap", colours = gray, prevalences = prevalences, detections = detections, min.prevalence = .5) + xlab("Detection Threshold (Relative Abundance (%))")
    plotly_build(coreplot)
  })
  output$coreHeatmap <- renderPlotly({
    coreHeatmapParams()
  })
  
  #Output code - core microbiota
  prevalenceAbsoluteCode <- reactive({
    paste0("as.data.frame(prevalence(microbiome::transform(",input$datasetSample,", type = 'compositional'), detection = ",input$detectionPrevalence,"/100, sort = TRUE, count = TRUE))")
  })
  prevalenceRelativeCode <- reactive({
    paste0("as.data.frame(prevalence(microbiome::transform(",input$datasetSample,", type = 'compositional'), detection = ",input$detectionPrevalence,"/100, sort = TRUE))")
  })
  
  corePhyloCode <- reactive({
    paste0("corePhylo <- core(microbiome::transform(",input$datasetSample,", type = 'compositional'), detection = ",input$detectionPrevalence2,", prevalence = ",input$prevalencePrevalence,")")
  })
  
  corePhyloSummaryCode <- reactive({
    paste0("summarize_phyloseq(corePhylo)")
  })
  
  coreTaxaCode <- reactive({
    paste0("taxa(corePhylo)")
  })
  
  #This is inefficient, but... #There's a bug where a vector of multiple values writes the whole paste0 again, rather than just c(value1,value2,...,valueN)
  coreLineplotCodeParams1 <- reactive({
    paste0("prevalences <- as.numeric(",input$prevalenceSelection,")")
  })
  coreLineplotCodeParams2 <- reactive({
    paste0("detections <- as.numeric(unlist(strsplit(",input$detectionForLineplot,", split = ',')))/100")
  })
  coreLineplotCodeParams3 <- reactive({
    paste0("lineplot <- plot_core(microbiome::transform(",input$datasetSample,", type = 'compositional'), prevalences = prevalences, detections = detections, plot.type = 'lineplot') + xlab('Relative Abundance (%)')")
  })
  coreLineplotCodeParams4 <- reactive({
    paste0("plotly_build(lineplot)")
  })
  
  coreHeatmapCodeParams1 <- reactive({
    paste0("prevalences <- as.numeric(",input$prevalenceSelectionHeat,")")
  })
  coreHeatmapCodeParams2 <- reactive({
    paste0("detections <- 10^seq(log10(as.numeric(",input$detectionMin,")), log10(as.numeric(",input$detectionMax,")), length = as.numeric(",input$maxLength,"))")
  })
  coreHeatmapCodeParams3 <- reactive({
    paste0("gray <- gray(seq(0,1,length=5))")
  })
  coreHeatmapCodeParams4 <- reactive({
    paste0("coreplot <- plot_core(microbiome::transform(",input$datasetSample,", type = 'compositional'), plot.type = 'heatmap', colours = gray, prevalences = prevalences, detections = detections, min.prevalence = .5) + xlab('Detection Threshold (Relative Abundance (%))')")
  })
  coreHeatmapCodeParams5 <- reactive({
    paste0("plotly_build(coreplot)")
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
  
  # Subset the database #
  datasetSubsetInput <- reactive({
    subset1 <- prune_samples(sample_data(datasetInput())[[input$z1]] == input$v1, datasetInput())
    subset1 <- prune_samples(sample_data(subset1)[[input$z2]] == input$v2, subset1)
    subset1 <- prune_samples(sample_data(subset1)[[input$z3]] == input$v3, subset1)
    subset2 <- prune_samples(sample_data(subset1)[[input$z1]] == input$v1, subset1) %>% aggregate_taxa(level = input$v4)
    microbiome::transform(subset2, "compositional")
  })
  
    
  # Make plots
  communityPlotParams <- reactive ({
    theme_set(theme_classic(21))
    communityplot <-
      datasetSubsetInput() %>% plot_composition(sample.sort = "Bacteroidetes", otu.sort = "abundance") +
      scale_fill_manual(values = default_colors("Phylum")[taxa(datasetSubsetInput())])
    print(communityplot)
    #plotly_build(communityplot)
  }) #otu.sort and and sample.sort need to be selectable
  
  output$communityPlot <- renderPlot({
      communityPlotParams()  
  })
  
  # Make plot
  communityPlotGenusParams <- reactive({
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
        subtitle = "Subtitle", #make it user input
        caption = "Caption text." #Ditto
      ) +
      theme_ipsum(grid = "Y")
    print(compositionplot)
    #plotly_build(compositionplot)
  })
  
  output$communityPlotGenus <- renderPlot({
    communityPlotGenusParams()
  })
  
  # Barplot, averaged by z1Average group
  communityBarplotParams <- reactive({
    compplot <- plot_composition(datasetSubsetInput(), average_by = input$z1Average, plot.type = input$plotTypeZ1Average )
    print(compplot)
    #plotly_build(compplot)
  })
  output$communityBarplot <- renderPlot({
    communityBarplotParams()
  })
  
  # And top it off with a taxa prevalence plot
  communityPrevalenceParams <- reactive({
    subset1 <- prune_samples(sample_data(datasetInput())[[input$z1]] == input$v1, datasetInput())
    subset1 <- prune_samples(sample_data(subset1)[[input$z2]] == input$v2, subset1)
    subset1 <- prune_samples(sample_data(subset1)[[input$z3]] == input$v3, subset1)
    subset2 <- prune_samples(sample_data(subset1)[[input$z1]] == input$v1, subset1) %>% aggregate_taxa(level = input$v4Plot)
    data <- microbiome::transform(subset2, "compositional")
    prevplot <- plot_taxa_prevalence(data, input$v4 ) #Can be changed to whatever taxonomic rank the input files have
    plotly_build(prevplot)
  })
  
  output$communityPrevalence <- renderPlotly({
    communityPrevalenceParams()
  })  
  
  #Output code - community composition
  
  datasetSubsetInputCode1 <- reactive({
    paste0("subset1 <- prune_samples(sample_data(datasetInput())[[",input$z1,"]] == ",input$v1,", ",input$datasetSample,")")
  })
  datasetSubsetInputCode2 <- reactive({
    paste0("subset1 <- prune_samples(sample_data(subset1)[[",input$z2,"]] == ",input$v2,", subset1)")
  })
  datasetSubsetInputCode3 <- reactive({
    paste0("subset1 <- prune_samples(sample_data(subset1)[[",input$z3,"]] == ",input$v3,", subset1)")
  })
  datasetSubsetInputCode4 <- reactive({
    paste0("subset2 <- prune_samples(sample_data(subset1)[[",input$z1,"]] == ",input$v1,", subset1) %>% aggregate_taxa(level = ",input$v4,")")
  })
  datasetSubsetInputCode5 <- reactive({
    paste0("subset2 <- microbiome::transform(subset2, 'compositional'')")
  })
  
  communityPlotCodeParams1 <- reactive({
    paste0("theme_set(theme_classic(21))")
  })
  communityPlotCodeParams2 <- reactive({
    paste0("communityplot <- subset2 %>% plot_composition(sample.sort = 'Bacteroidetes', otu.sort = 'abundance') + scale_fill_manual(values = default_colors('Phylum')[taxa(subset2)])")
  })
  communityPlotCodeParams3 <- reactive({
    paste0("print(communityplot)")
  })
  
  communityPlotGenusCodeParams <- reactive({
    paste0("compositionplot <- plot_composition( subset2, taxonomic.level = 'Genus', sample.sort = 'nationality', x.label = 'nationality' ) + guides(fill = guide_legend(ncol = 1)) + scale_y_percent() + labs( x = 'Samples', y = 'Relative abundance (%)', title = 'Relative abundance data') + theme_ipsum(grid = 'Y')")
  })
  
  communityBarplotCodeParams1 <- reactive({
    paste0("compplot <- plot_composition(subset2, average_by = ",input$z1Average,", plot.type = ",input$plotTypeZ1Average,")")
  })
  communityBarplotCodeParams2 <- reactive({
    paste0("print(compplot)")
  })
  communityPrevalenceCodeParams1 <- reactive({
    paste0("subset1 <- prune_samples(sample_data(",input$datasetSample,"))[[",input$z1,"]] == ",input$v1,", ",input$datasetSample,")")
  })
  communityPrevalenceCodeParams2 <- reactive({
    paste0("subset1 <- prune_samples(sample_data(",input$datasetSample,"))[[",input$z2,"]] == ",input$v2,", subset1)")
  })
  communityPrevalenceCodeParams3 <- reactive({
    paste0("subset1 <- prune_samples(sample_data(",input$datasetSample,"))[[",input$z3,"]] == ",input$v3,", subset1)")
  })
  communityPrevalenceCodeParams4 <- reactive({
    paste0("subset2 <- prune_samples(sample_data(subset1)[[",input$z1,"]] == input$v1, subset1) %>% aggregate_taxa(level = ",input$v4Plot,")")
  })
  communityPrevalenceCodeParams5 <- reactive({
    paste0("data <- microbiome::transform(subset2, 'compositional')")
  })
  communityPrevalenceCodeParams6 <- reactive({
    paste0("prevplot <- plot_taxa_prevalence(data, input$v4 )")
  })
  communityPrevalenceCodeParams7 <- reactive({
    paste0("print(prevplot)")
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
  
  abundanceParams <- reactive({
    datatable(abundances(datasetInput(), transform = "compositional"), options = list(scrollX = TRUE))
  })
  output$abundanceTable <- renderDataTable({
    abundanceParams()
  })
  
  # Updating SelectInputs when database changes #
  observe({
    updateSelectInput(session, "x",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "x2", choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "y",
                      choices = colnames(alpha(datasetInput())))
  })
  
  
  # Merged table - generate and output #
  mergedTable <- reactive({
    merge(meta(datasetInput()), alpha(datasetInput()), all.y = TRUE)
  })
  viewParams <-
    reactive({
      #Issue with search -> search by specific term (like "male") doesn't work -> causes problem when searching between female (which gives all female samples) and male (which give both male and female samples). Still, a slightly dysfunctional search function is better than none at all
      datatable(mergedTable(),options = list(scrollX = TRUE))
    })
  
  output$view <- DT::renderDataTable({
    viewParams()
  })
  
  # # Violin plot #
  # violinPlotParams <- reactive({
  #   violin <- ggviolin(
  #     mergedTable(),
  #     x = input$x,
  #     y = input$y,
  #     fill = input$x,
  #     palette = c("#a6cee3", "#b2df8a", "#fdbf6f"),
  #     add = "boxplot"
  #   )
  #   plotly_build(violin)
  # })
  # 
  # output$violinPlot <- renderPlotly({
  #   violinPlotParams()
  # })
  
  # Richness Plot #
  richnessPlotParams <- reactive({
    richnessplot <- plot_richness(
      datasetInput(),
      x = input$x2,
      measures = input$richnessChoices,
      color = input$x2
    )
    plotly_build(richnessplot)
  })
  
  output$richnessPlot <- renderPlotly({
    richnessPlotParams()
  })
  
  #Output code - Alpha diversity + summary
  
  #Summary code
  summaryCodeParams <- reactive({
    paste0("summarize_phyloseq(",input$datasetSample,")")
  })
  output$summaryCode <- renderPrint({
    summaryCodeParams()
  })
  
  #DT - Evenness
  evennessTableCodeParams <- reactive({
    paste0("datatable(evenness(", input$datasetSample ,"),options = list(scrollX = TRUE))")
  })
  output$evennessTableCode <- renderPrint({
    evennessTableCodeParams()
  })
  #DT - Abundance
  abundanceTableCodeParams <- reactive({
    paste0("datatable(abundances(", input$datasetSample ,", transform = 'compositional'), options = list(scrollX = TRUE))")
  })
  output$abundanceTableCode <- renderPrint({
    abundanceTableCodeParams()
  })
  #DT - Metadata, Diversity measures
  mergedTableCodeParams <- reactive({
    paste0("merge(meta(",input$datasetSample,"), alpha(",input$datasetSample,"), all.y = TRUE)")
  })
  mergedTableCode <- renderPrint({
    mergedTableCodeParams()
  })
  
  
  #Richness Plot
  richnessPlotCodeParams <- reactive({
      dataset <- input$datasetSample
      x <- input$x2
      richness <- input$richnessChoices
      paste0("richnessplot <- plot_richness( ", dataset , ", " , "x = " , x , ", " , "measures = " , richness , ", " , "color = " , x , " )")
  })
  
  output$richnessPlotCode <- renderPrint({
    richnessPlotCodeParams()
  })

  
  ## Beta Diversity ##
  
  # Updating SelectInputs #
  observe({
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
    p <- phyloseq::plot_ordination(datasetInput(), ordinateData(), color = input$xb, label = input$yb ) + geom_point(size = input$geom.size)
    plotly_build(p)
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
        color = input$yb,
        label = input$yb
      ) + geom_point(size = input$geom.size2)
    plotly_build(splitOrdplot)
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
      ) + geom_point(size = input$geom.size3)
    plotly_build(taxaOrdplot)
  })
  
  output$taxaOrd <- renderPlotly({
    taxaOrdParams()
  })
  
  #Output code - Beta diversity
  ordinateDataCodeParams <- reactive({
    paste0("ordinateData <- ordinate(",
      "microbiome::transform(",input$datasetSample,", 'compositional'),",
      "method = ", input$ordinate.method,",",
      "distance = ", input$ordinate.distance,")"
    )
  })
  ordinateDataSplitCodeParams <- reactive({
    paste0("ordinateDataSplit <- ordinate(",
           "microbiome::transform(",input$datasetSample,", 'compositional'),",
           "method = ", input$ordinate.method2,",",
           "distance = ", input$ordinate.distance2,")"
    )
  })
  ordinateDataTaxaCodeParams <- reactive({
    paste0("ordinateDataTaxa <- ordinate(",
           "microbiome::transform(",input$datasetSample,", 'compositional'),",
           "method = ", input$ordinate.method3,",",
           "distance = ", input$ordinate.distance3,")"
    )
  })
  ordinatePlotCodeParams <- reactive({
    paste0("p <- phyloseq::plot_ordination(",input$datasetSample,", ordinateData, color = ",input$xb,", label = ",input$yb," ) + geom_point(size = ",input$geom.size,")"," plotly_build(p)")
  })
  splitOrdCodeParams <- reactive({
    paste0( "splitOrdplot <- plot_ordination(",input$datasetSample,",ordinateDataSplit, type = 'split', shape = ",input$xb,", color = ",input$yb,", label = ",input$yb," ) + geom_point(size = ",input$geom.size2,")", " plotly_build(splitOrdplot)")
  })
  taxaOrdCodeParams <- reactive({
    paste0( "taxaOrdplot <- plot_ordination(",input$datasetSample,", ordinateDataTaxa, type = 'taxa', color = ",input$zb,", label = ",input$xb,") + geom_point(size = ",input$geom.size3,")"," plotly_build(taxaOrdplot)"
    )
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
    subset1 <- prune_samples(sample_data(datasetInput())[[input$metadataLandscape1]] == input$metadataLandscapeValue1, datasetInput())
    subset <- prune_samples(sample_data(subset1)[[input$metadataLandscape2]] == input$metadataLandscapeValue2, subset1)
    microbiome::transform(subset, "compositional")
  })
  
  #PCA plot
  landscapePCAParams <- reactive({ 
    p <- plot_landscape(datasetInput(), method = "PCA", transformation = "clr") +
      labs(title = paste("PCA / CLR"))
    plotly_build(p)
  })
  
  output$landscapePCA <- renderPlotly({
    landscapePCAParams()
  })
  
  #PCoA/MDS
  landscapePCoAMDSParams <- reactive({
    p <- plot_landscape(datasetSubsetLandscape,
                        method = "PCoA", distance = "bray") +
      labs(title = paste("PCoA / Compositional / Bray-Curtis"))
    plotly_build(p)
  })
  
  output$landscapePCoAMDS <- renderPlotly({
    landscapePCoAMDSParams()
  })
  
  #NMDS
  landscapeOrdinationParams <- reactive({
    x <- datasetInput()
    quiet(x.ord <- ordinate(x, input$ordinateMethodLandscape, input$ordinateDistanceLandscape))
    proj <- phyloseq::plot_ordination(x, x.ord, justDF=TRUE)
    names(proj)[1:2] <- paste("Comp", 1:2, sep=".")
    p <- plot_landscape(proj[, 1:2], col = proj[[input$metadataLandscape1]], add.points = TRUE, legend = TRUE)
    plotly_build(p)
  })
  
  output$landscapeOrdination <- renderPlotly({
    landscapeOrdinationParams()
  })

  
  #t-SNE
  landscapeTSneParams <- reactive({
    p <- plot_landscape(datasetInput(), "t-SNE",
                   distance = "euclidean", transformation = "hellinger") +
      labs(title = paste("t-SNE / Hellinger / Euclidean"))
    plotly_build(p)
  })
  output$landscapeTSne <- renderPlotly({
    landscapeTSneParams()
  })
  
  #Abundance histogram (Absolute/Relative)
  landscapeAbundanceHistAbsParams <- reactive({
    # Visualize population densities for specific taxa
    p <- plot_density(datasetInput(), "Dialister") + ggtitle("Absolute abundance")
    plotly_build(p)
  })
  output$landscapeAbundanceHistAbs <- renderPlotly({
    landscapeAbundanceHistAbsParams()
  })
  
  landscapeAbundanceHistRelParams <- reactive({
    # Visualize population densities for specific taxa
    p <- plot_density(compositionalInput(), "Dialister", log10 = TRUE) +
      ggtitle("Relative abundance") +
      xlab("Relative abundance (%)")
    plotly_build(p)
  })
  output$landscapeAbundanceHistRel <- renderPlotly({
    landscapeAbundanceHistRelParams()
  })
  
  #Output code - Landscape analysis
  
  landscapePCACodeParams <- reactive({
    paste0("p <- plot_landscape(",input$datasetSample,", method = 'PCA', transformation = 'clr') + labs(title = paste('PCA / CLR')) plotly_build(p)")
  })
  
  landscapeOrdinationCodeParams <- reactive({
    paste0("x <- ",input$datasetSample,"quiet(x.ord <- ordinate(x, ",input$ordinateMethodLandscape,", ",input$ordinateDistanceLandscape,")) proj <- phyloseq::plot_ordination(x, x.ord, justDF=TRUE) names(proj)[1:2] <- paste('Comp0, 1:2, sep='.') p <- plot_landscape(proj[, 1:2], col = proj[[",input$metadataLandscape1,"]], add.points = TRUE, legend = TRUE) plotly_build(p)")
  })
  
  landscapeTSneCodeParams <- reactive({
    paste0("p <- plot_landscape(",input$datasetSample,", 't-SNE', distance = 'euclidean', transformation = 'hellinger') + labs(title = paste('t-SNE / Hellinger / Euclidean')) plotly_build(p)")
  })
  
  landscapeAbundanceHistAbsCodeParams <- reactive({
    paste0("p <- plot_density(",input$datasetSample,", 'Dialister') + ggtitle('Absolute abundance') plotly_build(p)")
  })
  landscapeAbundanceHistRelCodeParams <- reactive({
    paste0("p <- plot_density(microbiome::transform(",input$datasetSample,", 'compositional'), 'Dialister', log10 = TRUE) + ggtitle('Relative abundance') + xlab('Relative abundance (%)') plotly_build(p)")
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
  
  dmmModelCheckParams <- reactive({
    lplc <- sapply(dmmModelFit(), laplace) # AIC / BIC / Laplace
    aic  <- sapply(dmmModelFit(), AIC) # AIC / BIC / Laplace
    bic  <- sapply(dmmModelFit(), BIC) # AIC / BIC / Laplace
    plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
    lines(aic, type="b", lty = 2)
    lines(bic, type="b", lty = 3)
  })
  
  output$dmmModelCheck <- renderPlot({
    dmmModelCheckParams()
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
  
  taxaContributionPerComponentParams <- reactive({
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

  output$taxaContributionPerComponent <- renderPlotly({
    taxaContributionPerComponentParams()
  })
  
  #Output code - DMM community typing
  dmmModelFitCode1 <- reactive({
    paste0("taxa <- core_members(microbiome::transform(",input$datasetSample,",'compositional'), detection = ",input$detectionDMM,"/100, prevalence = ",input$prevalenceDMM,"/100)")
  })
  dmmModelFitCode2 <- reactive({
    paste0("pseq <- prune_taxa(taxa, ",input$datasetSample,")")
  })
  dmmModelFitCode3 <- reactive({
    paste0("dat <- abundances(pseq)")
  })
  dmmModelFitCode4 <- reactive({
    paste0("count <- as.matrix(t(dat))")
  })
  dmmModelFitCode5 <- reactive({
    paste0("fit <- mclapply(1:",input$maxModelsDMM,", dmn, count = count, verbose=TRUE)")
  })

  dmmModelCheckCodeParams1 <- reactive({
    paste0("lplc <- sapply(fit, laplace)") # AIC / BIC / Laplace
  })
  
  dmmModelCheckCodeParams2 <- reactive({
    paste0("aic  <- sapply(fit, AIC)") # AIC / BIC / Laplace
  })
  dmmModelCheckCodeParams3 <- reactive({
    paste0("bic  <- sapply(fit, BIC)") # AIC / BIC / Laplace
  })
  dmmModelCheckCodeParams4 <- reactive({
    paste0("plot(lplc, type='b', xlab='Number of Dirichlet Components', ylab='Model Fit')")
  })
  dmmModelCheckCodeParams5 <- reactive({
    paste0("lines(aic, type='b', lty = 2)")
  })
  dmmModelCheckCodeParams6 <- reactive({
    paste0("lines(bic, type='b', lty = 3)")
  })
  
  dmmModelBestFitCode1 <- reactive({
    paste0("lplc <- sapply(fit, laplace)") # AIC / BIC / Laplace
  })
  
  dmmModelBestFitCode2 <- reactive({
    paste0("aic  <- sapply(fit, AIC)") # AIC / BIC / Laplace
  })
  dmmModelBestFitCode3 <- reactive({
    paste0("bic  <- sapply(fit, BIC)") # AIC / BIC / Laplace
    
  })
  dmmModelBestFitCode4 <- reactive({
    paste0("bestFit <- fit[[which.min(lplc)]]")
    
  })

  dmmParametersCode <- reactive({
    paste0("datatable(mixturewt(bestFit))")
  })
  sampleAssignmentsCode <- reactive({
    paste0("datatable(apply(mixture(bestFit), 1, which.max))")
  })
  
  taxaContributionPerComponentCodeParams1 <- reactive({
    paste0( "for (k in seq(ncol(fitted(bestFit)))) {")
  })
  
  taxaContributionPerComponentCodeParams2 <- reactive({
    paste0( "d <- melt(fitted(bestFit))")
  })
  taxaContributionPerComponentCodeParams3 <- reactive({
    paste0( "colnames(d) <- c('OTU', 'cluster', 'value')")
  })
  taxaContributionPerComponentCodeParams4 <- reactive({
    paste0("d <- subset(d, cluster == k) %>% arrange(value) %>% mutate(OTU = factor(OTU, levels = unique(OTU))) %>% filter(abs(value) > quantile(abs(value), 0.8))")
  })
  taxaContributionPerComponentCodeParams5 <- reactive({
    paste0("p <- ggplot(d, aes(x = OTU, y = value)) + geom_bar(stat = 'identity') + coord_flip() + labs(title = paste('Top drivers: community type', k))")
  })
  taxaContributionPerComponentCodeParams6 <- reactive({
    paste0("print(p)}")
  })
  
  ###########################
  ## Statistical analysis ###
  ###########################
  
  ## PERMANOVA ##
  #Update metadata column#
  observe({
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
                        data = meta, permutations = permnumber, method = "bray", parallel = getOption("mc.cores")
    )
  })
  
  densityPlotParams <- reactive({
    p <- plot_landscape(compositionalInput(), method = input$permanovaMethod, distance = input$permanovaDistanceMethod, col = input$permanovaColumn, size = input$permanovaPlotSize)
    plotly_build(p)
  })
  output$densityPlot <- renderPlotly({
    densityPlotParams()
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
    print(p)
  })
  
  output$topFactorPlot <- renderPlot({
    topFactorPlotParams()
  })
  
  netPlotParams <- reactive({
    n <- make_network(compositionalInput(), type = "otu", distance = input$permanovaDistanceMethodNet)
    p <- plot_network(n)
    plotly_build(p)  
  })
  
  output$netPlot <- renderPlotly({
    netPlotParams()
  })
  
  output$permaHeatmap <- renderPlotly({
    xp <- plot_heatmap(compositionalInput(), distance = ordinate(compositionalInput(), distance = input$permanovaDistanceMethodHeat), method = input$permanovaMethodHeat)
    plotly_build(p)
  })
  
  #Output code - PERMANOVA
  
  permanovaCode1 <- reactive({
    paste0("otu <- abundances(microbiome::transform(",input$datasetSample,",'compositional'))")
  })
  permanovaCode2 <- reactive({
    paste0("meta <- meta(microbiome::transform(",input$datasetSample,",'compositional'))")
  })
  permanovaCode3 <- reactive({
    paste0("permnumber <- ", input$permanovaPermutationsP)
  })
  permanovaCode4 <- reactive({
    paste0("metadata <- ", input$permanovaColumnP)
  })
  permanovaCode5 <- reactive({
    paste0("permanova <- adonis(t(otu) ~ meta[[metadata]], data = meta, permutations = permnumber, method = 'bray', parallel = getOption('mc.cores'))")
  })

  densityPlotCodeParams1 <- reactive({
    paste0("p <- plot_landscape(microbiome::transform(",input$datasetSample,",'compositional'), method = ",input$permanovaMethod,", distance = ",input$permanovaDistanceMethod,", col = ",input$permanovaColumn,", size = ",input$permanovaPlotSize,")")
  })
  densityPlotCodeParams2 <- reactive({
    paste0("plotly_build(p)")
  })
  
  pValueCode <- reactive({
    paste0("as.data.frame(permanova$aov.tab)")
  })
  
  homogenietyCodeParams1 <- reactive({
    paste0("otu <- abundances(microbiome::transform(",input$datasetSample,",'compositional'))")
  })
  
  homogenietyCodeParams2 <- reactive({
    paste0("meta <- meta(microbiome::transform(",input$datasetSample,",'compositional'))")
  })
  homogenietyCodeParams3 <- reactive({
    paste0("dist <- vegdist(t(otu))")
  })
  homogenietyCodeParams4 <- reactive({
    paste0("metadata <- ", input$permanovaColumnP)
  })
  homogenietyCodeParams5 <- reactive({
    paste0("homogeniety <- anova(betadisper(dist, meta[[metadata]]))")
  })
  
  homogenietyCodeParams6 <- reactive({
    paste0("datatable(homogeniety)")
  })
  
  topFactorPlotCodeParams1 <- reactive({
    paste0("p <- barplot(sort(top.coef), horiz = T, las = 1, main = 'Top taxa')")
  })
  topFactorPlotCodeParams2 <- reactive({
    paste0("print(p)")
  })
  
  netPlotCodeParams1 <- reactive({
    paste0("n <- make_network(microbiome::transform(",input$datasetSample,",'compositional'), type = 'otu', distance = ",input$permanovaDistanceMethodNet,")")
  })
  
  netPlotCodeParams2 <- reactive({
    paste0("p <- plot_network(n)")
  })
  
  netPlotCodeParams3 <- reactive({
    paste0("plotly_build(p)")  
  })
  
  permaHeatmapCode1 <- reactive({
    paste0("xp <- plot_heatmap(microbiome::transform(",input$datasetSample,",'compositional'), distance = ordinate(microbiome::transform(",input$datasetSample,",'compositional'), distance = ",input$permanovaDistanceMethodHeat,"), method = ",input$permanovaMethodHeat,")")
  })
  
  permaHeatmapCode2 <- reactive({
    paste0("plotly_build(p)")
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
  
# Alpha Diversity Results #
  # 
  # output$alphaCodeSummary <- renderPrint({
  #   paste0("library(\"phyloseq\")")
  #   paste0("library(\"microbiome\")")
  #   paste0("summarize_phyloseq(" + deparse(substitute(datasetInput())) + ")")
  # })


  output$downloadReportAlpha <- downloadHandler(
    filename = function() {
      paste('alpha-report', sep = '.', switch(
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
  
  #Community Landscape#
  output$downloadReportLandscape <- downloadHandler(
    filename = function() {
      paste('landscape-report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html'
      ))
    },
    content = function(file) {
      src <- normalizePath('landscape_report.Rmd')
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'landscape_report.Rmd', overwrite = TRUE)
      
      out <- rmarkdown::render('landscape_report.Rmd',
                               switch(input$format,
                                      PDF = pdf_document(), 
                                      HTML = html_document() 
                               ))
      file.rename(out, file)
    }
  )
  
  #DMM Community Typing#
  output$downloadReportDMM <- downloadHandler(
    filename = function() {
      paste('dmm-report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html'
      ))
    },
    content = function(file) {
      src <- normalizePath('dmm_report.Rmd')
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'dmm_report.Rmd', overwrite = TRUE)
      
      out <- rmarkdown::render('dmm_report.Rmd',
                               switch(input$format,
                                      PDF = pdf_document(), 
                                      HTML = html_document() 
                               ))
      file.rename(out, file)
    }
  )
  
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
