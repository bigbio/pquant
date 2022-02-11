## The code is written by Douer. 
## This is a Shiny application  for quantms pipeline data downstream analysis.
## You can run the application by clicking the 'Run App' button above.

library(shiny)
library(shinydashboard)
library(DT)
library(MSstats)
library(pheatmap)
library(rhandsontable)
library(shinyjs)
library(shinyWidgets)
library(ggplot2)
library(tidyverse)
library(IDPmisc)
library(data.table)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(proteus)
library(MSstatsTMT)
library(purrr)

library(org.Hs.eg.db)
library(org.Sc.sgd.db)
library(org.Rn.eg.db)
library(org.EcK12.eg.db)

# Source getPlots functions
source("./app/getSelector.R")
source("./app/volcano_live.R")
source("./app/FID_live.R")
source("./app/proteus_plots.R")
source("./app/TMT_plots.R")

  
ui <- dashboardPage(
    dashboardHeader(title = "pquantR"),
    dashboardSidebar(
                    useShinyjs(),
                    
                    # home page
                    conditionalPanel(condition = "input.main_tabs == 'home_condition'",
                                     sidebarMenu(
                                                 withTags({
                                                           div(
                                                               style = "text-align:center",br(),
                                                               h4(i(class = 'fab fa-github'), a(href="https://github.com/bigbio/pquant", "Github address")),br(),br(),
                                                               h5(i(class = 'fas fa-user-tie'),"Contact us"),
                                                               h5(i(class = 'fa fa-envelope-o'), "douerww@gmail.com"),
                                                               h5(i(class = 'fa fa-envelope-o'), "ypriverol@gmail.com"))
                                                          }))
                    ),
                    
                    # fileinput side bar menu
                    conditionalPanel(condition = "input.main_tabs == 'fileinput_condition'",
                                     sidebarMenu(
                                                 menuItem("Upload data",
                                                          startExpanded = TRUE, br(),
                                                          icon = icon("upload"),
                                                          h5('Choose the upload data type:'),
                                                          menuItem("URL",br(),
                                                                   actionButton(inputId = "inputData_URL",
                                                                                label = "Get data",
                                                                                width = 200,
                                                                                icon = icon("link"),
                                                                                style ="display: block; margin: 0 auto; width: 200px;color: black;"),br()),
                                                          menuItem("Local",
                                                                   fileInput(inputId = "inputData",
                                                                             label = "Upload .csv(.gz) / .Rdata",
                                                                             multiple = FALSE,
                                                                             accept = c(".csv",
                                                                                        ".gz",
                                                                                        ".RData"))),br(),
                                                          div(style = "text-align:left", 
                                                              "Pre: Proteins to Gene Accessions"),  
                                                          switchInput(inputId = "user_choose_pre_pro2gene",
                                                                      onLabel = "Yes", offLabel = "No",
                                                                      onStatus = "primary", offStatus = "default",
                                                                      value = FALSE),
                                                          helpText("You could perform gene expression analysis coming from the peptide info by choose 'Yes'."),
                                                          tags$hr(),
                                                          actionButton(inputId = "uploaddata_reset",
                                                                       label = "Reset all data",
                                                                       width = 200,
                                                                       icon = icon("redo")),br()
                    
                                                 ),
                                                 
                                                 menuItem("Parameters selection (LFQ)", #startExpanded = TRUE,
                                                          icon = icon("cogs"),br(),
                                                          radioButtons(inputId = "user_choose_logTrans",
                                                                       label = "Choose logarithm transformation",
                                                                       choices = c("2" = "2",
                                                                                   "10" = "10"),
                                                                       selected = "2"),
                                                         
                                                          radioButtons(inputId = "user_choose_normalization",
                                                                       label = "Choose normalization",
                                                                       choices = c("equalizeMedians" = "equalizeMedians",
                                                                                   "quantile" = "quantile",
                                                                                   "globalStandards" = "globalStandards",
                                                                                   "no normalization" = "FALSE"),
                                                                       selected = "equalizeMedians"),
                                                         
                                                          radioButtons(inputId = "user_choose_summaryMethod",
                                                                       label = "Choose summary method",
                                                                       choices = c("Tukey’s median polish(TMP)" = "TMP",
                                                                                   "linear mixed model" = "linear"),
                                                                       selected = "TMP"),
                                                                      
                                                          numericInput(inputId = "user_choose_maxQuantileforCensored",
                                                                       label = "Choose maxQuantile for deciding censored missing values",
                                                                       value = 0.999),br()
                                                 ),
                                                 
                                                 menuItem("Parameters selection (TMT)",
                                                          icon = icon("cogs"),br(),
                                                          radioButtons(inputId = "user_choose_TMT_method",
                                                                       label = "Choose logarithm transformation",
                                                                       choices = c("MSstats" = "msstats",
                                                                                   "MedianPolish" = "MedianPolish",
                                                                                   "Median" = "Median",
                                                                                   "LogSum" = "LogSum"),
                                                                       selected = "msstats"),
                                                          radioButtons(inputId = "user_choose_TMT_global_norm",
                                                                       label = "Global median normalization on peptide level data",
                                                                       choices = c("Yes" = "TRUE",
                                                                                   "No" = "FALSE"),
                                                                       selected = "TRUE"),
                                                          numericInput(inputId = "user_choose_TMT_maxQuantileforCensored",
                                                                       label = "Choose maxQuantile for deciding censored missing values",
                                                                       value = NULL),br()
                                                 ),
                                                 
                                                 
                                                 menuItem("Start preprocessing",
                                                          icon = icon("play"),br(),
                                                          div(style = "text-align:left", 
                                                              "Save preprocessed data to .RData?"),     
                                                          switchInput(inputId = "SaveProjectDataControl",
                                                                      onLabel = "Yes", offLabel = "No",
                                                                      onStatus = "primary", offStatus = "default",
                                                                      value = FALSE),
                                                          textInput(inputId = "SaveProjectDataName",
                                                                    label = "Project name"),
                                                          textInput(inputId = "SaveProjectFolderPath",
                                                                    label = "Folder path"),br(),
                                                          helpText('A 40Mb file takes about 6 mins.'),
                                                          actionButton(inputId = "start_preprocess",
                                                                       label = "Start Preprocessing",
                                                                       icon = icon("play-circle"),
                                                                       style ="display: block; margin: 0 auto; width: 200px;color: black;"),br()
                                                 ),
                                                 
                                                 
                                                 menuItem("Data process plots",
                                                          icon = icon("chart-bar"),br(),
                                                          menuItem("Profile & Condition plot",br(),
                                                                   uiOutput('initialData_profilecondition_selector'),br(),
                                                                   actionButton(inputId = "initialData_profilecondition_Render",
                                                                                label = "Render Plot",
                                                                                icon = icon("play-circle"),
                                                                                style ="display: block; margin: 0 auto; width: 200px;color: black;"),br()
                                                          ),br(),
                                                          menuItem("Quality control plot",br(),
                                                                   uiOutput('initialData_qualitycontrol_selector'),br(),
                                                                   actionButton(inputId = "initialData_qualitycontrol_Render",
                                                                                label = "Render Plot",
                                                                                icon = icon("play-circle"),
                                                                                style ="display: block; margin: 0 auto; width: 200px;color: black;"),br()
                                                          ),br()
                                                 ),
                                                 
                                                 menuItem("Download data",
                                                          icon = icon("download"),br(),
                                                          textInput(inputId = "download_name",
                                                                    label = "File name",
                                                                    placeholder = "Preprocessed data"),br(),
                                                          downloadButton(outputId = "download_Render",
                                                                         label = "Download",
                                                                         style="display: block; margin: 0 auto; width: 200px;color: black;"),br()
                                                 )
                                     )
                    ),
                    
                    
                    # default method sidebar
                    conditionalPanel(condition = "input.main_tabs == 'default_method_condition'",br(),
                                     sidebarMenu(
                                                 menuItem("Model-based QC plots",
                                                          icon = icon("chart-area"),br(),
                                                          menuItem("Residual & Normal Q-Q plot",br(),
                                                                   uiOutput('default_method_residual_qq_selector'),br(),
                                                                   actionButton(inputId = "default_method_residual_qq_Render",
                                                                                label = "Render Plot",
                                                                                icon = icon("play-circle"),
                                                                                style ="display: block; margin: 0 auto; width: 200px;color: black;"
                                                                   ),br()
                                                          ),br()
                                                  ),
                                           
                                                 menuItem("Group comparison plots",
                                                          icon = icon("chart-area"),br(),
                                                          menuItem("Volcano plot",br(),
                                                                   uiOutput('default_method_volcano_selector'),br(),
                                                                   actionButton(inputId = "default_method_volcano_Render",
                                                                                label = "Render Plot",
                                                                                icon = icon("play-circle"),
                                                                                style ="display: block; margin: 0 auto; width: 200px;color: black;"),br()
                                                          ),br(),
                                                          menuItem("Heatmap",br(),
                                                                   actionButton(inputId = "default_method_heatmap_Render",
                                                                                label = "Render Plot",
                                                                                icon = icon("play-circle"),
                                                                                style ="display: block; margin: 0 auto; width: 200px;color: black;"),br(),
                                                                   helpText("When there is only one set of comparisons, no heat map is generated.")
                                                          ),br(),
                                                          menuItem("Comparison plot",
                                                                   br(),
                                                                   uiOutput('default_method_comparison_selector'),
                                                                   br(),
                                                                   actionButton(inputId = "default_method_comparison_Render",
                                                                                label = "Render Plot",
                                                                                icon = icon("play-circle"),
                                                                                style ="display: block; margin: 0 auto; width: 200px;color: black;"),br()
                                                          ),br()
                                                 )
                                     )
                    ),
            
                    
                    # dynamic plots page
                    conditionalPanel(condition = "input.main_tabs == 'dynamic_volcano'",br(),
                                     sidebarMenu(
                                                 menuItem("Assign condition",br(),
                                                          icon = icon("clone"),
                                                          div(style = "text-align:center", "Click to render annotation table"),
                                                          actionButton(inputId = "start_anno",
                                                                       label = "Start",
                                                                       icon = icon("play-circle"),
                                                                       width = 200),
                                                          tags$hr(),
                                                          div(style = "text-align:center", "Your input in the second table",
                                                              br(), "will be used for stastical",
                                                              br(), "comparisons.",
                                                              br(),br(), "Example:"),
                                                          img(src = 'annotation_example.png', width = '100%'),
                                                          div(style = "text-align:center", "once done renaming press",br(),
                                                              "\'submit\' to lock in the annotation"),
                                                          actionButton(inputId = "submit_anno",
                                                                       label = "Submit",
                                                                       icon = icon("file-import"),
                                                                       style ="display: block; margin: 0 auto; width: 200px;color: black;"),
                                                          tags$hr(),
                                                          actionButton(inputId = "annoReset",
                                                                       label = "Reset annotation",
                                                                       width = 200,
                                                                       icon = icon("redo"))
                                                 ),
                                                 menuItem("Volcano plot",
                                                          icon = icon("chart-area"),br(),
                                                          uiOutput('dynamic_volcano_selector'),
                                                          actionButton(inputId = "dynamic_volcano_Render",
                                                                       label = "Render Plot",
                                                                       icon = icon("play-circle"),
                                                                       style ="display: block; margin: 0 auto; width: 200px;color: black;"),br()
                                                 ),
                                                 menuItem("Fold-change-intensity plot",
                                                          icon = icon("chart-bar"),br(),
                                                          uiOutput('dynamic_FID_selector'),
                                                          actionButton(inputId = "dynamic_FID_Render",
                                                                       label = "Render Plot",
                                                                       icon = icon("play-circle"),
                                                                       style ="display: block; margin: 0 auto; width: 200px;color: black;"),br()
                                                 )
                                     )
                    ),
                    
                    # Proteus page
                    conditionalPanel(condition = "input.main_tabs == 'proteus_condition'",br(),
                                     div(style = "text-align:center",
                                         "Use this part function must",br(),
                                         "submit annotation in 'Dynamic plots'",br()),
                                     tags$hr(),
                                     sidebarMenu(
                                                 menuItem("LFQ",br(),
                                                          menuItem("Peptide data",br(),
                                                                   icon = icon("images"),
                                                                   actionButton(inputId = "proteus_peptide_Render",
                                                                                label = "Render Plots",
                                                                                icon = icon("play-circle"),
                                                                                style ="display: block; margin: 0 auto; width: 200px;color: black;"),br()),br(),
                                                          menuItem("Protein data",br(),
                                                                   icon = icon("images"),
                                                                   actionButton(inputId = "proteus_protein_Render",
                                                                                label = "Render Plots",
                                                                                icon = icon("play-circle"),
                                                                                style ="display: block; margin: 0 auto; width: 200px;color: black;"),br()),br()
                                                 ),
                                                 menuItem("TMT",br(),
                                                          menuItem("Peptide data",br(),
                                                                   icon = icon("images"),
                                                                   actionButton(inputId = "proteus_TMT_peptide_Render",
                                                                                label = "Render Plots",
                                                                                icon = icon("play-circle"),
                                                                                style ="display: block; margin: 0 auto; width: 200px;color: black;"),br()),br(),
                                                          menuItem("Protein data",br(),
                                                                   icon = icon("images"),
                                                                   actionButton(inputId = "proteus_TMT_protein_Render",
                                                                                label = "Render Plots",
                                                                                icon = icon("play-circle"),
                                                                                style ="display: block; margin: 0 auto; width: 200px;color: black;"),br()),br()
                                                 )
                                                 
                                     )
                    )
    ),
  
    
    dashboardBody(
                  tabsetPanel(  # a layout
                              id = 'main_tabs',
                              
                              # home tab
                              tabPanel(title = 'Home',
                                       value = 'home_condition',
                                       icon = icon("home"),
                                       fluidPage(
                                                 fluidRow(h2("pquantR: a shiny-based application for downstream analysis of proteomics data processed"),br(),
                                                          p("As a powerful tool to analyze protein expression, proteomics becomes more and more important for biomedical researchers to convert large amounts of experimental data into biological knowledge. However, there is a major challenge in this field: data processing is rather complicated that may confuses researchers or raises a hard learning cost. To address this, we have created pquantR, a free open-source shiny-based analytics platform that allows users to analyze data from label-free quantification workflow, which is common in proteomics relative researches. The uploading data of this application comes from our 'quantms' pipeline, which used OpenMS and MSstats, with feature quantification, feature summarization, quality control and group-based statistical analysis. Besides, pquantR has a user-friendly interface, with automated data preprocessing pipeline which is easy to operate. A few clicks are needed to quickly analyze statistical data. Furthermore, multiple kinds of plots are supported and integrated in one application: visualization of processed data, potential systematic differences, protein differential expression and quality control plots. Users also could interact with volcano plot to select and inspect the detailed information."),
                                                          br(),br(),br(),
                                                          div(img(src = 'pipeline.png', width = '50%'), style="text-align: center;",),
                                                          br(),br(),br(),
                                                          h4("Data not ready?", a(href="./data/out_msstats_example.csv", download="out_msstats_example.csv", "Get demo data here"))),
                                       )
                              ),
                      
                      
                              # data tab
                              tabPanel(title = 'Data',
                                       value = 'fileinput_condition',
                                       icon = icon("table"),
                                       fluidPage(
                                                 tabBox(# No title
                                                        id = "initialData_tabbox", selected = "initialData_tabbox_data", width = 12,
                                                        tabPanel(
                                                                 title = ".csv data", value = "initialData_tabbox_data",
                                                                 fluidRow(
                                                                          column(11,
                                                                                 DT::DTOutput("contents", width = '90%')))
                                                        ),
                                                        tabPanel(
                                                                 title = "Data process plots", value = "initialData_tabbox_plot",
                                                                 fluidRow(
                                                                          column(6,
                                                                                 plotOutput("initialData_profilePlot_out")),
                                                                          column(6,
                                                                                 plotOutput("initialData_qualitycontrolPlot_out"))),br(),
                                                                 fluidRow(
                                                                          column(6,
                                                                                 plotOutput("initialData_conditionPlot_out")))
                                                        ),
                                                        tabPanel(
                                                                 title = "Preprocessed data", value ="initialData_tabbox_preprocesseddata",
                                                                 fluidRow(
                                                                          column(11,
                                                                                 DT::DTOutput("preprocessedData", width = '80%')))
                                                        )
                                                 )
                                       )
                                       
                              ),
                      
                      
                              # default method condition tab
                              tabPanel(title = 'MSstats method',
                                       value = 'default_method_condition',
                                       icon = icon("chart-area"),
                                       fluidPage(
                                                 tabBox(# No title
                                                        id = "msstats_tabbox", selected = "msstats_tabbox_modelbased", width = 12,
                                                        tabPanel(
                                                                 title = "Model-based QC plots", value = "msstats_tabbox_modelbased",
                                                                 fluidRow(
                                                                          column(6,
                                                                                 plotOutput("default_method_residualPlot_out")),
                                                                          column(6,
                                                                                 plotOutput("default_method_qqPlot_out")))
                                                        ),
                                                        tabPanel(
                                                                 title = "Group comparison plots", value = "msstats_tabbox_groupcomparison",
                                                                 fluidRow(
                                                                          column(6,
                                                                                 plotOutput("default_method_volcano_out")),
                                                                          column(6,
                                                                                 plotOutput("default_method_heat_out"))),br(),
                                                                 fluidRow(
                                                                          column(6,
                                                                                 plotOutput("default_method_comparisonPlot_out")))
                                                        )
                                                 )
                                       ),
                                       
                              ),
                      
                      
                      
                              # dynamic volcano tab
                              tabPanel(title = 'Dynamic Plots',
                                       value = 'dynamic_volcano',
                                       icon = icon("project-diagram"),
                                       fluidPage(
                                                 tabBox(# No title
                                                        id = "dynamic_tabbox1", selected = "dynamic_tabbox_data", width = 12,
                                                        tabPanel(
                                                                 title = "Data", value = "dynamic_tabbox_data",
                                                                 fluidRow(
                                                                          column(width = 9,
                                                                                 DT::DTOutput("dynamic_evidence_contents")),
                                                                          column(width = 3,
                                                                                 rHandsontableOutput("dynamic_define_metadata")))
                                                        ),
                                                        tabPanel(
                                                                title = "Volcano plot", value = "dynamic_tabbox_plot",
                                                               
                                                                ### from Proteus: live.R
                                                                fluidRow(
                                                                         column(6,
                                                                                plotOutput("dynamic_plotVolcano_out", height = "600px", width = "80%", brush = "plot_brush",hover="plot_hover"),
                                                                                DT::dataTableOutput("dynamic_significanceTable_out")),
                                                                         column(6,
                                                                                fluidRow(
                                                                                         column(4,
                                                                                                radioButtons("intensityScale","Intesity Scale:",choices = c("Linear scale" = "","Log scale"="Log"),inline = TRUE))),
                                                                                fluidRow(
                                                                                         column(4,
                                                                                                fluidRow(htmlOutput("dynamic_gap_out")),
                                                                                                fluidRow(plotOutput("dynamic_jitterPlot_out", height = "300px",width = "100%"))),
                                                                                         column(2,
                                                                                                fluidRow(tableOutput("dynamic_replicateTable_out")))))
                                                                ),
                                                                tags$hr(),
                                                                # Show main protein table
                                                                fluidRow(
                                                                         column(width = 12,
                                                                                DT::dataTableOutput("dynamic_allProteinTable_out")))
                                                       ),
                                                       tabPanel(title = "Fold-change-intensity plot (not support TMT data)", value = "dynamic_tabbox_plot2",
                                                                
                                                                ### from Proteus: live.R
                                                                fluidRow(
                                                                  column(6,
                                                                         plotOutput("dynamic_plotFID_out", height = "600px", width = "80%", brush = "plot_brush",hover="plot_hover"),
                                                                         DT::dataTableOutput("dynamic_FID_significanceTable_out")),
                                                                  column(6,
                                                                         fluidRow(
                                                                           column(4,
                                                                                  radioButtons("dynamic_FID_intensityScale","Intesity Scale:",choices = c("Linear scale" = "","Log scale"="Log"),inline = TRUE))),
                                                                         fluidRow(
                                                                           column(4,
                                                                                  fluidRow(htmlOutput("dynamic_FID_gap_out")),
                                                                                  fluidRow(plotOutput("dynamic_FID_jitterPlot_out", height = "300px",width = "100%"))),
                                                                           column(2,
                                                                                  fluidRow(tableOutput("dynamic_FID_replicateTable_out")))))
                                                                ),
                                                                tags$hr(),
                                                                # Show main protein table
                                                                fluidRow(
                                                                  column(width = 12,
                                                                         DT::dataTableOutput("dynamic_FID_allProteinTable_out")))
                                                       )
                                                   
                                                 )
                                       )
                              ),
                              
                              
                              # Proteus tab
                              tabPanel(title = "Proteus method",
                                       value = 'proteus_condition',
                                       icon = icon("chart-area"),
                                       fluidPage(
                                                 tabBox(# No title
                                                        id = "proteus_tabbox", selected = "proteus_tabbox_peptide", width = 12,
                                                        tabPanel(
                                                                 title = "Peptide data (LFQ)", value = "proteus_tabbox_peptide",
                                                                 fluidRow(
                                                                          column(6,
                                                                                 verbatimTextOutput("proteus_peptide_summary_out")),
                                                                          column(6,
                                                                                 plotOutput("proteus_peptide_number_out"))),br(),
                                                                 fluidRow(
                                                                          column(6,
                                                                                 plotOutput("proteus_peptide_distance_out")),
                                                                          column(6,
                                                                                 plotOutput("proteus_peptide_jaccard_out"))),br(),
                                                                 fluidRow(
                                                                          column(6,
                                                                                 plotOutput("proteus_peptide_pca_out")),
                                                                          column(6,
                                                                                 plotOutput("proteus_peptide_clustering_out")))
                                                        ),
                                                        tabPanel(
                                                                 title = "Protein data (LFQ)", value = "proteus_tabbox_protein",
                                                                 fluidRow(
                                                                          column(6,
                                                                                 verbatimTextOutput("proteus_protein_summary_out")),
                                                                          column(6,
                                                                                 plotOutput("proteus_protein_normalization_median_out"))),br(),
                                                                 fluidRow(
                                                                          column(6,
                                                                                 plotOutput("proteus_protein_mean_out")),
                                                                          column(6,
                                                                                 plotOutput("proteus_protein_clustering_out")))
                                                        ),
                                                        tabPanel(
                                                                 title = "Peptide data (TMT)", value = "proteus_TMT_tabbox_peptide",
                                                                 fluidRow(
                                                                          column(6,
                                                                                 verbatimTextOutput("proteus_TMT_peptide_summary_out")),
                                                                          column(6,
                                                                                 plotOutput("proteus_TMT_peptide_number_out")))
                                                        ),
                                                        tabPanel(
                                                                 title = "Protein data (TMT)", value = "proteus_TMT_tabbox_protein",
                                                                 fluidRow(
                                                                          column(6,
                                                                                 verbatimTextOutput("proteus_TMT_protein_summary_out")),
                                                                          column(6,
                                                                                 plotOutput("proteus_TMT_protein_normalization_median_out"))),br(),
                                                                 fluidRow(
                                                                          column(6,
                                                                                 plotOutput("proteus_TMT_protein_clustering_out")))
                                                        )
                                                 )
                                       )
                              )
                      
                              # enrichment analysis tab
                              #tabPanel(title = "Enrichment analysis",value = "enrichment_analsis")
                      
                      
                    )
                     
        )
      
)



# Define server logic ----
server <- function(input, output, session) {
    options(shiny.maxRequestSize=520*1024^2)
    options(ggrepel.max.overlaps = Inf)
    
    env <<- NULL
  
    prePquant <- reactiveValues(DDA2009.proposed = NULL,
                                DDA2009.comparisons = NULL,
                                inputData_rdata_out_msstats = NULL,
                                inputData_URL_out_msstats = NULL,
                                quant.msstats = NULL,
                                test.pairwise = NULL)
    
    renderCheck <- reactiveValues(initialData_qualitycontrolPlot = 0,
                                  initialData_profilecondition = 0,
                                  default_method_residual_qq = 0,
                                  default_method_volcano = 0,
                                  default_method_heatmap = 0,
                                  default_method_comparison = 0,
                                  dynamic_volcano = 0,
                                  dynamic_FID = 0,
                                  proteus_peptide = 0,
                                  proteus_protein = 0,
                                  proteus_pep_pro = 0,
                                  proteus_TMT_peptide = 0,
                                  proteus_TMT_protein = 0)
    
    volcanoLive <- reactiveValues(pdat = NULL,
                                  pepdat = NULL,
                                  prolfq_pepdat = NULL,
                                  prolfq_prodat.med = NULL,
                                  res = NULL,
                                  res_FID = NULL,
                                  max_points = NULL,
                                  protmt_pepd = NULL,
                                  protmt_prod = NULL)
    
    dataControl <- reactiveValues(annoStart = 0,
                                  annoSubmit = 0,
                                  inputData_state = NULL,
                                  inputData_rdata = NULL,
                                  inputData_URL_state = NULL,
                                  preprocess_judge = 0,
                                  db_type = NULL,
                                  data_type = "NULL")

    # -------------input---------------
    
    #### .csvfile, inputData part 
    
    observeEvent(input$inputData_URL,{
      dataControl$inputData_state <- "uploaded"
      dataControl$inputData_URL_state <- "uploaded"
      
      prePquant$inputData_URL_out_msstats <- data.table::fread("http://ftp.pride.ebi.ac.uk/pride/data/proteomes/proteogenomics/differential-expression/RPMID25238572.1-cell-lines/proteomics_lfq/out_msstats.csv") %>% as.data.frame
    })
    

    observeEvent(input$inputData,{
        dataControl$inputData_state <- "uploaded"
        if(input$inputData$type == ""){
            dataControl$inputData_rdata <- "uploaded"
        }
    })
    
    inputdf <- reactive({
        if(is.null(dataControl$inputData_state)){
            return(NULL)
        }
        else {
            if(is.null(dataControl$inputData_URL_state) != TRUE){
                fileData <- prePquant$inputData_URL_out_msstats
            }else{
                fileName <- input$inputData$name
                tmp <- strsplit(fileName, ".", fixed = TRUE)
                fileType <- tmp[[1]][length(tmp[[1]])]
                if(fileType == "csv" | fileType == "gz"){  #upload a csv file
                    fileData <- data.table::fread(input$inputData$datapath) %>% as.data.frame
                }
                else{
                    fileData <- prePquant$inputData_rdata_out_msstats}
            }
            return(fileData)
        }
    })
    
    ### load .rdata
    observeEvent(input$inputData, {
        if(is.null(dataControl$inputData_rdata) != TRUE){
            sendSweetAlert(
              session = session,
              title = "Success",
              text = "You are loading 'DDA2009.RData', please wait until done to the next step.",
              type = "success", 
              closeOnClickOutside = TRUE,
              width = 400
            )
            
            progress <- shiny::Progress$new()
            on.exit(progress$close())
            
            progress$set(message = "Begin to load data, please wait...", value = 0.5)
           
            n.env <- new.env()
            env <<- n.env
            load(input$inputData$datapath, envir = n.env)
            
            prePquant$inputData_rdata_out_msstats <- n.env$inputData_rdata_out_msstats
            if(sum("Mixture" == colnames(prePquant$inputData_rdata_out_msstats)) == 1){
                prePquant$test.pairwise <- n.env$test.pairwise
                prePquant$quant.msstats <- n.env$quant.msstats
                dataControl$data_type <- "TMT"
            }else{
                prePquant$DDA2009.proposed <- n.env$DDA2009.proposed
                prePquant$DDA2009.comparisons <- n.env$DDA2009.comparisons
                dataControl$data_type <- "LFQ"
            }
            
            progress$set(message = "Load over.", value = 1)
          
        }
    })

        
    output$contents <- renderDT({
        if(is.null(dataControl$inputData_state)){
            return(NULL)
        }
        else{
            fileData <- inputdf()
            datatable(fileData,options = list(scrollX = TRUE))
        }
    })
    
    
    ### enable, disable part
    observe({
        if(is.null(dataControl$inputData_state)){
            disable("start_preprocess")}
        else{enable("start_preprocess")}
    })
    
    observe({
        if(is.null(prePquant$DDA2009.comparisons) & is.null(prePquant$test.pairwise)){
            disable("download_Render")
        }else{
            enable("download_Render")
        }
    })
    
    # project name control
    observe({
        shinyjs::toggleState("SaveProjectDataName", input$SaveProjectDataControl == TRUE)
        shinyjs::toggleState("SaveProjectFolderPath", input$SaveProjectDataControl == TRUE)
    })
    
    observe({
      if(input$start_preprocess == 0 & is.null(dataControl$inputData_rdata)){
            disable("initialData_profilecondition_Render")
            disable("initialData_qualitycontrol_Render")
            disable("default_method_residual_qq_Render")
            disable("default_method_volcano_Render")
            disable("default_method_heatmap_Render")
            disable("default_method_comparison_Render")}
        else{
            enable("initialData_profilecondition_Render")
            enable("initialData_qualitycontrol_Render")
            enable("default_method_residual_qq_Render")
            enable("default_method_volcano_Render")
            enable("default_method_heatmap_Render")
            enable("default_method_comparison_Render")}
    })
    

    #dynamic volcano part
    observe({
        if(is.null(dataControl$inputData_state)) {
            disable("start_anno")}
        else{enable("start_anno")}
    })
    
    observe({
        if(dataControl$annoStart == 0) {
            disable("submit_anno")}
        else{enable("submit_anno")}
    })
    
    observe({
        if(is.null(dataControl$inputData_state) | dataControl$annoSubmit == 0){
            disable("dynamic_volcano_Render")
            disable("dynamic_FID_Render")}
            
        else{
            enable("dynamic_volcano_Render")
            enable("dynamic_FID_Render")}
    })
    
    
    ###start_preprocess 
    observeEvent(input$start_preprocess, {
      
        sendSweetAlert(
            session = session,
            title = "Judge species",
            text = "Start to judge the type of species, please wait...",
            type = "success", 
            closeOnClickOutside = TRUE,
            width = 400
        )
      
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        
        progress$set(message = "Begin to judge species, please wait...", value = 0.1)
      
        fileData <- inputdf()
        
        df <- fileData
        progress$set(message = "Begin to judge species, please wait...", value = 0.2)
        tmp_origin <- lapply(strsplit(as.character(df$ProteinName), "\\;"), "[")
        progress$set(message = "Begin to judge species, please wait....", value = 0.3)
        tmp_df <- as.data.frame(t(sapply(tmp_origin, "[", i = 1:max(sapply(tmp_origin, length)))))
        
        tmp <- as.data.frame(matrix(nrow=length(rownames(tmp_df)),ncol=length(colnames(tmp_df))))
        progress$set(message = "Begin to judge species, please wait...", value = 0.4)
        for(i in array(1:length(colnames(tmp_df)))){
          tmp[,i] <- unlist(lapply(strsplit(as.character(tmp_df[,i]), "\\|"), "[", 2))
        }
        
        progress$set(message = "Begin to judge species, please wait...", value = 0.6)
        s <- strsplit(as.character(tmp_df[,1]), "\\_")
        progress$set(message = "Begin to judge species, please wait...", value = 0.7)
        rs <- purrr::map2(s, lengths(s), pluck)
        progress$set(message = "Begin to judge species, please wait...", value = 0.8)
        rs_df <- t(as.data.frame(rs))
        species <- unique(unlist(rs))
        
        have_species <- c("HUMAN", "YEAST", "RAT", "ECOLI")
        sp_count <- 0
        for(sp in species){
            if(sum(grepl(sp, have_species))){
                sp_count <- sp_count + 1
            }
        }
        if(sp_count == length(species)){
            dataControl$db_type <- "yes"
        }
        
        progress$set(message = "Judging over.", value = 1)

        if(dataControl$db_type != "yes"){
            sendSweetAlert(
                session = session,
                title = "Warning",
                text = "Currently, only 'Human', 'Yeast', 'Rat' and 'Ecoli' species are supported。If the data contains other species, please contact developers to add.",
                type = "warning", 
                closeOnClickOutside = TRUE,
                width = 400
            )
        }else{
            if(sum("Mixture" == colnames(fileData)) == 1){
                dataControl$data_type <- "TMT"
            }else{
                dataControl$data_type <- "LFQ"
            }
          
            if(sum("Fraction" == colnames(fileData)) == 1 & sum("TechReplicate" == colnames(fileData)) == 0
               & dataControl$data_type == "LFQ")
            {
                sendSweetAlert(
                    session = session,
                    title = "Warning",
                    text = "For data with the [Fraction] column, please add the corresponding [TechReplicate] column then re-upload",
                    type = "warning", 
                    closeOnClickOutside = TRUE,
                    width = 400
                )
            }else{
                ### progress function
              
                sendSweetAlert(
                    session = session,
                    title = "Start",
                    text = "Your data is being preprocessed, please wait until done to the next step.",
                    type = "success", 
                    closeOnClickOutside = TRUE,
                    width = 400
                )
                
                progress <- shiny::Progress$new()
                on.exit(progress$close())
                
                progress$set(message = "Begin to preprocess data, please wait...", value = 0.1)
                
                if(input$user_choose_pre_pro2gene == TRUE){
                    for(sp in species){
                        if(sp == "HUMAN"){ db_type = org.Hs.eg.db
                        } else if(sp == "YEAST"){ db_type = org.Sc.sgd.db
                        } else if(sp == "RAT"){ db_type = org.Rn.eg.db 
                        } else if(sp == "ECOLI") { db_type = org.EcK12.eg.db 
                        }else{ db_type = "out" }
                        
                        df_s_proname <- as.data.frame(cbind(rs_df, df$ProteinName))
                        df_s_proname$V2[df_s_proname$V1 != sp] <- NA
                        
                        tmp_origin <- lapply(strsplit(as.character(df_s_proname$V2), "\\;"), "[")
                        tmp_df <- as.data.frame(t(sapply(tmp_origin, "[", i = 1:max(sapply(tmp_origin, length)))))
                        tmp_access <- as.data.frame(matrix(nrow=length(rownames(tmp_df)),ncol=length(colnames(tmp_df))))
                        for(i in array(1:length(colnames(tmp_df)))){
                            tmp_access[,i] <- unlist(lapply(strsplit(as.character(tmp_df[,i]), "\\|"), "[", 2))
                        }
                        
                        mapping_geneid <- data.frame(ProteinName = tmp_access[,1])
                        
                        for(i in array(1:length(colnames(tmp_df)))){
                            uniKeys <- as.character(tmp[,i])
                            
                            if('try-error' %in% class(try(AnnotationDbi::mapIds(db_type, keys=uniKeys, column="ENTREZID", keytype="UNIPROT")))){
                                tmp_geneid <- as.data.frame(matrix(nrow=length(rownames(tmp_df)),ncol=1))
                                colnames(tmp_geneid)[1] <- 'ENTREZID'
                                mapping_geneid <- cbind(mapping_geneid, tmp_geneid)
                            }else{
                                tmp_geneid <- AnnotationDbi::mapIds(db_type, keys=uniKeys, column="ENTREZID", keytype="UNIPROT")
                                tmp_geneid <- data.frame(matrix(lapply(tmp_geneid, as.character)))
                                tmp_geneid <- unlist(lapply(tmp_geneid[,1],function(x) if(identical(x,character(0))) NA else x))
                                tmp_geneid <- data.frame("ENTREZID"=tmp_geneid)
                                mapping_geneid <- cbind(mapping_geneid, tmp_geneid)
                            }
                        }
                    }
                    
                    mapping_geneid <- mapping_geneid[,-1]
                    cols_ENTREZID = "ENTREZID"
                    for(i in array(1:(length(colnames(mapping_geneid))-1))){
                        cols_ENTREZID <- cols_ENTREZID %>% append(paste("ENTREZID", ".", i, sep = ""))
                    }
                    
                    mapping_geneid <- tidyr::unite(mapping_geneid, "ENTREZID", cols_ENTREZID, sep=";", na.rm=TRUE)
                    #mapping_geneid <- mapping_geneid[,][nchar(mapping_geneid[,])>0]  # remove "" rows
                    
                    rows <- which(mapping_geneid[,1] == "", arr.ind = TRUE)
                    
                    mapping_geneid <- data.frame(mapping_geneid)
                    rs_df <- data.frame(rs_df)
                    
                    df_colnames <- array(1:length(row.names(rs_df)))
                    rownames(rs_df) <- df_colnames
                    rownames(mapping_geneid) <- df_colnames
                    
                    mapping_geneid <- mapping_geneid[!(row.names(mapping_geneid) %in% rows),]
                    rs_df <- rs_df[!(row.names(rs_df) %in% rows),]
                    df <- df[!(row.names(df) %in% rows),]
                    fileData <- fileData[!(row.names(fileData) %in% rows),]
                    
                    print(length(row.names(mapping_geneid)))
                    print(length(row.names(rs_df)))
                    print(length(row.names(df)))
                    print(length(row.names(fileData)))
                    
                    df$mapping_geneid <- mapping_geneid
                    df$rs_df <- rs_df
                    acc_sp <- tidyr::unite(df, "ProteinName", mapping_geneid, rs_df, sep="_")
                    fileData$ProteinName = acc_sp$ProteinName
                    
                } else {
                  
                    # Protein accession clean
                    df <- fileData
                    tmp_origin <- lapply(strsplit(as.character(df$ProteinName), "\\;"), "[")    # type: list
                    tmp_df <- as.data.frame(t(sapply(tmp_origin, "[", i = 1:max(sapply(tmp_origin, length)))))  # list to df
                    tmp_access <- as.data.frame(matrix(nrow=length(rownames(tmp_df)),ncol=length(colnames(tmp_df))))
                    for(i in array(1:length(colnames(tmp_df)))){
                        tmp_access[,i] <- unlist(lapply(strsplit(as.character(tmp_df[,i]), "\\|"), "[", 2))
                    }
                    
                    cols_proName = "ProteinName"
                    for(i in array(1:(length(colnames(tmp_access))-1))){
                        cols_proName <- cols_proName %>% append(paste("ProteinName", ".", i, sep = ""))
                    }
                    names(tmp_access) <- cols_proName
                    pro_accession <- tidyr::unite(tmp_access, "ProteinName", cols_proName, sep=";", na.rm=TRUE)
                    df$pro_accession <- pro_accession$ProteinName
                    df$rs_df <- rs_df
                    acc_sp <- tidyr::unite(df, "ProteinName", pro_accession, rs_df, sep="_")
                    fileData$ProteinName = acc_sp$ProteinName
                }
                
                progress$set(message = "Begin to preprocess data, please wait...", value = 0.2)
                # data type: LFQ
                if(dataControl$data_type == "LFQ"){
                    prePquant$DDA2009.proposed <- MSstats::dataProcess(raw = fileData,
                                                                       logTrans = as.numeric(input$user_choose_logTrans),
                                                                       normalization = input$user_choose_normalization,
                                                                       summaryMethod = input$user_choose_summaryMethod,
                                                                       maxQuantileforCensored = input$user_choose_maxQuantileforCensored,
                                                                       censoredInt = "NA",
                                                                       MBimpute = TRUE,
                                                                       use_log_file = FALSE)
                    
                    progress$set(message = "Begin to generate group comparison, please wait...", value = 0.5)
                    # Automatically create the manually created matrix in MSstats, user manual p23
                    len <- length(levels(prePquant$DDA2009.proposed$FeatureLevelData$GROUP))
                    
                    tmp <- t(combn(len,2))
                    matrix_len = length(t(combn(len,2))) / 2
                    
                    ourMatrix <- matrix(c(0:0),nrow=matrix_len,ncol=len)
                    
                    for(i in 1:matrix_len){
                        ourMatrix[i, tmp[i]] = -1
                        ourMatrix[i, tmp[i + matrix_len]] = 1
                    }
                    
                    ourCondition <- levels(prePquant$DDA2009.proposed$ProteinLevelData$GROUP)
                    tmp_name <- matrix(ourCondition, nr=len, nc=1)
                    name <- matrix(nr=matrix_len, nc=1)
                    for(i in 1:matrix_len){
                        name[i,1] <- sprintf('%s-%s', tmp_name[tmp[i+matrix_len]], tmp_name[tmp[i]])
                    }
                    row.names(ourMatrix) <- name
                    
                    #End of creation
                    colnames(ourMatrix) <- ourCondition
                    progress$set(message = "Begin to generate group comparison, please wait...", value = 0.6)
                    prePquant$DDA2009.comparisons <- MSstats::groupComparison(contrast.matrix = ourMatrix,
                                                                              data = prePquant$DDA2009.proposed,
                                                                              use_log_file = FALSE)
                }else{
                    # data type: TMT
                    input.om <- OpenMStoMSstatsTMTFormat(fileData, use_log_file = FALSE)
                    progress$set(message = "Begin to preprocess data, please wait...", value = 0.5)
                    if(is.na(input$user_choose_TMT_maxQuantileforCensored)){
                        prePquant$quant.msstats <- MSstatsTMT::proteinSummarization(input.om,
                                                                                    method = input$user_choose_TMT_method,
                                                                                    global_norm = as.logical(input$user_choose_TMT_global_norm),
                                                                                    reference_norm = TRUE,
                                                                                    remove_norm_channel = TRUE,
                                                                                    remove_empty_channel = TRUE,
                                                                                    maxQuantileforCensored = NULL,
                                                                                    use_log_file = FALSE) 
                    }else{
                        prePquant$quant.msstats <- MSstatsTMT::proteinSummarization(input.om,
                                                                                    method = input$user_choose_TMT_method,
                                                                                    global_norm = as.logical(input$user_choose_TMT_global_norm),
                                                                                    reference_norm = TRUE,
                                                                                    remove_norm_channel = TRUE,
                                                                                    remove_empty_channel = TRUE,
                                                                                    maxQuantileforCensored = input$user_choose_TMT_maxQuantileforCensored,
                                                                                    use_log_file = FALSE)}
                    
                      progress$set(message = "Begin to generate group comparison, please wait...", value = 0.6)
                      # test for all the possible pairs of conditions
                      prePquant$test.pairwise <- MSstatsTMT::groupComparisonTMT(prePquant$quant.msstats,
                                                                                moderated = TRUE,
                                                                                save_fitted_models = TRUE,
                                                                                use_log_file = FALSE)

                }
                
                progress$set(message = "Begin to convert ID, please wait...", value = 0.8)
                if(input$user_choose_pre_pro2gene == FALSE){
                    #map protein accsession to geneID + geneName

                    if(dataControl$data_type == "LFQ"){
                        df_protein <- as.character(prePquant$DDA2009.comparisons$ComparisonResult$Protein)
                    }else{
                        df_protein <- as.character(prePquant$test.pairwise$ComparisonResult$Protein)
                    }
                    
                    mapping_geneid <- data.frame(ProteinName = df_protein)
                    mapping_genename <- data.frame(ProteinName = df_protein)
                    
                    for(sp in species){
                        if(sp == "HUMAN"){ db_type = org.Hs.eg.db
                        } else if(sp == "YEAST"){ db_type = org.Sc.sgd.db
                        } else if(sp == "RAT"){ db_type = org.Rn.eg.db 
                        } else if(sp == "ECOLI") { db_type = org.EcK12.eg.db }

                        df_s <- unlist(lapply(strsplit(df_protein, "\\_"), "[", 2))
                        df_proname <- unlist(lapply(strsplit(df_protein, "\\_"), "[", 1))
                        df_s_proname <- data.frame(df_s, df_proname)
                        df_s_proname$df_proname[df_s_proname$df_s != sp] <- NA
                        
                        tmp_origin <- lapply(strsplit(as.character(df_s_proname$df_proname), "\\;"), "[")    # type: list
                        tmp_df <- as.data.frame(t(sapply(tmp_origin, "[", i = 1:max(sapply(tmp_origin, length)))))  # list to df
                        
                        tmp <- tmp_df
                        
                        for(i in array(1:length(colnames(tmp_df)))){
                            uniKeys <- as.character(tmp[,i])
                            
                            if('try-error' %in% class(try(AnnotationDbi::mapIds(db_type, keys=uniKeys, column="ENTREZID", keytype="UNIPROT")))){
                                tmp_geneid <- as.data.frame(matrix(nrow=length(rownames(tmp_df)),ncol=1))
                                colnames(tmp_geneid)[1] <- 'ENTREZID'
                                mapping_geneid <- cbind(mapping_geneid, tmp_geneid)
                            }else{
                                tmp_geneid <- AnnotationDbi::mapIds(db_type, keys=uniKeys, column="ENTREZID", keytype="UNIPROT")
                                tmp_geneid <- data.frame(matrix(lapply(tmp_geneid, as.character)))
                                tmp_geneid <- unlist(lapply(tmp_geneid[,1],function(x) if(identical(x,character(0))) NA else x))
                                tmp_geneid <- data.frame("ENTREZID"=tmp_geneid)
                                mapping_geneid <- cbind(mapping_geneid, tmp_geneid)
                            }
                            
                            if(species == "YEAST"){  #YEAST数据库没有SYMBOL列，对应的是GENENAME列
                                if('try-error' %in% class(try(AnnotationDbi::mapIds(db_type, keys=uniKeys, column="GENENAME", keytype="UNIPROT")))){
                                    tmp_genename <- as.data.frame(matrix(nrow=length(rownames(tmp_df)),ncol=1))
                                    colnames(tmp_genename)[1] <- 'GENENAME'
                                    mapping_genename <- cbind(mapping_genename, tmp_genename)
                                }else{
                                    tmp_genename <- AnnotationDbi::mapIds(db_type, keys=uniKeys, column="GENENAME", keytype="UNIPROT")
                                    tmp_genename <- data.frame(matrix(lapply(tmp_genename, as.character)))
                                    tmp_genename <- unlist(lapply(tmp_genename[,1],function(x) if(identical(x,character(0))) NA else x))
                                    tmp_genename <- data.frame("GENENAME"=tmp_genename)
                                    mapping_genename <- cbind(mapping_genename, tmp_genename)
                              }
                            } else {
                                if('try-error' %in% class(try(AnnotationDbi::mapIds(db_type, keys=uniKeys, column="SYMBOL", keytype="UNIPROT")))){
                                    tmp_genename <- as.data.frame(matrix(nrow=length(rownames(tmp_df)),ncol=1))
                                    colnames(tmp_genename)[1] <- 'GENENAME'
                                    mapping_genename <- cbind(mapping_genename, tmp_genename)
                                }else{
                                    tmp_genename <- AnnotationDbi::mapIds(db_type, keys=uniKeys, column="SYMBOL", keytype="UNIPROT")
                                    tmp_genename <- data.frame(matrix(lapply(tmp_genename, as.character)))
                                    tmp_genename <- unlist(lapply(tmp_genename[,1],function(x) if(identical(x,character(0))) NA else x))
                                    tmp_genename <- data.frame("GENENAME"=tmp_genename)
                                    mapping_genename <- cbind(mapping_genename, tmp_genename)
                                }
                            }
                        }
                    }
              
                    progress$set(message = "Begin to convert ID, please wait...", value = 0.9)
                    mapping_geneid <- mapping_geneid[,-1]
                    mapping_genename <- mapping_genename[,-1]
                    cols_ENTREZID = "ENTREZID"
                    cols_GENENAME = "GENENAME"
                    for(i in array(1:(length(colnames(mapping_geneid))-1))){
                        cols_ENTREZID <- cols_ENTREZID %>% append(paste("ENTREZID", ".", i, sep = ""))
                        cols_GENENAME <- cols_GENENAME %>% append(paste("GENENAME", ".", i, sep = ""))
                    }
            
                    mapping_geneid <- tidyr::unite(mapping_geneid, "ENTREZID", cols_ENTREZID, sep=";", na.rm=TRUE)
                    mapping_genename <- tidyr::unite(mapping_genename, "GENENAME", cols_GENENAME, sep=";", na.rm=TRUE)
                    
                    if(dataControl$data_type == "LFQ"){
                        prePquant$DDA2009.comparisons$ComparisonResult$ENTREZID <- mapping_geneid$ENTREZID
                        prePquant$DDA2009.comparisons$ComparisonResult$GENENAME <- mapping_genename$GENENAME
                    }else{
                        prePquant$test.pairwise$ComparisonResult$ENTREZID <- mapping_geneid$ENTREZID
                        prePquant$test.pairwise$ComparisonResult$GENENAME <- mapping_genename$GENENAME
                    }
                  
                } else {
                    # map geneID to geneName
                    if(dataControl$data_type == "LFQ"){
                        df_protein <- as.character(prePquant$DDA2009.comparisons$ComparisonResult$Protein)
                    }else{
                        df_protein <- as.character(prePquant$test.pairwise$ComparisonResult$Protein)
                    }
                    df_s <- unlist(lapply(strsplit(df_protein, "\\_"), "[", 2))
                    df_proname <- unlist(lapply(strsplit(df_protein, "\\_"), "[", 1))
                    df_s_proname <- data.frame(df_s, df_proname)
                    df_s_proname$df_proname[df_s_proname$df_s != sp] <- NA
                    
                    tmp_origin <- lapply(strsplit(as.character(df_s_proname$df_proname), "\\;"), "[")    # type: list
                    tmp_df <- as.data.frame(t(sapply(tmp_origin, "[", i = 1:max(sapply(tmp_origin, length)))))  # list to df

                    mapping_geneid <- tmp_df
                  
                    mapping_genename <- data.frame(ProteinName = mapping_geneid <- tmp_df[,1])
                    for(i in array(1:length(colnames(mapping_geneid <- tmp_df)))){
                        uniKeys <- as.character(mapping_geneid <- tmp_df[,i])
                        tmp_genename <- clusterProfiler::bitr(uniKeys, OrgDb = db_type, fromType = "ENTREZID", toType = "SYMBOL")
                        tmp_genename <- data.frame("GENENAME"=tmp_genename)
                        colnames(tmp_genename)[1] <- 'geneid'
                        
                        tmp <- data.frame(mapping_geneid <- tmp_df[,i])
                        colnames(tmp)[1] <- 'geneid'
                        tmp <- dplyr::left_join(tmp, tmp_genename, by = "geneid")
                        mapping_genename <- cbind(mapping_genename, tmp[2])
                    }
                    mapping_genename <- mapping_genename[,-1]
                    cols_GENENAME = "GENENAME.SYMBOL"
                    for(i in array(1:(length(colnames(mapping_geneid <- tmp_df))-1))){
                        cols_GENENAME <- cols_GENENAME %>% append(paste("GENENAME.SYMBOL", ".", i, sep = ""))
                    }
                    mapping_genename <- tidyr::unite(mapping_genename, "GENENAME", cols_GENENAME, sep=";", na.rm=TRUE)

                    if(dataControl$data_type == "LFQ"){
                        prePquant$DDA2009.comparisons$ComparisonResult$GENENAME <- mapping_genename$GENENAME
                    }else{
                        prePquant$test.pairwise$ComparisonResult$GENENAME <- mapping_genename$GENENAME
                    }
                }
                
                dataControl$preprocess_judge <- 1
                
                progress$set(message = "Preprocessing is over.", value = 1)
                
                if(input$SaveProjectDataControl == FALSE){
                    disable("start_preprocess")
                }
                
                
                observeEvent(input$start_preprocess, {
                    if((input$SaveProjectDataControl == TRUE) & dataControl$preprocess_judge == 1){
                        sendSweetAlert(
                            session = session,
                            title = "Start",
                            text = "Start to save data to local",
                            type = "success", 
                            closeOnClickOutside = TRUE,
                            width = 400
                        )
                        
                        progress <- shiny::Progress$new()
                        on.exit(progress$close())
                      
                        projectpathname <- paste(input$SaveProjectFolderPath, "/", input$SaveProjectDataName, sep = "")
                        dir.create(projectpathname)
                        
                        progress$set(message = "Begin to save data, please wait...", value = 0.5)
                        
                        inputData_rdata_out_msstats <- inputdf()
                        
                        if(dataControl$data_type == "LFQ"){
                            DDA2009.proposed <- isolate(prePquant$DDA2009.proposed)
                            DDA2009.comparisons <- isolate(prePquant$DDA2009.comparisons)
                            save(DDA2009.proposed, DDA2009.comparisons, inputData_rdata_out_msstats,
                                 file = paste(projectpathname,"/", "DDA2009.RData", sep = ""))
                        }else{
                            quant.msstats <- isolate(prePquant$quant.msstats)
                            test.pairwise <- isolate(prePquant$test.pairwise)
                            save(quant.msstats, test.pairwise, inputData_rdata_out_msstats,
                                 file = paste(projectpathname,"/", "DDA2009.RData", sep = ""))
                        }
                        
                        progress$set(message = "Save over.", value = 1)
                        
                        disable("start_preprocess")
                    }
                })
            }
        }
    })
   
     
    ## Reset uploaded data / reset all data
    observeEvent(input$uploaddata_reset, {
        enable("start_preprocess")
        
        reset("inputData")
        dataControl$inputData_URL_state <- NULL
        dataControl$inputData_state <- NULL
        dataControl$inputData_rdata <- NULL
        dataControl$db_type <- NULL
        dataControl$data_type <- "NULL"
        reset("download_name")
        reset("user_choose_pre_pro2gene")
        
        reset("SaveProjectDataControl")
        reset("SaveProjectDataName")
        reset("SaveProjectFolderPath")
        
        prePquant$DDA2009.proposed <- NULL
        prePquant$DDA2009.comparisons <- NULL
        prePquant$inputData_rdata_out_msstats <- NULL
        prePquant$inputData_URL_out_msstats <- NULL
        prePquant$quant.msstats <- NULL
        prePquant$test.pairwise <- NULL
        
        renderCheck$initialData_qualitycontrolPlot <- 0
        renderCheck$initialData_profilecondition <- 0
        renderCheck$default_method_residual_qq <- 0
        renderCheck$default_method_volcano <- 0
        renderCheck$default_method_heatmap <- 0
        renderCheck$default_method_comparison <- 0
        renderCheck$dynamic_volcano <- 0
        renderCheck$dynamic_FID <- 0
        
        dataControl$annoStart <- 0
        dataControl$annoSubmit <- 0
        dataControl$preprocess_judge <- 0
  
        volcanoLive$pdat <- NULL
        volcanoLive$res <- NULL
        volcanoLive$max_points <- NULL
        volcanoLive$pepdat <- NULL
        volcanoLive$prolfq_pepdat <- NULL
        volcanoLive$prolfq_prodat.med <- NULL
        volcanoLive$protmt_pepd <- NULL
        volcanoLive$protmt_prod <- NULL
        
        renderCheck$proteus_peptide <- 0
        renderCheck$proteus_protein <- 0
        renderCheck$proteus_pep_pro <- 0
        renderCheck$proteus_TMT_peptide <- 0
        renderCheck$proteus_TMT_protein <- 0
    })
    
    

    ##### ---initialData_tabbox_plot---
    initialData_profilecondition_select <- reactive({
        if(dataControl$data_type == "LFQ"){
            if(is.null(prePquant$DDA2009.proposed)) {
                selectInput(inputId = 'initialData_profilecondition_select_input',
                            label = 'Please upload data first',
                            choices = NULL)
            }
            else {
                initialData_profilecondition_selector <- levels(prePquant$DDA2009.proposed$ProteinLevelData$Protein)
                selectInput(inputId = 'initialData_profilecondition_select_input',
                            label = 'Options',
                            choices = as.list(initialData_profilecondition_selector))}
        }else{
            if(is.null(prePquant$quant.msstats)) {
                selectInput(inputId = 'initialData_profilecondition_select_input',
                            label = 'Please upload data first',
                            choices = NULL)
            }
            else {
                initialData_profilecondition_selector <- levels(prePquant$quant.msstats$ProteinLevelData$Protein)
                if(is.null(initialData_profilecondition_selector)){
                    initialData_profilecondition_selector <- unique(prePquant$quant.msstats$ProteinLevelData$Protein)
                }
                selectInput(inputId = 'initialData_profilecondition_select_input',
                            label = 'Options',
                            choices = as.list(initialData_profilecondition_selector))}}
    })
    
    initialData_qualitycontrol_select <- reactive({
        if(dataControl$data_type == "LFQ"){
            if(is.null(prePquant$DDA2009.proposed)) {
                selectInput(inputId = 'initialData_qualitycontrol_select_input',
                            label = 'Please upload data first',
                            choices = NULL)
            }
            else {
                tmp <- levels(prePquant$DDA2009.proposed$ProteinLevelData$Protein)
                initialData_qualitycontrol_selector <- append('allonly', tmp, 1)
                selectInput(inputId = 'initialData_qualitycontrol_select_input',
                            label = 'Options',
                            choices = as.list(initialData_qualitycontrol_selector))}
        }else{
            if(is.null(prePquant$quant.msstats)) {
                selectInput(inputId = 'initialData_qualitycontrol_select_input',
                            label = 'Please upload data first',
                            choices = NULL)
            }
            else {
                tmp <- levels(prePquant$quant.msstats$ProteinLevelData$Protein)
                if(is.null(tmp)){
                    tmp <- unique(prePquant$quant.msstats$ProteinLevelData$Protein)
                }
                initialData_qualitycontrol_selector <- append('allonly', tmp, 1)
                selectInput(inputId = 'initialData_qualitycontrol_select_input',
                            label = 'Options',
                            choices = as.list(initialData_qualitycontrol_selector))}
        }
        
    })
    
    output$initialData_profilecondition_selector <- renderUI({
        initialData_profilecondition_select()
    })
    
    output$initialData_qualitycontrol_selector <- renderUI({
        initialData_qualitycontrol_select()
    })
    
    
    #initialData_qualitycontrolPlot_out
    observeEvent(input$initialData_qualitycontrol_Render, {
        renderCheck$initialData_qualitycontrolPlot <- 1
    })
    
    output$initialData_qualitycontrolPlot_out <- renderPlot({
        if(renderCheck$initialData_qualitycontrolPlot > 0) {
            if(dataControl$data_type == "LFQ"){
                dataProcessPlots(data = prePquant$DDA2009.proposed, type="QCPlot",
                                 which.Protein=input$initialData_qualitycontrol_select_input,
                                 width=10, height=5, address=FALSE) 
            }else{
                dataProcessPlotsTMT(data = prePquant$quant.msstats, type='QCPlot',
                                    which.Protein=input$initialData_qualitycontrol_select_input,
                                    width=10, height=5, address=FALSE)}
        } else { return(NULL) }
    })
    
    #initialData_profilePlot_out
    observeEvent(input$initialData_profilecondition_Render, {
        renderCheck$initialData_profilecondition <- 1
    })
    
    output$initialData_profilePlot_out <- renderPlot({
        if(renderCheck$initialData_profilecondition > 0){
            if(dataControl$data_type == "LFQ"){
                dataProcessPlots(data = prePquant$DDA2009.proposed, type="ProfilePlot",
                                 which.Protein=input$initialData_profilecondition_select_input,
                                 width=10, height=5, address=FALSE)
            }else{
                dataProcessPlotsTMT(data = prePquant$quant.msstats, type = 'ProfilePlot',
                                    which.Protein=input$initialData_profilecondition_select_input,
                                    width=10, height=5, address=FALSE)
            }
        } else { return(NULL) }
        
    })
    
    #initialData_conditionPlot_out
    output$initialData_conditionPlot_out <- renderPlot({
        if(renderCheck$initialData_profilecondition > 0){
            if(dataControl$data_type == "LFQ"){
              dataProcessPlots(data = prePquant$DDA2009.proposed, type="ConditionPlot",
                               which.Protein=input$initialData_profilecondition_select_input,
                               width=10, height=5, address=FALSE)
            } else { return(NULL) }
        } else { return(NULL) }
        
    })
    
    
    ### preprocessed data  &  download
    output$preprocessedData <- renderDT({
        if(is.null(prePquant$DDA2009.comparisons) & is.null(prePquant$test.pairwise)){ return(NULL) }
        else{
            if(dataControl$data_type == "LFQ"){
                tmp <- prePquant$DDA2009.comparisons$ComparisonResult
                if(input$user_choose_pre_pro2gene == TRUE){ colnames(tmp)[1] <- 'Gene' }
                datatable(tmp,
                          options = list(scrollX = TRUE))
            }else{
                tmp <- prePquant$test.pairwise$ComparisonResult
                if(input$user_choose_pre_pro2gene == TRUE){ colnames(tmp)[1] <- 'Gene' }
                datatable(tmp,
                          options = list(scrollX = TRUE))
            }
        }
    })
    
    output$download_Render <- downloadHandler(
        filename = function() {
            paste(input$download_name, ".csv", sep = "")
        },
        content = function(file) {
            if(dataControl$data_type == "LFQ"){
                tmp <- prePquant$DDA2009.comparisons$ComparisonResult
            }else{
                tmp <- prePquant$test.pairwise$ComparisonResult
            }
            if(input$user_choose_pre_pro2gene == TRUE){ colnames(tmp)[1] <- 'Gene' }
            write.csv(tmp, file, row.names = FALSE)
        }
    )
    
    
    
    ##### -----default method-----  
    
    ### model-based QC plots
    default_method_residual_qq_select <- reactive({
        if(dataControl$data_type == "LFQ"){
            if(is.null(prePquant$DDA2009.proposed)) {
                selectInput(inputId = 'default_method_residual_qq_select_input',
                            label = 'Please upload data first',
                            choices = NULL)
            }
            else {
                default_method_residual_qq_selector <- levels(prePquant$DDA2009.proposed$ProteinLevelData$Protein)
                selectInput(inputId = 'default_method_residual_qq_select_input',
                            label = 'Options',
                            choices = as.list(default_method_residual_qq_selector))}
        }else{
            if(is.null(prePquant$quant.msstats)) {
                selectInput(inputId = 'default_method_residual_qq_select_input',
                            label = 'Please upload data first',
                            choices = NULL)
            }
            else {
                default_method_residual_qq_selector <- levels(prePquant$quant.msstats$ProteinLevelData$Protein)
                if(is.null(default_method_residual_qq_selector)){
                    default_method_residual_qq_selector <- unique(prePquant$quant.msstats$ProteinLevelData$Protein)
                }
                selectInput(inputId = 'default_method_residual_qq_select_input',
                            label = 'Options',
                            choices = as.list(default_method_residual_qq_selector))}
        }
        
    })
 
    output$default_method_residual_qq_selector <- renderUI({
        default_method_residual_qq_select()
    })
 
       
    #residual plot
    observeEvent(input$default_method_residual_qq_Render, {
        renderCheck$default_method_residual_qq <- 1
    })
    
    output$default_method_residualPlot_out <- renderPlot({
        if(renderCheck$default_method_residual_qq > 0){
            if(dataControl$data_type == "LFQ"){
                modelBasedQCPlots(data = prePquant$DDA2009.comparisons, type="ResidualPlots",
                                  which.Protein=input$default_method_residual_qq_select_input,
                                  width=10, height=5, address=FALSE)
            }else{
                modelBasedQCPlots(data = prePquant$test.pairwise, type="ResidualPlots",
                                  which.Protein=input$default_method_residual_qq_select_input,
                                  width=10, height=5, address=FALSE)}
        } else { return(NULL) }
    })
    
    #q-q plot (quantile-quantile)
    output$default_method_qqPlot_out <- renderPlot({
        if(renderCheck$default_method_residual_qq > 0){
            if(dataControl$data_type == "LFQ"){
                modelBasedQCPlots(data = prePquant$DDA2009.comparisons, type="QQPlots",
                                  which.Protein=input$default_method_residual_qq_select_input,
                                  width=10, height=5, address=FALSE)
            }else{
                modelBasedQCPlots(data = prePquant$test.pairwise, type="QQPlots",
                                  which.Protein=input$default_method_residual_qq_select_input,
                                  width=10, height=5, address=FALSE)}
        } else { return(NULL) }
    })
    
    
    ### group comparison plots
    # default_method volcano plot
    default_method_volcano_select <- reactive({
        if(is.null(prePquant$DDA2009.proposed) & is.null(prePquant$quant.msstats)) {
            selectInput(inputId = 'default_method_volcano_input',
                        label = 'Please upload data first',
                        choices = NULL)
        }
        else {
            if(dataControl$data_type == "LFQ"){
                default_method_vol_selector <- getSelector(inputdf(), flag = 'volcano',
                                                           prePquant$DDA2009.proposed)
            }else{
                default_method_vol_selector <- levels(factor(prePquant$test.pairwise$ComparisonResult$Label))
            }
            selectInput(inputId = 'default_method_volcano_input',
                        label = 'Options',
                        choices = as.list(default_method_vol_selector))
        }
    })
    
    output$default_method_volcano_selector <- renderUI({
        default_method_volcano_select()
    })
    
    observeEvent(input$default_method_volcano_Render, {
        renderCheck$default_method_volcano <- 1
    })
    
    output$default_method_volcano_out <- renderPlot({
        if(renderCheck$default_method_volcano > 0){
            if(dataControl$data_type == "LFQ"){
                groupComparisonPlots(data = prePquant$DDA2009.comparisons$ComparisonResult, type = 'VolcanoPlot',
                                     which.Comparison=input$default_method_volcano_input,
                                     width=5, height=5, address=FALSE)
            }else{
                groupComparisonPlots_TMT(data = prePquant$test.pairwise$ComparisonResult, type = 'VolcanoPlot',
                                         which.Comparison=input$default_method_volcano_input,
                                         width=5, height=5, address=FALSE)}
        } else { return(NULL) }
    })
    
    
    # default_method heatmap
    observeEvent(input$default_method_heatmap_Render, {
        renderCheck$default_method_heatmap <- 1
    })
    
    output$default_method_heat_out <- renderPlot({
        if(renderCheck$default_method_heatmap > 0){
            if(dataControl$data_type == "LFQ"){
                MS_output <- prePquant$DDA2009.comparisons$ComparisonResult
            }else{
                MS_output <- prePquant$test.pairwise$ComparisonResult
            }
          
            if(length(unique(na.omit(prePquant$test.pairwise$ComparisonResult$Label))) == 1 |
               length(unique(na.omit(prePquant$DDA2009.comparisons$ComparisonResult$Label))) == 1){
                return(NULL)
            }
          
            MS_output <- MS_output[,1:3]
            MS_output <- MS_output[is.finite(as.numeric(as.character(MS_output$log2FC))),]
            MS_output <- MS_output %>% filter(between(log2FC,-4,4)) %>% spread(key = Label, value = log2FC)
            MS_output <- IDPmisc::NaRV.omit(MS_output)
            
            heatmap <- MS_output[,-1]
            if(dataControl$data_type == "LFQ"){
                rownames(heatmap) <- MS_output$Protein
            }else{
                rownames(heatmap) <- levels(unlist(MS_output$Protein))
            }
            
            
            if(nrow(heatmap) > 50){
              pheatmap::pheatmap(heatmap, show_rownames = F)
            }else{ pheatmap::pheatmap(heatmap) }
            
            
        } else { return(NULL) }
    })
    
    
    # default_method comparison plot
    default_method_comparison_select <- reactive({
        if(is.null(prePquant$DDA2009.proposed) & is.null(prePquant$quant.msstats)) {
            selectInput(inputId = 'default_method_comparison_select_input',
                        label = 'Please upload data first',
                        choices = NULL)
        }
        else {
            if(dataControl$data_type == "LFQ"){
                default_method_comparison_selector <- levels(prePquant$DDA2009.proposed$ProteinLevelData$Protein)
            }else{
                default_method_comparison_selector <- levels(prePquant$quant.msstats$ProteinLevelData$Protein)
                if(is.null(default_method_comparison_selector)){
                    default_method_comparison_selector <- unique(prePquant$quant.msstats$ProteinLevelData$Protein)
                }
            }
            selectInput(inputId = 'default_method_comparison_select_input',
                        label = 'Options',
                        choices = as.list(default_method_comparison_selector))
        }
    })
  
    output$default_method_comparison_selector <- renderUI({
        default_method_comparison_select()
    })
    
    
    observeEvent(input$default_method_comparison_Render, {
        renderCheck$default_method_comparison <- 1
    })
    
    output$default_method_comparisonPlot_out <- renderPlot({
        if(renderCheck$default_method_comparison > 0){
            if(dataControl$data_type == "LFQ"){
                groupComparisonPlots(data = prePquant$DDA2009.comparisons$ComparisonResult, type="ComparisonPlot",
                                     which.Protein=input$default_method_comparison_select_input,
                                     width=10, height=5, address=FALSE)
            }else{
                groupComparisonPlots_TMT(data = prePquant$test.pairwise$ComparisonResult, type="ComparisonPlot",
                                         which.Protein=input$default_method_comparison_select_input,
                                         width=10, height=5, address=FALSE)}
        } else { return(NULL) }
    })
    
    
 
    
    ### ---dynamic volcano---
   
    # control data flow
    observeEvent(input$start_anno, {
       dataControl$annoStart <- 1
    })
    
    # Reset metadata settings
    observeEvent(input$annoReset, {
        dataControl$annoStart <- 0
        dataControl$annoSubmit <- 0
        
        renderCheck$dynamic_volcano <- 0
        renderCheck$dynamic_FID <- 0
        
        volcanoLive$pdat <- NULL
        volcanoLive$res <- NULL
        volcanoLive$max_points <- NULL
        volcanoLive$pepdat <- NULL
        volcanoLive$protmt_pepd <- NULL
        volcanoLive$protmt_prod <- NULL
    })
    
    ### start_convert_to_proteus
    
    output$dynamic_evidence_contents <- renderDT({
        if(is.null(dataControl$inputData_state)) {
            return(NULL)
        }
        else {
            fileData <- inputdf()
            dynamic_data <- unique(fileData[,grepl("Reference|Condition", colnames(fileData))])
            dynamic_data <- dynamic_data[order(dynamic_data$Condition),]

            datatable(dynamic_data,
                      rownames = FALSE,
                      options = list(scrollX = TRUE))
        }
    })

    
    ### annotate metadata for user
    categorial_anno <- reactive({
        fileData <- inputdf()
        if(dataControl$data_type == "LFQ"){
            dynamic_data <- unique(fileData[,grepl("Reference|Condition", colnames(fileData))])
            dynamic_data <- dynamic_data[order(dynamic_data$Condition),]
            colnames(dynamic_data)[1] <- 'condition'
            colnames(dynamic_data)[2] <- 'sample'
            metadata <- dynamic_data
            rownames(metadata) <- c(1:nrow(metadata))
        }else{
            meta_tmp <- unique(fileData[,grepl("Channel|Condition|TechRepMixture|Run", colnames(fileData))])
            meta_tmp <- tidyr::unite(meta_tmp, "sample", Run, Channel, TechRepMixture, sep="-", remove=FALSE)
            meta <- as.data.frame(matrix(nrow=length(rownames(meta_tmp)),ncol=2))
            meta_colnames <- c("condition", "sample")
            names(meta) <- meta_colnames
            meta$sample <- meta_tmp$sample
            meta$condition <- meta_tmp$Condition
            metadata <- meta
        }
        return(metadata)
    })
    
    
    output$dynamic_define_metadata <- renderRHandsontable({
        if (is.null(input$inputData)) { return(NULL) }
        else if (dataControl$annoStart > 0 & dataControl$inputData_state == "uploaded") {
            categorial_anno <- categorial_anno()
            categorial_anno$sample <- as.character(categorial_anno$sample)
            categorial_anno$condition <- as.character(categorial_anno$condition)
            tab <- rhandsontable(categorial_anno)
            return(tab)
        } else { return(NULL) }
    })
    
    
    #get data from user
    #displays need to signal data is submitted
    dynamic_metadata <- reactive({
        if (dataControl$annoSubmit > 0) {
            reps <- isolate(input$dynamic_define_metadata)
            repsOut <- hot_to_r(reps)
            return(repsOut)
        } else { return(NULL) }
    })
    
    
    #This controls enabling and disabling anno button
    observe({
        if (dataControl$annoSubmit > 0) {
            sendSweetAlert(
                session = session,
                title = "Success",
                text = "Your annotation is submitted",
                type = "success", 
                closeOnClickOutside = TRUE,
                width = 400
            )
            disable("submit_anno")
        }
    })
    
    # Proteus: create a peptide & protein dataset
    observeEvent(input$submit_anno, {
      if(input$start_preprocess == 0 & is.null(dataControl$inputData_rdata)){
          sendSweetAlert(
              session = session,
              title = "Fail",
              text = "Please preprocess '.csv(.gz)' file or  upload '.RData' file first",
              type = "error", 
              closeOnClickOutside = TRUE,
              width = 500)
      } else{
          dataControl$annoSubmit <- 1
          
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          
          progress$set(message = "Begin to preprocess data, please wait...", value = 0.1)
          fileData <- inputdf()
          
          if(dataControl$data_type == "LFQ"){
              meta <- dynamic_metadata()
              
              sequence.col <- "sequence"
              protein.col <- "protein"
              experiment.type <- "label-free"
              ncores = 4
              # mclapply doesn't work on Windows, so force 1 core
              if(Sys.info()[['sysname']] == "Windows") {
                  ncores <- 1
                  warning("Multicore processing not available in Windows. Using ncores=1")
              }
              
              progress$set(message = "Begin to preprocess data, please wait...", value = 0.3)
              measureColumns <- list(
                  Intensity = 'Intensity'
              )
              measures <- names(measureColumns)
              
              tabMelt <- fileData[c("PeptideSequence","Intensity","Reference")]
              names(tabMelt)[1] <- "sequence"
              names(tabMelt)[2] <- "value"
              names(tabMelt)[3] <- "sample"
              
              # create unique sequence names
              tabMelt$seqsam <- paste0(tabMelt$sequence, ".", tabMelt$sample)
              tabMelt <- tabMelt[order(tabMelt$seqsam),]
              s2s <- setNames(tabMelt$sequence, tabMelt$seqsam)
              r <- rle(tabMelt$seqsam)
              tabMelt$uniseq <- paste0(rep(s2s[r$values], times=r$lengths), ".", unlist(lapply(r$lengths, seq_len)))
              
              # cast into table: sample vs unique sequence
              # in this table there are multiple entries per peptide
              tab <- reshape2::dcast(tabMelt, uniseq ~ sample, sum, value.var="value")
              
              # original sequence
              u2s <- setNames(tabMelt$sequence, tabMelt$uniseq)
              sequences <- u2s[tab$uniseq]
              
              # extract numeric data as matrix
              tab <- as.matrix(tab[,2:ncol(tab)])
              tab[tab == 0] <- NA
              
              progress$set(message = "Begin to preprocess data, please wait...", value = 0.5)
              aggregateSum <- function(wp) {
                  row <- colSums(wp, na.rm=TRUE)
                  row[row==0] <- NA    # colSums puts zeroes where the column contains only NAs (!!!)
                  return(as.vector(row))
              }
              aggregate.fun=aggregateSum
              # aggregate peptides
              ptab <- parallel::mclapply(unique(sequences), function(s) {
                  wp <- tab[sequences == s,, drop=FALSE]
                  x <- aggregate.fun(wp)
                  row <- as.data.frame(t(as.vector(x)))
                  rownames(row) <- s
                  return(row)
              }, mc.cores=ncores)
              ptab <- as.matrix(do.call(rbind, ptab))
              colnames(ptab) <- colnames(tab)
              peptides <- row.names(ptab)
              
              # peptide to protein conversion
              pep2prot <- data.frame(sequence=fileData$PeptideSequence, protein=fileData$ProteinName)
              pep2prot <- unique(pep2prot)
              if(anyDuplicated(pep2prot$sequence) > 0) stop("Non-unique peptide-to-protein association found. Proteus requires that peptide sequence is uniquely associated with a protein or protein group.")
              rownames(pep2prot) <- pep2prot$sequence
              pep2prot <- pep2prot[peptides,]
              proteins <- levels(as.factor(pep2prot$protein))
              
              volcanoLive$pepdat <- ptab
              
              progress$set(message = "Begin to preprocess data, please wait...", value = 0.7)
              ### makeProteinTable
              tab <- ptab
              
              protlist <- list()
              for(cond in levels(factor(as.character(meta$condition)))) {
                  w <- tab[,which(meta$condition == cond), drop=FALSE]
                  samples <- colnames(w)
                  protcond <- parallel::mclapply(proteins, function(prot) {
                      sel <- which(pep2prot$protein == prot)
                      npep <- length(sel)
                      min.peptides=2
                      if(npep >= min.peptides) {
                          wp <- w[sel,, drop=FALSE]
                          x <- aggregate.fun(wp)
                          row <- as.data.frame(t(as.vector(x)))
                      } else {
                          row <- as.data.frame(t(rep(NA, length(samples))))
                      }
                      names(row) <- samples
                      row <- data.frame(protein=prot, row, check.names=FALSE)
                      return(row)
                  }, mc.cores=ncores)
                  protint <- do.call(rbind, protcond)
                  protlist[[cond]] <- protint
              }
            
              # dplyr join all tables
              protab <- Reduce(function(df1, df2) dplyr::full_join(df1,df2, by="protein"), protlist)
              
              # remove empty rows (happens when min.peptides > 1)
              protab <- protab[which(rowSums(!is.na(protab)) > 0), ]
              proteins <- protab$protein
              #protab <- as.matrix(protab[,as.character(unique(tabMelt$sample))])
              rownames(protab) <- protab[,1]
              protab <- protab[,-1]
              
              
              progress$set(message = "Begin to preprocess data, please wait...", value = 0.9)
              normalizeMedian <- function(tab) {
                  norm.fac <- apply(tab, 2, function(x) {median(x, na.rm=TRUE)})
                  norm.fac <- norm.fac / mean(norm.fac, na.rm=TRUE)
                  tab <- t(t(tab) / norm.fac)
                  return(tab)
              }
              protab <- normalizeMedian(protab)
              
              volcanoLive$pdat <- protab
          
          }else{
            
              df <- as.data.frame(matrix(nrow=length(rownames(fileData)),ncol=9))
              proteus_colnames <- c("sequence", "modified_sequence", "protein_group", "proteins", "protein",
                                    "experiment", "charge", "reverse", "contaminant")
              names(df) <- proteus_colnames
              
              df$sequence <- fileData$PeptideSequence
              df$proteins <- fileData$ProteinName
              df_protein <- lapply(strsplit(fileData$ProteinName, "\\;"), "[", 1)
              df$protein <- sapply(df_protein, "[", i = 1:max(sapply(df_protein, length)))
              df$charge <- fileData$Charge
              df$experiment <- fileData$Run
              
              df$remove1 <- fileData$Reference
              df$remove2 <- fileData$RetentionTime
              
              df$Intensity <- fileData$Intensity
              df$Channel <- fileData$Channel
              df <- df %>% spread(key = Channel, value = Intensity, fill = NA)
              df <- subset(df, select = -c(remove1, remove2))
              
              meta_tmp <- unique(fileData[,grepl("Channel|Condition|TechRepMixture|Run", colnames(fileData))])
              meta_tmp <- tidyr::unite(meta_tmp, "sample", Run, Channel, TechRepMixture, sep="-", remove=FALSE)
              
              meta <- as.data.frame(matrix(nrow=length(rownames(meta_tmp)),ncol=5))
              meta_colnames <- c("experiment", "measure", "sample", "condition", "replicate")
              names(meta) <- meta_colnames
              
              meta$experiment <- meta_tmp$Run
              meta$measure <- meta_tmp$Channel
              meta$sample <- meta_tmp$sample
              meta$condition <- meta_tmp$Condition
              meta$replicate <- meta_tmp$TechRepMixture
              
              measCols <- unique(meta$measure)
              names(measCols) <- unique(meta$measure)
              pepdat <- proteus::makePeptideTable(df, meta, measure.cols=measCols, aggregate.fun=aggregateMedian, experiment.type="TMT")
              prodat <- proteus::makeProteinTable(pepdat, aggregate.fun=aggregateHifly, hifly=3)
              prodat.norm <- proteus::normalizeTMT(prodat)
              
              volcanoLive$pdat <- prodat.norm$tab
              
              volcanoLive$protmt_pepd <- pepdat
              volcanoLive$protmt_prod <- prodat.norm
          }
          
          
          volcanoLive$max_points = 100
          
          progress$set(message = "Preprocessing is over.", value = 1)
        }
    })
    
    
    
    ### dynamic volcano Plot
    
    dynamic_volcano_select <- reactive({
        if(is.null(input$inputData)) {
            selectInput(inputId = 'dynamic_volcano_input',
                        label = 'Please upload data first',
                        choices = NULL)
        }
        else if(dataControl$annoSubmit == 0){
            selectInput(inputId = 'dynamic_volcano_input',
                        label = 'Please submit annotation first',
                        choices = NULL)
        }
        else {
            if(dataControl$data_type == "LFQ"){
                dynamic_vol_selector <- getSelector(inputdf(), flag = 'volcano',
                                                    prePquant$DDA2009.proposed)
            }else{
                dynamic_vol_selector <- levels(factor(prePquant$test.pairwise$ComparisonResult$Label))
            }
            selectInput(inputId = 'dynamic_volcano_input',
                        label = 'Options',
                        choices = as.list(dynamic_vol_selector))
        }
    })
    
    output$dynamic_volcano_selector <- renderUI({
        dynamic_volcano_select()
    })
    
    
    observeEvent(input$dynamic_volcano_Render, {
        renderCheck$dynamic_volcano <- 1
        
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        
        progress$set(message = "Begin to preprocess data, please wait...", value = 0.6)
        
        selector = input$dynamic_volcano_input
        if(dataControl$data_type == "LFQ"){
            comparisons <- subset(prePquant$DDA2009.comparisons$ComparisonResult, Label==selector)
        }else{
            comparisons <- subset(prePquant$test.pairwise$ComparisonResult, Label==selector)
        }
        volcanoLive$res <- comparisons
        
        volcanoLive$res$"-log10(pvalue)" <- -log10(volcanoLive$res$pvalue)
        
        rownames(volcanoLive$res) <- c(1:nrow(volcanoLive$res))
        
        progress$set(message = "Preprocessing is over.", value = 1)
        
        if(input$user_choose_pre_pro2gene == FALSE){
            dynamic_signif_flag = "n"
        } else {
            dynamic_signif_flag = "y"
        }

        output$dynamic_gap_out <- renderUI({HTML('<br/>')})
        
        output$dynamic_replicateTable_out <- dynamic_replicateTable(volcanoLive$res, input, volcanoLive$pdat, volcanoLive$max_points)
        output$dynamic_significanceTable_out <- dynamic_significanceTable(volcanoLive$res, volcanoLive$res, input, dynamic_signif_flag)
        output$dynamic_jitterPlot_out <- dynamic_jitterPlot(volcanoLive$res, input, volcanoLive$pdat, volcanoLive$max_points, dynamic_metadata())
        
        
        ## Volcano plot
        output$dynamic_plotVolcano_out <- renderPlot({
            if (renderCheck$dynamic_volcano > 0) {
                tab_idx <- as.numeric(input$allProteinTable_rows_selected)
                pVol <- dynamic_plotVolcano(volcanoLive$res, binhex=FALSE)
                if(length(tab_idx) > 0) {
                    pVol <- pVol + geom_point(data=volcanoLive$res[tab_idx,], size=3, color='red')
                }
                return(pVol)
            } else { return(NULL) }
        })
        
        
        output$dynamic_allProteinTable_out <- DT::renderDataTable({
            if (renderCheck$dynamic_volcano > 0) {
                # assume first column is id ("protein" or "peptide")
                idcol <- names(volcanoLive$res)[1]
                cols <- c(idcol, "log2FC", "pvalue", "adj.pvalue")
                if(dataControl$data_type == "LFQ"){
                    d <- volcanoLive$res[, cols]
                    d[, 2:ncol(d)] <- sapply(d[, 2:ncol(d)], function(x) signif(x, 3))
                }
                else{
                    d <- subset(volcanoLive$res, select = cols) 
                }
                d <- DT::datatable(d, class = 'cell-border strip hover')
                DT::formatStyle(d, 0, cursor = 'pointer')
            } else { return(NULL) }
        })
        
    })
    
    
    
    ### dynamic FID Plot
    
    dynamic_FID_select <- reactive({
        if(is.null(input$inputData)) {
            selectInput(inputId = 'dynamic_FID_input',
                        label = 'Please upload data first',
                        choices = NULL)
        }
        else if(dataControl$annoSubmit == 0){
            selectInput(inputId = 'dynamic_FID_input',
                        label = 'Please submit annotation first',
                        choices = NULL)
        }
        else {
            if(dataControl$data_type == "LFQ"){
                dynamic_fid_selector <- getSelector(inputdf(), flag = 'volcano',
                                                    prePquant$DDA2009.proposed)
            }else{
                dynamic_fid_selector <- levels(factor(prePquant$test.pairwise$ComparisonResult$Label))
            }
            selectInput(inputId = 'dynamic_FID_input',
                        label = 'Options',
                        choices = as.list(dynamic_fid_selector))
        }
    })
    
    output$dynamic_FID_selector <- renderUI({
        dynamic_FID_select()
    })
    
    
    observeEvent(input$dynamic_FID_Render, {
        renderCheck$dynamic_FID <- 1
        
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        
        progress$set(message = "Begin to preprocess data, please wait...", value = 0.4)
        
        selector = input$dynamic_FID_input
        if(dataControl$data_type == "LFQ"){
            comparisons <- subset(prePquant$DDA2009.comparisons$ComparisonResult, Label==selector)
        }else{
            comparisons <- subset(prePquant$test.pairwise$ComparisonResult, Label==selector)
        }
        volcanoLive$res_FID <- comparisons
        
        volcanoLive$res_FID$"-log10(pvalue)" <- -log10(volcanoLive$res_FID$pvalue)
        
        rownames(volcanoLive$res_FID) <- c(1:nrow(volcanoLive$res_FID))
        
        progress$set(message = "Begin to preprocess data, please wait...", value = 0.5)
        
        tmp <- dynamic_metadata()
        condition <- unique(tmp$condition)
        if(dataControl$data_type == "LFQ"){
            FID_pair <- strsplit(selector, "-")[[1]]
        }else{
            FID_pair <- strsplit(selector, "vs")[[1]]
        }
        
        
        # Generate the fold-change/intensity dataset. The same as in the FID plot.
        condMeans <- function(cond) {
            m <- rowMeans(log10(volcanoLive$pdat)[,which(condition == cond), drop=FALSE], na.rm=TRUE)
            m[is.nan(m)] <- NA
            m
        }
        m1 <- condMeans(FID_pair[1])
        m2 <- condMeans(FID_pair[2])
        good <- !is.na(m1) & !is.na(m2)
        fi <- data.frame(
            id = rownames(volcanoLive$pdat),
            x = (m1 + m2) / 2,
            y = m2 - m1,
            good = good
        )
        rownames(fi) <- 1:nrow(fi)
        
        mx <- 1.1 * max(abs(fi$y), na.rm=TRUE)
        m <- rbind(m1[!good], m2[!good])
        fi[!good, "x"] <- colSums(m, na.rm=TRUE)
        fi[!good, "y"] <- ifelse(is.na(m[1,]), mx, -mx)
        
        progress$set(message = "Preprocessing is over.", value = 1)
        
        if(input$user_choose_pre_pro2gene == FALSE){
          dynamic_signif_flag = "n"
        } else {
          dynamic_signif_flag = "y"
        }
        
        output$dynamic_FID_gap_out <- renderUI({HTML('<br/>')})
        
        output$dynamic_FID_replicateTable_out <- dynamic_fid_replicateTable(fi, input, volcanoLive$pdat, volcanoLive$max_points)
        output$dynamic_FID_significanceTable_out <- dynamic_fid_significanceTable(fi, volcanoLive$res_FID, input, dynamic_signif_flag)
        output$dynamic_FID_jitterPlot_out <- dynamic_fid_jitterPlot(fi, input, volcanoLive$pdat, volcanoLive$max_points, dynamic_metadata())
        
        
        ## FID plot
        output$dynamic_plotFID_out <- renderPlot({
            if (renderCheck$dynamic_FID > 0) {
                tab_idx <- as.numeric(input$allProteinTable_rows_selected)
                pFID <- dynamic_fid_plotFID(volcanoLive$pdat, condition, FID_pair, binhex=FALSE)
                if(length(tab_idx) > 0) {
                    pFID <- pFID + geom_point(data=fi[tab_idx,], size=3, color='red')
                }
                return(pFID)
            } else { return(NULL) }
        })
        
        
        output$dynamic_FID_allProteinTable_out <- DT::renderDataTable({
            if (renderCheck$dynamic_FID > 0) {
                # assume first column is id ("protein" or "peptide")
                idcol <- names(volcanoLive$res_FID)[1]
                cols <- c(idcol, "log2FC", "pvalue", "adj.pvalue")
                if(dataControl$data_type == "LFQ"){
                    d <- volcanoLive$res_FID[, cols]
                    d[, 2:ncol(d)] <- sapply(d[, 2:ncol(d)], function(x) signif(x, 3))
                }
                else{
                    d <- subset(volcanoLive$res_FID, select = cols) 
                }
                d <- DT::datatable(d, class = 'cell-border strip hover')
                DT::formatStyle(d, 0, cursor = 'pointer')
            } else { return(NULL) }
        })
      
    })
    
    
    
    ##### -----Proteus method-----
    
    ## annotation check control
    
    observe({
        shinyjs::toggleState("proteus_peptide_Render", dataControl$annoSubmit > 0 & dataControl$data_type == "LFQ")
        shinyjs::toggleState("proteus_protein_Render", dataControl$annoSubmit > 0 & dataControl$data_type == "LFQ")
        shinyjs::toggleState("proteus_TMT_peptide_Render", dataControl$annoSubmit > 0 & dataControl$data_type == "TMT")
        shinyjs::toggleState("proteus_TMT_protein_Render", dataControl$annoSubmit > 0 & dataControl$data_type == "TMT")
        #shinyjs::toggleState("proteus_TMT_protein_Render", input$submit_anno == TRUE)
    })
    
    ## proteus: peptide data
    
    observe({
        if(renderCheck$proteus_peptide > 0 | renderCheck$proteus_protein > 0){
            if(renderCheck$proteus_pep_pro == 0 & dataControl$data_type == "LFQ"){
                progress <- shiny::Progress$new()
                on.exit(progress$close())
                
                progress$set(message = "Begin to preprocess data, please wait...", value = 0.2)
                
                fileData <- inputdf()
                proteus_colnames <- c("PeptideSequence", "modified_sequence", "modifications", "protein_group", "protein",
                                      "experiment", "charge", "intensity", "sequence", "accession")
                
                df <- as.data.frame(matrix(nrow=length(rownames(fileData)),ncol=10))
                names(df) <- proteus_colnames
                
                df$protein_group <- fileData$ProteinName
                df$charge <- fileData$PrecursorCharge
                df$intensity <- fileData$Intensity
                
                df$experiment <- fileData$Reference
                
                df_protein <- lapply(strsplit(fileData$ProteinName, "\\;"), "[", 1)
                df_protein <- sapply(df_protein, "[", i = 1:max(sapply(df_protein, length)))
                df$protein <- df_protein
        
                df$sequence <- fileData$PeptideSequence
                
                progress$set(message = "Begin to preprocess data, please wait...", value = 0.4)
                
                
                dynamic_data <- unique(fileData[,grepl("Reference|Condition", colnames(fileData))])
                dynamic_data <- dynamic_data[order(dynamic_data$Condition),]
                colnames(dynamic_data)[1] <- 'condition'
                colnames(dynamic_data)[2] <- 'sample'
                metadata <- dynamic_data
                rownames(metadata) <- c(1:nrow(metadata))
                
                measure <- rep("Intensity", length(rownames(metadata)))
                metadata <- cbind(measure, metadata)
                
                experiment <- metadata$sample
                metadata <- cbind(experiment, metadata)
                
                progress$set(message = "Begin to preprocess data, please wait...", value = 0.6)
                volcanoLive$prolfq_pepdat <- proteus::makePeptideTable(df, metadata)
                progress$set(message = "Begin to preprocess data, please wait...", value = 0.8)
                prodat <- proteus::makeProteinTable(volcanoLive$prolfq_pepdat)
                progress$set(message = "Begin to preprocess data, please wait...", value = 0.9)
                volcanoLive$prolfq_prodat.med <- proteus::normalizeData(prodat)
                
                progress$set(message = "Preprocessing is over.", value = 1)
                
                renderCheck$proteus_pep_pro <- 1
            }
        }
        
        
        observeEvent(input$proteus_peptide_Render, {
            if(dataControl$data_type == "LFQ"){
                renderCheck$proteus_peptide <- 1
            } else { return(NULL) }
        })
        
        output$proteus_peptide_summary_out <- renderPrint({
            if(renderCheck$proteus_peptide > 0){
                summary(volcanoLive$prolfq_pepdat)
            } else { return(NULL) }
        })
        
        output$proteus_peptide_distance_out <- renderPlot({
            if(renderCheck$proteus_peptide > 0){
                if('try-error' %in% class(try(proteus::plotDistanceMatrix(volcanoLive$prolfq_pepdat)))){
                    return(NULL)
                }else{
                    proteus::plotDistanceMatrix(volcanoLive$prolfq_pepdat)
                }
            } else { return(NULL) }
        })
        
        output$proteus_peptide_number_out <- renderPlot({
            if(renderCheck$proteus_peptide > 0){
                #if('try-error' %in% class(try(proteus::plotCount(volcanoLive$prolfq_pepdat)))){
                #    return(NULL)
                #}else{
                    proteus::plotCount(volcanoLive$prolfq_pepdat)
                #}
            } else { return(NULL) }
        })
        
        output$proteus_peptide_jaccard_out <- renderPlot({
            if(renderCheck$proteus_peptide > 0){
                #if('try-error' %in% class(try(proteus::plotDetectionSimilarity(volcanoLive$prolfq_pepdat, bin.size = 0.02)))){
                #    return(NULL)
                #}else{
                    proteus::plotDetectionSimilarity(volcanoLive$prolfq_pepdat, bin.size = 0.02)
                #}
            } else { return(NULL) }
        })
        
        output$proteus_peptide_clustering_out <- renderPlot({
            if(renderCheck$proteus_peptide > 0){
                if('try-error' %in% class(try(proteus::plotClustering(volcanoLive$prolfq_pepdat)))){
                    return(NULL)
                }else{
                    proteus::plotClustering(volcanoLive$prolfq_pepdat)
                }
            } else { return(NULL) }
        })
        
        output$proteus_peptide_pca_out <- renderPlot({
            if(renderCheck$proteus_peptide > 0){
                if('try-error' %in% class(try(plotPCA_pquantr(volcanoLive$pepdat, dynamic_metadata())))){
                    return(NULL)
                }else{
                    plotPCA_pquantr(volcanoLive$pepdat, dynamic_metadata())
                }
            } else { return(NULL) }
        })
        
        
        ## proteus: protein data
        observeEvent(input$proteus_protein_Render, {
            if(dataControl$data_type == "LFQ"){
                renderCheck$proteus_protein <- 1
            } 
        })
        
        output$proteus_protein_summary_out <- renderPrint({
            if(renderCheck$proteus_protein > 0){
                summary(volcanoLive$prolfq_prodat.med)
            } else { return(NULL) }
        })
        
        output$proteus_protein_normalization_median_out <- renderPlot({
            if(renderCheck$proteus_protein > 0){
                proteus::plotSampleDistributions(volcanoLive$prolfq_prodat.med, title="Median normalization", fill="condition", method="violin")
            } else { return(NULL) }
        })
        
        output$proteus_protein_mean_out <- renderPlot({
            if(renderCheck$proteus_protein > 0){
                proteus::plotMV(volcanoLive$prolfq_prodat.med, with.loess=TRUE)
            } else { return(NULL) }
        })
        
        output$proteus_protein_clustering_out <- renderPlot({
            if(renderCheck$proteus_protein > 0){
                proteus::plotClustering(volcanoLive$prolfq_prodat.med)
            } else { return(NULL) }
        })
        
    })
    
    
    ## proteus TMT: peptide data
    observeEvent(input$proteus_TMT_peptide_Render, {
        if(dataControl$data_type == "TMT"){
            renderCheck$proteus_TMT_peptide <- 1
        }
    })
    
    output$proteus_TMT_peptide_summary_out <- renderPrint({
        if(renderCheck$proteus_TMT_peptide > 0){
            summary(volcanoLive$protmt_pepd)
        } else { return(NULL) }
    })
    
    output$proteus_TMT_peptide_number_out <- renderPlot({
        if(renderCheck$proteus_TMT_peptide > 0){
            proteus::plotCount(volcanoLive$protmt_pepd)
        } else { return(NULL) }
    })
    
    ## proteus TMT: protein data
    observeEvent(input$proteus_TMT_protein_Render, {
        if(dataControl$data_type == "TMT"){
            renderCheck$proteus_TMT_protein <- 1
        }
    })
    
    output$proteus_TMT_protein_summary_out <- renderPrint({
        if(renderCheck$proteus_TMT_protein > 0){
            summary(volcanoLive$protmt_prod)
        } else { return(NULL) }
    })
    
    output$proteus_TMT_protein_normalization_median_out <- renderPlot({
        if(renderCheck$proteus_TMT_protein > 0){
            proteus::plotSampleDistributions(volcanoLive$protmt_prod, log.scale=FALSE, fill="replicate") + labs(title="After")
        } else { return(NULL) }
    })
    
    output$proteus_TMT_protein_clustering_out <- renderPlot({
        if(renderCheck$proteus_TMT_protein > 0){
            proteus::plotClustering(volcanoLive$protmt_prod)
        } else { return(NULL) }
    })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
