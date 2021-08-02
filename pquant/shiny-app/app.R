library(shiny)
library(shinydashboard)
library(DT)
library('MSstats', warn.conflicts = F, quietly = T, verbose = F)
library(reticulate)
library(pheatmap)
library(rhandsontable)
library(MSstats)

setwd(getwd())

# Source getPlots functions
source("./app/getPlots.R")
source("./app/getSelector.R")

  
ui <- dashboardPage(
    dashboardHeader(title = "pQuantR"),
    dashboardSidebar(
        #useShinyjs(),
        
        # fileinput side bar menu
        conditionalPanel(condition = "input.main_tabs == 'fileinput_condition'",
                         sidebarMenu(
                             menuItem("Upload file", startExpanded = TRUE,
                                 # file selection box
                                 fileInput('csvFile', 'Choose the \'out_msstats.csv\'', multiple = FALSE, 
                                           accept=c('text/csv', 'text/comma-separated-values,text/plain')) # CSV text file
                                 #helpText(' Note: \'out_msstats.csv\' is ... ')
                             ),
                             
                             menuItem("Parameters selection", startExpanded = TRUE,
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
                                              value = 0.999
                                              )
                             ),
                             
                             
                             menuItem("Start preprocessing", startExpanded = TRUE,
                                 h5('Click the button to preprocess data,'),
                                 h5('and the progress bar will be'),
                                 h5('displayed in the lower right corner.'),
                                 br(),
                                 helpText('A 40Mb file takes about 6 mins.'),
                                 actionButton(inputId = "start_preprocess",
                                              label = "Start Preprocessing",
                                              icon = icon("play-circle"),
                                              style ="display: block; margin: 0 auto; width: 200px;color: black;"
                                              )   
                             )
                             #tags$hr()
                         )
                         
                         
                         
                         ),
        
        
        # default method sidebar
        conditionalPanel(condition = "input.main_tabs == 'default_method_condition'",
                         br(),
                         sidebarMenu(
                                # volcano sidebar
                               menuItem("Volcano plot",
                                         #tabName = 'default_method_volcano_show',
                                         br(),
                                         uiOutput('default_method_volcano_selector'),
                                         br(),
                                         actionButton(inputId = "default_method_volcano_Render",
                                                      label = "Render Plot",
                                                      icon = icon("play-circle"),
                                                      style ="display: block; margin: 0 auto; width: 200px;color: black;"
                                         ),
                                         br()
                                ),
                                
                                # heatmap sidebar
                               menuItem("Heatmap",
                                         #tabName = 'default_method_heatmap_show',
                                         br(),
                                         actionButton(inputId = "default_method_heatmap_Render",
                                                      label = "Render Plot",
                                                      icon = icon("play-circle"),
                                                      style ="display: block; margin: 0 auto; width: 200px;color: black;"
                                         ),
                                         br()
                                ),
                                
                                # qc sidebar
                               menuItem("QC plot",
                                         #tabName = 'default_method_qc_show',
                                         br(),
                                         uiOutput('default_method_qc_selector'),
                                         br(),
                                         actionButton(inputId = "default_method_qc_Render",
                                                      label = "Render Plot",
                                                      icon = icon("play-circle"),
                                                      style ="display: block; margin: 0 auto; width: 200px;color: black;"
                                         )
                                )
                         )
        ),
        
        # specific protein sidebar
        conditionalPanel(condition = "input.main_tabs == 'specific_protein_condition'",
                         br(),
                         h5('Select a specific protein: '),
                         br(),
                         uiOutput('specific_protein_selector'),
                         br(),
                         actionButton(inputId = "specific_protein_Render",
                                      label = "Render Plot",
                                      icon = icon("play-circle"),
                                      style ="display: block; margin: 0 auto; width: 200px;color: black;"
                         ),
                         br()

        )
),
  
    
    dashboardBody(
        tabsetPanel(  # a layout
            id = 'main_tabs',
            
            # data tab
            tabPanel(title = 'Data',
                     value = 'fileinput_condition',
                     fluidPage(
                        fluidRow(
                            column(10,
                                DT::DTOutput("contents", width = '80%')
                            )    
                        
                        )
                     )
            ),

            
            # default method condition tab
            tabPanel(title = 'MSstats method',
                     value = 'default_method_condition',
                     
                     fluidPage(
                          fluidRow(
                              column(6,
                                     plotOutput('default_method_volcano_out')
                                     #downloadButton(outputId = "default_method_volcano_downloader", 
                                     #                label = "Download Volcano plot",
                                     #                style="color: black;")
                              ),
                              column(6,
                                     plotOutput('default_method_heat_out')
                                     #downloadButton(outputId = "default_method_heatmap_downloader", 
                                     #                label = "Download Heatmap",
                                     #                style="color: black;")
                              ),
                          ),
                          
                          br(),
                            
                          fluidRow(
                              column(9, 
                                     #align="center",   # doesn't work
                                     plotOutput('default_method_qc_out')
                                     #downloadButton(outputId = "default_method_qc_downloader", 
                                     #                label = "Download QC plot",
                                     #                style="color: black;")
                              )
                          )

                   
                     )
            ),
            
            # specific protein tab
            tabPanel(title = 'Specific protein',
                     value = 'specific_protein_condition',
                     
                     fluidPage(
                       fluidRow(
                         column(6,
                                plotOutput('specific_protein_qcplot_out')
                         ),
                         column(6,
                                plotOutput('specific_protein_profileplot_out')
                         ),
                       ),
                       
                       br(),
                       
                       fluidRow(
                         column(9, 
                                #align="center",   # doesn't work
                                plotOutput('specific_protein_conditionplot_out')
                         )
                       )                      
                     )
                    )
            )
                     
        )
      
)



# Define server logic ----
server <- function(input, output, session) {
    options(shiny.maxRequestSize=500*1024^2)
    

    # -------------input---------------
    inputdf <- reactive({
        inFile <- input$csvFile
        if(is.null(inFile)) # the initialization should be NULL
            return(NULL)
        
        fileData <- read.csv(inFile$datapath)
    })
    
    
    output$contents <- renderDT({
        datatable(inputdf(),
                  options = list(scrollX = TRUE)
        )
    })
    
    
    ### progress function
    prePquant <- reactiveValues()
    
    observeEvent(input$start_preprocess, {
      
      fileData <- read.csv(input$csvFile$datapath)
      
      print("test")
      
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      
      progress$set(message = "Begin to preprocess data, please wait...", value = 0.1)
      setwd("../data/")
      
      prePquant$DDA2009.proposed <- MSstats::dataProcess(raw = fileData,
                                               logTrans = as.numeric(input$user_choose_logTrans),
                                               normalization = input$user_choose_normalization,
                                               summaryMethod = input$user_choose_summaryMethod,
                                               maxQuantileforCensored = input$user_choose_maxQuantileforCensored,
                                               censoredInt = "NA",
                                               MBimpute = TRUE)
      
      progress$set(message = "Begin to preprocess data, please wait...", value = 0.4)
      prePquant$DDA2009.TMP <- MSstats::dataProcess(raw = fileData,
                                          logTrans = as.numeric(input$user_choose_logTrans),
                                          normalization = input$user_choose_normalization,
                                          summaryMethod = input$user_choose_summaryMethod,
                                          maxQuantileforCensored = input$user_choose_maxQuantileforCensored,
                                          censoredInt = NULL,
                                          MBimpute = FALSE)
      
      progress$set(message = "Begin to generate group comparison, please wait...", value = 0.8)
      # Automatically create the manually created matrix in MSstats, user manual p23
      len <- length(levels(prePquant$DDA2009.TMP$FeatureLevelData$GROUP))
      
      ourMatrix <- matrix(c(0:0),nrow=len,ncol=len)
      diag(ourMatrix) = -1
      for(i in 1:len-1){
        ourMatrix[i,i+1] = 1
      }
      ourMatrix[len,1] = 1
      
      ourCondition <- levels(prePquant$DDA2009.TMP$ProteinLevelData$GROUP)
      len2 <- length(ourCondition)
      tmp <- matrix(ourCondition, nr=len2, nc=1)
      name <- matrix(nr=len2, nc=1)
      for(i in 1:len2-1){
        name[i,1] <- sprintf('%s-%s', tmp[i+1,1], tmp[i,1])
      }
      name[len2,1] <- sprintf('%s-%s', tmp[1,1], tmp[len2,1])
      
      row.names(ourMatrix) <- name
      #----------End of creation-----------
      colnames(ourMatrix) <- ourCondition
      prePquant$DDA2009.comparisons <- groupComparison(contrast.matrix = ourMatrix,
                                             data = prePquant$DDA2009.proposed)
 
      write.csv(prePquant$DDA2009.comparisons$ComparisonResult, file="MSstats_output.csv")
      
      #! /usr/bin/python
      #conda_install(packages = 'pandas') # If you are using it for the first time, you need to install the pandas package
      
      py_run_file('../py/MSstatas to pheatmap.py')
      
    
      progress$set(message = "Preprocessing is over.", value = 1)
    })
    
    
    
    # default_method volcano plot
    default_method_volcano_plot <- reactive({
        if(input$default_method_volcano_Render == 0) {
            return(NULL)
        }
        else {
            getPlot(inputdf(), flag = 'volcano', selector = input$default_method_volcano_input,
                    prePquant$DDA2009.proposed, prePquant$DDA2009.TMP, prePquant$DDA2009.comparisons)
        }
    })
    
    default_method_volcano_select <- reactive({
        if(input$csvFile == 0) {
            return(NULL)
        }
        else {
            default_method_vol_selector <- getSelector(inputdf(), flag = 'volcano',
                                                       prePquant$DDA2009.proposed, prePquant$DDA2009.TMP)
            selectInput(inputId = 'default_method_volcano_input',
                        label = 'Options',
                        choices = as.list(default_method_vol_selector)
            )
        }
    })
    
    output$default_method_volcano_selector <- renderUI({
        default_method_volcano_select()
    })
    
    output$default_method_volcano_out <- renderPlot({
        default_method_volcano_plot()
    })
    
    
    # default_method heatmap
    default_method_heatmap_plot <- reactive({
        if(input$default_method_heatmap_Render == 0) {
            return(NULL)
        }
        else {
            getPlot(inputdf(), flag = 'heat',
                    prePquant$DDA2009.proposed, prePquant$DDA2009.TMP, prePquant$DDA2009.comparisons)
        }
    })
    
    output$default_method_heat_out <- renderPlot({
        default_method_heatmap_plot()
    })
    
    
    # default_method qc plot
    default_method_qc_plot <- reactive({
        if(input$default_method_qc_Render == 0) {
            return(NULL)
        }
        else {
            getPlot(inputdf(), flag = 'qc', selector = input$default_method_qc_input,
                    prePquant$DDA2009.proposed, prePquant$DDA2009.TMP, prePquant$DDA2009.comparisons)
        }
    })
    
    default_method_qc_select <- reactive({
        if(input$csvFile == 0) {
            return(NULL)
        }
        else {
            default_method_qc_selector <- getSelector(inputdf(), flag = 'qc',
                                                      prePquant$DDA2009.proposed, prePquant$DDA2009.TMP)
            selectInput(inputId = 'default_method_qc_input',
                        label = 'Options',
                        choices = as.list(default_method_qc_selector)
            )
            }
    })
    
    output$default_method_qc_selector <- renderUI({
        default_method_qc_select()
    })
    
    output$default_method_qc_out <- renderPlot({
        default_method_qc_plot()        
    })
    

        
    ### specific protein
    output$specific_protein_selector <- renderUI({
      specific_protein_select()
    })
    
    specific_protein_select <- reactive({
      if(input$csvFile == 0) {
        return(NULL)
      }
      else {
        specific_protein_selector <- levels(prePquant$DDA2009.proposed$ProteinLevelData$Protein)
        selectInput(inputId = 'specific_protein_select_input',
                    label = 'Options',
                    choices = as.list(specific_protein_selector)
        )
      }
    })
    
    
    #specific_protein_qcplot_out
    specific_protein_qc_plot <- reactive({
      if(input$specific_protein_Render == 0) {
        return(NULL)
      }
      else {
        dataProcessPlots(data = prePquant$DDA2009.proposed, type="QCPlot",
                         which.Protein=input$specific_protein_select_input,
                         width=10, height=5, address=FALSE)
      }
    })
    
    output$specific_protein_qcplot_out <- renderPlot({
      specific_protein_qc_plot()
    })
    
    #specific_protein_profileplot_out
    specific_protein_profileplot_plot <- reactive({
      if(input$specific_protein_Render == 0) {
        return(NULL)
      }
      else {
        dataProcessPlots(data = prePquant$DDA2009.proposed, type="ProfilePlot",
                         which.Protein=input$specific_protein_select_input,
                         width=10, height=5, address=FALSE)
      }
    })
    
    output$specific_protein_profileplot_out <- renderPlot({
      specific_protein_profileplot_plot()
    })
    
    #specific_protein_conditionplot_out
    specific_protein_conditionplot_plot <- reactive({
      if(input$specific_protein_Render == 0) {
        return(NULL)
      }
      else {
        dataProcessPlots(data = prePquant$DDA2009.proposed, type="ConditionPlot",
                         which.Protein=input$specific_protein_select_input,
                         width=10, height=5, address=FALSE)
      }
    })
    
    output$specific_protein_conditionplot_out <- renderPlot({
      specific_protein_conditionplot_plot()
    })
     
}

# Run the app ----
shinyApp(ui = ui, server = server)
