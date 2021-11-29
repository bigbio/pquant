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

library(org.Hs.eg.db)
library(org.Sc.sgd.db)
library(org.Rn.eg.db)

# Source getPlots functions
source("./app/getSelector.R")
source("./app/volcano_live.R")


LoadToEnvironment <- function(RData, env=new.env()) {
    load(RData, env)
    return(env)
}

  
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
                                                                                        ".RData"))),
                                                          tags$hr(),
                                                          actionButton(inputId = "uploaddata_reset",
                                                                       label = "Reset all data",
                                                                       width = 200,
                                                                       icon = icon("redo")),br()
                    
                                                 ),
                                                 
                                                 menuItem("Parameters selection", #startExpanded = TRUE,
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
                                                 
                                                 
                                                 menuItem("Start preprocessing",
                                                          icon = icon("play"),br(),
                                                          div(style = "text-align:left", 
                                                              "Save preprocessed data to local?"),     
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
                                                                                style ="display: block; margin: 0 auto; width: 200px;color: black;"),br()
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
            
                    
                    # dynamic volcano
                    conditionalPanel(condition = "input.main_tabs == 'dynamic_volcano'",br(),
                                     div(style = "text-align:center",
                                         "Use this part function",br(),
                                         "must upload a \'.csv\' file",br()),br(),
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
                                                          div(img(src = 'pipeline.svg', width = '50%'), style="text-align: center;",),
                                                          br(),br(),br(),
                                                          h4("Data not ready?", a(href="./data/out_msstats_example.csv", download="out_msstats_example", "Get demo data here"))),
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
                                                                                 DT::DTOutput("contents", width = '80%')))
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
                                                                 title = "Preprocessed data", vaule ="initialData_tabbox_preprocesseddata",
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
                      
                      
                      
                              # dynamic volcano tab
                              tabPanel(title = 'Dynamic volcano',
                                       value = 'dynamic_volcano',
                                       icon = icon("project-diagram"),
                                       fluidPage(
                                                 tabBox(# No title
                                                        id = "dynamic_tabbox1", selected = "dynamic_tabbox_data", width = 12,
                                                        tabPanel(
                                                                 title = "Data", value = "dynamic_tabbox_data",
                                                                 fluidRow(
                                                                          column(width = 8,
                                                                                 DT::DTOutput("dynamic_evidence_contents")),
                                                                          column(width = 4,
                                                                                 rHandsontableOutput("dynamic_define_metadata")))
                                                        ),
                                                        tabPanel(
                                                                title = "Plot", value = "dynamic_tabbox_plot",
                                                               
                                                                ### from Proteus: live.R
                                                                fluidRow(
                                                                         column(5,
                                                                                plotOutput("dynamic_plotVolcano_out", height = "600px", width = "100%", brush = "plot_brush",hover="plot_hover"),
                                                                                tableOutput("dynamic_significanceTable_out")),
                                                                         column(7,
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
                                                               
                                                                # Show main protein table
                                                                fluidRow(
                                                                         column(width = 12,
                                                                                DT::dataTableOutput("dynamic_allProteinTable_out")))
                                                       )
                                                   
                                                 )
                                       )
                              )#,
                      
                      
                              # enrichment analysis tab
                              #tabPanel(title = "Enrichment analysis",value = "enrichment_analsis")
                      
                      
                    )
                     
        )
      
)



# Define server logic ----
server <- function(input, output, session) {
    options(shiny.maxRequestSize=520*1024^2)
  
    prePquant <- reactiveValues(DDA2009.proposed = NULL,
                                DDA2009.comparisons = NULL,
                                inputData_rdata_out_msstats = NULL,
                                inputData_URL_out_msstats = NULL)
    
    renderCheck <- reactiveValues(initialData_qualitycontrolPlot = 0,
                                  initialData_profilecondition = 0,
                                  default_method_residual_qq = 0,
                                  default_method_volcano = 0,
                                  default_method_heatmap = 0,
                                  default_method_comparison = 0,
                                  dynamic_volcano = 0)
    
    volcanoLive <- reactiveValues(pdat = NULL,
                                  res = NULL,
                                  max_points = NULL)
    
    dataControl <- reactiveValues(annoStart = 0,
                                  annoSubmit = 0,
                                  inputData_state = NULL,
                                  inputData_rdata = NULL,
                                  inputData_URL_state = NULL,
                                  preprocess_judge = 0,
                                  db_type = NULL)

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
                if(input$inputData$type == "application/vnd.ms-excel" | input$inputData$type == "application/x-gzip"){  #upload a csv file
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
        
        env <- reactiveFileReader(1000, session, input$inputData$datapath, LoadToEnvironment)
        #print(names(env()))  [1] "inputData_rdata_out_msstats" "DDA2009.comparisons"         "DDA2009.TMP"                 "DDA2009.proposed"
        prePquant$inputData_rdata_out_msstats <- env()[[names(env())[1]]]
        prePquant$DDA2009.proposed <- env()[[names(env())[3]]]
        prePquant$DDA2009.comparisons <- env()[[names(env())[2]]]
        
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
    
    output$preprocessedData <- renderDT({
        if(is.null(prePquant$DDA2009.comparisons)){
            return(NULL)
        }
        else{
            datatable(prePquant$DDA2009.comparisons$ComparisonResult,
                      options = list(scrollX = TRUE))
        }
    })
    
    
    ### enable, disable part
    observe({
        if(is.null(dataControl$inputData_state)){
            disable("start_preprocess")}
        else{enable("start_preprocess")}
    })
    
    observe({
        if(is.null(prePquant$DDA2009.comparisons)){
            disable("download_Render")}
        else{enable("download_Render")}
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
            disable("dynamic_volcano_Render")}
        else{
            enable("dynamic_volcano_Render")}
    })
    
    
    ###start_preprocess 
    observeEvent(input$start_preprocess, {
      
        fileData <- inputdf()
        
        species <- unlist(lapply(strsplit(as.character(fileData[1,1]), "\\_"), "[", 2))
        if(species == "HUMAN"){
            db_type = org.Hs.eg.db
        }else if(species == "YEAST"){
            db_type = org.Sc.sgd.db
        }else if(species == "RAT"){
            db_type = org.Rn.eg.db
        }else{
            dataControl$db_type = "out"
        }

        if(!(is.null(dataControl$db_type))){
            sendSweetAlert(
                session = session,
                title = "Warning",
                text = "Currently, only 'Human', 'Yeast', and 'Rat' species are supported.",
                type = "warning", 
                closeOnClickOutside = TRUE,
                width = 400
            )
        }else{
            if(sum("Fraction" == colnames(fileData)) == 1 & input$user_choose_summaryMethod == 'TMP')
            {
                sendSweetAlert(
                    session = session,
                    title = "Warning",
                    text = "Your data is not suitable for TMP summary method, please consider switching to 'linear method'(low precision) or changing data.",
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
                
                progress$set(message = "Begin to preprocess data, please wait...", value = 0.2)
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
                prePquant$DDA2009.comparisons <- groupComparison(contrast.matrix = ourMatrix,
                                                                 data = prePquant$DDA2009.proposed,
                                                                 use_log_file = FALSE)
                       
                progress$set(message = "Begin to convert ID, please wait...", value = 0.8)
                tmp_origin <- lapply(strsplit(as.character(prePquant$DDA2009.comparisons$ComparisonResult$Protein), "\\;"), "[")    # type: list
                tmp_df <- as.data.frame(t(sapply(tmp_origin, "[", i = 1:max(sapply(tmp_origin, length)))))  # list to df
                tmp <- as.data.frame(matrix(nrow=length(rownames(tmp_df)),ncol=length(colnames(tmp_df))))
                for(i in array(1:length(colnames(tmp_df)))){
                    tmp[,i] <- unlist(lapply(strsplit(as.character(tmp_df[,i]), "\\|"), "[", 2))
                }
        
                mapping_geneid <- data.frame(ProteinName = tmp[,1])
                mapping_genename <- data.frame(ProteinName = tmp[,1])
    
                for(i in array(1:length(colnames(tmp_df)))){
                    uniKeys <- as.character(tmp[,i])
                    
                    tmp_geneid <- AnnotationDbi::mapIds(db_type, keys=uniKeys, column="ENTREZID", keytype="UNIPROT")
                    tmp_geneid <- data.frame(matrix(lapply(tmp_geneid, as.character)))
                    tmp_geneid <- unlist(lapply(tmp_geneid[,1],function(x) if(identical(x,character(0))) NA else x))
                    tmp_geneid <- data.frame("ENTREZID"=tmp_geneid)
                    mapping_geneid <- cbind(mapping_geneid, tmp_geneid)
                    
                    tmp_genename <- AnnotationDbi::mapIds(db_type, keys=uniKeys, column="ENTREZID", keytype="UNIPROT")
                    tmp_genename <- data.frame(matrix(lapply(tmp_genename, as.character)))
                    tmp_genename <- unlist(lapply(tmp_genename[,1],function(x) if(identical(x,character(0))) NA else x))
                    tmp_genename <- data.frame("GENENAME"=tmp_genename)
                    mapping_genename <- cbind(mapping_genename, tmp_genename)
                }
        
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
        
                prePquant$DDA2009.comparisons$ComparisonResult$ENTREZID <- mapping_geneid$ENTREZID
                prePquant$DDA2009.comparisons$ComparisonResult$GENENAME <- mapping_genename$GENENAME
        
                ###
                
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
                        
                        DDA2009.proposed <- isolate(prePquant$DDA2009.proposed)
                        DDA2009.comparisons <- isolate(prePquant$DDA2009.comparisons)
                        
                      
                        projectpathname <- paste(input$SaveProjectFolderPath, "/", input$SaveProjectDataName, sep = "")
                        dir.create(projectpathname)
                        
                        progress$set(message = "Begin to save data, please wait...", value = 0.5)
                        
                        inputData_rdata_out_msstats <- inputdf()
                        
                        save(DDA2009.proposed, DDA2009.comparisons, inputData_rdata_out_msstats,
                             file = paste(projectpathname,"/", "DDA2009.RData", sep = ""))
                        
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
        reset("download_name")
        
        reset("SaveProjectDataControl")
        reset("SaveProjectDataName")
        reset("SaveProjectFolderPath")
        
        prePquant$DDA2009.proposed <- NULL
        prePquant$DDA2009.comparisons <- NULL
        prePquant$inputData_rdata_out_msstats <- NULL
        prePquant$inputData_URL_out_msstats <- NULL
        
        renderCheck$initialData_qualitycontrolPlot <- 0
        renderCheck$initialData_profilecondition <- 0
        renderCheck$default_method_residual_qq <- 0
        renderCheck$default_method_volcano <- 0
        renderCheck$default_method_heatmap <- 0
        renderCheck$default_method_comparison <- 0
        renderCheck$dynamic_volcano <- 0
        
        dataControl$annoStart <- 0
        dataControl$annoSubmit <- 0
        dataControl$preprocess_judge <- 0
  
        volcanoLive$pdat <- NULL
        volcanoLive$res <- NULL
        volcanoLive$max_points <- NULL
    })
    
    

    ##### ---initialData_tabbox_plot---
    initialData_profilecondition_select <- reactive({
        if(is.null(prePquant$DDA2009.proposed)) {
            selectInput(inputId = 'initialData_profilecondition_select_input',
                        label = 'Please upload data first',
                        choices = NULL)
        }
        else {
            initialData_profilecondition_selector <- levels(prePquant$DDA2009.proposed$ProteinLevelData$Protein)
            selectInput(inputId = 'initialData_profilecondition_select_input',
                        label = 'Options',
                        choices = as.list(initialData_profilecondition_selector))
        }
    })
    
    initialData_qualitycontrol_select <- reactive({
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
                        choices = as.list(initialData_qualitycontrol_selector))
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
            dataProcessPlots(data = prePquant$DDA2009.proposed, type="QCPlot",
                             which.Protein=input$initialData_qualitycontrol_select_input,
                             width=10, height=5, address=FALSE) 
        } else { return(NULL) }
    })
    
    #initialData_profilePlot_out
    observeEvent(input$initialData_profilecondition_Render, {
        renderCheck$initialData_profilecondition <- 1
    })
    
    output$initialData_profilePlot_out <- renderPlot({
        if(renderCheck$initialData_profilecondition > 0){
            dataProcessPlots(data = prePquant$DDA2009.proposed, type="ProfilePlot",
                             which.Protein=input$initialData_profilecondition_select_input,
                             width=10, height=5, address=FALSE)
        } else { return(NULL) }
    })
    
    #initialData_conditionPlot_out
    output$initialData_conditionPlot_out <- renderPlot({
        if(renderCheck$initialData_profilecondition > 0){
            dataProcessPlots(data = prePquant$DDA2009.proposed, type="ConditionPlot",
                             which.Protein=input$initialData_profilecondition_select_input,
                             width=10, height=5, address=FALSE)
        } else { return(NULL) }
    })
    
    
    ### preprocessed data download
    output$download_Render <- downloadHandler(
        filename = function() {
            paste(input$download_name, ".csv", sep = "")
        },
        content = function(file) {
            write.csv(prePquant$DDA2009.comparisons$ComparisonResult, file, row.names = FALSE)
        }
    )
    
    
    
    ##### -----default method-----  
    
    ### model-based QC plots
    default_method_residual_qq_select <- reactive({
        if(is.null(prePquant$DDA2009.proposed)) {
            selectInput(inputId = 'default_method_residual_qq_select_input',
                        label = 'Please upload data first',
                        choices = NULL)
        }
        else {
            default_method_residual_qq_selector <- levels(prePquant$DDA2009.proposed$ProteinLevelData$Protein)
            selectInput(inputId = 'default_method_residual_qq_select_input',
                        label = 'Options',
                        choices = as.list(default_method_residual_qq_selector))
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
            modelBasedQCPlots(data = prePquant$DDA2009.comparisons, type="ResidualPlots",
                              which.Protein=input$default_method_residual_qq_select_input,
                              width=10, height=5, address=FALSE)
        } else { return(NULL) }
    })
    
    #q-q plot (quantile-quantile)
    output$default_method_qqPlot_out <- renderPlot({
        if(renderCheck$default_method_residual_qq > 0){
            modelBasedQCPlots(data = prePquant$DDA2009.comparisons, type="QQPlots",
                              which.Protein=input$default_method_residual_qq_select_input,
                              width=10, height=5, address=FALSE)
        } else { return(NULL) }
    })
    
    
    ### group comparison plots
    # default_method volcano plot
    default_method_volcano_select <- reactive({
        if(is.null(prePquant$DDA2009.proposed)) {
            selectInput(inputId = 'default_method_volcano_input',
                        label = 'Please upload data first',
                        choices = NULL)
        }
        else {
            default_method_vol_selector <- getSelector(inputdf(), flag = 'volcano',
                                                       prePquant$DDA2009.proposed)
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
            groupComparisonPlots(data = prePquant$DDA2009.comparisons$ComparisonResult, type = 'VolcanoPlot',
                                 which.Comparison=input$default_method_volcano_input,
                                 width=5, height=5, address=FALSE)
        } else { return(NULL) }
    })
    
    
    # default_method heatmap
    observeEvent(input$default_method_heatmap_Render, {
        renderCheck$default_method_heatmap <- 1
    })
    
    output$default_method_heat_out <- renderPlot({
        if(renderCheck$default_method_heatmap > 0){
            MS_output <- prePquant$DDA2009.comparisons$ComparisonResult
            MS_output <- MS_output[,1:3]
            MS_output <- MS_output[is.finite(as.numeric(as.character(MS_output$log2FC))),]
            MS_output <- MS_output %>% filter(between(log2FC,-4,4)) %>% spread(key = Label, value = log2FC)
            MS_output <- IDPmisc::NaRV.omit(MS_output)
            
            heatmap <- MS_output[,-1]
            rownames(heatmap) <- MS_output[,1]
            
            if (nrow(heatmap) > 50) {pheatmap(heatmap, show_rownames = F)}
            else {pheatmap(heatmap)}
        } else { return(NULL) }
    })
    
    
    # default_method comparison plot
    default_method_comparison_select <- reactive({
        if(is.null(prePquant$DDA2009.proposed)) {
            selectInput(inputId = 'default_method_comparison_select_input',
                        label = 'Please upload data first',
                        choices = NULL)
        }
        else {
            default_method_comparison_selector <- levels(prePquant$DDA2009.proposed$ProteinLevelData$Protein)
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
            groupComparisonPlots(data = prePquant$DDA2009.comparisons$ComparisonResult, type="ComparisonPlot",
                                 which.Protein=input$default_method_comparison_select_input,
                                 width=10, height=5, address=FALSE)
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
        
        volcanoLive$pdat <- NULL
        volcanoLive$res <- NULL
        volcanoLive$max_points <- NULL
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
        dynamic_data <- unique(fileData[,grepl("Reference|Condition", colnames(fileData))])
        dynamic_data <- dynamic_data[order(dynamic_data$Condition),]
        dynamic_condition <- sort(dynamic_data[,grepl("Condition", colnames(dynamic_data))])
        metadata <- data.frame(condition = dynamic_condition)
        return(metadata)
    })
    
    
    output$dynamic_define_metadata <- renderRHandsontable({
        if (is.null(input$inputData)) { return(NULL) }
        else if ((dataControl$annoStart > 0) & (input$inputData != 0)) {
          categorial_anno <- categorial_anno()
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
          
          volcanoLive$max_points = 100
          
          progress$set(message = "Preprocessing is over.", value = 1)
        }
    })
    
    
    
    ### dynamic volcano Plot
    
    dynamic_volcano_select <- reactive({
        if(is.null(input$inputData)) {
            selectInput(inputId = 'dynamic_volcano_input',
                        label = 'Please upload .csv data first',
                        choices = NULL)
        }
        else if(dataControl$annoSubmit == 0){
            selectInput(inputId = 'dynamic_volcano_input',
                        label = 'Please submit annotation first',
                        choices = NULL)
        }
        else {
            dynamic_vol_selector <- getSelector(inputdf(), flag = 'volcano',
                                                       prePquant$DDA2009.proposed)
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
        comparisons <- subset(prePquant$DDA2009.comparisons$ComparisonResult, Label==selector)
        volcanoLive$res <- comparisons
        
        volcanoLive$res$"-log10(pvalue)" <- -log10(volcanoLive$res$pvalue)
        
        rownames(volcanoLive$res) <- c(1:nrow(volcanoLive$res))
        
        progress$set(message = "Preprocessing is over.", value = 1)
    })
  
    
    # else
    output$dynamic_gap_out <- renderUI({HTML('<br/>')})
    
    output$dynamic_replicateTable_out <- dynamic_replicateTable(volcanoLive$res, input, volcanoLive$pdat, volcanoLive$max_points)
    output$dynamic_significanceTable_out <- dynamic_significanceTable(volcanoLive$res, volcanoLive$res, input)
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
        if(renderCheck$dynamic_volcano > 0) {
            # assume first column is id ("protein" or "peptide")
            idcol <- names(volcanoLive$res)[1]
            cols <- c(idcol, "log2FC", "pvalue", "adj.pvalue")
            d <- volcanoLive$res[, cols]
            d[, 2:ncol(d)] <- sapply(d[, 2:ncol(d)], function(x) signif(x, 3))
            d <- DT::datatable(d, class = 'cell-border strip hover')
            DT::formatStyle(d, 0, cursor = 'pointer')
        } else { return(NULL) }
    })

    
}

# Run the app ----
shinyApp(ui = ui, server = server)
