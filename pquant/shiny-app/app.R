library(shiny)
library(shinydashboard)
library(DT)
library('MSstats', warn.conflicts = F, quietly = T, verbose = F)
library(reticulate)
library(pheatmap)
library(rhandsontable)
library(MSstats)
library(shinyjs)
library(shinyWidgets)
library(ggplot2)

setwd(getwd())

# Source getPlots functions
source("./app/getPlots.R")
source("./app/getSelector.R")
source("./app/volcano_live.R")


  
ui <- dashboardPage(
    dashboardHeader(title = "pQuantR"),
    dashboardSidebar(
        #useShinyjs(),
        
        # fileinput side bar menu
        conditionalPanel(condition = "input.main_tabs == 'fileinput_condition'",
                         sidebarMenu(
                             menuItem("Upload file",
                                 # file selection box
                                 fileInput('csvFile', 'Choose the \'out_msstats.csv\'', multiple = FALSE, 
                                           accept=c('text/csv', 'text/comma-separated-values,text/plain')), # CSV text file
                                 helpText('If the preprocessing error occurred,'),
                                 helpText('please change summary method.')
                             ),
                             
                             menuItem("Parameters selection", 
                                      #startExpanded = TRUE,
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
                             
                             
                             menuItem("Start preprocessing",
                                 h5('Click the button to preprocess data,'),
                                 h5('and the progress bar will be'),
                                 h5('displayed in the lower right corner.'),
                                 br(),
                                 helpText('A 40Mb file takes about 6 mins.'),
                                 actionButton(inputId = "start_preprocess",
                                              label = "Start Preprocessing",
                                              icon = icon("play-circle"),
                                              style ="display: block; margin: 0 auto; width: 200px;color: black;"
                                              ),
                                 br()
                             ),
                             
                             
                             menuItem("Data process plots",
                                      br(),
                                      helpText("Please wait until the data is preprocessed."),
                                      br(),
                                      menuItem("Profile & Condition plot",
                                               br(),
                                               uiOutput('initialData_profilecondition_selector'),
                                               br(),
                                               actionButton(inputId = "initialData_profilecondition_Render",
                                                            label = "Render Plot",
                                                            icon = icon("play-circle"),
                                                            style ="display: block; margin: 0 auto; width: 200px;color: black;"
                                               ),
                                               br()
                                      ),
                                      br(),
                                      menuItem("Quality control plot",
                                               br(),
                                               uiOutput('initialData_qualitycontrol_selector'),
                                               br(),
                                               actionButton(inputId = "initialData_qualitycontrol_Render",
                                                            label = "Render Plot",
                                                            icon = icon("play-circle"),
                                                            style ="display: block; margin: 0 auto; width: 200px;color: black;"
                                               ),
                                               br()
                                      ),
                                      br()
                             )
                             
                         )
                         
                         
                         
        ),
        
        
        # default method sidebar
        conditionalPanel(condition = "input.main_tabs == 'default_method_condition'",
                         br(),
                         sidebarMenu(
                               menuItem("Group comparison plots",
                                        br(),
                                        menuItem("Residual & Normal Q-Q plot",
                                                 br(),
                                                 uiOutput('default_method_residual_qq_selector'),
                                                 br(),
                                                 actionButton(inputId = "default_method_residual_qq_Render",
                                                              label = "Render Plot",
                                                              icon = icon("play-circle"),
                                                              style ="display: block; margin: 0 auto; width: 200px;color: black;"
                                                 ),
                                                 br()
                                        ),
                                        br()
                                ),
                               
                               menuItem("Model-based QC plots",
                                        br(),
                                        menuItem("Volcano plot",
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
                                        br(),
                                        menuItem("Heatmap",
                                                 br(),
                                                 actionButton(inputId = "default_method_heatmap_Render",
                                                              label = "Render Plot",
                                                              icon = icon("play-circle"),
                                                              style ="display: block; margin: 0 auto; width: 200px;color: black;"
                                                 ),
                                                 br()
                                        ),
                                        br(),
                                        menuItem("Comparison plot",
                                                 br(),
                                                 uiOutput('default_method_comparison_selector'),
                                                 br(),
                                                 actionButton(inputId = "default_method_comparison_Render",
                                                              label = "Render Plot",
                                                              icon = icon("play-circle"),
                                                              style ="display: block; margin: 0 auto; width: 200px;color: black;"
                                                 ),
                                                 br()
                                        ),
                                        br()
                               )
                         )
                         
                         
        ),

        
        # dynamic volcano
        conditionalPanel(condition = "input.main_tabs == 'dynamic_volcano'",
                         br(),
                         sidebarMenu(
                           menuItem("Assign condition",
                                    div(style = "text-align:center", "Click to render annotation table"),
                                    actionButton(inputId = "start_anno",
                                                 label = "Start",
                                                 icon = icon("play-circle"),
                                                 width = 200),
                                    br(),
                                    div(style = "text-align:center", "Your input in the second table",
                                        br(), "will be used for stastical",
                                        br(), "comparisons.",
                                        br(),br(), "Example:"),
                                    img(src = 'annotation_example.png', width = '100%'),
                                    div(style = "text-align:center", "once done renaming press",
                                        br(), "\'submit\' to lock in the annotation"),
                                    actionButton(inputId = "submit_anno",
                                                 label = ">  Submit",
                                                 width = 200),
                                    br(),br(),
                                    
                                    actionButton(inputId = "annoReset",
                                                 label = "Reset annotation",
                                                 width = 200,
                                                 icon = icon("redo"))
                           ),
                           menuItem("Volcano plot",
                                    br(),
                                    uiOutput('dynamic_volcano_selector'),
                                    actionButton(inputId = "dynamic_volcano_Render",
                                                 label = "Render Plot",
                                                 icon = icon("play-circle"),
                                                 style ="display: block; margin: 0 auto; width: 200px;color: black;"
                                    ),
                                    br()
                           )
                         )
        )
),
  
    
    dashboardBody(
        tabsetPanel(  # a layout
            id = 'main_tabs',
            
            # data tab
            tabPanel(title = 'Data',
                     value = 'fileinput_condition',
                     fluidPage(
                       tabBox(# No title
                         id = "initialData_tabbox", selected = "initialData_tabbox_data", width = 12,
                         tabPanel(
                           title = "Data", value = "initialData_tabbox_data",
                           fluidRow(
                             column(11,
                                    DT::DTOutput("contents", width = '80%')
                             )
                           )),
                         tabPanel(
                           title = "Data process plots", value = "initialData_tabbox_plot",
                           fluidRow(
                             column(6,
                                    plotOutput("initialData_profilePlot_out")),
                             column(6,
                                    plotOutput("initialData_qualitycontrolPlot_out"))
                           ),
                           br(),
                           fluidRow(
                             column(6,
                                    plotOutput("initialData_conditionPlot_out"))
                           )
                         )
                       )
                     )
                     
            ),
            
            
            # default method condition tab
            tabPanel(title = 'MSstats method',
                     value = 'default_method_condition',
                     tabBox(# No title
                       id = "msstats_tabbox", selected = "msstats_tabbox_modelbased", width = 12,
                       tabPanel(
                         title = "Model-based QC plots", value = "msstats_tabbox_modelbased",
                         fluidRow(
                           column(6,
                                  plotOutput("default_method_residualPlot_out")),
                           column(6,
                                  plotOutput("default_method_qqPlot_out"))
                         )),
                       tabPanel(
                         title = "Group comparison plots", value = "msstats_tabbox_groupcomparison",
                         fluidRow(
                           column(6,
                                  plotOutput("default_method_volcano_out")),
                           column(6,
                                  plotOutput("default_method_heat_out"))
                         ),
                         br(),
                         fluidRow(
                           column(6,
                                  plotOutput("default_method_comparisonPlot_out"))
                         )
                       )
                     )
            ),
            
            
            
            # dynamic volcano tab
            tabPanel(title = 'Dynamic volcano',
                     value = 'dynamic_volcano',
                     fluidPage(
                       tabBox(# No title
                         id = "dynamic_tabbox1", selected = "dynamic_tabbox_data", width = 12,
                         tabPanel(
                           title = "Data", value = "dynamic_tabbox_data",
                           fluidRow(
                             column(width = 8,
                                    DT::DTOutput("dynamic_evidence_contents")),
                             column(width = 4,
                                    rHandsontableOutput("dynamic_define_metadata"))
                           )),
                         tabPanel(
                           title = "Plot", value = "dynamic_tabbox_plot",
                           
                           ### from Proteus: live.R
                           fluidRow(
                             column(5, plotOutput("dynamic_plotVolcano_out", height = "600px", width = "100%", brush = "plot_brush",hover="plot_hover")),
                             column(7,
                                    fluidRow(
                                      column(4,
                                             radioButtons("intensityScale","Intesity Scale:",choices = c("Linear scale" = "","Log scale"="Log"),inline = TRUE)
                                      )
                                    ),
                                    fluidRow(
                                      column(4,
                                             fluidRow(plotOutput("dynamic_jitterPlot_out", height = "300px",width = "100%")),
                                             fluidRow(htmlOutput("dynamic_gap_out")),
                                             fluidRow(tableOutput("dynamic_significanceTable_out1")),
                                             fluidRow(tableOutput("dynamic_significanceTable_out2"))
                                      ),
                                      column(2,
                                             fluidRow(tableOutput("dynamic_replicateTable_out"))
                                      )
                                    )
                             )
                           ),
                           
                           # Show main protein table
                           fluidRow(
                             column(width = 12,
                                    DT::dataTableOutput("dynamic_allProteinTable_out"))
                           )
                         )
                         
                       )
                     )
            )
            
            
          )
                     
        )
      
)



# Define server logic ----
server <- function(input, output, session) {
    options(shiny.maxRequestSize=520*1024^2)
    

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
      
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      
      progress$set(message = "Begin to preprocess data, please wait...", value = 0.1)
      
      prePquant$DDA2009.proposed <- MSstats::dataProcess(raw = fileData,
                                               logTrans = as.numeric(input$user_choose_logTrans),
                                               normalization = input$user_choose_normalization,
                                               summaryMethod = input$user_choose_summaryMethod,
                                               maxQuantileforCensored = input$user_choose_maxQuantileforCensored,
                                               censoredInt = "NA",
                                               MBimpute = TRUE,
                                               use_log_file = FALSE)
      
      progress$set(message = "Begin to preprocess data, please wait...", value = 0.4)
      prePquant$DDA2009.TMP <- MSstats::dataProcess(raw = fileData,
                                          logTrans = as.numeric(input$user_choose_logTrans),
                                          normalization = input$user_choose_normalization,
                                          summaryMethod = input$user_choose_summaryMethod,
                                          maxQuantileforCensored = input$user_choose_maxQuantileforCensored,
                                          censoredInt = NULL,
                                          MBimpute = FALSE,
                                          use_log_file = FALSE)
      
      progress$set(message = "Begin to generate group comparison, please wait...", value = 0.7)
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
      #End of creation
      colnames(ourMatrix) <- ourCondition
      prePquant$DDA2009.comparisons <- groupComparison(contrast.matrix = ourMatrix,
                                             data = prePquant$DDA2009.proposed,
                                             use_log_file = FALSE)
 
      progress$set(message = "Preprocessing is over.", value = 1)
    })
    
    
    

    
    ### ---initialData_tabbox_plot---
    initialData_profilecondition_select <- reactive({
      if(input$csvFile == 0) {
        return(NULL)
      }
      else {
        initialData_profilecondition_selector <- levels(prePquant$DDA2009.proposed$ProteinLevelData$Protein)
        selectInput(inputId = 'initialData_profilecondition_select_input',
                    label = 'Options',
                    choices = as.list(initialData_profilecondition_selector)
        )
      }
    })
    
    initialData_qualitycontrol_select <- reactive({
      if(input$csvFile == 0) {
        return(NULL)
      }
      else {
        tmp <- levels(prePquant$DDA2009.proposed$ProteinLevelData$Protein)
        initialData_qualitycontrol_selector <- append('allonly', tmp, 1)
        selectInput(inputId = 'initialData_qualitycontrol_select_input',
                    label = 'Options',
                    choices = as.list(initialData_qualitycontrol_selector)
        )
      }
    })
    
    output$initialData_profilecondition_selector <- renderUI({
      initialData_profilecondition_select()
    })
    
    output$initialData_qualitycontrol_selector <- renderUI({
      initialData_qualitycontrol_select()
    })
    
    
    #initialData_qualitycontrolPlot_out
    initialData_qualitycontrolPlot <- reactive({
      if(input$initialData_qualitycontrol_Render == 0) {
        return(NULL)
      }
      else {
        dataProcessPlots(data = prePquant$DDA2009.proposed, type="QCPlot",
                         which.Protein=input$initialData_qualitycontrol_select_input,
                         width=10, height=5, address=FALSE)
      }
    })
    
    output$initialData_qualitycontrolPlot_out <- renderPlot({
      initialData_qualitycontrolPlot()
    })
    
    #initialData_profilePlot_out
    initialData_profilePlot <- reactive({
      if(input$initialData_profilecondition_Render == 0) {
        return(NULL)
      }
      else {
        dataProcessPlots(data = prePquant$DDA2009.proposed, type="ProfilePlot",
                         which.Protein=input$initialData_profilecondition_select_input,
                         width=10, height=5, address=FALSE)
      }
    })
    
    output$initialData_profilePlot_out <- renderPlot({
      initialData_profilePlot()
    })
    
    #initialData_conditionPlot_out
    initialData_conditionPlot <- reactive({
      if(input$initialData_profilecondition_Render == 0) {
        return(NULL)
      }
      else {
        dataProcessPlots(data = prePquant$DDA2009.proposed, type="ConditionPlot",
                         which.Protein=input$initialData_profilecondition_select_input,
                         width=10, height=5, address=FALSE)
      }
    })
    
    output$initialData_conditionPlot_out <- renderPlot({
      initialData_conditionPlot()
    })
    
    
    
    ##### -----default method-----  
    
    ### model-based QC plots
    default_method_residual_qq_select <- reactive({
      if(input$csvFile == 0) {
        return(NULL)
      }
      else {
        default_method_residual_qq_selector <- levels(prePquant$DDA2009.proposed$ProteinLevelData$Protein)
        selectInput(inputId = 'default_method_residual_qq_select_input',
                    label = 'Options',
                    choices = as.list(default_method_residual_qq_selector)
        )
      }
    })
 
    output$default_method_residual_qq_selector <- renderUI({
      default_method_residual_qq_select()
    })
 
       
    #residual plot
    default_method_residualPlot <- reactive({
      if(input$default_method_residual_qq_Render == 0) {
        return(NULL)
      }
      else {
        modelBasedQCPlots(data = prePquant$DDA2009.comparisons, type="ResidualPlots",
                          which.Protein=input$default_method_residual_qq_select_input,
                          width=10, height=5, address=FALSE)
      }
    })
    
    output$default_method_residualPlot_out <- renderPlot({
      default_method_residualPlot()
    })
    
    #q-q plot (quantile-quantile)
    default_method_qqPlot <- reactive({
      if(input$default_method_residual_qq_Render == 0) {
        return(NULL)
      }
      else {
        modelBasedQCPlots(data = prePquant$DDA2009.comparisons, type="QQPlots",
                          which.Protein=input$default_method_residual_qq_select_input,
                          width=10, height=5, address=FALSE)
      }
    })
    
    output$default_method_qqPlot_out <- renderPlot({
      default_method_qqPlot()
    })
    
    
    ### group comparison plots
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
            getPlot(inputdf(), flag = 'heat', selector = input$default_method_qc_input,
                    prePquant$DDA2009.proposed, prePquant$DDA2009.TMP, prePquant$DDA2009.comparisons)
        }
    })
    
    output$default_method_heat_out <- renderPlot({
        default_method_heatmap_plot()
    })
    
    #comparison plot
    default_method_comparison_select <- reactive({
      if(input$csvFile == 0) {
        return(NULL)
      }
      else {
        default_method_comparison_selector <- levels(prePquant$DDA2009.proposed$ProteinLevelData$Protein)
        selectInput(inputId = 'default_method_comparison_select_input',
                    label = 'Options',
                    choices = as.list(default_method_comparison_selector)
        )
      }
    })
  
    output$default_method_comparison_selector <- renderUI({
      default_method_comparison_select()
    })
    
    
    default_method_comparisonPlot <- reactive({
      if(input$default_method_comparison_Render == 0) {
        return(NULL)
      }
      else {
        groupComparisonPlots(data = prePquant$DDA2009.comparisons$ComparisonResult, type="ComparisonPlot",
                             which.Protein=input$default_method_comparison_select_input,
                             width=10, height=5, address=FALSE)
      }
    })
    
    output$default_method_comparisonPlot_out <- renderPlot({
      default_method_comparisonPlot()
    })
    
    
 
    
    ### ---dynamic volcano---
    
    
    # control data flow
    dataControl = reactiveValues(annoStart = 0,
                                 annoSubmit = 0
    )
    
    observeEvent(input$start_anno, {
      dataControl$annoStart <-  dataControl$annoStart + 1
    })
    
    observeEvent(input$submit_anno, {
      dataControl$annoSubmit <-  dataControl$annoSubmit + 1
    })
    
    
    # Reset metadata settings
    observeEvent(input$annoReset, {
      dataControl$annoStart <- 0
      dataControl$annoSubmit <- 0
    })
    
    
    ### start_convert_to_proteus
    dynamic_method_start_convert <- reactive({
      inFile <- input$csvFile
      if(is.null(inFile)) {
        return(NULL)
      }
      else {
        fileData <- read.csv(input$csvFile$datapath)
        dynamic_data <- sort(unique(fileData[, ncol(fileData)]))
        dynamic_data <- cbind(Reference=dynamic_data, measure="Intensity")
        dynamic_data
      }
    })
    
    output$dynamic_evidence_contents <- renderDT({
      datatable(dynamic_method_start_convert(),
                options = list(scrollX = TRUE)
      )
    })

    
    ### annotate metadata for user
    categorial_anno <- reactive({
      fileData <- read.csv(input$csvFile$datapath)
      dynamic_null <- rep('NULL', times = length(unique(fileData[, ncol(fileData)])))
      
      metadata <- data.frame(condition = dynamic_null)
      return(metadata)
    })
    
    
    output$dynamic_define_metadata <- renderRHandsontable({
      if (dataControl$annoStart > 0) {
        categorial_anno <- categorial_anno()
        categorial_anno$condition <- as.character(categorial_anno$condition)
        tab <- rhandsontable(categorial_anno)
        return(tab)
      } else {
        return(NULL)
      }
    })
    
    
    #get data from user
    #displays need to signal data is submitted
    dynamic_metadata <- reactive({
      if (dataControl$annoSubmit > 0) {
        reps <- isolate(input$dynamic_define_metadata)
        repsOut <- hot_to_r(reps)
        return(repsOut)
      } else {
        return(NULL)
      }
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
        disable("submit_anno")  # hide("dynamic_define_metadata")
      } else {
        enable("submit_anno")   # show("dynamic_define_metadata")
      }
    })
    
    # Proteus: create a peptide & protein dataset
    volcanoLive <- reactiveValues()
    
    observeEvent(input$submit_anno, {

      progress <- shiny::Progress$new()
      on.exit(progress$close())
      
      progress$set(message = "Begin to preprocess data, please wait...", value = 0.1)
      fileData <- read.csv(input$csvFile$datapath)
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
          if(npep >= min.peptides)
          {
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
    })
    
    observeEvent(input$dynamic_volcano_Render, {
      
      
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
    
    
    ### dynamic volcano Plot
    
    dynamic_volcano_select <- reactive({
      if(input$csvFile == 0) {
        return(NULL)
      }
      else {
        dynamic_vol_selector <- getSelector(inputdf(), flag = 'volcano',
                                                   prePquant$DDA2009.proposed, prePquant$DDA2009.TMP)
        selectInput(inputId = 'dynamic_volcano_input',
                    label = 'Options',
                    choices = as.list(dynamic_vol_selector)
        )
      }
    })
    
    output$dynamic_volcano_selector <- renderUI({
      dynamic_volcano_select()
    })
    
    
    
    # else
    output$dynamic_gap_out <- renderUI({HTML('<br/>')})
    output$dynamic_replicateTable_out <- dynamic_replicateTable(volcanoLive$res, input, volcanoLive$pdat, volcanoLive$max_points)
    output$dynamic_significanceTable_out1 <- dynamic_significanceTable1(volcanoLive$res, volcanoLive$res, input)
    output$dynamic_significanceTable_out2 <- dynamic_significanceTable2(volcanoLive$res, volcanoLive$res, input)
    output$dynamic_jitterPlot_out <- dynamic_jitterPlot(volcanoLive$res, input, volcanoLive$pdat, volcanoLive$max_points, dynamic_metadata())
    
    # Volcano plot
    output$dynamic_plotVolcano_out <- renderPlot({
      tab_idx <- as.numeric(input$allProteinTable_rows_selected)
      pVol <- dynamic_plotVolcano(volcanoLive$res, binhex=FALSE)
      if(length(tab_idx) > 0) {
        pVol <- pVol + geom_point(data=volcanoLive$res[tab_idx,], size=3, color='red')
      }
      pVol
    })
    
    output$dynamic_allProteinTable_out <- dynamic_allProteinTable(volcanoLive$res)
    
}

# Run the app ----
shinyApp(ui = ui, server = server)
