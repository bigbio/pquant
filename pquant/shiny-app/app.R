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
library(proteus)
library(ggplot2)

setwd(getwd())

# Source getPlots functions
source("./app/getPlots.R")
source("./app/getSelector.R")
source("./app/proteus_live.R")


  
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
                                           accept=c('text/csv', 'text/comma-separated-values,text/plain')) # CSV text file
                                 #helpText(' Note: \'out_msstats.csv\' is ... ')
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

        ),
        
        # proteus volcano sidebar
        conditionalPanel(condition = "input.main_tabs == 'proteus_volcano'",
                         br(),
                         sidebarMenu(
                           menuItem("Data handling",
                                    helpText('Click to convert \'out.csv\' to a'),
                                    helpText('Proteus-compliant input format.'),
                                    br(),
                                    helpText('Make sure you have \'out_mztab.csv\''),
                                    helpText('in your data file folder.'),
                                    actionButton(inputId = "proteus_method_start_Render",
                                                 label = "Start",
                                                 icon = icon("play-circle"),
                                                 width = 200),
                                    br()
                           ),
                           menuItem("Assign metadata",
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
                                    #imageOutput('proteus_annotation_example'),
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
                                    uiOutput('proteus_volcano_selector'),
                                    actionButton(inputId = "proteus_volcano_Render",
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
                              column(5,
                                     plotOutput('default_method_volcano_out')
                                     #downloadButton(outputId = "default_method_volcano_downloader", 
                                     #                label = "Download Volcano plot",
                                     #                style="color: black;")
                              )
                          ),
                          
                          br(),
                          
                          fluidRow(
                              column(5,
                                     plotOutput('default_method_heat_out')
                                     #downloadButton(outputId = "default_method_heatmap_downloader", 
                                     #                label = "Download Heatmap",
                                     #                style="color: black;")
                              )
                          ),
                          
                          br(),
                            
                          fluidRow(
                              column(5, 
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
                                plotOutput('specific_protein_conditionplot_out')
                         )
                       )                      
                     )
            ),
            
            # proteus volcano tab
            tabPanel(title = 'Proteus volcano',
                     value = 'proteus_volcano',
                     fluidPage(
                       tabBox(# No title
                         id = "proteus_tabbox1", selected = "proteus_tabbox_data", width = 12,
                         tabPanel(
                           title = "Data", value = "proteus_tabbox_data",
                           fluidRow(
                             column(width = 8,
                                    DT::DTOutput("proteus_evidence_contents")),
                             column(width = 4,
                                    rHandsontableOutput("proteus_define_metadata"))
                           )),
                         tabPanel(
                           title = "Plot", value = "proteus_tabbox_plot",
                           
                           ### from Proteus: live.R
                           fluidRow(
                             column(5, plotOutput("proteus_plotVolcano_out", height = "600px", width = "100%", brush = "plot_brush",hover="plot_hover")),
                             column(7,
                                    fluidRow(tableOutput("proteus_proteinInfo_out")),
                                    fluidRow(
                                      column(4,
                                             radioButtons("intensityScale","Intesity Scale:",choices = c("Linear scale" = "","Log scale"="Log"),inline = TRUE)
                                      )
                                    ),
                                    fluidRow(
                                      column(4,
                                             fluidRow(plotOutput("proteus_jitterPlot_out", height = "300px",width = "100%")),
                                             fluidRow(htmlOutput("proteus_gap_out")),
                                             fluidRow(tableOutput("proteus_significanceTable_out"))
                                      ),
                                      column(3,
                                             fluidRow(tableOutput("proteus_replicateTable_out"))
                                      )
                                    )
                             )
                           ),
                           
                           # Show main protein table
                           fluidRow(
                             column(width = 12,
                                    DT::dataTableOutput("proteus_allProteinTable_out"))
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
      
      progress$set(message = "Begin to preprocess data, please wait...", value = 0.3)
      prePquant$DDA2009.TMP <- MSstats::dataProcess(raw = fileData,
                                          logTrans = as.numeric(input$user_choose_logTrans),
                                          normalization = input$user_choose_normalization,
                                          summaryMethod = input$user_choose_summaryMethod,
                                          maxQuantileforCensored = input$user_choose_maxQuantileforCensored,
                                          censoredInt = NULL,
                                          MBimpute = FALSE,
                                          use_log_file = FALSE)
      
      progress$set(message = "Begin to generate group comparison, please wait...", value = 0.6)
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
                                             data = prePquant$DDA2009.proposed,
                                             use_log_file = FALSE)
 
      
      #! /usr/bin/python
      #conda_install(packages = 'pandas') # If you are using it for the first time, you need to install the pandas package

      progress$set(message = "Begin to preprocess proteus data, please wait...", value = 0.8)
      
      setwd('../data/')
      reticulate::py_run_file('../py/get_proteus_evidence.py')
      proteus_evidence_file <- read.csv('out_proteus.csv', row.names = NULL)
      # Prevent the occurrence of ''(null value)
      proteus_evidence_file <- unique(proteus_evidence_file)
      proteus_evidence_file[proteus_evidence_file == ''] <- NA
      proteus_evidence_file <- na.omit(proteus_evidence_file)
      write.csv(proteus_evidence_file, 'out_proteus.csv', row.names = F)
      
      prePquant$out_proteus <- proteus_evidence_file
    
      progress$set(message = "Preprocessing is over.", value = 1)
      setwd(getwd())
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
            getPlot(inputdf(), flag = 'heat', selector = input$default_method_qc_input,
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
    
    
    
    
    ### ---Proteus volcano---
    
    
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
    proteus_method_start_convert <- reactive({
      if(input$proteus_method_start_Render == 0) {
        return(NULL)
      }
      else {
        proteus_evidence_file <- prePquant$out_proteus
      }
    })
    
    output$proteus_evidence_contents <- renderDT({
      datatable(proteus_method_start_convert(),
                options = list(scrollX = TRUE)
      )
    })
    
    
    ### annotate metadata for user
    categorial_anno <- reactive({
      proteus_data <- prePquant$out_proteus
      proteus_experiment <- unique(proteus_data$experiment)
      proteus_intensity <- rep('Intensity', times = length(proteus_experiment))
      proteus_null <- rep('NULL', times = length(proteus_experiment))
      
      metadata <- data.frame(experiment = proteus_experiment,
                             measure = proteus_intensity,
                             sample = proteus_experiment,
                             condition = proteus_null,
                             replicate_notNecessary = proteus_null
      )
      return(metadata)
    })
    
    
    output$proteus_define_metadata <- renderRHandsontable({
      if (dataControl$annoStart > 0) {
        categorial_anno <- categorial_anno()
        categorial_anno$experiment <- as.character(categorial_anno$experiment)
        categorial_anno$measure <- as.character(categorial_anno$measure)
        categorial_anno$sample <- as.character(categorial_anno$experiment)
        categorial_anno$condition <- as.character(categorial_anno$condition)
        categorial_anno$replicate_notNecessary <- as.character(categorial_anno$replicate_notNecessary)
        tab <- rhandsontable(categorial_anno) %>%
          hot_col("experiment", readOnly = T)
        return(tab)
      } else {
        return(NULL)
      }
    })
    
    
    #get data from user
    #displays need to signal data is submitted
    proteus_metadata <- reactive({
      if (dataControl$annoSubmit > 0) {
        reps <- isolate(input$proteus_define_metadata)
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
        disable("submit_anno")  # hide("proteus_define_metadata")
      } else {
        enable("submit_anno")   # show("proteus_define_metadata")
      }
    })
    
    # Proteus: create a peptide & protein dataset
    proteusLive <- reactiveValues()
    
    observeEvent(input$submit_anno, {
      
      evi <- prePquant$out_proteus
      
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      
      progress$set(message = "Begin to preprocess data, please wait...", value = 0.1)
      meta <- proteus_metadata()
      
      progress$set(message = "Begin to preprocess data, please wait...", value = 0.3)
      proteus_pepdat <- makePeptideTable(evi, meta)
      
      progress$set(message = "Begin to preprocess data, please wait...", value = 0.5)
      proteus_prodat <- makeProteinTable(proteus_pepdat)

      progress$set(message = "Begin to preprocess data, please wait...", value = 0.7)
      proteusLive$pdat <- normalizeData(proteus_prodat)
      
      proteusLive$max_points = 100
    })
    
    
    observeEvent(input$proteus_volcano_Render, {
      
      
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      
      progress$set(message = "Begin to preprocess data, please wait...", value = 0.9)
      tmp_selector = input$proteus_volcano_input
      selector = unlist(strsplit(tmp_selector, split=' '))
      print(selector)
      proteusLive$res <- limmaDE(proteusLive$pdat, conditions=selector)
      
      proteusLive$res$"-log10(P.Value)" <- -log10(proteusLive$res$P.Value)
      progress$set(message = "Preprocessing is over.", value = 1)
    })
    
    
    
    ### Proteus Plot
    
    proteus_volcano_select <- reactive({
      if(input$csvFile == 0) {
        return(NULL)
      }
      else {
        con_list <- proteusLive$pdat[["conditions"]]
        tmp_selector <- c()
        x = 1
        for (i in 1:(length(con_list)-1)){
          for (j in (i+1):length(con_list)){
            tmp <- paste(con_list[i],con_list[j],collapse=' ')
            tmp_selector[x] <- tmp
            x=x+1
          }
        }
        proteus_vol_selector <- tmp_selector
        selectInput(inputId = 'proteus_volcano_input',
                    label = 'Options',
                    choices = as.list(proteus_vol_selector)
        )
      }
    })
    
    output$proteus_volcano_selector <- renderUI({
      proteus_volcano_select()
    })
    
    
    
    # else
    output$proteus_proteinInfo_out <- proteus_proteinInfo(proteusLive$res, input, proteusLive$pdat, proteusLive$max_points)
    output$proteus_gap_out <- renderUI({HTML('<br/>')})
    output$proteus_replicateTable_out <- proteus_replicateTable(proteusLive$res, input, proteusLive$pdat, proteusLive$max_points)
    output$proteus_significanceTable_out <- proteus_significanceTable(proteusLive$res, proteusLive$res, input)
    output$proteus_jitterPlot_out <- proteus_jitterPlot(proteusLive$res, input, proteusLive$pdat, proteusLive$max_points)
    
    # Volcano plot
    output$proteus_plotVolcano_out <- renderPlot({
      tab_idx <- as.numeric(input$allProteinTable_rows_selected)
      pVol <- plotVolcano(proteusLive$res, binhex=FALSE)
      if(length(tab_idx) > 0) {
        pVol <- pVol + geom_point(data=proteusLive$res[tab_idx,], size=3, color='red')
      }
      pVol
    })
    
    output$proteus_allProteinTable_out <- proteus_allProteinTable(proteusLive$res)
    
}

# Run the app ----
shinyApp(ui = ui, server = server)
