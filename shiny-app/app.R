library(shiny)
library(shinydashboard)
library(DT)
library('MSstats', warn.conflicts = F, quietly = T, verbose = F)
library(reticulate)
library(pheatmap)

# Source getPlots functions
source("app/getPlots.R")
source("app/getSelector.R")


ui <- dashboardPage(
    dashboardHeader(title = "pQuantR"),
    dashboardSidebar(
        #useShinyjs(),
        
        # fileinput side bar menu
        conditionalPanel(condition = "input.main_tabs == 'fileinput_condition'",
                         h5('Note: this shiny app ...'),
                         tags$hr(), 
                         
                         # file selection box
                         fileInput('csvFile', 'Choose a \'out_msstats.csv\' file', multiple = FALSE, 
                                   accept=c('text/csv', 'text/comma-separated-values,text/plain')), # CSV text file
                         
                         helpText('Note: \'.mzTab\' is ... ')
        ),
        
        # volcano sidebar
        conditionalPanel(condition = "input.main_tabs == 'volcano_condition'",
                         br(),
                         helpText('Please wait a second until the \"Options\" appear, then click \"Render Plot\".'),
                         uiOutput('volcano_selector'),
                         br(),
                         actionButton(inputId = "volcano_Render",
                                      label = "Render Plot",
                                      icon = icon("play-circle"),
                                      style ="display: block; margin: 0 auto; width: 200px;color: black;"
                         )
        ),
        
        # heatmap sidebar
        conditionalPanel(condition = "input.main_tabs == 'heatmap_condition'",
                         br(),
                         actionButton(inputId = "heatmap_Render",
                                      label = "Render Plot",
                                      icon = icon("play-circle"),
                                      style ="display: block; margin: 0 auto; width: 200px;color: black;"
                         )
        ),
        
        # qc sidebar
        conditionalPanel(condition = "input.main_tabs == 'qc_condition'",
                         br(),
                         #div(style = "text: align-left; color: white,", tags$b("Current comparison")),
                         #verbatimTextOutput(outputId = "currentCompareText"),
                         #selectInput(inputId = 'qc_input',
                         #            label = 'Options',
                         #            choices = qc_input_options,
                         #            selectize = FALSE
                         #            )
                         helpText('Please wait a second until the \"Options\" appear, then click \"Render Plot\".'),
                         uiOutput('qc_selector'),
                         br(),
                         actionButton(inputId = "qc_Render",
                                      label = "Render Plot",
                                      icon = icon("play-circle"),
                                      style ="display: block; margin: 0 auto; width: 200px;color: black;"
                         )
        )
    ),
    
    
    dashboardBody(
        tabsetPanel(  # a layout
            id = 'main_tabs',
            
            # data tab
            tabPanel(title = 'data',
                     value = 'fileinput_condition',
                     box(DT::DTOutput("contents"))), # output as DT control
            
            # volcano tab
            tabPanel(title = 'Volcano plot',
                     value = 'volcano_condition',
                     plotOutput('volcano')),
            
            # heatmap tab
            tabPanel(title = 'Heatmap',
                     value = 'heatmap_condition',
                     plotOutput('heat')),
            
            # qc tab
            tabPanel(title = 'QC plot',
                     value = 'qc_condition',
                     plotOutput('qc'))
        )
    )
)

# Define server logic ----
server <- function(input, output) {
    
    inputdf <- reactive({
        inFile <- input$csvFile
        if(is.null(inFile)) # the initialization should be NULL
            return(NULL)
        
        fileData <- read.csv(inFile$datapath)
    })
    
    
    output$contents <- renderDT({
        inputdf()
    })
    
    
    # volcano plot
    volcano_plot <- reactive({
        if(input$volcano_Render == 0) {
            return(NULL)
        }
        else {
            getPlot(inputdf(), flag = 'volcano', selector = input$volcano_input)
        }
    })
    
    volcano_select <- reactive({
        if(input$csvFile == 0) {
            return(NULL)
        }
        else {
            vol_selector <- getSelector(inputdf(), flag = 'volcano')
            selectInput(inputId = 'volcano_input',
                        label = 'Options',
                        choices = as.list(vol_selector)
            )
        }
    })
    
    output$volcano_selector <- renderUI({
        volcano_select()
    })
    
    output$volcano <- renderPlot({
        volcano_plot()
    })
    
    
    # heatmap
    heatmap_plot <- reactive({
        if(input$heatmap_Render == 0) {
            return(NULL)
        }
        else {
            getPlot(inputdf(), flag = 'heat')
        }
    })
    
    output$heat <- renderPlot({
        heatmap_plot()
    })
    
    
    # qc plot
    qc_plot <- reactive({
        if(input$qc_Render == 0) {
            return(NULL)
        }
        else {
            getPlot(inputdf(), flag = 'qc', selector = input$qc_input)
        }
    })
    
    qc_select <- reactive({
        if(input$csvFile == 0) {
            return(NULL)
        }
        else {
            #output$currentCompareText <- renderPrint(input$qc_input)
            selector <- getSelector(inputdf(), flag = 'qc')
            selectInput(inputId = 'qc_input',
                        label = 'Options',
                        choices = as.list(selector)
            )
            }
    })
    
    output$qc_selector <- renderUI({
        qc_select()
    })
    
    output$qc <- renderPlot({
        qc_plot()        
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
