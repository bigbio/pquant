library(shiny)
library(shinydashboard)
library(DT)
library('MSstats', warn.conflicts = F, quietly = T, verbose = F)

# Source getPlots functions
source("app/getPlots.R")


ui <- dashboardPage(
    dashboardHeader(title = "pQuantR"),
    dashboardSidebar(
        #描述shiny app
        h5('Note: this shiny app ...'),
        
        # 水平线条
        tags$hr(), 
        
        # 文件选择框
        fileInput('csvFile', 'Choose a \'out_msstats.csv\' file', multiple = FALSE, 
                  accept=c('text/csv', 'text/comma-separated-values,text/plain')), # CSV文本文件
        
        #描述文件
        helpText('Note: \'.mzTab\' is ... ')
    ),
    
    dashboardBody(
        tabsetPanel(  # 一种布局
            tabPanel('data',
                     box(DT::DTOutput("contents"))), # 以DT控件输出
            tabPanel('Volcano plot',
                     plotOutput('volcano')),
            tabPanel('Heatmap',
                     plotOutput('heat')),
            tabPanel('QC plot',
                     plotOutput('qc'))
    )
)
)

# Define server logic ----
server <- function(input, output) {
    
    inputdf <- reactive({
        inFile <- input$csvFile
    if(is.null(inFile)) # 初始应该为 NULL
        return(NULL)
    
    fileData <- read.csv(inFile$datapath)
    })
    
    
    output$contents <- renderDT({
        inputdf()
    })
    
    output$volcano <- renderPlot({
        getPlot(inputdf(), flag = 'volcano')
    })
    
    output$heat <- renderPlot({
        getPlot(inputdf(), flag = 'heat')
    })
    
    output$qc <- renderPlot({
        getPlot(inputdf(), flag = 'qc')         
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
