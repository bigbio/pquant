library(shiny)
library(shinydashboard)
library(DT)
library('MSstats', warn.conflicts = F, quietly = T, verbose = F)
library(reticulate)
library(pheatmap)
library(rhandsontable)


# Source getPlots functions
source("../shiny-app/app/getPlots.R")
source("../shiny-app/app/getSelector.R")

setwd(getwd())


pquant_shiny <- function(DDA2009.proposed, DDA2009.TMP, DDA2009.comparisons){
  
ui <- dashboardPage(
    dashboardHeader(title = "pQuantR"),
    dashboardSidebar(
        #useShinyjs(),
        
        # fileinput side bar menu
        conditionalPanel(condition = "input.main_tabs == 'fileinput_condition'",
                         h5('Note: this shiny app ...'),
                         tags$hr(), 
                         
                         # file selection box
                         fileInput('csvFile', 'Choose the same \'out_msstats.csv\'', multiple = FALSE, 
                                   accept=c('text/csv', 'text/comma-separated-values,text/plain')), # CSV text file
                         br(),
                         
                         helpText('Note: \'out_msstats.csv\' is ... ')
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
    
    
    # default_method volcano plot
    default_method_volcano_plot <- reactive({
        if(input$default_method_volcano_Render == 0) {
            return(NULL)
        }
        else {
            getPlot(inputdf(), flag = 'volcano', selector = input$default_method_volcano_input,
                    DDA2009.proposed, DDA2009.TMP, DDA2009.comparisons)
        }
    })
    
    default_method_volcano_select <- reactive({
        if(input$csvFile == 0) {
            return(NULL)
        }
        else {
            default_method_vol_selector <- getSelector(inputdf(), flag = 'volcano',
                                                       DDA2009.proposed, DDA2009.TMP)
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
            getPlot(inputdf(), flag = 'heat', DDA2009.proposed, DDA2009.TMP, DDA2009.comparisons)
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
                    DDA2009.proposed, DDA2009.TMP, DDA2009.comparisons)
        }
    })
    
    default_method_qc_select <- reactive({
        if(input$csvFile == 0) {
            return(NULL)
        }
        else {
            default_method_qc_selector <- getSelector(inputdf(), flag = 'qc',
                                                      DDA2009.proposed, DDA2009.TMP)
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
     
}

# Run the app ----
shinyApp(ui = ui, server = server)
}
