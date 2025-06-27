library(shiny)
library(shinyBS)
library(shinyjs)
library(shinythemes)
library(DiscrimOD)

source("shiny_demo_Fidalgo_ui.R", local = TRUE)
source("shiny_demo_Arrhenius_ui.R", local = TRUE)
source("shiny_demo_Meeker_ui.R", local = TRUE)

source("shiny_demo_Fidalgo_server.R", local = TRUE)
source("shiny_demo_Arrhenius_server.R", local = TRUE)
source("shiny_demo_Meeker_server.R", local = TRUE)

ui <- fluidPage(
  theme = shinytheme("united"),
  useShinyjs(),
  
  tags$head(
    tags$style(HTML("
      body, label, input, button, select, h1, h2, h3, h4, h5, h6, table, td, th, .shiny-output-error {
        font-family: 'Times New Roman', Times, serif !important;
      }
    "))
  ),
  
  div(
    style = "background-color: #e0f9e4; padding: 30px 20px; text-align: center; 
             border-bottom: 2px solid #ccc; margin-bottom: 10px;",
    HTML("
      <h1 style='font-weight: bold; font-size: 32px; margin-bottom: 10px;'>
        Model Discrimination Design Generation for Accelerated Life Testing Experiments via Hybridized Optimization Algorithms
      </h1>
      <h4 style='margin: 5px 0;'>Kuan-Yuan Lin</h4>
      <h5 style='color: #555;'>Department of Statistics, National Taipei University, Taiwan</h5>
    ")
  ),
  tabsetPanel(
    tabPanel("Fidalgo", shiny_demo_Fidalgo_ui),
    tabPanel("Arrhenius", shiny_demo_Arrhenius_ui),
    tabPanel("Meeker", shiny_demo_Meeker_ui)
  ),
  
  tags$footer(
    style = "background-color: #e0f9e4; padding: 10px; text-align: right; font-size: 13px; color: #333;",
    HTML("&#9675; Maintainer: Kuan-Yuan Lin (<a href='mailto:a0921129003@gmail.com'>a0921129003@gmail.com</a>)")
  )
)

server <- function(input, output, session) {
  shiny_demo_Fidalgo_server(input, output, session)
  shiny_demo_Arrhenius_server(input, output, session)
  shiny_demo_Meeker_server(input, output, session)
}

shinyApp(ui, server)