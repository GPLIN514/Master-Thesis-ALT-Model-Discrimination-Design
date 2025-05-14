library(shiny)
library(shinyjs)
library(shinythemes)

source("shiny_demo_Fidalgo_ui.R", local = TRUE)
source("shiny_demo_Arrhenius_ui.R", local = TRUE)
source("shiny_demo_Meeker_ui.R", local = TRUE)

source("shiny_demo_Fidalgo_server.R", local = TRUE)
source("shiny_demo_Arrhenius_server.R", local = TRUE)
source("shiny_demo_Meeker_server.R", local = TRUE)

ui <- fluidPage(
  theme = shinytheme("united"),
  titlePanel("Interactive Model Discrimination Design Platform"),
  tabsetPanel(
    tabPanel("Fidalgo", shiny_demo_Fidalgo_ui),
    tabPanel("Arrhenius", shiny_demo_Arrhenius_ui),
    tabPanel("Meeker", shiny_demo_Meeker_ui)
  )
)

server <- function(input, output, session) {
  shiny_demo_Fidalgo_server(input, output, session)
  shiny_demo_Arrhenius_server(input, output, session)
  shiny_demo_Meeker_server(input, output, session)
}

shinyApp(ui, server)