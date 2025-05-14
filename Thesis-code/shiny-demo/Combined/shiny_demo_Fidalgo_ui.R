# Define UI
shiny_demo_Fidalgo_ui <- fluidPage(
  theme = shinytheme("united"),
  
  titlePanel("MMM vs. MM under pharmacokinetic model (López-Fidalgo et al., 2007)"),
  
  sidebarLayout(
    sidebarPanel(
      tags$h4("Model Discrimination"),
      withMathJax(tags$p("True Model: $$\\frac{\\zeta_1 \\cdot x} {\\zeta_2 + x} + \\zeta_3\ x$$")),
      withMathJax(tags$p("Rival Model: $$\\frac{\\delta_1 \\cdot x} {\\delta_2 + x}$$")),
      
      tags$h4("Select Calculate Options"),
      selectInput("calculate_option_Fidalgo", "Divergence Method:",
                  choices = c("Close Form" = "close-form",
                              "Integration" = "integration"),
                  selected = "close-form"),
      
      selectInput("distribution_Fidalgo", "Select Distribution:",
                  choices = c("Log-Normal" = "lognormal", "Weibull" = "weibull"),
                  selected = "lognormal"),
      
      tags$h4("Design Space Bound Settings"),
      sliderInput("nSupp_Fidalgo", "Support Points:", min = 1, max = 10, value = 3, step = 1),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("dsLower_Fidalgo", "Lower Bound:", value = 0.1, step = 0.1, width = "45%"),
        numericInput("dsUpper_Fidalgo", "Upper Bound:", value = 5, step = 0.1, width = "45%")
      ),
      
      h4("True Model Mean Response Parameters"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("p1_Fidalgo", withMathJax("\\(\\zeta_1\\):"), value = 1, step = 0.1, width = "30%"),
        numericInput("p2_Fidalgo", withMathJax("\\(\\zeta_2\\):"), value = 1, step = 0.1, width = "30%"),
        numericInput("p3_Fidalgo", withMathJax("\\(\\zeta_3\\):"), value = 1, step = 0.1, width = "30%")
      ),
      
      h4("Dispersion Parameter"),
      numericInput("disp_Fidalgo", "Dispersion (\\(\\sigma\\)):", value = 1, step = 0.01),
      
      tags$h4("Rival Model Mean Response Parameters Bound"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("lower_delta_1_Fidalgo", "\\(\\delta_1\\) Lower:", value = 0.1, step = 0.1, width = "45%"),
        numericInput("upper_delta_1_Fidalgo", "\\(\\delta_1\\) Upper:", value = 100, step = 0.1, width = "45%")
      ),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("lower_delta_2_Fidalgo", "\\(\\delta_2\\) Lower:", value = 0.1, step = 0.1, width = "45%"),
        numericInput("upper_delta_2_Fidalgo", "\\(\\delta_2\\) Upper:", value = 100, step = 0.1, width = "45%")
      ),
      
      tags$h4("Algorithm Settings"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("nSwarm_outer_Fidalgo", "Outer PSO Particles:", value = 64, min = 1, step = 1, width = "45%"),
        numericInput("maxIter_outer_Fidalgo", "Outer PSO Iterations:", value = 5, min = 1, step = 5, width = "45%")
      ),
      checkboxInput("use_inner_pso_Fidalgo", "Use Inner PSO", value = FALSE),
      conditionalPanel(
        condition = "input.use_inner_pso_Fidalgo == true",
        tags$div(
          style = "display: flex; justify-content: space-between;",
          numericInput("nSwarm_inner_Fidalgo", "Inner PSO Particles:", value = 32, min = 1, step = 1, width = "45%"),
          numericInput("maxIter_inner_Fidalgo", "Inner PSO Iterations:", value = 100, min = 1, step = 1, width = "45%")
        )
      ),
      conditionalPanel(
        condition = "input.use_inner_pso_Fidalgo == false",
        numericInput("LBFGS_RETRY_Fidalgo", "BFGS Retries:", value = 2, min = 1, step = 1)
      ),
      
      actionButton("run_Fidalgo", "Run Optimization", style = "background-color:#FF5733; color:white; font-weight:bold;")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Result",
                 fluidRow(
                   column(6, 
                          h4("Support Points and Weights:"),
                          tableOutput("bestDesign_Fidalgo")),
                   column(6,
                          h4("Criterion Value and CPU Time:"),
                          tableOutput("optimizationSummary_Fidalgo"))
                 ),
                 h4("Model 1 Parameters:"),
                 tableOutput("model1_params_Fidalgo"),
                 h4("Model 2 Parameters:"),
                 tableOutput("model2_params_Fidalgo"),
                 h4("Directional Derivative Plot:"),
                 plotOutput("derivativePlot_Fidalgo")
        ),
        tabPanel("Setting",
                 h4("Selected Settings:"),
                 tableOutput("selectedSettings_Fidalgo")
        )
      )
    )
  )
)