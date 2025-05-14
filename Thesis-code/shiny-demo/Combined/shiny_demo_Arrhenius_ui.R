# Define UI
shiny_demo_Arrhenius_ui <- fluidPage(
  theme = shinytheme("united"),
  
  titlePanel("Quadratic vs. Linear under Arrhenius model with Type I censored data"),
  
  sidebarLayout(
    sidebarPanel(
      tags$h4("Model Discrimination"),
      withMathJax(tags$p("True Model: $$\\zeta_1 + \\zeta_2 \\cdot \\frac{11605}{x + 273.15} + \\zeta_3\ \\cdot \\left( \\frac{11605}{x + 273.15} \\right)^2$$")),
      withMathJax(tags$p("Rival Model: $$\\delta_1 + \\delta_2 \\cdot \\frac{11605}{x + 273.15}$$")),
      
      tags$h4("Select Divergence Measure"),
      selectInput("divergence_Arrhenius", "Divergence Method:",
                  choices = c("KL Divergence" = "kl",
                              "LW Divergence" = "lw",
                              "Bhattacharyya Distance" = "b",
                              "Chi-square Distance" = "chi"),
                  selected = "kl"),
      
      selectInput("distribution_Arrhenius", "Select Distribution:",
                  choices = c("Log-Normal" = "lognormal", "Weibull" = "weibull"),
                  selected = "lognormal"),
      
      tags$h4("Design Space Bound Settings"),
      sliderInput("nSupp_Arrhenius", "Support Points:", min = 1, max = 10, value = 3, step = 1),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("dsLower_Arrhenius", "Lower Bound:", value = 10, step = 1, width = "45%"),
        numericInput("dsUpper_Arrhenius", "Upper Bound:", value = 80, step = 1, width = "45%")
      ),
      
      h4("True Model Mean Response Parameters"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("p1_Arrhenius", withMathJax("\\(\\zeta_1\\):"), value = -5, step = 0.1, width = "30%"),
        numericInput("p2_Arrhenius", withMathJax("\\(\\zeta_2\\):"), value = -1.5, step = 0.1, width = "30%"),
        numericInput("p3_Arrhenius", withMathJax("\\(\\zeta_3\\):"), value = 0.05, step = 0.01, width = "30%")
      ),
      
      h4("True Model Dispersion Parameter"),
      numericInput("disp1_Arrhenius", "Dispersion (\\(\\sigma_1\\)):", value = 0.9780103, step = 0.01),
      
      tags$h4("Rival Model Mean Response Parameters Bound"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("lower_delta_1_Arrhenius", "\\(\\delta_1\\) Lower:", value = -100, step = 0.1, width = "45%"),
        numericInput("upper_delta_1_Arrhenius", "\\(\\delta_1\\) Upper:", value = -10, step = 0.1, width = "45%")
      ),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("lower_delta_2_Arrhenius", "\\(\\delta_2\\) Lower:", value = 0.1, step = 0.1, width = "45%"),
        numericInput("upper_delta_2_Arrhenius", "\\(\\delta_2\\) Upper:", value = 5, step = 0.1, width = "45%")
      ),
      
      h4("Rival Model Dispersion Parameter Bound"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("lower_disp2_Arrhenius", "\\(\\sigma_2\\) Lower:", value = 0.9780103, step = 0.01, width = "45%"),
        numericInput("upper_disp2_Arrhenius", "\\(\\sigma_2\\) Upper:", value = 0.9780103, step = 0.01, width = "45%")
      ),
      
      h4("Censoring Threshold Setting"),
      numericInput("tc_Arrhenius", "Censoring Threshold (tc):", value = 5000, min = 1, step = 100),
      
      tags$h4("Algorithm Settings"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("nSwarm_outer_Arrhenius", "Outer PSO Particles:", value = 64, min = 1, step = 1, width = "45%"),
        numericInput("maxIter_outer_Arrhenius", "Outer PSO Iterations:", value = 5, min = 1, step = 5, width = "45%")
      ),
      checkboxInput("use_inner_pso_Arrhenius", "Use Inner PSO", value = FALSE),
      conditionalPanel(
        condition = "input.use_inner_pso_Arrhenius == true",
        tags$div(
          style = "display: flex; justify-content: space-between;",
          numericInput("nSwarm_inner_Arrhenius", "Inner PSO Particles:", value = 32, min = 1, step = 1, width = "45%"),
          numericInput("maxIter_inner_Arrhenius", "Inner PSO Iterations:", value = 100, min = 1, step = 1, width = "45%")
        )
      ),
      conditionalPanel(
        condition = "input.use_inner_pso_Arrhenius == false",
        numericInput("LBFGS_RETRY_Arrhenius", "BFGS Retries:", value = 2, min = 1, step = 1)
      ),
      
      actionButton("run_Arrhenius", "Run Optimization", style = "background-color:#FF5733; color:white; font-weight:bold;")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Result",
                 fluidRow(
                   column(6, 
                          h4("Support Points and Weights:"),
                          tableOutput("bestDesign_Arrhenius")),
                   column(6,
                          h4("Criterion Value and CPU Time:"),
                          tableOutput("optimizationSummary_Arrhenius"))
                 ),
                 h4("Model 1 Parameters:"),
                 tableOutput("model1_params_Arrhenius"),
                 h4("Model 2 Parameters:"),
                 tableOutput("model2_params_Arrhenius"),
                 h4("Directional Derivative Plot:"),
                 plotOutput("derivativePlot_Arrhenius")
        ),
        tabPanel("Setting",
                 h4("Selected Settings:"),
                 tableOutput("selectedSettings_Arrhenius")
        )
      )
    )
  )
)