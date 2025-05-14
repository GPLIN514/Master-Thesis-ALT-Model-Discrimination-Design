# Define UI
shiny_demo_Meeker_ui <- fluidPage(
  theme = shinytheme("united"),
  useShinyjs(),
  
  titlePanel("Stress-Dependent Variance case under fatigue life model(Pascual and Meeker, 1997)"),
  
  sidebarLayout(
    sidebarPanel(
      tags$h4("Model Discrimination"),
      withMathJax(tags$p("True Model: $$\\zeta_1 + \\zeta_2 \\cdot \\log(x - 75.71)$$")),
      withMathJax(tags$p("Rival Model: $$\\delta_1 + \\delta_2 \\cdot \\log(x - 75.71)$$")),
      withMathJax(tags$p("True Model Dispersion: $$\\exp\\left\\{\\phi_1 + \\phi_2 \\cdot \\log(x - 75.71)\\right\\}$$")),
      withMathJax(tags$p("Rival Model Dispersion: $$\\exp\\left\\{\\kappa_1 + \\kappa_2 \\cdot \\log(x - 75.71)\\right\\}$$")),      
      
      
      tags$h4("Select Model Distributions"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        selectInput("distribution_M1_Meeker", "Model 1:",
                    choices = c("Log-Normal" = "lognormal", "Weibull" = "weibull"),
                    selected = "lognormal", width = "48%"),
        
        selectInput("distribution_M2_Meeker", "Model 2:",
                    choices = c("Log-Normal" = "lognormal", "Weibull" = "weibull"),
                    selected = "weibull", width = "48%")
      ),
      
      tags$h4("Design Space Bound Settings"),
      sliderInput("nSupp_Meeker", "Support Points:", min = 1, max = 10, value = 4, step = 1),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("dsLower_Meeker", "Lower Bound:", value = 76, step = 1, width = "45%"),
        numericInput("dsUpper_Meeker", "Upper Bound:", value = 150, step = 1, width = "45%")
      ),
      
      h4("True Model Mean Response Parameters"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("p1_Meeker", withMathJax("\\(\\zeta_1\\):"), value = 14.75, step = 0.01, width = "45%"),
        numericInput("p2_Meeker", withMathJax("\\(\\zeta_2\\):"), value = -1.39, step = 0.01, width = "45%")
      ),
      
      h4("True Model Dispersion Parameter"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("disp_phi1_Meeker", "Dispersion (\\(\\phi_1\\)):", value = 10.97, step = 0.01, width = "45%"),
        numericInput("disp_phi2_Meeker", "Dispersion (\\(\\phi_2\\)):", value = -2.50, step = 0.01, width = "45%")
      ),
      
      tags$h4("Rival Model Mean Response Parameters Bound"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("lower_delta_1_Meeker", "\\(\\delta_1\\) Lower:", value = 12.06, step = 0.01, width = "45%"),
        numericInput("upper_delta_1_Meeker", "\\(\\delta_1\\) Upper:", value = 17.44, step = 0.01, width = "45%")
      ),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("lower_delta_2_Meeker", "\\(\\delta_2\\) Lower:", value = -2.02, step = 0.01, width = "45%"),
        numericInput("upper_delta_2_Meeker", "\\(\\delta_2\\) Upper:", value = -0.76, step = 0.01, width = "45%")
      ),
      
      h4("Rival Model Dispersion Parameter Bound"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("lower_disp_kappa1_Meeker", "(\\(\\kappa_1\\)) Lower:", value = 10.00, step = 0.01, width = "45%"),
        numericInput("upper_disp_kappa1_Meeker", "(\\(\\kappa_1\\)) Upper:", value = 20.00, step = 0.01, width = "45%")
      ),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("lower_disp_kappa2_Meeker", "(\\(\\kappa_2\\)) Lower:", value = -3.00, step = 0.01, width = "45%"),
        numericInput("upper_disp_kappa2_Meeker", "(\\(\\kappa_2\\)) Upper:", value = -0.01, step = 0.01, width = "45%")
      ),
      
      h4("Censoring Threshold Setting"),
      numericInput("tc_Meeker", "Censoring Threshold (tc):", value = 1000, min = 1, step = 100),
      
      tags$h4("Algorithm Settings"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("nSwarm_outer_Meeker", "Outer PSO Particles:", value = 64, min = 1, step = 1, width = "45%"),
        numericInput("maxIter_outer_Meeker", "Outer PSO Iterations:", value = 5, min = 1, step = 5, width = "45%")
      ),
      checkboxInput("use_inner_pso_Meeker", "Use Inner PSO", value = FALSE),
      conditionalPanel(
        condition = "input.use_inner_pso_DeviceA == true",
        tags$div(
          style = "display: flex; justify-content: space-between;",
          numericInput("nSwarm_inner_Meeker", "Inner PSO Particles:", value = 32, min = 1, step = 1, width = "45%"),
          numericInput("maxIter_inner_Meeker", "Inner PSO Iterations:", value = 100, min = 1, step = 1, width = "45%")
        )
      ),
      conditionalPanel(
        condition = "input.use_inner_pso_Meeker == false",
        numericInput("LBFGS_RETRY_Meeker", "BFGS Retries:", value = 2, min = 1, step = 1)
      ),
      
      actionButton("run_Meeker", "Run Optimization", style = "background-color:#FF5733; color:white; font-weight:bold;")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Result",
                 fluidRow(
                   column(6, 
                          h4("Support Points and Weights:"),
                          tableOutput("bestDesign_Meeker")),
                   column(6,
                          h4("Criterion Value and CPU Time:"),
                          tableOutput("optimizationSummary_Meeker"))
                 ),
                 h4("Model 1 Parameters:"),
                 tableOutput("model1_params_Meeker"),
                 h4("Model 2 Parameters:"),
                 tableOutput("model2_params_Meeker"),
                 h4("Directional Derivative Plot:"),
                 plotOutput("derivativePlot_Meeker")
        ),
        tabPanel("Setting",
                 h4("Selected Settings:"),
                 tableOutput("selectedSettings_Meeker")
        )
      )
    )
  )
)