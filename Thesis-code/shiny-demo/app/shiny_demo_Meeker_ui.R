# Define UI
shiny_demo_Meeker_ui <- fluidPage(
  theme = shinytheme("united"),
  useShinyjs(),
  
  tags$head(
    tags$style(HTML("
      .loader {
        border: 6px solid #f3f3f3;
        border-top: 6px solid #FF5733;
        border-radius: 50%;
        width: 32px;
        height: 32px;
        animation: spin 1s linear infinite;
      }

      @keyframes spin {
        0% { transform: rotate(0deg); }
        100% { transform: rotate(360deg); }
      }
    "))
  ),
  
  titlePanel(HTML("Stress-Dependent Variance case under fatigue life model(Pascual and Meeker, 1997)")),
  
  sidebarLayout(
    sidebarPanel(
      style = "background-color: #fffef3;
         border: 1px solid #ffc107;
         border-radius: 10px;
         padding: 15px;
         margin-bottom: 15px;
         box-shadow: 2px 2px 6px rgba(0,0,0,0.1);
         font-size: 16px;",
      
      tags$h4("🧠 Select Model Distributions"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        selectInput("distribution_M1_Meeker", "True Model:",
                    choices = c("Log-Normal" = "lognormal", "Weibull" = "weibull"),
                    selected = "lognormal", width = "48%"),
        
        selectInput("distribution_M2_Meeker", "Rival Model:",
                    choices = c("Log-Normal" = "lognormal", "Weibull" = "weibull"),
                    selected = "weibull", width = "48%")
      ),
      
      tags$h4("🧭 Design Space Bound Settings"),
      sliderInput("nSupp_Meeker", "Support Points:", min = 1, max = 10, value = 4, step = 1),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("dsLower_Meeker", "Lower Bound:", value = 76, step = 1, width = "45%"),
        numericInput("dsUpper_Meeker", "Upper Bound:", value = 150, step = 1, width = "45%")
      ),
      
      h4("🧪 True Model Mean Response Parameters"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("p1_Meeker", withMathJax("\\(\\zeta_1\\):"), value = 14.75, step = 0.01, width = "45%"),
        numericInput("p2_Meeker", withMathJax("\\(\\zeta_2\\):"), value = -1.39, step = 0.01, width = "45%")
      ),
      
      h4("🎯 True Model Dispersion Parameter"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("disp_phi1_Meeker", "Dispersion (\\(\\phi_1\\)):", value = 10.97, step = 0.01, width = "45%"),
        numericInput("disp_phi2_Meeker", "Dispersion (\\(\\phi_2\\)):", value = -2.50, step = 0.01, width = "45%")
      ),
      
      tags$h4("⚔️ Rival Model Mean Response Parameters Bound"),
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
      
      h4("🎯 Rival Model Dispersion Parameter Bound"),
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
      
      h4("🚦 Censoring Threshold Setting"),
      numericInput("tc_Meeker", "Censoring Threshold (tc):", value = 1000, min = 1, step = 100),
      
      tags$h4("⚙️ Algorithm Settings"),
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
      
      div(
        style = "display: flex; justify-content: center; margin-top: 20px;",
        actionButton(
          "run_Meeker", 
          "🚀 Run Optimization", 
          style = "background-color:#FF5733; color:white; font-weight:bold; width: 80%; max-width: 300px;"
        )
      ),
      
      div(
        id = "loading_spinner_Meeker",
        style = "margin-top: 15px; display: none; text-align: center;",
        tags$div(
          style = "display: inline-block;",
          tags$div(
            class = "loader",  # CSS spinner
            style = "margin-bottom: 5px;"
          ),
          tags$div(
            span("Running Optimization...", style = "color: #FF5733; font-weight: bold; font-size: 15px;")
          )
        )
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Result",
                 withMathJax(
                   div(
                     style = "background-color: #fff8f0; border: 1px solid #ffc107; border-radius: 10px; padding: 15px; margin-bottom: 15px;
                box-shadow: 2px 2px 6px rgba(0,0,0,0.1); font-size:16px;",
                     HTML("
                      <b>🔍 Model Discrimination</b><br/>
                      <div style='display: flex; justify-content: space-between;'>
                        <div style='width: 48%;'>
                          <b>True Model:</b><br/>
                          $$\\zeta_1 + \\zeta_2 \\cdot \\log(x - 75.71)$$
                          <b>True Model Dispersion:</b><br/>
                          $$\\exp\\left\\{\\phi_1 + \\phi_2 \\cdot \\log(x - 75.71)\\right\\}$$
                       </div>
                       <div style='width: 48%;'>
                          <b>Rival Model:</b><br/>
                          $$\\delta_1 + \\delta_2 \\cdot \\log(x - 75.71)$$
                          <b>Rival Model Dispersion:</b><br/>
                          $$\\exp\\left\\{\\kappa_1 + \\kappa_2 \\cdot \\log(x - 75.71)\\right\\}$$
                       </div>
                       </div>
                      The goal of this design is to discriminate between the true model and a competing rival model by maximizing the divergence between them over the design space.<br/>
                      The general objective function is given by:
                      $$\\max_{\\xi\\in \\Xi} \\text{CKL}_{r,tr}(\\xi)=\\max_{\\xi\\in \\Xi} \\min_{\\substack{\\theta_r\\in \\Theta_r}}\\int_{X}D_{\\text{CKL}}(M_{tr},M_r,x,\\theta_r) \\xi(dx)$$
                      Formally, the optimal design criterion value is defined as:<br/>
                      $$C^* = \\max_{\\xi \\in \\Xi} \\min_{\\theta_r \\in \\Theta_r} \\left\\{ \\int_X D_{\\text{CKL}} \\left( M_{tr}(x, \\theta_{tr}), M_r(x, \\theta_r) \\right) \\xi(dx) \\right\\}$$
                    ")
                   )
                 ),
                 bsCollapse(id = "resultCollapse", multiple = TRUE,
                            bsCollapsePanel("📌 Support Points and Weights",
                                            tableOutput("bestDesign_Meeker"), style = "info"
                            ),
                            bsCollapsePanel("📊 Criterion Value and CPU Time",
                                            tableOutput("optimizationSummary_Meeker"), style = "info"
                            ),
                            bsCollapsePanel("📐 Model Parameters",
                                            fluidRow(
                                              column(6,
                                                     tags$h5("True Model Parameters"),
                                                     tableOutput("model1_params_Meeker")
                                              ),
                                              column(6,
                                                     tags$h5("Rival Model Parameters"),
                                                     tableOutput("model2_params_Meeker")
                                              )
                                            ),
                                            style = "success"
                            ),
                            bsCollapsePanel("📈 Directional Derivative Plot",
                                            plotOutput("derivativePlot_Meeker"), style = "primary"
                            )
                 )
        ),
        tabPanel("Setting",
                 div(
                   style = "background-color: #fff8f0;
                    border: 1px solid #ffc107;
                    border-radius: 10px;
                    padding: 15px;
                    margin-bottom: 15px;
                    box-shadow: 2px 2px 6px rgba(0,0,0,0.1);
                    font-size: 16px;",
                   h4("📋 Selected Settings:"),
                   tableOutput("selectedSettings_Meeker")
                 )
        )
      )
    )
  )
)