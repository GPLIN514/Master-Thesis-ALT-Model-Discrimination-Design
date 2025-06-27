# Define UI
shiny_demo_Fidalgo_ui <- fluidPage(
  useShinyjs(),
  theme = shinytheme("united"),
  
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
  
  titlePanel(HTML("MMM vs. MM under pharmacokinetic model (López-Fidalgo et al., 2007)")),
  
  sidebarLayout(
    sidebarPanel(
      style = "background-color: #fffef3;
         border: 1px solid #ffc107;
         border-radius: 10px;
         padding: 15px;
         margin-bottom: 15px;
         box-shadow: 2px 2px 6px rgba(0,0,0,0.1);
         font-size: 16px;",
      
      tags$h4("🧠 Select Calculate Options"),
      selectInput("calculate_option_Fidalgo", "Divergence Method:",
                  choices = c("Closed-Form" = "close-form",
                              "Integration" = "integration"),
                  selected = "close-form"),
      
      selectInput("distribution_Fidalgo", "Select Distribution:",
                  choices = c("Log-Normal" = "lognormal", "Weibull" = "weibull"),
                  selected = "lognormal"),
      
      tags$h4("🧭 Design Space Bound Settings"),
      sliderInput("nSupp_Fidalgo", "Support Points:", min = 1, max = 10, value = 3, step = 1),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("dsLower_Fidalgo", "Lower Bound:", value = 0.1, step = 0.1, width = "45%"),
        numericInput("dsUpper_Fidalgo", "Upper Bound:", value = 5, step = 0.1, width = "45%")
      ),
      
      h4("🧪 True Model Mean Response Parameters"),
      tags$div(
        style = "display: flex; justify-content: space-between;",
        numericInput("p1_Fidalgo", withMathJax("\\(\\zeta_1\\):"), value = 1, step = 0.1, width = "30%"),
        numericInput("p2_Fidalgo", withMathJax("\\(\\zeta_2\\):"), value = 1, step = 0.1, width = "30%"),
        numericInput("p3_Fidalgo", withMathJax("\\(\\zeta_3\\):"), value = 1, step = 0.1, width = "30%")
      ),
      
      h4("🎯 Dispersion Parameter"),
      numericInput("disp_Fidalgo", "Dispersion (\\(\\sigma\\)):", value = 1, step = 0.01),
      
      tags$h4("⚔️ Rival Model Mean Response Parameters Bound"),
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
      
      tags$h4("⚙️ Algorithm Settings"),
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
      
      div(
        style = "display: flex; justify-content: center; margin-top: 20px;",
        actionButton(
          "run_Fidalgo", 
          "🚀 Run Optimization", 
          style = "background-color:#FF5733; color:white; font-weight:bold; width: 80%; max-width: 300px;"
        )
      ),
      
      div(
        id = "loading_spinner_Fidalgo",
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
                          $$\\frac{\\zeta_1 \\cdot x} {\\zeta_2 + x} + \\zeta_3\\ x$$
                       </div>
                       <div style='width: 48%;'>
                          <b>Rival Model:</b><br/>
                          $$\\frac{\\delta_1 \\cdot x} {\\delta_2 + x}$$
                       </div>
                       </div>
                      The goal of this design is to discriminate between the true model and a competing rival model by maximizing the divergence between them over the design space.<br/>
                      The general objective function is given by:
                      $$\\max_{\\xi\\in \\Xi} \\text{KL}_{r,tr}(\\xi)=\\max_{\\xi\\in \\Xi} \\min_{\\substack{\\theta_r\\in \\Theta_r}}\\int_{X}D_{\\text{KL}}(M_{tr},M_r,x,\\theta_r) \\xi(dx)$$
                      Formally, the optimal design criterion value is defined as:<br/>
                      $$C^* = \\max_{\\xi \\in \\Xi} \\min_{\\theta_r \\in \\Theta_r} \\left\\{ \\int_X D_{\\text{KL}}\\left( M_{tr}(x, \\theta_{tr}), M_r(x, \\theta_r) \\right) \\xi(dx) \\right\\}$$
                    ")
                   )
                 ),
                 bsCollapse(id = "resultCollapse", multiple = TRUE,
                            bsCollapsePanel("📌 Approximation Design",
                                            tableOutput("bestDesign_Fidalgo"), style = "info"
                            ),
                            bsCollapsePanel("📊 Criterion Value and CPU Time",
                                            tableOutput("optimizationSummary_Fidalgo"), style = "info"
                            ),
                            bsCollapsePanel("📐 Model Parameters",
                                            fluidRow(
                                              column(6,
                                                     tags$h5("True Model Parameters"),
                                                     tableOutput("model1_params_Fidalgo")
                                              ),
                                              column(6,
                                                     tags$h5("Rival Model Parameters"),
                                                     tableOutput("model2_params_Fidalgo")
                                              )
                                            ),
                                            style = "success"
                            ),
                            bsCollapsePanel("📈 Directional Derivative Plot",
                                            plotOutput("derivativePlot_Fidalgo"), style = "primary"
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
                   tableOutput("selectedSettings_Fidalgo")
                 )
        )
      )
    )
  )
)