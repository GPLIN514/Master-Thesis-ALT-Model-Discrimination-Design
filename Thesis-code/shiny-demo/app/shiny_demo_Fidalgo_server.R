# Define Server
shiny_demo_Fidalgo_server <- function(input, output, session) {
  
  observeEvent(input$run_Fidalgo, {
    shinyjs::show("loading_spinner_Fidalgo")
    
    af1_mean_Fidalgo <- function(x, p) (p[1]*x)/(p[2] + x) + p[3]*x
    af2_mean_Fidalgo <- function(x, p) (p[1]*x)/(p[2] + x)
    af1_disp_Fidalgo <- function(x, p) rep(p[1], length(x))
    af2_disp_Fidalgo <- function(x, p) rep(p[1], length(x))
    
    # calculate options based on user selection
    if (input$calculate_option_Fidalgo == "close-form") {
      
      if (input$distribution_Fidalgo == "lognormal") {
        observed_func <- function(xt, xr, vt, vr) {
          mu_t <- log(xt) - 0.5 * log(1 + (vt / (xt^2)))
          mu_r <- log(xr) - 0.5 * log(1 + (vr / (xr^2)))
          var_t <- log(1 + (vt / (xt^2)))
          var_r <- log(1 + (vr / (xr^2)))
          
          log(sqrt(var_r)) - log(sqrt(var_t)) - ((var_r - var_t - (mu_r - mu_t)^2) / (2 * var_r))
        }
      } else if (input$distribution_Fidalgo == "weibull") {
        observed_func <- function(xt, xr, vt, vr) {
          gamma_const <- 0.57721566490153286060  # Euler-Mascheroni constant
          k1 <- 1/sqrt(vt)
          k2 <- 1/sqrt(vr)
          lambda1 <- exp(xt)
          lambda2 <- exp(xr)
          
          term1 <- log(k1 / k2)
          term2 <- log(lambda2 / lambda1)
          term3 <- -((k1 - 1) / k1) * gamma_const
          term4 <- -1
          term5 <- -(k2 - 1) * log(lambda1 / lambda2)
          term6 <- ((k2 - 1) / k1) * gamma_const
          term7 <- (lambda1 / lambda2)^k2 * gamma((k2 / k1) + 1)
          
          term1 + term2 + term3 + term4 + term5 + term6 + term7
        }
      }
      
    } else if (input$calculate_option_Fidalgo == "integration") {
      
      if (input$distribution_Fidalgo == "lognormal") {
        observed_func <- function(y, m1, m2, s1, s2) {
          mu_t <- log(m1) - 0.5 * log(1 + (s1 / m1^2))
          mu_r <- log(m2) - 0.5 * log(1 + (s2 / m2^2))
          var_t <- log(1 + (s1 / m1^2))
          var_r <- log(1 + (s2 / m2^2))
          
          lpdf1 <- dlnorm(y, mu_t, sqrt(var_t), log = TRUE)
          lpdf2 <- dlnorm(y, mu_r, sqrt(var_r), log = TRUE)
          pdf1 <- exp(lpdf1)
          pdf1 * (lpdf1 - lpdf2)
        }
      } else if (input$distribution_Fidalgo == "weibull") {
        observed_func <- function(y, m1, m2, s1, s2) {
          lpdf1 <- log(dweibull(y, 1/s1, exp(m1)) + 1e-12)
          lpdf2 <- log(dweibull(y, 1/s2, exp(m2)) + 1e-12)
          pdf1 <- exp(lpdf1)
          pdf1 * (lpdf1 - lpdf2)
        }
      }
    }
    
    Fidalgo_div_dynamic <- function(xt, xr, st, sr) {
      kl_divs <- rep(0, length(xt))
      
      for (i in 1:length(xt)) {
        if (input$calculate_option_Fidalgo == "close-form") {
          kl_divs[i] <- observed_func(xt[i], xr[i], st[i], sr[i])
          
        } else if (input$calculate_option_Fidalgo == "integration") {
          
          if (input$distribution_Fidalgo == "lognormal") {
            upper <- 10000
          } else if (input$distribution_Fidalgo == "weibull") {
            upper <- qweibull(0.99, shape = 1/st[i], scale = exp(xt[i]))
          }
          
          kl_divs[i] <- integrate(observed_func, 0, upper,
                                  m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                                  subdivisions = 1000, stop.on.error = FALSE)$value
        }
      }
      
      return(kl_divs)
    }
    
    
    # Model parameters from input
    model_info_Fidalgo <- list(
      list(
        mean = af1_mean_Fidalgo, disp = af1_disp_Fidalgo,
        meanPara = c(input$p1_Fidalgo, input$p2_Fidalgo, input$p3_Fidalgo),
        dispPara = input$disp_Fidalgo
      ),
      list(
        mean = af2_mean_Fidalgo, disp = af2_disp_Fidalgo,
        meanParaLower = c(input$lower_delta_1_Fidalgo, input$lower_delta_2_Fidalgo),
        meanParaUpper = c(input$upper_delta_1_Fidalgo, input$upper_delta_2_Fidalgo),
        dispParaLower = input$disp_Fidalgo,
        dispParaUpper = input$disp_Fidalgo
      )
    )
    
    # PSO and BFGS options
    if (input$use_inner_pso_Fidalgo) {
      PSO_INFO <- getPSOInfo(nSwarm = c(input$nSwarm_outer_Fidalgo, input$nSwarm_inner_Fidalgo), 
                             maxIter = c(input$maxIter_outer_Fidalgo, input$maxIter_inner_Fidalgo))
      LBFGS_INFO <- getLBFGSInfo(IF_INNER_LBFGS = FALSE)
    } else {
      PSO_INFO <- getPSOInfo(nSwarm = input$nSwarm_outer_Fidalgo, maxIter = input$maxIter_outer_Fidalgo)
      LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = input$LBFGS_RETRY_Fidalgo)
    }
    
    # Run optimization
    af_res_Fidalgo <- DiscrimOD(MODEL_INFO = model_info_Fidalgo, DISTANCE = Fidalgo_div_dynamic,
                                nSupp = input$nSupp_Fidalgo, dsLower = input$dsLower_Fidalgo, dsUpper = input$dsUpper_Fidalgo,
                                crit_type = "pair_fixed_true", PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                                seed = 100, verbose = TRUE)
    
    # Store results
    best_design_Fidalgo <- af_res_Fidalgo$BESTDESIGN
    best_val_Fidalgo <- af_res_Fidalgo$BESTVAL
    cpu_time_Fidalgo <- af_res_Fidalgo$CPUTIME
    
    # Equivalence test
    equivalence_res_Fidalgo <- equivalence(
      ngrid = 100,
      PSO_RESULT = af_res_Fidalgo,
      MODEL_INFO = model_info_Fidalgo,
      DISTANCE = Fidalgo_div_dynamic,
      dsLower = input$dsLower_Fidalgo,
      dsUpper = input$dsUpper_Fidalgo,
      crit_type = "pair_fixed_true",
      PSO_INFO = PSO_INFO,
      LBFGS_INFO = LBFGS_INFO
    )
    
    # Calculate design criterion for the optimal design
    design_criterion_result_Fidalgo <- designCriterion(
      af_res_Fidalgo$BESTDESIGN,
      MODEL_INFO = model_info_Fidalgo,
      DISTANCE = Fidalgo_div_dynamic,
      dsLower = input$dsLower_Fidalgo,
      dsUpper = input$dsUpper_Fidalgo,
      crit_type = "pair_fixed_true",
      MaxMinStdVals = NULL,
      PSO_INFO = PSO_INFO,
      LBFGS_INFO = LBFGS_INFO
    )
    
    # Extract only the theta2 part from the result
    theta2_Fidalgo <- design_criterion_result_Fidalgo$theta2
    
    # Update UI outputs
    output$bestDesign_Fidalgo <- renderTable({
      colnames(best_design_Fidalgo) <- c("Support Point", "Weight")
      round(best_design_Fidalgo, 3)
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c", rownames = FALSE)
    
    output$designCriterionResult_Fidalgo <- renderTable({
      data.frame(
        Parameter = c("Design Criterion"),
        Value = c(round(design_criterion_result_Fidalgo, 6))
      )
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    # Model 1 Parameters
    output$model1_params_Fidalgo <- renderTable({
      model1_df <- data.frame(
        Zeta1 = round(theta2_Fidalgo[1, 1], 4),
        Zeta2 = round(theta2_Fidalgo[1, 2], 4),
        Zeta3 = round(theta2_Fidalgo[1, 3], 4),
        Sigma1 = round(theta2_Fidalgo[1, 4], 4)
      )
      colnames(model1_df) <- c("ζ₁", "ζ₂", "ζ₃", "σ")
      model1_df
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    # Model 2 Parameters
    output$model2_params_Fidalgo <- renderTable({
      model2_df <- data.frame(
        Delta1 = round(theta2_Fidalgo[2, 1], 4),
        Delta2 = round(theta2_Fidalgo[2, 2], 4),
        Sigma2 = round(theta2_Fidalgo[2, 3], 4)
      )
      colnames(model2_df) <- c("δ₁", "δ₂", "σ") 
      model2_df
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    output$optimizationSummary_Fidalgo <- renderTable({
      data.frame(
        Criterion_Value = format(round(best_val_Fidalgo, 6), nsmall = 6), 
        CPU_Time = round(cpu_time_Fidalgo, 3) 
      )
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    output$derivativePlot_Fidalgo <- renderPlot({
      plot(equivalence_res_Fidalgo$Grid_1, equivalence_res_Fidalgo$DirDeriv, type = "l", col = "blue",
           main = "Directional Derivative", xlab = "x", ylab = "Directional Derivative")
      abline(h = 0)
      points(best_design_Fidalgo[,1], rep(0, nrow(best_design_Fidalgo)), pch = 16, col = "red")
    })
    output$selectedSettings_Fidalgo <- renderTable({
      settings <- data.frame(
        Setting = c(
          "Outer PSO Particles",
          "Outer PSO Iterations",
          "Use Inner PSO",
          "Inner PSO Particles",
          "Inner PSO Iterations",
          "BFGS Retries",
          "Support Points",
          "Design Bounds",
          "True Model Parameters",
          "Model Dispersion"
        ),
        Value = c(
          input$nSwarm_outer_Fidalgo,
          input$maxIter_outer_Fidalgo,
          ifelse(input$use_inner_pso_Fidalgo, "TRUE", "FALSE"),
          ifelse(input$use_inner_pso_Fidalgo, input$nSwarm_inner_Fidalgo, "N/A"),
          ifelse(input$use_inner_pso_Fidalgo, input$maxIter_inner_Fidalgo, "N/A"),
          ifelse(!input$use_inner_pso_Fidalgo, input$LBFGS_RETRY_Fidalgo, "N/A"),
          input$nSupp_Fidalgo,
          paste0("[", input$dsLower_Fidalgo, ", ", input$dsUpper_Fidalgo, "]"),
          paste0("(", input$p1_Fidalgo, ", ", input$p2_Fidalgo, ", ", input$p3_Fidalgo, ")"),
          paste0(input$disp_Fidalgo)
        )
      )
      settings
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    shinyjs::hide("loading_spinner_Fidalgo")
  })
}