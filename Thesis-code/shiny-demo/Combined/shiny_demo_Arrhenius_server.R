# Define Server
shiny_demo_Arrhenius_server <- function(input, output, session) {
  
  observeEvent(input$run_Arrhenius, {
    
    # Fixed model structures
    af1_mean_Arrhenius <- function(x, p) p[1] + p[2] * (11605/(x+273.15)) + p[3] * (11605/(x+273.15))^2
    af2_mean_Arrhenius <- function(x, p) p[1] + p[2] * (11605/(x+273.15))
    af1_disp_Arrhenius <- function(x, p) rep(p[1], length(x))
    af2_disp_Arrhenius <- function(x, p) rep(p[1], length(x))
    
    # Divergence functions based on user selection
    if (input$distribution_Arrhenius == "lognormal") {
      
      if (input$divergence_Arrhenius == "kl") {
        observed_func <- function(y, m1, m2, s1, s2) {
          lpdf1 <- dlnorm(y, m1, s1, log = TRUE)
          lpdf2 <- dlnorm(y, m2, s2, log = TRUE)
          pdf1 <- exp(lpdf1)
          
          val <- pdf1*(lpdf1 - lpdf2)
          return(val)
        }
        censored_func <- function(y, m1, m2, s1, s2) {
          lcdf1 <- log(1 - plnorm(y, m1, s1) + 1e-12)
          lcdf2 <- log(1 - plnorm(y, m2, s2) + 1e-12)
          cdf1 <- exp(lcdf1)
          #
          val <- cdf1*(lcdf1 - lcdf2)
          return(val)
        }
        
      } else if (input$divergence_Arrhenius == "lw") {
        observed_func <- function(y, m1, m2, s1, s2) {
          pdf1 <- dlnorm(y, m1, s1)
          pdf2 <- dlnorm(y, m2, s2)
          
          val <- pdf1*(log(2*pdf1) - log(pdf1+pdf2))
          return(val)
        }
        censored_func <- function(y, m1, m2, s1, s2) {
          cdf1 <- (1 - plnorm(y, m1, s1) + 1e-12)
          cdf2 <- (1 - plnorm(y, m2, s2) + 1e-12)
          #
          val <- cdf1*(log(2*cdf1) - log(cdf1+cdf2))
        }
        
      } else if (input$divergence_Arrhenius == "b") {
        observed_func <- function(y, m1, m2, s1, s2) {
          pdf1 <- dlnorm(y, m1, s1)
          pdf2 <- dlnorm(y, m2, s2)
          
          val <- sqrt(pdf1*pdf2)
          return(val)
        }
        censored_func <- function(y, m1, m2, s1, s2) {
          cdf1 <- (1 - plnorm(y, m1, s1) + 1e-12)
          cdf2 <- (1 - plnorm(y, m2, s2) + 1e-12)
          #
          val <- sqrt(cdf1*cdf2)
          return(val)
        }
        
      } else if (input$divergence_Arrhenius == "chi") {
        observed_func <- function(y, m1, m2, s1, s2) {
          pdf1 <- dlnorm(y, m1, s1)
          pdf2 <- dlnorm(y, m2, s2)
          
          val <- ((pdf1^2)/(pdf2+1e-12))
          return(val)
        }
        censored_func <- function(y, m1, m2, s1, s2) {
          cdf1 <- (1 - plnorm(y, m1, s1) + 1e-12)
          cdf2 <- (1 - plnorm(y, m2, s2) + 1e-12)
          #
          val <- ((cdf1^2/cdf2))
          return(val)
        }
      }
      
    } else if (input$distribution_Arrhenius == "weibull") {
      
      if (input$divergence_Arrhenius == "kl") {
        observed_func <- function(y, m1, m2, s1, s2) {
          lpdf1 <- dweibull(y, 1/s1, exp(m1), log = TRUE)
          lpdf2 <- dweibull(y, 1/s2, exp(m2), log = TRUE)
          pdf1 <- exp(lpdf1)
          
          val <- pdf1*(lpdf1 - lpdf2)
          return(val)
        }
        censored_func <- function(y, m1, m2, s1, s2) {
          lcdf1 <- log(1 - pweibull(y, 1/s1, exp(m1)) + 1e-12)
          lcdf2 <- log(1 - pweibull(y, 1/s2, exp(m2)) + 1e-12)
          cdf1 <- exp(lcdf1)
          #
          val <- cdf1*(lcdf1 - lcdf2)
          return(val)
        }
        
      } else if (input$divergence_Arrhenius == "lw") {
        observed_func <- function(y, m1, m2, s1, s2) {
          pdf1 <- dweibull(y, 1/s1, exp(m1))
          pdf2 <- dweibull(y, 1/s2, exp(m2))
          
          val <- pdf1*(log(2*pdf1) - log(pdf1+pdf2))
          return(val)
        }
        censored_func <- function(y, m1, m2, s1, s2) {
          cdf1 <- (1 - pweibull(y, 1/s1, exp(m1)) + 1e-12)
          cdf2 <- (1 - pweibull(y, 1/s2, exp(m2)) + 1e-12)
          #
          val <- cdf1*(log(2*cdf1) - log(cdf1+cdf2))
          return(val)
        }
        
      } else if (input$divergence_Arrhenius == "b") {
        observed_func <- function(y, m1, m2, s1, s2) {
          pdf1 <- dweibull(y, 1/s1, exp(m1))
          pdf2 <- dweibull(y, 1/s2, exp(m2))
          
          val <- sqrt(pdf1*pdf2)
          return(val)
        }
        censored_func <- function(y, m1, m2, s1, s2) {
          cdf1 <- (1 - pweibull(y, 1/s1, exp(m1)) + 1e-12)
          cdf2 <- (1 - pweibull(y, 1/s2, exp(m2)) + 1e-12)
          #
          val <- sqrt(cdf1*cdf2)
          return(val)
        }
        
      } else if (input$divergence_Arrhenius == "chi") {
        observed_func <- function(y, m1, m2, s1, s2) {
          pdf1 <- dweibull(y, 1/s1, exp(m1))
          pdf2 <- dweibull(y, 1/s2, exp(m2))
          
          val <- ((pdf1^2)/pdf2)
          return(val)
        }
        censored_func <- function(y, m1, m2, s1, s2) {
          cdf1 <- (1 - pweibull(y, 1/s1, exp(m1)) + 1e-12)
          cdf2 <- (1 - pweibull(y, 1/s2, exp(m2)) + 1e-12)
          #
          val <- ((cdf1^2/cdf2))
          return(val)
        }
      }
    }
    
    Arrhenius_censored_div_dynamic <- function(xt, xr, st, sr, tc) {
      intVec <- rep(0, length(xt))
      for (i in 1:length(xt)) {
        intg_part <- integrate(observed_func, 0, tc,
                               m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                               subdivisions = 100, stop.on.error = FALSE)$value
        cens_part <- censored_func(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
        intVec[i] <- intg_part + cens_part
      }
      intVec
    }
    
    # Wrapper for fixed tc
    Arrhenius_censored_div_wrapper <- function(xt, xr, st, sr) {
      Arrhenius_censored_div_dynamic(xt, xr, st, sr, tc = input$tc_Arrhenius)
    }
    
    # Model parameters from input
    model_info_Arrhenius <- list(
      list(
        mean = af1_mean_Arrhenius, disp = af1_disp_Arrhenius,
        meanPara = c(input$p1_Arrhenius, input$p2_Arrhenius, input$p3_Arrhenius),
        dispPara = input$disp1_Arrhenius
      ),
      list(
        mean = af2_mean_Arrhenius, disp = af2_disp_Arrhenius,
        meanParaLower = c(input$lower_delta_1_Arrhenius, input$lower_delta_2_Arrhenius),
        meanParaUpper = c(input$upper_delta_1_Arrhenius, input$upper_delta_2_Arrhenius),
        dispParaLower = input$lower_disp2_Arrhenius,
        dispParaUpper = input$upper_disp2_Arrhenius
      )
    )
    
    # PSO and BFGS options
    if (input$use_inner_pso_Arrhenius) {
      PSO_INFO <- getPSOInfo(nSwarm = c(input$nSwarm_outer_Arrhenius, input$nSwarm_inner_Arrhenius), 
                             maxIter = c(input$maxIter_outer_Arrhenius, input$maxIter_inner_Arrhenius))
      LBFGS_INFO <- getLBFGSInfo(IF_INNER_LBFGS = FALSE)
    } else {
      PSO_INFO <- getPSOInfo(nSwarm = input$nSwarm_outer_Arrhenius, maxIter = input$maxIter_outer_Arrhenius)
      LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = input$LBFGS_RETRY_Arrhenius)
    }
    
    # Run optimization
    af_res_Arrhenius <- DiscrimOD(MODEL_INFO = model_info_Arrhenius, DISTANCE = Arrhenius_censored_div_wrapper,
                                  nSupp = input$nSupp_Arrhenius, dsLower = input$dsLower_Arrhenius, dsUpper = input$dsUpper_Arrhenius,
                                  crit_type = "pair_fixed_true", PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                                  seed = 100, verbose = TRUE)
    
    # Store results
    best_design_Arrhenius <- af_res_Arrhenius$BESTDESIGN
    best_val_Arrhenius <- af_res_Arrhenius$BESTVAL
    cpu_time_Arrhenius <- af_res_Arrhenius$CPUTIME
    
    # Equivalence test
    equivalence_res_Arrhenius <- equivalence(
      ngrid = 100,
      PSO_RESULT = af_res_Arrhenius,
      MODEL_INFO = model_info_Arrhenius,
      DISTANCE = Arrhenius_censored_div_wrapper,
      dsLower = input$dsLower_Arrhenius,
      dsUpper = input$dsUpper_Arrhenius,
      crit_type = "pair_fixed_true",
      PSO_INFO = PSO_INFO,
      LBFGS_INFO = LBFGS_INFO
    )
    
    # Calculate design criterion for the optimal design
    design_criterion_result_Arrhenius <- designCriterion(
      af_res_Arrhenius$BESTDESIGN,
      MODEL_INFO = model_info_Arrhenius,
      DISTANCE = Arrhenius_censored_div_wrapper,
      dsLower = input$dsLower_Arrhenius,
      dsUpper = input$dsUpper_Arrhenius,
      crit_type = "pair_fixed_true",
      MaxMinStdVals = NULL,
      PSO_INFO = PSO_INFO,
      LBFGS_INFO = LBFGS_INFO
    )
    
    # Extract only the theta2 part from the result
    theta2_Arrhenius <- design_criterion_result_Arrhenius$theta2
    
    # Update UI outputs
    output$bestDesign_Arrhenius <- renderTable({
      round(best_design_Arrhenius, 3)
    })
    
    output$designCriterionResult_Arrhenius <- renderTable({
      data.frame(
        Parameter = c("Design Criterion"),
        Value = c(round(design_criterion_result_Arrhenius, 6))
      )
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    # Model 1 Parameters
    output$model1_params_Arrhenius <- renderTable({
      model1_df <- data.frame(
        Zeta1 = round(theta2_Arrhenius[1, 1], 4),
        Zeta2 = round(theta2_Arrhenius[1, 2], 4),
        Zeta3 = round(theta2_Arrhenius[1, 3], 4),
        Sigma1 = round(theta2_Arrhenius[1, 4], 4)
      )
      colnames(model1_df) <- c("ζ₁", "ζ₂", "ζ₃", "σ₁")
      model1_df
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    # Model 2 Parameters
    output$model2_params_Arrhenius <- renderTable({
      model2_df <- data.frame(
        Delta1 = round(theta2_Arrhenius[2, 1], 4),
        Delta2 = round(theta2_Arrhenius[2, 2], 4),
        Sigma2 = round(theta2_Arrhenius[2, 3], 4)
      )
      colnames(model2_df) <- c("δ₁", "δ₂", "σ₂") 
      model2_df
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    output$optimizationSummary_Arrhenius <- renderTable({
      data.frame(
        Criterion_Value = format(round(best_val_Arrhenius, 6), nsmall = 6), 
        CPU_Time = round(cpu_time_Arrhenius, 3) 
      )
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    output$derivativePlot_Arrhenius <- renderPlot({
      plot(equivalence_res_Arrhenius$Grid_1, equivalence_res_Arrhenius$DirDeriv, type = "l", col = "blue",
           main = "Directional Derivative", xlab = "x", ylab = "Directional Derivative")
      abline(h = 0)
      points(best_design_Arrhenius[,1], rep(0, nrow(best_design_Arrhenius)), pch = 16, col = "red")
      output$selectedSettings_Arrhenius <- renderTable({
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
            "Model 1 Parameters",
            "Model 1 Dispersion Parameter",
            "Model 2 first Parameter Bound",
            "Model 2 second Parameter Bound",
            "Model 2 Dispersion Parameter Bound"
          ),
          Value = c(
            input$nSwarm_outer_Arrhenius,
            input$maxIter_outer_Arrhenius,
            ifelse(input$use_inner_pso_Arrhenius, "TRUE", "FALSE"),
            ifelse(input$use_inner_pso_Arrhenius, input$nSwarm_inner_Arrhenius, "N/A"),
            ifelse(input$use_inner_pso_Arrhenius, input$maxIter_inner_Arrhenius, "N/A"),
            ifelse(!input$use_inner_pso_Arrhenius, input$LBFGS_RETRY_Arrhenius, "N/A"),
            input$nSupp_Arrhenius,
            paste0("[", input$dsLower_Arrhenius, ", ", input$dsUpper_Arrhenius, "]"),
            paste0("(", input$p1_Arrhenius, ", ", input$p2_Arrhenius, ", ", input$p3_Arrhenius, ")"),
            paste0(input$disp1_Arrhenius),
            paste0("[", input$lower_delta_1_Arrhenius, ", ", input$upper_delta_1_Arrhenius, "]"),
            paste0("[", input$lower_delta_2_Arrhenius, ", ", input$upper_delta_2_Arrhenius, "]"),
            paste0("[", input$lower_disp2_Arrhenius, ", ", input$upper_disp2_Arrhenius, "]")
          )
        )
        settings
      }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    })
  })
}