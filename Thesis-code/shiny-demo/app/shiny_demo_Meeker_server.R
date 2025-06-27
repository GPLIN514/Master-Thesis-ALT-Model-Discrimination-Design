# Define Server
shiny_demo_Meeker_server <- function(input, output, session) {
  
  observe({
    updateSelectInput(session, "distribution_M2_Meeker",
                      selected = if (input$distribution_M1_Meeker == "lognormal") "weibull" else "lognormal")
    
    disable("distribution_M2_Meeker")
  })
  
  observeEvent(input$run_Meeker, {
    shinyjs::show("loading_spinner_Meeker")
    
    # Fixed model structures
    af1_mean_Meeker <- function(x, p) p[1] + p[2] * (11605/(x+273.15)) + p[3] * (11605/(x+273.15))^2
    af2_mean_Meeker <- function(x, p) p[1] + p[2] * (11605/(x+273.15))
    af1_disp_Meeker <- function(x, p) exp(p[1] + p[2]*log(x + 1e-12))
    af2_disp_Meeker <- function(x, p) exp(p[1] + p[2]*log(x + 1e-12))
    
    # Divergence functions based on user selection
    if (input$distribution_M1_Meeker == "lognormal") {
      observed_func <- function(y, m1, m2, s1, s2) {
        if (!all(is.finite(c(y, m1, m2, s1, s2))) || any(c(s1, s2) <= 0)) return(1e6)
        tryCatch({
          lpdf1 <- log(dlnorm(y, m1, s1) + 1e-12)
          lpdf2 <- log(dweibull(y, 1/s2, exp(m2)) + 1e-12)
          pdf1 <- exp(lpdf1)
          pdf1 * (lpdf1 - lpdf2)
        }, error = function(e) return(1e6))
      }
      
      censored_func <- function(y, m1, m2, s1, s2) {
        #
        lcdf1 <- log(1 - plnorm(y, m1, s1) + 1e-12)
        lcdf2 <- log(1 - pweibull(y, 1/s2, exp(m2)) + 1e-12)
        cdf1 <- exp(lcdf1)
        #
        val <- cdf1*(lcdf1 - lcdf2)
        return(val)
      }
      
    } else if (input$distribution_M1_Meeker == "weibull") {
      observed_func <- function(y, m1, m2, s1, s2) {
        lpdf1 <- log(dweibull(y, 1/s1, exp(m1)) + 1e-12)
        lpdf2 <- log(dlnorm(y, m2, s2) + 1e-12)
        pdf1 <- exp(lpdf1)
        
        val <- pdf1*(lpdf1 - lpdf2)
        return(val)
      }
      
      censored_func <- function(y, m1, m2, s1, s2) {
        lcdf1 <- log(1 - pweibull(y, 1/s1, exp(m1)) + 1e-12)
        lcdf2 <- log(1 - plnorm(y, m2, s2) + 1e-12)
        cdf1 <- exp(lcdf1)
        
        val <- cdf1*(lcdf1 - lcdf2)
        return(val)
      }
    }
    
    Meeker_censored_div_dynamic <- function(xt, xr, st, sr, tc) {
      intVec <- numeric(length(xt))
      
      for (i in seq_along(xt)) {
        m1 <- xt[i]
        m2 <- xr[i]
        s1 <- st[i]
        s2 <- sr[i]
        
        if (!all(is.finite(c(m1, m2, s1, s2, tc))) || any(c(s1, s2) <= 0)) {
          intVec[i] <- NA
          next
        }
        
        intg_part <- tryCatch({
          res <- integrate(function(y) {
            val <- tryCatch({
              out <- observed_func(y, m1, m2, s1, s2)
              if (length(out) != 1 || !is.finite(out)) stop("Invalid value")
              out
            }, error = function(e) {
              warning(sprintf("Observed function error at i=%d, y=%.4f: %s", i, y, e$message))
              return(1e6)
            })
            return(out)
          }, lower = 0, upper = tc, subdivisions = 100, stop.on.error = FALSE)
          res$value
        }, error = function(e) {
          warning(sprintf("Integration failed at i=%d: %s", i, e$message))
          return(1e6)
        })
        
        cens_part <- tryCatch({
          val <- censored_func(tc, m1, m2, s1, s2)
          if (!is.finite(val)) return(1e6)
          val
        }, error = function(e) {
          warning(sprintf("Censored function failed at i=%d: %s", i, e$message))
          return(1e6)
        })
        
        intVec[i] <- intg_part + cens_part
      }
      
      return(intVec)
    }
    
    # Wrapper for fixed tc
    Meeker_censored_div_wrapper <- function(xt, xr, st, sr) {
      Meeker_censored_div_dynamic(xt, xr, st, sr, tc = input$tc_Meeker)
    }
    
    # Model parameters from input
    model_info_Meeker <- list(
      list(
        mean = af1_mean_Meeker, disp = af1_disp_Meeker,
        meanPara = c(input$p1_Meeker, input$p2_Meeker),
        dispPara = c(input$disp_phi1_Meeker,input$disp_phi2_Meeker)
      ),
      list(
        mean = af2_mean_Meeker, disp = af2_disp_Meeker,
        meanParaLower = c(input$lower_delta_1_Meeker, input$lower_delta_2_Meeker),
        meanParaUpper = c(input$upper_delta_1_Meeker, input$upper_delta_2_Meeker),
        dispParaLower = c(input$lower_disp_kappa1_Meeker,input$lower_disp_kappa2_Meeker),
        dispParaUpper = c(input$upper_disp_kappa1_Meeker,input$upper_disp_kappa2_Meeker)
      )
    )
    
    # PSO and BFGS options
    if (input$use_inner_pso_Meeker) {
      PSO_INFO <- getPSOInfo(nSwarm = c(input$nSwarm_outer_Meeker, input$nSwarm_inner_Meeker), 
                             maxIter = c(input$maxIter_outer_Meeker, input$maxIter_inner_Meeker))
      LBFGS_INFO <- getLBFGSInfo(IF_INNER_LBFGS = FALSE)
    } else {
      PSO_INFO <- getPSOInfo(nSwarm = input$nSwarm_outer_Meeker, maxIter = input$maxIter_outer_Meeker)
      LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = input$LBFGS_RETRY_Meeker)
    }
    
    # Run optimization
    af_res_Meeker <- DiscrimOD(MODEL_INFO = model_info_Meeker, DISTANCE = Meeker_censored_div_wrapper,
                               nSupp = input$nSupp_Meeker, dsLower = input$dsLower_Meeker, dsUpper = input$dsUpper_Meeker,
                               crit_type = "pair_fixed_true", PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                               seed = 100, verbose = TRUE)
    
    # Store results
    best_design_Meeker <- af_res_Meeker$BESTDESIGN
    best_val_Meeker <- af_res_Meeker$BESTVAL
    cpu_time_Meeker <- af_res_Meeker$CPUTIME
    
    # Equivalence test
    equivalence_res_Meeker <- equivalence(
      ngrid = 100,
      PSO_RESULT = af_res_Meeker,
      MODEL_INFO = model_info_Meeker,
      DISTANCE = Meeker_censored_div_wrapper,
      dsLower = input$dsLower_Meeker,
      dsUpper = input$dsUpper_Meeker,
      crit_type = "pair_fixed_true",
      PSO_INFO = PSO_INFO,
      LBFGS_INFO = LBFGS_INFO
    )
    
    # Calculate design criterion for the optimal design
    design_criterion_result_Meeker <- designCriterion(
      af_res_Meeker$BESTDESIGN,
      MODEL_INFO = model_info_Meeker,
      DISTANCE = Meeker_censored_div_wrapper,
      dsLower = input$dsLower_Meeker,
      dsUpper = input$dsUpper_Meeker,
      crit_type = "pair_fixed_true",
      MaxMinStdVals = NULL,
      PSO_INFO = PSO_INFO,
      LBFGS_INFO = LBFGS_INFO
    )
    
    # Extract only the theta2 part from the result
    theta2_Meeker <- design_criterion_result_Meeker$theta2
    
    # Update UI outputs
    output$bestDesign_Meeker <- renderTable({
      colnames(best_design_Meeker) <- c("Support Point", "Weight")
      round(best_design_Meeker, 3)
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c", rownames = FALSE)
    
    output$designCriterionResult_Meeker <- renderTable({
      data.frame(
        Parameter = c("Design Criterion"),
        Value = c(round(design_criterion_result_Meeker, 6))
      )
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    # Model 1 Parameters
    output$model1_params_Meeker <- renderTable({
      model1_df <- data.frame(
        Zeta1 = round(theta2_Meeker[1, 1], 4),
        Zeta2 = round(theta2_Meeker[1, 2], 4),
        Zeta3 = round(theta2_Meeker[1, 3], 4),
        Sigma1 = round(theta2_Meeker[1, 4], 4)
      )
      colnames(model1_df) <- c("ζ₁", "ζ₂", "ζ₃", "σ₁")
      model1_df
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    # Model 2 Parameters
    output$model2_params_Meeker <- renderTable({
      model2_df <- data.frame(
        Delta1 = round(theta2_Meeker[2, 1], 4),
        Delta2 = round(theta2_Meeker[2, 2], 4),
        Sigma2 = round(theta2_Meeker[2, 3], 4)
      )
      colnames(model2_df) <- c("δ₁", "δ₂", "σ₂") 
      model2_df
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    output$optimizationSummary_Meeker <- renderTable({
      data.frame(
        Criterion_Value = format(round(best_val_Meeker, 6), nsmall = 6), 
        CPU_Time = round(cpu_time_Meeker, 3) 
      )
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    output$derivativePlot_Meeker <- renderPlot({
      plot(equivalence_res_Meeker$Grid_1, equivalence_res_Meeker$DirDeriv, type = "l", col = "blue",
           main = "Directional Derivative", xlab = "x", ylab = "Directional Derivative")
      abline(h = 0)
      points(best_design_Meeker[,1], rep(0, nrow(best_design_Meeker)), pch = 16, col = "red")
    })
    output$selectedSettings_Meeker <- renderTable({
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
          "True Model Dispersion Parameter",
          "Rival Model first Parameter Bound",
          "Rival Model second Parameter Bound",
          "Rival Model Dispersion first Parameter Bound",
          "Rival Model Dispersion second Parameter Bound"
        ),
        Value = c(
          input$nSwarm_outer_Meeker,
          input$maxIter_outer_Meeker,
          ifelse(input$use_inner_pso_Meeker, "TRUE", "FALSE"),
          ifelse(input$use_inner_pso_Meeker, input$nSwarm_inner_Meeker, "N/A"),
          ifelse(input$use_inner_pso_Meeker, input$maxIter_inner_Meeker, "N/A"),
          ifelse(!input$use_inner_pso_Meeker, input$LBFGS_RETRY_Meeker, "N/A"),
          input$nSupp_Meeker,
          paste0("[", input$dsLower_Meeker, ", ", input$dsUpper_Meeker, "]"),
          paste0("(", input$p1_Meeker, ", ", input$p2_Meeker, ", ", input$p3_Meeker, ")"),
          paste0("(", input$disp_phi1_Meeker,", ", input$disp_phi2_Meeker, ")"),
          paste0("[", input$lower_delta_1_Meeker, ", ", input$upper_delta_1_Meeker, "]"),
          paste0("[", input$lower_delta_2_Meeker, ", ", input$upper_delta_2_Meeker, "]"),
          paste0("[", input$lower_disp_kappa1_Meeker, ", ", input$upper_disp_kappa1_Meeker, "]"),
          paste0("[", input$lower_disp_kappa2_Meeker, ", ", input$upper_disp_kappa2_Meeker, "]")
        )
      )
      settings
    }, striped = TRUE, bordered = TRUE, spacing = "m", align = "c")
    
    shinyjs::hide("loading_spinner_Meeker")
  })
}