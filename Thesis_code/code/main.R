install.packages("devtools")
devtools::install_github("PingYangChen/DiscrimOD")
library(DiscrimOD)
source("/home/kylin/Rcode/Paper.code/divFunctions.R")

#
# Fidalgo-example-settings
#

af1_mean_Fidalgo <- function(x, p) (p[1]*x)/(p[2] + x) + p[3]*x
af2_mean_Fidalgo <- function(x, p) (p[1]*x)/(p[2] + x)
af1_disp_Fidalgo <- function(x, p) rep(p[1], length(x))
af2_disp_Fidalgo <- function(x, p) rep(p[1], length(x))
#
af1_para_Fidalgo <- c(1, 1, 1)
#
model_info_Fidalgo <- list(
  # The first list should be the true model and the specified nominal values
  list(mean = af1_mean_Fidalgo, disp = af1_disp_Fidalgo, meanPara = af1_para_Fidalgo, dispPara = 1),
  # Then the rival models are listed accordingly. We also need to specify the model space.
  list(mean = af2_mean_Fidalgo, disp = af2_disp_Fidalgo,
       meanParaLower = rep(0.1, 2), meanParaUpper = rep(100, 2),
       dispParaLower = c(1), dispParaUpper = c(1) )
)

nSupp_Fidalgo <- 3
dsRange_Fidalgo <- c(0.1, 5)
#
#
# ----------------------------------------------------------------------
#
#
#
# Fidalgo-example-lognormal-closeform

# Initialize PSO and BFGS options
PSO_INFO <- getPSOInfo(nSwarm = 64, maxIter = 200)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 2)

# Find KL-optimal design for models
res_Fidalgo <- DiscrimOD(MODEL_INFO = model_info_Fidalgo, DISTANCE = kldiv_ln_closeform,
                 nSupp = nSupp_Fidalgo, dsLower = dsRange_Fidalgo[1], dsUpper = dsRange_Fidalgo[2], crit_type = "pair_fixed_true",
                 PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Fidalgo$BESTDESIGN, 3) # The resulting design
res_Fidalgo$BESTVAL # The KL-optimal criterion value
res_Fidalgo$CPUTIME # CPU time

designCriterion(res_Fidalgo$BESTDESIGN, MODEL_INFO = model_info_Fidalgo, DISTANCE = kldiv_ln_closeform,
                dsLower = dsRange_Fidalgo[1], dsUpper = dsRange_Fidalgo[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

# Test optimality by equivalence theorem
eqv_Fidalgo <- equivalence(ngrid = 100, PSO_RESULT = res_Fidalgo, MODEL_INFO = model_info_Fidalgo,
                   DISTANCE = kldiv_ln_closeform,
                   dsLower = dsRange_Fidalgo[1], dsUpper = dsRange_Fidalgo[2], crit_type = "pair_fixed_true",
                   PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

# Draw the directional derivative curve
plot(eqv_Fidalgo$Grid_1, eqv_Fidalgo$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Fidalgo$BESTDESIGN[,1], rep(0, nrow(res_Fidalgo$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Fidalgo-example-lognormal-integrate

PSO_INFO <- getPSOInfo(nSwarm = 64, maxIter = 200)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 2)

res_Fidalgo <- DiscrimOD(MODEL_INFO = model_info_Fidalgo, DISTANCE = kldiv_ln_integrate,
                         nSupp = nSupp_Fidalgo, dsLower = dsRange_Fidalgo[1], dsUpper = dsRange_Fidalgo[2], crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Fidalgo$BESTDESIGN, 3)
res_Fidalgo$BESTVAL
res_Fidalgo$CPUTIME

designCriterion(res_Fidalgo$BESTDESIGN, MODEL_INFO = model_info_Fidalgo, DISTANCE = kldiv_ln_integrate,
                dsLower = dsRange_Fidalgo[1], dsUpper = dsRange_Fidalgo[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Fidalgo <- equivalence(ngrid = 100, PSO_RESULT = res_Fidalgo, MODEL_INFO = model_info_Fidalgo,
                           DISTANCE = kldiv_ln_integrate,
                           dsLower = dsRange_Fidalgo[1], dsUpper = dsRange_Fidalgo[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Fidalgo$Grid_1, eqv_Fidalgo$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Fidalgo$BESTDESIGN[,1], rep(0, nrow(res_Fidalgo$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Fidalgo-example-weibull-closeform
#

PSO_INFO <- getPSOInfo(nSwarm = 64, maxIter = 200)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 5)

res_Fidalgo <- DiscrimOD(MODEL_INFO = model_info_Fidalgo, DISTANCE = kldiv_wei_closeform,
                         nSupp = nSupp_Fidalgo, dsLower = dsRange_Fidalgo[1], dsUpper = dsRange_Fidalgo[2], crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Fidalgo$BESTDESIGN, 3)
res_Fidalgo$BESTVAL
res_Fidalgo$CPUTIME

designCriterion(res_Fidalgo$BESTDESIGN, MODEL_INFO = model_info_Fidalgo, DISTANCE = kldiv_wei_closeform,
                dsLower = dsRange_Fidalgo[1], dsUpper = dsRange_Fidalgo[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Fidalgo <- equivalence(ngrid = 100, PSO_RESULT = res_Fidalgo, MODEL_INFO = model_info_Fidalgo,
                           DISTANCE = kldiv_wei_closeform,
                           dsLower = dsRange_Fidalgo[1], dsUpper = dsRange_Fidalgo[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Fidalgo$Grid_1, eqv_Fidalgo$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Fidalgo$BESTDESIGN[,1], rep(0, nrow(res_Fidalgo$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Fidalgo-example-weibull-integrate
#

PSO_INFO <- getPSOInfo(nSwarm = 64, maxIter = 200)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 5)

res_Fidalgo <- DiscrimOD(MODEL_INFO = model_info_Fidalgo, DISTANCE = kldiv_wei_integrate,
                         nSupp = nSupp_Fidalgo, dsLower = dsRange_Fidalgo[1], dsUpper = dsRange_Fidalgo[2], crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Fidalgo$BESTDESIGN, 3)
res_Fidalgo$BESTVAL
res_Fidalgo$CPUTIME

designCriterion(res_Fidalgo$BESTDESIGN, MODEL_INFO = model_info_Fidalgo, DISTANCE = kldiv_wei_integrate,
                dsLower = dsRange_Fidalgo[1], dsUpper = dsRange_Fidalgo[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Fidalgo <- equivalence(ngrid = 100, PSO_RESULT = res_Fidalgo, MODEL_INFO = model_info_Fidalgo,
                           DISTANCE = kldiv_wei_integrate,
                           dsLower = dsRange_Fidalgo[1], dsUpper = dsRange_Fidalgo[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Fidalgo$Grid_1, eqv_Fidalgo$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Fidalgo$BESTDESIGN[,1], rep(0, nrow(res_Fidalgo$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Arrhenius-example-settings-lognormal-given both variance 0.9780103
#
af1_mean_Arrhenius <- function(x, p) p[1] + p[2] * (11605/(x+273.15)) + p[3] * (11605/(x+273.15))^2
af2_mean_Arrhenius <- function(x, p) p[1] + p[2] * (11605/(x+273.15))
af1_disp_Arrhenius <- function(x, p) rep(p[1], length(x))
af2_disp_Arrhenius <- function(x, p) rep(p[1], length(x))
#
af1_para_Arrhenius <- c(-5, -1.5, 0.05)
#
model_info_Arrhenius <- list(
  # The first list should be the true model and the specified nominal values
  list(mean = af1_mean_Arrhenius, disp = af1_disp_Arrhenius, meanPara = af1_para_Arrhenius, dispPara = 0.9780103),
  # Then the rival models are listed accordingly. We also need to specify the model space.
  list(mean = af2_mean_Arrhenius, disp = af2_disp_Arrhenius,
       meanParaLower = c(-100, 0.1), meanParaUpper = c(-10, 5),
       dispParaLower = c(0.9780103), dispParaUpper = c(0.9780103) )
)

nSupp_Arrhenius <- 3
dsRange_Arrhenius <- c(10, 80)
PSO_INFO <- getPSOInfo(nSwarm = 64, maxIter = 200)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 50)

#
#
# ----------------------------------------------------------------------
#
#
#
# Find CKL-optimal design for models(c=5000)

res_Arrhenius <- DiscrimOD(MODEL_INFO = model_info_Arrhenius, DISTANCE = kldiv_lnln_censored5000,
                         nSupp = nSupp_Arrhenius, dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Arrhenius$BESTDESIGN, 3)
res_Arrhenius$BESTVAL
res_Arrhenius$CPUTIME

designCriterion(res_Arrhenius$BESTDESIGN, MODEL_INFO = model_info_Arrhenius, DISTANCE = kldiv_lnln_censored5000,
                dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Arrhenius <- equivalence(ngrid = 100, PSO_RESULT = res_Arrhenius, MODEL_INFO = model_info_Arrhenius,
                           DISTANCE = kldiv_lnln_censored5000,
                           dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Arrhenius$Grid_1, eqv_Arrhenius$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Arrhenius$BESTDESIGN[,1], rep(0, nrow(res_Arrhenius$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Find CLW-optimal design for models(c=5000)

res_Arrhenius <- DiscrimOD(MODEL_INFO = model_info_Arrhenius, DISTANCE = lwdiv_lnln_censored5000,
                           nSupp = nSupp_Arrhenius, dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Arrhenius$BESTDESIGN, 3)
res_Arrhenius$BESTVAL
res_Arrhenius$CPUTIME

designCriterion(res_Arrhenius$BESTDESIGN, MODEL_INFO = model_info_Arrhenius, DISTANCE = lwdiv_lnln_censored5000,
                dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Arrhenius <- equivalence(ngrid = 100, PSO_RESULT = res_Arrhenius, MODEL_INFO = model_info_Arrhenius,
                             DISTANCE = lwdiv_lnln_censored5000,
                             dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Arrhenius$Grid_1, eqv_Arrhenius$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Arrhenius$BESTDESIGN[,1], rep(0, nrow(res_Arrhenius$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Find CB-optimal design for models(c=5000)

res_Arrhenius <- DiscrimOD(MODEL_INFO = model_info_Arrhenius, DISTANCE = bdiv_lnln_censored5000,
                           nSupp = nSupp_Arrhenius, dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Arrhenius$BESTDESIGN, 3)
res_Arrhenius$BESTVAL
res_Arrhenius$CPUTIME

designCriterion(res_Arrhenius$BESTDESIGN, MODEL_INFO = model_info_Arrhenius, DISTANCE = bdiv_lnln_censored5000,
                dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Arrhenius <- equivalence(ngrid = 100, PSO_RESULT = res_Arrhenius, MODEL_INFO = model_info_Arrhenius,
                             DISTANCE = bdiv_lnln_censored5000,
                             dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Arrhenius$Grid_1, eqv_Arrhenius$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Arrhenius$BESTDESIGN[,1], rep(0, nrow(res_Arrhenius$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Find CCS-optimal design for models(c=5000)

res_Arrhenius <- DiscrimOD(MODEL_INFO = model_info_Arrhenius, DISTANCE = csdiv_lnln_censored5000,
                           nSupp = nSupp_Arrhenius, dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Arrhenius$BESTDESIGN, 3)
res_Arrhenius$BESTVAL
res_Arrhenius$CPUTIME

designCriterion(res_Arrhenius$BESTDESIGN, MODEL_INFO = model_info_Arrhenius, DISTANCE = csdiv_lnln_censored5000,
                dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Arrhenius <- equivalence(ngrid = 100, PSO_RESULT = res_Arrhenius, MODEL_INFO = model_info_Arrhenius,
                             DISTANCE = csdiv_lnln_censored5000,
                             dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Arrhenius$Grid_1, eqv_Arrhenius$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Arrhenius$BESTDESIGN[,1], rep(0, nrow(res_Arrhenius$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Arrhenius-example-settings-weibull-given both variance 0.9780103
#

af1_mean_Arrhenius <- function(x, p) p[1] + p[2] * (11605/(x+273.15)) + p[3] * (11605/(x+273.15))^2
af2_mean_Arrhenius <- function(x, p) p[1] + p[2] * (11605/(x+273.15))
af1_disp_Arrhenius <- function(x, p) rep(p[1], length(x))
af2_disp_Arrhenius <- function(x, p) rep(p[1], length(x))
#
af1_para_Arrhenius <- c(-5, -1.5, 0.05)
#
model_info_Arrhenius <- list(
  # The first list should be the true model and the specified nominal values
  list(mean = af1_mean_Arrhenius, disp = af1_disp_Arrhenius, meanPara = af1_para_Arrhenius, dispPara = 0.9780103),
  # Then the rival models are listed accordingly. We also need to specify the model space.
  list(mean = af2_mean_Arrhenius, disp = af2_disp_Arrhenius,
       meanParaLower = c(-100, 0.1), meanParaUpper = c(-10, 5),
       dispParaLower = c(0.9780103), dispParaUpper = c(0.9780103) )
)

nSupp_Arrhenius <- 3
dsRange_Arrhenius <- c(10, 80)
PSO_INFO <- getPSOInfo(nSwarm = 64, maxIter = 200)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 50)

#
#
# ----------------------------------------------------------------------
#
#
#
# Find CKL-optimal design for models(c=5000)

res_Arrhenius <- DiscrimOD(MODEL_INFO = model_info_Arrhenius, DISTANCE = kldiv_wewe_censored5000,
                           nSupp = nSupp_Arrhenius, dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Arrhenius$BESTDESIGN, 3)
res_Arrhenius$BESTVAL
res_Arrhenius$CPUTIME

designCriterion(res_Arrhenius$BESTDESIGN, MODEL_INFO = model_info_Arrhenius, DISTANCE = kldiv_wewe_censored5000,
                dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Arrhenius <- equivalence(ngrid = 100, PSO_RESULT = res_Arrhenius, MODEL_INFO = model_info_Arrhenius,
                             DISTANCE = kldiv_wewe_censored5000,
                             dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Arrhenius$Grid_1, eqv_Arrhenius$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Arrhenius$BESTDESIGN[,1], rep(0, nrow(res_Arrhenius$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Find CLW-optimal design for models(c=5000)

res_Arrhenius <- DiscrimOD(MODEL_INFO = model_info_Arrhenius, DISTANCE = lwdiv_wewe_censored5000,
                           nSupp = nSupp_Arrhenius, dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Arrhenius$BESTDESIGN, 3)
res_Arrhenius$BESTVAL
res_Arrhenius$CPUTIME

designCriterion(res_Arrhenius$BESTDESIGN, MODEL_INFO = model_info_Arrhenius, DISTANCE = lwdiv_wewe_censored5000,
                dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Arrhenius <- equivalence(ngrid = 100, PSO_RESULT = res_Arrhenius, MODEL_INFO = model_info_Arrhenius,
                             DISTANCE = lwdiv_wewe_censored5000,
                             dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Arrhenius$Grid_1, eqv_Arrhenius$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Arrhenius$BESTDESIGN[,1], rep(0, nrow(res_Arrhenius$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Find CB-optimal design for models(c=5000)

res_Arrhenius <- DiscrimOD(MODEL_INFO = model_info_Arrhenius, DISTANCE = bdiv_wewe_censored5000,
                           nSupp = nSupp_Arrhenius, dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Arrhenius$BESTDESIGN, 3)
res_Arrhenius$BESTVAL
res_Arrhenius$CPUTIME

designCriterion(res_Arrhenius$BESTDESIGN, MODEL_INFO = model_info_Arrhenius, DISTANCE = bdiv_wewe_censored5000,
                dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Arrhenius <- equivalence(ngrid = 100, PSO_RESULT = res_Arrhenius, MODEL_INFO = model_info_Arrhenius,
                             DISTANCE = bdiv_wewe_censored5000,
                             dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Arrhenius$Grid_1, eqv_Arrhenius$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Arrhenius$BESTDESIGN[,1], rep(0, nrow(res_Arrhenius$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Find CCS-optimal design for models(c=5000)

res_Arrhenius <- DiscrimOD(MODEL_INFO = model_info_Arrhenius, DISTANCE = csdiv_wewe_censored5000,
                           nSupp = nSupp_Arrhenius, dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Arrhenius$BESTDESIGN, 3)
res_Arrhenius$BESTVAL
res_Arrhenius$CPUTIME

designCriterion(res_Arrhenius$BESTDESIGN, MODEL_INFO = model_info_Arrhenius, DISTANCE = csdiv_wewe_censored5000,
                dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Arrhenius <- equivalence(ngrid = 100, PSO_RESULT = res_Arrhenius, MODEL_INFO = model_info_Arrhenius,
                             DISTANCE = csdiv_wewe_censored5000,
                             dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Arrhenius$Grid_1, eqv_Arrhenius$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Arrhenius$BESTDESIGN[,1], rep(0, nrow(res_Arrhenius$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Arrhenius-example-settings
#

af1_mean_Arrhenius <- function(x, p) p[1] + p[2] * (11605/(x+273.15)) + p[3] * (11605/(x+273.15))^2
af2_mean_Arrhenius <- function(x, p) p[1] + p[2] * (11605/(x+273.15))
af1_disp_Arrhenius <- function(x, p) rep(p[1], length(x))
af2_disp_Arrhenius <- function(x, p) rep(p[1], length(x))
#
af1_para_Arrhenius <- c(-5, -1.5, 0.05)
#
nSupp_Arrhenius <- 3
dsRange_Arrhenius <- c(10, 80)
PSO_INFO <- getPSOInfo(nSwarm = 64, maxIter = 200)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 50)


#
#
# ----------------------------------------------------------------------
#
#
#
# Arrhenius-lognormal-given first model variance 0.9780103,another one [0.47,4.97] (c=5000)
#

model_info_Arrhenius <- list(
  list(mean = af1_mean_Arrhenius, disp = af1_disp_Arrhenius, meanPara = af1_para_Arrhenius, dispPara = 0.9780103),
  list(mean = af2_mean_Arrhenius, disp = af2_disp_Arrhenius,
       meanParaLower = c(-100, 0.1), meanParaUpper = c(-10, 5),
       dispParaLower = c(0.4780103), dispParaUpper = c(4.9780103) )
)

res_Arrhenius <- DiscrimOD(MODEL_INFO = model_info_Arrhenius, DISTANCE = kldiv_lnln_censored5000,
                           nSupp = nSupp_Arrhenius, dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Arrhenius$BESTDESIGN, 3)
res_Arrhenius$BESTVAL
res_Arrhenius$CPUTIME

designCriterion(res_Arrhenius$BESTDESIGN, MODEL_INFO = model_info_Arrhenius, DISTANCE = kldiv_lnln_censored5000,
                dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Arrhenius <- equivalence(ngrid = 100, PSO_RESULT = res_Arrhenius, MODEL_INFO = model_info_Arrhenius,
                             DISTANCE = kldiv_lnln_censored5000,
                             dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Arrhenius$Grid_1, eqv_Arrhenius$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Arrhenius$BESTDESIGN[,1], rep(0, nrow(res_Arrhenius$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Arrhenius-lognormal-given first model variance 1.4780103,another one [0.47,4.97](c=5000)
#

model_info_Arrhenius <- list(
  list(mean = af1_mean_Arrhenius, disp = af1_disp_Arrhenius, meanPara = af1_para_Arrhenius, dispPara = 1.4780103),
  list(mean = af2_mean_Arrhenius, disp = af2_disp_Arrhenius,
       meanParaLower = c(-100, 0.1), meanParaUpper = c(-10, 5),
       dispParaLower = c(0.4780103), dispParaUpper = c(4.9780103) )
)

res_Arrhenius <- DiscrimOD(MODEL_INFO = model_info_Arrhenius, DISTANCE = kldiv_lnln_censored5000,
                           nSupp = nSupp_Arrhenius, dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Arrhenius$BESTDESIGN, 3)
res_Arrhenius$BESTVAL
res_Arrhenius$CPUTIME

designCriterion(res_Arrhenius$BESTDESIGN, MODEL_INFO = model_info_Arrhenius, DISTANCE = kldiv_lnln_censored5000,
                dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Arrhenius <- equivalence(ngrid = 100, PSO_RESULT = res_Arrhenius, MODEL_INFO = model_info_Arrhenius,
                             DISTANCE = kldiv_lnln_censored5000,
                             dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Arrhenius$Grid_1, eqv_Arrhenius$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Arrhenius$BESTDESIGN[,1], rep(0, nrow(res_Arrhenius$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Arrhenius-weibull-given first model variance 0.9780103,another one [0.47,4.97] (c=5000)
#

model_info_Arrhenius <- list(
  list(mean = af1_mean_Arrhenius, disp = af1_disp_Arrhenius, meanPara = af1_para_Arrhenius, dispPara = 0.9780103),
  list(mean = af2_mean_Arrhenius, disp = af2_disp_Arrhenius,
       meanParaLower = c(-100, 0.1), meanParaUpper = c(-10, 5),
       dispParaLower = c(0.4780103), dispParaUpper = c(4.9780103) )
)

res_Arrhenius <- DiscrimOD(MODEL_INFO = model_info_Arrhenius, DISTANCE = kldiv_wewe_censored5000,
                           nSupp = nSupp_Arrhenius, dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Arrhenius$BESTDESIGN, 3)
res_Arrhenius$BESTVAL
res_Arrhenius$CPUTIME

designCriterion(res_Arrhenius$BESTDESIGN, MODEL_INFO = model_info_Arrhenius, DISTANCE = kldiv_wewe_censored5000,
                dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Arrhenius <- equivalence(ngrid = 100, PSO_RESULT = res_Arrhenius, MODEL_INFO = model_info_Arrhenius,
                             DISTANCE = kldiv_wewe_censored5000,
                             dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Arrhenius$Grid_1, eqv_Arrhenius$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Arrhenius$BESTDESIGN[,1], rep(0, nrow(res_Arrhenius$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Arrhenius-weibull-given first model variance 1.4780103,another one [0.47,4.97](c=5000)
#

model_info_Arrhenius <- list(
  list(mean = af1_mean_Arrhenius, disp = af1_disp_Arrhenius, meanPara = af1_para_Arrhenius, dispPara = 1.4780103),
  list(mean = af2_mean_Arrhenius, disp = af2_disp_Arrhenius,
       meanParaLower = c(-100, 0.1), meanParaUpper = c(-10, 5),
       dispParaLower = c(0.4780103), dispParaUpper = c(4.9780103) )
)

res_Arrhenius <- DiscrimOD(MODEL_INFO = model_info_Arrhenius, DISTANCE = kldiv_wewe_censored5000,
                           nSupp = nSupp_Arrhenius, dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Arrhenius$BESTDESIGN, 3)
res_Arrhenius$BESTVAL
res_Arrhenius$CPUTIME

designCriterion(res_Arrhenius$BESTDESIGN, MODEL_INFO = model_info_Arrhenius, DISTANCE = kldiv_wewe_censored5000,
                dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Arrhenius <- equivalence(ngrid = 100, PSO_RESULT = res_Arrhenius, MODEL_INFO = model_info_Arrhenius,
                             DISTANCE = kldiv_wewe_censored5000,
                             dsLower = dsRange_Arrhenius[1], dsUpper = dsRange_Arrhenius[2], crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Arrhenius$Grid_1, eqv_Arrhenius$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Arrhenius$BESTDESIGN[,1], rep(0, nrow(res_Arrhenius$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Meeker-Stress-Dependent Variance(true model follow lognormal,rival model follow weibull)(c=1000)

af1_mean_Meeker <- function(x, p) p[1] + p[2]*log(x - 75.71 + 1e-12)
af2_mean_Meeker <- function(x, p) p[1] + p[2]*log(x - 75.71 + 1e-12)
af1_disp_Meeker <- function(x, p) exp(p[1] + p[2]*log(x + 1e-12))
af2_disp_Meeker <- function(x, p) exp(p[1] + p[2]*log(x + 1e-12))
#
af1_mean_para_Meeker <- c(14.75, -1.39)
af1_disp_para_Meeker <- c(10.97, -2.50)
#
model_info <- list(
  list(mean = af1_mean_Meeker, disp = af1_disp_Meeker, meanPara = af1_mean_para_Meeker, dispPara = af1_disp_para_Meeker),
  list(mean = af2_mean_Meeker, disp = af2_disp_Meeker,
       meanParaLower = c(12.06, -2.02), meanParaUpper = c(17.44, -0.76),
       dispParaLower = c(10, -3.0), dispParaUpper = c(20, -0.01) )
)

nSupp_Meeker <- 4
dsRange_Meeker <- c(76, 150)
PSO_INFO <- getPSOInfo(nSwarm = 64, maxIter = 200)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 50)

#
#
# ----------------------------------------------------------------------
#
#
#
# Find CKL-optimal design for models(lognormal-integrate)

res_Meeker <- DiscrimOD(MODEL_INFO = model_info_Meeker, DISTANCE = kldiv_censored_lnIsTrue_tc_1000,
                 nSupp = nSupp_Meeker, dsLower = dsRange_Meeker[1], dsUpper = dsRange_Meeker[2], crit_type = "pair_fixed_true",
                 PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Meeker$BESTDESIGN, 3)
res_Meeker$BESTVAL
res_Meeker$CPUTIME

designCriterion(res_Meeker$BESTDESIGN, MODEL_INFO = model_info_Meeker, DISTANCE = kldiv_censored_lnIsTrue_tc_1000,
                dsLower = dsRange_Meeker[1], dsUpper = dsRange_Meeker[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Meeker <- equivalence(ngrid = 100, PSO_RESULT = res_Meeker, MODEL_INFO = model_info_Meeker,
                   DISTANCE = kldiv_censored_lnIsTrue_tc_1000,
                   dsLower = dsRange_Meeker[1], dsUpper = dsRange_Meeker[2], crit_type = "pair_fixed_true",
                   PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Meeker$Grid_1, eqv_Meeker$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Meeker$BESTDESIGN[,1], rep(0, nrow(res_Meeker$BESTDESIGN)), pch = 16)

#
#
# ----------------------------------------------------------------------
#
#
#
# Meeker-Stress-Dependent Variance(true model follow weibull,rival model follow lognormal)(c=1000)

af1_mean_Meeker <- function(x, p) p[1] + p[2]*log(x - 75.71 + 1e-12)
af2_mean_Meeker <- function(x, p) p[1] + p[2]*log(x - 75.71 + 1e-12)
af1_disp_Meeker <- function(x, p) exp(p[1] + p[2]*log(x + 1e-12))
af2_disp_Meeker <- function(x, p) exp(p[1] + p[2]*log(x + 1e-12))
#
af1_mean_para_Meeker <- c(14.75, -1.39)
af1_disp_para_Meeker <- c(10.97, -2.50)
#
model_info <- list(
  list(mean = af1_mean_Meeker, disp = af1_disp_Meeker, meanPara = af1_mean_para_Meeker, dispPara = af1_disp_para_Meeker),
  list(mean = af2_mean_Meeker, disp = af2_disp_Meeker,
       meanParaLower = c(12.06, -2.02), meanParaUpper = c(17.44, -0.76),
       dispParaLower = c(10, -3.0), dispParaUpper = c(20, -0.01) )
)

nSupp_Meeker <- 4
dsRange_Meeker <- c(76, 150)
PSO_INFO <- getPSOInfo(nSwarm = 64, maxIter = 200)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 50)

#
#
# ----------------------------------------------------------------------
#
#
#
# Find CKL-optimal design for models

res_Meeker <- DiscrimOD(MODEL_INFO = model_info_Meeker, DISTANCE = kldiv_censored_wbIsTrue_tc_1000,
                        nSupp = nSupp_Meeker, dsLower = dsRange_Meeker[1], dsUpper = dsRange_Meeker[2], crit_type = "pair_fixed_true",
                        PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)

round(res_Meeker$BESTDESIGN, 3)
res_Meeker$BESTVAL
res_Meeker$CPUTIME

designCriterion(res_Meeker$BESTDESIGN, MODEL_INFO = model_info_Meeker, DISTANCE = kldiv_censored_wbIsTrue_tc_1000,
                dsLower = dsRange_Meeker[1], dsUpper = dsRange_Meeker[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_Meeker <- equivalence(ngrid = 100, PSO_RESULT = res_Meeker, MODEL_INFO = model_info_Meeker,
                          DISTANCE = kldiv_censored_wbIsTrue_tc_1000,
                          dsLower = dsRange_Meeker[1], dsUpper = dsRange_Meeker[2], crit_type = "pair_fixed_true",
                          PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv_Meeker$Grid_1, eqv_Meeker$DirDeriv, type = "l", col = "blue",
     main = "", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res_Meeker$BESTDESIGN[,1], rep(0, nrow(res_Meeker$BESTDESIGN)), pch = 16)
