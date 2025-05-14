install.packages(c("devtools", "Rcpp", "RcppArmadillo"))
devtools::install_github("PingYangChen/DiscrimOD")
library(DiscrimOD)

### Define the mean response and dispersion functions for both the true model and the rival model
af1_mean <- function(x, p) p[1] + p[2] * (11605/(x+273.15)) + p[3] * (11605/(x+273.15))^2
af2_mean <- function(x, p) p[1] + p[2] * (11605/(x+273.15))
af1_disp <- function(x, p) rep(p[1], length(x))
af2_disp <- function(x, p) rep(p[1], length(x))

### Set the nominal values for the true model
af1_para <- c(-5, -1.5, 0.05)
model_info <- list(
  ### The first list should be the true model and the specified nominal values
  list(mean = af1_mean, disp = af1_disp, meanPara = af1_para, dispPara = 0.9780103),
  ### Then the rival models are listed accordingly. We also need to specify the model space.
  list(mean = af2_mean, disp = af2_disp,
       meanParaLower = c(-100, 0.1), meanParaUpper = c(-10, 5),
       dispParaLower = c(0.9780103), dispParaUpper = c(0.9780103) )
)

### Define the CKL divergence function
# xt is the mean values of the true model
# xr is the mean values of the rival model
kl_lnln_observed <- function(y, m1, m2, s1, s2) {
  lpdf1 <- dlnorm(y, m1, s1, log = TRUE)
  lpdf2 <- dlnorm(y, m2, s2, log = TRUE)
  pdf1 <- exp(lpdf1)
  val <- pdf1*(lpdf1 - lpdf2)
  return(val)
}

kl_lnln_censored <- function(y, m1, m2, s1, s2) {
  lcdf1 <- log(1 - plnorm(y, m1, s1) + 1e-12)
  lcdf2 <- log(1 - plnorm(y, m2, s2) + 1e-12)
  cdf1 <- exp(lcdf1)
  val <- cdf1*(lcdf1 - lcdf2)
  return(val)
}

kldiv_lnln_censored5000 <- function(xt, xr, st, sr) {
  tc <- 5000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    intg_part <- integrate(kl_lnln_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], 
                           s1 = st[i], s2 = sr[i],
                           subdivisions = 100,
                           stop.on.error = FALSE)$value
    cens_part <- kl_lnln_censored(tc, m1 = xt[i], m2 = xr[i],
                                  s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}

### Setting up the algorithm parameters
PSO_INFO <- getPSOInfo(nSwarm = 64, maxIter = 200)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 50)

### The number of support points and the range of the support points
nSupp <- 3
dsRange <- c(10, 80)

### The output of the DiscrimOD function
res <- DiscrimOD(MODEL_INFO = model_info, DISTANCE = kldiv_lnln_censored5000,
                 nSupp = nSupp, dsLower = dsRange[1], dsUpper = dsRange[2],
                 crit_type = "pair_fixed_true",
                 PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                 seed = 100, verbose = TRUE)

### The best design points
round(res$BESTDESIGN, 3)

### The best design criterion value
res$BESTVAL

### The elapsed CPU time (in seconds) during the optimization process
res$CPUTIME

### The estimated parameters of the rival model and recompute the criterion value
designCriterion(res$BESTDESIGN, MODEL_INFO = model_info,
                DISTANCE = kldiv_lnln_censored5000,
                dsLower = dsRange[1], dsUpper = dsRange[2],
                crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

### Compute the directional derivative
eqv <- equivalence(ngrid = 100, PSO_RESULT = res,
                   MODEL_INFO = model_info,
                   DISTANCE = kldiv_lnln_censored5000,
                   dsLower = dsRange[1], dsUpper = dsRange[2],
                   crit_type = "pair_fixed_true",
                   PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

# Draw the directional derivative curve
plot(eqv$Grid_1, eqv_Arrhenius$DirDeriv, type = "l",
     col = "blue", main = "",
     xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res$BESTDESIGN[,1], rep(0, nrow(res$BESTDESIGN)), pch = 16)