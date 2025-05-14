# Non censored data and two model follow lognormal(close form)

kldiv_ln_closeform <- function(xt, xr, vt ,vr) {
  
  mu_t <- log(xt) - 0.5 * log(1 + (vt / (xt^2)))
  mu_r <- log(xr) - 0.5 * log(1 + (vr / (xr^2)))
  var_t <- log(1 + (vt / (xt^2)))
  var_r <- log(1 + (vr / (xr^2)))
  
  log(sqrt(var_r)) - log(sqrt(var_t)) - ((var_r - var_t - (mu_r - mu_t)^2) / (2 * var_r))
}

#
# --------------------------------------------------
#
# Non censored data and two model follow lognormal(integrate)

kl_ln_int <- function(y, m1, m2, s1, s2) {
  #
  mu_t <- log(m1) - 0.5 * log(1 + (s1 / m1^2))
  mu_r <- log(m2) - 0.5 * log(1 + (s2 / m2^2))
  var_t <- log(1 + (s1 / m1^2))
  var_r <- log(1 + (s2 / m2^2))
  
  lpdf1 <- dlnorm(y, mu_t, sqrt(var_t), log = TRUE)
  lpdf2 <- dlnorm(y, mu_r, sqrt(var_r), log = TRUE)
  pdf1 <- exp(lpdf1)
  kl_int <- pdf1 * (lpdf1 - lpdf2)
  return(kl_int)
}

kldiv_ln_integrate <- function(xt, xr, st, sr) {
  kl_divs <- rep(0, length(xt))
  
  for (i in 1:length(xt)) {
    kl_divs[i] <- integrate(kl_ln_int,0, 100,
                            m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                            subdivisions = 1000, stop.on.error = FALSE)$value
  }
  return(kl_divs)
}

#
# --------------------------------------------------
#
# Non censored data and two model follow weibull(close form)

kldiv_wei_closeform <- function(xt, xr, vt, vr) {
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
  term7 <- (lambda1 / lambda2)^k2 * gamma( (k2 / k1) + 1 )
  
  kl_div <- term1 + term2 + term3 + term4 + term5 + term6 + term7
  return(kl_div)
}

#
# --------------------------------------------------
#
# Non censored data and two model follow weibull(integrate)

kl_int <- function(y, m1, m2, s1, s2) {
  #
  lpdf1 <- log(dweibull(y, 1/s1, exp(m1)) + 1e-12)
  lpdf2 <- log(dweibull(y, 1/s2, exp(m2)) + 1e-12)
  pdf1 <- exp(lpdf1)
  
  kl_int <- pdf1*(lpdf1 - lpdf2)
  return(kl_int)
}

kldiv_wei_integrate <- function(xt, xr, st, sr) {
  kl_divs <- rep(0, length(xt))
  
  for (i in 1:length(xt)) {
    upper <- qweibull(0.99, shape = 1/st[i], scale = exp(xt[i]))
    
    kl_divs[i] <- integrate(kl_int, 0, upper,
                            m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                            subdivisions = 1000000, stop.on.error = FALSE)$value
  }
  return(kl_divs)
}

#
# --------------------------------------------------
#
#
# --------------------------------------------------
#
#
# --------------------------------------------------
#
# Censored data and two model follow lognormal(c=5000)(CKL)

kl_lnln_observed <- function(y, m1, m2, s1, s2) {
  #
  lpdf1 <- dlnorm(y, m1, s1, log = TRUE)
  lpdf2 <- dlnorm(y, m2, s2, log = TRUE)
  pdf1 <- exp(lpdf1)
  
  val <- pdf1*(lpdf1 - lpdf2)
  return(val)
}

kl_lnln_censored <- function(y, m1, m2, s1, s2) {
  #
  lcdf1 <- log(1 - plnorm(y, m1, s1) + 1e-12)
  lcdf2 <- log(1 - plnorm(y, m2, s2) + 1e-12)
  cdf1 <- exp(lcdf1)
  #
  val <- cdf1*(lcdf1 - lcdf2)
  return(val)
}

kldiv_lnln_censored5000 <- function(xt, xr, st, sr) {
  tc <- 5000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    intg_part <- integrate(kl_lnln_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- kl_lnln_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}

#
# --------------------------------------------------
#
# Censored data and two model follow lognormal(c=5000)(CLW)

lw_lnln_observed <- function(y, m1, m2, s1, s2) {
  #
  pdf1 <- dlnorm(y, m1, s1)
  pdf2 <- dlnorm(y, m2, s2)
  
  val <- pdf1*(log(2*pdf1) - log(pdf1+pdf2))
  return(val)
}

lw_lnln_censored <- function(y, m1, m2, s1, s2) {
  #
  cdf1 <- (1 - plnorm(y, m1, s1) + 1e-12)
  cdf2 <- (1 - plnorm(y, m2, s2) + 1e-12)
  #
  val <- cdf1*(log(2*cdf1) - log(cdf1+cdf2))
  return(val)
}

lwdiv_lnln_censored5000 <- function(xt, xr, st, sr) {
  tc <- 5000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    intg_part <- integrate(lw_lnln_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- lw_lnln_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}

#
# --------------------------------------------------
#
# Censored data and two model follow lognormal(c=5000)(CB)

b_lnln_observed <- function(y, m1, m2, s1, s2) {
  #
  pdf1 <- dlnorm(y, m1, s1)
  pdf2 <- dlnorm(y, m2, s2)
  
  val <- sqrt(pdf1*pdf2)
  return(val)
}

b_lnln_censored <- function(y, m1, m2, s1, s2) {
  #
  cdf1 <- (1 - plnorm(y, m1, s1) + 1e-12)
  cdf2 <- (1 - plnorm(y, m2, s2) + 1e-12)
  #
  val <- sqrt(cdf1*cdf2)
  return(val)
}

bdiv_lnln_censored5000 <- function(xt, xr, st, sr) {
  tc <- 5000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    intg_part <- integrate(b_lnln_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- b_lnln_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}

#
# --------------------------------------------------
#
# Censored data and two model follow lognormal(c=5000)(CCS)

cs_lnln_observed <- function(y, m1, m2, s1, s2) {
  #
  pdf1 <- dlnorm(y, m1, s1)
  pdf2 <- dlnorm(y, m2, s2)
  
  val <- ((pdf1^2)/(pdf2+1e-12))
  return(val)
}

cs_lnln_censored <- function(y, m1, m2, s1, s2) {
  #
  cdf1 <- (1 - plnorm(y, m1, s1) + 1e-12)
  cdf2 <- (1 - plnorm(y, m2, s2) + 1e-12)
  #
  val <- ((cdf1^2/cdf2))
  return(val)
}

csdiv_lnln_censored5000 <- function(xt, xr, st, sr) {
  tc <- 5000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    intg_part <- integrate(cs_lnln_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- cs_lnln_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part -1
  }
  return(intVec)
}

#
# --------------------------------------------------
#
#
# --------------------------------------------------
#
#
# --------------------------------------------------
#
# Censored data and two model follow weibull(c=5000)(CKL)

kl_wewe_observed <- function(y, m1, m2, s1, s2) {
  #
  lpdf1 <- dweibull(y, 1/s1, exp(m1), log = TRUE)
  lpdf2 <- dweibull(y, 1/s2, exp(m2), log = TRUE)
  pdf1 <- exp(lpdf1)
  
  val <- pdf1*(lpdf1 - lpdf2)
  return(val)
}

kl_wewe_censored <- function(y, m1, m2, s1, s2) {
  #
  lcdf1 <- log(1 - pweibull(y, 1/s1, exp(m1)) + 1e-12)
  lcdf2 <- log(1 - pweibull(y, 1/s2, exp(m2)) + 1e-12)
  cdf1 <- exp(lcdf1)
  #
  val <- cdf1*(lcdf1 - lcdf2)
  return(val)
}

kldiv_wewe_censored5000 <- function(xt, xr, st, sr) {
  tc <- 5000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    intg_part <- integrate(kl_wewe_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- kl_wewe_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}

#
# --------------------------------------------------
#
# Censored data and two model follow weibull(c=5000)(CLW)

lw_wewe_observed <- function(y, m1, m2, s1, s2) {
  #
  pdf1 <- dweibull(y, 1/s1, exp(m1))
  pdf2 <- dweibull(y, 1/s2, exp(m2))
  
  val <- pdf1*(log(2*pdf1) - log(pdf1+pdf2))
  return(val)
}

lw_wewe_censored <- function(y, m1, m2, s1, s2) {
  #
  cdf1 <- (1 - pweibull(y, 1/s1, exp(m1)) + 1e-12)
  cdf2 <- (1 - pweibull(y, 1/s2, exp(m2)) + 1e-12)
  #
  val <- cdf1*(log(2*cdf1) - log(cdf1+cdf2))
  return(val)
}

lwdiv_wewe_censored5000 <- function(xt, xr, st, sr) {
  tc <- 5000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    intg_part <- integrate(lw_wewe_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- lw_wewe_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}

#
# --------------------------------------------------
#
# Censored data and two model follow weibull(c=5000)(CB)

b_wewe_observed <- function(y, m1, m2, s1, s2) {
  #
  pdf1 <- dweibull(y, 1/s1, exp(m1))
  pdf2 <- dweibull(y, 1/s2, exp(m2))
  
  val <- sqrt(pdf1*pdf2)
  return(val)
}

b_wewe_censored <- function(y, m1, m2, s1, s2) {
  #
  cdf1 <- (1 - pweibull(y, 1/s1, exp(m1)) + 1e-12)
  cdf2 <- (1 - pweibull(y, 1/s2, exp(m2)) + 1e-12)
  #
  val <- sqrt(cdf1*cdf2)
  return(val)
}

bdiv_wewe_censored5000 <- function(xt, xr, st, sr) {
  tc <- 5000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    intg_part <- integrate(b_wewe_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- b_wewe_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}

#
# --------------------------------------------------
#
# Censored data and two model follow weibull(c=5000)(CCS)

cs_wewe_observed <- function(y, m1, m2, s1, s2) {
  #
  pdf1 <- dweibull(y, 1/s1, exp(m1))
  pdf2 <- dweibull(y, 1/s2, exp(m2))
  
  val <- ((pdf1^2)/pdf2)
  return(val)
}

cs_wewe_censored <- function(y, m1, m2, s1, s2) {
  #
  cdf1 <- (1 - pweibull(y, 1/s1, exp(m1)) + 1e-12)
  cdf2 <- (1 - pweibull(y, 1/s2, exp(m2)) + 1e-12)
  #
  val <- ((cdf1^2/cdf2))
  return(val)
}

csdiv_wewe_censored5000 <- function(xt, xr, st, sr) {
  tc <- 5000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    intg_part <- integrate(cs_wewe_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- cs_wewe_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part -1
  }
  return(intVec)
}

#
# --------------------------------------------------
#
#
# --------------------------------------------------
#
#
# --------------------------------------------------
# two model follow different distribution

lnIsTrue_observed <- function(y, m1, m2, s1, s2) {
  #
  lpdf1 <- dlnorm(y, m1, s1, log = TRUE)
  lpdf2 <- log(dweibull(y, 1/s2, exp(m2)) + 1e-12)
  pdf1 <- exp(lpdf1)
  #
  val <- pdf1*(lpdf1 - lpdf2)
  return(val)
}

lnIsTrue_censored <- function(y, m1, m2, s1, s2) {
  #
  lcdf1 <- log(1 - plnorm(y, m1, s1) + 1e-12)
  lcdf2 <- log(1 - pweibull(y, 1/s2, exp(m2)) + 1e-12)
  cdf1 <- exp(lcdf1)
  #
  val <- cdf1*(lcdf1 - lcdf2)
  return(val)
}

wbIsTrue_observed <- function(y, m1, m2, s1, s2) {
  #
  lpdf1 <- dweibull(y, 1/s1, exp(m1), log = TRUE)
  lpdf2 <- dlnorm(y, m2, s2, log = TRUE)
  pdf1 <- exp(lpdf1)
  #
  val <- pdf1*(lpdf1 - lpdf2)
  return(val)
}

wbIsTrue_censored <- function(y, m1, m2, s1, s2) {
  #
  lcdf1 <- log(1 - pweibull(y, 1/s1, exp(m1)) + 1e-12)
  lcdf2 <- log(1 - plnorm(y, m2, s2) + 1e-12)
  cdf1 <- exp(lcdf1)
  #
  val <- cdf1*(lcdf1 - lcdf2)
  return(val)
}

kldiv_censored_lnIsTrue_tc_5000 <- function(xt, xr, st, sr) {
  tc <- 5000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    intg_part <- integrate(lnIsTrue_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- lnIsTrue_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}

kldiv_censored_wbIsTrue_tc_5000  <- function(xt, xr, st, sr) {
  tc <- 5000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    intg_part <- integrate(wbIsTrue_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- wbIsTrue_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}

kldiv_censored_lnIsTrue_tc_1000 <- function(xt, xr, st, sr) {
  tc <- 1000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    intg_part <- integrate(lnIsTrue_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- lnIsTrue_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}

kldiv_censored_wbIsTrue_tc_1000  <- function(xt, xr, st, sr) {
  tc <- 1000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    intg_part <- integrate(wbIsTrue_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- wbIsTrue_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}
