## March, 2022
## Function library
#
# Parameters: quantile (kappa) and a
# mu: kappa 
# sigma: a
# q0: prob. for the quantile
#
library(gamlss)
library(gamlss.dist)
#
## Distributions generated from the Gompertz distribution
# Distribution 1
# Transformation: exp(-x)
#
# Probability density function
# mu = kappa (quantile) and sigma = a
dLG <- function (y, mu = 0.5, sigma = 1, log = FALSE) 
{
    if (any(mu <= 0) | (any(mu >= 1)))
        stop(paste("mu must be in (0, 1) ", "\n", ""))    
    if (any(y < 0) | any(y > 1))
        stop(paste("y must be in [0, 1] ", "\n", ""))
    ly <- max(length(y), length(mu), length(sigma))
    y <- rep(y, length = ly)
    sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)
    b <- sigma * log(q0) / (1 - mu^(-sigma))

    fy <- b * y^(-sigma - 1) * exp(-b * (y^(-sigma) - 1) / sigma)
    if (log == TRUE) fy <- log(fy)
    fy   
}
#
# Cumulative distribution function
# mu = kappa (quantile) and sigma = a
pLG <- function (q, mu = 0.5, sigma = 1, lower.tail = TRUE, 
                    log.p = FALSE) {
   if (any(mu <= 0) | (any(mu >= 1)))
      stop(paste("mu must be in [0, 1] ", "\n", ""))    
   if (any(q < 0) | any(q > 1))
      stop(paste("y must be in [0, 1] ", "\n", ""))
    lq <- max(length(q), length(mu), length(sigma))
    q <- rep(q, length = lq)
    sigma <- rep(sigma, length = lq)
    mu <- rep(mu, length = lq)
    b <- sigma * log(q0) / (1 - mu^(-sigma))
    
    cdf <- exp(-b * (q^(-sigma) - 1) / sigma)
    if(lower.tail == FALSE) cdf <- 1 - cdf
    if (log.p == TRUE)  cdf <- log(cdf)   
    cdf    
}
#
# Random number generation
# mu = kappa (quantile) and sigma = a
# (sigma > 0)
#
rLG <- function (n, mu = 0.5, sigma = 1 , q0 = 0.5) {
  if (any(mu <= 0) | (any(mu >= 1)))
    stop(paste("mu must be in [0, 1] ", "\n", ""))    
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  b <- sigma * log(q0) / (1 - mu^(-sigma))
  n <- ceiling(n)
  r <- (1 - (sigma * log(runif(n))/b) )^(-1 / sigma)
  r <- ifelse(is.na(r)==TRUE,1,r)
} 
#
# Family
# mu = kappa (quantile) and sigma = a
LG = function (mu.link = "logit", sigma.link = "identity")
{
    mstats = checklink("mu.link", "Distribution1", 
        substitute(mu.link), c("logit", "identity"))
    dstats = checklink("sigma.link", "Distribution1", 
        substitute(sigma.link), c("identity"))
    structure(
        list(family = c("LG", "Distribution1"), 
        parameters = list(mu = TRUE, sigma = TRUE),
        nopar = 2, 
        type = "Continuous",
        mu.link = as.character(substitute(mu.link)),
        sigma.link = as.character(substitute(sigma.link)),
        mu.linkfun = mstats$linkfun,  
        sigma.linkfun = dstats$linkfun,
        mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
        mu.dr = mstats$mu.eta, 
        sigma.dr = dstats$mu.eta,

        dldm = function(y, mu, sigma, nu) {
            nd = gamlss:::numeric.deriv(dLG(y, mu, sigma, 
                 log = TRUE), "mu", delta = 1e-04)
            dldm = as.vector(attr(nd, "gradient"))
            dldm },     
        
        dldd = function(y, mu, sigma, nu) {
            nd = gamlss:::numeric.deriv(dLG(y, mu, sigma, 
                log = TRUE), "sigma", delta = 1e-04)
            dldd = as.vector(attr(nd, "gradient"))
            dldd },     
             
        d2ldm2 = function(y, mu, sigma, nu) {
            nd = gamlss:::numeric.deriv(dLG(y, mu, sigma,
                log = TRUE), "mu", delta = 1e-04)
            dldm = as.vector(attr(nd, "gradient"))
            d2ldm2 = -dldm * dldm
            d2ldm2 },
             
        d2ldmdd = function(y, mu, sigma, nu) {
            nd = gamlss:::numeric.deriv(dLG(y, mu, sigma, 
                log = TRUE), "mu", delta = 1e-04)
            dldm = as.vector(attr(nd, "gradient"))
            nd = gamlss:::numeric.deriv(dLG(y, mu, sigma,
                log = TRUE), "sigma", delta = 1e-04)
            dldd = as.vector(attr(nd, "gradient"))           
            d2ldmdd = -dldm * dldd
            d2ldmdd },

        d2ldd2 = function(y, mu, sigma, nu) {
            nd = gamlss:::numeric.deriv(dLG(y, mu, sigma, 
               log = TRUE), "sigma", delta = 1e-04)
            dldd = as.vector(attr(nd, "gradient"))
            d2ldd2 = -dldd * dldd
            d2ldd2 },

        G.dev.incr = function(y, mu, sigma,...)
          -2 * dLG(y, mu = mu, sigma = sigma, log = TRUE),
        rqres = expression(rqres(pfun = "pLG", type = "Continuous",
            ymin = 0, y = y, mu = mu, sigma = sigma)), 
        mu.initial = expression(mu = rep(0.5, length(y))), 
        sigma.initial = expression(sigma = rep(1, length(y))),
        mu.valid = function(mu) all(mu > 0) | all(mu < 1),
        sigma.valid = function(sigma) TRUE,
        y.valid = function(y) all(y >= 0) | all(y <= 1)), 
        class = c("gamlss.family", "family"))
}
#
## Distribution 2
# Transformation: 1 - exp(-x)
#
# Probability density function
# mu = kappa (quantile) and sigma = a
dCLG <- function (y, mu = 0.5, sigma = 1, log = FALSE) 
{
   if (any(mu <= 0) | (any(mu >= 1)))
      stop(paste("mu must be in (0, 1) ", "\n", ""))    
   if (any(y < 0) | any(y > 1))
      stop(paste("y must be in [0, 1] ", "\n", ""))
   ly <- max(length(y), length(mu), length(sigma))
   y <- rep(y, length = ly)
   sigma <- rep(sigma, length = ly)
   mu <- rep(mu, length = ly)
   b <- sigma * log(1 - q0) / (1 - (1 - mu)^(-sigma))
   
   fy <- b * (1 - y)^(-sigma - 1) * exp(-b * ((1 - y)^(-sigma) - 1) / 
            sigma)
   if (log == TRUE) fy <- log(fy)
   fy   
}
#
# Cumulative distribution function
# mu = kappa (quantile) and sigma = a
pCLG <- function (q, mu = 0.5, sigma = 1, lower.tail = TRUE, 
                    log.p = FALSE) {
   if (any(mu <= 0) | (any(mu >= 1)))
      stop(paste("mu must be in [0, 1] ", "\n", ""))    
   if (any(q < 0) | any(q > 1))
      stop(paste("y must be in [0, 1] ", "\n", ""))
   lq <- max(length(q), length(mu), length(sigma))
   q <- rep(q, length = lq)
   sigma <- rep(sigma, length = lq)
   mu <- rep(mu, length = lq)
   b <- sigma * log(1 - q0) / (1 - (1 - mu)^(-sigma))
   
   cdf <- 1 - exp(-b * ((1 - q)^(-sigma) - 1) / sigma)
   if(lower.tail == FALSE) cdf <- 1 - cdf
   if (log.p == TRUE)  cdf <- log(cdf)   
   cdf    
}
#
# Random number generation
# mu = kappa (quantile) and sigma = a
# (sigma > 0)
#
rCLG <- function (n, mu = 0.5, sigma = 1 , q0 = 0.5) {
  if (any(mu <= 0) | (any(mu >= 1)))
    stop(paste("mu must be in [0, 1] ", "\n", ""))    
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  b <- sigma * log(1 - q0) / (1 - (1 - mu)^(-sigma))
  n <- ceiling(n)
  r <- 1-(1- (sigma * log(1-runif(n)) / b) )^(-1 / sigma)
  #r <- min(r,0.98)
} 
#
# Family
# mu = kappa (quantile) and sigma = a
CLG = function (mu.link = "logit", sigma.link = "identity")
{
   mstats = checklink("mu.link", "Distribution2", 
                      substitute(mu.link), c("logit", "identity"))
   dstats = checklink("sigma.link", "Distribution2", 
                      substitute(sigma.link), c("identity"))
   structure(
      list(family = c("CLG", "Distribution2"), 
           parameters = list(mu = TRUE, sigma = TRUE),
           nopar = 2, 
           type = "Continuous",
           mu.link = as.character(substitute(mu.link)),
           sigma.link = as.character(substitute(sigma.link)),
           mu.linkfun = mstats$linkfun,  
           sigma.linkfun = dstats$linkfun,
           mu.linkinv = mstats$linkinv, 
           sigma.linkinv = dstats$linkinv,
           mu.dr = mstats$mu.eta, 
           sigma.dr = dstats$mu.eta,
           
           dldm = function(y, mu, sigma, nu) {
              nd = gamlss:::numeric.deriv(dCLG(y, mu, sigma, 
                     log = TRUE), "mu", delta = 1e-04)
              dldm = as.vector(attr(nd, "gradient"))
              dldm },     
           
           dldd = function(y, mu, sigma, nu) {
              nd = gamlss:::numeric.deriv(dCLG(y, mu, sigma, 
                     log = TRUE), "sigma", delta = 1e-04)
              dldd = as.vector(attr(nd, "gradient"))
              dldd },     
           
           d2ldm2 = function(y, mu, sigma, nu) {
              nd = gamlss:::numeric.deriv(dCLG(y, mu, sigma,
                     log = TRUE), "mu", delta = 1e-04)
              dldm = as.vector(attr(nd, "gradient"))
              d2ldm2 = -dldm * dldm
              d2ldm2 },
           
           d2ldmdd = function(y, mu, sigma, nu) {
              nd = gamlss:::numeric.deriv(dCLG(y, mu, sigma, 
                     log = TRUE), "mu", delta = 1e-04)
              dldm = as.vector(attr(nd, "gradient"))
              nd = gamlss:::numeric.deriv(dCLG(y, mu, sigma,
                     log = TRUE), "sigma", delta = 1e-04)
              dldd = as.vector(attr(nd, "gradient"))           
              d2ldmdd = -dldm * dldd
              d2ldmdd },
           
           d2ldd2 = function(y, mu, sigma, nu) {
              nd = gamlss:::numeric.deriv(dCLG(y, mu, sigma, 
                     log = TRUE), "sigma", delta = 1e-04)
              dldd = as.vector(attr(nd, "gradient"))
              d2ldd2 = -dldd * dldd
              d2ldd2 },
           
           G.dev.incr = function(y, mu, sigma,...)
              -2 * dCLG(y, mu = mu, sigma = sigma, log = TRUE),
           rqres = expression(rqres(pfun = "pCLG", type = "Continuous",
                     ymin = 0, y = y, mu = mu, sigma = sigma)), 
           mu.initial = expression(mu = rep(0.5, length(y))), 
           sigma.initial = expression(sigma = rep(1, length(y))),
           mu.valid = function(mu) all(mu > 0) | all(mu < 1),
           sigma.valid = function(sigma) TRUE,
           y.valid = function(y) all(y >= 0) | all(y <= 1)), 
      class = c("gamlss.family", "family"))
}