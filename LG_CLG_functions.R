pdfLGq<-function(y,k,a)
{
  b <- (a * log(0.5)) / (1-(k)^-a)
  b*((y)^(-a-1))*exp(b/a)*exp(-(b*(y)^(-a))/a)
}

pdfLGq0<-function(y,k,a,q)
{
  b <- (a * log(q)) / (1-(k)^-a)
  b*((y)^(-a-1))*exp(b/a)*exp(-(b*(y)^(-a))/a)
}

pdfLG<-function(y,a,b)
{
  b*((y)^(-a-1))*exp(b/a)*exp(-(b*(y)^(-a))/a)
}

pdfCLG<-function(y,a,b)
{
  b*((1-y)^(-a-1))*exp(b/a)*exp(-(b*(1-y)^(-a))/a)
}

pdfCLGq0<-function(y,k,a,q)
{
  b <- (a*log(1-q)) / (1-(1-k)^(-a))
  b*((1-y)^(-a-1))*exp(b/a)*exp(-(b*(1-y)^(-a))/a)
}

LLLG <- function(a, b) {
  pdf <- b*((y)^(-a-1))*exp(b/a)*exp(-(b*(y)^(-a))/a)
  #
  -sum(log(pdf))
}

LLLGquantile <- function(k,a) {
  b <- (a * log(0.5)) / (1-(k)^-a)
  pdf <- b*((y)^(-a-1))*exp(b/a)*exp(-(b*(y)^(-a))/a)
  #
  -sum(log(pdf))
}

LLLGquantile0 <- function(k,a,q) {
  b <- (a * log(q)) / (1-(k)^-a)
  pdf <- b*((y)^(-a-1))*exp(b/a)*exp(-(b*(y)^(-a))/a)
  #
  -sum(log(pdf))
}

LLCLG <- function(a, b) {
  pdf <- b*((1-y)^(-a-1))*exp(b/a)*exp(-(b*(1-y)^(-a))/a)
  #
  -sum(log(pdf))
}

LLCLGquantile <- function(k, a) {
  b <- (a*log(1-0.5)) / (1-(1-k)^(-a))
  pdf <- b*((1-y)^(-a-1))*exp(b/a)*exp(-(b*(1-y)^(-a))/a)
  #
  -sum(log(pdf))
}

LLCLGquantile0 <- function(k, a, q) {
  b <- (a*log(1-q)) / (1-(1-k)^(-a))
  pdf <- b*((1-y)^(-a-1))*exp(b/a)*exp(-(b*(1-y)^(-a))/a)
  #
  -sum(log(pdf))
}

LLBeta <- function(shape1, shape2) {
  pdf <- dbeta(y, shape1, shape2)
  #
  -sum(log(pdf))
}

LLBetaOne <- function(mu, sigma, nu) {
  pdf <- dBEOI(y, mu, sigma, nu, log = FALSE)
  #
  -sum(log(pdf))
}
