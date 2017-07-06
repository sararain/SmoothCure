################################################################################
# Generate Shared Frailty Survival Data
# Author: Yu Liu
# Date: June 2017
################################################################################
#------------------------------------------------------------------------------#
# Samuli(2000) Session 5
# Estimation of Multivariate Frailty Models Using Penalized Partial Likelihood
#------------------------------------------------------------------------------#

library(dplyr)
rfrailty.base <- function(sigma2, # Log Normal
                          size = 100, rep = 2, #
                          lambda = 0.1, # Exponential Distribution rate lambda
                          censor.rate = 0.9, # percentage of censored data
                          theta = c(1, -0.7, 0.5),
                          seed = NULL){

  if(!is.null(seed)) set.seed(seed)
  num.samples <- size * rep
  b <- rep(rnorm(size, sd = sqrt(sigma2)), each = rep)
  errors <- rexp(num.samples, rate = lambda)

  # Three Predictors
  X1 <- rnorm(num.samples)
  X2 <- rep(rnorm(size), each = rep)
  X3 <- rep(rbinom(n = size, size = 1, prob = 0.5), each = rep)
  X <- cbind(X1,X2,X3)

  time <- exp(- b - X %*% theta)*errors
  status <- rbinom(num.samples, 1, censor.rate)

  group <- rep(1:size, each = rep)
  df <- data.frame(time = time, status = status,
                   group = group, X1 = X1, X2 = X2, X3 = X3)
  df$group = as.factor(df$group)
  return(df)
}

rfrailty.transform <- function(sigma2, size = 100, rep = 2, ...){

  df = rfrailty.base(sigma2 = sigma2, size = size, rep = rep, ...)

  Z <- matrix(
    rep(diag(rep(1,size)), each = rep),
    nrow = size * rep,
    ncol = size
  )

  cov.m = matrix(rep(1,(size-1)^2), nrow = size -1)
  diag(cov.m) = 2

  eigen.cov.m = eigen(cov.m)
  sqrt.cov.m = eigen.cov.m$vector %*%
    diag(sqrt(eigen.cov.m$value)) %*%
    t(eigen.cov.m$vector)

  #sum(sqrt.cov.m%*%sqrt.cov.m - cov.m)
  transformed.Z = as.data.frame(Z[,-1] %*% sqrt.cov.m)

  df = cbind(df, transformed.Z)
  return(df)
}
