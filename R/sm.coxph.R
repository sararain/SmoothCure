################################################################################
# sm.coxph
# Author: Yu Liu
# Date: June 2017
################################################################################

library(dplyr)
sm.coxph <- function(formula, data = NULL){
  # Sort Data

  mf <- model.frame(formula = formula, data = data)
  time <- model.response(mf)[,1]
  data = data[order(time),]

  XZ <- model.matrix(object = formula, data = data)
  XZ <- XZ[,colnames(XZ) %>% setdiff(.,"(Intercept)")]

  pen.terms <- regexpr("ridge", colnames(XZ)) == 1
  Z <- XZ[, colnames(XZ)[pen.terms]] %>% as.matrix
  X <- XZ[, colnames(XZ)[!pen.terms]] %>% as.matrix

  if(length(colnames(XZ)[!pen.terms]) > 0){
    fix.terms.str <- paste0(' + ', paste0(colnames(XZ)[!pen.terms], collapse = ' + '))
  } else{
    fix.terms.str <- ''
  }
  # Penalized Terms for ridge formula
  pen.terms.str <- regmatches(
    colnames(XZ),
    regexpr('\\).+', colnames(XZ))
  ) %>%
    sapply(., function(x) substring(x, 2), USE.NAMES = FALSE) %>%
    paste0(., collapse = ', ')

  rm(XZ)

  ##### Precalculate Constant Value #####
  q = sum(pen.terms) # Number of random effects
  ZxtZ <- apply(Z,1, function(x) x%*%t(x))
  Ik <- diag(rep(1,q))

  ##### Initial Value #####
  fit <- coxph(formula, data = data, model = TRUE)
  b = summary(fit)$coefficient[,1]
  beta = b[names(b)[!pen.terms]] %>% as.matrix
  b = b[names(b)[pen.terms]] %>% as.matrix

  unique.time.index <- cumsum(c(TRUE, diff(data$time) != 0))
  hazard <- basehaz(fit, centered=FALSE)$hazard[unique.time.index] %>% as.matrix

  old_sigma = sigma = 10
  error = 1
  n = 1
  ##### Optimization Iteration #####
  while(error > 0.0001){

    hessian = matrix(ZxtZ %*% (exp(X%*%beta + Z%*%b)*hazard), q, q) + 1/sigma*Ik
    sigma = sum(b^2)/q + sum((1/eigen(hessian)$value))/q
    surv_formula = formula(
      paste0('Surv(time,status) ~ ridge(',
             pen.terms.str,
             ', scale = FALSE, theta = ',
             1/sigma,')',
             fix.terms.str)
    )
    surv.fit <- coxph(surv_formula, data = data, model = TRUE)
    b = summary(surv.fit)$coefficient[,1]
    beta = b[names(b)[!pen.terms]] %>% as.matrix
    b = b[names(b)[pen.terms]] %>% as.matrix
    hazard <- basehaz(surv.fit, centered=FALSE)$hazard[unique.time.index] %>% as.matrix
    error = abs(sigma - old_sigma)
    old_sigma = sigma
    n = n + 1
    if(n%%200 == 0){
      break
    }
  }
  return(list(sigma = sigma, beta = beta, surv.fit = surv.fit))
}









