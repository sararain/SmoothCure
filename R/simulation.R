################################################################################
# Simulation
# Samuli(2000) Session 5
# Estimation of Multivariate Frailty Models Using Penalized Partial Likelihood
# Author: Yu Liu
# Date: June 2017
################################################################################

#------------------------------------------------------------------------------#
# The first simulation study looks at the performance of frailty variance
# estimates in the absence of fixed effects.
#------------------------------------------------------------------------------#

.simulation1 <- function(sigma2 = 1, size = 50, rep = 5){
  #sigma2 = 1; size = 50; rep = 5
  res = c()
  data = rfrailty.transform(
    sigma2 = sigma2, size = size, rep = rep, theta = c(0,0,0)
  )
  data = data[order(data$time), ]
  ########
  parfm.fit <- parfm(Surv(time, status) ~ 1, cluster = 'group', data = data,
                     dist = "exponential", frailty = 'lognormal')
  res = rbind(res, parfm.fit[,1])
  ########
  # fP.fit <- frailtyPenal(
  #   Surv(time, status) ~ cluster(group),
  #   RandDist = 'LogN',
  #   data = data,
  #   hazard = "Weibull"
  # )
  # fP.pred <- c(fP.fit$sigma2, 1/fP.fit$scale.weib[1], fP.fit$coef)
  # res = rbind(res, fP.pred)
  ########
  sigma = 10
  var_names = paste0(paste0('V',1:(size - 1)), collapse = ',')
  initial.coxph.formula = formula(
    paste0('Surv(time,status) ~ ridge(',
           var_names,
           ', scale = FALSE, theta = ',
           1/sigma,')')
  )

  sm.pred <- sm.coxph(initial.coxph.formula, data)
  res <- rbind(res, sm.pred)
  rownames(res) <- NULL
  return(res)
}

simulation1 <- function(){
  library(parfm)
  library(dplyr)
  library(data.table)
  res <- c()
  for(i in 1:20){
    if(i%%10 == 0) print(i)
    res <- rbind(res, .simulation1(1, 50, 5))
  }

  res = data.table(res)
  res$alg <- rep(c('parfm','smcoxph'), 20)
  mean_res = res %>% group_by(alg) %>%
    summarise(
      m.sigma2 = mean(sigma2)
    )
  print(mean_res)
  return(res)
}

#------------------------------------------------------------------------------#
# The second simulation study looks at the performance of frailty variance
# estimates with three fixed effects.
#------------------------------------------------------------------------------#

.simulation2 <- function(sigma2 = 1, size = 50, rep = 5){
  #sigma2 = 1; size = 50; rep = 5
  res = c()
  data = rfrailty.transform(sigma2 = sigma2, size = size, rep = rep)
  data = data[order(data$time), ]
  ########
  parfm.fit <- parfm(
    Surv(time, status) ~ X1 + X2 + X3,
    cluster = 'group',
    data = data,
    dist = "exponential",
    frailty = 'lognormal'
  )
  res = rbind(res, parfm.fit[,1])
  ########
  # fP.fit <- frailtyPenal(
  #   Surv(time, status) ~ cluster(group),
  #   RandDist = 'LogN',
  #   data = data,
  #   hazard = "Weibull"
  # )
  # fP.pred <- c(fP.fit$sigma2, 1/fP.fit$scale.weib[1], fP.fit$coef)
  # res = rbind(res, fP.pred)
  ########
  sigma = 10
  var_names = paste0(paste0('V',1:(size - 1)), collapse = ',')
  initial.coxph.formula = formula(
    paste0('Surv(time,status) ~ ridge(',
           var_names,
           ', scale = FALSE, theta = ',
           1/sigma,')',' + X1 + X2 + X3')
  )

  sm.pred <- sm.coxph(initial.coxph.formula, data)
  res <- rbind(res, sm.pred)
  rownames(res) <- NULL
  return(res)
}

simulation2 <- function(){
  library(parfm)
  library(dplyr)
  library(data.table)
  res <- c()
  for(i in 1:20){
    if(i%%10 == 0) print(i)
    res <- rbind(res, .simulation2(1, 50, 5))
  }
  res = data.table(res)
  res$alg <- rep(c('parfm','smcoxph'), 20)
  mean_res = res %>% group_by(alg) %>%
    summarise(
      m.sigma2 = mean(sigma2),
      m.X1 = mean(X1),
      m.X2 = mean(X2),
      m.X3 = mean(X3)
    )
  print(mean_res)
  return(res)
}
