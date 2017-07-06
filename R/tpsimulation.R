################################################################################
# Simulation
# sm.coxph on Truncated Piecewise Coxph Data
# Author: Yu Liu
# Date: June 2017
################################################################################
simulation3 <-function(){
  library(data.table)
  library(dplyr)
  library(survival)
  m = c()
  for(j in seq(0.4,2,0.1)){
    print(paste0("Running True Sigma2 = ", j))
    for(i in 1:200){
      data = rweibull.age(j)
      formula1 = Surv(time, status) ~ ridge(age1,age2,age3,age4,age5,age6,age7,age8,age9,age10, scale = FALSE, theta = 1/10) + age0
      res <- c(sm.coxph(formula = formula1, data = data),j)
      names(res) <- c('sigma2','lambda','age','true_sigma2')
      m <- rbind(m, res)
    }
  }
  m <- data.table(m)

  summary.m <- m %>% group_by(true_sigma2) %>%
    summarise(
      m.sigma2 = mean(sigma2),
      m.age = mean(age)
    )
  return(summary.m)
}
