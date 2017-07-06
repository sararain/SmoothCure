################################################################################
# Generate Truncated Piecewised Coxph Data
# Author: Yu Liu
# Date: June 2017
################################################################################
rweibull.age <- function(sigma2 = 1, size = 1000, slope = 1,
                         nknots = 10, knots = NULL,
                         lambda = 0.1, cencor.rate = 0.8){
  age <- sample(40:90, size = size, replace = TRUE)/5
  theta <- c(slope, rnorm(nknots, sd = sqrt(sigma2)))
  XZ <- tp(x = age, knots = knots, nknots = nknots, include_raw = TRUE)
  surv.data <- as.data.frame(XZ)
  errors <- rexp(size, rate = lambda)
  surv.data$time <- exp(- XZ %*% theta)*errors
  surv.data$status <- rbinom(size, 1, cencor.rate)
  return(surv.data)
}
