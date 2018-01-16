library(dplyr)
library(survival)
################################################################################
# Data Generator
################################################################################
generate_data <- function(n, f, c_f){
  x <- runif(n = n, min = 0, max = 1)
  rexp_1 <- rexp(n = n, rate = 1)
  uncensored_time <- exp(-f(x)) * rexp_1
  censor <- sapply(x, c_f)
  data <- data.frame(x = x)
  data$time <- pmin(uncensored_time, censor)
  data$status <- as.integer(uncensored_time <= censor)
  data$sex <- rbinom(n = n, size = 1, prob = 0.5)
  return(data)
}

################################################################################
# Single simulation Wrapper Function
################################################################################
single_simulation <- function(n, f, c_f, return_raw_fit = TRUE){
  data <- generate_data(
    n = 1000, f = f, c_f = c_f)
  print(paste0(
    "Censor Rate is ",
    as.character(1 - sum(data$status) / nrow(data)
    )))
  data$x <- data$x
  formula <- Surv(time, status) ~ tp(x) + sex
  sm.fit <- sm.coxph2(formula, data)

  if(!return_raw_fit){
    loglik_ratio <- sm.fit$loglik_ratio
    sigma <- sm.fit$sigma
    beta_x <- sm.fit$fit$beta['x',]
    beta_sex <- sm.fit$fit$beta['sex', ]

    sm.fit <- c(loglik_ratio, sigma, beta_x, beta_sex)
    names(sm.fit) <- c('loglik_ratio', 'sigma', 'x', 'sex')
  } else {
    sm.fit$x <- data$x
  }
  return(sm.fit)
}

###############################################################################
# Single simulation f(x) = 2x + sin(30 * x)
################################################################################
f <- function(x) 2 * x + sin(30 * x)
c_f <- function(x) 2
set.seed(1)
fit <- single_simulation(
  n = 1000, f = f, c_f = c_f)

# num_iterations = 16
fit$num_iterations

# sigma = 796.3493
fit$sigma

# slope = 29.00818
fit$fit$beta['x',]

# plot tp(age)
x <- sort(unique(fit$x))
tp_matrix <- tp(x = x, name = 'x') %>% as.matrix
y <- tp_matrix %*% fit$fit$b + x * fit$fit$beta['x',]
plot(x, f(x), cex = 0.2, col = 'blue')
lines(x, y, col = 'red')

# log likelihood
fit$loglik_ratio # 337.0645

###############################################################################
# Single simulation f(x) = 2x
################################################################################
f <- function(x) 2 * x
c_f <- function(x) 2
set.seed(1)
fit <- single_simulation(
  n = 1000, f = f, c_f = c_f)

# num_iterations = 500
fit$num_iterations

# sigma = 9.999992e-09
fit$sigma

# slope = 2.140525
fit$fit$beta['x',]

plot(log(fit$sigma_iterations))

# plot tp(age)
x <- sort(unique(fit$x))
tp_matrix <- tp(x = x, name = 'x') %>% as.matrix
y <- tp_matrix %*% fit$fit$b + x * fit$fit$beta['x',]
plot(x, f(x), cex = 0.2, col = 'blue')
lines(x, y, col = 'red')

# log likelihood
fit$loglik_ratio # 6.542523e-07


################################################################################
# Monte Carlo Simulation Alternative f(x) = 2x + sin(30x)
################################################################################
file_name <- "./simulation/Alter_2x_sin30x_c2.csv"
if(!file.exists(file_name)){
  write('loglik_ratio,sigma,x,sex',
        file = file_name,
        append = TRUE)
}

f <- function(x) 2 * x + sin(30 * x)
c_f <- function(x) 2
for(i in 1:100){
  simulation_null <- single_simulation(
    n = 1000, f = f, c_f = c_f, return_raw_fit = FALSE)
  text_line <- paste0(simulation_null, collapse = ',')
  write(text_line,
        file = file_name,
        append = TRUE)
}

simulation_result <- read.csv(file_name)
simulation_result$loglik_ratio %>% round(5)
plot(density(simulation_result$sigma))

################################################################################
# Monte Carlo Simulation Null f(x) = 2x
################################################################################
file_name <- "./simulation/Null_2x_c2.csv"
if(!file.exists(file_name)){
  write('loglik_ratio,sigma,x,sex',
        file = file_name,
        append = TRUE)
}

f <- function(x) 2 * x
c_f <- function(x) 2
for(i in 1:1000){
  simulation_null <- single_simulation(
    n = 1000, f = f, c_f = c_f, return_raw_fit = FALSE)
  text_line <- paste0(simulation_null, collapse = ',')
  write(text_line,
        file = file_name,
        append = TRUE)
}

simulation_result <- read.csv(file_name)
simulation_result$loglik_ratio %>% round(5)
plot(density(simulation_result$sigma))

################################################################################
# Monte Carlo Simulation Null f(x) = 2x with c = rnorm(2, 0.5)
################################################################################
file_name <- "./simulation/Null_2x_norm2_0.5.csv"
if(!file.exists(file_name)){
  write('loglik_ratio,sigma,x,sex',
        file = file_name,
        append = TRUE)
}

f <- function(x) 2 * x
c_f <- function(x) rnorm(1, 2, 0.5)
for(i in 1:1000){
  simulation_null <- single_simulation(
    n = 1000, f = f, c_f = c_f, return_raw_fit = FALSE)
  text_line <- paste0(simulation_null, collapse = ',')
  write(text_line,
        file = file_name,
        append = TRUE)
}

simulation_result <- read.csv(file_name)
simulation_result$loglik_ratio %>% round(5)
plot(density(simulation_result$sigma))

################################################################################
# Monte Carlo Simulation Null f(x) = 2x with c = rnorm(0.4, 0.5)
################################################################################
file_name <- "./simulation/Null_2x_norm0.4_0.5.csv"
if(!file.exists(file_name)){
  write('loglik_ratio,sigma,x,sex',
        file = file_name,
        append = TRUE)
}

f <- function(x) 2 * x
c_f <- function(x) rnorm(1, 0.4, 0.5)
for(i in 1:1000){
  simulation_null <- single_simulation(
    n = 1000, f = f, c_f = c_f, return_raw_fit = FALSE)
  text_line <- paste0(simulation_null, collapse = ',')
  write(text_line,
        file = file_name,
        append = TRUE)
}

simulation_result <- read.csv(file_name)
simulation_result$loglik_ratio %>% round(5)
plot(density(simulation_result$sigma))


################################################################################
# Monte Carlo Simulation Alternative f(x) = 2x + sin(30x)
# with c = rnorm(1.4, 0.5)
################################################################################
file_name <- "./simulation/Alter_2x_sin30x_norm1.4_0.5.csv"
if(!file.exists(file_name)){
  write('loglik_ratio,sigma,x,sex',
        file = file_name,
        append = TRUE)
}

f <- function(x) 2 * x + sin(30 * x)
c_f <- function(x) rnorm(1, 1.4, 0.5)
for(i in 1:100){
  simulation_null <- single_simulation(
    n = 1000, f = f, c_f = c_f, return_raw_fit = FALSE)
  text_line <- paste0(simulation_null, collapse = ',')
  write(text_line,
        file = file_name,
        append = TRUE)
}

simulation_result <- read.csv(file_name)
simulation_result$loglik_ratio %>% round(5)
plot(density(simulation_result$sigma))

################################################################################
# Monte Carlo Simulation Alternative f(x) = 2x + sin(30x)
# with c = rnorm(0.8, 0.5)
################################################################################
file_name <- "./simulation/Alter_2x_sin30x_norm0.8_0.5.csv"
if(!file.exists(file_name)){
  write('loglik_ratio,sigma,x,sex',
        file = file_name,
        append = TRUE)
}

f <- function(x) 2 * x + sin(30 * x)
c_f <- function(x) rnorm(1, 0.8, 0.5)
for(i in 1:100){
  simulation_null <- single_simulation(
    n = 1000, f = f, c_f = c_f, return_raw_fit = FALSE)
  text_line <- paste0(simulation_null, collapse = ',')
  write(text_line,
        file = file_name,
        append = TRUE)
}

simulation_result <- read.csv(file_name)
simulation_result$loglik_ratio %>% round(5)
plot(density(simulation_result$sigma))

################################################################################
# Monte Carlo Simulation Alternative f(x) = 2x + sin(30x)
# with c = rnorm(0.4, 0.5)
################################################################################
file_name <- "./simulation/Alter_2x_sin30x_norm0.4_0.5.csv"
if(!file.exists(file_name)){
  write('loglik_ratio,sigma,x,sex',
        file = file_name,
        append = TRUE)
}

f <- function(x) 2 * x + sin(30 * x)
c_f <- function(x) rnorm(1, 0.4, 0.5)
for(i in 1:100){
  simulation_null <- single_simulation(
    n = 1000, f = f, c_f = c_f, return_raw_fit = FALSE)
  text_line <- paste0(simulation_null, collapse = ',')
  write(text_line,
        file = file_name,
        append = TRUE)
}

simulation_result <- read.csv(file_name)
simulation_result$loglik_ratio %>% round(5)
plot(density(simulation_result$sigma))
