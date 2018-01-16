################################################################################
# sm.coxph
# Author: Yu Liu
# Date: Jan 2018
################################################################################


.parse_formula <- function(formula){
  if(class(formula) != 'formula'){
    stop("The formula should be a formula object.")
  }
  term_labels <- terms(formula) %>% attr(., which = "term.labels")
  is_tp <- startsWith(term_labels, prefix = 'tp(')
  tp_term <- term_labels[is_tp]
  if(length(tp_term) > 1) stop("Only one variable allowed in tp() format.")
  length_tp_term <- nchar(tp_term)
  tp_term_label <- substr(tp_term, 4, length_tp_term - 1)
  term_labels[is_tp] <- tp_term_label
  time_label <- all.vars(formula)[1]
  return(list(tp_term_label = tp_term_label,
              term_labels = term_labels,
              time_label = time_label))
}


.get_hazard_rate <- function(time, fit){
  unique_time_index <- cumsum(c(TRUE, diff(time) != 0))
  sfit <- survfit(fit, se.fit = FALSE)
  zcoef <- ifelse(is.na(coef(fit)), 0, coef(fit))
  offset <- sum(fit$means * zcoef)
  chaz <- sfit$cumhaz * exp(-offset)
  return(chaz[unique_time_index])
}


.create_surv_formula <- function(formula, tp_vars, sigma){
  tp_var_add_str <- paste(tp_vars, collapse = ' , ')
  tp_vars_formula_str <- paste0(
    "ridge(", tp_var_add_str, ", scale = FALSE, theta = ", 1/sigma, ')')
  term_labels <- SmoothCure:::.parse_formula(formula)$term_labels
  term_labels <- c(tp_vars_formula_str, term_labels)
  coxph_formula <- reformulate(
    termlabels = term_labels,
    response = formula[[2]]
  )
  return(coxph_formula)
}


.penalized_coxph <- function(formula, data, labels){
  fit <- coxph(formula, data = data, model = TRUE)
  coeff <- summary(fit)$coefficient[,1, drop = FALSE]
  b <- coeff[labels$ridge_tp_vars, , drop = FALSE]
  beta <- coeff[labels$term_labels, , drop = FALSE]
  hazard <- SmoothCure:::.get_hazard_rate(
    time = data[[labels$time_label]], fit = fit
  )
  loglik <- fit$loglik[2]
  return(list(b = b, beta = beta, hazard = hazard, loglik = loglik))
}


sm.coxph2 <- function(formula, data){

  labels <- SmoothCure:::.parse_formula(formula)
  data <- data[order(data[[labels$time_label]]),]

  tp_data <- tp(
    x = data[[labels$tp_term_label]],
    name = labels$tp_term_label,
    raw = FALSE
  )

  labels$tp_vars <- colnames(tp_data)
  labels$ridge_tp_vars <- paste0('ridge(', labels$tp_vars, ')')

  data <- cbind(data, tp_data)
  rm(tp_data)

  X <- data[,labels$term_labels] %>% as.matrix
  Z <- data[,labels$tp_vars] %>% as.matrix

  q <- length(labels$tp_vars) # Number of random effects
  ZxtZ <- apply(Z,1, function(x) x%*%t(x))
  Ik <- diag(rep(1,q))

  # Null Hypotheis Coxph
  coxph_formula <- reformulate(
    termlabels = labels$term_labels,
    response = formula[[2]]
  )
  coxph.fit <- coxph(coxph_formula, data)
  null_loglik <- coxph.fit$loglik[2]

  # Grid search to find the best initial value of sigma
  sigma_list <- c(1e-8, 1e-6, 1e-4, 1e-2, 1, 100, 1e4, 1e6)
  loglik_list <- c()
  for(sigma in sigma_list){
    surv_formula <- SmoothCure:::.create_surv_formula(
      formula = formula,
      tp_vars = labels$tp_vars,
      sigma = sigma
    )
    pfit <- SmoothCure:::.penalized_coxph(
      formula = surv_formula,
      data = data, labels = labels)
    hessian <- matrix(
      ZxtZ %*% (exp(X%*%pfit$beta + Z%*%pfit$b)*pfit$hazard),
      nrow = q, ncol = q
    ) + 1/sigma*Ik
    loglik_value <- pfit$loglik - log(det(hessian)) / 2 - q * log(sigma) /2
    loglik_list <- c(loglik_list, loglik_value)
  }

  sigma <- sigma_list[which.max(loglik_list)]
  print(paste0('Initial Sigma is ', as.character(sigma)))

  #
  old_sigma <- sigma
  relative_error = 1
  num_iterations = 1
  sigma_iterations = c(sigma)
  x_iterations = c()

  ##### Optimization Iteration #####
  while(relative_error > 1e-6){
    surv_formula <- SmoothCure:::.create_surv_formula(
      formula = formula,
      tp_vars = labels$tp_vars,
      sigma = sigma
    )
    pfit <- SmoothCure:::.penalized_coxph(
      formula = surv_formula,
      data = data, labels = labels)

    hessian <- matrix(
      ZxtZ %*% (exp(X%*%pfit$beta + Z%*%pfit$b)*pfit$hazard),
      nrow = q, ncol = q
    ) + 1/sigma*Ik

    sigma <- sum(pfit$b^2)/q + sum((1/eigen(hessian)$value))/q
    sigma_iterations <- c(sigma_iterations, sigma)

    relative_error <- abs(sqrt(sigma) - sqrt(old_sigma)) / abs(sqrt(old_sigma))
    old_sigma <- sigma
    num_iterations <- num_iterations + 1

    if(num_iterations %% 500 == 0){
      print("Reach 500 itertions")
      break
    }
  }

  loglik <- pfit$loglik - log(det(hessian)) / 2 - q * log(sigma) /2
  loglik_ratio <- max(0, 2 * (loglik - null_loglik))

  return(list(sigma = sigma,
              fit = pfit,
              loglik_ratio = loglik_ratio,
              sigma_iterations = sigma_iterations,
              num_iterations = num_iterations))
}

