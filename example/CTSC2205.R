library(dplyr)
library(survival)

################################################################################
# Read Raw Data
################################################################################
raw.data <- read.csv(
  '~/Downloads/STSwActualSizeCTSC2205.csv',
  stringsAsFactors = FALSE)

################################################################################
# Prepare Data Frame with All numeric variables
################################################################################
data <- data.frame(age = raw.data$ageatdiagnosis / 200)
data$status <- as.integer(raw.data$vitalstatus == 'Alive')
data$time <- raw.data$survivalmonths
data$radiotherapy <- as.integer(raw.data$radiationyesno == 'Yes')
data$sex <- as.integer(raw.data$sex == 'Male')
data$racewhite <- as.integer(raw.data$race == 'White')
# site
data$siteHeadNeck <- as.integer(raw.data$site == 'Head/Neck')
data$siteTrunk <- as.integer(raw.data$site == 'Trunk')
data$siteOther <- as.integer(
  raw.data$site == 'Visceral' |
  raw.data$site == 'Thoracic' |
  raw.data$site == 'GU/GYN')
# tumor size
data$size10_15CM <- as.integer(raw.data$nsize == "10-15 CM")
data$size15_CM <- as.integer(raw.data$nsize == '15+ CM')
data$size5_10CM <- as.integer(raw.data$nsize == '5-10 CM')

################################################################################
# Fit Smooth Coxph
################################################################################
formula <- Surv(time, status) ~
  tp(age) + sex + radiotherapy + racewhite + siteHeadNeck +
  siteOther + siteTrunk + size5_10CM + size15_CM + size10_15CM

sm.fit <- sm.coxph2(formula, data)

# num_iterations = 24
sm.fit$num_iterations

#
plot(sm.fit$sigma_iterations, type = 'l')

# sigma = 1.822226
sm.fit$sigma

# log likelihood ratio = 68.65589
sm.fit$loglik_ratio

################################################################################
# plot tp(age)
################################################################################
x <- sort(unique(data$age))
tp_matrix <- tp(x = x, name = 'age') %>% as.matrix

y <- tp_matrix %*% sm.fit$fit$b + x * sm.fit$fit$beta['age',]
plot(x, y, type = 'l')





