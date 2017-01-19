library("dplyr")

# 3.21.1
data <- read.table("/Users/garidor/Desktop/cps09mar.txt")
data <- dplyr::rename(data, age = V1, female = V2, hisp = V3, educ = V4, earnings = V5,
       hours_worked = V6, weeks_worked = V7, union = V8, uncov = V9, region = V10, race = V11, marital = V12)
data2 <- filter(data, race == 1, hisp == 1, female == 0)
experience <- data2$age - data2$educ - 6
exp_2 <- experience^2
lwage <- log(data2$earnings / (data2$hours_worked * data2$weeks_worked))
dummy_married <- as.numeric(data2$marital <= 3)
dummy_widowed <- as.numeric(data2$marital == 4) + as.numeric(data2$marital == 5)
dummy_widowed <- as.numeric(dummy_widowed > 0)
dummy_separated <- as.numeric(data2$marital == 6)
dummy_ne <- as.numeric(data2$region == 1)
dummy_s <- as.numeric(data2$region == 3)
dummy_w <- as.numeric(data2$region == 4)
model1 <- lm(lwage ~ data2$educ + experience + exp_2 + dummy_married + dummy_widowed + dummy_separated + dummy_ne + dummy_s + dummy_w)
print(summary(model1))

calculate_v_bar <- function(X, model) {
  XX <- solve(t(X)%*%X)
  k <- ncol(X)
  out <- lm.influence(model)
  ehat <- model$residuals
  hii<-out$hat
  V_bar<-XX%*%t(X)%*%diag(ehat^2/(1-hii))%*%X%*%XX
  return(V_bar)
}

#4.12
X <- cbind(matrix(1,length(data2$educ),1), data2$educ, experience, exp_2, dummy_married, dummy_widowed, dummy_separated, dummy_ne, dummy_s, dummy_w)
cat("Calculate standard error using Horn-Horn-Duncan")
a <- sqrt(diag(calculate_v_bar(X, model1)))
print(a)
cat("\n")

#6.11
exp_2 <- experience^2/100
model2 <- lm(lwage ~ data2$educ + experience + exp_2)
print(summary(model2))
X <- cbind(matrix(1, length(data2$educ), 1), data2$educ, experience, exp_2)
V_bar <- calculate_v_bar(X, model2)
cat("Standard errors ")
print(sqrt(diag(V_bar)))
beta_1 <- summary(model2)$coefficients[2]
beta_2 <- summary(model2)$coefficients[3]
beta_3 <- summary(model2)$coefficients[4]
val <- beta_1 / (beta_2 + (beta_3 * experience / 50))
plot(experience, val, ylim=c(-50, 50))

RT <- cbind(1/(beta_2 + (experience *  beta_3 / 50)), -1 * beta_1 / (beta_2 + (experience *  beta_3 / 50)^2), -50 * beta_1 / (50*beta_2 + (experience *  beta_3)^2))
R <- t(RT)
s_vec <- c()
for (i in 1:length(experience)) {
  s_vec <- c(s_vec, sqrt(R[, i] %*% V_bar[2:4, 2:4] %*% R[, i]))
}
plot(experience, s_vec, ylim=c(0, 10))

m_e <- mean(experience)
RT <- cbind(1/(beta_2 + (m_e *  beta_3 / 50)), -1 * beta_1 / (beta_2 + (m_e *  beta_3 / 50)^2), -50 * beta_1 / (50*beta_2 + (m_e *  beta_3)^2))
R <- t(RT)
s_bar <- sqrt(RT %*% V_bar[2:4, 2:4] %*% R)
cat("s_hat(theta_hat) ", s_bar, "\n")

# Evaluate the CI at the mean level of experience
theta_bar <- beta_1 / (beta_2 + beta_3 * m_e / 50)
ci_lower <- theta_bar - 1.645 * s_bar
ci_upper <- theta_bar + 1.645 * s_bar
cat("90% Confidence interval ", ci_lower, " to ", ci_upper, "\n")

zt <- as.matrix(c(1, 12, 20, 4))
z <- t(zt)
ri <- 1.96 * sqrt(z %*% V_bar %*% zt)
z_beta_hat <- z %*% as.matrix(model2$coefficients)
cat("95% Confidence interval for regression function ", z_beta_hat - ri, " to ", z_beta_hat + ri, "\n")

zt <- as.matrix(c(1,16,5,1/4))
z <- t(zt)
z_beta_hat <- z %*% as.matrix(model2$coefficients)
zi <- 1.28 * sqrt(z %*% V_bar %*% zt)
cat("80% Forecast interval for log wage", z_beta_hat - zi, " to ", z_beta_hat + zi, "\n")
cat("80% Forecast interval for wage", z_beta_hat - exp(zi), " to ", z_beta_hat + exp(zi))