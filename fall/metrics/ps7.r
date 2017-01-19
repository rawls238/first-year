library("dplyr")

min_distance <- function(beta_hat, mat, R, c) {
  tR = t(R)
  constraint = tR %*% beta_hat - c
  mid <- solve(tR %*% mat %*% R)
  return (beta_hat - mat %*% R %*% mid %*% constraint)
}

calculate_v <- function(X, XX, e, n, k) {
  mid_term <- t(X) %*%diag(e^2)%*% X
  return(n/(n-k)*XX %*% mid_term %*% XX)
}

data <- read.table("/Users/garidor/Desktop/cps09mar.txt")
data <- dplyr::rename(data, age = V1, female = V2, hisp = V3, educ = V4, earnings = V5,
                      hours_worked = V6, weeks_worked = V7, union = V8, uncov = V9, region = V10, race = V11, marital = V12)
data2 <- filter(data, race == 1, hisp == 1, female == 0)
educ <- data2$educ
experience <- data2$age - data2$educ - 6
exp_2 <- experience^2/100
lwage <- log(data2$earnings / (data2$hours_worked * data2$weeks_worked))
dummy_married_1 <- as.numeric(data2$marital == 1)
dummy_married_2 <- as.numeric(data2$marital == 2)
dummy_married_3 <- as.numeric(data2$marital == 3)
dummy_widowed <- as.numeric(data2$marital == 4)
dummy_divorced <- as.numeric(data2$marital == 5)
dummy_separated <- as.numeric(data2$marital == 6)
model1 <- lm(lwage ~ 1 + educ + experience + exp_2 + dummy_married_1 + dummy_married_2 + dummy_married_3 + dummy_widowed + dummy_divorced + dummy_separated)
print(summary(model1))

X <- cbind(matrix(1,length(data2$educ),1), educ, experience, exp_2, dummy_married_1, dummy_married_2, dummy_married_3, dummy_widowed, dummy_divorced, dummy_separated)
constraint <- cbind(c(0, 0, 0, 0, 1, 0, 0, -1, 0, 0), c(0, 0, 0, 0, 0, 0, 0, 0, 1, -1))
c <- c(0, 0)
XX <- solve(t(X)%*%X)
cls_est <- matrix(ncol=2, nrow=10)
cls_est[, 1] <- min_distance(model1$coefficients, XX, constraint, c)
e_tilda <- lwage - X %*% as.matrix(cls_est[, 1])
n <- length(data2$educ)
k <- 10
cls_est[, 2] <- sqrt(diag(calculate_v(X, XX, as.vector(e_tilda), n, k)))
print(cls_est)

emd_est <- matrix(ncol=2, nrow=10)
V_hat <- calculate_v(X, XX, as.vector(model1$residuals), n, k)
emd_est[,1] <- min_distance(model1$coefficients, V_hat, constraint, c)
covar_matrix <- V_hat - V_hat %*% constraint %*% solve(t(constraint) %*% V_hat %*% constraint) %*% t(constraint) %*% V_hat
emd_est[,2] <- sqrt(diag(covar_matrix))
print(emd_est)