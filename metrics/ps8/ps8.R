library("dplyr")
library("sandwich")

calculate_wald_stat <- function(R, V_hat, theta_0, coefficients) {
  mid <- solve(t(R) %*% V_hat %*% R)
  side <- t(R) %*% coefficients - theta_0
  return (t(side) %*% mid %*% side)
}

min_distance <- function(beta_hat, mat, R, c) {
  tR = t(R)
  constraint = tR %*% beta_hat - c
  mid <- solve(tR %*% mat %*% R)
  return (beta_hat - mat %*% R %*% mid %*% constraint)
}

#8.3
invest_data <- read.table("/Users/garidor/Desktop/first-year/metrics/ps8/invest.dat")
invest_data <- dplyr::rename(invest_data, I = V1, Q = V2, C = V3, D = V4);

#8.3.a
model1 <- lm(invest_data$I ~ 0 + invest_data$Q + invest_data$C + invest_data$D)
print(summary(model1))

#8.3.b
coefficients <- coef(summary(model1))[,1]
X <- cbind(invest_data$Q, invest_data$I, invest_data$C)
XX <- solve(t(X)%*%X)
n <- length(invest_data$I)
V_hat <- vcovHC(model1,type=c("HC0"))
std_errors <- sqrt(diag(V_hat))
ci <- matrix(nrow = length(coefficients), ncol=2)
ci[,1] <- coefficients - 1.96 * std_errors
ci[,2] <- coefficients + 1.96 * std_errors
print("95% confidence interval ")
print(ci)

#8.3.c
R <- cbind(c(0, 1, 0), c(0, 0, 1))

# test joint hypothesis of C and I = 0
wald <- calculate_wald_stat(R, V_hat, c(0, 0), coefficients)
cv <- qchisq(0.95, df=2)
cat("8.3.c Wald stat ", abs(wald), " with critical value of ", cv, "\n")
t <- coefficients[2] / std_errors[2]
cv <- 1.96
cat("8.3.c T stat ", abs(t), " with critical value of ", cv, "\n")

#8.3.d
I <- invest_data$I
Q <- invest_data$Q
D <- invest_data$D
C <- invest_data$C
n <- length(invest_data$I)
model2 <- lm(I ~ 0 + Q + C + D + I(Q^2) + I(C^2) + I(D^2) + I(Q*C) + I(Q*D) + I(C*D))
print(summary(model2))
R1 <- c(0, 0, 0, 1, 0, 0, 0, 0, 0)
R2 <- c(0, 0, 0, 0, 1, 0, 0, 0, 0)
R3 <- c(0, 0, 0, 0, 0, 1, 0, 0, 0)
R4 <- c(0, 0, 0, 0, 0, 0, 1, 0, 0)
R5 <- c(0, 0, 0, 0, 0, 0, 0, 1, 0)
R6 <- c(0, 0, 0, 0, 0, 0, 0, 0, 1)
R <- cbind(R1, R2, R3, R4, R5, R6)
V_hat <- vcovHC(model2,type=c("HC0"))
wald <- calculate_wald_stat(R, V_hat, c(0, 0, 0, 0, 0, 0), model2$coefficients)
cv <- qchisq(0.95, df=6)
cat("Wald stat ", wald, " with rejection value ", cv, "\n")

#8.4.a
nerlov_data <- read.table("/Users/garidor/Desktop/first-year/metrics/ps8/nerlov.dat")
nerlov_data <- nerlov_data[-c(146)] # some random data import issue
nerlov_data <- dplyr::rename(nerlov_data, TC = V1, Q = V2, PL = V3, PF = V4, PK = V5)
log_tc <- log(nerlov_data$TC)
log_q <- log(nerlov_data$Q)
log_pl <- log(nerlov_data$PL)
log_pk <- log(nerlov_data$PK)
log_pf <- log(nerlov_data$PF)
n <- length(log_tc)
model3 <- lm(log_tc ~ 1 + log_q + log_pl + log_pk + log_pf)
print(summary(model3))
V_hat_ols <- vcovHC(model3,type=c("HC0"))

#8.4.c
X <- cbind(matrix(1,n,1), log_q, log_pl, log_pk, log_pf)
XX <- solve(t(X) %*% X)
k <- 5
constraint <- as.matrix(c(0, 0, 1, 1, 1))
c_val <- c(1)
cls_est <- matrix(ncol=2, nrow=length(model3$coefficients))
cls_est[,1] <- min_distance(model3$coefficients, XX, constraint, c_val)
tmp <- diag(k) - XX %*% constraint %*% solve(t(constraint) %*% XX %*% constraint) %*% t(constraint)
cls_est[,2] <- sqrt(diag(tmp %*% V_hat_ols %*% t(tmp)))
print("CLS estimates and standard errors 8.4.c")
print(cls_est)

#8.4.d
emd_est <- matrix(ncol=2, nrow=length(model3$coefficients))
emd_est[,1] <- min_distance(model3$coefficients, V_hat_ols, constraint, c_val)
covar_matrix <- V_hat_ols - V_hat_ols %*% constraint %*% solve(t(constraint) %*% V_hat_ols %*% constraint) %*% t(constraint) %*% V_hat_ols
emd_est[,2] <- sqrt(diag(covar_matrix))
print("EMD estimates and standard errors 8.4.d")
print(emd_est)

#8.4.e
wald_stat <- calculate_wald_stat(constraint, V_hat_ols, c_val, model3$coefficients)
cv <- qchisq(0.95, df=3)
cat("Wald stat for 8.4.e is ", wald_stat, " with critical value ", cv, "\n")

#8.4.f
min_distance_stat <- t(as.matrix(model3$coefficients) - emd_est[, 1]) %*% solve(V_hat_ols) %*% (model3$coefficients - emd_est[, 1])
cat("Min distance stat for 8.4.f is ", min_distance_stat)