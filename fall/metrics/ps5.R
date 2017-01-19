x1 <- rnorm(100, 0, 1)
x2 <- rnorm(100, 0, 1)
e <- rnorm(100, 0, 1)

generate_y <- function(x1, x2, e) {
  beta0 <- 5
  beta1 <- 1
  beta2 <- 2
  return(beta0 + beta1 * x1 + beta2 * x2 + e)
}

calculated_r <- function(model) {
  e_sum <- sum(model$residuals^2)
  y1 <- model$fitted + model$residuals
  mean <- mean(y1)
  denom <- sum((y1 - mean)^2)
  return(1 - (e_sum/denom))
}

y <- generate_y(x1, x2, e)
model1 <- lm(y ~ x1 + x2)
model2 <- lm(y ~ 0 + x1 + x2)
print(summary(model1))
print(summary(model2))
cat("Calculated R^2 for the non-intercept model is ", calculated_r(model2), "\n")
cat("Calculated R^2 for the intercept model is ", calculated_r(model1), "\n")
k1 <- sum(x1 * model2$residuals)
k2 <- sum(x2 * model2$residuals)
cat("1e for k = 1", k1, "\n")
cat("1e for k = 2", k2, "\n")

model3 <- lm(y ~ x1 + x2)
cat("2a.1: ", sum(model3$residuals), "\n")
cat("2a.2: ", sum(x1 * model3$residuals), "\n")
cat("2a.3: ", sum(x2 * model3$residuals), "\n")
cat("2a.4, k = 1: ", sum(x1^2 * model3$residuals), "\n")
cat("2a.4, k = 2: ", sum(x2^2 * model3$residuals), "\n")
cat("2a.5: ", sum(model3$fitted * model3$residuals), "\n")

out <- lm.influence(model3)
cat("sum of hii is 3 is", all.equal(sum(out$hat), 3), "\n")
ei_tilda <- (1/(1-out$hat))*model3$residuals
ei_tilda_mean <- mean(ei_tilda)
cat("2c.1: ", ei_tilda_mean, "\n")
ei_tilda_hat <- mean(model3$residuals)
cat("2c.2: ", ei_tilda_hat, "\n")


X <- matrix(nrow=3, ncol=length(x1))
X[1, ] <- rep(1, 100)
X[2, ] <- x1
X[3, ] <- x2
XT <- t(X)
XXT <- X %*% XT
XXT_inv <- solve(XXT)
e_hat <- model3$residuals
k <- sum(out$hat)
n <- length(x1)
s_squared <- (1/(n-k))*sum(e_hat^2)
V0_hat <- XXT_inv * s_squared
cat("2d.1: ")
print(V0_hat)
cat("\n")
mid_term <- X %*%diag(e_hat^2)%*% XT
V_hat <- n/(n-k)*XXT_inv %*% mid_term %*% XXT_inv
cat("2d.2: ")
print(V_hat)
cat("\n")
e_bar <- (1-out$hat)^(-1/2)*e_hat
mid_term <- X %*%diag(e_bar^2)%*% XT
V_bar <- XXT_inv %*% mid_term %*% XXT_inv
cat("2d.3: ")
print(V_bar)
cat("\n")
mid_term <- X %*%diag(ei_tilda^2)%*% XT
V_tilda <- XXT_inv %*% mid_term %*% XXT_inv
cat("2d.4: ")
print(V_tilda)
cat("\n")
  
x3 <- rnorm(100, 0, 1)
model4 <- lm(y ~ x1 + x2 + x3)
print(summary(model4))
x2[length(x2)] = 10
model5 <- lm(y ~ x1 + x2 + x3)
out5 <- lm.influence(model5)
ei_tilda <- (1/(1-out5$hat))*model5$residuals
cat("2f: The maximum is attained at observation", which.max(abs(out5$hat * ei_tilda)))