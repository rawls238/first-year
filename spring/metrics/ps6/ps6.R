library("quantreg")

data("CPS1985", package = "AER")
cps <- CPS1985
model <- rq(log(wage) ~ education + experience, data=cps,tau=seq(0.01, 0.99, by=0.01))
plot(summary(model))


#Q2
#modified from SO
top <- function(x, d, n=100){
  result <- numeric()
  for(i in 1:n){
    j <- which.max(x)
    result[i] <- d[j]
    x[j] <- -Inf
  }
  result
}

num_val <- 1000
top_10 <- num_val / 10;
d <- rnorm(num_val)
y <- c()
for (i in seq(0, num_val)) {
  dev <- sqrt((1 + d[i])^2)
  y <- c(y, rnorm(1, 0, dev))
}
vals <- tail(sort(y), top_10)
x_vals <- top(y, d, top_10)

X <- cbind(matrix(1,length(x_vals),1), x_vals)
XX <- solve(t(X)%*%X)
ols_est <- XX %*% (t(X) %*% vals)

quantile_model <- rq(vals ~ x_vals, tau=0.9)


