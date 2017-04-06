library(sem)
library(foreign)

data <- read.dta("/Users/garidor/Desktop/first-year/spring/metrics/ps9/fertility.dta")
ols <- lm(samesex ~ morekids, data = data)

t <- tsls(weeksm1 ~ morekids, instruments=~samesex, data = data)

weekly_earning <- 835.92
intercept <- coef(t)['(Intercept)']
res <- coef(t)['morekids']
two_kids <- intercept * weekly_earning
three_kids <- (intercept + res) * weekly_earning
percent_change <- 100 * (1 - (three_kids / two_kids))

ols_res <- lm(weeksm1 ~ morekids, data = data)
print(summary(t))
print(summary(ols_res))


set.seed(1)
n <- 200
l_list <- c(1, 2, 3, 10, 20, 50, 150, 175, 200)
num_datasets <- 1000
average_bias <- c()
for(l in l_list) {
  bias <- c()
  for (i in seq(1, num_datasets)) {
    gamma <- rep(1, l)
    gamma[1] <- 1
    v <- rnorm(n)
    eta <- rnorm(n)
    u <- v + eta
    Z <- matrix(rnorm(n * l), n, l)
    x <- Z %*% gamma + v
    y <- as.matrix(u)
    x <- as.matrix(x)
    ols <- lm(y ~ x)
    ols_est <- coef(ols)['x']
    P <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
    tsls_est <- solve(t(x) %*% P %*% x) %*% t(x) %*% P %*% y
    bias <- c(bias, tsls_est - ols_est)
  }
  average_bias <- c(average_bias, mean(bias))
}
plot(x=l_list,y=average_bias)