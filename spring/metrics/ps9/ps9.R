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
gamma <- matrix(0, nrow = n, ncol = n) #????
gamma[1,] <- rep(1, n)
num_datasets <- 1000
for(l in l_list) {
  for (i in seq(1, num_datasets)) {
    v <- rnorm(l)
    eta <- rnorm(l)
    u <- v + eta
    Z <- matrix(rnorm(n * l), n, l)
    x <- t(Z) %*% gamma + v
    y <- u
    ols <- lm(y ~ x)
    nodata <- data.frame(x= numeric(0), y= numeric(0), z = numeric(0))
    tsls_est <- tsls(y ~ x, instruments=Z)
  }
}