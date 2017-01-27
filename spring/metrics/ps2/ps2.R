get_mean <- function(n, mu) {
  return(sqrt(n) * mu)
}

simulate <- function(num_iter, n, mu, var, c) {
  pass_stat_1 <- 0
  pass_stat_2 <- 0
  for (year in 1:num_iter) {
    x <- rnorm(n, mu, var)
    z <- rnorm(n, mu, var)
    stat_1 <- sqrt(n) * (mean(z) / mean(x) - c)
    stat_2 <- sqrt(n) * (mean(z) - c * mean(x))
    err_1 <- sqrt((var(x)*mean(z)^2)/mean(x)^4+var(z)/(c^2 + mean(x)^2))
    err_2 <- sqrt(var(z) + c^2 * var(x))
    t_stat_1 <- stat_1 / err_1
    t_stat_2 <- stat_2 / err_2
    if (abs(t_stat_1) > 1.96) {
      pass_stat_1 <- pass_stat_1 + 1
    }
    if (abs(t_stat_2) > 1.96) {
      pass_stat_2 <- pass_stat_2 + 1
    }
  }
  cat(pass_stat_1 / num_iter, pass_stat_2 / num_iter, '\n')
}

n <- 100
num_iter <- 10000
var <- 1
c <- 1
cat('mu = 1, c = 1, ')
simulate(n, num_iter, 1, var, c)
cat('mu = .1, c = 1, ')
simulate(n, num_iter, .1, var, c)
cat('mu = .01, c = 1, ')
simulate(n, num_iter, .01, var, c)

c <- 2
cat('mu = 1, c = 2, ')
simulate(n, num_iter, 1, var, c)
cat('mu = .1, c = 2, ')
simulate(n, num_iter, .1, var, c)
cat('mu = .01, c = 2, ')
simulate(n, num_iter, .01, var, c)