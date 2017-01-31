get_mean <- function(n, mu) {
  return(sqrt(n) * mu)
}

simulate <- function(num_iter, m, n, mu=0.01, var=1.0, c=1.0) {
  s1 <- c()
  s2 <- c()
  for (year in 1:num_iter) {
    for(i in seq(n, 10000, by = 1000)) {
      x <- rnorm(i, m/sqrt(i), var)
      z <- rnorm(i, mu, var)
      c_n <- sqrt(i) * mu / m
      stat_1 <- sqrt(i) * (mean(z) / mean(x) - c_n)
      stat_2 <- sqrt(i) * (mean(z) - c_n * mean(x))
      err_1 <- sqrt((var(x)*mean(z)^2)/mean(x)^4+var(z)/(c_n^2 + mean(x)^2))
      err_2 <- sqrt(var(z) + c_n^2 * var(x))
      t_stat_1 <- stat_1 / err_1
      t_stat_2 <- stat_2 / err_2
      if (abs(t_stat_1) > 1.96) {
        if (is.null(s1) || is.na(s1[i])) {
          s1[i] <- 1
        } else {
          s1[i] <- s1[i] + 1
        }
      }
      if (abs(t_stat_2) > 1.96) {
        if (is.null(s2) || is.na(s2[i])) {
          s2[i] <- 1
        } else {
          s2[i] <- s2[i] + 1
        }
      }
    }
  }
  return(c(s1, s2))
}

num_iter <- 10000
var <- 1
c <- 1
n <- 1000
mu <- 0.01
m <- get_mean(n, mu)
b <- simulate(num_iter, m, n, mu)
