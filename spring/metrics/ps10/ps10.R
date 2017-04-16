library(gmm)
library(sem)
library(parallel)
library(ivmodel)

g1 <- function(beta,data) {
  d <- dim(data)[2]
  y <- data[,d]
  x <- data[,d-1]
  z <- data[,1:(d-2)]
  return(z * (y - x*beta))
}

compute_vals <- function(x, n, l, cv) {
  gamma <- rep(0, l)
  gamma[1] <- 1
  v <- rnorm(n)
  eta <- rnorm(n)
  u <- v + eta
  Z <- matrix(rnorm(n * l), n, l)
  x <- Z %*% gamma + v
  y <- as.matrix(u)
  x <- as.matrix(x)
  
  mod <- ivmodel(y, x, Z, intercept=FALSE, beta0=0)
  cof <- coef(mod)
  ols_est <- cof[1,2]
  ols_t_stat_0 <- abs(cof[1,4]) < cv
  tsls_est <- cof[2,2]
  tsls_t_stat_0 <- abs(cof[2,4]) < cv
  ar_test_0 <- mod$AR$p.value < 0.05
  
  mod2 <- ivmodel(y, x, Z, intercept=FALSE, beta0=1)
  cof <- coef(mod2)
  ols_t_stat_1 <- abs(cof[1,4]) < cv
  tsls_t_stat_1 <- abs(cof[2,4]) < cv
  ar_test_1 <- mod2$AR$p.value < 0.05

  el_est <- gel(g1,cbind(Z,x,y),tet=0.0,type="EL")
  el_coef <- coef(el_est)[1]
  el_se <- sqrt(vcov.gel(el_est)[1,1])
  el_t_stat_0 <- abs((el_coef / el_se)) < cv
  el_t_stat_1 <- abs(((el_coef - 1) / el_se)) < cv
  
  return(c(ols_est, tsls_est, el_coef, ols_t_stat_0, ols_t_stat_1, tsls_t_stat_0, tsls_t_stat_1, el_t_stat_0, el_t_stat_1, ar_test_0, ar_test_1))
}

set.seed(1)
n <- 200
l_list <- c(2, 3, 10, 20, 50, 150, 175, 198)
num_datasets <- 75
tsls_bias <- c()
ols_bias <- c()
el_bias <- c()
ols_t_acceptance_tot_0 <- c()
ols_t_acceptance_tot_1 <- c()
tsls_t_acceptance_tot_0 <- c()
tsls_t_acceptance_tot_1 <- c()
ar_acceptance_tot_0 <- c()
ar_acceptance_tot_1 <- c()
el_acceptance_tot_0 <- c()
el_acceptance_tot_1 <- c()
for(l in l_list) {
  ols_est_list <- c()
  tsls_est_list <- c()
  el_est_list <- c()
  ols_t_acceptance_0 <- c()
  ols_t_acceptance_1 <- c()
  tsls_t_acceptance_0 <- c()
  tsls_t_acceptance_1 <- c()
  ar_acceptance_0 <- c()
  ar_acceptance_1 <- c()
  el_acceptance_0 <- c()
  el_acceptance_1 <- c()
  res <- lapply(seq(1, num_datasets), compute_vals, n=n, l=l, cv=1.96)
  for(i in seq(1, num_datasets)) {
    ols_est_list <- c(ols_est_list, res[[i]][1])
    tsls_est_list <- c(tsls_est_list, res[[i]][2])
    el_est_list <- c(el_est_list, res[[i]][3])
    ols_t_acceptance_0 <- c(ols_t_acceptance_0, res[[i]][4])
    ols_t_acceptance_1 <- c(ols_t_acceptance_1, res[[i]][5])
    tsls_t_acceptance_0 <- c(tsls_t_acceptance_0, res[[i]][6])
    tsls_t_acceptance_1 <- c(tsls_t_acceptance_1, res[[i]][7])
    el_acceptance_0 <- c(el_acceptance_0, res[[i]][8])
    el_acceptance_1 <- c(el_acceptance_1, res[[i]][9])
    ar_acceptance_0 <- c(ar_acceptance_0, res[[i]][10])
    ar_acceptance_1 <- c(ar_acceptance_1, res[[i]][11])
  }
  tsls_bias <- c(tsls_bias, mean(tsls_est_list))
  ols_bias <- c(ols_bias, mean(ols_est_list))
  el_bias <- c(el_bias, mean(el_est_list))
  ols_t_acceptance_tot_0 <- c(ols_t_acceptance_tot_0, mean(ols_t_acceptance_0))
  ols_t_acceptance_tot_1 <- c(ols_t_acceptance_tot_1, mean(ols_t_acceptance_1))
  tsls_t_acceptance_tot_0 <- c(tsls_t_acceptance_tot_0, mean(tsls_t_acceptance_0))
  tsls_t_acceptance_tot_1 <- c(tsls_t_acceptance_tot_1, mean(tsls_t_acceptance_1))
  el_acceptance_tot_0 <- c(el_acceptance_tot_0, mean(el_acceptance_0))
  el_acceptance_tot_1 <- c(el_acceptance_tot_1, mean(el_acceptance_1))
  ar_acceptance_tot_0 <- c(ar_acceptance_tot_0, mean(ar_acceptance_0))
  ar_acceptance_tot_1 <- c(ar_acceptance_tot_1, mean(ar_acceptance_1))
}
plot(x=l_list,y=tsls_bias, type='b', lty=1,lwd=1, xlab="Number of instruments", ylab="Bias", ylim=range( c(min(tsls_bias) - .2, max(ols_bias) + .5)))
lines(l_list,el_bias,type='b', lty=2, lwd=2)
lines(l_list,ols_bias, type='b', lty=3,lwd=3)
legend("bottomright",legend=c('2SLS bias', 'EL bias', 'OLS bias'),lty=c(1,2,3),bg="white",lwd=2)
