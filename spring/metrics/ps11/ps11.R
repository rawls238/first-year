library(parallel)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
base <- 4
clusterExport(cl, "base")


bootstrap_1 <- function(row, B, n, q) {
  t <- mean(row)
  v <- var(row)
  m <- unlist(lapply(1:B, function(i, data, n) { 
    s <- sample(data, n, replace=TRUE)
    mean(s)
  }, data = row, n = n))
  bootstrap_est <- quantile(m, .95)
  clt_est <- mean(row) + q * sqrt(v)
  return(c(t, bootstrap_est, clt_est))
}


bootstrap_2 <- function(row, B, n, q) {
  v <- var(row)
  t <- sqrt(n)*(mean(row) - 1) / sqrt(v)
  m <- unlist(lapply(1:B, function(i, data, n) { 
    s <- sample(data, n, replace=TRUE)
    v2 <- var(s)
    return(sqrt(n)*(mean(s) - 1) / sqrt(v2))
  }, data = row, n = n))
  bootstrap_est <- quantile(m, .95)
  return(c(t, bootstrap_est, q))
}

set.seed(1)
n <- 20
B <- 999
num_samples <- 10000
q_95 <- qnorm(.95)
data <- matrix(rexp(n * num_samples), num_samples, n)
d <- parApply(cl, data, 1, bootstrap_1, B=B, n = n, q=q_95)

t_95_a <- quantile(d[1,], .95)
bootstrap_mean_a <- mean(d[2,])
bootstrap_var_a <- var(d[2,])
clt_mean_a <- mean(d[3,])
clt_var_a <- var(d[3,])

d2 <- parApply(cl, data, 1, bootstrap_2, B=B, n = n, q=q_95)
t_95_b <- quantile(d2[1,], .95)
bootstrap_mean_b <- mean(d2[2,])
bootstrap_var_b <- var(d2[2,])
clt_mean_b <- mean(d2[3,])
clt_var_b <- var(d2[3,])

a <- matrix(c(t_95_a, 0.00, clt_mean_a, clt_var_a, bootstrap_mean_a, bootstrap_var_a), ncol=3)
colnames(a) <- c("t", "CLT", "Bootstrap")
rownames(a) <- c("mean", "variance")

b <- matrix(c(t_95_b, 0.00, clt_mean_b, clt_var_b, bootstrap_mean_b, bootstrap_var_b), ncol=3)
colnames(b) <- c("t", "CLT", "Bootstrap")
rownames(b) <- c("mean", "variance")


stopCluster(cl)