library(foreign)
library(parallel)
library(sandwich)
library(lmtest)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
base <- 4
clusterEvalQ(cl, library(sandwich))
clusterEvalQ(cl, library(lmtest))
clusterExport(cl, "base")

wild_sample <- function(u_val) {
  const <- (sqrt(5) - 1) / (2*sqrt(5))
  s <- unlist(lapply(runif(length(u_val)), function(u, u_val) { 
    if (u < const) { 
      return((1 + sqrt(5)) / 2)
    } else { 
      return ((1 - sqrt(5)) / 2)
    }
  }))
  return (s * u_val)
}

clusterExport(cl, c('wild_sample'))


nonparametric_bootstrap <- function(sample, B, n, betas, cl) {
  m <- parLapply(cl, 1:B, function(i, data, n, betas) {
    s <- data[sample(nrow(data),n,replace=TRUE),]
    new_model <- lm(lwage1 ~ educ + marr + nonwhite + covered + exper, data=s)
    errs <- sqrt(diag(vcovHC(new_model)))
    t <- (coef(new_model) - betas) / errs
    return(t)
  }, data = sample, n = n, betas = betas)
  m <- matrix(unlist(m), ncol=length(betas), nrow=B)
  cv <- apply(m, 2, function(col) {
    return(c(quantile(col, .05), quantile(col, .95), quantile(abs(col), .9)))
  })
  return(cv)
}

wild_bootstrap <- function(sample, B, n, sample_beta, sample_residuals, cl) {
  ones <- rep(1,n)
  m <- parLapply(cl, 1:B, function(i, data, n, sample_betas, u, ones) {
    u_star <- wild_sample(u)
    x_star <- cbind(ones, data['educ'], data['marr'], data['nonwhite'], data['covered'], data['exper'])
    data$ystar <- as.matrix(x_star) %*% sample_betas + u_star
    new_model <- lm(ystar ~ educ + marr + nonwhite + covered + exper, data=data)
    errs <- sqrt(diag(vcovHC(new_model)))
    t <- (coef(new_model) - sample_betas) / errs
    return(t)
  }, data = sample, n = n, sample_betas=sample_beta, u=sample_residuals, ones=ones)
  m <- matrix(unlist(m), ncol=length(sample_betas), nrow=B)
  cv <- apply(m, 2, function(col) {
    return(c(quantile(col, .05), quantile(col, .95), quantile(abs(col), .9)))
  })
  return(cv)
}

reject_equivalent <- function(i, t_stat, cvs) {
  return(t_stat[i] < cvs[2,i] && t_stat[i] > cvs[1,i])
}

reject_symmetric <- function(i, t_stat, cvs) {
  return(abs(t_stat[i]) < cvs[3,i])
}

data <- read.dta("/Users/garidor/Desktop/first-year/spring/metrics/ps11/USMenWages.dta")
model <- lm(lwage1 ~ educ + marr + nonwhite + covered + exper, data=data)
true_beta <- coef(model)
n <- 500
B <- 1000
num_iterations <- 500
non_parametric_equivalent_reject_sum <- rep(0,length(true_beta))
non_parametric_symmetric_reject_sum <- rep(0,length(true_beta))
wild_equivalent_reject_sum <- rep(0,length(true_beta))
wild_symmetric_reject_sum <- rep(0,length(true_beta))
clt_reject_sum <- rep(0,length(true_beta))
for (i in seq(1, num_iterations)) {
  sample_data <- data[sample(nrow(data), n),]
  sample_model <- lm(lwage1 ~ educ + marr + nonwhite + covered + exper, data=sample_data)
  u <- sample_model$residuals
  sample_betas <- coef(sample_model)
  std_errors <- coeftest(sample_model)
  t_stat <- (sample_betas - true_beta) / std_errors[,2]
  clt_cv <- qnorm(.95)
  non_parametric_cvs <- nonparametric_bootstrap(sample_data,B,n,as.matrix(sample_betas), cl)
  wild_cvs <- wild_bootstrap(sample_data,B,n,as.matrix(sample_betas),u, cl)
  non_parametric_equivalent_reject <- unlist(lapply(1:6, reject_equivalent, t_stat = t_stat, cvs = non_parametric_cvs))
  non_parametric_equivalent_reject_sum <- non_parametric_equivalent_reject_sum + non_parametric_equivalent_reject
  non_parametric_symmetric_reject <- unlist(lapply(1:6, reject_symmetric, t_stat = t_stat, cvs = non_parametric_cvs))
  non_parametric_symmetric_reject_sum <- non_parametric_symmetric_reject_sum + non_parametric_symmetric_reject
  wild_equivalent_reject <- unlist(lapply(1:6, reject_equivalent, t_stat = t_stat, cvs = wild_cvs))
  wild_equivalent_reject_sum <- wild_equivalent_reject_sum + wild_equivalent_reject
  wild_symmetric_reject <- unlist(lapply(1:6, reject_symmetric, t_stat = t_stat, cvs = wild_cvs))
  wild_symmetric_reject_sum <- wild_symmetric_reject_sum + wild_symmetric_reject
  clt_reject <- abs(t_stat) < clt_cv
  clt_reject_sum <- clt_reject_sum + clt_reject
}

cov_t <- function(est, cov, num_iterations) {
  return((est-cov)/sqrt(cov*(1-cov)/num_iterations))
}
cov <- .9
clt_cov_est <- clt_reject_sum / num_iterations
clt_t <- cov_t(clt_cov_est, cov, num_iterations)
wild_cov_sym_est <- wild_symmetric_reject_sum / num_iterations
wild_sym_t <- cov_t(wild_cov_sym_est, cov, num_iterations)
wild_cov_eq_est <- wild_equivalent_reject_sum / num_iterations
wild_eq_t <- cov_t(wild_cov_eq_est, cov, num_iterations)

np_cov_sym_est <- non_parametric_symmetric_reject_sum / num_iterations
np_sym_t <- cov_t(np_cov_sym_est, cov, num_iterations)
np_cov_eq_est <- non_parametric_equivalent_reject_sum / num_iterations
np_eq_t <- cov_t(np_cov_eq_est, cov, num_iterations)

results <- matrix(c(clt_cov_est, wild_cov_eq_est, wild_cov_sym_est, np_cov_eq_est, np_cov_sym_est), ncol=5, nrow=6)
colnames(results) <- c("CLT", "Wild Eq", "Wild Sym", "NP Eq", "NP Sym")
rownames(results) <- c("(Intercept)", "educ", "marr", "nonwhite", "covered", "exper")

t_results <- abs(results) < 1.96
colnames(t_results) <- c("CLT", "Wild Eq", "Wild Sym", "NP Eq", "NP Sym")
rownames(t_results) <- c("(Intercept)", "educ", "marr", "nonwhite", "covered", "exper")

stopCluster(cl)