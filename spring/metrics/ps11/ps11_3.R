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
    new_model <- coeftest(new_model)
    t <- (new_model[,1] - betas) / new_model[,2]
    return(t)
  }, data = sample, n = n, betas = betas)
  m <- matrix(unlist(m), ncol=length(betas), nrow=B)
  cv <- apply(m, 2, function(col) {
    return(c(quantile(col, .05), quantile(col, .95), quantile(abs(col), .9)))
  })
  return(cv)
}

wild_bootstrap <- function(sample, B, n, true_betas, sample_beta, sample_residuals, cl) {
  ones <- rep(1,n)
  m <- parLapply(cl, 1:B, function(i, data, n, true_betas, sample_betas, u, ones) {
    u_star <- wild_sample(u)
    x_star <- cbind(ones, data['educ'], data['marr'], data['nonwhite'], data['covered'], data['exper'])
    data$ystar <- as.matrix(x_star) %*% sample_betas + u_star
    new_model <- lm(ystar ~ educ + marr + nonwhite + covered + exper, data=data)
    model_summary <- coeftest(new_model)
    t <- (model_summary[,1] - true_betas) / model_summary[,2]
    return(t)
  }, data = sample, n = n, true_betas = true_betas, sample_betas=sample_beta, u=sample_residuals, ones=ones)
  m <- matrix(unlist(m), ncol=length(true_betas), nrow=B)
  cv <- apply(m, 2, function(col) {
    return(c(quantile(col, .05), quantile(col, .95), quantile(abs(col), .9)))
  })
  return(cv)
}

reject_equivalent <- function(i, t_stat, cvs) {
  return(t_stat[i] > cvs[2,i] || t_stat[i] < cvs[1,i])
}

reject_symmetric <- function(i, t_stat, cvs) {
  return(abs(t_stat[i]) > cvs[3,i])
}

data <- read.dta("/Users/garidor/Desktop/first-year/spring/metrics/ps11/USMenWages.dta")
model <- lm(lwage1 ~ educ + marr + nonwhite + covered + exper, data=data)
true_beta <- coef(model)
n <- 500
B <- 10
num_iterations <- 5
non_parametric_equivalent_reject_sum <- rep(0,length(true_beta))
non_parametric_symmetric_reject_sum <- rep(0,length(true_beta))
wild_equivalent_reject_sum <- rep(0,length(true_beta))
wild_symmetric_reject_sum <- rep(0,length(true_beta))
clt_reject_sum <- rep(0,length(true_beta))
for (i in seq(1, 500)) {
  sample_data <- data[sample(nrow(data), 500),]
  sample_model <- lm(lwage1 ~ educ + marr + nonwhite + covered + exper, data=sample_data)
  u <- sample_model$residuals
  sample_betas <- coef(sample_model)
  std_errors <- coeftest(sample_model)
  t_stat <- (sample_betas - true_beta) / std_errors[,2]
  clt_cv <- qnorm(.95)
  non_parametric_cvs <- nonparametric_bootstrap(sample_data,B,n,true_beta, cl)
  wild_cvs <- wild_bootstrap(sample_data,B,n,true_beta,as.matrix(sample_betas),u, cl)
  non_parametric_equivalent_reject <- unlist(lapply(1:6, reject_equivalent, t_stat = t_stat, cvs = non_parametric_cvs))
  non_parametric_equivalent_reject_sum <- non_parametric_equivalent_reject_sum + non_parametric_equivalent_reject
  non_parametric_symmetric_reject <- unlist(lapply(1:6, reject_symmetric, t_stat = t_stat, cvs = non_parametric_cvs))
  non_parametric_symmetric_reject_sum <- non_parametric_symmetric_reject_sum + non_parametric_symmetric_reject
  wild_equivalent_reject <- unlist(lapply(1:6, reject_equivalent, t_stat = t_stat, cvs = wild_cvs))
  wild_equivalent_reject_sum <- wild_equivalent_reject_sum + wild_equivalent_reject
  wild_symmetric_reject <- unlist(lapply(1:6, reject_symmetric, t_stat = t_stat, cvs = wild_cvs))
  wild_symmetric_reject_sum <- wild_symmetric_reject_sum + wild_symmetric_reject
  clt_reject <- t_stat > clt_cv
  clt_reject_sum <- clt_reject_sum + clt_reject
}

stopCluster(cl)