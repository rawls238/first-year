library(foreign)
library(parallel)
library(sandwich)
library(lmtest)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
base <- 4
clusterExport(cl, "base")

nonparametric_bootstrap <- function(sample, B, n, betas) {
  m <- lapply(1:B, function(i, data, n, betas) {
    s <- data[sample(nrow(data),n,replace=TRUE),]
    new_model <- regress(s)
    t <- (new_model[,1] - betas) / new_model[,2]
    return(t)
  }, data = sample, n = n, betas = betas)
  m <- matrix(unlist(m), ncol=length(betas), nrow=B)
  cv <- apply(m, 2, function(col) {
    return(c(quantile(col, .05), quantile(col, .95), quantile(abs(col), .9)))
  })
  return(cv)
}

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

wild_bootstrap <- function(sample, B, n, true_betas, sample_beta, sample_residuals) {
  ones <- rep(1,n)
  m <- lapply(1:B, function(i, data, n, true_betas, sample_betas, u, ones) {
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
 
regress <- function(data) {
  model <- lm(lwage1 ~ educ + marr + nonwhite + covered + exper, data=data)
  vars <- coeftest(model)
  return(vars)
}

data <- read.dta("/Users/garidor/Desktop/first-year/spring/metrics/ps11/USMenWages.dta")
model <- lm(lwage1 ~ educ + marr + nonwhite + covered + exper, data=data)
true_beta <- coef(model)
n <- 500
B <- 10
num_iterations <- 5
for (i in seq(1, 500)) {
  sample_data <- data[sample(nrow(data), 500),]
  sample_model <- lm(lwage1 ~ educ + marr + nonwhite + covered + exper, data=sample_data)
  u <- sample_model$residuals
  sample_betas <- coef(sample_model)
  std_errors <- coeftest(sample_model)
  t_stat <- (sample_betas - true_beta) / std_errors[,2]
  clt_cv <- qnorm(.95)
  non_parametric_cvs <- nonparametric_bootstrap(sample_data,B,n,true_beta)
  wild_cvs <- wild_bootstrap(sample_data,B,n,true_beta,as.matrix(sample_betas),u)
  
}

stopCluster(cl)