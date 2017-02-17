library("foreign")
library("sandwich")
library("lmtest")
library("MASS")
data <- read.dta("/Users/garidor/Desktop/first-year/spring/metrics/ps5/lakisha_aer.dta")

#b
ratio <- xtabs(~ race + call, data)
ratio_table <- prop.table(ratio)
c <- chisq.test(ratio)
cat("Reject the null hypothesis that they are equal", c$p.value < 0.05)

#d

d_probit <- glm(call ~ 0 + race + sex + education + yearsexp + military + col + h, data=data, family=binomial(link="probit"))
d_logit <- glm(call ~ 0 + race + sex + education + yearsexp + military + col + h, data=data, family=binomial(link="logit"))
d_ols <- lm(call ~ 0 + race + sex + education + yearsexp + military + col + h, data=data)
d_ols_se <- coeftest(d_ols, vcov. = vcovHC)

#e

e_probit <- glm(call ~  eoe + race*eoe + sex*eoe + race + sex + education + yearsexp + military + col + h, data=data, family=binomial(link="probit"))
e_logit <- glm(call ~ eoe + race*eoe + sex*eoe + race + sex + education + yearsexp + military + col + h, data=data, family=binomial(link="logit"))
e_ols <- lm(call ~ eoe + race*eoe + sex*eoe + race + sex + education + yearsexp + military + col + h, data=data)
e_ols_se <- coeftest(d_ols, vcov. = vcovHC)

#f
probit_coeff <- coef(e_probit)
logit_coeff <- coef(e_logit)
mean_educ <- mean(!is.na(data["education"]))
mean_exp <- mean(!is.na(data["yearsexp"]))
predicted_1_w_p <-  probit_coeff['racew'] + probit_coeff['sexm'] + probit_coeff['education'] * mean_educ + probit_coeff['yearsexp'] * mean_exp + probit_coeff['h'] + probit_coeff['col']
predicted_1_b_p <- probit_coeff['sexm'] + probit_coeff['education'] * mean_educ + probit_coeff["yearsexp"] * mean_exp + probit_coeff['h'] + probit_coeff['col']
predicted_1_w_l <- (logit_coeff['racew'] + logit_coeff['sexm'] + logit_coeff['education'] * mean_educ + logit_coeff['yearsexp'] * mean_exp + logit_coeff['h'] + logit_coeff['col'])
predicted_1_b_l <-  (logit_coeff['sexm'] + logit_coeff['education'] * mean_educ + logit_coeff["yearsexp"] * mean_exp + logit_coeff['h'] + logit_coeff['col'])

predicted_2_w_p <- probit_coeff['racew'] + probit_coeff['education'] * mean_educ + probit_coeff['yearsexp'] * mean_exp + probit_coeff['h'] + probit_coeff['col']
predicted_2_b_p <-probit_coeff['education'] * mean_educ + probit_coeff["yearsexp"] * mean_exp + probit_coeff['h'] + probit_coeff['col']
predicted_2_w_l <-  logit_coeff['racew'] + logit_coeff['education'] * mean_educ + logit_coeff['yearsexp'] * mean_exp + logit_coeff['h'] + logit_coeff['col']
predicted_2_b_l <-  logit_coeff['education'] * mean_educ + logit_coeff["yearsexp"] * mean_exp + logit_coeff['h'] + logit_coeff['col']

prev_probit_wf <- 0.0
prev_logit_wf <- 0.0
prev_probit_bf <- 0.0
prev_logit_bf <- 0.0
diff_pwf <- c()
diff_pbf <- c()
diff_lwf <- c()
diff_lbf <- c()
for(i in seq(1,4)) {
  a <- i * probit_coeff['education'] + probit_coeff['racew'] + probit_coeff['yearsexp'] * mean_exp + probit_coeff['h'] + probit_coeff['col']
  b <- predicted_3_bf_p <- i * probit_coeff['education'] + probit_coeff['yearsexp'] * mean_exp + probit_coeff['h'] + probit_coeff['col']
  d <- i * logit_coeff['education'] + logit_coeff['racew'] + logit_coeff['yearsexp'] * mean_exp + logit_coeff['h'] + logit_coeff['col']
  e <- logit_coeff['education'] * i + logit_coeff["yearsexp"] * mean_exp + logit_coeff['h'] + logit_coeff['col']
  if (i >= 2) {
    diff_pwf <- c(diff_pwf, a - prev_probit_wf)
    diff_pbf <- c(diff_pbf, b - prev_probit_bf)
    diff_lwf <- c(diff_lwf, d - prev_logit_wf)
    diff_lbf <- c(diff_lbf, e - prev_logit_bf)
  }
  prev_probit_wf <- a
  prev_logit_wf <- d
  prev_probit_bf <- b
  prev_logit_bf <- e
}


#q2
data <- read.dta("/Users/garidor/Desktop/first-year/spring/metrics/ps5/mus10data.dta")
data <- subset(data,year01==1)
model_poisson <- glm(docvis ~ private + chronic + female + income, data=data,family=poisson(link = "log"))
hist(data['docvis'][!is.na(data['docvis'])])
hist(predict(model_poisson))
model_nb <- glm.nb(docvis ~ private + chronic + female + income, data=data)
hist(predict(model_nb))

r = seq(round(min(data['income'])), round(max(data['income'])))
vals_chronic = 1 - exp(-1 * (coef(model_poisson)['income'] * r + coef(model_poisson)['chronic'] + coef(model_poisson)['private'] + coef(model_poisson)['female']))
vals_no_chronic = 1 - exp(-1 * (coef(model_poisson)['income'] * r + coef(model_poisson)['private'] + coef(model_poisson)['female']))
plot(vals_chronic)
plot(vals_no_chronic)

