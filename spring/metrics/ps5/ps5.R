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
predicted_1_w_p <-  predict(e_probit, data.frame(sex='m',race='w',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=0,military=0))
predicted_1_b_p <- predict(e_probit, data.frame(sex='m',race='b',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=0,military=0))
predicted_1_w_l <- predict(e_logit, data.frame(sex='m',race='w',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=0,military=0))
predicted_1_b_l <- predict(e_logit, data.frame(sex='m',race='b',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=0,military=0))
diff_p_1 <- predicted_1_w_p - predicted_1_b_p
diff_l_1 <- predicted_1_w_l - predicted_1_b_l 

predicted_2_w_p <- predict(e_probit, data.frame(sex='f',race='w',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=0,military=0))
predicted_2_b_p <- predict(e_probit, data.frame(sex='f',race='b',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=0,military=0))
predicted_2_w_l <-  predict(e_logit, data.frame(sex='f',race='w',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=0,military=0))
predicted_2_b_l <- predict(e_logit, data.frame(sex='f',race='b',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=0,military=0))
diff_p_2 <- predicted_2_w_p - predicted_2_b_p
diff_l_2 <- predicted_2_w_l - predicted_2_b_l 


prev_probit_wf <- 0.0
prev_logit_wf <- 0.0
prev_probit_bf <- 0.0
prev_logit_bf <- 0.0
diff_pwf <- c()
diff_pbf <- c()
diff_lwf <- c()
diff_lbf <- c()
for(i in seq(1,4)) {
  a <- predict(e_probit, data.frame(sex='f',race='w',education=i, yearsexp=mean_exp,h=1,col=1,eoe=0,military=0))
  b <- predict(e_probit, data.frame(sex='f',race='b',education=i, yearsexp=mean_exp,h=1,col=1,eoe=0,military=0))
  d <- predict(e_logit, data.frame(sex='f',race='w',education=i, yearsexp=mean_exp,h=1,col=1,eoe=0,military=0))
  e <- predict(e_logit, data.frame(sex='f',race='b',education=i, yearsexp=mean_exp,h=1,col=1,eoe=0,military=0))
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

