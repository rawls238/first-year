library("foreign")
library("sandwich")
library("lmtest")
library("MASS")
data <- read.dta("/Users/garidor/Desktop/first-year/spring/metrics/ps5/lakisha_aer.dta")
data1 <- data
#b
ratio <- xtabs(~ race + call, data)
ratio_table <- prop.table(ratio)
ch <- chisq.test(ratio)
print(ratio_table)
cat("Reject the null hypothesis that they are equal", ch$p.value < 0.05)

#d
d_probit <- glm(call ~ race + sex + education + yearsexp + military + col + h, data=data, family=binomial(link="probit"))
d_logit <- glm(call ~  race + sex + education + yearsexp + military + col + h, data=data, family=binomial(link="logit"))
d_ols <- lm(call ~ race + sex + education + yearsexp + military + col + h, data=data)
d_ols_se <- coeftest(d_ols, vcov. = vcovHC)

#e

e_probit <- glm(call ~  eoe + race*eoe + sex*eoe + race + sex + education + yearsexp + military + col + h, data=data, family=binomial(link="probit"))
e_logit <- glm(call ~ eoe + race*eoe + sex*eoe + race + sex + education + yearsexp + military + col + h, data=data, family=binomial(link="logit"))
e_ols <- lm(call ~ eoe + race*eoe + sex*eoe + race + sex + education + yearsexp + military + col + h, data=data)
e_ols_se <- coeftest(d_ols, vcov. = vcovHC)

#f
probit_coeff <- coef(e_probit)
logit_coeff <- coef(e_logit)
mean_educ <- mean(data$education)
mean_exp <- mean(data$yearsexp)
predicted_1_w_p <-  predict(e_probit, data.frame(sex='m',race='w',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=1,military=0), type="response")
predicted_1_b_p <- predict(e_probit, data.frame(sex='m',race='b',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=1,military=0), type="response")
predicted_1_w_l <- predict(e_logit, data.frame(sex='m',race='w',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=1,military=0), type="response")
predicted_1_b_l <- predict(e_logit, data.frame(sex='m',race='b',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=1,military=0), type="response")
diff_p_1 <- predicted_1_w_p - predicted_1_b_p
diff_l_1 <- predicted_1_w_l - predicted_1_b_l 

predicted_2_w_p <- predict(e_probit, data.frame(sex='f',race='w',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=1,military=0), type="response")
predicted_2_b_p <- predict(e_probit, data.frame(sex='f',race='b',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=1,military=0), type="response")
predicted_2_w_l <-  predict(e_logit, data.frame(sex='f',race='w',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=1,military=0), type="response")
predicted_2_b_l <- predict(e_logit, data.frame(sex='f',race='b',education=mean_educ, yearsexp=mean_exp,h=1,col=1,eoe=1,military=0), type="response")
diff_p_2 <- predicted_2_w_p - predicted_2_b_p
diff_l_2 <- predicted_2_w_l - predicted_2_b_l 

data$dummy_educ0 <- as.numeric(data$education == 0)
data$dummy_educ1 <- as.numeric(data$education == 1)
data$dummy_educ2 <- as.numeric(data$education == 2)
data$dummy_educ3 <- as.numeric(data$education == 3)
data$dummy_educ4 <- as.numeric(data$education == 4)
f_probit <- glm(call ~  eoe + race*eoe + sex*eoe + race + sex + dummy_educ0 + dummy_educ1 + dummy_educ2 + dummy_educ3 + dummy_educ4 + yearsexp + military + col + h, data=data, family=binomial(link="probit"))
f_logit <- glm(call ~ eoe + race*eoe + sex*eoe + race + sex + dummy_educ0 + dummy_educ1 + dummy_educ2 + dummy_educ3 + dummy_educ4  + yearsexp + military + col + h, data=data, family=binomial(link="logit"))


data$lo <- f_logit$coefficients['education']*dlogis(predict(f_logit,data))
data$pr <- f_probit$coefficients['education']*dnorm(predict(f_probit,data))

par_effect_logit_w = c()
par_effect_probit_w = c()
par_effect_logit_b = c()
par_effect_probit_b = c()
par_effect_probit_w = c(par_effect_probit_w, mean(data[data$race=='w' & data$sex=='f',]$pr,na.rm=TRUE))
par_effect_logit_w = c(par_effect_logit_w, mean(data[data$race=='w' & data$sex=='f',]$lo,na.rm=TRUE))
par_effect_probit_b = c(par_effect_probit_b, mean(data[data$race=='b' & data$sex=='f',]$pr,na.rm=TRUE))
par_effect_logit_b = c(par_effect_logit_b, mean(data[data$race=='b' & data$sex=='f',]$lo,na.rm=TRUE))


#q2
data <- read.dta("/Users/garidor/Desktop/first-year/spring/metrics/ps5/mus10data.dta")
data <- subset(data,year01==1)
model_poisson <- glm(docvis ~ private + chronic + female + income, data=data,family=poisson(link = "log"))
hist(data['docvis'][!is.na(data['docvis'])], xlim=c(0, 25), breaks=150, main="Data")
hist(predict(model_poisson), main="Poisson")
model_nb <- glm.nb(docvis ~ private + chronic + female + income, data=data)
hist(predict(model_nb), main="Negative Binomial")

r = seq(round(min(data['income'])), round(max(data['income'])))
vals_chronic = 1 - exp(-1 * (coef(model_poisson)['income'] * r + coef(model_poisson)['chronic'] + coef(model_poisson)['private'] + coef(model_poisson)['female']))
vals_no_chronic = 1 - exp(-1 * (coef(model_poisson)['income'] * r + coef(model_poisson)['private'] + coef(model_poisson)['female']))
plot(vals_chronic, main="Chronic")
plot(vals_no_chronic, main="Not Chronic")

