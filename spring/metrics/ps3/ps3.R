# First, we load the data set, which is taken from the 1998 CPS
data("CPS1985", package = "AER")
cps <- CPS1985

# To illustrate the usage of NLS, we run a nonlinear regression of wages on the
# exponential of a linear index of education and experience

mod <-nls(log(wage) ~ a1 + a2*education + a3*experience + a4*(experience / (1 + exp(-(experience - a5) / a6))), data = cps,
    start = list(a1 = 2.0, a2 = 32.0, a3=34.0, a4=6.0, a5=7.0,a6=8.0))

ceduc <- cps["education"]
educ_mean <- mean(ceduc[!is.na(ceduc)])
c <- coef(mod)
r <- min(cps["experience"]):max(cps["experience"])
predicted <- c["a1"] + c["a2"] * educ_mean + c["a3"] * r + c["a4"] * (r / (1 + exp(-(r - c["a5"]) / c["a6"])))
plot(predicted)