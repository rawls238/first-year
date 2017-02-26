library("quantreg")

data("CPS1985", package = "AER")
cps <- CPS1985
model <- rq(log(wage) ~ education + experience, data=cps,tau=seq(0.01, 0.99, by=0.01))
plot(summary(model))