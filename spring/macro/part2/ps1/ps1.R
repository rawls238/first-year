library(vars)

data <- read.csv("/Users/garidor/Desktop/first-year/spring/macro/part2/ps1/ps1_data.csv")

gdp <- data[,3][3:194]
gdp <- as.numeric(as.character(gdp))
gdp_deflator <- data[,2][3:194]
gdp_deflator <- as.numeric(as.character(gdp_deflator))
fed_funds <- data[,6]
fed_funds <- fed_funds[3:length(fed_funds)]
fed_funds <- as.numeric(as.character(fed_funds))
pmi <- data[,7]
pmi <- pmi[3:length(pmi)]
pmi <- as.numeric(as.character(pmi))


inflation <- c()
growth <- c()
for(i in seq(1, length(gdp))) {
  inflation <- c(inflation, 400*log(gdp_deflator[i]/gdp_deflator[i-1]))
  growth <- c(growth, 400 * log(gdp[i] / gdp[i-1]))
}

quarterly_rate <- c()
quarterly_pmi <- c()
for (i in seq(2, 192)) {
  b <- (i-1) * 3
  quarterly_pmi <- c(quarterly_pmi, (pmi[b+1] + pmi[b+2] + pmi[b+3]) / 3)
  quarterly_rate <- c(quarterly_rate, (fed_funds[b+1] + fed_funds[b+2] + fed_funds[b+3]) / 3)
}

# question B
b_dat <- cbind(inflation, growth, quarterly_rate)
b_model <- VAR(b_dat, p=4)
b_mat <- cbind(c(NA,NA,NA), c(0,NA,NA), c(0,0,NA))
b_model <- SVAR(b_model, Bmat=b_mat)
print(b_model)

b_irfs <- irf(b_model, n.ahead=20)
plot(b_irfs)


fevd_b <- fevd(b_model, n.ahead=20)
inf <- fevd_b$inflation * 100
gr <- fevd_b$growth * 100
ffr <- fevd_b$quarterly_rate * 100
print("Numbers below are Inflation, Growth, FFR (Percentage Points)")
print("Inflation FEVD")
cat("1 period ahead", inf[1,], "\n")
cat("8 periods ahead", inf[8,], "\n")
cat("20 periods ahead", inf[20,], "\n")
print("Growth Rate FEVD")
cat("1 period ahead", gr[1,], "\n")
cat("8 periods ahead", gr[8,], "\n")
cat("20 periods ahead", gr[20,], "\n")
print("FFR FEVD")
cat("1 period ahead", ffr[1,], "\n")
cat("8 periods ahead", ffr[8,], "\n")
cat("20 periods ahead", ffr[20,], "\n")
  
  
# question C
c_dat <- cbind(inflation, growth, quarterly_pmi, quarterly_rate)
c_model <- VAR(c_dat, p=4)
b_mat <- cbind(c(NA,NA,NA,NA), c(0,NA,NA,NA), c(0,0,NA,NA), c(0,0,0,NA))
c_model <- SVAR(c_model, Bmat=b_mat)
print(c_model)

c_irfs <- irf(c_model, n.ahead=20)
plot(c_irfs)

fevd_c <- fevd(c_model, n.ahead=20)
inf <- fevd_c$inflation * 100
gr <- fevd_c$growth * 100
ffr <- fevd_c$quarterly_rate * 100
fpmi <- fevd_c$quarterly_pmi * 100
print("Numbers below are Inflation, Growth, PMI, FFR (Percentage Points) including PMI")
print("Inflation FEVD")
cat("1 period ahead", inf[1,], "\n")
cat("8 periods ahead", inf[8,], "\n")
cat("20 periods ahead", inf[20,], "\n")
print("Growth Rate FEVD")
cat("1 period ahead", gr[1,], "\n")
cat("8 periods ahead", gr[8,], "\n")
cat("20 periods ahead", gr[20,], "\n")
print("PMI FEVD")
cat("1 period ahead", fpmi[1,], "\n")
cat("8 periods ahead", fpmi[8,], "\n")
cat("20 periods ahead", fpmi[20,], "\n")
print("FFR FEVD")
cat("1 period ahead", ffr[1,], "\n")
cat("8 periods ahead", ffr[8,], "\n")
cat("20 periods ahead", ffr[20,], "\n")


#question D
cutoff_point <- 78
pre_dat <- b_dat[1:cutoff_point,]
post_dat <- b_dat[(cutoff_point+1):nrow(b_dat),]

#pre 1979
pre_model <- VAR(pre_dat, p=4)
b_mat <- cbind(c(NA,NA,NA), c(0,NA,NA), c(0,0,NA))
pre_model <- SVAR(pre_model, Bmat=b_mat)
print(pre_model)

pre_irfs <- irf(pre_model, n.ahead=20)
plot(pre_irfs)

fevd_pre <- fevd(pre_model, n.ahead=20)
inf <- fevd_pre$inflation * 100
gr <- fevd_pre$growth * 100
ffr <- fevd_pre$quarterly_rate * 100
print("Numbers below are Inflation, Growth, FFR (Percentage Points) for pre-1979")
print("Inflation FEVD")
cat("1 period ahead", inf[1,], "\n")
cat("8 periods ahead", inf[8,], "\n")
cat("20 periods ahead", inf[20,], "\n")
print("Growth Rate FEVD")
cat("1 period ahead", gr[1,], "\n")
cat("8 periods ahead", gr[8,], "\n")
cat("20 periods ahead", gr[20,], "\n")
print("FFR FEVD")
cat("1 period ahead", ffr[1,], "\n")
cat("8 periods ahead", ffr[8,], "\n")
cat("20 periods ahead", ffr[20,], "\n")

#post 1979
post_model <- VAR(post_dat, p=4)
b_mat <- cbind(c(NA,NA,NA), c(0,NA,NA), c(0,0,NA))
post_model <- SVAR(post_model, Bmat=b_mat)
print(post_model)

post_irfs <- irf(post_model, n.ahead=20)
plot(post_irfs)

fevd_post <- fevd(post_model, n.ahead=20)
inf <- fevd_post$inflation * 100
gr <- fevd_post$growth * 100
ffr <- fevd_post$quarterly_rate * 100
print("Numbers below are Inflation, Growth, FFR (Percentage Points) for post-1979")
print("Inflation FEVD")
cat("1 period ahead", inf[1,], "\n")
cat("8 periods ahead", inf[8,], "\n")
cat("20 periods ahead", inf[20,], "\n")
print("Growth Rate FEVD")
cat("1 period ahead", gr[1,], "\n")
cat("8 periods ahead", gr[8,], "\n")
cat("20 periods ahead", gr[20,], "\n")
print("FFR FEVD")
cat("1 period ahead", ffr[1,], "\n")
cat("8 periods ahead", ffr[8,], "\n")
cat("20 periods ahead", ffr[20,], "\n")

  