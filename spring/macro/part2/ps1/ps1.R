library(vars)

gdp <- read.csv("/Users/garidor/Desktop/first-year/spring/macro/part2/ps1/GDP.csv")[2]
gdp_deflator <- read.csv("/Users/garidor/Desktop/first-year/spring/macro/part2/ps1/GDPDEF.csv")[2]
fed_funds <- read.csv("/Users/garidor/Desktop/first-year/spring/macro/part2/ps1/FEDFUNDS.csv")
pmi <- read.csv("/Users/garidor/Desktop/first-year/spring/macro/part2/ps1/ISM-MAN_PMI.csv")
pmi <- pmi[order(nrow(pmi):1),]

inflation <- c()
growth <- c()
for(i in seq(2, nrow(gdp))) {
  inflation <- c(inflation, 400*log(gdp_deflator[i,]/gdp_deflator[i-1,]))
  growth <- c(growth, 400 * log(gdp[i,] / gdp[i-1,]))
}

# oh boy
cur_year <- as.numeric(unlist(strsplit(toString(fed_funds[1, 1]), '-'))[1])
cur_rate_sum <- 0
cur_pmi_sum <- 0
num_periods <- 0
quarterly_rate <- c()
quarterly_pmi <- c()
for (i in seq(1, nrow(fed_funds[2]))) {
  date <- unlist(strsplit(toString(fed_funds[i, 1]), '-'))
  month <- as.numeric(date[2])
  year <- as.numeric(date[1])
  if (month == 4 || month == 7 || month == 10 || year > cur_year) {
    quarterly_rate <- c(quarterly_rate, cur_rate_sum / num_periods)
    quarterly_pmi <- c(quarterly_pmi, cur_pmi_sum / num_periods)
    cur_year <- year
    num_periods <- 0
    cur_rate_sum <- 0
    cur_pmi_sum <- 0
  }
  num_periods <- num_periods + 1
  cur_pmi_sum <- cur_pmi_sum + pmi[i, 2]
  cur_rate_sum <- cur_rate_sum + fed_funds[i, 2]
}

# question B
b_dat <- cbind(inflation, growth, quarterly_rate)
b_model <- VAR(b_dat, p=4, type="none")
b_mat <- cbind(c(NA,NA,NA), c(0,NA,NA), c(0,0,NA))
b_model <- SVAR(b_model, Bmat=b_mat)


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
dat <- cbind(inflation, growth, quarterly_pmi, quarterly_rate)
c_model <- VAR(dat, p=4, type="none")
b_mat <- cbind(c(NA,NA,NA,NA), c(0,NA,NA,NA), c(0,0,NA,NA), c(0,0,0,NA))
c_model <- SVAR(c_model, Bmat=b_mat)

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
cutoff_point <- 79
pre_dat <- b_dat[1:cutoff_point,]
post_dat <- b_dat[(cutoff_point+1):nrow(b_dat),]

#pre 1979
pre_model <- VAR(b_dat, p=4, type="none")
b_mat <- cbind(c(NA,NA,NA), c(0,NA,NA), c(0,0,NA))
pre_model <- SVAR(pre_model, Bmat=b_mat)

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
post_model <- VAR(post_dat, p=4, type="none")
b_mat <- cbind(c(NA,NA,NA), c(0,NA,NA), c(0,0,NA))
post_model <- SVAR(post_model, Bmat=b_mat)

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

  