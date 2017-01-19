ar1.sim<-arima.sim(model=list(ar=c(.7)),n=200)
ts.plot(ar1.sim)
ar.acf<-acf(ar1.sim,20,type="correlation",plot=T)


arma.sim<-arima.sim(model=list(ar=c(.7,-.3),ma=c(.8,.4)),n=200)
ts.plot(arma.sim)
arma.acf<-acf(arma.sim,20,type="correlation",plot=T)