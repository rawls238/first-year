library(gmm)

processFile <- function(filepath) {
  con = file(filepath, "r")
  dat <- as.data.frame(matrix(0, ncol = 3))
  count <- 1
  while ( TRUE ) {
    line = readLines(con, n = 1)
    res <- as.numeric(trimws(strsplit(line, split="       ")[[1]]))
    if ( length(res) == 0 ) {
      break
    }
    dat[count,] <-res
    count <- count + 1
  }
  
  close(con)
  return(dat)
}

capmMomentCondition <- function(dat, X) {
  forward <- X[2:nrow(X),]
  X <- X[1:(nrow(X)-1), ]
  z <- cbind(X[,1], X[,3])
  return((dat[1] * (forward[,1]^(dat[2]-1)) * forward[,3] - 1) * z)
}

capmMomentConditionWithOneLag <- function(dat, X) {
  forward <- X[3:nrow(X),]
  laggedGrowth <- X[1:(nrow(X)-2),1] 
  X <- X[2:(nrow(X)-1),]
  z <- cbind(X[,1], X[,3], laggedGrowth)
  return((dat[1] * (forward[,1]^(dat[2]-1)) * forward[,3] - 1) * z)
}

capmMomentConditionWithTwoLags <- function(dat, X) {
  forward <- X[4:nrow(X),]
  laggedGrowth <- X[1:(nrow(X)-3),2]
  X <- X[3:(nrow(X)-1),]
  z <- cbind(X[,1], X[,3], laggedGrowth)
  return((dat[1] * (forward[,1]^(dat[2]-1)) * forward[,3] - 1) * z)
}

filename <- "/Users/garidor/Desktop/first-year/spring/metrics/ps8/capm-data.txt"
data <- processFile(filename)

cat("No lags", "\n")
identGmm <- gmm(capmMomentCondition,data,t0=c(.99,2.0),vcov="iid",wmatrix="ident",type="iterative")
optGmm <- gmm(capmMomentCondition,data,t0=c(.99,2.0),vcov="iid",wmatrix="optimal",type="twoStep")
elRes <- gel(capmMomentCondition,data,tet=c(.99,2.0),type="EL")
print(summary(identGmm))
print(summary(optGmm))
print(summary(elRes))


cat("One Lag", "\n")
identGmmOneLag <- gmm(capmMomentConditionWithOneLag,data,t0=c(.99,2.0),vcov="iid",wmatrix="ident",type="iterative")
optGmmOneLag <- gmm(capmMomentConditionWithOneLag,data,t0=c(.99,2.0),vcov="iid",wmatrix="optimal",type="twoStep")
elResOneLag <- gel(capmMomentConditionWithOneLag,data,tet=c(.99,20.0),type="EL")
print(summary(identGmmOneLag))
print(summary(optGmmOneLag))
print(summary(elResOneLag))

cat("Two lags", "\n")
identGmmTwoLag <- gmm(capmMomentConditionWithTwoLags,data,t0=c(.99,2.0),vcov="iid",wmatrix="ident",type="iterative")
optGmmTwoLag <- gmm(capmMomentConditionWithTwoLags,data,t0=c(.99,2.0),vcov="iid",wmatrix="optimal",type="twoStep")
elResTwoLag <- gel(capmMomentConditionWithTwoLags,data,tet=c(.99,2.0),type="EL")
print(summary(identGmmTwoLag))
print(summary(optGmmTwoLag))
print(summary(elResTwoLag))


capmMomentCondition <- function(dat, X) {
  forward <- X[2:nrow(X),]
  X <- X[1:(nrow(X)-1), ]
  z <- cbind(X[,1], X[,2])
  return((dat[1] * (forward[,1]^(dat[2]-1)) * forward[,3] - 1) * z)
}

capmMomentConditionWithOneLag <- function(dat, X) {
  forward <- X[3:nrow(X),]
  laggedGrowth <- X[1:(nrow(X)-2),1] 
  X <- X[2:(nrow(X)-1),]
  z <- cbind(X[,1], X[,2], laggedGrowth)
  return((dat[1] * (forward[,1]^(dat[2]-1)) * forward[,3] - 1) * z)
}

capmMomentConditionWithTwoLags <- function(dat, X) {
  forward <- X[4:nrow(X),]
  laggedGrowth <- X[1:(nrow(X)-3),2]
  X <- X[3:(nrow(X)-1),]
  z <- cbind(X[,1], X[,2], laggedGrowth)
  return((dat[1] * (forward[,1]^(dat[2]-1)) * forward[,3] - 1) * z)
}


cat("No lags", "\n")
identGmm <- gmm(capmMomentCondition,data,t0=c(.99,2.0),vcov="iid",wmatrix="ident",type="iterative")
optGmm <- gmm(capmMomentCondition,data,t0=c(.99,2.0),vcov="iid",wmatrix="optimal",type="twoStep")
elRes <- gel(capmMomentCondition,data,tet=c(.99,2.0),type="EL")
print(summary(identGmm))
print(summary(optGmm))
print(summary(elRes))


cat("One Lag", "\n")
identGmmOneLag <- gmm(capmMomentConditionWithOneLag,data,t0=c(.99,2.0),vcov="iid",wmatrix="ident",type="iterative")
optGmmOneLag <- gmm(capmMomentConditionWithOneLag,data,t0=c(.99,2.0),vcov="iid",wmatrix="optimal",type="twoStep")
elResOneLag <- gel(capmMomentConditionWithOneLag,data,tet=c(.99,20.0),type="EL")
print(summary(identGmmOneLag))
print(summary(optGmmOneLag))
print(summary(elResOneLag))

cat("Two lags", "\n")
identGmmTwoLag <- gmm(capmMomentConditionWithTwoLags,data,t0=c(.99,2.0),vcov="iid",wmatrix="ident",type="iterative")
optGmmTwoLag <- gmm(capmMomentConditionWithTwoLags,data,t0=c(.99,2.0),vcov="iid",wmatrix="optimal",type="twoStep")
elResTwoLag <- gel(capmMomentConditionWithTwoLags,data,tet=c(.99,2.0),type="EL")
print(summary(identGmmTwoLag))
print(summary(optGmmTwoLag))
print(summary(elResTwoLag))
