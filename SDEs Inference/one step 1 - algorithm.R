library(MASS)
library(mvtnorm)
library(coda)

###Data###
data <- cbind(time = seq(0,4,0.5)[-8], 
              st = c(254,235,201,153,121,110,97,83),
              it = c(7,14,22,29,20,8,8,0))

par <- log(c(exp(-4.5), 1.5))


###SIR drift and variance functions###
alphaSIR <- function(x, theta) {
  c(-theta[1]*x[1]*x[2],
    theta[1]*x[1]*x[2]-theta[2]*x[2])
}

betaSIR <- function(x, theta) {
  cov <- -theta[1]*x[1]*x[2]
  var1 <- theta[1]*x[1]*x[2]
  var2 <- theta[1]*x[1]*x[2] + theta[2]*x[2]
  matrix(c(var1,cov,cov,var2), ncol=2, byrow=T)
}


###Log likelihood/ joint density of data###
loglike <- function(data, par, afun, bfun) {
  val <- rep(0, nrow(data)-1)
  for(i in 2:nrow(data)) {
    deltat <- data[i,1] - data[i-1,1]
    val[i-1] <- dmvnorm(x = data[i,2:3],
                        mean = data[i-1,2:3]+afun(data[i-1,2:3],exp(par))*deltat,
                        sigma = bfun(data[i-1,2:3],exp(par))*deltat, log = TRUE)
  }
  sum(val)
}


###Log posterior density###
lpost <- function(data, par, afun, bfun) {
  lprior <- sum(dnorm(par, log=TRUE))
  llike <- loglike(data, par, afun, bfun)
  lprior + llike
}


###RWM###
RWM <- function(N, Vtune, data, par, afun, bfun) {
  mat <- matrix(0, nrow=N, ncol=2)
  psi <- par
  mat[1,] <- psi
  count <- 0 
  for (i in 2:N) {
    can <- psi + mvrnorm(1, rep(0,2), Vtune) #proposed candidate
    laprob <- lpost(data,can,afun,bfun) - lpost(data,psi,afun,bfun)
    if (log(runif(1)) < laprob) {
      psi <- can 
      count <- count+1
    } 
    mat[i,] <- psi
  }
  print(count/(N-1)) #empirical overall acceptance rate
  return(mat) 
}


###Measuring the CPU time used###
wrapper <- function(code){
  start <- Sys.time()
  code
  end <- Sys.time()
  end-start
}


###Executing RWM###
set.seed(1417)
wrapper(
  out <- RWM(N=5000, Vtune=0.12^2*diag(2), data, par, alphaSIR, betaSIR))
# 0.3784757, Time difference of 22.61751 secs

