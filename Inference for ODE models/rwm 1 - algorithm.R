library(deSolve)
library(MASS)

###Eyam plague data###
it <- c(7,14,22,29,20,8,8,NA,0)
yt <- it #infectious data
psi0 <- log(c(exp(-4.5), 1.5, 3)) #log(beta0,gamma0,sigma0)
x0 <- c(254,7) #c(S0,I0)


###Solve ODEs###
euler <- function(x0, par, T, deltat) {
  fct <- function(t, x0, par) {
    with(as.list(c(x0, par)), {
      dS <- -par[1]*x0[1]*x0[2]
      dI <- par[1]*x0[1]*x0[2] - par[2]*x0[2]
      list(c(dS, dI))
    })
  }
  times <- seq(0, T, by=deltat)
  out <- ode(y=x0, times=times, func=fct, parms=par)
  R <- c(rep(sum(x0), nrow(out)) - out[,2] - out[,3])
  mat <- cbind(out, R)
}


###Evaluate log likelihood###
loglike <- function(x0, psi, yt, T, deltat) {
  sim <- euler(x0, exp(psi[1:2]), T, deltat)
  It <- sim[c(seq(1,nrow(sim),0.5/deltat)), 3]
  sum(dnorm(yt,It,exp(psi[3]),log=TRUE), na.rm=TRUE)    
}


###Evaluate the log of posterior distribution (up to normalising constant)###
lpost <- function(x0, psi, yt, T, deltat) {
  lprior <- sum(dnorm(psi, log=TRUE))
  llike <- loglike(x0, psi, yt, T, deltat) 
  lprior + llike
}


###RWM Function###
RWM <- function(N, Vtune, psi0, x0, yt, T, deltat) {
  mat <- matrix(0, nrow=N, ncol=3)
  psi <- psi0
  mat[1,] <- psi
  count <- 0 
  for (i in 2:N) {
    innov <- mvrnorm(1, rep(0,3), Vtune) #innovation
    can <- psi + innov #proposed candidate
    laprob <- lpost(x0,can,yt,T,deltat) - lpost(x0,psi,yt,T,deltat)
    if (log(runif(1)) < laprob) { #accept 
      psi <- can 
      count <- count+1
    } 
    mat[i,] <- psi #reject
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


###Executing RWN###
set.seed(0105)
wrapper(
  out <- RWM(N=2000, Vtune=0.07^2*diag(3), psi0, x0, yt, 4, 0.005))
# pilot run: 0.3121561, Time difference of 3.030774 mins

set.seed(2112)
wrapper(
  out2 <- RWM(N=5000, Vtune=1.2^2*var(out[-(1:400),]), psi0, x0, yt, 4, 0.005))
# main run: 0.3080616, Time difference of 2.944661  mins
