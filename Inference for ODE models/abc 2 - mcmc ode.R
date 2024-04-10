library(MASS)
library(coda)

###Solve ODEs###
euler <- function(x0, par, T, deltat) {
  n <- T/deltat + 1
  mat <- matrix(0, nrow=n, ncol=2)
  colnames(mat) <- c("S", "I")
  mat[1,] <- c(x0[1], x0[2])
  for(i in 2:n) {
    St <- mat[i-1,1]
    It <- mat[i-1,2]
    mat[i,1] <- St - par[1]*St*It*deltat
    mat[i,2] <- It + par[1]*St*It*deltat - par[2]*It*deltat
  }
  mat
}


###ABC MCMC###
abcmcmc <- function(N, epsilon, x0, psi0, T, deltat, yt, Vtune) {
  mat <- matrix(0, nrow=N, ncol=2) #store values
  psi <- psi0
  mat[1,] <- psi
  for (i in 2:N) {
    can <- psi + mvrnorm(1, rep(0,2), Vtune) #proposed candidate
    simdata <- euler(x0, exp(can), T, deltat)
    zdata <- simdata[c(seq(1,nrow(simdata),0.5/deltat)), 2]
    rho <- sqrt(sum(na.omit(zdata-yt)^2)) #Euclidean distance
    if (rho <= epsilon) {
       laprob <- sum(dnorm(can, log=TRUE)) - sum(dnorm(psi, log=TRUE))
       if (log(runif(1)) < laprob) {
         psi <- can
       } 
    }
    mat[i,] <- psi
  }
  mat
}


###Measuring the CPU time used###
wrapper <- function(code){
  start <- Sys.time()
  code
  end <- Sys.time()
  end-start
}


###Executing ABC-MCMC###
x0<-c(254,7); psi<-log(c(0.019, 3.12)); T<-4; deltat<-0.005 
st <- c(254,235,201,153,121,110,97,NA,83)
it <- c(7,14,22,29,20,8,8,NA,0)

set.seed(1329)
wrapper( 
  abc <- abcmcmc(N=5000, epsilon=20, x0=x0, psi0=psi, 4, 0.005, it, Vtune=0.2^2*diag(2)))
# Time difference of 31.99775 secs

set.seed(2216)
wrapper( 
  abc2 <- abcmcmc(N=11200, epsilon=20, x0=x0, psi0=psi, 4, 0.005, it, Vtune=1.3^2*var(abc)))
# Time difference of 1.186351 mins


###Posterior samples - abc mcmc###
out2 <- exp(abc2)

# ESS
effectiveSize(out2) #1003.124 1093.615

# summaries
mean(out2[,1]) #0.0211898
sd(out2[,1]) #0.003393583
quantile(out2[,1],c(0.025,0.975)) #0.01432931 0.02647355 

mean(out2[,2]) #3.413172
sd(out2[,2]) #0.5497504
quantile(out2[,2],c(0.025,0.975)) #2.396206 4.360882 

# mixing
par(mfrow=c(2,2))
plot(ts(out2[,1]), xlab="Iteration", ylab="beta", main="")
plot(ts(out2[,2]), xlab="Iteration", ylab="gamma", main="")

acf(out2[,1], main="beta")
acf(out2[,2], main="gamma")

