library(MASS)
library(coda)

###E-M approximation - regular time grid###
euler.m <- function(x0, par, T, deltat) {
  n <- T/deltat + 1
  mat <- matrix(0, nrow=n, ncol=2)
  colnames(mat) <- c("S", "I")
  mat[1,] <- x0
  for(i in 2:n) {
    h1 <- par[1]*mat[i-1,1]*mat[i-1,2]
    h2 <- par[2]*mat[i-1,2]
    infmean <- c(-h1, h1-h2)
    infsd <- t(-chol(matrix(c(h1,-h1,-h1,h1+h2), ncol=2, byrow=T)))
    mat[i,] <- (mat[i-1,]+infmean*deltat) + infsd%*%rnorm(2,0,sqrt(deltat))
    for(j in 1:2) {
      if(mat[i,j] < 0.001) {
        mat[i,j] <- 0.001
      }
    }
  }
  mat
}


###ABC MCMC###
abcmcmc <- function(N, epsilon, x0, psi0, T, deltat, yt, Vtune) {
  mat <- matrix(0, nrow=N, ncol=2) #store values
  psi <- psi0
  mat[1,] <- psi
  count <- 0 
  rho <- rep(0,N)
  for (i in 2:N) {
    innov <- mvrnorm(1, rep(0,2), Vtune) #innovation
    can <- psi + innov #proposed candidate
    simdata <- euler.m(x0, exp(can), T, deltat)
    zdata <- simdata[c(seq(1,nrow(simdata),0.5/deltat)),2]
    rho[i] <- sqrt(sum(na.omit(zdata-yt)^2)) #Euclidean distance
    if (rho[i] <= epsilon) {
      laprob <- sum(dnorm(can, log=TRUE)) - sum(dnorm(psi, log=TRUE))
      if (log(runif(1)) < laprob) { #accept 
        psi <- can
        count <- count+1
      } 
    }
    mat[i,] <- psi
  }
  mat
}


###Measuring the CPU time used###
wrapper <- function(code) {
  start <- Sys.time()
  code
  end <- Sys.time()
  end-start
}


###Executing ABC-MCMC###
x0<-c(254,7); psi<-log(c(0.019, 3.12)); T<-4; deltat<-0.005 
it <- c(7,14,22,29,20,8,8,NA,0)

set.seed(2307)
wrapper(
  abc <- abcmcmc(N=2000, epsilon=25, x0=x0, psi0=psi, 4, 0.005, it, Vtune=0.3^2*diag(2)))
# Time difference of 55.35868 secs

set.seed(1332)
wrapper(
  abc2 <- abcmcmc(N=10000, epsilon=25, x0=x0, psi0=psi, 4, 0.005, it, Vtune=1.5^2*var(abc)))
# Time difference of 5.738656 mins


###Posterior samples - abc mcmc###
out2 <- exp(abc2)

# ESS
effectiveSize(exp(abc)) #31.23650 41.32112
effectiveSize(out2) #142.9984 232.7711 

# summaries
mean(out2[,1]) #0.02779104
sd(out2[,1]) #0.008777211
quantile(out2[,1],c(0.025,0.975)) #0.01373195 0.05775493 

mean(out2[,2]) #4.025225
sd(out2[,2]) #1.088616
quantile(out2[,2],c(0.025,0.975)) #2.269341 6.790416 

# mixing
par(mfrow=c(2,2))
plot(ts(out2[,1]), xlab="Iteration", ylab="beta", main="")
plot(ts(out2[,2]), xlab="Iteration", ylab="gamma", main="")

acf(out2[,1], main="beta")
acf(out2[,2], main="gamma")



