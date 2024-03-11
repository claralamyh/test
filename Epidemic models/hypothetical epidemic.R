library(MASS)
library(mvtnorm)
library(coda)
library(deSolve)

###Data###
data <- cbind(time = seq(0,4,0.5)[-8], 
              st = c(254,235,201,153,121,110,97,83),
              it = c(7,14,22,29,20,8,8,0))

data.seen <- data[1:6,]
par <- log(c(0.01916658, 3.074714))


###MCMC with data augmentation### - "eyam sim 1 - algorithm.R"
set.seed(1740)
out2 <- gibbs(N=5000, m=2, Vtune=0.12^2*diag(2), data.seen, par, alphaSIR, betaSIR)
post2 <- out2[[1]][-(1:100),]
effectiveSize(post2)
plot(ts(post2))

new.par <- c(mean(post2[,1]),mean(post2[,2])) #0.01813936 3.02533917


###Prediction###
euler <- function(x0, par, T, deltat) {
  fct <- function(t, x0, par) {
    with(as.list(c(x0, par)), {
      dS <- -par[1]*x0[1]*x0[2]
      dI <- par[1]*x0[1]*x0[2] - par[2]*x0[2]
      list(c(dS, dI))
    })
  }
  times <- seq(0, T, by=deltat)
  ode(y=x0, times=times, func=fct, parms=par)
}

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

x0.unseen <- data[6,2:3]
(R_0 <- new.par[1]/new.par[2]*261) #1.564906
1/R_0 #0.6390159
x0.unseen[1]/261 #0.4214559  

ode <- euler(x0.unseen, new.par, 1.5, 0.005)

m <- 50
outarray <- array(0, c(nrow(ode), m, 2))

set.seed(1538)
for (i in 1:m) {
  sim <- euler.m(x0.unseen, new.par, 1.5, 0.005)
  outarray[,i,1] <- sim[,1]
  outarray[,i,2] <- sim[,2]
}
S.mat <- outarray[,,1]; I.mat <- outarray[,,2]

T <- 2.5
plot(ts(S.mat[,1], start=T, deltat=0.005), ylim=c(0,260), xlim=c(0,4),
     xlab = "Time (month)", ylab = "Population", col="grey")
for(i in 2:m) {
  lines(ts(S.mat[,i], start=T, deltat=0.005), col="grey")
}
for(i in 1:m) {
  lines(ts(I.mat[,i], start=T, deltat=0.005), col="grey")
  lines(ts(rep(261,nrow(ode))-S.mat[,i]-I.mat[,i], start=T, deltat=0.005), col="grey")
}

points(data[,1], data[,2])
points(data[,1], data[,3], pch=2)

lines(ts(apply(S.mat, 1, mean), start=T, deltat=0.005), col="green")
lines(ts(apply(I.mat, 1, mean), start=T, deltat=0.005), col="red")
lines(ts(rep(261,nrow(ode))-apply(S.mat,1,mean)-apply(I.mat,1,mean),
         start=T, deltat=0.005), col="blue")

lines(ts(ode[,2], start=T, deltat=0.005), col="green", lty=5, lwd=1.5)
lines(ts(ode[,3], start=T, deltat=0.005), col="red", lty=5, lwd=1.5)
lines(ts(rep(261,nrow(ode))-ode[,2]-ode[,3], start=T, deltat=0.005),
      col="blue", lty=5, lwd=1.5)
