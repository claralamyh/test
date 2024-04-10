# E-M approximation - regular time grid
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


abcrej.sde <- function(N, epsilon, x0, par, T, deltat, yt) {
  theta <- matrix(0, nrow=N, ncol=2)
  theta[1,] <- log(par)
  theta[2:N,1] <- rep(log(par[1]), N-1) + rnorm(N-1) #log(beta)
  theta[2:N,2] <- rep(log(par[2]), N-1) + rnorm(N-1) #log(gamma)
  simdata <- matrix(0, ncol=N, nrow=T/deltat+1)
  zdata <- matrix(0, ncol=N, nrow=length(yt))
  rho <- rep(0,N)
  for(i in 1:N) {
    simdata[,i] <- euler.m(x0, exp(theta[i,]), T, deltat)[,2]
    zdata[,i] <- simdata[c(seq(1,nrow(simdata),0.5/deltat)), i]
    rho[i] <- sqrt(sum(na.omit(zdata[,i]-yt)^2)) #Euclidean distance
  }
  print(mean(rho))
  out <- theta[rho <= epsilon,]
  out.filtered <- out[out[,1]<=0 & out[,2]<=5,]
  return(out.filtered)
}

x0<-c(254,7); par<-c(0.019, 3.11); T<-4; deltat<-0.005 
it <- c(7,14,22,29,20,8,8,NA,0)

set.seed(2212)
abc <- abcrej.sde(N=2000, epsilon=45, x0, par, T, deltat, it) #mean(rho)=77.41851
abc2 <- abcrej.sde(N=2000, epsilon=40, x0, par, T, deltat, it) #80.1617
abc3 <- abcrej.sde(N=2000, epsilon=35, x0, par, T, deltat, it) #81.36512

out <- exp(abc); out2 <- exp(abc2); out3 <- exp(abc3)

par(mfrow=c(3,2))
hist(out[,1], freq=F, main="epsilon=45", xlab="beta")
hist(out[,2], freq=F, main="epsilon=45", xlab="gamma")

hist(out2[,1], freq=F, main="epsilon=40", xlab="beta")
hist(out2[,2], freq=F, main="epsilon=40", xlab="gamma")

hist(out3[,1], freq=F, main="epsilon=35", xlab="beta")
hist(out3[,2], freq=F, main="epsilon=35", xlab="gamma")
par(mfrow=c(1,1))


