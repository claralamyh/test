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

abcrej <- function(N, epsilon, x0, par, T, deltat, yt) {
  theta <- matrix(0, nrow=N, ncol=2)
  theta[1,] <- log(par)
  theta[2:N,1] <- rep(log(par[1]), N-1) + rnorm(N-1) #log(beta)
  theta[2:N,2] <- rep(log(par[2]), N-1) + rnorm(N-1) #log(gamma)
  simdata <- matrix(0, ncol=N, nrow=T/deltat+1)
  zdata <- matrix(0, ncol=N, nrow=length(yt))
  rho <- rep(0,N)
  for(i in 1:N) {
    simdata[,i] <- euler(x0, exp(theta[i,]), T, deltat)[,2]
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

set.seed(2002)
abc <- abcrej(N=2000, epsilon=45, x0, par, T, deltat, it) #mean(rho)=75.43174
abc2 <- abcrej(N=2000, epsilon=40, x0, par, T, deltat, it) #78.88682
abc3 <- abcrej(N=2000, epsilon=35, x0, par, T, deltat, it) #79.91357

out <- exp(abc); out2 <- exp(abc2); out3 <- exp(abc3)

par(mfrow=c(3,2))
hist(out[,1], freq=F, main="epsilon=45", xlab="beta")
hist(out[,2], freq=F, main="epsilon=45", xlab="gamma")

hist(out2[,1], freq=F, main="epsilon=40", xlab="beta")
hist(out2[,2], freq=F, main="epsilon=40", xlab="gamma")

hist(out3[,1], freq=F, main="epsilon=35", xlab="beta")
hist(out3[,2], freq=F, main="epsilon=35", xlab="gamma")


N <- nrow(out2)
endT <- 4
dt <- 0.005
pred <- array(0, dim=c(endT/dt+1, N, 2))
for(i in 1:N) {
  sim <- euler(c(254,7), out2[i,1:2], endT, dt)
  pred[,i,1] <- sim[,1]
  pred[,i,2] <- sim[,2]
}

S.mat <- pred[,,1]
S.mean <- apply(S.mat, 1, mean)
S.lq <- apply(S.mat, 1, quantile, 0.025)
S.uq <- apply(S.mat, 1, quantile, 0.975)
st <- c(254,235,201,153,121,110,97,NA,83)

I.mat <- pred[,,2]
I.mean <- apply(I.mat, 1, mean)
I.lq <- apply(I.mat, 1, quantile,0.025)
I.uq <- apply(I.mat, 1, quantile,0.975)
it <- c(7,14,22,29,20,8,8,NA,0)

par(mfrow=c(1,2))
plot(ts(S.mean, start=0, deltat=0.005), ylim=c(30,260),
     xlab = "Time (month)", ylab = "Susceptible", main="")
lines(ts(S.lq, start=0, deltat=0.005))
lines(ts(S.uq, start=0, deltat=0.005))
lines(seq(0,4,0.5)[-8], na.omit(st), type="l", col="red")

plot(ts(I.mean, start=0, deltat=0.005), ylim=c(0,60),
     xlab = "Time (month)", ylab = "Infectious", main="")
lines(ts(I.lq, start=0, deltat=0.005))
lines(ts(I.uq, start=0, deltat=0.005))
lines(seq(0,4,0.5)[-8], na.omit(it), type="l", col="red")

