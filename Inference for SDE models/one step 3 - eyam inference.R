###Basic reproduction number###
R0 <- post[,1]/post[,2]*261
par(mfrow=c(1,3))
hist(R0, xlab="R_0", main="")
abline(v=mean(R0), lwd=1.3) 
# 1.64642
abline(v=quantile(R0, c(0.025,0.975)), col="blue", lty="dashed", lwd=1.3)
# 1.351404 2.021509 


###E-M irregular time-grid###
emsim <- function(x0, time, par, afun, bfun) {
  n <- length(time)
  mat <- matrix(0, nrow=n, ncol=2)
  colnames(mat) <- c("S", "I")
  mat[1,] <- x0
  for(i in 2:n) {
    deltat <- time[i] - time[i-1]
    mat[i,] <- (mat[i-1,]+afun(mat[i-1,],par)*deltat) +
      t(-chol(bfun(mat[i-1,],par)))%*%rnorm(2,0,sqrt(deltat))
    for(j in 1:2) {
      if(mat[i,j] < 0.001) {
        mat[i,j] <- 0.001
      }
    }
  }
  mat
}

# within sample predictive distribution
N <- nrow(post)
pred <- array(0, dim=c(nrow(data), N, 2))
for(i in 1:N) {
  sim <- emsim(data[1,2:3], data[,1], post[i,1:2], alphaSIR, betaSIR)
  pred[,i,1] <- sim[,1]
  pred[,i,2] <- sim[,2]
}

S.mat <- pred[,,1]
S.mean <- apply(S.mat, 1, mean)
S.lq <- apply(S.mat, 1, quantile, 0.025)
S.uq <- apply(S.mat, 1, quantile, 0.975)

I.mat <- pred[,,2]
I.mean <- apply(I.mat, 1, mean)
I.lq <- apply(I.mat, 1, quantile,0.025)
I.uq <- apply(I.mat, 1, quantile,0.975)

plot(data[,1], S.mat[,1], type="l", ylim=c(0,260), col="grey",
     xlab = "Time (month)", ylab = "Susceptible", main="")
for(i in 2:N) {
  lines(data[,1], S.mat[,i], col="grey")}
lines(data[,1], S.mean, lwd=1.3)
lines(data[,1], S.lq, col="blue", lwd=1.3)
lines(data[,1], S.uq, col="blue", lwd=1.3)
lines(data[,1], data[,2], col="red", lwd=1.3)

plot(data[,1], I.mat[,1], type="l", ylim=c(0,100), col="grey",
     xlab = "Time (month)", ylab = "Infectious", main="")
for(i in 2:N) {
  lines(data[,1], I.mat[,i], col="grey")}
lines(data[,1], I.mean, lwd=1.3)
lines(data[,1], I.lq, col="blue", lwd=1.3)
lines(data[,1], I.uq, col="blue", lwd=1.3)
lines(data[,1], data[,3], col="red", lwd=1.3)

