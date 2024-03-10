###Posterior samples###
post <- exp(out2[-(1:100),])

# basic reproduction number
R0 <- post[,1]/post[,2]*261
par(mfrow=c(1,1))
hist(R0, xlab="Basic reproduction number", main="")
abline(v=mean(R0), lwd=1.5)
abline(v=quantile(R0, c(0.025,0.975)), col="blue", lty="dashed", lwd=1.5)

# within sample predictive distribution
N <- nrow(post)
endT <- 4
dt <- 0.005
pred <- array(0, dim=c(endT/dt+1, N, 2))
for(i in 1:N) {
  sim <- euler(c(254,7), post[i,1:2], endT, dt)
  pred[,i,1] <- sim[,2]
  pred[,i,2] <- sim[,3]
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

par(mfrow=c(1,2))
plot(ts(S.mean, start=0, deltat=dt), ylim=c(50,260),
     xlab = "Time (month)", ylab = "Susceptible")
lines(ts(S.lq, start=0, deltat=dt), col="blue", lty="dashed")
lines(ts(S.uq, start=0, deltat=dt), col="blue", lty="dashed")
lines(seq(0,4,0.5)[-8], na.omit(st), type="l", col="red")

plot(ts(I.mean, start=0, deltat=dt), ylim=c(0,32),
     xlab = "Time (month)", ylab = "Infectious")
lines(ts(I.lq, start=0, deltat=dt), col="blue", lty="dashed")
lines(ts(I.uq, start=0, deltat=dt), col="blue", lty="dashed")
lines(seq(0,4,0.5)[-8], na.omit(it), type="l", col="red")