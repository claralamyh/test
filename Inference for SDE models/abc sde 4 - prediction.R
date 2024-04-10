par(mfrow=c(1,2))

##out
# basic reproduction number
R0 <- out[,1]/out[,2]*261
hist(R0, xlab="R_0", main="ABC rejection", xlim=c(0.5,3))
abline(v=mean(R0), lwd=1.3) # 1.520872
abline(v=quantile(R0, c(0.025,0.975)), 
       col="blue", lty="dashed", lwd=1.3) #0.9068752 2.1034270 

# within sample predictive distribution
endT <- 4
dt <- 0.005

pred <- array(0, dim=c(endT/dt+1, nrow(out), 2))
for(i in 1:nrow(out)) {
  sim <- euler.m(c(254,7), out[i,], endT, dt)
  pred[,i,1] <- sim[,1]
  pred[,i,2] <- sim[,2]
}

thinned <- out2[seq(1,nrow(out2),10),]
pred2 <- array(0, dim=c(endT/dt+1, nrow(thinned), 2))
for(i in 1:nrow(thinned)) {
  sim <- euler.m(c(254,7), thinned[i,], endT, dt)
  pred2[,i,1] <- sim[,1]
  pred2[,i,2] <- sim[,2]
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

plot(ts(S.mat[,1], start=0, deltat=dt), type="l", ylim=c(0,260), col="grey",
     xlab="Time (month)", ylab="Susceptible", main="ABC rejection")
for(i in 2:nrow(out)) {
  lines(ts(S.mat[,i], start=0, deltat=dt), col="grey")}
lines(ts(S.mean, start=0, deltat=dt), lwd=1.3)
lines(ts(S.lq, start=0, deltat=dt), col="blue", lwd=1.3)
lines(ts(S.uq, start=0, deltat=dt), col="blue", lwd=1.3)
lines(seq(0,4,0.5)[-8], na.omit(st), col="red", lwd=1.3)

plot(ts(I.mat[,1], start=0, deltat=dt), type="l", ylim=c(0,100), col="grey",
     xlab="Time (month)", ylab="Infectious", main="ABC rejection")
for(i in 2:nrow(out)) {
  lines(ts(I.mat[,i], start=0, deltat=dt), col="grey")}
lines(ts(I.mean, start=0, deltat=dt), lwd=1.3)
lines(ts(I.lq, start=0, deltat=dt), col="blue", lwd=1.3)
lines(ts(I.uq, start=0, deltat=dt), col="blue", lwd=1.3)
lines(seq(0,4,0.5)[-8], na.omit(it), col="red", lwd=1.3)

##out2
R0 <- thinned[,1]/thinned[,2]*261
hist(R0, xlab="R_0", main="ABC-MCMC", xlim=c(0.5,3))
abline(v=mean(R0), lwd=1.3) #1.810282
abline(v=quantile(R0, c(0.025,0.975)), 
       col="blue", lty="dashed", lwd=1.3) #1.273021 2.358940 

S.mat <- pred2[,,1]
S.mean <- apply(S.mat, 1, mean)
S.lq <- apply(S.mat, 1, quantile, 0.025)
S.uq <- apply(S.mat, 1, quantile, 0.975)

I.mat <- pred2[,,2]
I.mean <- apply(I.mat, 1, mean)
I.lq <- apply(I.mat, 1, quantile,0.025)
I.uq <- apply(I.mat, 1, quantile,0.975)

plot(ts(S.mat[,1], start=0, deltat=dt), type="l", ylim=c(0,260), col="grey",
     xlab="Time (month)", ylab="Susceptible", main="ABC-MCMC")
for(i in 2:nrow(pred2)) {
  lines(ts(S.mat[,i], start=0, deltat=dt), col="grey")}
lines(ts(S.mean, start=0, deltat=dt), lwd=1.3)
lines(ts(S.lq, start=0, deltat=dt), col="blue", lwd=1.3)
lines(ts(S.uq, start=0, deltat=dt), col="blue", lwd=1.3)
lines(seq(0,4,0.5)[-8], na.omit(st), col="red", lwd=1.3)

plot(ts(I.mat[,1], start=0, deltat=dt), type="l", ylim=c(0,100), col="grey",
     xlab="Time (month)", ylab="Infectious", main="ABC-MCMC")
for(i in 2:nrow(pred2)) {
  lines(ts(I.mat[,i], start=0, deltat=dt), col="grey")}
lines(ts(I.mean, start=0, deltat=dt), lwd=1.3)
lines(ts(I.lq, start=0, deltat=dt), col="blue", lwd=1.3)
lines(ts(I.uq, start=0, deltat=dt), col="blue", lwd=1.3)
lines(seq(0,4,0.5)[-8], na.omit(it), col="red", lwd=1.3)
