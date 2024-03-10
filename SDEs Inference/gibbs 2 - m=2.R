library(coda)

# m=2
set.seed(0225)
out2 <- gibbs(N=5000, m=2, Vtune=0.12^2*diag(2), data, par, alphaSIR, betaSIR)

post2 <- out2[[1]][-(1:100),]
effectiveSize(post2) #492.3801 413.1365 


# summaries
post <- post2
mean(post[,1]) #0.01916658
sd(post[,1]) #0.001674287
quantile(post[,1],c(0.025,0.975)) #0.01602585 0.02268682 

mean(post[,2]) #3.074714
sd(post[,2]) #0.2630541
quantile(post[,2],c(0.025,0.975)) #2.591748 3.622998


# kernel density estimate (marginal density)
par(mfrow=c(2,2))
plot(density(post2[,1]), xlab="beta", main="m=2")
abline(v=mean(post2[,1]))
abline(v=quantile(post2[,1],c(0.025,0.975)), col="blue", lty="dashed")

plot(density(post2[,2]), xlab="gamma", main="m=2")
abline(v=mean(post2[,2]))
abline(v=quantile(post2[,2],c(0.025,0.975)), col="blue", lty="dashed")


# kernel density estimate (joint density)
dens <- kde2d(post2[,1],post2[,2])
par(mfrow=c(1,1))
filled.contour(dens, nlevels=6, main="m=2", xlab="beta", ylab="gamma", col=rev(hcl.colors(9, "YlGnBu")),
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(dens, nlevels=6, add=TRUE, lwd=1.8, labcex=0.8)
               })


#within-sample predictive distribution
st2 <- out2[[2]][-(1:100),]; it2 <- out2[[3]][-(1:100),]
par(mfrow=c(2,2))
plot(out2[[5]], st2[1,], type="l", ylim=c(80,260),
     xlab="Time (month)", ylab="Susceptible", main="m=2")
for(i in 2:nrow(st2)) {
  lines(out2[[5]], st2[i,], col=i)
}
plot(out2[[5]], it2[1,], type="l", ylim=c(0,50),
     xlab="Time (month)", ylab="Infectious", main="m=2")
for(i in 2:nrow(it2)) {
  lines(out2[[5]], it2[i,], col=i)
}


###Basic reproduction number###
par(mfrow=c(2,3))
R0 <- post2[,1]/post2[,2]*261
hist(R0, xlab="R_0", main="m=2")
abline(v=mean(R0), lwd=1.3) 
# 1.635866
abline(v=quantile(R0, c(0.025,0.975)), col="blue", lty="dashed", lwd=1.3)


#out of-sample predictive distribution
N <- nrow(post2)
pred2 <- array(0, dim=c(length(out2[[5]]), N, 2))
for(i in 1:N) {
  sim <- emsim(data[1,2:3], out2[[5]], post2[i,1:2], alphaSIR, betaSIR)
  pred2[,i,1] <- sim[,2]
  pred2[,i,2] <- sim[,3]
}

S.mat <- pred2[,,1]
S.mean <- apply(S.mat, 1, mean)
S.lq <- apply(S.mat, 1, quantile, 0.025)
S.uq <- apply(S.mat, 1, quantile, 0.975)

I.mat <- pred2[,,2]
I.mean <- apply(I.mat, 1, mean)
I.lq <- apply(I.mat, 1, quantile,0.025)
I.uq <- apply(I.mat, 1, quantile,0.975)

plot(out2[[5]], S.mat[,1], type="l", ylim=c(0,260), col="grey",
     xlab = "Time (month)", ylab = "Susceptible", main="m=2")
for(i in 2:N) {
  lines(out2[[5]], S.mat[,i], col="grey")}
lines(out2[[5]], S.mean, lwd=1.3)
lines(out2[[5]], S.lq, col="blue", lwd=1.3)
lines(out2[[5]], S.uq, col="blue", lwd=1.3)
lines(data[,1], data[,2], col="red", lwd=1.3)

plot(out2[[5]], I.mat[,1], type="l", ylim=c(0,100), col="grey",
     xlab = "Time (month)", ylab = "Infectious", main="m=2")
for(i in 2:N) {
  lines(out2[[5]], I.mat[,i], col="grey")}
lines(out2[[5]], I.mean, lwd=1.3)
lines(out2[[5]], I.lq, col="blue", lwd=1.3)
lines(out2[[5]], I.uq, col="blue", lwd=1.3)
lines(data[,1], data[,3], col="red", lwd=1.3)

