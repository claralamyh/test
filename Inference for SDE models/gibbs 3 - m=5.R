library(coda)

# m=5
set.seed(2338)
out5 <- gibbs(N=5000, m=5, Vtune=0.12^2*diag(2), data, par, alphaSIR, betaSIR)

post5 <- out5[[1]][-(1:100),]
effectiveSize(post5) #280.6367 279.2175 


# summaries
post <- post5
mean(post[,1]) #0.01976993
sd(post[,1]) #0.001820084
quantile(post[,1],c(0.025,0.975)) #0.01642087 0.02356089

mean(post[,2]) #3.124241
sd(post[,2]) #0.2763269
quantile(post[,2],c(0.025,0.975)) #2.592993 3.698658


# kernel density estimate (marginal density)
plot(density(post5[,1]), xlab="beta", main="m=5", ylim=c(0,250))
abline(v=mean(post5[,1]))
abline(v=quantile(post5[,1],c(0.025,0.975)), col="blue", lty="dashed")

plot(density(post5[,2]), xlab="gamma", main="m=5", ylim=c(0,1.5))
abline(v=mean(post5[,2]))
abline(v=quantile(post5[,2],c(0.025,0.975)), col="blue", lty="dashed")


# kernel density estimate (joint density)
dens <- kde2d(post5[,1],post5[,2])
par(mfrow=c(1,1))
filled.contour(dens, nlevels=6, main="m=5", xlab="beta", ylab="gamma", col=rev(hcl.colors(9, "YlGnBu")),
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(dens, nlevels=6, add=TRUE, lwd=1.8, labcex=0.8)
               })


#within-sample predictive distribution
st5 <- out5[[2]][-(1:5),]; it5 <- out5[[3]][-(1:5),]
plot(out5[[5]], st5[1,], type="l", ylim=c(80,260),
     xlab="Time (month)", ylab="Susceptible", main="m=5")
for(i in 2:nrow(st5)) {
  lines(out5[[5]], st5[i,], col=i)
}
plot(out5[[5]], it5[1,], type="l", ylim=c(0,50),
     xlab="Time (month)", ylab="Infectious", main="m=5")
for(i in 2:nrow(it5)) {
  lines(out5[[5]], it5[i,], col=i)
}



###Basic reproduction number###
R0 <- post5[,1]/post5[,2]*261
hist(R0, xlab="R_0", main="m=5")
abline(v=mean(R0), lwd=1.3) 
# 1.661333
abline(v=quantile(R0, c(0.025,0.975)), col="blue", lty="dashed", lwd=1.3)


# latent process marginal
N <- nrow(post5)
pred5 <- array(0, dim=c(length(out5[[5]]), N, 2))
for(i in 1:N) {
  sim <- emsim(data[1,2:3], out5[[5]], post5[i,1:2], alphaSIR, betaSIR)
  pred5[,i,1] <- sim[,2]
  pred5[,i,2] <- sim[,3]
}

S.mat <- pred5[,,1]
S.mean <- apply(S.mat, 1, mean)
S.lq <- apply(S.mat, 1, quantile, 0.025)
S.uq <- apply(S.mat, 1, quantile, 0.975)

I.mat <- pred5[,,2]
I.mean <- apply(I.mat, 1, mean)
I.lq <- apply(I.mat, 1, quantile,0.025)
I.uq <- apply(I.mat, 1, quantile,0.975)

plot(out5[[5]], S.mat[,1], type="l", ylim=c(0,260), col="grey",
     xlab = "Time (month)", ylab = "Susceptible", main="m=5")
for(i in 2:N) {
  lines(out5[[5]], S.mat[,i], col="grey")}
lines(out5[[5]], S.mean, lwd=1.3)
lines(out5[[5]], S.lq, col="blue", lwd=1.3)
lines(out5[[5]], S.uq, col="blue", lwd=1.3)
lines(data[,1], data[,2], col="red", lwd=1.3)

plot(out5[[5]], I.mat[,1], type="l", ylim=c(0,100), col="grey",
     xlab = "Time (month)", ylab = "Infectious", main="m=5")
for(i in 2:N) {
  lines(out5[[5]], I.mat[,i], col="grey")}
lines(out5[[5]], I.mean, lwd=1.3)
lines(out5[[5]], I.lq, col="blue", lwd=1.3)
lines(out5[[5]], I.uq, col="blue", lwd=1.3)
lines(data[,1], data[,3], col="red", lwd=1.3)


