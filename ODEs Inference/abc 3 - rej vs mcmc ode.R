###Posterior samples - abc rejection###
out <- exp(abc)

###Posterior samples - abc mcmc###
out2 <- exp(abc2)

# marginal densities
par(mfrow=c(2,2))
plot(density(out[,1]), xlab="beta", main="ABC rejection", ylim=c(0,100))
abline(v=mean(out[,1]))
abline(v=quantile(out[,1],c(0.025,0.975)), col="blue", lty="dashed")

plot(density(out[,2]), xlab="gamma", main="ABC rejection")
abline(v=mean(out[,2]))
abline(v=quantile(out[,2],c(0.025,0.975)), col="blue", lty="dashed")

plot(density(out2[,1]), xlab="beta", main="ABC-MCMC", xlim=c(0.01,0.03))
abline(v=mean(out2[,1]))
abline(v=quantile(out2[,1],c(0.025,0.975)), col="blue", lty="dashed")

plot(density(out2[,2]), xlab="gamma", main="ABC-MCMC", xlim=c(2,5))
abline(v=mean(out2[,2]))
abline(v=quantile(out2[,2],c(0.025,0.975)), col="blue", lty="dashed")

# marginal joint densities
dens1 <- kde2d(out[,1],out[,2])
dens2 <- kde2d(out2[,1],out2[,2])

cowplot::plot_grid(~filled.contour(dens1, nlevels = 6, xlab="beta", ylab="gamma", main="ABC rejection", col=rev(hcl.colors(9, "YlGnBu")),
                                   plot.axes = {
                                     axis(1)
                                     axis(2)
                                     contour(dens1, nlevels = 6, add = TRUE, lwd = 2, labcex=0.8)
                                   }),
                   ~filled.contour(dens2, nlevels = 6, xlab="beta", ylab="gamma", main="ABC-MCMC", col=rev(hcl.colors(9, "YlGnBu")),
                                   plot.axes = {
                                     axis(1)
                                     axis(2)
                                     contour(dens2, nlevels = 6, add = TRUE, lwd = 2, labcex=0.8)
                                   }))

# within sample predictive distribution
endT <- 4
dt <- 0.005

N <- nrow(out)
pred <- array(0, dim=c(endT/dt+1, N, 2))
for(i in 1:N) {
  sim <- euler(c(254,7), out[i,1:2], endT, dt)
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

par(mfrow=c(2,2))
plot(ts(S.mean, start=0, deltat=dt), ylim=c(50,260),
     xlab = "Time (month)", ylab = "Susceptible", main="ABC rejection")
lines(ts(S.lq, start=0, deltat=dt), col="blue", lty="dashed")
lines(ts(S.uq, start=0, deltat=dt),  col="blue",lty="dashed")
lines(seq(0,4,0.5)[-8], na.omit(st), type="l", col="red")

plot(ts(I.mean, start=0, deltat=dt), ylim=c(0,50),
     xlab = "Time (month)", ylab = "Infectious", main="ABC rejection")
lines(ts(I.lq, start=0, deltat=dt), col="blue", lty="dashed")
lines(ts(I.uq, start=0, deltat=dt), col="blue", lty="dashed")
lines(seq(0,4,0.5)[-8], na.omit(it), type="l", col="red")

N <- nrow(out2)
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

I.mat <- pred[,,2]
I.mean <- apply(I.mat, 1, mean)
I.lq <- apply(I.mat, 1, quantile,0.025)
I.uq <- apply(I.mat, 1, quantile,0.975)

plot(ts(S.mean, start=0, deltat=dt), ylim=c(50,260),
     xlab = "Time (month)", ylab = "Susceptible", main="ABC-MCMC")
lines(ts(S.lq, start=0, deltat=dt), col="blue", lty="dashed")
lines(ts(S.uq, start=0, deltat=dt), col="blue", lty="dashed")
lines(seq(0,4,0.5)[-8], na.omit(st), type="l", col="red")

plot(ts(I.mean, start=0, deltat=dt), ylim=c(0,50),
     xlab = "Time (month)", ylab = "Infectious", main="ABC-MCMC")
lines(ts(I.lq, start=0, deltat=dt), col="blue", lty="dashed")
lines(ts(I.uq, start=0, deltat=dt), col="blue", lty="dashed")
lines(seq(0,4,0.5)[-8], na.omit(it), type="l", col="red")
