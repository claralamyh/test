###Posterior samples###
post <- exp(out2[-(1:100),])

# marginal summaries
mean(post[,1]) #0.01934984
sd(post[,1]) #0.0008983461
quantile(post[,1],c(0.025,0.975)) #0.01751695 0.02121653

mean(post[,2]) #3.119765
sd(post[,2]) #0.14217
quantile(post[,2],c(0.025,0.975)) #2.850083 3.427365

mean(post[,3]) #2.584337
sd(post[,3]) #0.8160047
quantile(post[,3],c(0.025,0.975)) #1.516564 4.579346

# marginal posterior densities
par(mfrow=c(1,3))
plot(density(post[,1]), xlab="beta", main="")
abline(v=mean(post[,1]))
abline(v=quantile(post[,1],c(0.025,0.975)), col="blue", lty="dashed")

plot(density(post[,2]), xlab="gamma", main="")
abline(v=mean(post[,2]))
abline(v=quantile(post[,2],c(0.025,0.975)), col="blue", lty="dashed")

plot(density(post[,3]), xlab="sigma", main="")
abline(v=mean(post[,3]))
abline(v=quantile(post[,3],c(0.025,0.975)), col="blue", lty="dashed")

# marginal posterior joint density for beta and gamma
dens <- kde2d(post[,1],post[,2])
par(mfrow=c(1,1))
filled.contour(dens, nlevels=6, xlab="beta", ylab="gamma", col=rev(hcl.colors(9, "YlGnBu")),
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(dens, nlevels=6, add=TRUE, lwd=1.8, labcex=0.8)
               })

