###Inference###
post <- exp(out[101:5000,])
effectiveSize(post) #728.1088 754.5061

# summaries
mean(post[,1]) #0.01903765
sd(post[,1]) #0.001389289
quantile(post[,1],c(0.025,0.975)) #0.01646219 0.02193448

mean(post[,2]) #3.033178
sd(post[,2]) #0.2162525
quantile(post[,2],c(0.025,0.975)) #2.640886 3.456923

# sampler mixing and summaries
par(mfrow=c(2,3))
plot(ts(post[,1], start=0), xlab="Iteration", ylab="beta")
acf(post[,1], main="beta")
plot(density(post[,1]), xlab="beta", main="")
abline(v=mean(post[,1]))
abline(v=quantile(post[,1],c(0.025,0.975)), col="blue", lty="dashed")

plot(ts(post[,2], start=0), xlab="Iteration", ylab="gamma")
acf(post[,2], main="gamma") 
plot(density(post[,2]), xlab="gamma", main="")
abline(v=mean(post[,2]))
abline(v=quantile(post[,2],c(0.025,0.975)), col="blue", lty="dashed")

# kernel density estimate (joint density)
dens <- kde2d(post[,1],post[,2])
par(mfrow=c(1,1))
filled.contour(dens, nlevels=6, xlab="beta", ylab="gamma", col=rev(hcl.colors(9, "YlGnBu")),
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(dens, nlevels=6, add=TRUE, lwd=1.8, labcex=0.8)
               })


