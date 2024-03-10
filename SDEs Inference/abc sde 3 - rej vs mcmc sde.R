###Posterior samples - abc rejection###
out <- exp(abc)

###Posterior samples - abc mcmc###
out2 <- exp(abc2)


# marginal densities
par(mfrow=c(2,2))
plot(density(out[,1]), xlab="beta", main="ABC rejection", xlim=c(0,0.07))
abline(v=mean(out[,1]))
abline(v=quantile(out[,1],c(0.025,0.975)), col="blue", lty="dashed")

plot(density(out[,2]), xlab="gamma", main="ABC rejection", xlim=c(0,8))
abline(v=mean(out[,2]))
abline(v=quantile(out[,2],c(0.025,0.975)), col="blue", lty="dashed")

plot(density(out2[,1]), xlab="beta", main="ABC-MCMC", xlim=c(0,0.07), ylim=c(0,60))
abline(v=mean(out2[,1]))
abline(v=quantile(out2[,1],c(0.025,0.975)), col="blue", lty="dashed")

plot(density(out2[,2]), xlab="gamma", main="ABC-MCMC", xlim=c(0,8))
abline(v=mean(out2[,2]))
abline(v=quantile(out2[,2],c(0.025,0.975)), col="blue", lty="dashed")


# marginal joint densities
dens1 <- kde2d(out[,1],out[,2])
dens2 <- kde2d(out2[,1],out2[,2])

cowplot::plot_grid(~filled.contour(dens1, nlevels = 6, xlab="beta", ylab="gamma", main="ABC rejection", col=rev(hcl.colors(9, "YlGnBu")),
                                   plot.axes = {
                                     axis(1)
                                     axis(2)
                                     contour(dens1, nlevels = 6, add = TRUE, lwd = 1.8, labcex=0.8)
                                   }),
                   ~filled.contour(dens2, nlevels = 6, xlab="beta", ylab="gamma", main="ABC-MCMC", col=rev(hcl.colors(9, "YlGnBu")),
                                   plot.axes = {
                                     axis(1)
                                     axis(2)
                                     contour(dens2, nlevels = 6, add = TRUE, lwd = 1.8, labcex=0.8)
                                   }))
