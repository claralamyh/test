library(coda)

###Posterior samples###
post <- exp(out2[-(1:100),])

# mixing
par(mfrow=c(2,3))
plot(ts(post[,1]), xlab="Iteration", ylab="beta", main="beta")
plot(ts(post[,2]), xlab="Iteration", ylab="gamma", main="gamma")
plot(ts(post[,3]), xlab="Iteration", ylab="sigma", main="sigma")
acf(post[,1], main="")
acf(post[,2], main="")
acf(post[,3], main="")

# ESS
effectiveSize(post) #454.8073 465.0061 290.8386