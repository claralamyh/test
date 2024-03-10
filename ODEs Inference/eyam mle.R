library(deSolve)

###Solve ODEs###
euler <- function(x0, par, T, deltat) {
  fct <- function(t, x0, par) {
    with(as.list(c(x0, par)), {
      dS <- -par[1]*x0[1]*x0[2]
      dI <- par[1]*x0[1]*x0[2] - par[2]*x0[2]
      list(c(dS, dI))
    })
  }
  times <- seq(0, T, by=deltat)
  out <- ode(y=x0, times=times, func=fct, parms=par)
  R <- c(rep(sum(x0), nrow(out)) - out[,2] - out[,3])
  mat <- cbind(out, R)
}


###Evaluate NEGATIVE log likelihood###
loglikeNeg <- function(par, yt, x0, T, deltat, i) {
  sim <- euler(x0, exp(par[1:2]), T, deltat)
  #i=2 gives St, i=3 gives It
  StIt <- sim[c(seq(1,nrow(sim),0.5/deltat)), i]
  val <- sum(dnorm(yt,StIt,exp(par[3]),log=TRUE), na.rm=TRUE)    
  -val #want mle, hence need -val
}


###Eyam plague data###
st <- c(254,235,201,153,121,110,97,NA,83)
it <- c(7,14,22,29,20,8,8,NA,0)
par <- log(c(exp(-4.5), 1.5, 3)) #log(beta,gamma,sigma)


###Inference###
# Susceptible
out.s <- optim(par=par, fn=loglikeNeg, y=st, x0=c(254,7), T=4, deltat=0.005, i=2)
(mle.s <- exp(out.s$par)) #mles: 0.01713254 2.71342790 3.38861275
out.s$value #minimised -l(\theta|y) value: 21.11476

# Infectious
out.i <- optim(par=par, fn=loglikeNeg, y=it, x0=c(254,7), T=4, deltat=0.005, i=3)
(mle.i <- exp(out.i$par)) #mles: 0.01920681 3.09469807 2.15539531
out.i$value #minimised -l(\theta|y) value: 17.49746

# Producing simulations
sim.par <- euler(c(254,7),c(exp(par[1:2])),4,0.005)
sim.mle.s <- euler(c(254,7),c(mle.s[1:2]),4,0.005)
sim.mle.i <- euler(c(254,7),c(mle.i[1:2]),4,0.005)

plotsir <- function(time, sir) {
  matplot(x=time, y=sir,
          type = "l", lty = "solid",
          col = c("green", "red", "blue"),
          xlab = "Time (month)", ylab = "Population")
  lines(seq(0,4,0.5), st, type="p")
  lines(seq(0,4,0.5), it, type="p", pch=2)
} 

par(mfrow=c(3,1))
plotsir(sim.par[,1], sim.par[,2:4])
title("Simulation using initial parameters")
plotsir(sim.mle.s[,1], sim.mle.s[,2:4])
title("Simulation using MLEs from Susceptible data")
plotsir(sim.mle.i[,1], sim.mle.i[,2:4])
title("Simulation using MLEs from Infectious data")
par(mfrow=c(1,1))

