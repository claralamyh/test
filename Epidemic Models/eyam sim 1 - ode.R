library(deSolve)

###Eyam plague data###
st <- c(254,235,201,153,121,110,97,NA,83)
it <- c(7,14,22,29,20,8,8,NA,0)


###Eyam plague initial condition###
x0 <- c(254, 7)
par <- c(exp(-4.5), 1.5)


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

out <- euler(x0, par, 4, 0.005)


###Plot sir function###
plotsir <- function(time, sir) {
  matplot(x=time, y=sir,
          type = "l", lty = "solid",
          col = c("green", "red", "blue"),
          xlab = "Time (month)", ylab = "Population")
  lines(seq(0,4,0.5), st, type="p")
  lines(seq(0,4,0.5), it, type="p", pch=2)
} 

plotsir(out[,1], out[,2:4])


###Inference on R_0###
(R_0 <- unname(par[1]/par[2]*sum(x0)))
max <- which(out[,2]/sum(x0) < 1/R_0)[1]
(maxt <- out[max,1]) #2.285 
abline(v=maxt, lty="dashed")
