library(deSolve)

#x0=c(Sr,Er,Ir)
x0_r <- c(Sr=1074103, Er=1074103, Ir=1074103)
#par=c(1-DELTA,2-mu,3-nu,4-phi,5-betarr)
par_r <- c(DELTA=0.9, mu=0.002, nu=0.0001, phi=0.007, betarr=0.0057)


euler_r <- function(x0, par, T, deltat) {
  fct <- function(t, x0, par) {
    with(as.list(c(x0, par)), {
      lambdar <- betarr*Ir/(Sr+Er+Ir)
      dSr <- DELTA - (lambdar + mu)*Sr
      dEr <- lambdar*Sr -(phi + mu)*Er
      dIr <- phi*Er - (mu + nu)*Ir
      list(c(dSr, dEr, dIr))
    })
  }
  times <- seq(0, T, by=deltat)
  ode(y=x0, times=times, func=fct, parms=par)
}

out_r <- euler_r(x0_r, par_r, 70, 0.1)
# plot(ts(out_r[,2:4], start=0, deltat=0.1))




