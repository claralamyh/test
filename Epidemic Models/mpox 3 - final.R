library(deSolve)

###Rodents###
x0_r <- c(Sr=1074103, Er=1074103, Ir=1074103)
par_r <- c(LAMBDA=0.9, mu=0.002, nu=0.0001, phi=0.007, betarr=0.0057)

euler_r <- function(x0, par, T, deltat) {
  fct <- function(t, x0, par) {
    with(as.list(c(x0, par)), {
      lambdar <- betarr*Ir/(Sr+Er+Ir)
      dSr <- LAMBDA - (lambdar + mu)*Sr
      dEr <- lambdar*Sr -(phi + mu)*Er
      dIr <- phi*Er - (mu + nu)*Ir
      list(c(dSr, dEr, dIr))
    })
  }
  times <- seq(0, T, by=deltat)
  ode(y=x0, times=times, func=fct, parms=par)
}

out_r <- euler_r(x0_r, par_r, 70, 0.1)
plot(ts(out_r[,2:4], start=0, deltat=0.1))


###Human###
x0_h <- c(S=68530739, E=0, I=31412, Q=0, R=0)
par_h <- c(LAMBDA=8644, mu=0.025, nu=0.00008, phi=0.007, xi=0.00001, theta=0.029, gamma=0.0081, psi=0.025, delta=0.012, betarh=0.000009, betahh=0.03)

euler_h <- function(x0, par, T, deltat) {
  Ir <- function(t) {
    out_r[t/deltat+1, 4]
  }
  fct <- function(t, x0, par, input) {
    with(as.list(c(x0, par)), {
      Ir <- input(t)
      lambda <- (betarh*Ir + betahh*I) / (S+E+I+Q+R)
      dS <- LAMBDA + xi*R + theta*Q - (lambda + mu)*S
      dE <- lambda*S -(gamma + phi + mu)*E
      dI <- phi*E - (psi + mu + nu)*I
      dQ <- gamma*E - (theta + delta + mu + nu)*Q
      dR <- psi*I + delta*Q - (xi + mu)*R
      list(c(dS, dE, dI, dQ, dR))
    })
  }
  times <- seq(0, T, by=deltat)
  ode(y=x0, times=times, func=fct, parms=par, input=Ir)
}

out_h <- euler_h(x0_h, par_h, 70, 0.1)
plot(ts(out_h[,2:6], start=0, deltat=0.1))

changeI <- par_h["phi"]*(out_h[,3])
plot(ts(changeI, start=0, deltat=0.1)) 

# basic reproduction number
n <- par_h["betahh"]*par_h["LAMBDA"]*par_h["phi"]
d <- par_h["mu"]*(par_h["gamma"]+par_h["phi"]+par_h["mu"])*(par_h["psi"]+par_h["mu"]+par_h["nu"])
unname(R0 <- n/d)
