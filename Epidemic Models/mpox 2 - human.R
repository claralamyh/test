library(deSolve)

#x0=c(S,E,I,Q,R)
x0_h <- c(S=68530739, E=0, I=31412, Q=0, R=0)
#par=c(1-DELTA,2-xi,3-theta,4-mu,5-nu,6-phi,7-psi,8-gamma,9-delta,10-betahh,11-betarh)
par_h <- c(DELTA=8644, xi=0.00001, theta=0.029, mu=0.05, nu=0.00008, phi=0.007, psi=0.056, gamma=0.0081, delta=0.012, betahh=0.00008, betarh=0.000009)


euler_h <- function(x0, par, T, deltat) {
  Ir <- function(t) {
    out_r[t/deltat+1, 4]
  }
  fct <- function(t, x0, par, input) {
    with(as.list(c(x0, par)), {
      Ir <- input(t)
      lambda <- (betarh*Ir + betahh*I) / (S+E+I+Q+R)
      dS <- DELTA + xi*R + theta*Q - (lambda + mu)*S
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


# euler_h <- function(x0, par, Ir, T, deltat) {
#   S <- x0[1]; E <- x0[2]; I <- x0[3]; Q <- x0[4]; R <- x0[5]
#   n <- T/deltat + 1
#   mat <- matrix(0, nrow=n, ncol=6)
#   colnames(mat) <- c("Time", "S", "E", "I", "Q", "R")
#   mat[1,] <- c(0, S, E, I, Q, R)
#   for(i in 2:n) {
#     lambda <- (par[11]*Ir[i-1]+par[10]*I) / sum(mat[i-1,2:6])
#     S <- S + (par[1]+par[2]*R+par[3]*Q-lambda*S-par[4]*S)*deltat
#     E <- E + (lambda*S-(par[8]+par[6]+par[4])*E)*deltat
#     I <- I + (par[6]*E-(par[7]+par[4]+par[5])*I)*deltat
#     Q <- Q + (par[8]*E-(par[3]+par[9]+par[4]+par[5])*Q)*deltat
#     R <- R + (par[7]*I+par[9]*Q-par[2]*R-par[4]*R)*deltat
#     mat[i,] <- c((i-1)*deltat, S, E, I, Q, R)
#   }
#   mat
# }
# 
# out_h <- euler_h(x0_h, par_h, out_r[,4], 70, 0.1)

changeI <- 0.007*out_h[,3]
plot(ts(changeI, start=0, deltat=0.1))
