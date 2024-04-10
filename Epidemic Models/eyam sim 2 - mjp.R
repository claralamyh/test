library(devtools)
library(MASS)
library(issb)

###Eyam plague initial condition###
x0 <- c(S=254, I=7)
st <- c(254,235,201,153,121,110,97,NA,83)
it <- c(7,14,22,29,20,8,8,NA,0)
par <- c(beta=exp(-4.5), gamma=1.5)
dt <- 0.005


###Stochastic setup###
# hazard functions
h <- function(x, pars) {
  hazs <- rep(0, length(pars))
  hazs[1] <- pars[1]*x[1]*x[2]
  hazs[2] <- pars[2]*x[2]
  hazs
}

# stoichiometry matrix
smat <- matrix(c(-1, 0, 1, -1), 
               byrow = TRUE, nrow = 2,
               dimnames = list(c("S", "I")))

###Gillespie###
model <- create_model(smat, h, x0, par)
set.seed(1227)
g <- gillespie(model, maxtime=4)
R_sim <- rep(sum(x0),nrow(g)) - g[,1] - g[,2]
g.out <- cbind(g, R_sim)


###Poisson leap###
poissonleap <- function(x0, par, smat, T, deltat) {
  n <- T/deltat + 1
  mat <- matrix(0, nrow=n, ncol=2)
  colnames(mat) <- c("S", "I")
  mat[1,] <- x0
  for(t in 2:n) {
    deltaR1 <- rpois(1, h(mat[t-1,], par)[1]*deltat)
    deltaR2 <- rpois(1, h(mat[t-1,], par)[2]*deltat)
    deltaR <- c(deltaR1, deltaR2)
    mat[t,] <- mat[t-1,] + t(smat %*% deltaR)
  }
  R <- c(rep(sum(x0), nrow(mat)) - mat[,1] - mat[,2])
  out <- cbind(mat, R)
}

set.seed(3421)
p.out <- poissonleap(x0, par, smat, 4, dt)


###Plot sir function###
plotsir <- function(time, sir) {
  matplot(x=time, y=sir,
          type = "l", lty = "solid",
          col = c("green", "red", "blue"),
          xlab = "Time (month)", ylab = "Population")
  lines(seq(0,4,0.5), st, type="p")
  lines(seq(0,4,0.5), it, type="p", pch=2)
}


# single realisation
par(mfrow=c(1,2))
plotsir(g.out[,1], g.out[,2:4])
title(main="Gillespie algorithm")

plotsir(seq(0,4,dt), p.out[,1:3])
title(main="Poisson leap")


# multiple realisations
m <- 50

#gillespie
outarray <- array(0, c(4/dt+1, m, 2))
set.seed(1335)
for (i in 1:m) {
  g <- gillespie(model, maxtime=4, tstep=0.005)
  outarray[,i,1] <- g[,2]
  outarray[,i,2] <- g[,3]
}

S.mat <- outarray[,,1]
I.mat <- outarray[,,2]

smean <- apply(S.mat, 1, mean)
slq <- apply(S.mat, 1, quantile, 0.025)
suq <- apply(S.mat, 1, quantile, 0.975)

imean <- apply(I.mat, 1, mean)
ilq <- apply(I.mat, 1, quantile, 0.025)
iuq <- apply(I.mat, 1, quantile, 0.975)

plot(ts(S.mat[,1], start=0, deltat=0.005), ylim=c(40,260), col="grey",
     xlab = "Time (month)", ylab = "Susceptible", main = "Gillespie algorithm")
for(i in 2:50) {
  lines(ts(S.mat[,i], start=0, deltat=0.005), col="grey")
}
lines(ts(smean, start=0, deltat=0.005), lwd=1.3)
lines(ts(slq, start=0, deltat=0.005), col="blue", lwd=1.3)
lines(ts(suq, start=0, deltat=0.005), col="blue", lwd=1.3)
lines(seq(0,4,0.5)[-8], na.omit(st), col="red", lwd=1.3)

plot(ts(I.mat[,1], start=0, deltat=0.005), ylim=c(0,100), col="grey",
     xlab = "Time (month)", ylab = "Infectious", main = "Gillespie algorithm")
for(i in 2:50) {
  lines(ts(I.mat[,i], start=0, deltat=0.005), col="grey")
}
lines(ts(imean, start=0, deltat=0.005), lwd=1.3)
lines(ts(ilq, start=0, deltat=0.005), col="blue", lwd=1.3)
lines(ts(iuq, start=0, deltat=0.005), col="blue", lwd=1.3)
lines(seq(0,4,0.5)[-8], na.omit(it), col="red", lwd=1.3)

#poisson
set.seed(1156)
for (i in 1:m) {
  p <- poissonleap(x0, par, smat, 4, dt)
  outarray[,i,1] <- p[,1]
  outarray[,i,2] <- p[,2]
}

S.mat <- outarray[,,1]
I.mat <- outarray[,,2]

smean <- apply(S.mat, 1, mean)
slq <- apply(S.mat, 1, quantile, 0.025)
suq <- apply(S.mat, 1, quantile, 0.975)

imean <- apply(I.mat, 1, mean)
ilq <- apply(I.mat, 1, quantile, 0.025)
iuq <- apply(I.mat, 1, quantile, 0.975)

plot(ts(S.mat[,1], start=0, deltat=0.005), ylim=c(40,260), col="grey",
     xlab = "Time (month)", ylab = "Susceptible", main = "Poisson leap")
for(i in 2:50) {
  lines(ts(S.mat[,i], start=0, deltat=0.005), col="grey")
}
lines(ts(smean, start=0, deltat=0.005), lwd=1.3)
lines(ts(slq, start=0, deltat=0.005), col="blue", lwd=1.3)
lines(ts(suq, start=0, deltat=0.005), col="blue", lwd=1.3)
lines(seq(0,4,0.5)[-8], na.omit(st), col="red", lwd=1.3)

plot(ts(I.mat[,1], start=0, deltat=0.005), ylim=c(0,100), col="grey",
     xlab = "Time (month)", ylab = "Infectious", main = "Poisson leap")
for(i in 2:50) {
  lines(ts(I.mat[,i], start=0, deltat=0.005), col="grey")
}
lines(ts(imean, start=0, deltat=0.005), lwd=1.3)
lines(ts(ilq, start=0, deltat=0.005), col="blue", lwd=1.3)
lines(ts(iuq, start=0, deltat=0.005), col="blue", lwd=1.3)
lines(seq(0,4,0.5)[-8], na.omit(it), col="red", lwd=1.3)
