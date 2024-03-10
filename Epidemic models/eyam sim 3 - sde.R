library(MASS)

# Euler-Maruyama - regular time grid
euler.m <- function(x0, par, T, deltat) {
  n <- T/deltat + 1
  mat <- matrix(0, nrow=n, ncol=2)
  colnames(mat) <- c("S", "I")
  mat[1,] <- x0
  for(i in 2:n) {
    h1 <- par[1]*mat[i-1,1]*mat[i-1,2]
    h2 <- par[2]*mat[i-1,2]
    infmean <- c(-h1, h1-h2)
    infsd <- t(-chol(matrix(c(h1,-h1,-h1,h1+h2), ncol=2, byrow=T)))
    mat[i,] <- (mat[i-1,]+infmean*deltat) + infsd%*%rnorm(2,0,sqrt(deltat))
    for(j in 1:2) {
      if(mat[i,j] < 0.001) {
        mat[i,j] <- 0.001
      }
    }
  }
  R <- c(rep(sum(x0),nrow(mat)) - mat[,1] - mat[,2])
  out <- cbind(mat, R)
}


###Eyam plague###
# data and initial condition
x0 <- c(254, 7)
par <- c(exp(-4.5), 1.5)
st <- c(254,235,201,153,121,110,97,NA,83)
it <- c(7,14,22,29,20,8,8,NA,0)

# single realisation
set.seed(1529)
out <- euler.m(x0, par, 4, 0.005)
plotsir(seq(0,4,0.005), out[,1:3])
(R_0 <- unname(par[1]/par[2]*sum(x0))) #1.93
max <- which(out[,1]/sum(x0) < 1/R_0)[1]
(maxt <- seq(0,4,0.005)[max]) #1.95
abline(v=maxt, lty="dashed")

# multiple realisations
m <- 50
outarray <- array(0, c(nrow(out), m, 2))

set.seed(1538)
for (i in 1:m) {
  sim <- euler.m(x0, par, 4, 0.005)
  outarray[,i,1] <- sim[,1]
  outarray[,i,2] <- sim[,2]
}

S.mat <- outarray[,,1]
I.mat <- outarray[,,2]

smean <- apply(S.mat, 1, mean)
slq <- apply(S.mat, 1, quantile, 0.025)
suq <- apply(S.mat, 1, quantile, 0.975)

imean <- apply(I.mat, 1, mean)
ilq <- apply(I.mat, 1, quantile, 0.025)
iuq <- apply(I.mat, 1, quantile, 0.975)

par(mfrow=c(1,2))
plot(ts(S.mat[,1], start=0, deltat=0.005), ylim=c(40,260),
     xlab = "Time (month)", ylab = "Susceptible", col="grey")
for(i in 2:50) {
  lines(ts(S.mat[,i], start=0, deltat=0.005), col="grey")
}
lines(ts(smean, start=0, deltat=0.005), lwd=1.3)
lines(ts(slq, start=0, deltat=0.005), col="blue", lwd=1.3)
lines(ts(suq, start=0, deltat=0.005), col="blue", lwd=1.3)
lines(seq(0,4,0.5)[-8], na.omit(st), col="red", lwd=1.3)

plot(ts(I.mat[,1], start=0, deltat=0.005), ylim=c(0,80),
     xlab = "Time (month)", ylab = "Infectious", col="grey")
for(i in 2:50) {
  lines(ts(I.mat[,i], start=0, deltat=0.005), col="grey")
}
lines(ts(imean, start=0, deltat=0.005), lwd=1.3)
lines(ts(ilq, start=0, deltat=0.005), col="blue", lwd=1.3)
lines(ts(iuq, start=0, deltat=0.005), col="blue", lwd=1.3)
lines(seq(0,4,0.5)[-8], na.omit(it), col="red", lwd=1.3)





###Plot sir function###
plotsir <- function(time, sir) {
  matplot(x=time, y=sir,
          type = "l", lty = "solid",
          col = c("green", "red", "blue"),
          xlab = "Time (month)", ylab = "Population")
  lines(seq(0,4,0.5), st, type="p")
  lines(seq(0,4,0.5), it, type="p", pch=2)
}
