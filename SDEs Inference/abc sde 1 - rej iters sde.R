###E-M approximation - regular time grid###
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
  mat
}


###ABC-Rejection###
abcrej.iters <- function(N, epsilon, x0, par, T, deltat, yt) {
  mat <- matrix(0, nrow=N, ncol=2) #store values
  mat[1,] <- par
  for (i in 2:N) {
    rho <- epsilon + 1
    while (rho > epsilon) {
      can <- c(par[1]+rnorm(1), par[2]+rnorm(1)) #proposed candidate
      simdata <- euler.m(x0, exp(can), T, deltat)
      zdata <- simdata[c(seq(1,nrow(simdata),0.5/deltat)),2]
      rho <- sqrt(sum(na.omit(zdata-yt)^2)) #Euclidean distance
      if (rho <= epsilon) {
        mat[i,] <- can
        break
      }
    }
  }
  mat
}


###Measuring the CPU time used###
wrapper <- function(code) {
  start <- Sys.time()
  code
  end <- Sys.time()
  end-start
}


###Executing ABC-Rejection###
x0<-c(254,7); par<-log(c(0.019, 3.12)); T<-4; deltat<-0.005 
st <- c(254,235,201,153,121,110,97,NA,83)
it <- c(7,14,22,29,20,8,8,NA,0)

set.seed(1126)
wrapper(
  abc <- abcrej.iters(N=140, epsilon=25, x0, par, T, deltat, it))
# Time difference of 1.327933 mins


###Posterior samples - abc rejection###
out <- exp(abc)

# summaries
mean(out[,1]) #0.01758186
sd(out[,1]) #0.007139055
quantile(out[,1],c(0.025,0.975)) #0.005468097 0.034007212

mean(out[,2]) #3.00196
sd(out[,2]) #1.014942
quantile(out[,2],c(0.025,0.975)) #1.204525 5.403973

