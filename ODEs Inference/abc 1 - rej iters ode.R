###Solve ODEs###
euler <- function(x0, par, T, deltat) {
  n <- T/deltat + 1
  mat <- matrix(0, nrow=n, ncol=2)
  colnames(mat) <- c("S", "I")
  mat[1,] <- c(x0[1], x0[2])
  for(i in 2:n) {
    St <- mat[i-1,1]
    It <- mat[i-1,2]
    mat[i,1] <- St - par[1]*St*It*deltat
    mat[i,2] <- It + par[1]*St*It*deltat - par[2]*It*deltat
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
      simdata <- euler(x0, exp(can[1:2]), T, deltat)
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
wrapper <- function(code){
  start <- Sys.time()
  code
  end <- Sys.time()
  end-start
}


###Executing ABC-Rejection###
x0<-c(254,7); par<-log(c(0.019, 3.12)); T<-4; deltat<-0.005 
st <- c(254,235,201,153,121,110,97,NA,83)
it <- c(7,14,22,29,20,8,8,NA,0)

set.seed(2002)
wrapper(
abc <- abcrej.iters(N=1000, epsilon=20, x0, par, T, deltat, it))
# Time difference of 1.79585 mins


###Posterior samples###
out <- exp(abc)


# summaries
mean(out[,1]) #0.01920807
sd(out[,1]) #0.003568882
quantile(out[,1],c(0.025,0.975)) #0.01335317 0.02596277

mean(out[,2]) #3.174201
sd(out[,2]) #0.5541015
quantile(out[,2],c(0.025,0.975)) #2.330247 4.293062 

