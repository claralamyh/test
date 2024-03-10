library(MASS)
library(mvtnorm)

###Data###
data <- cbind(time = seq(0,4,0.5)[-8], 
              st = c(254,235,201,153,121,110,97,83),
              it = c(7,14,22,29,20,8,8,0))

par <- log(c(exp(-4.5), 1.5))


###SIR drift and variance functions###
alphaSIR <- function(x, theta) {
  c(-theta[1]*x[1]*x[2],
    theta[1]*x[1]*x[2]-theta[2]*x[2])
}

betaSIR <- function(x, theta) {
  cov <- -theta[1]*x[1]*x[2]
  var1 <- theta[1]*x[1]*x[2]
  var2 <- theta[1]*x[1]*x[2] + theta[2]*x[2]
  matrix(c(var1,cov,cov,var2), ncol=2, byrow=T)
}


###E-M approximation - irregular time grid###
emsim <- function(x0, time, par, afun, bfun) {
  n <- length(time)
  mat <- matrix(0, nrow=n, ncol=2)
  colnames(mat) <- c("S", "I")
  mat[1,] <- x0
  for(i in 2:n) {
    deltat <- time[i] - time[i-1]
    mat[i,] <- (mat[i-1,]+afun(mat[i-1,],par)*deltat) +
      t(-chol(bfun(mat[i-1,],par)))%*%rnorm(2,0,sqrt(deltat))
    for(j in 1:2) {
      if(mat[i,j] < 0.001) {
        mat[i,j] <- 0.001
      }
    }
  }
  out <- cbind(time, mat)
}


###One-step euler transition density###
pe <- function(to, from, par, deltat, afun, bfun) {
  dmvnorm(x = to,
          mean = from + afun(from,par)*deltat,
          sigma = bfun(from,par)*deltat, log = TRUE)
}


###Log likelihood/ joint density of all x###
loglike <- function(xdata, time, psi, afun, bfun) {
  val <- rep(0, ncol(xdata)-1)
  for(i in 2:ncol(xdata)) {
    deltat <- time[i] - time[i-1]
    val[i-1] <- pe(xdata[1:2,i], xdata[1:2,i-1], exp(psi), deltat, afun, bfun)
  }
  sum(val)
}


###Log posterior###
lpost <- function(xdata, time, psi, afun, bfun) {
  lprior <- sum(dnorm(psi, log=TRUE))
  llike <- loglike(xdata, time, psi, afun, bfun)
  lprior + llike
}


###Latent process segment update###
latent <- function(i, m, x.all, time, par, afun, bfun, count2) {
  t <- time[(1+(i-1)*m):(i*m)]
  deltat <- t[2]-t[1]
  can.x <- matrix(emsim(x.all[1:2,1+(i-1)*m], t, par, afun, bfun)[2:m,], nrow=m-1)
  x.unobs <- t(x.all[1:2,(2+(i-1)*m):(i*m)])
  laprob <- pe(x.all[1:2,1+i*m], can.x[m-1,2:3], par, deltat, afun, bfun) -
    pe(x.all[1:2,1+i*m], x.unobs[m-1,], par, deltat, afun, bfun)
  if (log(runif(1)) < laprob) {
    x.unobs <- matrix(can.x[,2:3], ncol=2) 
    count2[i] <- count2[i] + 1
  } 
  return(list(x.unobs, count2))
}


###Gibbs###
gibbs <- function(N, m, Vtune, data, par, afun, bfun) {
  n <- nrow(data)-1
  mat.psi <- matrix(0, nrow=N, ncol=2)   #store psi = log theta
  mat.s <- matrix(0, nrow=N, ncol=n*m+1) #store all (unobs+obs) st
  mat.i <- matrix(0, nrow=N, ncol=n*m+1) #store all (unobs+obs) it
  
  #initialise parameters
  psi <- par
  mat.psi[1,] <- psi
  
  #initialise the process
  for(j in 1:N) { #observed x
    mat.s[j,seq(1,n*m+1,m)] <- data[,2]
    mat.i[j,seq(1,n*m+1,m)] <- data[,3]
  }
  
  tau <<- as.numeric(rep(data[n+1,1], n*m+1)) #new time step
  
  for(i in 1:n) { #unobserved x
    t <- seq(data[i,1], data[i+1,1], length.out=m+1)[-(m+1)]
    x.unobs <- emsim(data[i,2:3], t, exp(par), afun, bfun)
    tau[(1+(i-1)*m):(i*m)] <<- t
    mat.s[1,(2+(i-1)*m):(m+(i-1)*m)] <- x.unobs[2:m,2]
    mat.i[1,(2+(i-1)*m):(m+(i-1)*m)] <- x.unobs[2:m,3]
  }
  
  #count acceptance
  count1 <- 0 #theta 
  count2 <- rep(0,n) #unobs
  
  for(j in 2:N) {
    # update parameters
    x.all <- matrix(c(mat.s[j-1,], mat.i[j-1,]), nrow=2, byrow=TRUE)
    can.psi <- psi + mvrnorm(1, rep(0,2), Vtune)
    laprob <- lpost(x.all, tau, can.psi, afun, bfun) -
      lpost(x.all, tau, psi, afun, bfun)
    if (log(runif(1)) < laprob) {
      psi <- can.psi 
      count1 <- count1+1
    }
    mat.psi[j,] <- psi
    
    # update unobserved process
    for(i in 1:n) {
      latent <- latent(i, m, x.all, tau, exp(mat.psi[j,]), afun, bfun, count2)
      x.unobs <- latent[[1]]
      mat.s[j,(2+(i-1)*m):(m+(i-1)*m)] <- x.unobs[,1]
      mat.i[j,(2+(i-1)*m):(m+(i-1)*m)] <- x.unobs[,2]
      count2 <- latent[[2]]
    }
  }
  
  print(count1/(N-1))
  return(list(exp(mat.psi), mat.s, mat.i, count2/rep(N-1,n), tau))
}


