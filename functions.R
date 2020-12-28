## Computes the likelihood for one given set of observations for (S1,S2,W).
## Since the likelihood involves nested integrals, each integral is approximated
## by a Riemann sum, the quality of the approximation can be controlled with N.
## Input:
## x1, x2        - Observed values for S1 and S2 (one value each)
## w             - Observed value for W (can be equal to x1)
## theta1_polar, - Parameters to test given in polar coordinates, theta1_polar 
## theta2_polar    is used for the incubation distribution and theta2_polar for
##                 the generation time distribution.
## N             - Length of the Riemann sum in the integral approximation. De-
##                 fault value is 500.
## Output: Value of the likelihood
inc_likelihood <- function(x1,x2,w,theta1_polar,theta2_polar,N=500) {
  theta1 <- polar2rect(1,theta1_polar)
  theta2 <- polar2rect(1,theta2_polar)
  
  ## Approximating sequence for y
  y     <- seq(from=0,to=x2,length.out=N)
  Delta <- y[2]-y[1]
  
  ## Approximating sequence for x1 and x2
  ax1 <- x1-(0:(N-1))*Delta
  ax2 <- x2-(0:(N-1))*Delta
  
  ## Compute Laguerre polynomials on grid
  fIx1 <- laguerre_density_pos(ax1,theta1)*(ax1>=x1-min(c(x1,w)))
  fIx2 <- laguerre_density_pos(ax2,theta1)
  fG   <- laguerre_density_pos(y  ,theta2)
  
  ## Produce Matrix for inner integral and the inner integral
  A <- toeplitz(fIx1)
  A[lower.tri(A,diag=FALSE)] <- 0
  A <- t(t(A)*fIx2)
  
  I <- rowSums((A[,1:(N-1)]+A[,2:N])/2)*Delta
  
  ## Compute Integral
  oI <- fG*I
  log_lik <- log(sum((oI[1:(N-1)]+oI[2:N])/2)*Delta)
  
  return(log_lik)
}

## Computes f_{theta}(z), i.e.,  the Laguerre density on the positives as speci-
## fied in the paper.
## Input:
##  z     - Vector of points at which the density shall be evaluated
##  theta - Parameter vector provided in Cartesian coordinates (must have norm 1)
## Output: A vector of the same length as z which contains the corresponding
##         values of the density.
laguerre_density_pos <- function(z,theta) {
  f <- rep(0,length(z))
  
  ## Prepare Data matrices
  m <- length(theta)-1
  ind <- which(z>=0)
  n <- length(ind)
  z <- z[ind]
  
  ## Compute the matrix of binomial coefficients
  B <- matrix(0,ncol=m+1,nrow=m+1)
  for(k in 1:(m+1)) {
    B[k,1:k] <- choose(k-1,0:(k-1))
  }
  
  ## Compute the powers of z
  Zp <- matrix(1,ncol=m+1,nrow=n)
  if(m+1>=2) {
    for(i in 2:(m+1)) {
      Zp[,i] <- Zp[,i-1]*(-z)/(i-1)
    }
  }
  
  ## Compute Laguerre Polynomials
  Lp <- Zp%*%t(B)
  
  ## Compute the density
  f[ind] <- exp(-z)*(Lp%*%theta)^2
  
  return(f)
}

## Compute the Laguerre Polynom of given order.
## Input:
##  z - points at which the Laguerre polynomial shall be evaluated (must be non-
##      negative).
##  k - Order of the Laguerre-polynomial to compute (natural number, k>=0)
## Output: Vector of the same length as z containing the respective evaluations.
laguerre_polynomial <- function(z,k) {
  ## Compute binomial coefficients
  B <- choose(k,0:k)

  ## Compute the powers of z
  n <- length(z)
  Zp <- matrix(1,ncol=k+1,nrow=n)
  if(k>=1) {
    for(i in 2:(k+1)) {
      Zp[,i] <- Zp[,i-1]*(-z)/(i-1)
    }
  }
  
  ## Compute Laguerre Polynomials
  return(Zp%*%B)
}


## Computes the negative log-likelihood for a given parameter and observation
## set. This function needs to be minimized in order to obtain the sieve MLE.
## Input:
## theta   - Vector of length m1+m2 with the polar coordinates of the two para-
##           meters theta1 (entries 1:m1, for the incubation period) and theta2
##           (entries (m1+1):(m1+m2), for the generation time).
## x1,x2   - Two vectors of the same length, every entry corresponds to one ob-
##           servation. The entries in x1 are the observations for S1 and the
##           entries in x2 are the observations for S2.
## w       - Vector of the same length as x1 and x2. Each entry corresponds to
##           an observation of min(S1,W).
## m1,m2   - Two integers (both at least 1) which specify the size of the fitted
##           model. See also the description of theta.
## N       - Same as for inc_likelihood
## cluster - Either FALSE (the default), then normal serial computations are
##           carried out, or an cluster identifier (returned from makeCluster)
##           to do parallel computations by using the foreach framework (requi-
##           res packages doParallel and foreach).
## Output: The un-normalized (i.e. not divided by the number of observations),
##         negative log-likelihood corresponding to the given observations
##         and models.
whole_inc_likelihood <- function(theta,x1,x2,w,m1,m2,N=500,cluster=FALSE) {
  n <- length(x1)
  
  theta1_polar <- theta[1:m1]
  theta2_polar <- theta[(m1+1):(m1+m2)]
  
  if(isFALSE(cluster)) {
    ## DO Serial Computation
    out <- 0
    for(i in 1:n) {
      out <- out+inc_likelihood(x1[i],x2[i],w[i],theta1_polar,theta2_polar,N)
    }
  } else {
    clusterCall(cluster,source,'functions.R')
    out <- foreach(i=1:n,.combine=cbind,.packages='SphericalCubature') %dopar% {
      inc_likelihood(x1[i],x2[i],w[i],theta1_polar,theta2_polar,N)
    }
    out <- sum(out)
  }
  
  return(-out)
}

## Computes the best Laguerre approximation to a density for a given degree.
## Input
##  fdens - Function which shall be approximated. The first argument of fdens
##          must be the vector of evaluation points and the output of fdens must
##          be a vector of the same langth containing the respective evaluations.
##  m     - Order of the Laguerre Approximation
##  ...   - Additional arguments which are passed to fdens.
## Output: Vector of length m+1 which contains the normalized Cartesian coordi-
##         nates of the best Laguerre approximation of fdens of order m.
laguerre_approx <- function(fdens,m,...) {
  ## Compute the un-scaled coefficients
  theta <- rep(0,m+1)
  for(k in 1:(m+1)) {
    theta[k] <- cubintegrate(function(x,...) {return(exp(-x/2)*sqrt(fdens(x,...))*laguerre_polynomial(x,k-1))},lower=0,upper=Inf,...)$integral
  }
  
  return(theta/sqrt(sum(theta^2)))
}