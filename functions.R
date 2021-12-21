## Fits Laguerre approximating densities of given degrees to a given data-set.
## Input:
##  S1           - Symptom onset times of primary cases
##  S2           - Symptom onset times of secondary cases (must have the same
##                 length as S1)
##  w            - End of exposure windows (must have the same length as S2)
##  rC           - Exponential growth rates to be applied for corresponding
##                 primary case infection (must have the same length as S1)
##  x0           - Starting value of the first global optimisation, if NULL
##                 (default) a random starting value is chosen.
##  m1           - Degree of approximation of the incubation time.
##  m2           - Degree of approximation of the generation time.
##  glob_repeats - Number of global optimisation tasks to be carried out from
##                 random starting points. The first optimisation is started
##                 from x0 if specified. glob_repeats must be specified and it
##                 must be a natural number larger than or equal to 1.
##  xtol_rel     - Value of xtol_rel for the local optimisation (default is
##                 0.0001), see nloptr for details. The global optimisation is
##                 carried out with 1000*xtol_rel.
##  maxeval      - Value of maxeval (see nloptr for details) for both, the local
##                 and global optimisation. Default is 10000.
##  int_prec     - Length of Riemann sum for approximating the integrals in the
##                 likelihood (default is 500).
##  print        - Logical, if TRUE (default) status messages are printed.
## Output: A list of four elements
##  theta1   - Optimal parameter values for the incubation time.
##  theta2   - Optimal parameter values for the generation time.
##  glob_opt - Output of nloptr of that global optimisation which yielded the
##             best result.
##  loc_opt  - Output of nloptr of the local optimisation.
fit_laguerre <- function(S1,S2,w,rC,x0=NULL,m1,m2,glob_repeats,xtol_rel=0.0001,maxeval=10000,int_prec=500,print=TRUE) {
  ## Do estimation
  if(print) {
    print_level=0
  } else {
    print_level=0
  }
  
  ## Set Random Starting Point if necesseary
  if(is.null(x0)==TRUE) {
    if(print) {
      cat("Set random starting Value\n")
    }
    x0 <- runif(m1+m2,min=0,max=pi)
  }
  
  ## Optimize
  currmin <- Inf
  for(k in 1:glob_repeats) {
    if(print) {
      cat("Do Global Optimisation ",k,"/",glob_repeats,"\n")
    }
    opts  <- list(algorithm="NLOPT_GN_CRS2_LM",print_level=print_level,xtol_rel=1000*xtol_rel,maxeval=maxeval)
    glob_opt_temp <- nloptr(x0=x0,eval_f=whole_inc_likelihood,lb=rep(0,m1+m2),ub=rep(pi,m1+m2),opts=opts,x1=S1,x2=S2,w=w,rC=rC,m1=m1,m2=m2,N=int_prec)
    x0 <- runif(m1+m2,min=0,max=pi)
    
    if(glob_opt_temp$objective<currmin) {
      if(print) {
        cat("Update to new minimum",glob_opt_temp$objective,"\n")
      }
      glob_opt <- glob_opt_temp
      currmin <- glob_opt_temp$objective
    }
  }
  
  
  ## Refine By local Optimisation
  if(print) {
    cat("Do Local Optimisation\n")
  }
  opts  <- list(algorithm="NLOPT_LN_SBPLX",print_level=print_level,xtol_rel=xtol_rel,maxeval=maxeval)
  loc_opt <- nloptr(x0=glob_opt$solution,eval_f=whole_inc_likelihood,lb=rep(0,m1+m2),ub=rep(pi,m1+m2),opts=opts,x1=S1,x2=S2,w=w,rC=rC,m1=m1,m2=m2,N=int_prec)
  
  ## Convert solution to rectangular coordinates
  theta1_est <- polar2rect(1,loc_opt$solution[1:m1])
  theta2_est <- polar2rect(1,loc_opt$solution[(m1+1):(m1+m2)])
  
  return(list(theta1=theta1_est,theta2=theta2_est,glob_opt=glob_opt,loc_opt=loc_opt))
}

## This function computes the optimal likelihood for a range of degrees. The
## output is not guaranteed to be monotonic due to numerical differences. This
## can be achieved by calling force_monotonicity in a second stage.
## Input:
##  S1, S2, w, RC - Input data, see the description in fit_laguerre for details.
##  m1_max        - Maximal degree for the approximation of the incubation time
##  m2_max        - Maximal degree for the approximation of the generation time 
##  print         - Logical, if TRUE (default) status information is printed
##  glob_repeats  - Same as for fit_Laguerre
##  xtol_rel      - Same as for fit_Laguerre
##  maxeval       - Same as for fit_Laguerre
##  int_prec      - Same as for fit_Laguerre
##  cluster       - If FALSE (default) computations are carried out serially. If
##                  this is a cluster object as returned by makeCluster, the
##                  computations for different approximating degrees are carried
##                  in parallel.
## Output: List of six elements
##  m1_max       - The provided value of m1_max
##  m2_max       - The provided value of m2_max
##  L_matrix     - Matrix of optimal negative log-likelihood values.
##                 L_matrix[m1,m2] contains the optimal value corresponding to a
##                 degree of m1 in the incubation time and m2 in the generation
##                 time.
##  fit_results  - A list of the return values from fit_laguerre for each
##                 degree. The i-th element of fit_results contains the results
##                 belonging to the degrees provided in dim_list[i,] (see next
##                 line).
##  dim_list     - A matrix with two columns: dim_list[i,]=c(m1,m2) and
##                 fit_results[[i]] contains the fir for degrees (m1,m2).
##  index_matrix - m1_max x m2_max matrix: index_matrix[m1,m2]=1 if and only if
##                 dim_list[i,]
model_selection_laguerre <- function(S1,S2,w,rC,m1_max,m2_max,print=TRUE,glob_repeats=5,xtol_rel=0.0001,maxeval=10000,int_prec=500,cluster=FALSE) {
  ## Write the matrix with dimensions to test
  test_dims <- matrix(0,ncol=2,nrow=m1_max*m2_max)
  i <- 1
  for(m1 in 1:m1_max) {
    for(m2 in 1:m2_max) {
      test_dims[i,] <- c(m1,m2)
      i <- i+1
    }
  }
  
  ## Compute the best estimate for each model level to test (possibly in parallel)
  if(isFALSE(cluster)==TRUE) {
    out <- vector(mode="list",length=m1_max*m2_max)
    for(i in 1:(m1_max*m2_max)) {
      set.seed(i)
      m1 <- test_dims[i,1]
      m2 <- test_dims[i,2]
      cat("Consider m1=",m1,"/",m1_max,", m2=",m2,"/",m2_max,".\n")
      
      ## Do the actual model fit
      Lfit <- fit_laguerre(S1,S2,w,rC,x0=NULL,m1,m2,glob_repeats,xtol_rel=xtol_rel,maxeval=maxeval,int_prec=int_prec,print=print)
      
      out[[i]] <- list(Lfit=Lfit,m1=m1,m2=m2)
    }
  } else {
    ## Do the same thing but in parallel
    clusterCall(cluster,source,'functions.R')
    clusterCall(cluster,dyn.load,"c_function.dll")
    
    out <- foreach(i=1:(m1_max*m2_max),.packages=c('SphericalCubature','nloptr')) %dopar% {
      set.seed(i)
      m1 <- test_dims[i,1]
      m2 <- test_dims[i,2]
      cat("Consider m1=",m1,"/",m1_max,", m2=",m2,"/",m2_max,".\n")
      
      ## Do the actual model fit
      Lfit <- fit_laguerre(S1,S2,w,rC,x0=NULL,m1,m2,glob_repeats,xtol_rel=xtol_rel,maxeval=maxeval,int_prec=int_prec,print=print)

      list(Lfit=Lfit,m1=m1,m2=m2)
    }
  }
  
  ## Create outputs
  index_matrix <- matrix(0,nrow=m1_max,ncol=m2_max)
      L_matrix <- matrix(0,nrow=m1_max,ncol=m2_max)
  fit_results <- vector(mode="list",length=m1_max*m2_max)
  dim_list <- matrix(0,ncol=2,nrow=m1_max*m2_max)
  for(i in 1:(m1_max*m2_max)) {
    m1 <- out[[i]]$m1
    m2 <- out[[i]]$m2
    
    index_matrix[m1,m2] <- i
    dim_list[i,] <- c(m1,m2)
    L_matrix[m1,m2] <- out[[i]]$Lfit$loc_opt$objective
    fit_results[[i]] <- out[[i]]$Lfit
  }

  return(list(m1_max=m1_max,m2_max=m2_max,L_matrix=L_matrix,dim_list=dim_list,fit_results=fit_results,index_matrix=index_matrix))
}

## This function can be called subsequently to model_selection_laguerre in order
## to guarantee that larger models yield a larger likelihood.
## Input:
##  MSL        - Output of a previous call of model_selection_laguerre.
##  S1,S2,w,rC - Dataset, provided in the same way as for fit_laguerre. This
##               dataset must be the same as the one with which
##               model_selection_laguerre has been called before.
##  xtol_rel   - Same as for fit_laguerre
##  maxeval    - Same as for fit_laguerre
##  int_prec   - Same as for fit_laguerre
##  print      - Logical, if TRUE (default) status messages will be printed.
## Output: A list of the exact same format as for model_selection_laguerre but
##         now the likelihoods in the matrix L_matrix are guaranteed to be
##         decreasing.
force_monotonicity <- function(MSL,S1,S2,w,rC,xtol_rel=0.0001,maxeval=10000,int_prec=500,print=TRUE) {
  ## Control for monotonicity
  for(m1 in 1:MSL$m1_max) {
    for(m2 in 1:MSL$m2_max) {
      rerun <- 0
      m1obj <- Inf
      m2obj <- Inf
      m1start <- 0
      m2start <- 0

      if(m1>=2) {
        if(MSL$L_matrix[m1,m2]>=MSL$L_matrix[m1-1,m2]) {
          rerun <- 1
          i <- MSL$index_matrix[m1-1,m2]

          m1start <- c(my_rect2polar(c(MSL$fit_results[[i]]$theta1,0)),my_rect2polar(MSL$fit_results[[i]]$theta2))
          m1obj   <- MSL$fit_results[[i]]$loc_opt$objective
        }
      }
      if(m2>=2) {
        if(MSL$L_matrix[m1,m2]>=MSL$L_matrix[m1,m2-1]) {
          rerun <- 1
          i <- MSL$index_matrix[m1,m2-1]

          m2start <- c(my_rect2polar(MSL$fit_results[[i]]$theta1),my_rect2polar(c(MSL$fit_results[[i]]$theta2,0)))
          m2obj   <- MSL$fit_results[[i]]$loc_opt$objective
        }
      }
      if(rerun==1) {
        if(print) {
          cat("Do a re-run for m1=",m1,"/",MSL$m1_max," and m2=",m2,"/",MSL$m2_max,"\n")
        }
	  ## Do a rerun
        ## Find best starting value
        if(m1obj<m2obj) {
          x0 <- m1start
        } else {
          x0 <- m2start
        }

        ## Do the re-run
        rfit <- fit_laguerre(S1,S2,w,rC,x0=x0,m1,m2,1,xtol_rel=xtol_rel,maxeval=maxeval,int_prec=int_prec,print=print)

        ## Save the output
        i <- MSL$index_matrix[m1,m2]
        MSL$L_matrix[m1,m2] <- rfit$loc_opt$objective
        MSL$fit_results[[i]] <- rfit
      }
    }
  }

  return(MSL)
}

## Computes f_{theta}(z), i.e.,  the Laguerre density on the positives as speci-
## fied in the paper.
## Input:
##  z     - Vector of points at which the density shall be evaluated
##  theta - Parameter vector provided in Cartesian coordinates (must have norm
##          1)
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

## Computes \int_0^zf_{theta}(x)dx, i.e.,  the distribution function of the cor-
## responding Laguerre density on the positives as specified in the paper.
## Input:
##  z     - Vector of points at which the distribution shall be evaluated
##  theta - Parameter vector provided in Cartesian coordinates (must have norm
##           1)
## Output: A vector of the same length as z which contains the corresponding
##         values of the density.
laguerre_distr_pos <- function(z,theta) {
  m <- length(theta)-1
  n <- length(z)
  
  ## Compute integrals \int_0^z x^ke^(-x)dx for k=0,...,2m
  pint <- matrix(0,ncol=2*m+1,nrow=n)
  pint[,1] <- 1-exp(-z)
  for(k in 2:(2*m+1)) {
    pint[,k] <- -z^(k-1)*exp(-z)+(k-1)*pint[,k-1]
  }
  
  ## Compute distribution function
  out <- rep(0,n)
  for(k1 in 0:m) {
    for(k2 in 0:m) {
      for(i1 in 0:k1) {
        for(i2 in 0:k2) {
          out <- out+theta[k1+1]*theta[k2+1]*choose(k1,i1)*choose(k2,i2)*(-1)^(i1+i2)/factorial(i1)/factorial(i2)*pint[,i1+i2+1]
        }
      }
    }
  }
  
  return(out)
}

## Computes the basic reproduction number belonging to a certain generation time
## density specified through a Laguerre density on the positives as specified in
## the paper.
## Input:
##  theta - Parameter vector provided in Cartesian coordinates (must have norm
##          1)
##  r     - Exponential growth parameter of the pandemic.
## Output: The corresponding basic reproduction number.
laguerre_R0 <- function(theta,r) {
  m <- length(theta)-1

  ## Compute integrals \int_0^infinity x^ke^(-(1+r)x)dx for k=0,...,2m
  pint <- rep(0,2*m+1)
  pint[1] <- 1/(1+r)
  for(k in 2:(2*m+1)) {
    pint[k] <- (k-1)/(1+r)*pint[k-1]
  }
  
  ## Compute FR
  FR <- 0
  for(k1 in 0:m) {
    for(k2 in 0:m) {
      for(i1 in 0:k1) {
        for(i2 in 0:k2) {
          FR <- FR+theta[k1+1]*theta[k2+1]*choose(k1,i1)*choose(k2,i2)*(-1)^(i1+i2)/factorial(i1)/factorial(i2)*pint[i1+i2+1]
        }
      }
    }
  }
  
  return(1/FR)
}

## Computes the quantile of a given Laguerre density by numerically evaluating
## the density between 0 and a user specified upper bound ub on a equidistant
## grid of length L. The returned value equals the average of the two closest
## quantiles on the grid.
## Input:
##  theta - Parameter of Laguerre density, cf. laguerre_density_pos
##  tau   - Quantile to be found
##  ub    - Upper bound of the evaluating grid (default is 20)
##  L     - Number of points in the evaluating grid (default is 10000)
## Output: An estimate of the tau-th quantile.
laguerre_quantile <- function(theta,tau,ub=20,L=10000) {
  z <- seq(from=0,to=ub,length.out=L)
  distr <- laguerre_distr_pos(z,theta)
  
  q1 <- z[min(which(distr>=tau))]
  q2 <- z[max(which(distr<=tau))]
  return((q1+q2)/2)
}

## Compute the Laguerre Polynomial of given order.
## Input:
##  z - Points at which the Laguerre polynomial shall be evaluated (must be non-
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
## rC      - Vector of the same length as x1 and x2. Each entry corresponds to
##           an observation of the exponential parameter at the observed
##           location.
## m1,m2   - Two integers (both at least 1) which specify the size of the fitted
##           model. See also the description of theta.
## N       - This is the same as int_prec from fit_laguerre.
## Output: The un-normalized (i.e. not divided by the number of observations)
##         negative log-likelihood corresponding to the given observations
##         and models.
whole_inc_likelihood <- function(theta,x1,x2,w,rC,m1,m2,N=500) {
  n <- length(x1)
  
  theta1_polar <- theta[1:m1]
  theta2_polar <- theta[(m1+1):(m1+m2)]
  
  theta1 <- polar2rect(1,theta1_polar)
  theta2 <- polar2rect(1,theta2_polar)
  
  out <- .Call("whole_inc_likelihood",theta1,theta2,x1,x2,w,rC,as.integer(N))
  
  return(-out)
}

## Computes the best Laguerre approximation to a density for a given degree.
## Input
##  fdens - Function which shall be approximated. The first argument of fdens
##          must be the vector of evaluation points and the output of fdens must
##          be a vector of the same length containing the respective evaluations.
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

## This function computes the Hellinger distance between a Laguerre density and
## another density td.
## Input:
##  theta - Parameters of the Laguerre density to use
##  td    - Density to which the Laguerre density will be compared. It must
##          provided in the same way as fdens in laguerre_approx.
##  ...   - Further arguments which are passed to td
## Output: Value of the Hellinger distance
hellinger_distance <- function(theta,td,...) {
  return(cubintegrate(function(x,...) { return((sqrt(laguerre_density_pos(x,theta))-sqrt(td(x,...)))^2) },lower=0,upper=Inf,...)$integral)
}

## This function computes the Hellinger distance between two Laguerre densities.
## Input:
##  theta1, theta2 - Parameter vectors of the two Laguerre polynomials to
##                   compare.
## Output: Value of the Hellinger distance
hellinger_distance_lag <- function(theta1,theta2) {
  return(cubintegrate(function(x) { return((sqrt(laguerre_density_pos(x,theta1))-sqrt(laguerre_density_pos(x,theta2)))^2) },lower=0,upper=Inf)$integral)
}


################################################################################
## Internal functions which are normally not direclty called by the user.     ##
################################################################################

my_rect2polar <- function(x) {
  phi0 <- rect2polar(x)$phi
  if(phi0[length(phi0)]>pi) {
    phi <- rect2polar(-x)$phi
  } else {
    phi <- phi0
  }
  
  return(phi)
}