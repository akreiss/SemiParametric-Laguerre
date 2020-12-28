## Libraries
library(doParallel)
library(foreach)
library(nloptr)
library(SphericalCubature)
set.seed(2020)


## Set information
n <- 10 # Number of observations

## Generate Data
w <- rexp(n,rate=0.3820225)

first_infection  <-                 runif(n,min=0,w)
second_infection <- first_infection+rweibull(n,shape=2.826,scale=5.665)
  
S1 <- first_infection +rlnorm(n,meanlog=1.644,sd=0.363)
S2 <- second_infection+rlnorm(n,meanlog=1.644,sd=0.363)

w <- pmin(w,S1)



################################################################################
## Estimation ##################################################################
################################################################################
m1 <- 2
m2 <- 2

opts <- list(algorithm="NLOPT_LN_SBPLX",print_level=0,xtol_rel=0.001,maxeval=10000)

noc <- detectCores()-1
cl <- makeCluster(noc,outfile="progress")
registerDoParallel(cl)

## Prepare output
trials <- 5
out <- vector(mode="list",length=trials)

## Do Estimation
currmin <- Inf
lowest <- 1
for(i in 1:trials) {
  out[[i]] <- nloptr(x0=runif(m1+m2,min=0,max=pi),eval_f=whole_inc_likelihood,lb=rep(0,m1+m2),ub=rep(pi,m1+m2),opts=opts,x1=S1,x2=S2,w=w,m1=m1,m2=m2,N=500,cluster=cl)
  if(out[[i]]$objective<currmin) {
    lowest <- i
    currmin <- out[[i]]$objective
  }
  cat("Trial ",i,".\n")
}
stopCluster(cl)

## Extract Estimate
estimate <-  out[[lowest]]$solution
theta1_est <- polar2rect(1,estimate[     1:(m1   )])
theta2_est <- polar2rect(1,estimate[(m1+1):(m1+m2)])

## Compute estimated densities
x <- seq(from=0,to=15,length.out=1000)
inc_est <- laguerre_density_pos(x,theta1_est)
gen_est <- laguerre_density_pos(x,theta2_est)

## Compute true density
inc <- dlnorm(x,meanlog=1.644,sd=0.363)
gen <- dweibull(x,shape=2.826,scale=5.665)

## Compute the best Laguerre approximation to the truth
best_theta_inc <- laguerre_approx(dlnorm,m=m1,meanlog=1.644,sd=0.363)
best_theta_gen <- laguerre_approx(dweibull,m=m2,shape=2.826,scale=5.665)

best_inc <- laguerre_density_pos(x,best_theta_inc)
best_gen <- laguerre_density_pos(x,best_theta_gen)

## Plot
dev.new()
par(mfrow=c(1,2))
plot(x,inc,main="Incubation Period",type="l")
lines(x,inc_est,lty=2)
lines(x,best_inc,lty=3)
plot(x,gen,main="Generation Time",type="l")
lines(x,gen_est,lty=2)
lines(x,best_gen,lty=3)

################################################################################
## Model Selection #############################################################
################################################################################
trials <- 5
m1_max <- 4
m2_max <- 4
out_matrix <- matrix(0,nrow=m1_max,ncol=m2_max)
opts <- list(algorithm="NLOPT_LN_SBPLX",print_level=3,xtol_rel=0.001,maxeval=10000)
out <- vector(mode="list",length=trials)

noc <- detectCores()-1
cl <- makeCluster(noc,outfile="progress")
registerDoParallel(cl)

for(m1 in 1:m1_max) {
  for(m2 in 1:m2_max) {
    for(i in 1:trials) {
      currmin <- Inf
      lowest <- 0
      out[[i]] <- nloptr(x0=runif(m1+m2,min=0,max=pi),eval_f=whole_inc_likelihood,lb=rep(0,m1+m2),ub=rep(pi,m1+m2),opts=opts,x1=S1,x2=S2,w=w,m1=m1,m2=m2,N=500,cluster=cl)
      if(out[[i]]$objective<currmin) {
        lowest <- i
        currmin <- out[[i]]$objective
      }
      cat("Trial ",i,".\n")
    }
    out_matrix[m1,m2] <- out[[lowest]]$objective
  }
}
stopCluster(cl)

BIC <- matrix(0,ncol=m2_max,nrow=m1_max)
AIC <- matrix(0,ncol=m2_max,nrow=m1_max)
for(m1 in 1:m1_max) {
  for(m2 in 1:m2_max) {
    BIC[m1,m2] <- (m1+m2+2)*log(n)-2*(-out_matrix[m1,m2])
    AIC[m1,m2] <- (m1+m2+2)       -2*(-out_matrix[m1,m2])
  }
}
