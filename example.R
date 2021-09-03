## Libraries
library(doParallel)
library(foreach)
library(nloptr)
library(SphericalCubature)
set.seed(2020)


## Set information
n <- 10 # Number of observations
r <- log(2)/5 # Exponential growth rate

## Generate Data
w <- rexp(n,rate=0.3820225)
location<- sample(c(0,1),n,replace=TRUE,prob=c(26/40,14/40))
rC <- rep(0,n)
rC[location==1] <- r

no_exp_growth <- which(location==0)
exp_growth    <- which(location==1)

first_infection  <- rep(0,n)
first_infection[no_exp_growth] <- runif(length(no_exp_growth),min=0,max=w[no_exp_growth])
first_infection[   exp_growth] <- w[exp_growth]-(rexp(length(exp_growth),rate=r) %% w[exp_growth])

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
  out[[i]] <- nloptr(x0=runif(m1+m2,min=0,max=pi),eval_f=whole_inc_likelihood,lb=rep(0,m1+m2),ub=rep(pi,m1+m2),opts=opts,x1=S1,x2=S2,w=w,rC=rC,m1=m1,m2=m2,N=500,cluster=cl)
  if(out[[i]]$objective<currmin) {
    lowest <- i
    currmin <- out[[i]]$objective
  }
  cat("Trial ",i,".\n")
}
stopCluster(cl)

################################################################################
## Visualize density estimate ##################################################
################################################################################

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
## Compute true and estimated quantiles and true and estimated basic repro-   ##
## duction numbers.                                                           ##
################################################################################

## Quantiles to test
tau <- c(0.3,0.5,0.7,0.9)

## True quantiles
true_quant_inc <- qlnorm(tau,meanlog=1.644,sd=0.363)
true_quant_gen <- qweibull(tau,shape=2.826,scale=5.665)

## Quantiles of the best Laguerre Approximation
best_quant_inc <- rep(0,length(tau))
best_quant_gen <- rep(0,length(tau))
for(k in 1:length(tau)) {
  best_quant_inc[k] <- lag_quantile(best_theta_inc,tau[k])
  best_quant_gen[k] <- lag_quantile(best_theta_gen,tau[k])
}

## Estimated Quantiles
est_quant_inc <- rep(0,length(tau))
est_quant_gen <- rep(0,length(tau))

for(k in 1:length(tau)) {
  est_quant_inc[k] <- lag_quantile(theta1_est,tau[k])
  est_quant_gen[k] <- lag_quantile(theta2_est,tau[k])
}

## Print result
cat("Quantiles of interest:\n")
print(tau)
cat("True quantiles of Incubation Time distribution:\n")
print(true_quant_inc)
cat("Estimated Quantiles of the Incubation Time distribution:\n")
print(est_quant_inc)
cat("Quantiles of the best Laguerre Approximation to the Incubation Time Dsitribution:\n")
print(best_quant_inc)
cat("True quantiles of Generation Time distribution:\n")
print(true_quant_gen)
cat("Estimated Quantiles of the Generation Time distribution:\n")
print(est_quant_gen)
cat("Quantiles of the best Laguerre Approximation to the Generation Time Dsitribution:\n")
print(best_quant_gen)


## Compute true R0
helpf <- function(x) {
  return(exp(-r*x)*dweibull(x,shape=2.826,scale=5.665))
}
R0 <- 1/integrate(helpf,lower=0,upper=Inf)$value

## Compute estimated R0
R0_est <- laguerre_R0(theta2_est,r)

## Print result
cat("True R0 is ",R0," while the estimate is ",R0_est,".\n")


################################################################################
## Model Selection #############################################################
################################################################################
trials <- 5
m1_max <- 4
m2_max <- 4
out_matrix <- matrix(0,nrow=m1_max,ncol=m2_max)
opts <- list(algorithm="NLOPT_LN_SBPLX",print_level=0,xtol_rel=0.001,maxeval=10000)
out <- vector(mode="list",length=trials)

noc <- detectCores()-1
cl <- makeCluster(noc,outfile="progress")
registerDoParallel(cl)

for(m1 in 1:m1_max) {
  for(m2 in 1:m2_max) {
    for(i in 1:trials) {
      currmin <- Inf
      lowest <- 0
      out[[i]] <- nloptr(x0=runif(m1+m2,min=0,max=pi),eval_f=whole_inc_likelihood,lb=rep(0,m1+m2),ub=rep(pi,m1+m2),opts=opts,x1=S1,x2=S2,w=w,rC=rC,m1=m1,m2=m2,N=500,cluster=cl)
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
