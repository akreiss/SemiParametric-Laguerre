## Load Code
source('./functions.R')
dyn.load("c_function.dll")

## Libraries
library(doParallel)
library(foreach)
library(nloptr)
library(SphericalCubature)


## Set information
n <- 40  # Number of patients
N <- 1000 # Number of repetitions
trials <- 10 # Number of optimisation tasks


################################################################################
## Generate Data                                                              ##
################################################################################
set.seed(2021)

## Prepare output
symptom_onset_infector <- matrix(0,nrow=n,ncol=N)
symptom_onset_infected <- matrix(0,nrow=n,ncol=N)
exposure_window        <- matrix(0,nrow=n,ncol=N)
location               <- matrix(0,nrow=n,ncol=N)

## Generate Data
for(k in 1:N) {
  exposure_window[,k] <- rexp(n,rate=0.3820225)
  location[,k]        <- sample(c(0,1),n,replace=TRUE,prob=c(26/40,14/40))
  
  no_exp_growth <- which(location[,k]==0)
  exp_growth    <- which(location[,k]==1)
  
  first_infection  <- rep(0,n)
  first_infection[no_exp_growth] <- runif(length(no_exp_growth),min=0,max=exposure_window[no_exp_growth,k])
  first_infection[   exp_growth] <- exposure_window[exp_growth,k]-(rexp( length(   exp_growth),rate=log(2)/5) %% exposure_window[exp_growth,k])
  second_infection <- first_infection+rweibull(n,shape=2.826,scale=5.665)

  symptom_onset_infector[,k] <- first_infection +rlnorm(n,meanlog=1.644,sd=0.363)
  symptom_onset_infected[,k] <- second_infection+rlnorm(n,meanlog=1.644,sd=0.363)
}

################################################################################
## Model Selection #############################################################
################################################################################
m1_max <- 4
m2_max <- 4

S1 <- symptom_onset_infector[,1]
S2 <- symptom_onset_infected[,1]
w <- exposure_window[,1]
rC <- rep(0,n)
rC[location[,1]==1] <- log(2)/5

noc <- detectCores()
cl <- makeCluster(noc,outfile="progress")
registerDoParallel(cl)

MSL <- model_selection_laguerre(S1,S2,w,rC,m1_max,m2_max,TRUE,trials,0.0001,10000,500,cl)
MSL <- force_monotonicity(MSL,S1,S2,w,rC,0.0001,10000,500,TRUE)

stopCluster(cl)

BIC <- matrix(0,ncol=m2_max,nrow=m1_max)
AIC <- matrix(0,ncol=m2_max,nrow=m1_max)
for(m1 in 1:m1_max) {
  for(m2 in 1:m2_max) {
    BIC[m1,m2] <- (m1+m2+2)*log(n)-2*(-MSL$L_matrix[m1,m2])
    AIC[m1,m2] <- (m1+m2+2)       -2*(-MSL$L_matrix[m1,m2])
  }
}

save(BIC,AIC,symptom_onset_infector,symptom_onset_infected,exposure_window,location,n,N,trials,file="n40_ModelSelection.RData")

BIC==min(BIC)
AIC==min(AIC)

################################################################################
## Do Estimation ###############################################################
################################################################################
m1 <- 3
m2 <- 2

noc <- detectCores()
cl <- makeCluster(noc,outfile="progress")
registerDoParallel(cl)

clusterCall(cl,source,'functions.R')
clusterCall(cl,dyn.load,"c_function.dll")

out <- foreach(k=1:N,.packages=c('SphericalCubature','nloptr')) %dopar% {
  set.seed(2021+k)

  ## Write starting Time
  if(k==1) {
    starting_time <- Sys.time()
    save(starting_time,file="starting_time.RData")
  }

  ## Print Status
  if(k %% 50 == 0) {
    current_time <- Sys.time()
    load("starting_time.RData")
    dt <- as.numeric(difftime(current_time,starting_time,units="secs"))
    cat("Step ",k," of ",N,". Estimated remaining time: ",dt/k*(N-k)/60,"min\n")
  }

  ## Fit the model  
  rC <- rep(0,n)
  rC[location[,k]==1] <- log(2)/5
  
  fit_laguerre(S1=symptom_onset_infector[,k],S2=symptom_onset_infected[,k],w=pmin(exposure_window[,k],symptom_onset_infector[,k]),rC=rC,m1=m1,m2=m2,glob_repeats=trials)
}
stopCluster(cl)

## Format output
estimates <- matrix(0,nrow=N,ncol=m1+m2+2)

for(k in 1:N) {
  estimates[k,     1:(m1   +1)] <- out[[k]]$theta1
  estimates[k,(m1+2):(m1+m2+2)] <- out[[k]]$theta2
}

## Save Simulation Result
save(symptom_onset_infected,symptom_onset_infector,exposure_window,location,N,n,m1,m2,estimates,file='SimulationStudy_40Patients_32.RData')

################################################################################
## Visualisation ###############################################################
################################################################################

## Print BIC values ############################################################
library(xtable)
print.xtable(xtable(BIC))


## Visualize closest Laguerre Polynomials ######################################
## Compute best approximation to the true densities for various m
m1_max <- 6
m2_max <- 4
x <- seq(from=0,to=15,length.out=1000)

inc_best_approx <- matrix(0,nrow=m1_max+1,ncol=length(x))
gen_best_approx <- matrix(0,nrow=m2_max+1,ncol=length(x))

inc <- dlnorm(x,meanlog=1.644,sd=0.363)
gen <- dweibull(x,shape=2.826,scale=5.665)

for(i in 0:m1_max) {
  best_theta_inc        <- laguerre_approx(dlnorm,m=i,meanlog=1.644,sd=0.363)
  inc_best_approx[i+1,] <- laguerre_density_pos(x,best_theta_inc)
}
for(i in 0:m2_max) {
  best_theta_gen        <- laguerre_approx(dweibull,m=i,shape=2.826,scale=5.665)
  gen_best_approx[i+1,] <- laguerre_density_pos(x,best_theta_gen)
}

## Plot comparison
yRange_inc <- max(inc_best_approx)
yRange_gen <- max(gen_best_approx)

dev.new(width=1920,height=1080,unit="px",noRStudioGD = TRUE)
par(mfrow=c(1,2))
plot(x,inc,main="Incubation Period",xlab="Time in Days",ylab="Density",type="l",ylim=c(0,0.3),lwd=3)
for(i in 0:m1_max) {
  lines(x,inc_best_approx[i+1,],lty=i+2)
}
legend("topright",lty=1:(m1_max+2),lwd=c(3,rep(1,m1_max+1)),legend=c("Truth",0:m1_max))

plot(x,gen,main="Generation Time",xlab="Time in Days",ylab="Density",type="l",ylim=c(0,0.3),lwd=3)
for(i in 0:m2_max) {
  lines(x,gen_best_approx[i+1,],lty=i+2)
}
legend("topright",lty=1:(m2_max+2),lwd=c(3,rep(1,m2_max+1)),legend=c("Truth",0:m2_max))


## Visualize the fitted curves #################################################
## Compute all estimated functions
lower <- 0
upper <- 15
leng <- 1000
x <- seq(from=lower,to=upper,length.out=leng)
incubation <- matrix(0,nrow=N,ncol=leng)
generation <- matrix(0,nrow=N,ncol=leng)

for(k in 1:N) {
  theta1_est <- estimates[k,     1:(m1   +1)]
  theta2_est <- estimates[k,(m1+2):(m1+m2+2)]
  
  incubation[k,] <- laguerre_density_pos(x,theta1_est)
  generation[k,] <- laguerre_density_pos(x,theta2_est)
}


## Compute point-wise quantiles
probs <- c(0.5,2.5,5,7.5,10,12.5,15,17.5,20,80,82.5,85,87.5,90,92.5,95,97.5,99.5)/100
colours <- c("gray90","gray80","gray70","gray60","gray50","gray40","gray30","gray20","gray10")
quantiles_incubation <- apply(incubation,2,quantile,simplify=TRUE,probs=probs)
quantiles_generation <- apply(generation,2,quantile,simplify=TRUE,probs=probs)

## Compute Optimal Laguerre Fits
best_theta_inc <- laguerre_approx(dlnorm,m=m1,meanlog=1.644,sd=0.363)
best_theta_gen <- laguerre_approx(dweibull,m=m2,shape=2.826,scale=5.665)

best_inc <- laguerre_density_pos(x,best_theta_inc)
best_gen <- laguerre_density_pos(x,best_theta_gen)

## Compute true densities
inc <- dlnorm(x,meanlog=1.644,sd=0.363)
gen <- dweibull(x,shape=2.826,scale=5.665)

## Plot results
dev.new(width=1920,height=1080,unit="px",noRStudioGD = TRUE)
par(mfrow=c(1,2))

## Incubation Period
plot(x,inc,ylim=c(0,0.4),main="Incubation Period",xlab="Time in days",ylab="Density",type="n")
K <- dim(quantiles_incubation)[1]
for(k in 1:(K/2)) {
  polygon(c(x,rev(x)),c(quantiles_incubation[k,],rev(quantiles_incubation[K-k,])),col=colours[k],border=NA)
}
lines(x,inc,lwd=2,col="royalblue3")
lines(x,best_inc,lty=4,lwd=2,col="royalblue3")
legend("topright",lty=c(1,4),col=c("royalblue3","royalblue3"),legend=c("Truth","Closest Laguerre"))

## Generation Time
plot(x,gen,ylim=c(0,0.4),main="Generation Time",xlab="Time in days",ylab="Density",type="n")
K <- dim(quantiles_generation)[1]
for(k in 1:(K/2)) {
  polygon(c(x,rev(x)),c(quantiles_generation[k,],rev(quantiles_generation[K-k,])),col=colours[k],border=NA)
}
lines(x,gen,lwd=2,col="royalblue3")
lines(x,best_gen,lty=4,lwd=2,col="royalblue3")
legend("topright",lty=c(1,4),col=c("royalblue3","royalblue3"),legend=c("Truth","Closest-Laguerre"))


## Compare Hellinger Distances to truth and to closest Laguerre, lognormal #####
## Compute differences in Hellinger distance to the truth
hellinger_to_truth <- matrix(0,ncol=2,nrow=N)
for(i in 1:N) {
  hellinger_to_truth[i,1] <- hellinger_distance(estimates[i,     1:(m1   +1)],dlnorm  ,meanlog=1.644,sd   =0.363)
  hellinger_to_truth[i,2] <- hellinger_distance(estimates[i,(m1+2):(m1+m2+2)],dweibull,  shape=2.826,scale=5.665)
}

dev.new(width=1920,height=1080,unit="px",noRStudioGD = TRUE)
par(mfrow=c(1,2))
hist(hellinger_to_truth[,1],main="Incubation Period",xlab="Squared Hellinger Distance",xlim=c(0,0.35))
hist(hellinger_to_truth[,2],main="Generation Time"  ,xlab="Squared Hellinger Distance",xlim=c(0,0.35))

quantile(hellinger_to_truth[,1],probs=c(10:90,95,99)/100)
quantile(hellinger_to_truth[,2],probs=c(10:90,95,99)/100)

## Compare estimated and true R0 ###############################################
r <- log(2)/5 # Growth of pandemic

## Compute true R0
helpf <- function(x) {
  return(exp(-r*x)*dweibull(x,shape=2.826,scale=5.665))
}
R0 <- 1/integrate(helpf,lower=0,upper=Inf)$value

## Compute estimated R0
R0_est <- rep(0,N)
for(k in 1:N) {
  R0_est[k] <- laguerre_R0(estimates[k,(m1+2):(m1+m2+2)],r)
}

## Visualize
dev.new()
x <- seq(from=min(R0_est)-1,to=max(R0_est)+1,length.out=1000)
y <- dnorm(x,mean=mean(R0_est),sd=sd(R0_est))
hist(R0_est,freq=FALSE,main="Basic Reproduction Number R0",xlab="R0")
abline(v=R0,lty=2)
lines(x,y)


## Compare estimated and true quantiles ########################################
## Quantiles to test
tau <- c(0.3,0.5,0.7,0.9)

## True quantiles
true_quant_inc <- qlnorm(tau,meanlog=1.644,sd=0.363)
true_quant_gen <- qweibull(tau,shape=2.826,scale=5.665)

## Quantiles of the best Laguerre Approximation
best_quant_inc <- rep(0,length(tau))
best_quant_gen <- rep(0,length(tau))
for(k in 1:length(tau)) {
  best_quant_inc[k] <- laguerre_quantile(best_theta_inc,tau[k])
  best_quant_gen[k] <- laguerre_quantile(best_theta_gen,tau[k])
}

## Estimated Quantiles
est_quant_inc <- matrix(0,nrow=N,ncol=length(tau))
est_quant_gen <- matrix(0,nrow=N,ncol=length(tau))

for(i in 1:N) {
  for(k in 1:length(tau)) {
    est_quant_inc[i,k] <- laguerre_quantile(estimates[i,     1:(m1   +1)],tau[k])
    est_quant_gen[i,k] <- laguerre_quantile(estimates[i,(m1+2):(m1+m2+2)],tau[k])
  }
}

## Visualize
x <- seq(from=0,to=20,length.out=1000)
for(k in 1:length(tau)) {
  dev.new(width=1920,height=1080,unit="px",noRStudioGD = TRUE)
  par(mfrow=c(1,2))
  hist(est_quant_inc[,k],freq = FALSE,main=sprintf("Incubation Time - tau=%.2f",tau[k]),xlab="Quantile",ylim=c(0,1.6))
  abline(v=true_quant_inc[k],lty=2)
  abline(v=best_quant_inc[k],lty=3)
  legend("topleft",legend=c("Normal Approx.","True Quantile","Best Quantile"),lty=c(1,2,3))
  lines(x,dnorm(x,mean=mean(est_quant_inc[,k]),sd=sd(est_quant_inc[,k])))
  
  hist(est_quant_gen[,k],freq = FALSE,main=sprintf("Generation Time - tau=%.2f",tau[k]),xlab="Quantile")
  abline(v=true_quant_gen[k],lty=2)
  abline(v=best_quant_gen[k],lty=3)
  legend("topleft",legend=c("Normal Approx.","True Quantile","Best Quantile"),lty=c(1,2,3))
  lines(x,dnorm(x,mean=mean(est_quant_gen[,k]),sd=sd(est_quant_gen[,k])))
}