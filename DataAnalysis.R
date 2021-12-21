## Libraries
library(doParallel)
library(foreach)
library(cubature)
library(nloptr)
library(SphericalCubature)

dyn.load("c_function.dll")
source('./functions.R')


################################################################################
## Load and Prepare Dataset ####################################################
################################################################################
set.seed(2021)

## Load Data
data <- read.csv('extended_data.csv',sep=';')
n <- dim(data)[1]

## Read Information
first_exposure         <- data[,4]
last_exposure          <- data[,5]
infected_last_exposure <- data[,7]
symptom_onset_infector <- data[,2]
symptom_onset_infected <- data[,3]

## Replace NAs by natural bounds
last_exposure <- pmin(last_exposure,symptom_onset_infector,symptom_onset_infected,infected_last_exposure,na.rm=TRUE)
fe_na_ind <- which(is.na(first_exposure))
first_exposure[fe_na_ind] <- symptom_onset_infector[fe_na_ind]-2*30

## Compute data for estimation
S1 <- symptom_onset_infector-first_exposure+runif(n)
S2 <- symptom_onset_infected-first_exposure+runif(n)
w  <-          last_exposure-first_exposure+runif(n)
rC <- rep(0,n)
rC_ind <- which(data[,8]=="YES")
rC[rC_ind] <- log(2)/5


################################################################################
## Model Selection #############################################################
################################################################################
set.seed(2021)

trials <- 10
m1_max <- 6
m2_max <- 6
out_matrix <- matrix(0,nrow=m1_max,ncol=m2_max)
optimal_outputs <- vector(mode="list",length=m1_max*m2_max)


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

save(MSL,BIC,AIC,file="DataAnalysis_ModelSelection.RData")




################################################################################
## Make Simulations for the test.                                             ##
################################################################################
set.seed(2021)

## Set information
n <- length(symptom_onset_infector)  # Number of patients
N <- 1000 # Number of repetitions

## Prepare output
symptom_onset_infector_sim <- matrix(0,nrow=n,ncol=N)
symptom_onset_infected_sim <- matrix(0,nrow=n,ncol=N)
exposure_window_sim        <- matrix(0,nrow=n,ncol=N)
location_sim               <- matrix(0,nrow=n,ncol=N)

## Generate Data
for(k in 1:N) {
  exposure_window_sim[,k] <- rexp(n,rate=0.3820225)
  location_sim[,k]        <- sample(c(0,1),n,replace=TRUE,prob=c(26/40,14/40))
  
  no_exp_growth <- which(location_sim[,k]==0)
  exp_growth    <- which(location_sim[,k]==1)
  
  first_infection  <- rep(0,n)
  first_infection[no_exp_growth] <- runif(length(no_exp_growth),min=0,max=exposure_window_sim[no_exp_growth,k])
  first_infection[   exp_growth] <- exposure_window_sim[exp_growth,k]-(rexp( length(   exp_growth),rate=log(2)/5) %% exposure_window_sim[exp_growth,k])
  second_infection <- first_infection+rweibull(n,shape=2.826,scale=5.665)
  
  symptom_onset_infector_sim[,k] <- floor(first_infection +rlnorm(n,meanlog=1.644,sd=0.363))+runif(n)
  symptom_onset_infected_sim[,k] <- floor(second_infection+rlnorm(n,meanlog=1.644,sd=0.363))+runif(n)
  exposure_window_sim[,k] <- floor(exposure_window_sim[,k])+runif(n)
}

## Do Estimation 
m1 <- 4
m2 <- 3
trials <- 10

noc <- detectCores()
cl <- makeCluster(noc,outfile="progress")
registerDoParallel(cl)

clusterCall(cl,source,'functions.R')
clusterCall(cl,dyn.load,"c_function.dll")

out <- foreach(k=1:N,.packages=c('SphericalCubature','nloptr')) %dopar% {
  set.seed(2021+k)
  
  if(k %% 50 ==0) {
    cat("Do Step ",k," of ",N,".\n")
    
  }
  
  rC_sim <- rep(0,n)
  rC_sim[location_sim[,k]==1] <- log(2)/5
  
  fit_laguerre(S1=symptom_onset_infector_sim[,k],S2=symptom_onset_infected_sim[,k],w=pmin(exposure_window_sim[,k],symptom_onset_infector_sim[,k]),rC=rC_sim,m1=m1,m2=m2,glob_repeats=trials)
}
stopCluster(cl)

## Format output
estimates <- matrix(0,nrow=N,ncol=m1+m2+2)

for(k in 1:N) {
  estimates[k,     1:(m1   +1)] <- out[[k]]$theta1
  estimates[k,(m1+2):(m1+m2+2)] <- out[[k]]$theta2
}

## Find closest Laguerre poynomial
best_theta_inc <- laguerre_approx(dlnorm  ,m=m1,meanlog=1.644,sd=0.363)
best_theta_gen <- laguerre_approx(dweibull,m=m2,shape=2.826,scale=5.665)

## Compute differences in Hellinger distance to closest Laguerre density
hellinger_inc <- rep(0,N)
hellinger_gen <- rep(0,N)

for(i in 1:N) {
  hellinger_inc[i] <- hellinger_distance_lag(best_theta_inc,estimates[i,     1:(m1   +1)])
  hellinger_gen[i] <- hellinger_distance_lag(best_theta_gen,estimates[i,(m1+2):(m1+m2+2)])
}
save(estimates,hellinger_gen,hellinger_inc,best_theta_inc,best_theta_gen,file="simulated_hellinger.RData")


################################################################################
## Visualize ###################################################################
################################################################################

## Print BIC values ############################################################
library(xtable)
print.xtable(xtable(BIC))

## Show estimates ##############################################################
m1 <- 5
m2 <- 3

ind <- MSL$index_matrix[m1,m2]

theta1_est <- MSL$fit_results[[ind]]$theta1
theta2_est <- MSL$fit_results[[ind]]$theta1

x <- seq(from=0,to=15,length.out=1000)
inc_est <- laguerre_density_pos(x,theta1_est)
gen_est <- laguerre_density_pos(x,theta2_est)

inc <- dlnorm(x,meanlog=1.644,sd=0.363)
gen <- dweibull(x,shape=2.826,scale=5.665)

best_inc <- laguerre_density_pos(x,best_theta_inc)
best_gen <- laguerre_density_pos(x,best_theta_gen)

dev.new(width=1920,height=1080,unit="px",noRStudioGD = TRUE)
par(mfrow=c(1,2))
plot(x,inc,main="Incubation Period",type="l",ylim=c(0,0.3))
lines(x,inc_est,lty=2)
lines(x,best_inc,lty=3)
legend("topright",lty=c(1,2,3),legend=c("Parametric","Laguerre","Best-Laguerre"))
plot(x,gen,main="Generation Time",type="l",ylim=c(0,0.4))
lines(x,gen_est,lty=2)
lines(x,best_gen,lty=3)
legend("topright",lty=c(1,2,3),legend=c("Parametric","Laguerre","Best-Laguerre"))

## Compute differences in Hellinger distance to best density
inc_real_dist <- hellinger_distance_lag(best_theta_inc,theta1_est)
gen_real_dist <- hellinger_distance_lag(best_theta_gen,theta2_est)

N <- length(hellinger_gen)
p1 <- sum(hellinger_inc>inc_real_dist)/N*100
p2 <- sum(hellinger_gen>gen_real_dist)/N*100
p3 <- sum(hellinger_gen>gen_real_dist & hellinger_inc>inc_real_dist)/N*100

dev.new(width=1920,height=1080,unit="px",noRStudioGD = TRUE)
par(mfrow=c(1,2))
hist(hellinger_inc,main="Incubation Period",xlab="Hellinger Distance",ylab="",xlim=c(0,0.06))
abline(v=inc_real_dist)
hist(hellinger_gen,main="Generation Time",xlab="Hellinger Distance",ylab="")
abline(v=gen_real_dist)

## p-values
sum(hellinger_inc>=inc_real_dist)/N*100
sum(hellinger_gen>=gen_real_dist)/N*100
sum(hellinger_inc>=inc_real_dist & hellinger_gen>=gen_real_dist)/N*100

