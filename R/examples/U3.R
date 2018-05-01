# Model = linear drift, CEV diffusion (Univariate Model U3), dX = b(a - X) dt + c X^d dW
# 4 parameters to be estimated: (a,b,c,d)
#
# Works for all values of d > 1/2 using two different forms of the loglik for d=1 and d~=1 in file closedform_expansion_cev_loglik
# starting values for MLE algorithm and simulated series
source('R/inc.R')
#source('ModelU3.R')



eps<-1e-3
args<-list()
args$plot <- TRUE


args$nloptr<-list(maxiter=100,
                  method = "NLOPT_LD_MMA", #NLOPT_LD_TNEWTON", #"NLOPT_LN_COBYLA",
                  l = c(eps,eps,0.1+eps,0.5+eps),         #a,b,c,d
                  u = c(0.5-eps,1.0-eps,1.0-eps,1.1-eps), #a,b,c,d
                  eval_g_ineq = NULL,
                  eval_jac_g_ineq = NULL,
                  check_derivatives = TRUE,
                  check_derivatives_print="all",
                  ftol_abs=1e-16,
                  xtol_rel=1e-16,
                  ftol_rel=1e-16,
                  print_level=3)


args$DEoptim$maxiter    <- 50
args$DEoptim$population <- 100
args$DEoptim$strategy   <- 2
args$mode <- 'direct'
##Step 1: Simulating CEV data

a<-0.2
b<-0.5
c<-0.9
d<-0.9

param_0 <-c(a,b,c,d)


# daily data: del = 1/252
del <- 1/252

# read your data here: create a vector x where x(1) = first observations, x(n) = last observations
# instead, for illustrations purposes, let's just simulate a series from the model

set.seed(11)
x_0 <- a
factor <- 10
n <- 5000
n.burnin <- 0
nsimul <- factor*n
delta <- del/factor
param <- param_0

exact <- logdensity2loglik(ModelU3,x,del,param_0,args)
print(paste("Maximum log likelihood is ", exact))

rndn <- rnorm(n.burnin + nsimul)
xsimul <- rep(nsimul,0)
xsimul[1] <- x_0
for (i in 2:nsimul){
  xsimul[i] <- xsimul[i-1] + b*(a-xsimul[i-1])*delta + c*xsimul[i-1]^d*sqrt(delta)*rndn[n.burnin+i] + 0.5*c^2*d*xsimul[i-1]^(2*d-1)*(rndn[n.burnin+i]*rndn[n.burnin+i]-1)*delta
}
x <- rep(n,0)

for (i in 1:n){
  x[i] <- xsimul[1+(i-1)*factor]
}


summary(ModelU3,x,del,param_0,args)

load('U3.rda')
## Step 2: optionally change the initial parameter to observe how stable the calibration is when the initial value for the optimization
## is not the solution parameter

m <- 1 #1000  # number of simulations
params <-matrix(0,m,length(param_0))
set.seed(11)
for (k in 1:m){

  # randomize the initial condition
  for (j in 1:length(param_0)){
    param[j]<-runif(1,args$nloptr$l[j], args$nloptr$u[j])
  }
  print(i)

  ## Step 3: estimate the MLE parameters
  output <- mle(ModelU3,x,del,param,args)
  params[k,] <- output$solution
}
## Step 4: compute diagnostic information and plot the log likelihood function

res <- summary(ModelU3,x,del,output$solution,args)

print(res$se)
print(res$se_robust)



