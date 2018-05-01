#' An example of using maximum likelihood estimation to fit the geometric brownian motions to simulated stock prices
#'
#' Univariate Model U2: dX = muX dt + sigma X dW
#'  2 parameters to be estimated: (b=mu,d=sigma)
#
#' Works for all values of d > 0
#' @keywords Geometric Brownian Motion

source('R/inc.R')
#source('R/BlackScholes.R')    # exact likelihood function for diagnostics
#source('ModelU2.R') # approximate likelihood function

set.seed(22)
eps<-1e-3
args<-list()
args$plot <- TRUE


args$nloptr<-list(maxiter=100,
                  method = "NLOPT_LD_MMA", #'NLOPT_LD_TNEWTON', 'NLOPT_LN_COBYLA'
                  l = c(eps,0.1+ eps),         #b,d
                  u = c(1.0-eps,1.0-eps),      #b,d
                  eval_g_ineq = NULL,
                  eval_jac_g_ineq = NULL,
                  check_derivatives = TRUE,
                  check_derivatives_print="all",
                  ftol_abs=1e-16,
                  xtol_rel=1e-16,
                  ftol_rel=1e-16,
                  print_level=3)

args$DEoptim$maxiter    <- 0
args$DEoptim$population <- 100
args$DEoptim$strategy   <- 2

args$mode <- 'direct'

##Step 1: Simulating GBM data
b<-0.5
d<-0.2

x_0<-10
param_0 <-c(b,d)
del <- 1/252 # weekly data: del = 1/52

# read your data here: create a vector x where x(1) = first observations, x(n) = last observations
# instead, for illustrations purposes, let's just simulate a series from the model


n <- 10000
n.burnin <- 500
factor <- 100
nsimul <- factor*(n+n.burnin)
delta <- del/factor
param <- param_0

rndn <- rnorm(nsimul)
xsimul <- rep(nsimul,0)
xsimul[1] <- x_0
for (i in 2:nsimul){
  xsimul[i] <- xsimul[i-1] + (b*xsimul[i-1])*delta + d*xsimul[i-1]*sqrt(delta)*rndn[i]
}
x <- rep(n,0)
for (i in 1:n){
  x[i] <- xsimul[n.burnin*factor+1+(i-1)*factor]
}


## Step 2: optionally change the initial parameter to observe how stable the calibration is when the initial value for the optimization
## is not the solution parameter
m <- 1 #1000  # number of simulations
params <-matrix(0,m,length(param))

for (i in 1:m){
  for (j in 1:length(param)){
    param[j]<-runif(1,args$nloptr$l[j], args$nloptr$u[j])
  }
  print(i)
## Step 3: estimate the MLE parameters
  output <- mle(ModelU2,x,del,param,args)
  params[i,] <- output$solution
}
## Step 4: compute diagnostic information and plot the log likelihood function
res <- summary(ModelU2,x,del,output$solution,args)


# Step 5: Perform additional diagnostics comparing the approximate likelihood
## function with the exact value
objfun <- function(param){
  objfun <- logdensity2loglik(ModelU2,x,del,param,args)$llk
}

exact <- objfun(param_0)
grad(objfun,param_0)
hessian(objfun, param_0)

print(paste("Maximum log likelihood is ", exact))
print(paste("Standard Error Estimate ", res$se))
print(paste("Huber Sandwich Error Estimate ", res$se_robust))

exact<-list()
exact$score <- exactscore(x,del,output$solution)
exact$InfoMatrix <- exactinformationmatrix(x,del,output$solution)
exact$H <- exacthessian(x,del,output$solution)
Variance <- solve(exact$InfoMatrix)
invH <- solve(exact$H)
Variance_Robust <- invH %*% exact$InfoMatrix %*% t(invH)

exact$se <-sqrt(diag(Variance))
exact$se_robust <- sqrt(diag(Variance_Robust))
print(paste("Exact Standard Error ", res$se))
print(paste("Exact Huber Sandwich Error ", res$se_robust))

print(paste("Error in S.E. Estimate ", norm(as.matrix(res$se-exact$se))))
print(paste("Error in H.S.E. Estimate ", norm(as.matrix(res$se_robust-exact$se_robust))))
print(paste("L2 Norm of Score Error ", norm(as.matrix(res$score-exact$score))))
print(paste("L2 Norm of Hessian Error ", norm(as.matrix(res$H-exact$H))))
print(paste("L2 Norm of Information matrix Error ", norm(as.matrix(res$InfoMatrix-exact$InfoMatrix))))


