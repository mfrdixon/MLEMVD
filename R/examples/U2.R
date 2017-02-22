#' An example of using maximum likelihood estimation to fit the geometric brownian motions to simulated stock prices
#' 
#' Univariate Model U2: dX = muX dt + sigma X dW
#'  2 parameters to be estimated: (b=mu,d=sigma)
#
#' Works for all values of d > 0 
#' @keywords Geometric Brownian Motion

source('inc.R')
source('utilities/BlackScholes.R')    # exact likelihood function for diagnostics
source('models/univariate/ModelU2.R') # approximate likelihood function

eps<-1e-3
args<-list()
args$plot <- TRUE
args$nloptr<-list(maxiter=100,
                  method = 'NLOPT_LN_COBYLA',
                  l = c(-1.0+eps,0.1+ eps),    #b,d
                  u = c(1.0-eps,1.0-eps),      #b,d
                  eval_g_ineq = NULL,
                  ftol_abs=1e-14, 
                  xtol_rel=1e-14, 
                  ftol_rel=1e-11, 
                  print_level=3)

args$DEoptim$maxiter    <- 10
args$DEoptim$population <- 100
args$DEoptim$strategy   <- 2

args$mode <- 'direct'

##Step 1: Simulating CEV data
rate<-0.01
a<-rate
b<-0.5
d<-0.2
del<-1/52 
x_0<-10
param_0 <-c(b,d)
del <- 1/52 # weekly data: del = 1/52

# read your data here: create a vector x where x(1) = first observations, x(n) = last observations
# instead, for illustrations purposes, let's just simulate a series from the model

n <- 500
factor <- 10
nsimul <- factor*n
delta <- del/factor

rndn <- rnorm(nsimul)
xsimul <- rep(nsimul,0)
xsimul[1] <- x_0
for (i in 2:nsimul){
  xsimul[i] <- xsimul[i-1] + (a + b*xsimul[i-1])*delta + d*xsimul[i-1]*sqrt(delta)*rndn[i]
}
x <- rep(n,0)
x[1] <- x_0
for (i in 2:n){
  x[i] <- xsimul[1+(i-1)*factor]
}


## Step 2: optionally change the initial parameter to observe how stable the calibration is when the initial value for the optimization 
## is not the solution parameter
param <- c(0.1,0.4)

## Step 3: estimate the MLE parameters 
output <- mle(ModelU2,x,del,param,args)
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


