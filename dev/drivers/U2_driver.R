# Model =  geometric brownian motion (Univariate Model U1), dX = muX dt + sigma X dW
# 2 parameters to be estimated: (b=mu,d=sigma)
#
# Works for all values of d > 0 
source('inc.R')
source('models/univariate/ModelU2.R')

eps<-1e-3
args<-list()
args$plot <- TRUE
args$nloptr<-list(maxiter=1000,
                  method = 'NLOPT_GN_ISRES',
                  l = c(-1.0+eps,0.1+ eps),    #b,d
                  u = c(1.0-eps,1.0-eps),      #b,d
                  eval_g_ineq = NULL,
                  ftol_abs=1e-14, 
                  xtol_rel=1e-14, 
                  ftol_rel=1e-11, 
                  print_level=3)

# NLOPT_GN_DIRECT
# NLOPT_GN_DIRECT_L
# NLOPT_GN_DIRECT_L_RAND
# NLOPT_GN_DIRECT_NOSCAL
# NLOPT_GN_DIRECT_L_NOSCAL
# NLOPT_GN_DIRECT_L_RAND_NOSCAL
# NLOPT_GN_ORIG_DIRECT
# NLOPT_GN_ORIG_DIRECT_L
# NLOPT_GD_STOGO
# NLOPT_GD_STOGO_RAND
# NLOPT_LD_SLSQP
# NLOPT_LD_LBFGS_NOCEDAL
# NLOPT_LD_LBFGS
# NLOPT_LN_PRAXIS
# NLOPT_LD_VAR1
# NLOPT_LD_VAR2
# NLOPT_LD_TNEWTON
# NLOPT_LD_TNEWTON_RESTART
# NLOPT_LD_TNEWTON_PRECOND
# NLOPT_LD_TNEWTON_PRECOND_RESTART
# NLOPT_GN_CRS2_LM
# NLOPT_GN_MLSL
# NLOPT_GD_MLSL
# NLOPT_GN_MLSL_LDS
# NLOPT_GD_MLSL_LDS
# NLOPT_LD_MMA
# NLOPT_LN_COBYLA
# NLOPT_LN_NEWUOA
# NLOPT_LN_NEWUOA_BOUND
# NLOPT_LN_NELDERMEAD
# NLOPT_LN_SBPLX
# NLOPT_LN_AUGLAG
# NLOPT_LD_AUGLAG
# NLOPT_LN_AUGLAG_EQ
# NLOPT_LD_AUGLAG_EQ
# NLOPT_LN_BOBYQA
# NLOPT_GN_ISRES


args$DEoptim$maxiter    <- 10
args$DEoptim$population <- 100
args$DEoptim$strategy   <- 2
##Step 1: Simulating CEV data

b<-0.5
d<-0.2
del<-1/52 
x_0<-10
param_0 <-c(b,d)
#set.seed(177)
# weekly data: del = 1/52
del <- 1/52

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



objfun <- function(param){
  objfun <- logdensity2loglik(ModelU2,x,del,param,args)
} 



param <- c(0.1,0.4)
output <- mymle(ModelU2,x,del,param,args)
res <- summary(output$solution,ModelU2,x,del,args)

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
print(paste("L2 Norm of Information matrix ", norm(as.matrix(res$InfoMatrix-exact$InfoMatrix))))

exactloglik<-function(x,del,param){
  
  mu    <- param[1]
  sigma <- param[2]
  n<-length(x)-1
  output <-0 
  for (i in 1:n){
    output<-output+log(2*pi*del*sigma^2*x[i+1]^2)+ (log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)^2/(sigma^2*del)
  }
  output <- -0.5*output
  return(output)
}
exactscore<-function(x,del,param){
  
  mu    <- param[1]
  sigma <- param[2]
  n<-length(x)-1
  output <-rep(0,2) 
  output[2] <- -n/sigma
  for (i in 1:n){
    output[1]<-output[1] + (log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)/sigma^2
    output[2]<-output[2] -((log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)/(sigma*sqrt(del)))*((-(log(x[i+1]/x[i]) -mu*del)/(sqrt(del)*sigma^2)) + sqrt(del)/2)
  }
  return(output)
}

exacthessian<-function(x,del,param){
  
  mu    <- param[1]
  sigma <- param[2]
  n<-length(x)-1
  output <-matrix(rep(0,4),2,2) 
  output[2,2] <- n/(sigma^2)
  for (i in 1:n){
    output[1,1]<-output[1,1] -del/sigma^2
    output[1,2]<-output[1,2] -2*(log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)/sigma^3 + del/sigma
    g<- (log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)/(sigma*sqrt(del))
    h<- (-(log(x[i+1]/x[i]) -mu*del)/(sqrt(del)*sigma^2)) + sqrt(del)/2
    output[2,2]<-output[2,2] -h^2 +2*g*(mu*del-log(x[i+1]/x[i]))/(sqrt(del)*sigma^3)
  }
  output[2,1]<-output[1,2]
  return(output)
}

exactinformationmatrix<-function(x,del,param){

  mu    <- param[1]
  sigma <- param[2]
  output <- 0  
  n <- length(x)-1
  score_i <- rep(0,2)
  for (i in 1:n)
  {
    score_i[1] <- (log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)/sigma^2
    score_i[2] <- -1/sigma -((log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)/(sigma*sqrt(del)))*((-(log(x[i+1]/x[i]) -mu*del)/(sqrt(del)*sigma^2)) + sqrt(del)/2)
    output <- output + kronecker(t(score_i),score_i)
  }
  return(output)
}
