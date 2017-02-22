# Model = linear drift, CEV diffusion (Univariate Model U3), dX = b(a - X) dt + c X^d dW
# 4 parameters to be estimated: (a,b,c,d)
#
# Works for all values of d > 1/2 using two different forms of the loglik for d=1 and d~=1 in file closedform_expansion_cev_loglik
# starting values for MLE algorithm and simulated series
source('inc.R')
source('models/univariate/ModelU3.R')

args<-list()
args$nloptr<-list(maxiter=1000,
                  method = 'NLOPT_LN_COBYLA',
                  l = c(eps,eps,0.1+eps,0.5+eps),         #a,b,c,d
                  u = c(0.5-eps,1.0-eps,1.0-eps,1.1-eps), #a,b,c,d
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
##Step 1: Simulating CEV data

a<-0.08 
b<-0.5
c<-0.7
d<-0.7
del<-1/52 
x_0<-0.08
param_0 <-c(a,b,c,d)
set.seed(177)
# weekly data: del = 1/52
del <- 1/52

# read your data here: create a vector x where x(1) = first observations, x(n) = last observations
# instead, for illustrations purposes, let's just simulate a series from the model

n <- 500
x_0 <- a
factor <- 10
nsimul <- factor*n
delta <- del/factor

rndn <- rnorm(nsimul)
xsimul <- rep(nsimul,0)
xsimul[1] <- x_0
for (i in 2:nsimul){
   xsimul[i] <- xsimul[i-1] + b*(a-xsimul[i-1])*delta + c*xsimul[i-1]^d*sqrt(delta)*rndn[i]
}
x <- rep(n,0)
x[1] <- x_0
for (i in 2:n){
   x[i] <- xsimul[1+(i-1)*factor]
}

exact <- logdensity2loglik(ModelU3,x,del,param_0,args)

print(paste("Maximum log likelihood is ", exact))

#param_1 <- c(0.2,0.1,0.4,0.6)
output <- mymle(ModelU3,x,del,param_0,args)
output$solution<-c(0.0786,0.5533,0.62,0.6823)
res <- summary(output$solution,ModelU3,x,del,args)

print(res$se)
print(res$se_robust)



