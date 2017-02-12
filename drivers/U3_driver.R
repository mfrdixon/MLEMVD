# Model = linear drift, CEV diffusion (Univariate Model U3), dX = b(a - X) dt + c X^d dW
# 4 parameters to be estimated: (a,b,c,d)
#
# Works for all values of d > 1/2 using two different forms of the loglik for d=1 and d~=1 in file closedform_expansion_cev_loglik
# starting values for MLE algorithm and simulated series
source('inc.R')
source('models/univariate/ModelU3.R')

args<-list(maxiter=200, eps=1e-9, print_level=3)
eps <-1e-3
#a,b,c,d
args$l <- c(eps,eps,eps,0.5+eps)
args$u <- c(2.0-eps,2.0-eps,2.0-eps,2.0-eps)
args$eval_g_ineq <- NULL
##Step 1: Simulating CEV data

a<-0.08 
b<-0.5
c<-0.7
d<-1.5
del<-1/52 
x_0<-0.08
param_0 <-c(a,b,c,d)

# weekly data: del = 1/52
del <- 1/52

# read your data here: create a vector x where x(1) = first observations, x(n) = last observations
# instead, for illustrations purposes, let's just simulate a series from the model

n <- 5000
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
param_1 <- c(0.2,0.1,0.4,0.6)
output <- mymle(ModelU3,x,del,param_1,args)

res <- summary(output$par,ModelU3,x,del,args)
print(res$se)
print(res$se_robust)

