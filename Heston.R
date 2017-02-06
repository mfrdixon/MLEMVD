
source('logdensity2loglik.R')
source('mymle.R')
source('models/ModelB4.R')
require(nloptr)

args<-list(maxiter=500, eps=1e-8, print_level=3)
eps <-1e-8
#b,rho, kappa, theta, sigma
args$l <- c(eps, -1.0 + eps,eps,eps, eps)
args$u <- c(2.0-eps,1.0-eps,4.0-eps,1.0-eps,2.0-eps)
#b,rho, kappa, theta, sigma
eval_g_ineq <- function (x) {
  grad <- c(0,0, -2.0*x[4],-2.0*x[3],2.0*x[5])
  return(list("constraints"=c(x[5]*x[5] - 2.0*x[3]*x[4]), "jacobian"=grad))  
}
# Model = B4;
# dx1 = (a + b*x2)*dt + sqrt(x2)*(sqrt(1 - rho^2)*dW1 + rho*dW2)
# dx2 = kappa*(theta - x2)*dt + sigma*sqrt(x2)*dW2
# 5 parameters to be estimated: (b,rho,kappa,theta,sigma)

# Heston model:
# dln(S_t) = \mu dt + \sqrt{V_t}dW_t^{1}
# dV_t = \kappa(\theta - V_t)dt  + \sigma \sqrt{V_t}dW_t^{2}

  
  
# starting values for MLE algorithm and simulated series
# Step 1: Simulating Heston data

S_0      <-50 
v_0      <-0.1  # the initial variance is assumed to be unknown
rate     <-0.1  # the risk free rate is assumed to be known
q        <-0

a_0      <- rate -q
b_0      <- 1.8
rho_0    <- -0.8
kappa_0  <-3
theta_0  <-0.2
sigma_0  <-0.25 

param_0<-c(b_0,rho_0,kappa_0,theta_0,sigma_0)
 
args$mode = 'option' # calibration to ATM option prices 
args$callput = 'C'
# daily data: del = 1/52
delta <- 1/252
n     <- 50
# s_t=ln(S_t)
x1_0 <- log(S_0) #s_0 <- ln S_0;
x2_0 <- v_0
factor  <- 1
nsimul  <- factor*n
delsimul<- delta/factor

rndn1 <- rnorm(nsimul)
rndn2 <- rnorm(nsimul)

x1simul <- rep(0, nsimul)   
x2simul <- rep(0, nsimul)  

x1simul[1] <- x1_0 
x2simul[1] <- x2_0


if(args$mode=='option'){
  T_0  <- 0.1
  callput <-'C'
  x3simul <- rep(0, nsimul)  
  x3_0 <- HestonCOS(S_0,S_0,T_0,rate,q,sigma_0,kappa_0,theta_0,v_0,rho_0,'C')
  x3simul[1] <- x3_0
}


for (i in 2:nsimul){
  # dx2 = kappa*(theta - x2)*dt + sigma*sqrt(x2)*dW2
  x2simul[i] <- x2simul[i-1]+kappa_0*(theta_0 - x2simul[i-1])*delsimul + sigma_0*sqrt(x2simul[i-1])*sqrt(delsimul)*rndn2[i]
  # dx1 = (a_0 + b_0*x2)*dt + sqrt(x2)*(sqrt(1 - rho^2)*dW1 + rho*dW2)
  x1simul[i] <- x1simul[i-1]+(a_0 + b_0*x2simul[i-1])*delsimul + sqrt(x2simul[i-1])*(sqrt(1 - rho_0^2)*sqrt(delsimul)*rndn1[i] + rho_0*sqrt(delsimul)*rndn2[i])    
  
  if(args$mode=='option'){ # calibrate to option prices
   S<-exp(x1simul[i])
   v<-x2simul[i]
   x3simul[i] <- HestonCOS(S,S,T_0,rate,q,sigma_0,kappa_0,theta_0,v,rho_0,'C') 
  }
}


if (args$mode=='state'){  
  x   <- cbind(rep(0,n),rep(0,n))
  x_0  <- c(x1_0,x2_0)
  x[1,]<- x_0
  for (i in 2:n){
    x[i,1] <- x1simul[1+(i-1)*factor] #prices
    x[i,2] <- x2simul[1+(i-1)*factor] #volatilities
  }
} else if (args$mode == 'option')
{
  x   <- cbind(rep(0,n),rep(0,n),rep(0,n))
  x_0  <- c(x1_0,x2_0,x3_0)
  x[1,]<- x_0
  for (i in 2:n){
    x[i,1] <- x1simul[1+(i-1)*factor] # prices
    x[i,2] <- x2simul[1+(i-1)*factor] # volatilties
    x[i,3] <- x3simul[1+(i-1)*factor] # option prices
  }
}

# change the parameter and determine whether the calibration can converge to the correct values
param_0<-c(0.2,-0.3,2.0,0.1,0.2)

output <- mymle(ModelHeston, x, delta, param_0, args)
print(output)
