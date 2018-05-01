#' An example of using maximum likelihood estimation to fit the Heston model to simulated stock and option prices
#'
#' Heston model:
#' \eqn{dln(S_t) = (\mu - \sigma^2/2) dt + \sqrt{V_t}dW_t^{1}}
#' \eqn{dV_t = \kappa(\theta - V_t)dt  + \sigma \sqrt{V_t}dW_t^{2}}
#' 4 parameters to be estimated: (rho, kappa, theta, sigma)
#'
#' Works for all values of sigma,kappa,theta > 0
#' @keywords Heston model

source('R/inc.R')
#source('ModelB6.R')
#source('ModelHeston.R')

set.seed(20)
eps<-1e-3
#rho, kappa, theta, sigma
eval_g_ineq <- function (x) {
  grad <- c(0,-2.0*x[3],-2.0*x[2],2.0*x[4])
  return(list("constraints"=c(x[4]*x[4] - 2.0*x[2]*x[3]), "jacobian"=grad))
}

eval_g0 <- function (x) {
  return(x[4]*x[4] - 2.0*x[2]*x[3])
}

eval_jac_g_ineq <- function(x){
  return(c(0,-2.0*x[3],-2.0*x[2],2.0*x[4]))
}


args<-list()
args$N <- 128
args$plot <- TRUE
args$nloptr<-list(maxiter=500,
                  method = "NLOPT_LD_MMA", #'NLOPT_LD_TNEWTON', 'NLOPT_LN_COBYLA'
                  #rho, kappa, theta, sigma
                  l = c(-1.0 + eps,eps, 0.1 + eps, eps),
                  u = c(0.0 - eps,4.0-eps,0.5-eps,0.5-eps),
                  eval_g_ineq = eval_g0,
                  eval_jac_g_ineq = eval_jac_g_ineq,
                  check_derivatives = TRUE,
                  check_derivatives_print="all",
                  ftol_abs=1e-16,
                  xtol_rel=1e-16,
                  ftol_rel=1e-16,
                  print_level=3)

args$DEoptim$maxiter    <- 0
args$DEoptim$population <- 100
args$DEoptim$strategy   <- 5




## Step 1: Simulating Heston data

# starting values for MLE algorithm and simulated series
S_0       <-100
args$v_0  <-0.1  # the initial variance is assumed to be unknown
args$rate <-0.1  # the risk free rate is assumed to be known
args$q    <-0
mu_0      <-args$rate-args$q

rho_0    <- -0.5
kappa_0  <- 3.0
theta_0  <- 0.3
sigma_0  <- 0.15

param_0<-c(rho_0,kappa_0,theta_0,sigma_0)

#set.seed(99)

args$mode <- 'direct' # implied = calibration to ATM option prices
args$callput <- 'C'
delta <- 1/252 # daily data
n <- 500
n.burnin <- 500
factor <- 10
nsimul <- factor*(n+n.burnin)
delsimul<- delta/factor


# s_t=ln(S_t)
x1_0 <- log(S_0) #s_0 <- ln S_0;
x2_0 <- args$v_0


rndn1 <- rnorm(nsimul)
rndn2 <- rnorm(nsimul)

x1simul <- rep(0, nsimul)
x2simul <- rep(0, nsimul)

x1simul[1] <- x1_0
x2simul[1] <- x2_0


if(args$mode=='implied'){
  args$T_0  <- 0.1
  K_0  <- S_0
  callput <-'C'
  x3simul <- rep(0, nsimul)
  x3_0 <- HestonCOS(S_0,K_0,args$T_0,args$rate,args$q,sigma_0,kappa_0,theta_0,args$v_0,rho_0,'C',args$N)
  x3simul[1] <- x3_0
}


for (i in 2:nsimul){
  #dx2 = kappa*(theta - x2)*dt + sigma*sqrt(x2)*dW2
  #x2simul[i] <- x2simul[i-1]+kappa_0*(theta_0 - x2simul[i-1])*delsimul + sigma_0*sqrt(x2simul[i-1])*sqrt(delsimul)*rndn2[i]
  #dx1 = (a_0 + b_0*x2)*dt + sqrt(x2)*(sqrt(1 - rho^2)*dW1 + rho*dW2)
  #x1simul[i] <- x1simul[i-1]+(mu_0 -0.5*x2simul[i-1])*delsimul + sqrt(x2simul[i-1])*(sqrt(1 - rho_0^2)*sqrt(delsimul)*rndn1[i] + rho_0*sqrt(delsimul)*rndn2[i])

  # dx2 = kappa*(theta - x2)*dt + sigma*sqrt(x2)*(r*dW1 + sqrt(1-r^2)*dW2)
   x2simul[i] <- x2simul[i-1]+kappa_0*(theta_0 - x2simul[i-1])*delsimul + sigma_0*sqrt(x2simul[i-1]*delsimul)*(rho_0*rndn1[i] + sqrt(1-rho_0^2)*rndn2[i])
  # dx1 = (a_0 + b_0*x2)*dt + sqrt(x2)*dW1
   x1simul[i] <- x1simul[i-1]+(mu_0 -0.5*x2simul[i-1])*delsimul + sqrt(x2simul[i-1]*delsimul)*rndn1[i]



  if(args$mode=='implied'){ # calibrate to option prices
   S<-exp(x1simul[i])
   v<-x2simul[i]
   K<-S
   x3simul[i] <- HestonCOS(S,K,args$T_0,args$rate,args$q,sigma_0,kappa_0,theta_0,v,rho_0,'C',args$N)
  }
}



if (args$mode=='direct'){
  x   <- cbind(rep(0,n),rep(0,n))
  #x_0  <- c(x1_0,x2_0)
  #x[1,]<- x_0
  for (i in 1:n){
    x[i,1] <- x1simul[n.burnin*factor+1+(i-1)*factor] #prices
    x[i,2] <- x2simul[n.burnin*factor+1+(i-1)*factor] #volatilities
  }
}else if (args$mode == 'implied')
{
  x   <- cbind(rep(0,n),rep(0,n),rep(0,n))
  #x_0  <- c(x1_0,x2_0,x3_0)
  #x[1,]<- x_0
  for (i in 1:n){
    x[i,1] <- x1simul[n.burnin*factor+1+(i-1)*factor] # prices
    x[i,2] <- x2simul[n.burnin*factor+1+(i-1)*factor] # volatilties
    x[i,3] <- x3simul[n.burnin*factor+1+(i-1)*factor] # option prices
  }
}

# change the parameter and determine whether the calibration can converge to the correct values
#param_0<-c(-0.3,2.0,0.1,0.2)

exact <- logdensity2loglik(ModelHeston,x,delta,param_0,args)

print(paste("Maximum log likelihood is ", exact))

output <- mle(ModelHeston, x, delta, param_0, args)

rho.est    <- output$solution[1]
kappa.est  <- output$solution[2]
theta.est  <- output$solution[3]
sigma.est  <- output$solution[4]

v <- logdensity2loglik(ModelHeston,x,delta,output$solution,args)$v
res <- summary(ModelHeston,x,delta,output$solution,args)



print(paste("Standard Error Estimate ", res$se))
print(paste("Huber Sandwich Error Estimate ", res$se_robust))

require("quantmod")
require("xts")


st<-"2016-06-15"
end<-"2016-12-31"
v[1]<- getImpliedVolatility(S_0,x[1,3],S_0,args$T_0,args$rate,args$q,sigma.est,kappa.est,theta.est,args$v_0,rho.est,args$callput, 0.01,100, args$N) # the implied volatility

dates<-seq(as.Date(st), as.Date(end), "days")
vol.sim <-as.xts(x[,2],order.by=dates,frequency = NULL)
vol.model<-as.xts(v,order.by=dates,frequency = NULL)

chartSeries(vol.sim-vol.model)

price.sim <-as.xts(x[,3],order.by=dates,frequency = NULL)
prices <- rep(0,n)

for (i in 1:n){
  prices[i]<-HestonCOS(exp(x[i,1]),exp(x[i,1]),args$T_0,args$rate,args$q,sigma.est,kappa.est,theta.est,v[i],rho.est,'C',args$N)
}
price.model<-as.xts(prices,order.by=dates,frequency = NULL)

chartSeries(price.model-price.sim)
#addTA(price.sim, on=1, col='red')

