
objfun<-function(v){
  objfun<- HestonCOS(S,K,T_0,rate,q,sigma_0,kappa_0,theta_0,v,rho_0,'C',args$N)
  return(objfun)
}

st<-proc.time()
logdensity2loglik(ModelHeston,x,delta,param_0,args)
end<-proc.time() -st

st<-proc.time()
vega<-HestonCOS_vega(S,K,T_0,rate,q,sigma_0,kappa_0,theta_0,v_0,rho_0,'C',args$N)
end<-proc.time() -st
print(vega)
print(end)
st<-proc.time()
vega<-HestonVega(S,K,T_0,rate,q,sigma_0,kappa_0,theta_0,v_0,rho_0,'C',args$N)
end<-proc.time() -st
print(vega)
print(end)

# Black-Scholes Option Value
# Call value is returned in values[1], vega in values[2]
blackscholes <- function(S, X, rf, T, sigma) {
  values <- c(2)
  
  d1 <- (log(S/X)+(rf+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  values[1] <- S*pnorm(d1) - X*exp(-rf*T)*pnorm(d2) # call
  values[2] <- S*dnorm(d1)*sqrt(T)#vega
  
  values
}

sigma_0 <- 0.00001
kappa_0 <- 0.00001
theta_0 <- v_0
rho_0 <-0
T_0 <- 1
print("Limiting condition \sigma->0, \kappa->0, v_0=\theta->\sigma_{BS}^2")
HC<-HestonCOS(S,K,T_0,rate,q,sigma_0,kappa_0,theta_0,v_0,rho_0,'C',args$N)
HCv<-HestonVega(S,K,T_0,rate,q,sigma_0,kappa_0,theta_0,v_0,rho_0,'C',args$N)
print(paste('Heston Call:', HC))
print(paste('Heston vega:', HCv))
BS<-blackscholes(S,K,rate,T_0,sqrt(v_0))
print(paste('BS Call:', BS[1]))
print(paste('BS vega:', BS[2]))


