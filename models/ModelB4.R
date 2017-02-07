

ModelHeston<- function(x,x0,del,param,args){
  
  output <- 0
  if (args$mode=='option')
  {
    e<- 1e-5
    S<-exp(x[1])
    rho   <- param[1]
    kappa <- param[2]
    theta <- param[3]
    sigma <- param[4]
    
    
    objfun<-function(vega){
      objfun<- abs(HestonCOS(S,S,T_0,rate,q,sigma,kappa,theta,vega,rho,args$callput)-x[3]) 
    }
    # Infer v_0 from observed option prices using a root finding method.
    res<- nloptr( x0=v_0, 
                  eval_f=objfun, 
                  lb = c(0.01), 
                  ub = c(5), 
                  opts = list("algorithm"="NLOPT_LN_COBYLA", "maxeval" = 50, "xtol_rel" = args$eps, "print_level"=0))
    
    x[2] < - res$solution # the implied volatility
   
    # calculate the vega to obtain the Jacobian
    dVdv0 <- HestonCOS_vega(S,S,T_0,rate,q,sigma,kappa,theta,x[2],rho,args$callput)
    
    J <- dVdv0 
    if (is.nan(log(J))){
      print('Warning: NAN occured')
      print(J)
      print(x[2])
    }
    else
      output <- -log(J)
  }
  
  
  a_0      <- rate -q   # risk free rate - annualized dividend yield
  a_1      <- -0.5      # b
  a        <- param[3]  # theta
  b        <- param[2]  # kappa
  L2       <- 0
  g        <- param[4]  # sigma
  beta     <- 0.5
  f        <- 1
  r        <- param[1]  # rho
  
  param_prime <- c(a_0,a_1,a,b,L2,g,beta,f,r)
  output <- output + ModelB4(x,x0,del,param_prime)
  
  return(output)
}

ModelB4 <- function(x,x0,del,param){

  # Model:
  # mu1[x1,x2] <- a0+a1*x2 
  # mu2[x1,x2] <- b*(a-x2)+L2*g*x2^beta*sqrt(a+f*(x2-a)) 
  # sigma11[x1,x2] <- sqrt(1-r^2)*sqrt(a+f*(x2-a)) 
  # sigma22[x1,x2] <- g*x2^beta 
  # sigma12[x1,x2] <- r*sqrt(a+f*(x2-a)) 
  # sigma21[x1,x2] <- 0 
  
  x1 <- x[1] 
  x2 <- x[2] 
  x10 <- x0[1] 
  x20 <- x0[2] 
    
  a0 <- param[1]  
  a1 <- param[2] 
  a <- param[3] 
  b <- param[4] 
  L2 <- param[5] 
  g <- param[6] 
  beta <- param[7] 
  f <- param[8] 
  r <- param[9] 
  
  Dv <- log(g)+(1/2)*log(1-r^2)+beta*log(x2)+(1/2)*log(a+f*(x2-a)) 
  cm1 <- -((beta*(x2 - x20)^3*x20^(-1 - 2*beta))/(2*g^2*(-1 + r^2))) - (f*(x1 - x10)^2*(x2 - x20))/(4*(-1 + r^2)*(a*(-1 + f) - f*x20)^2) + (f^2*g^2*(x1 - x10)^4*x20^(2*beta))/ (96*(-1 + r^2)^2*(a - a*f + f*x20)^4) - (f^2*g*r*(x1 - x10)^3*(x2 - x20)*x20^beta)/(24*(-1 + r^2)^2*(a - a*f + f*x20)^(7/2)) + (r*(x1 - x10)*(x2 - x20)^2*x20^(-1 - beta)*(2*a*beta*(-1 + f) - (1 + 2*beta)*f*x20))/(4*g*(-1 + r^2)*(a*(-1 + f) - f*x20)* sqrt(a - a*f + f*x20)) - (f*(x1 - x10)^2*(x2 - x20)^2*(2*a*beta*(-1 + f)*(-1 + r^2) + f*(6 - 9*r^2 - 2*beta*(-1 + r^2))*x20))/ (48*(-1 + r^2)^2*x20*(a - a*f + f*x20)^3) + (r*(x1 - x10)*(x2 - x20)^3*x20^(-2 - beta)*(-4*a^2*beta*(1 + beta)*(-1 + f)^2*(-1 + r^2) + 4*a*beta*(3 + 2*beta)*(-1 + f)*f*(-1 + r^2)*x20 + f^2*(2 - 3*r^2 - 8*beta*(-1 + r^2) - 4*beta^2*(-1 + r^2))*x20^2))/ (24*g*(-1 + r^2)^2*(a*(-1 + f) - f*x20)^2*sqrt(a - a*f + f*x20)) + ((x2 - x20)^4*(4*a^2*beta*(4 + 7*beta)*(-1 + f)^2*(-1 + r^2) - 8*a*beta*(4 + 7*beta)*(-1 + f)*f*(-1 + r^2)*x20 + f^2*(r^2 + 16*beta*(-1 + r^2) + 28*beta^2*(-1 + r^2))*x20^2))/x20^(2*(1 + beta))/(96*g^2*(-1 + r^2)^2*(a - a*f + f*x20)^2) + (a*(-1 + f)*(x2 - x20)^2 - f*(x2 - x20)^2*x20 - g*(x1 - x10)*x20^beta*(g*(x1 - x10)*x20^beta + 2*r*(-x2 + x20)*sqrt(a - a*f + f*x20)))/ x20^(2*beta)/(2*g^2*(-1 + r^2)*(a*(-1 + f) - f*x20)) 
  c0 <- (f*g^2*(x1 - x10)^2*x20^(-1 + 2*beta)*(2*a*beta*(-1 + f) + (3 - 2*beta)*f*x20))/ (48*(-1 + r^2)*(a - a*f + f*x20)^3) + ((x1 - x10)*(x2 - x20)*x20^(-2 - beta)*(-24*a^3*b*beta*(-1 + f)^2*r*x20 + 12*a^2*(-1 + f)*r*(b*(2 - f + beta*(-2 + 6*f))*x20^2 - (-1 + beta)*beta*(-1 + f)*g^2*x20^(2*beta)) + f*x20^2*(12*b*(-1 + 2*beta)*f*r*x20^2 - g*x20^beta*((15 - 28*beta + 12*beta^2)*f*g*r*x20^beta - 24*a0*sqrt(a - a*f + f*x20))) - 4*a*x20*(3*b*f*(3 - 2*f + beta*(-4 + 6*f))*r*x20^2 - 2*(-1 + f)*g*x20^beta*(-5*beta*f*g*r*x20^beta + 3*beta^2*f*g*r*x20^beta + 3*a1*x20*sqrt(a - a*f + f*x20)))))/(48*g*(-1 + r^2)*(a*(-1 + f) - f*x20)^2*sqrt(a - a*f + f*x20)) + (1/(4*g^2*(-1 + r^2)*(a - a*f + f*x20)))*(x2 - x20)*x20^(-1 - 2*beta)*(4*a^2*b*(-1 + f)*x20 + a*(4*b*(1 - 2*f)*x20^2 - 2*(-1 + f)*g*x20^beta*(beta*g*x20^beta - 2*L2*x20*sqrt(a - a*f + f*x20))) + x20*(4*b*f*x20^2 + g*x20^beta*(4*r*(a0 + a1*x20)*sqrt(a - a*f + f*x20) + f*((-1 + 2*beta)*g*x20^beta - 4*L2*x20*sqrt(a - a*f + f*x20))))) + (1/(4*g*(-1 + r^2)*(a*(-1 + f) - f*x20)*sqrt(a - a*f + f*x20)))*(x1 - x10)*x20^(-1 - beta)* (4*a^2*b*(-1 + f)*r*x20 - 2*a*r*(2*b*(-1 + 2*f)*x20^2 + (-1 + f)*g*x20^beta*(beta*g*x20^beta - 2*L2*x20*sqrt(a - a*f + f*x20))) + x20*(4*b*f*r*x20^2 + g*x20^beta*(4*(a0 + a1*x20)*sqrt(a - a*f + f*x20) + f*r*((-1 + 2*beta)*g*x20^beta - 4*L2*x20*sqrt(a - a*f + f*x20))))) + ((1/(48*g^2*(-1 + r^2)*(a - a*f + f*x20)^2))*(x2 - x20)^2* (48*a^3*b*beta*(-1 + f)^2*x20 - 12*a^2*(-1 + f)*(2*b*(1 - f + beta*(-2 + 6*f))*x20^2 + beta*(-1 + f)*g*x20^beta* (g*x20^beta - 2*L2*x20*sqrt(a - a*f + f*x20))) - f*x20^2*(24*b*(-1 + 2*beta)*f*x20^2 + g*x20^beta*(12*r*(a0 + 2*a0*beta + a1*(-1 + 2*beta)*x20)*sqrt(a - a*f + f*x20) + f*((-9 + 14*beta)*g*x20^beta + 12*(1 - 2*beta)*L2*x20*sqrt(a - a*f + f*x20)))) + 2*a*x20*(24*b*f*(1 - f + beta*(-2 + 3*f))*x20^2 + (-1 + f)*g*x20^beta*(6*(f*L2 - 2*a1*r)*x20*sqrt(a - a*f + f*x20) + beta*(12*r*(a0 + a1*x20)*sqrt(a - a*f + f*x20) + f*(13*g*x20^beta - 24*L2*x20*sqrt(a - a*f + f*x20)))))))/x20^(2*(1 + beta)) 
  c1 <- ((1/(96*g^2*(-1 + r^2)*(a - a*f + f*x20)^2))*(48*a^4*b^2*(-1 + f)^2*x20^2 - 48*a^3*(-1 + f)*x20*(2*b^2*(-1 + 2*f)*x20^2 + (-1 + f)^2*g^2*L2^2*x20^(1 + 2*beta) - b*(-1 + f)*g*x20^beta* (beta*g*(-2 + r^2)*x20^beta + 2*L2*x20*sqrt(a - a*f + f*x20))) + 12*a^2*(4*b^2*(1 - 6*f + 6*f^2)*x20^4 - (-1 + f)^2*g^2*x20^(2*beta)*(8*a0*L2*r*x20^2 - 12*f*L2^2*x20^3 + 8*a1*L2*r*x20^3 - 2*beta*g^2*x20^(2*beta) + beta^2*g^2*x20^(2*beta) + 2*beta*g^2*r^2*x20^(2*beta) - 2*beta^2*g^2*r^2*x20^(2*beta) + 4*beta*g*L2*x20^(1 + beta)*sqrt(a - a*f + f*x20)) - 2*b*(-1 + f)*g*x20^(2 + beta)* (g*(-2 + 2*f + 2*r^2 - f*r^2 + 2*beta*(-1 + 3*f)*(-2 + r^2))*x20^beta + 4*sqrt(a - a*f + f*x20)* ((-a0)*r - (L2 - 3*f*L2 + a1*r)*x20))) + f*x20^2*(48*b^2*f*x20^4 - 24*b*g*x20^(2 + beta)* (-4*r*(a0 + a1*x20)*sqrt(a - a*f + f*x20) + f*((-1 + 2*beta)*g*(-2 + r^2)*x20^beta + 4*L2*x20*sqrt(a - a*f + f*x20))) + g^2*x20^(2*beta)*(48*a0^2*x20 + 48*a1^2*x20^3 - 24*a1*r*x20*(4*f*L2*x20^2 - (-1 + 2*beta)*g*x20^beta*sqrt(a - a*f + f*x20)) + f*(48*f*L2^2*x20^3 + g*x20^beta*(g*(-3 + 6*r^2 + beta*(28 - 40*r^2) + 12*beta^2*(-1 + 2*r^2))*x20^beta + 24*(1 - 2*beta)*L2*x20*sqrt(a - a*f + f*x20))) + 24*a0*(4*a1*x20^2 - r*(4*f*L2*x20^2 - (-1 + 2*beta)*g*x20^beta*sqrt( a - a*f + f*x20))))) - 4*a*x20*(24*b^2*f*(-1 + 2*f)*x20^4 - 6*b*g*x20^(2 + beta)*(4*r*(a0 + a1*x20)*sqrt(a - a*f + f*x20) + 2*f^2*((-1 + 3*beta)*g*(-2 + r^2)*x20^beta + 6*L2*x20*sqrt(a - a*f + f*x20)) - f*(g*(4 - 3*r^2 + 4*beta*(-2 + r^2))*x20^beta + 8*sqrt(a - a*f + f*x20)*(a0*r + (L2 + a1*r)*x20))) + (-1 + f)*g^2*x20^(2*beta)*(12*a0^2*x20 + 12*a1^2*x20^3 - 12*a1*r*x20*(4*f*L2*x20^2 - beta*g*x20^beta*sqrt(a - a*f + f*x20)) + 12*a0*(2*a1*x20^2 - 4*f*L2*r*x20^2 + beta*g*r*x20^beta*sqrt(a - a*f + f*x20)) + f*(36*f*L2^2*x20^3 + g*x20^beta*(6*beta^2*g*(-1 + 2*r^2)*x20^beta + 6*L2*x20*sqrt(a - a*f + f*x20) + beta*((-g)*(-13 + 16*r^2)*x20^beta - 24*L2*x20*sqrt(a - a*f + f*x20))))))))/x20^(2*(1 + beta)) 
  output <- -log(2*pi*del) - Dv + cm1/del + c0 + c1*del 
  return(output)
}
