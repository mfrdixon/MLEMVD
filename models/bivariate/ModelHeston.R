ModelHeston<- function(x,x0,del,param,args){
  
  output <- 0
  rho   <- param[1]
  kappa <- param[2]
  theta <- param[3]
  sigma <- param[4]
  
  if (args$mode=='option')
  {

    S<-exp(x[1])
    K<-S
    objfun<-function(v){
      objfun<- HestonCOS(S,K,T,rate,q,param[4],param[2],param[3],v,param[1],args$callput, args$N)
    }
    
    x[2] <- getImpliedVolatility(x[3],K,T_0,param,v_0,args) # the implied volatility
    #print(paste('implied vol estimate: ',x[2]))
    # calculate the vega to obtain the Jacobian
    
    dVdv0 <- grad(objfun,x[2])
    #HestonVega(S,K,T_0,rate,q,sigma,kappa,theta,x[2],rho,args$callput,args$N)
    #print(paste('Jacobian: ',dVdv0))
    J <- dVdv0 
    if (is.nan(log(J))){
      # print('Warning: NAN occured')
      # print(paste("T: ", T_0))
      # print(paste("K: ", K))
      # print(paste("v: ", x[2]))
      # print(paste("V: ",x[3]))
      # print(paste("param: ",param))
    }
    else
      output <- -log(J)
  }

  # m,a,b,s,r
  param_prime <- c(rate, kappa*theta, kappa, sigma, rho)
  output <- output + ModelB6(x,x0,del,param_prime)
  
  return(output)
}

