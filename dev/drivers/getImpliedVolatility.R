getImpliedVolatility<-function(price,K,T,param,v_0,args){
  
  # rho   <- param[1]
  # kappa <- param[2]
  # theta <- param[3]
  # sigma <- param[4]
  objfun<-function(v){
    objfun<- abs(HestonCOS(S,K,T,rate,q,param[4],param[2],param[3],v,param[1],args$callput, args$N)-price)
  }
  jac<-function(v){
     #jac<-grad(objfun,v)
    jac<-HestonVega(S,K,T,rate,q,param[4],param[2],param[3],v,param[1],args$callput, args$N)
  }
  # Infer v_0 from observed option prices using a root finding method.
  # res<- nloptr( x0=v_0,
  #               eval_f=objfun,
  #               eval_grad_f=jac,
  #               lb = 0.05,
  #               ub = 1.0,
  #               opts = list("algorithm"='NLOPT_LN_COBYLA', "maxeval" = 10, "xtol_rel" = 1e-3, "ftol_rel"= 1e-3, "ftol_abs"=1e-3, "print_level"=1)) #args$nloptr$print_level
  # 
  res<-list()
  res$solution <-v_0
  return(res$solution)
}