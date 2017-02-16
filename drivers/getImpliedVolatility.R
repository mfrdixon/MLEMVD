getImpliedVolatility<-function(price,K,T,param,v_0,args){
  
  # rho   <- param[1]
  # kappa <- param[2]
  # theta <- param[3]
  # sigma <- param[4]
  objfun<-function(v){
    objfun<- abs(HestonCOS(S,K,T,rate,q,param[4],param[2],param[3],v,param[1],args$callput, args$N)-price)
  }
  jac<-function(v){
    jac<-grad(objfun,v)
      #HestonVega(S,K,T,rate,q,param[4],param[2],param[3],v,param[1],args$callput, args$N)
  }
  # Infer v_0 from observed option prices using a root finding method.
  res<- nloptr( x0=v_0,
                eval_f=objfun,
                eval_grad_f=jac,
                lb = 0.1,
                ub = 2.0,
                opts = list("algorithm"='NLOPT_LD_LBFGS', "maxeval" = args$nloptr$maxiter, "xtol_rel" = args$nloptr$xtol_rel, "ftol_rel"= args$nloptr$ftol_rel, "ftol_abs"=args$nloptr$ftol_abs, "print_level"=0)) #args$nloptr$print_level
  
  return(res$solution)
}