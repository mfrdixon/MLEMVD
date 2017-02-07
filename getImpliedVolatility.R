getImpliedVolatility<-function(price, v_0, args){
  
  objfun<-function(v){
    objfun<- abs(HestonCOS(S,K_0,T_0,rate,q,sigma,kappa,theta,v,rho,args$callput)-price) 
  }
  # Infer v_0 from observed option prices using a root finding method.
  res<- nloptr( x0=v_0, 
                eval_f=objfun, 
                lb = 0.01, 
                ub = 5, 
                opts = list("algorithm"="NLOPT_LN_COBYLA", "maxeval" = 50, "xtol_rel" = args$eps, "print_level"=0))
  
  return(res$solution)
}