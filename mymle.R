library('nloptr')
library('pracma')


mymle<-function(logdensity,x_prime,del,param0,args){

    objfun <- function(param){
      # Set the objective function
      objfun <- -logdensity2loglik(logdensity,x_prime,del,param,args)
    }  
  
  # Optimize
  # nloptr.print.options()  
  
  
  res<- nloptr( x0=param0, 
                eval_f=objfun, 
                eval_g_ineq=eval_g_ineq,
                lb = args$l, 
                ub = args$u, 
                opts = list("algorithm"="NLOPT_LN_COBYLA", "maxeval" = args$maxiter, "xtol_rel" = args$eps, "print_level"=args$print_level))
  
  
  #  Compute the var-cov matrix
  #  The derivatives are computed numerically. The Hessian may not be
  #  positive definite. We report the inverse[Infomation], as well as the
  #  robust sandwich matrix.
   H <- hessian(objfun,res$solution)
   InfoMatrix <- logdensity2info(logdensity,x_prime,del,res$solution,args)
   Variance <- solve(InfoMatrix)
   invH <- solve(H)
   Variance_Robust <- invH %*% InfoMatrix %*% t(invH)
   print(InfoMatrix)
   
   res$variance <- Variance
   res$variance_robust <- Variance_Robust
   res$se <-sqrt(diag(Variance))
   res$se_robust <- sqrt(diag(Variance_Robust))
  return(res)
}

