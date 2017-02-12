summary<-function(res,logdensity,x_prime,del,args){

  #  Compute the var-cov matrix
  #  The derivatives are computed numerically. The Hessian may not be
  #  positive definite. We report the inverse[Infomation], as well as the
  #  robust sandwich matrix.
  objfun <- function(param){
      objfun <- -logdensity2loglik(logdensity,x_prime,del,param,args)
  }  
    
  H <- hessian(objfun,res)
  InfoMatrix <- logdensity2info(logdensity,x_prime,del,res,args)
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