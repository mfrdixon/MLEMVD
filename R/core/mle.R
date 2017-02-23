#' Entry point for the maximum likelihood estimation
#'
#' 
#' @param logdensity: the model log likelihood function
#' @param x: Time series of the observed state variables
#' @param del: The uniform time step between observations
#' @param param0: The parameter vector
#' @param args: Arguments passed to the numerical optimization method 
#' @keywords maximum likelihood estimator
#' @return list containing the MLE parameter and diagnostics
#' @examples
#'  mle(ModelU1,x,0.1,c(0.2,0.3,0.1,0.5,0.9),args)
#' 
#' @export
mle<-function(logdensity,x,del,param0,args){
    
    #Set the objective function for the negative log likelihood funciton 
    #with a one dimensional argument
    objfun <- function(param){
      m <- -logdensity2loglik(logdensity,x,del,param,args)$llk
      return(m)
    }  
    
    jac <- function(param){
      jac<- grad(objfun, param)
    }
    #Perform Differential Evolution to reduce dependency of solution on the initial condition
  
    if (args$DEoptim$maxiter >0){
      DEres<-DEoptim(fn=objfun, lower=args$nloptr$l, upper=args$nloptr$u,
                   control=list(NP=args$DEoptim$population, itermax=args$DEoptim$maxiter, strategy=args$DEoptim$strategy))
      param0<-as.numeric(DEres$optim$bestmem) 
    }
   
   # minimize the negative log likelihood function
   # nloptr.print.options()  
   print('NLOPTR')
   res<- nloptr( x0=param0,
                 eval_f=objfun,
                 eval_grad_f=jac,
                 eval_g_ineq=args$nloptr$eval_g_ineq,
                 lb = args$nloptr$l,
                 ub = args$nloptr$u,
                 opts = list("algorithm"=args$nloptr$method, "maxeval" = args$nloptr$maxiter, "xtol_rel" = args$nloptr$xtol_rel, "ftol_rel"= args$nloptr$ftol_rel, "ftol_abs"=args$nloptr$ftol_abs, "print_level"=args$nloptr$print_level))


  #res<-mlsl(x0 = param0,fn = objfun, gr = NULL, lower = args$l, upper = args$u,local.method = "LBFGS", low.discrepancy = TRUE,
  #      nl.info = TRUE, control = list(ftol_abs=args$ftol_abs, ftol_rel=args$ftol_rel, xtol_rel=args$xtol_rel, maxeval=args$maxiter))
  #res<-lbfgs(x0 = param0,fn = objfun, gr = NULL, lower = args$l, upper = args$u,
  #       nl.info = TRUE, control = list(ftol_abs=args$ftol_abs, ftol_rel=args$ftol_rel, xtol_rel=args$xtol_rel, maxeval=args$maxiter))
 
 
  return(res)
}

