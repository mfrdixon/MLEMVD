#' Entry point for the maximum likelihood estimation
#'
#' 
#' @param logdensity: is the model log likelihood function
#' @param x_prime: Observation of the state variable or the price series
#' @param del: The time step between the current and previous observation
#' @param param0: Initial set of parameters for the calibration 
#' @param args: Arguments passed to the numerical optimization method 
#' @keywords maximum likelihood estimator
#' @export
#' @examples
#' mymle(ModelHeston,x, 0.1,c(0.2,0.3,0.1,0.5,0.9),args)
#' 

mymle<-function(logdensity,x_prime,del,param0,args){

    #Set the objective function for the negative log likelihood funciton 
    #with a one dimensional argument
    objfun <- function(param){
      objfun <- -logdensity2loglik(logdensity,x_prime,del,param,args)
    }  
  
  # Optimize the negative log likelihood function
  # nloptr.print.options()  
   # res<- nloptr( x0=param0, 
   #               eval_f=objfun, 
   #               eval_g_ineq=args$eval_g_ineq,
   #               lb = args$l, 
   #               ub = args$u, 
   #               opts = list("algorithm"="NLOPT_LN_COBYLA", "maxeval" = args$maxiter, "xtol_rel" = args$eps, "print_level"=args$print_level))
   # 
   
  #res<-mlsl(x0 = param0,fn = objfun, gr = NULL, lower = args$l, upper = args$u,local.method = "LBFGS", low.discrepancy = TRUE,
  #      nl.info = TRUE, control = list(xtol_rel=args$eps, ftol_rel=args$eps, maxeval=args$maxiter))
   res<-lbfgs(x0 = param0,fn = objfun, gr = NULL, lower = args$l, upper = args$u,
         nl.info = TRUE, control = list(xtol_rel=args$eps, ftol_rel=args$eps, maxeval=args$maxiter))
   
  
  return(res)
}

