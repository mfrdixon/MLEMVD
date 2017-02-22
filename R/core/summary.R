#' Compute the summary information of the maximum likelihood estimate 
#'
#' @param logdensity: the model log likelihood function
#' @param x: Time series of the observed state variables
#' @param del: The uniform time step between observations
#' @param param: The maximum likelihood estimated parameter
#' @export output: the summary information on the MLE
#' @examples
#'  summary(ModelU1,x,0.1,c(0.01,0.2))
#'
summary<-function(logdensity,x,del,param0,args=NULL){

  #  Compute the var-cov matrix
  #  The derivatives are computed numerically. The Hessian may not be
  #  positive definite. We report the inverse[Infomation], as well as the
  #  robust sandwich matrix.
  objfun <- function(param){
      objfun <- logdensity2loglik(logdensity,x,del,param,args)$llk
  }  
  info <- list()  
  opt<- objfun(param0)
  info$score <- grad(objfun,param0)
  info$H <- hessian(objfun,param0)
  info$InfoMatrix <- logdensity2info(logdensity,x,del,param0,args)
  Variance <- solve(info$InfoMatrix)
  invH <- solve(info$H)
  Variance_Robust <- invH %*% info$InfoMatrix %*% t(invH)
  
  
  info$variance <- Variance
  info$variance_robust <- Variance_Robust
  info$se <-sqrt(diag(Variance))
  info$se_robust <- sqrt(diag(Variance_Robust))
  
  if (args$plot==TRUE){
    pdf('MLE.pdf')
    n.param <- length(param0)
    par(mfrow=c(ceil(n.param/2),2))
    
    eps.grid <- 0.01
    for (j in 1:n.param){
     grid <- seq(args$nloptr$l[j], args$nloptr$u[j],eps.grid)
     n.grid <- length(grid)
     llk  <- rep(0,n.grid)
     for (i in 1:n.grid){
       param_prime  <- param0
       param_prime[j] <- grid[i]
       llk[i] <- objfun(param_prime)
     }
     plot(grid,llk, xlab=paste('param ',j), ylab="log likelihood", type='l')
     points(x=param0[j],y=opt, type='p', pch=19, col='red', cex=1)
    }
    dev.off()
  }
  
  return(info)
}
