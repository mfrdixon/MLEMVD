summary<-function(res,logdensity,x_prime,del,args){

  #  Compute the var-cov matrix
  #  The derivatives are computed numerically. The Hessian may not be
  #  positive definite. We report the inverse[Infomation], as well as the
  #  robust sandwich matrix.
  objfun <- function(param){
      objfun <- logdensity2loglik(logdensity,x_prime,del,param,args)
  }  
  info <- list()  
  opt<- objfun(res)
  info$score <- grad(objfun,res)
  info$H <- hessian(objfun,res)
  info$InfoMatrix <- logdensity2info(logdensity,x_prime,del,res,args)
  Variance <- solve(info$InfoMatrix)
  invH <- solve(info$H)
  Variance_Robust <- invH %*% info$InfoMatrix %*% t(invH)
  
  
  info$variance <- Variance
  info$variance_robust <- Variance_Robust
  info$se <-sqrt(diag(Variance))
  info$se_robust <- sqrt(diag(Variance_Robust))
  
  if (args$plot==TRUE){
    pdf('MLE.pdf')
    n.param <- length(res)
    par(mfrow=c(ceil(n.param/2),2))
    
    eps.grid <- 0.01
    for (j in 1:n.param){
     grid <- seq(args$nloptr$l[j], args$nloptr$u[j],eps.grid)
     n.grid <- length(grid)
     llk  <- rep(0,n.grid)
     for (i in 1:n.grid){
       res_prime  <- res
       res_prime[j] <- grid[i]
       #print(res_prime)
       llk[i] <- objfun(res_prime)
     }
     plot(grid,llk, xlab=paste('param ',j), ylab="log likelihood", type='l')
     points(x=res[j],y=opt, type='p', pch=19, col='red', cex=1)
    }
    dev.off()
  }
  
  return(info)
}