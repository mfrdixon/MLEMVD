# Author: Matthew Dixon
# Email:  mdixon7@stuart.iit.edu
# Date:   2/22/2017
# Notes:  Adapted from the Matlab library 'ClosedFormMLE'
#         https://www.princeton.edu/~yacine/ClosedFormMLEUsersGuide.pdf

logdensity2loglik<-function(logdensity,x,del,param,args){
# Inputs:
#
# logdensity: is the function handle of the transition density with the
# structure such as logdensity(x,x0,del,param)
#
# x: the price series nxk, n is the number of observations and k is the
# dimension of the multivariate series
#
# del: sampling interval
# param: parameter vector
  
# example input of x 
# x<-cbind(c(1,2,3), c(3,4,5))  
  y <- as.matrix(x)
  n <- dim(y)[1] - 1
  output <- list()
  output$llk <- 0
  output$v <- rep(0,n+1)
  output$v[1] <-v_0
  for (i in 1:n){
      m <- logdensity(y[i+1,],y[i,],del,param,args)
      output$llk <- output$llk + m$llk  # update the log likelihood function
      output$v[i+1] <- m$v              # update the implied volatility
  }

  return(output)
}
