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

  n <- dim(x)[1] - 1
  output <- 0
  for (i in 1:n){
      output <- output + logdensity(x[i+1,],x[i,],del,param, args)
  }
  #output <- output/n

  return(output)
}
