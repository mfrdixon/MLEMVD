

logdensity2info <- function(logdensity,x,del,param,args){
  # Compute the information matrix = Sum_{i=1}^n Score_i' Score_i
  n <- dim(x)[1] - 1
  output <- 0
  lgd <- function(param){
    theta<-logdensity(x_ip1,x_i,del,param,args)
  }
  for (i in 1:n)
  {
    x_ip1 <- x[i+1,]
    x_i  <- x[i,]
    score_i <- grad(lgd,param) # This is a row vector
    output <- output + kronecker(t(score_i),score_i)
  }
  output <- output/n
  return(output)
}