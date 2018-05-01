#' Compute the information matrix = Sum_{i=1}^n Score_i' Score_i
#'
#' @param logdensity the model log likelihood function
#' @param x Time series of the observed state variables
#' @param del The uniform time step between observations
#' @param param The parameter vector
#' @param args Not currently used
#' @export output the score of the MLE
#' @examples
#'  logdensity2info(ModelU1,c(0.1,0.2,0.13,0.14),0.1,c(0.01,0.2))
#'
logdensity2info <- function(logdensity,x,del,param,args=NULL){
  y <- as.matrix(x)
  n <- dim(y)[1] - 1
  output <- 0
  lgd <- function(param){
    theta<-logdensity(x_ip1,x_i,del,param,args)$llk
  }
  for (i in 1:n)
  {
    x_ip1 <- y[i+1,]
    x_i  <- y[i,]
    score_i <- pracma::grad(lgd,param) # This is a row vector
    output <- output + kronecker(t(score_i),score_i)
  }
  return(output)
}
