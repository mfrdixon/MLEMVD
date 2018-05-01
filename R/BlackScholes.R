#' Compute the exact log likelihood estimate of GBM
#'
#'
#' @param x Observation of the state variable at time t
#' @param del The time step between the current and previous observation
#' @param param The parameter 2-vector (mu,sigma)
#' @export output a list with a llk variable storing the result of the exact log likelihood
#' @examples
#' exactloglik(c(0.1,0.2,0.3,0.4,0.5,0.6),0.1,c(0.1,0.2))
#'

exactloglik<-function(x,del,param){

  mu    <- param[1]
  sigma <- param[2]
  n<-length(x)-1
  output <-0
  for (i in 1:n){
    output<-output+log(2*pi*del*sigma^2*x[i+1]^2)+ (log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)^2/(sigma^2*del)
  }
  output <- -0.5*output
  return(output)
}

#' Compute the exact score of GBM
#'
#'
#' @param x Observation of the state variable at time t
#' @param del The time step between the current and previous observation
#' @param param The parameter 2-vector (mu,sigma)
#' @export output a 2-list of scores using the exact likelihood estimator
#' @examples
#' exactscore(c(0.1,0.2,0.3,0.4,0.5,0.6),0.1,c(0.1,0.2))
#'
exactscore<-function(x,del,param){

  mu    <- param[1]
  sigma <- param[2]
  n<-length(x)-1
  output <-rep(0,2)
  output[2] <- -n/sigma
  for (i in 1:n){
    output[1]<-output[1] + (log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)/sigma^2
    output[2]<-output[2] -((log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)/(sigma*sqrt(del)))*((-(log(x[i+1]/x[i]) -mu*del)/(sqrt(del)*sigma^2)) + sqrt(del)/2)
  }
  return(output)
}
#' Compute the exact Hessian of GBM
#'
#'
#' @param x Observation of the state variable at time t
#' @param del The time step between the current and previous observation
#' @param param The parameter 2-vector (mu,sigma)
#' @export output a 2-by-2 hessian matrix using the exact likelihood estimator
#' @examples
#' exacthessian(c(0.1,0.2,0.3,0.4,0.5,0.6),0.1,c(0.1,0.2))
#'
exacthessian<-function(x,del,param){

  mu    <- param[1]
  sigma <- param[2]
  n<-length(x)-1
  output <-matrix(rep(0,4),2,2)
  output[2,2] <- n/(sigma^2)
  for (i in 1:n){
    output[1,1]<-output[1,1] -del/sigma^2
    output[1,2]<-output[1,2] -2*(log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)/sigma^3 + del/sigma
    g<- (log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)/(sigma*sqrt(del))
    h<- (-(log(x[i+1]/x[i]) -mu*del)/(sqrt(del)*sigma^2)) + sqrt(del)/2
    output[2,2]<-output[2,2] -h^2 +2*g*(mu*del-log(x[i+1]/x[i]))/(sqrt(del)*sigma^3)
  }
  output[2,1]<-output[1,2]
  return(output)
}
#' Compute the exact information matrix of GBM
#'
#' @param x Observation of the state variable at time t
#' @param del The time step between the current and previous observation
#' @param param The parameter 2-vector (mu,sigma)
#' @export output an exact information matrix using the exact likelihood estimator
#' @examples
#' exactinformationmatrix(c(0.1,0.2,0.3,0.4,0.5,0.6),0.1,c(0.1,0.2))
#'
exactinformationmatrix<-function(x,del,param){

  mu    <- param[1]
  sigma <- param[2]
  output <- 0
  n <- length(x)-1
  score_i <- rep(0,2)
  for (i in 1:n)
  {
    score_i[1] <- (log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)/sigma^2
    score_i[2] <- -1/sigma -((log(x[i+1]/x[i]) - (mu-sigma^2/2)*del)/(sigma*sqrt(del)))*((-(log(x[i+1]/x[i]) -mu*del)/(sqrt(del)*sigma^2)) + sqrt(del)/2)
    output <- output + kronecker(t(score_i),score_i)
  }
  return(output)
}
