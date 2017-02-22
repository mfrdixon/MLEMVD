#' Compute the maximum likelihood estimate of Model U1
#'
#' 
#' @param x: Observation of the state variable at time t 
#' @param x0: Observation of the state variable at time t-1 
#' @param del: The time step between the current and previous observation
#' @param param: The parameter 3-vector 
#' @export output: a list with a llk variable storing the result of the log likelihood calculation
#' @examples
#' ModelU1(0.4,0.3,0.1,c(0.1,0.2,0.3))
#'
ModelU1 <- function(x,x0,del,param,args=NULL){
  a <- param[1] 
  b <- param[2] 
  sigma <- param[3] 
  
  y <- 2/(sqrt(x)*sigma) 
  y0 <- 2/(sqrt(x0)*sigma) 
  sx <- sigma*x^(3/2) 
  output <- list()  
  output$llk <- -log(2*pi*del)/2 - log(sx)   
  + (-(1/2))*(y - y0)^2 / del   
  + (a*(-y^2 + y0^2)* sigma^2 + (-8*b + 6* sigma^2)*log(y/y0))/(4* sigma^2)  
  + ( -((1/(24*y*y0* sigma^4))* (48*b^2 + 24*b*(-2 + a*y*y0)* sigma^2 + (9 - 24*a*y*y0 + a^2*y*y0*(y^2 + y*y0 + y0^2))* sigma^4)) )*del   
  + ( -((48*b^2 - 48*b* sigma^2 + (9 + a^2*y^2*y0^2)* sigma^4)/ (24*y^2*y0^2* sigma^4)) )*(del^2/2) 
  
  return(output) 
}

