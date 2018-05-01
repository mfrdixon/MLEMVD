#' Compute the maximum likelihood estimate of Model U2
#'
#'
#' @param x Observation of the state variable at time t
#' @param x0 Observation of the state variable at time t-1
#' @param del The time step between the current and previous observation
#' @param param The parameter 3-vector
#' @param args Not currently used
#' @export output a list with a llk variable storing the result of the log likelihood calculation
#' @examples
#' ModelU2(0.4,0.3,0.1,c(0.1,0.2,0.3))
#'
ModelU2 <- function(x,x0,del,param,args=NULL)
{
  if (length(param)==3){
   a <- param[1]
   b <- param[2]
   d <- param[3]
  }
  else if (length(param)==2){
   a <- 0
   b <- param[1]
   d <- param[2]
  }
  y <- log(x)/d
  y0 <- log(x0)/d


  E <- exp(1)
  sx <- d*x

  cYm1 <- (-(1/2))*(y - y0)^2
  cY0 <- (E^((-d)*y) - E^((-d)*y0))*(-(a/d^2)) + (y - y0)*(b/d - d/2)

  if (y != y0)
  {
    cY1 <- (a^2/(4*d^3))*((E^(-2*d*y) - E^(-2*d*y0))/(y - y0)) + ((a*b)/d^3 - a/d)*((E^((-d)*y) - E^((-d)*y0))/(y - y0)) -  (2*b - d^2)^2/(8*d^2)
    cY2 <- (-(a^2/(2*d^3)))*((E^(-2*d*y) - E^(-2*d*y0))/(y - y0)^3) +
      ((2*a)/d - (2*a*b)/d^3)*((E^((-d)*y) - E^((-d)*y0))/(y - y0)^3) + (-(a^2/(2*d^2)))*((E^(-2*d*y) + E^(-2*d*y0))/(y - y0)^2) +
      (a - (a*b)/d^2)*((E^((-d)*y) + E^((-d)*y0))/(y - y0)^2)
  }
  else
  {
    cY1 <- (-4*a^2 - 8*a*(b - d^2)*E^(d*y) - (-2*b + d^2)^2*E^(2*d*y))/(E^(2*d*y)*(8*d^2))
    cY2 <- ((1/6)*a*(-2*a + (-b + d^2)*E^(d*y)))/E^(2*d*y)
  }

  output <- list()
  output$llk <- (-(1/2))*log(2*pi*del) - log(sx) + cYm1/del + cY0 +  cY1*del + cY2*(del^2/2)

  return(output)
}
