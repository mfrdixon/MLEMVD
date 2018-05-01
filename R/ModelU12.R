#' Compute the maximum likelihood estimate of Model U12
#'
#' @param x Observation of the state variable at time t
#' @param x0 Observation of the state variable at time t-1
#' @param del The time step between the current and previous observation
#' @param param The parameter 3-vector (alpha,beta,gamma)
#' @export output a list with a llk variable storing the result of the log likelihood calculation
#' @examples
#' ModelU12(0.4,0.3,0.1,c(0.1,0.3,0.2))
#'

ModelU12 <- function(x,x0,del,param)
{
  beta  <- param[1]
  alpha <- param[2]
  gamma <- param[3]

  s_x <- gamma * x ^ (1/2)

  y <- (2*sqrt(x))/gamma
  y0 <- (2*sqrt(x0))/gamma

  cm1 <- (-(1/2))*(y - y0)^2
  c0 <- -((4*beta)/(gamma^4*y^2)) - (1/192)*alpha*gamma^4*y^6 + (4*beta)/(gamma^4*y0^2) + (1/192)*alpha*gamma^4*y0^6 - (1/2)*log(y/y0)
  c1 <- (1/12)*alpha*beta*y^2 + (1/80)*alpha*gamma^4*y^4 - (alpha^2*gamma^8*y^10)/22528 - (32*beta^2)/(5*gamma^8*y*y0^5) - (32*beta^2)/(5*gamma^8*y^2*y0^4) - (32*beta^2)/(5*gamma^8*y^3*y0^3) + (16*beta)/(3*gamma^4*y*y0^3) - (32*beta^2)/(5*gamma^8*y^4*y0^2) + (16*beta)/(3*gamma^4*y^2*y0^2) - (32*beta^2)/(5*gamma^8*y^5*y0) + (16*beta)/(3*gamma^4*y^3*y0) - 3/(8*y*y0) + (1/12)*alpha*beta*y*y0 + (1/80)*alpha*gamma^4*y^3*y0 - (alpha^2*gamma^8*y^9*y0)/22528 + (1/12)*alpha*beta*y0^2 + (1/80)*alpha*gamma^4*y^2*y0^2 - (alpha^2*gamma^8*y^8*y0^2)/22528 + (1/80)*alpha*gamma^4*y*y0^3 - (alpha^2*gamma^8*y^7*y0^3)/22528 + (1/80)*alpha*gamma^4*y0^4 - (alpha^2*gamma^8*y^6*y0^4)/22528 - (alpha^2*gamma^8*y^5*y0^5)/22528 - (alpha^2*gamma^8*y^4*y0^6)/22528 - (alpha^2*gamma^8*y^3*y0^7)/22528 - (alpha^2*gamma^8*y^2*y0^8)/22528 - (alpha^2*gamma^8*y*y0^9)/22528 - (alpha^2*gamma^8*y0^10)/22528
  c2 <- (alpha*beta)/12 + (3/80)*alpha*gamma^4*y^2 - (9*alpha^2*gamma^8*y^8)/22528 - (32*beta^2)/(gamma^8*y^2*y0^6) - (256*beta^2)/(5*gamma^8*y^3*y0^5) - (288*beta^2)/(5*gamma^8*y^4*y0^4) + (16*beta)/(gamma^4*y^2*y0^4) - (256*beta^2)/(5*gamma^8*y^5*y0^3) + (64*beta)/(3*gamma^4*y^3*y0^3) - (32*beta^2)/(gamma^8*y^6*y0^2) + (16*beta)/(gamma^4*y^4*y0^2) - 3/(8*y^2*y0^2) + (1/20)*alpha*gamma^4*y*y0 - (alpha^2*gamma^8*y^7*y0)/1408 + (3/80)*alpha*gamma^4*y0^2 - (21*alpha^2*gamma^8*y^6*y0^2)/22528 - (3*alpha^2*gamma^8*y^5*y0^3)/2816 - (25*alpha^2*gamma^8*y^4*y0^4)/22528 - (3*alpha^2*gamma^8*y^3*y0^5)/2816 - (21*alpha^2*gamma^8*y^2*y0^6)/22528 - (alpha^2*gamma^8*y*y0^7)/1408 - (9*alpha^2*gamma^8*y0^8)/22528
  output <- list()
  output$llk <- (-(1/2))*log(2*pi*del) - log(s_x) + cm1/del + c0 + c1 * del + c2 * (del^2/2)

  return(output)
}


# ModelU12(3,4,1/52,c(0.1,0.2,0.3))
