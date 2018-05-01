#' Compute the maximum likelihood estimate of Model U4
#'
#' CIR model: dX = kappa (alpha - X) dt + sigma sqrt(X) dW
#' @param x Observation of the state variable at time t
#' @param x0 Observation of the state variable at time t-1
#' @param del The time step between the current and previous observation
#' @param param The parameter 3-vector (kappa,alpha,sigma)
#' @export output a list with a llk variable storing the result of the log likelihood calculation
#' @examples
#' ModelU4(0.4,0.3,0.1,c(0.1,0.3,0.2))
#'

ModelU4 <- function(x,x0,del,param)
{
#
#
#
#
# pcirX(del,x,x0,3)=(1/(26542080*sqrt(del)*sqrt(2*Pi)*sigma^6*x^2*x0^2))*(E^(((-(2+del*kappa))*x+4*sqrt(x)*sqrt(x0)+(-2+del*kappa)*x0)/(del*sigma^2))*(sqrt (x)/sigma)^(-(1/2)+(2*alpha*kappa)/sigma^2)*(sqrt (x0)/sigma)^(3/2-(2*alpha*kappa)/sigma^2)*(26542080*sigma^6*x^(3/2)*x0^(3/2)-276480*del*sigma^4*x*x0*(48*alpha^2*kappa^2-48*alpha*kappa*sigma^2+9*sigma^4+16*kappa^2*x^(3/2)*sqrt(x0)+16*kappa^2*x*x0+16*kappa^2*sqrt(x)*sqrt(x0)*(-6*alpha+x0))+1440*del^2*sigma^2*sqrt(x)*sqrt(x0)*(9*(256*alpha^4*kappa^4-512*alpha^3*kappa^3*sigma^2+224*alpha^2*kappa^2*sigma^4+32*alpha*kappa*sigma^6-15*sigma^8)+256*kappa^4*x^3*x0+512*kappa^4*x^(5/2)*x0^(3/2)+96*kappa^2*(16*alpha^2*kappa^2-16*alpha*kappa*sigma^2+3*sigma^4)*sqrt(x)*sqrt(x0)*(-6*alpha+x0)+768*kappa^4*x^2*x0*(-4*alpha+x0)+32*kappa^2*x*x0*(336*alpha^2*kappa^2-48*alpha*kappa*sigma^2-3*sigma^4-96*alpha*kappa^2*x0+8*kappa^2*x0^2)+32*kappa^2*x^(3/2)*sqrt(x0)*(48*alpha^2*kappa^2-48*alpha*kappa*sigma^2+9*sigma^4-96*alpha*kappa^2*x0+16*kappa^2*x0^2))-del^3*(135*(4096*alpha^6*kappa^6-12288*alpha^5*kappa^5*sigma^2+6400*alpha^4*kappa^4*sigma^4+7680*alpha^3*kappa^3*sigma^6-5456*alpha^2*kappa^2*sigma^8-432*alpha*kappa*sigma^10+315*sigma^12)+20480*kappa^6*x^(9/2)*x0^(3/2)+61440*kappa^6*x^4*x0^2+2160*kappa^2*(256*alpha^4*kappa^4-512*alpha^3*kappa^3*sigma^2+224*alpha^2*kappa^2*sigma^4+32*alpha*kappa*sigma^6-15*sigma^8)*sqrt(x)*sqrt(x0)*(-6*alpha+x0)+122880*kappa^6*x^(7/2)*x0^(3/2)*(-3*alpha+x0)+720*kappa^2*(16*alpha^2*kappa^2-16*alpha*kappa*sigma^2+3*sigma^4)*x*x0*(624*alpha^2*kappa^2-48*alpha*kappa*sigma^2+9*sigma^4-192*alpha*kappa^2*x0+16*kappa^2*x0^2)+1536*kappa^4*x^(5/2)*x0^(3/2)*(1680*alpha^2*kappa^2-240*alpha*kappa*sigma^2-63*sigma^4-720*alpha*kappa^2*x0+80*kappa^2*x0^2)+1280*kappa^4*x^3*x0*(144*alpha^2*kappa^2-144*alpha*kappa*sigma^2+27*sigma^4-576*alpha*kappa^2*x0+112*kappa^2*x0^2)+768*kappa^4*x^2*x0*(-180*alpha*(16*alpha^2*kappa^2-16*alpha*kappa*sigma^2+3*sigma^4)+9*(400*alpha^2*kappa^2-80*alpha*kappa*sigma^2-17*sigma^4)*x0-960*alpha*kappa^2*x0^2+80*kappa^2*x0^3)+16*kappa^2*x^(3/2)*sqrt(x0)*(135*(256*alpha^4*kappa^4-512*alpha^3*kappa^3*sigma^2+224*alpha^2*kappa^2*sigma^4+32*alpha*kappa*sigma^6-15*sigma^8)-8640*alpha*kappa^2*(48*alpha^2*kappa^2-16*alpha*kappa*sigma^2-sigma^4)*x0+288*kappa^2*(560*alpha^2*kappa^2-80*alpha*kappa*sigma^2-21*sigma^4)*x0^2-23040*alpha*kappa^4*x0^3+1280*kappa^4*x0^4))))
#
#
# ---------------------------
  #
# numerical example to try to check the translation:
  # kappa = 0.5  alpha = 0.07  sigma = 0.15  del=1/12
  #
  #  at x0=0.06 and x=0.08, the value of the density is:
  #
  # closed form = 6.6937743258
  # order 2 approx = 6.693777568
  # order 3 approx (the formula above) = 6.69377439438

  kappa <- param[1]
  alpha <- param[2]
  sigma <- param[3]
  output <- list()
  output$llk <- (1/(26542080*sqrt(del)*sqrt(2*pi)*sigma^6*x^2*x0^2))*
  (exp(((-(2+del*kappa))*x+4*sqrt(x)*sqrt(x0)+(-2+del*kappa)*x0)/(del*sigma^2))*(sqrt (x)/sigma)^(-(1/2)+(2*alpha*kappa)/sigma^2)*(sqrt (x0)/sigma)^(3/2-(2*alpha*kappa)/sigma^2)*(26542080*sigma^6*x^(3/2)*x0^(3/2)-276480*del*sigma^4*x*x0*(48*alpha^2*kappa^2-48*alpha*kappa*sigma^2+9*sigma^4+16*kappa^2*x^(3/2)*sqrt(x0)+16*kappa^2*x*x0+16*kappa^2*sqrt(x)*sqrt(x0)*(-6*alpha+x0))+1440*del^2*sigma^2*sqrt(x)*sqrt(x0)*(9*(256*alpha^4*kappa^4-512*alpha^3*kappa^3*sigma^2+224*alpha^2*kappa^2*sigma^4+32*alpha*kappa*sigma^6-15*sigma^8)+256*kappa^4*x^3*x0+512*kappa^4*x^(5/2)*x0^(3/2)+96*kappa^2*(16*alpha^2*kappa^2-16*alpha*kappa*sigma^2+3*sigma^4)*sqrt(x)*sqrt(x0)*(-6*alpha+x0)+768*kappa^4*x^2*x0*(-4*alpha+x0)+32*kappa^2*x*x0*(336*alpha^2*kappa^2-48*alpha*kappa*sigma^2-3*sigma^4-96*alpha*kappa^2*x0+8*kappa^2*x0^2)+32*kappa^2*x^(3/2)*sqrt(x0)*(48*alpha^2*kappa^2-48*alpha*kappa*sigma^2+9*sigma^4-96*alpha*kappa^2*x0+16*kappa^2*x0^2))-del^3*(135*(4096*alpha^6*kappa^6-12288*alpha^5*kappa^5*sigma^2+6400*alpha^4*kappa^4*sigma^4+7680*alpha^3*kappa^3*sigma^6-5456*alpha^2*kappa^2*sigma^8-432*alpha*kappa*sigma^10+315*sigma^12)+20480*kappa^6*x^(9/2)*x0^(3/2)+61440*kappa^6*x^4*x0^2+2160*kappa^2*(256*alpha^4*kappa^4-512*alpha^3*kappa^3*sigma^2+224*alpha^2*kappa^2*sigma^4+32*alpha*kappa*sigma^6-15*sigma^8)*sqrt(x)*sqrt(x0)*(-6*alpha+x0)+122880*kappa^6*x^(7/2)*x0^(3/2)*(-3*alpha+x0)+720*kappa^2*(16*alpha^2*kappa^2-16*alpha*kappa*sigma^2+3*sigma^4)*x*x0*(624*alpha^2*kappa^2-48*alpha*kappa*sigma^2+9*sigma^4-192*alpha*kappa^2*x0+16*kappa^2*x0^2)+1536*kappa^4*x^(5/2)*x0^(3/2)*(1680*alpha^2*kappa^2-240*alpha*kappa*sigma^2-63*sigma^4-720*alpha*kappa^2*x0+80*kappa^2*x0^2)+1280*kappa^4*x^3*x0*(144*alpha^2*kappa^2-144*alpha*kappa*sigma^2+27*sigma^4-576*alpha*kappa^2*x0+112*kappa^2*x0^2)+768*kappa^4*x^2*x0*(-180*alpha*(16*alpha^2*kappa^2-16*alpha*kappa*sigma^2+3*sigma^4)+9*(400*alpha^2*kappa^2-80*alpha*kappa*sigma^2-17*sigma^4)*x0-960*alpha*kappa^2*x0^2+80*kappa^2*x0^3)+16*kappa^2*x^(3/2)*sqrt(x0)*(135*(256*alpha^4*kappa^4-512*alpha^3*kappa^3*sigma^2+224*alpha^2*kappa^2*sigma^4+32*alpha*kappa*sigma^6-15*sigma^8)-8640*alpha*kappa^2*(48*alpha^2*kappa^2-16*alpha*kappa*sigma^2-sigma^4)*x0+288*kappa^2*(560*alpha^2*kappa^2-80*alpha*kappa*sigma^2-21*sigma^4)*x0^2-23040*alpha*kappa^4*x0^3+1280*kappa^4*x0^4))))
  output$llk <- log(output$llk)

  return(output)

}

# ModelU4(3,4,1/23,c(0.1,0.2,0.3))
