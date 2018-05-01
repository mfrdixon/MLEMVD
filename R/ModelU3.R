#' Compute the maximum likelihood estimate of Model U3
#'
#' linear drift, CEV diffusion (Univariate Model U3), dX = b(a - X) dt + c X^d dW
#' @param x Observation of the state variable at time t
#' @param x0 Observation of the state variable at time t-1
#' @param del The time step between the current and previous observation
#' @param param The parameter 4-vector (a,b,c,d)
#' @param args Not currently used
#' @export output a list with a llk variable storing the result of the log likelihood calculation
#' @examples
#' ModelU3(0.4,0.3,0.1,c(0.1,0.2,0.3,0.1))
#'

ModelU3 <- function(x,x0,del,param,args)
{
  # formula = log(pX(del,x,x0)) with K=2
  #
  # Model = linear drift, CEV diffusion, dX = b(a - X) dt + c X^d dW
  #
  # separate formulae for 1/2 < d < 1 and d >1 on the one hand, and d = 1 on the other
  # also separate formulae at x = x0 based on l'Hopital's Rule

  a <- param[1]
  b <- param[2]
  c <- param[3]
  d <- param[4]
  output <- list()

  # Check each parameter to see if it is out of bounds
  # Set the flag to 1 if so
  # If flag = 0, evaluate likelihood,
  # otherwise return infinity so this will never be a minimum


  if ((a>0) & (b>0) & (c>0) & (1/2<d) & (d!=1)){

    # case where d is >1/2 but not equal to 1
    if (x != x0){

      output$llk <- -((x^(1 - d)/(c*(-1 + d)) - x0^(1 - d)/(c*(-1 + d)))^2/(2*del)) + (1/(24*c^2*(2 - 9*d + 9*d^2)*((-x^d)*x0 + x*x0^d)))* (del*x^(-1 + 2*d)*x0^(-1 + 2*d)*(3*c^4*(-2 + d)*d*(2 + 9*(-1 + d)*d)*x^(1 - d) + 4*b^2*(2 + 9*(-1 + d)*d)*x^(1 - d)* x0^(4 - 4*d) + 12*b*c^2*(-1 + 2*d)*(2 + 9*(-1 + d)*d)*x^(1 - d)*x0^(2 - 2*d) - 3*c^4*(-2 + d)*d*(2 + 9*(-1 + d)*d)* x0^(1 - d) - 4*b^2*(2 + 9*(-1 + d)*d)*x^(4 - 4*d)*x0^(1 - d) - 12*b*c^2*(-1 + 2*d)*(2 + 9*(-1 + d)*d)*x^(2 - 2*d)* x0^(1 - d) + 12*a^2*b^2*(-1 + d)*(-2 + 3*d)*x^(1 - 4*d)*x0^(1 - 4*d)*(x^(3*d)*x0 - x*x0^(3*d)) - 24*a*b*(-1 + d)*(-1 + 3*d)*x^(1 - 4*d)*x0^(1 - 4*d)*((-b)*x^2*x0^(3*d) - c^2*(-2 + 3*d)*x^(2*d)*x0^(3*d) + x^(3*d)*(b*x0^2 + c^2*(-2 + 3*d)*x0^(2*d))))) + ((c - c*d)^3*del^2*x^(-2 - 3*d)*x0^(-2 - 3*d)*((2 - 9*d + 9*d^2)*(x^d*x0 - x*x0^d)^3* (-4*b^2*x^2*x0^2 + 3*c^4*(-2 + d)*d*x^(2*d)*x0^(2*d)) - 12*a^2*b^2*(-2 + 3*d)*x^2*x0^2* ((1 + d)*x^(3*d)*x0 + (1 - 3*d)*x^(1 + 2*d)*x0^d - (1 + d)*x*x0^(3*d) + (-1 + 3*d)*x^d*x0^(1 + 2*d)) - 24*a*b*(-1 + 3*d)*x*x0*((-c^2)*d*(-2 + 3*d)*x^(3*d)*x0^(2 + 2*d) - b*(-2 + 3*d)*x^(2 + d)*x0^(2 + 2*d) + b*d*x^3*x0^(1 + 3*d) - c^2*(4 - 8*d + 3*d^2)*x^(1 + 2*d)*x0^(1 + 3*d) + (-2 + 3*d)*x^(2 + 2*d)*x0^d* (b*x0^2 + c^2*d*x0^(2*d)) + x^(1 + 3*d)*x0*((-b)*d*x0^2 + c^2*(4 - 8*d + 3*d^2)*x0^(2*d)))))/ (48*c^3*(-1 + d)*(2 + 9*(-1 + d)*d)*(x^(1 - d) - x0^(1 - d))^3) - (1/2)*log(2*del*pi) - log(c*x^d) + (1/(2*c^2*(1 - 3*d + 2*d^2)))* ((b*(-2*a*(-1 + d)*x*x0^(2*d) + (-1 + 2*d)*x^2*x0^(2*d) - x^(2*d)*x0*(2*a - 2*a*d - x0 + 2*d*x0)) - c^2*d*(1 - 3*d + 2*d^2)*x^(2*d)*x0^(2*d)*log(x) + c^2*d*(1 - 3*d + 2*d^2)*x^(2*d)*x0^(2*d)*log(x0))/(x^(2*d)*x0^(2*d)))
    }
    else # case x=x0
    {
      output$llk <- (1/(48*x0^4))*(del^2*(-4*a^2*b^2*d*(1+d)*x0^2+4*a*b^2*d*(-1+2*d)*x0^3-4*b^2*(-1+d)^2*x0^4+3*c^4*(-2+d)*(-1+d)^2*d*x0^(4*d)-4*a*b*c^2*(-2+d)*d*x0^(1+2*d)))+(1/(8*c^2))*((del*(-4*a^2*b^2*x0^2+8*a*b^2*x0^3-4*b^2*x0^4+c^4*(-2+d)*d*x0^(4*d)+8*a*b*c^2*d*x0^(1+2*d)-4*b*c^2*(-1+2*d)*x0^(2+2*d)))/x0^(2*(1+d)))-(1/2)*log(2*del*pi)-log(c*x0^d);
    }
  }

  else if ((a>0) & (b>0) & (c>0) & (d==1)){
  # case where d=1
   if (x != x0){
      output$llk <- (-(1/2))*log(2*del*pi) - log(c*x) - (log(x)/c - log(x0)/c)^2/(2*del) - (2*a*b*(1/x - 1/x0) + (2*b + c^2)*log(x) - (2*b + c^2)*log(x0))/(2*c^2) + (del*(2*a*b*(x - x0)*((-a)*b*x0 + x*((-a)*b + 4*(b + c^2)*x0)) - (2*b + c^2)^2*x^2*x0^2*log(x) + (2*b + c^2)^2*x^2*x0^2*log(x0)))/(8*c^2*x^2*x0^2*(log(x) - log(x0))) + (1/(4*x^2*x0^2*(log(x) - log(x0))^3))* (a*b*del^2*(2*(b + c^2)*x^2*x0*(-2 + log(x) - log(x0)) - a*b*x^2*(-1 + log(x) - log(x0)) + 2*(b + c^2)*x*x0^2*(2 + log(x) - log(x0)) + a*b*x0^2*(-1 - log(x) + log(x0))))
    }
    else # case x=x0
    {
      output$llk <- (a*b*del^2*(-2*a*b+(b+c^2)*x0))/(12*x0^2)-(del*(4*a^2*b^2-8*a*b*(b+c^2)*x0+(2*b+c^2)^2*x0^2))/(8*c^2*x0^2)-(1/2)*log(2*del*pi)-log(c*x0)
    }
  }

  else{
    output$llk <- -Inf # if parameter restrictions are not verified, send f to +infinity
    }


  return(output)
}

