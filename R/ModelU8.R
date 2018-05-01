#' Compute the maximum likelihood estimate of Model U8
#'
#' @param x Observation of the state variable at time t
#' @param x0 Observation of the state variable at time t-1
#' @param del The time step between the current and previous observation
#' @param param The parameter 6-vector (am1,a0,a1,a2,sigma,rho>=1)
#' @export output a list with a llk variable storing the result of the log likelihood calculation
#' @examples
#' ModelU8(0.4,0.3,0.1,c(0.1,0.3,0.2,0.1,0.2,1.3))
#'

ModelU8 <- function(x,x0,del,param)

# (* DO NOT LET THE PARAMETER rho GO BELOW 1! THE EXPANSION IS DIFFERENT THERE.
   #
   # If you use this code, please acknowledge the source:Ait-Sahalia,Yacine,1999, Transition Densities for Interest Rate and Other Nonlinear Densities, Journal of Finance 54,1361-1395.
   # Ait-Sahalia,Yacine,2002, Maximum-Likelihood Estimation of Discretely Sampled Diffusions: A Closed-Form Approximation Approach, Econometrica 70,223-262. *)

{
  am1 <- param[1]
  a0 <- param[2]
  a1 <- param[3]
  a2 <- param[4]
  sigma <- param[5]
  rho <- param[6]

  if (rho<1)
  # If rho<1, set the log density to -Infinity for estimation purpose.
  {
    output <- -Inf
  }


  y <- x^(1 - rho)/(sigma*(rho - 1))
  y0 <- x0^(1 - rho)/(sigma*(rho - 1))
  sx <- sigma*x^rho

  cYm1 <- -(1/2)*(y - y0)^2
  cY0 <- -(1/2) * a1 * (-1 + rho) * (y - y0) * (y + y0) - ( a2 * (-1 + rho)^((-3 + 2 * rho)/(-1 + rho)) *
                                                           sigma^(1/(  1 - rho)) * (y^(2 + 1/(1 - rho)) - y0^(2 + 1/(1 - rho))))/(-3 + 2 * rho) - ( a0 * (-1 + rho)^(1 + rho/(-1 + rho)) *
                                                                                                                                                    sigma^(  1/(-1 + rho)) * (y^(2 + 1/(-1 + rho)) - y0^(2 + 1/(-1 + rho))))/(-1 +   2 * rho)
  - ( am1 * (-1 + rho)^((2 * rho)/(-1 + rho))* sigma^(  2/(-1 + rho)) * (y^((2 * rho)/(-1 + rho)) - y0^((2 * rho)/(-1 + rho))))/( 2 * rho)
  + (rho * log(y/y0))/(-2 + 2 * rho)
  cY1 <- -(1/(y - y0)) * (( am1^2 * (-1 + rho)^((1 + 3 * rho)/(-1 + rho)) * sigma^( 4/(-1 + rho)) *
                           (y^(3 + 4/(-1 + rho)) - y0^(3 + 4/(-1 + rho))))/( 2 + 6 * rho) + 1/6 * (a1^2 * (-1 + rho)^2 * (y^3 - y0^3) +
                                                                                                   ( 3 * a0^2 * (-1 + rho)^(1 + (2 * rho)/(-1 + rho)) * sigma^( 2/(-1 + rho)) * (y^(3 + 2/(-1 + rho)) - y0^(3 + 2/(-1 + rho))))/(-1 + 3 * rho)
                                                                                                   + ( 3 * a1 * (1 - rho) * (((4 - 11 * rho + 6 * rho^2) * (y - y0))/(-1 + rho) - 2 * a2 * (-1 + rho)^((-3 + 2 * rho)/(-1 + rho)) *
                                                                                                                             sigma^(1/( 1 - rho)) * (y^(3 + 1/(1 - rho)) - y0^(3 + 1/(1 - rho)))))/(-4 + 3 * rho)
                                                                                                   + 3 * (-(((-2 + rho) * rho * (y - y0))/(4 * (-1 + rho)^2 * y * y0))
                                                                                                          + ( a2^2 * (-1 + rho)^((-5 + 3 * rho)/(-1 + rho)) * sigma^(-( 2/(-1 + rho))) * (y^(3 - 2/(-1 + rho)) - y0^( 3 - 2/(-1 + rho))))/(-5 + 3 * rho)
                                                                                                          - ( 2 * a2 * (-1 + rho)^((-3 + 2 * rho)/(-1 + rho)) * sigma^(1/( 1 - rho)) * (y^((-2 + rho)/(-1 + rho))
                                                                                                                                                                                        - y0^((-2 + rho)/(-1 + rho))))/(-2 + rho)) - (1/(rho * (-2 + 3 * rho))) * a0 * (-1 + rho)^(rho/(-1 + rho)) * sigma^( 1/(-1 + rho)) *
                                                                                                   (-6 * a1 * (-1 + rho)^2 * rho * (y^(3 + 1/(-1 + rho)) - y0^( 3 + 1/(-1 + rho))) + (-2 + 3 * rho) * (-2 * a2 *(-1 + rho)^((-2 + rho)/(-1 + rho)) *
                                                                                                                                                                                                       rho * sigma^(1/( 1 - rho)) * (y^3 - y0^3) + 6 * rho * (y^(rho/(-1 + rho)) - y0^(rho/(-1 + rho))))))
                         - (1/( 6 * rho * (2 - 7 * rho + 9 * rho^3))) * am1 * (-1 + rho)^((1 + rho)/(-1 + rho)) * sigma^( 2/(-1 + rho)) *
                         (-2 * a0 * (-1 + rho)^( rho/(-1 + rho)) * (-2 + 9 * rho - 7 * rho^2 - 9 * rho^3 + 9 * rho^4) * sigma^( 1/(-1 + rho)) *
                           (y^((3 * rho)/(-1 + rho)) - y0^((3 * rho)/(-1 + rho))) + 3 * rho * (-2 * a1 * (-1 + rho)^2 * (-2 + rho + 3 * rho^2) * (y^(3 + 2/(-1 + rho))
                                                                                                                                                  - y0^(3 + 2/(-1 + rho))) + (-1 + 3 * rho) * (-2 * a2 * (-1 + rho)^((-2 + rho)/(-1 + rho)) * (-1 + rho^2)
                                                                                                                                                                                               * sigma^(1/( 1 - rho)) * (y^(3 + 1/(-1 + rho)) - y0^(3 + 1/(-1 + rho))) + (1 + 2 * rho) * (-2 + 3 * rho) * (y^((1 + rho)/(-1 + rho))
                                                                                                                                                                                                                                                                                                           - y0^(( 1 + rho)/(-1 + rho)))))))
  output <- list()
  output$llk <- (-(1/2))*log(2*pi*del) - log(sx) + cYm1/del + cY0 + cY1*del

  return(output)
}

# ModelU8(4,5,1/32,c(0.1,0.2,0.3,0.4,0.5,4))
