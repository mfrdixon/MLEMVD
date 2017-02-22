#' Compute the maximum likelihood estimate of Model U5
#'
#' mu(x) <- theta0 + theta1*x + theta2*x^2 + theta3*x^3 
#' s(x) <- gamma*x^rho 
#' @param x: Observation of the state variable at time t 
#' @param x0: Observation of the state variable at time t-1 
#' @param del: The time step between the current and previous observation
#' @param param: The parameter 6-vector (gamma,rho>=1,theta0,theta1,theta2,theta3)
#' @export output: a list with a llk variable storing the result of the log likelihood calculation
#' @examples
#' ModelU5(0.4,0.3,0.1,c(0.1,1.3,0.2,0.3,0.2,0.1))
#'

ModelU5 <- function(x,x0,del,param)
{

  gamma <- param[1]
  rho <- param[2] 
  theta0 <- param[3]
  theta1 <- param[4]
  theta2 <- param[5]
  theta3 <- param[6]
  
  alpham1 <- rho/(2*(-1 + rho)) 
  alpha0 <- ((-gamma^(1/(-1 + rho)))*theta0)/(-1 + rho)^(rho/(1 - rho)) 
  alpha1 <- theta1*(1 - rho) 
  alpha2 <- (-gamma^(-(1/(-1 + rho))))*(-1 + rho)^((-2 + rho)/(-1 + rho))*theta2 
  alpha3 <- (-gamma^(-(2/(-1 + rho))))*(-1 + rho)^((-3 + rho)/(-1 + rho))*theta3 
  
  sx <- gamma*x^rho 
  y <- x^(1 - rho)/(gamma*(rho - 1)) 
  y0 <- x0^(1 - rho)/(gamma*(rho - 1)) 
  cYm1 <- (-(1/2))*(y - y0)^2 
  cY0 <- -((6*alpha1*y^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3)) + (19*alpha1*rho*y^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (16*alpha1*rho^2*y^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (4*alpha1*rho^3*y^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (12*alpha0*y^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (26*alpha0*rho*y^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (18*alpha0*rho^2*y^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (4*alpha0*rho^3*y^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (3*alpha3*y^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (11*alpha3*rho*y^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (12*alpha3*rho^2*y^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (4*alpha3*rho^3*y^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (4*alpha2*y^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (14*alpha2*rho*y^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (14*alpha2*rho^2*y^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) 
         + (4*alpha2*rho^3*y^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (6*alpha1*y0^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (19*alpha1*rho*y0^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (16*alpha1*rho^2*y0^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (4*alpha1*rho^3*y0^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (12*alpha0*y0^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (26*alpha0*rho*y0^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (18*alpha0*rho^2*y0^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (4*alpha0*rho^3*y0^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (3*alpha3*y0^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (11*alpha3*rho*y0^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (12*alpha3*rho^2*y0^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (4*alpha3*rho^3*y0^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (4*alpha2*y0^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (14*alpha2*rho*y0^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (14*alpha2*rho^2*y0^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (4*alpha2*rho^3*y0^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (12*alpham1*log(y/y0))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (38*alpham1*rho*log(y/y0))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (32*alpham1*rho^2*log(y/y0))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (8*alpham1*rho^3*log(y/y0))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) 
  cY1 <- -(alpha1/2) - alpha1*alpham1 - alpham1/(2*y*(y - y0)) + alpham1^2/(2*y*(y - y0)) - (alpha1^2*y^3)/(6*(y - y0)) - (alpha0*alpha2*y^3)/(3*(y - y0)) + alpham1/(2*(y - y0)*y0) - alpham1^2/(2*(y - y0)*y0) + (alpha1^2*y0^3)/(6*(y - y0)) + (alpha0*alpha2*y0^3)/(3*(y - y0)) + (alpha2*alpham1*(y^(1 + 1/(1 - rho)) - y0^(1 + 1/(1 - rho))))/((-2 + rho)*(y - y0)) - (alpha2*rho*(y^(1 + 1/(1 - rho)) - y0^(1 + 1/(1 - rho))))/(2*(-2 + rho)*(y - y0)) - (alpha2*alpham1*rho*(y^(1 + 1/(1 - rho)) - y0^(1 + 1/(1 - rho))))/((-2 + rho)*(y - y0)) - (alpha2*(-y^(1 + 1/(1 - rho)) + y0^(1 + 1/(1 - rho))))/((-2 + rho)*(y - y0)) - (alpha1*alpha2*(-1 + rho)*(y^(3 + 1/(1 - rho)) - y0^(3 + 1/(1 - rho))))/((-4 + 3*rho)*(y - y0)) - (alpha0*alpha3*(-1 + rho)*(y^(3 + 1/(1 - rho)) - y0^(3 + 1/(1 - rho))))/ ((-4 + 3*rho)*(y - y0)) - (alpha3^2*(-1 + rho)*(y^(3 - 4/(-1 + rho)) - y0^(3 - 4/(-1 + rho))))/(2*(-7 + 3*rho)*(y - y0)) - (alpha2*alpha3*(-1 + rho)*(y^(3 - 3/(-1 + rho)) - y0^(3 - 3/(-1 + rho))))/(3*(-2 + rho)*(y - y0)) - (alpha2^2*(-1 + rho)*(y^(3 - 2/(-1 + rho)) - y0^(3 - 2/(-1 + rho))))/ (2*(-5 + 3*rho)*(y - y0)) - (alpha1*alpha3*(-1 + rho)*(y^(3 - 2/(-1 + rho)) - y0^(3 - 2/(-1 + rho))))/((-5 + 3*rho)*(y - y0)) - (alpha0*alpha1*(-1 + rho)*(y^(3 + 1/(-1 + rho)) - y0^(3 + 1/(-1 + rho))))/((-2 + 3*rho)*(y - y0)) - (alpha0^2*(-1 + rho)*(y^(3 + 2/(-1 + rho)) - y0^(3 + 2/(-1 + rho))))/ (2*(-1 + 3*rho)*(y - y0)) + (3*alpha3*(y^((-3 + rho)/(-1 + rho)) - y0^((-3 + rho)/(-1 + rho))))/(2*(-3 + rho)*(y - y0)) + 
          (alpha3*alpham1*(y^((-3 + rho)/(-1 + rho)) - y0^((-3 + rho)/(-1 + rho))))/((-3 + rho)*(y - y0)) - (alpha3*rho*(y^((-3 + rho)/(-1 + rho)) - y0^((-3 + rho)/(-1 + rho))))/ (2*(-3 + rho)*(y - y0)) 
        - (alpha3*alpham1*rho*(y^((-3 + rho)/(-1 + rho)) - y0^((-3 + rho)/(-1 + rho))))/((-3 + rho)*(y - y0)) - (alpha0*(y^(rho/(-1 + rho)) - y0^(rho/(-1 + rho))))/(2*(y - y0)) - (alpha0*alpham1*(y^(rho/(-1 + rho)) - y0^(rho/(-1 + rho))))/(y - y0) - (alpha0*alpham1*(-y^(rho/(-1 + rho)) + y0^(rho/(-1 + rho))))/(rho*(y - y0)) output <- -(1/2)*log(2*pi*del) - log(sx) + cYm1/del + cY0 + cY1*del 
  output <-list()
  output$llk <- -(1/2)*log(2*pi*del) - log(sx) + cYm1/del + cY0 + cY1*del 
  
  
  return(output)
}


# Produced by file Nonlinear Drift with CEV.nb
# 
# 
# If you use this code,you must acknowledge the source:
  # 
# Ait-Sahalia,Yacine,1999, Transition Densities for Interest Rate and Other Nonlinear Densities, Journal of Finance 54,1361-1395. 
# 
# Ait-Sahalia,Yacine,2002, Maximum-Likelihood Estimation of Discretely Sampled Diffusions: A Closed-Form Approximation Approach, Econometrica 70,223-262.
# 
# 
# 
# 
# DO NOT LET THE PARAMETER rho GO BELOW 1! 
  # 
# 
# Transformation:
  # 
# g(x) <- x^(1 - rho)/(gamma*(rho - 1)) 
# 
# 
# New parameters:
  # 
# alpham1 <- rho/(2*(-1 + rho)) 
# 
# alpha0 <- ((-gamma^(1/(-1 + rho)))*theta0)/(-1 + rho)^(rho/(1 - rho)) 
# 
# alpha1 <- theta1*(1 - rho) 
# 
# alpha2 <- (-gamma^(-(1/(-1 + rho))))*(-1 + rho)^((-2 + rho)/(-1 + rho))*theta2 
# 
# alpha3 <- (-gamma^(-(2/(-1 + rho))))*(-1 + rho)^((-3 + rho)/(-1 + rho))*theta3 
# 
# log-density function at order 1 in del:
  # 
# lnpX(del, x, x0, 1, exact) <- -(1/2)*log(2*Pi*del) - log(s(x)) + cY(g(x), g(x0), -1, exact)/del + cY(g(x), g(x0), 0, exact) + cY(g(x), g(x0), 1, exact)*del 
# 
# 
# 
# Coefficients:
  # 
# cY(y, y0, -1, exact) <- (-(1/2))*(y - y0)^2 
# 
# cY(y, y0, 0, exact) <- -((6*alpha1*y^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3)) + (19*alpha1*rho*y^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (16*alpha1*rho^2*y^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (4*alpha1*rho^3*y^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (12*alpha0*y^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (26*alpha0*rho*y^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (18*alpha0*rho^2*y^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (4*alpha0*rho^3*y^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (3*alpha3*y^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (11*alpha3*rho*y^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (12*alpha3*rho^2*y^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (4*alpha3*rho^3*y^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (4*alpha2*y^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (14*alpha2*rho*y^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (14*alpha2*rho^2*y^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (4*alpha2*rho^3*y^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (6*alpha1*y0^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (19*alpha1*rho*y0^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (16*alpha1*rho^2*y0^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (4*alpha1*rho^3*y0^2)/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (12*alpha0*y0^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (26*alpha0*rho*y0^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (18*alpha0*rho^2*y0^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (4*alpha0*rho^3*y0^(2 + 1/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (3*alpha3*y0^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (11*alpha3*rho*y0^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (12*alpha3*rho^2*y0^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (4*alpha3*rho^3*y0^((-4 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (4*alpha2*y0^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (14*alpha2*rho*y0^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (14*alpha2*rho^2*y0^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (4*alpha2*rho^3*y0^((-3 + 2*rho)/(-1 + rho)))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (12*alpham1*log(y/y0))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (38*alpham1*rho*log(y/y0))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) - (32*alpham1*rho^2*log(y/y0))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) + (8*alpham1*rho^3*log(y/y0))/(-12 + 38*rho - 32*rho^2 + 8*rho^3) 
# 
# cY(y, y0, 1, exact) <- -(alpha1/2) - alpha1*alpham1 - alpham1/(2*y*(y - y0)) + alpham1^2/(2*y*(y - y0)) - (alpha1^2*y^3)/(6*(y - y0)) - (alpha0*alpha2*y^3)/(3*(y - y0)) + alpham1/(2*(y - y0)*y0) - alpham1^2/(2*(y - y0)*y0) + (alpha1^2*y0^3)/(6*(y - y0)) + (alpha0*alpha2*y0^3)/(3*(y - y0)) + (alpha2*alpham1*(y^(1 + 1/(1 - rho)) - y0^(1 + 1/(1 - rho))))/((-2 + rho)*(y - y0)) - (alpha2*rho*(y^(1 + 1/(1 - rho)) - y0^(1 + 1/(1 - rho))))/(2*(-2 + rho)*(y - y0)) - (alpha2*alpham1*rho*(y^(1 + 1/(1 - rho)) - y0^(1 + 1/(1 - rho))))/((-2 + rho)*(y - y0)) - (alpha2*(-y^(1 + 1/(1 - rho)) + y0^(1 + 1/(1 - rho))))/((-2 + rho)*(y - y0)) - (alpha1*alpha2*(-1 + rho)*(y^(3 + 1/(1 - rho)) - y0^(3 + 1/(1 - rho))))/((-4 + 3*rho)*(y - y0)) - (alpha0*alpha3*(-1 + rho)*(y^(3 + 1/(1 - rho)) - y0^(3 + 1/(1 - rho))))/ ((-4 + 3*rho)*(y - y0)) - (alpha3^2*(-1 + rho)*(y^(3 - 4/(-1 + rho)) - y0^(3 - 4/(-1 + rho))))/(2*(-7 + 3*rho)*(y - y0)) - (alpha2*alpha3*(-1 + rho)*(y^(3 - 3/(-1 + rho)) - y0^(3 - 3/(-1 + rho))))/(3*(-2 + rho)*(y - y0)) - (alpha2^2*(-1 + rho)*(y^(3 - 2/(-1 + rho)) - y0^(3 - 2/(-1 + rho))))/ (2*(-5 + 3*rho)*(y - y0)) - (alpha1*alpha3*(-1 + rho)*(y^(3 - 2/(-1 + rho)) - y0^(3 - 2/(-1 + rho))))/((-5 + 3*rho)*(y - y0)) - (alpha0*alpha1*(-1 + rho)*(y^(3 + 1/(-1 + rho)) - y0^(3 + 1/(-1 + rho))))/((-2 + 3*rho)*(y - y0)) - (alpha0^2*(-1 + rho)*(y^(3 + 2/(-1 + rho)) - y0^(3 + 2/(-1 + rho))))/ (2*(-1 + 3*rho)*(y - y0)) + (3*alpha3*(y^((-3 + rho)/(-1 + rho)) - y0^((-3 + rho)/(-1 + rho))))/(2*(-3 + rho)*(y - y0)) + (alpha3*alpham1*(y^((-3 + rho)/(-1 + rho)) - y0^((-3 + rho)/(-1 + rho))))/((-3 + rho)*(y - y0)) - (alpha3*rho*(y^((-3 + rho)/(-1 + rho)) - y0^((-3 + rho)/(-1 + rho))))/ (2*(-3 + rho)*(y - y0)) - (alpha3*alpham1*rho*(y^((-3 + rho)/(-1 + rho)) - y0^((-3 + rho)/(-1 + rho))))/((-3 + rho)*(y - y0)) - (alpha0*(y^(rho/(-1 + rho)) - y0^(rho/(-1 + rho))))/(2*(y - y0)) - (alpha0*alpham1*(y^(rho/(-1 + rho)) - y0^(rho/(-1 + rho))))/(y - y0) - (alpha0*alpham1*(-y^(rho/(-1 + rho)) + y0^(rho/(-1 + rho))))/(rho*(y - y0)) 
