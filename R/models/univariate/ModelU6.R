#' Compute the maximum likelihood estimate of Model U6
#'
#' mu[x_] = a + b*x + c*x^2 + d*x^3  
#' s[x_] = f  
#' @param x: Observation of the state variable at time t 
#' @param x0: Observation of the state variable at time t-1 
#' @param del: The time step between the current and previous observation
#' @param param: The parameter 5-vector (a,b,c,d,f)
#' @export output: a list with a llk variable storing the result of the log likelihood calculation
#' @examples
#' ModelU6(0.4,0.3,0.1,c(0.1,0.3,0.2,0.3,0.1))
#'

ModelU6 <- function(x,x0,del,param)
{
  # (* If you use this code,you must acknowledge the source:
       # Ait-Sahalia,Yacine,1996 , 
     #   Testing Continuous-Time Models of the Spot Interest Rate, 
     #   Review of Financial Studies,9,385-426.
     # Ait-Sahalia,Yacine,1999, 
     #   Transition Densities for Interest Rate and Other Nonlinear Densities, 
     #   Journal of Finance 54,1361-1395.
     # Ait-Sahalia,Yacine,2002, 
     #   Maximum-Likelihood Estimation of Discretely Sampled Diffusions: 
       #         A Closed-Form Approximation Approach, Econometrica 70,223-262. *)
  # 
  # mu[x_] = a + b*x + c*x^2 + d*x^3  
  # s[x_] = f 
  # delorder = 2 
  
  a = param[1]
  b = param[2]
  c = param[3]
  d = param[4]
  f = param[5]
  
  sx = f 
  
  cm1 = -(x - x0)^2/(2 * f^2) 
  c0 = (4* c* x^3 + 3 *d *x^4 + 12 *a *(x - x0) - 4 *c *x0^3 - 3 *d *x0^4 + 6 *b *(x^2 - x0^2))/(12* f^2) 
  c1 = -1/(420 * f^2)* (210* a^2 + 70* b^2 *(x^2 + x *x0 + x0^2) +  
                        35* a *(6 *b *(x + x0) + 4 *c *(x^2 + x *x0 + x0^2) +  
                                3 *d *(x^3 + x^2 *x0 + x *x0^2 + x0^3)) +  
                        21 *b *(10 *f^2 + 5 *c *(x^3 + x^2 *x0 + x *x0^2 + x0^3) +  
                                4 *d *(x^4 + x^3 *x0 + x^2 *x0^2 + x *x0^3 + x0^4)) +  
                        2 *(21 *c^2 *(x^4 + x^3 *x0 + x^2 *x0^2 + x *x0^3 + x0^4) +  
                            35* c *(x + x0) *(3 *f^2 + d *(x^4 + x^2 *x0^2 + x0^4)) +  
                            15 *d *(7 *f^2 *(x^2 + x *x0 + x0^2) +  
                                    d *(x^6 + x^5 *x0 + x^4 *x0^2 + x^3 *x0^3 + x^2 *x0^4 + x *x0^5 +  
                                        x0^6)))) 
  c2 = 1/210 * (-35 *b^2 - 105* d* f^2 - 63* c^2* x^2 - 140* c* d* x^3 - 75 *d^2 *x^4 -  
                84 *c^2* x *x0 - 210 *c *d *x^2 *x0 - 120* d^2* x^3* x0 - 63* c^2* x0^2 -  
                210 *c *d *x *x0^2 - 135* d^2* x^2* x0^2 - 140* c* d* x0^3 - 120* d^2* x* x0^3 -  
                75* d^2* x0^4 - 35* a* (2* c + 3* d* (x + x0)) -  
                21 *b *(5* c *(x + x0) + 2 *d *(3 *x^2 + 4 *x *x0 + 3 *x0^2))) 
  output <- list()
  output$llk <- -log(2* pi * del)/2 - log(sx) + cm1/del + c0 + c1*del + c2*del^2/2 
  
  
  return(output)

}

# ModelU6(3,4,1/23,c(0.1,0.2,0.3,0.4,0.5))
