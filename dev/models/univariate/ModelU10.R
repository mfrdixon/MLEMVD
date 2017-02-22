ModelU10 <- function(x,x0,del,param)
{
# (* If you use this code, you must acknowledge the source: 
     # 
   # Ait-Sahalia,Yacine,1996 , Testing Continuous-Time Models of the Spot Interest Rate, Review of Financial Studies,9,385-426. 
   # 
   # Ait-Sahalia,Yacine,1999, Transition Densities for Interest Rate and Other Nonlinear Densities, Journal of Finance 54,1361-1395. 
   # 
   # Ait-Sahalia,Yacine,2002, Maximum-Likelihood Estimation of Discretely Sampled Diffusions: A Closed-Form Approximation Approach, Econometrica 70,223-262. *) 

  am1 <- param[1] 
  a0 <- param[2] 
  a1 <- param[3] 
  a2 <- param[4] 
  b0 <- param[5] 
  b1 <- param[6] 
  b2 <- param[7] 
  b3 <- param[8] 
  
  sx <- b0 + b1*x + b2*x^b3 
  cm1 <- -((x - x0)^2/(2*(b0 + b1*x0 + b2*x0^b3)^2)) + ((x - x0)^3*(b1 + b2*b3*x0^(-1 + b3)))/(2*(b0 + b1*x0 + b2*x0^b3)^3) + ((x - x0)^4*(-11*(b1 + b2*b3*x0^(-1 + b3))^2 + 4*b2*(-1 + b3)*b3* x0^(-2 + b3)*(b0 + b1*x0 + b2*x0^b3)))/(24*(b0 + b1*x0 + b2*x0^b3)^4)  
  c0 <- ((x - x0)*((-(b1 + b2*b3*x0^(-1 + b3)))*(b0 + b1*x0 + b2*x0^b3) + 2*(a0 + am1/x0 + x0*(a1 + a2*x0))))/(2*(b0 + b1*x0 + b2*x0^b3)^2) + ((x - x0)^2*((-b2)*(-1 + b3)*b3*x0^(-2 + b3)*(b0 + b1*x0 + b2*x0^b3)^2 - 4*(b1 + b2*b3*x0^(-1 + b3))*(a0 + am1/x0 + x0*(a1 + a2*x0)) + (b0 + b1*x0 + b2*x0^b3)*(2*(a1 - am1/x0^2 + 2*a2*x0) + (b1 + b2*b3*x0^(-1 + b3))^2)))/(4*(b0 + b1*x0 + b2*x0^b3)^3) 
  c1 <- (-(1/(8*(b0 + b1*x0 + b2*x0^b3)^2)))*(-8*(b1 + b2*b3*x0^(-1 + b3))*(b0 + b1*x0 + b2*x0^b3)* (a0 + am1/x0 + x0*(a1 + a2*x0)) + 4*(a0 + am1/x0 + x0*(a1 + a2*x0))^2 + (b0 + b1*x0 + b2*x0^b3)^2*(4*(a1 - am1/x0^2 + 2*a2*x0) + (b1 + b2*b3*x0^(-1 + b3))^2 - 2*b2*(-1 + b3)*b3*x0^(-2 + b3)*(b0 + b1*x0 + b2*x0^b3)))  
  output <- list()
  output$llk <- -(1/2)*log(2*pi*del) - log(sx) + cm1/del + c0 + c1*del  

  return(output)
  
}

# ModelU10(3,4,1/52,c(0.1,0.2,0.3,0.1,0.2,0.3,0.1,0.2))
