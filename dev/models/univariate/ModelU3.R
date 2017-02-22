# linear drift, CEV diffusion (Univariate Model U3), dX = b(a - X) dt + c X^d dW
# 4 parameters to be estimated: (a,b,c,d)
ModelU3 <- function(z,z0,del,param,args)
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
    if (z!=z0){ 
      
      output$llk = -((z^(1 - d)/(c*(-1 + d)) - z0^(1 - d)/(c*(-1 + d)))^2/(2*del)) + (1/(24*c^2*(2 - 9*d + 9*d^2)*((-z^d)*z0 + z*z0^d)))* (del*z^(-1 + 2*d)*z0^(-1 + 2*d)*(3*c^4*(-2 + d)*d*(2 + 9*(-1 + d)*d)*z^(1 - d) + 4*b^2*(2 + 9*(-1 + d)*d)*z^(1 - d)* z0^(4 - 4*d) + 12*b*c^2*(-1 + 2*d)*(2 + 9*(-1 + d)*d)*z^(1 - d)*z0^(2 - 2*d) - 3*c^4*(-2 + d)*d*(2 + 9*(-1 + d)*d)* z0^(1 - d) - 4*b^2*(2 + 9*(-1 + d)*d)*z^(4 - 4*d)*z0^(1 - d) - 12*b*c^2*(-1 + 2*d)*(2 + 9*(-1 + d)*d)*z^(2 - 2*d)* z0^(1 - d) + 12*a^2*b^2*(-1 + d)*(-2 + 3*d)*z^(1 - 4*d)*z0^(1 - 4*d)*(z^(3*d)*z0 - z*z0^(3*d)) - 24*a*b*(-1 + d)*(-1 + 3*d)*z^(1 - 4*d)*z0^(1 - 4*d)*((-b)*z^2*z0^(3*d) - c^2*(-2 + 3*d)*z^(2*d)*z0^(3*d) + z^(3*d)*(b*z0^2 + c^2*(-2 + 3*d)*z0^(2*d))))) + ((c - c*d)^3*del^2*z^(-2 - 3*d)*z0^(-2 - 3*d)*((2 - 9*d + 9*d^2)*(z^d*z0 - z*z0^d)^3* (-4*b^2*z^2*z0^2 + 3*c^4*(-2 + d)*d*z^(2*d)*z0^(2*d)) - 12*a^2*b^2*(-2 + 3*d)*z^2*z0^2* ((1 + d)*z^(3*d)*z0 + (1 - 3*d)*z^(1 + 2*d)*z0^d - (1 + d)*z*z0^(3*d) + (-1 + 3*d)*z^d*z0^(1 + 2*d)) - 24*a*b*(-1 + 3*d)*z*z0*((-c^2)*d*(-2 + 3*d)*z^(3*d)*z0^(2 + 2*d) - b*(-2 + 3*d)*z^(2 + d)*z0^(2 + 2*d) + b*d*z^3*z0^(1 + 3*d) - c^2*(4 - 8*d + 3*d^2)*z^(1 + 2*d)*z0^(1 + 3*d) + (-2 + 3*d)*z^(2 + 2*d)*z0^d* (b*z0^2 + c^2*d*z0^(2*d)) + z^(1 + 3*d)*z0*((-b)*d*z0^2 + c^2*(4 - 8*d + 3*d^2)*z0^(2*d)))))/ (48*c^3*(-1 + d)*(2 + 9*(-1 + d)*d)*(z^(1 - d) - z0^(1 - d))^3) - (1/2)*log(2*del*pi) - log(c*z^d) + (1/(2*c^2*(1 - 3*d + 2*d^2)))* ((b*(-2*a*(-1 + d)*z*z0^(2*d) + (-1 + 2*d)*z^2*z0^(2*d) - z^(2*d)*z0*(2*a - 2*a*d - z0 + 2*d*z0)) - c^2*d*(1 - 3*d + 2*d^2)*z^(2*d)*z0^(2*d)*log(z) + c^2*d*(1 - 3*d + 2*d^2)*z^(2*d)*z0^(2*d)*log(z0))/(z^(2*d)*z0^(2*d)))
    } 
    else # case x=x0
    {
      output$llk = (1/(48*z0^4))*(del^2*(-4*a^2*b^2*d*(1+d)*z0^2+4*a*b^2*d*(-1+2*d)*z0^3-4*b^2*(-1+d)^2*z0^4+3*c^4*(-2+d)*(-1+d)^2*d*z0^(4*d)-4*a*b*c^2*(-2+d)*d*z0^(1+2*d)))+(1/(8*c^2))*((del*(-4*a^2*b^2*z0^2+8*a*b^2*z0^3-4*b^2*z0^4+c^4*(-2+d)*d*z0^(4*d)+8*a*b*c^2*d*z0^(1+2*d)-4*b*c^2*(-1+2*d)*z0^(2+2*d)))/z0^(2*(1+d)))-(1/2)*log(2*del*pi)-log(c*z0^d);
    }
  }
  
  else if ((a>0) & (b>0) & (c>0) & (d==1)){ 
  # case where d=1
   if (z!=z0){
      output$llk = (-(1/2))*log(2*del*pi) - log(c*z) - (log(z)/c - log(z0)/c)^2/(2*del) - (2*a*b*(1/z - 1/z0) + (2*b + c^2)*log(z) - (2*b + c^2)*log(z0))/(2*c^2) + (del*(2*a*b*(z - z0)*((-a)*b*z0 + z*((-a)*b + 4*(b + c^2)*z0)) - (2*b + c^2)^2*z^2*z0^2*log(z) + (2*b + c^2)^2*z^2*z0^2*log(z0)))/(8*c^2*z^2*z0^2*(log(z) - log(z0))) + (1/(4*z^2*z0^2*(log(z) - log(z0))^3))* (a*b*del^2*(2*(b + c^2)*z^2*z0*(-2 + log(z) - log(z0)) - a*b*z^2*(-1 + log(z) - log(z0)) + 2*(b + c^2)*z*z0^2*(2 + log(z) - log(z0)) + a*b*z0^2*(-1 - log(z) + log(z0))));
    }
    else # case x=x0
    {  
      output$llk = (a*b*del^2*(-2*a*b+(b+c^2)*z0))/(12*z0^2)-(del*(4*a^2*b^2-8*a*b*(b+c^2)*z0+(2*b+c^2)^2*z0^2))/(8*c^2*z0^2)-(1/2)*log(2*del*pi)-log(c*z0);
    }
  }
  
  else{
    output$llk = -inf # if parameter restrictions are not verified, send f to +infinity
    }
   
  
  return(output)
}

