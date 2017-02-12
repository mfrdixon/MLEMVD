ModelU11 <- function(x,x0,del,param)
{

  a <- param[1]
  b <- param[2]
  f <- param[3]
  d <- param[4]
  
  sx <- f + d*x 
  y <- log(1 + (d*x)/f)/d 
  y0 <- log(1 + (d*x0)/f)/d 
  
  E <- exp(1) 
  
  cYm1 <- (-(1/2))*(y - y0)^2  
  
  cY0 <- (E^((-d)*y) - E^((-d)*y0))*((b*f - a*d)/(d^2*f)) + (y - y0)*((2*b - d^2)/(2*d))  
  
  if (y != y0)
  {
    cY1 <- (1/(2*d))*(b^2/(2*d^2) + a^2/(2*f^2) - (a*b)/(d*f))*((E^(-2*d*y) - E^(-2*d*y0))/(y - y0)) + 
      ((a*b)/(d^2*f) - b^2/d^3 + b/d - a/f)*((E^((-d)*y) - E^((-d)*y0))/(y - y0)) - (2*b - d^2)^2/(8*d^2) 
  }
   
  else
  {
    cY1 <- -((a*d - b*f)^2/(E^(2*d*y)*(2*f^2*d^2))) + 
     ((b - d^2)*((-a)*d + b*f))/(E^(d*y)*(f*d^2)) - (2*b - d^2)^2/(8*d^2)  
  }
   
  
  output <- (-(1/2))*log(2*pi*del) - log(sx) + cYm1/del + cY0 + cY1*del  
  
  print(output)
  return(output)
}

# ModelU11(3,4,1/52,c(0.3,0.4,0.5,0.6))
