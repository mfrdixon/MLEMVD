ModelB11<-function(x,x0,del,param){

# dx1 <- (k1 + k2*x2)*dt + sqrt(x2)*(sqrt(1 - rho^2)*dW1 + rho*dW2)
# dx2 <- kappa*(theta - x2)*dt + sigma*x2*dW2

  x1 <- x[1]
  x2 <- x[2]
  x10 <- x0[1]
  x20 <- x0[2]
  
  k1    <- param[1]
  k2    <- param[2] 
  rho   <- param[3]
  kappa <- param[4]
  theta <- param[5]
  sigma <- param[6]
  
  Dv <- (1/2)*log(x2^3*(1 - rho^2)*sigma^2)
  cm1 <- ((x1 - x10)^2*(x2 - x20)^3*(4 - 7*rho^2))/(32*x20^4*(-1 + rho^2)^2) - 
    ((x1 - x10)^2*(x2 - x20))/(4*x20^2*(-1 + rho^2)) + ((x1 - x10)^2*(x2 - x20)^2*
                        (-8 + 11*rho^2))/(48*x20^3*(-1 + rho^2)^2) + 
  ((x1 - x10)^2*(x2 - x20)^4*(384 - 1179*rho^2 + 815*rho^4))/(3840*x20^5*(-1 + rho^2)^3) + ((x2 - x20)^5*(20 - 21*rho^2))/(48*x20^5*(-1 + rho^2)^2*sigma^2) - 
    (x2 - x20)^3/(2*x20^3*(-1 + rho^2)*sigma^2) + ((x2 - x20)^4*(-44 + 45*rho^2))/(96*x20^4*(-1 + rho^2)^2*sigma^2) + 
  ((x2 - x20)^6*(4384 - 9105*rho^2 + 4725*rho^4))/(11520*x20^6*(-1 + rho^2)^3*sigma^2) + 
  (7*(x1 - x10)*(x2 - x20)^4*rho*(-13 + 15*rho^2))/(192*x20^(9/2)*(-1 + rho^2)^2*sigma) + 
  ((x1 - x10)*(x2 - x20)^3*(14*rho - 15*rho^3))/(24*x20^(7/2)*(-1 + rho^2)^2*sigma) +
  ((x1 - x10)*(x2 - x20)^5*(-764*rho + 1705*rho^3 - 945*rho^5))/(1920*x20^(11/2)*(-1 + rho^2)^3*sigma) + 
  ((x1 - x10)^3*(x2 - x20)^3*rho*(153 - 173*rho^2)*sigma)/(2880*x20^(9/2)*(-1 + rho^2)^3) + (5*(x1 - x10)^3*(x2 - x20)^2*rho*sigma)/(96*x20^(7/2)*
                                           (-1 + rho^2)^2) - ((x1 - x10)^3*(x2 - x20)*rho*sigma)/(24*x20^(5/2)*(-1 + rho^2)^2) - ((x1 - x10)^4*(x2 - x20)*
                                  sigma^2)/(96*x20^3*(-1 + rho^2)^2) + ((x1 - x10)^4*sigma^2)/
  (96*x20^2*(-1 + rho^2)^2) + ((x1 - x10)^4*(x2 - x20)^2*
                               (-9 + 14*rho^2)*sigma^2)/(960*x20^4*(-1 + rho^2)^3) - 
  ((x1 - x10)^5*(x2 - x20)*rho*sigma^3)/(480*x20^(7/2)*
                                         (-1 + rho^2)^3) + ((x1 - x10)^6*sigma^4)/(2880*x20^3*(-1 + rho^2)^3) + (3*(x1 - x10)*(x2 - x20)^2*
                                 rho)/(x20^(5/2)*(-4*sigma + 4*rho^2*sigma)) + (1/(2*x20^2*(-1 + rho^2)*sigma^2))*(x2^2 - 2*x2*sqrt(x20)*(sqrt(x20) + 
  (x1 - x10)*rho*sigma) + x20*(x20 + 
  2*(x1 - x10)*sqrt(x20)*rho*sigma + (x1 - x10)^2*sigma^2)) 
  
  
  c0 <- ((x1 - x10)^2*sigma^2)/
  (-48*x20 + 48*x20*rho^2) + (7*(x1 - x10)^4*sigma^4)/
  (11520*x20^2*(-1 + rho^2)^2) - 
  ((x2 - x20)^2*(24*x20*kappa - 48*theta*kappa + 
                 36*k1*sqrt(x20)*rho*sigma + 
                 12*k2*x20^(3/2)*rho*sigma + 5*x20*sigma^2))/(48*
                                                              x20^3*(-1 + rho^2)*
                                                              sigma^2) + ((x2 - x20)*(-4*theta*kappa + 
                                                                                      4*k1*sqrt(x20)*rho*sigma + 
                                                                                      4*k2*x20^(3/2)*rho*sigma + x20*(4*kappa + sigma^2)))/
  (4*x20^2*(-1 + rho^2)*sigma^2) - 
  ((x1 - x10)^3*
  sigma*(-12*theta*kappa*rho + 8*k1*sqrt(x20)*sigma + 
         x20*rho*(4*kappa + sigma^2)))/(192*
                                        x20^(5/2)*(-1 + rho^2)^2) - 
  ((x1 - x10)*(-4*theta*kappa*rho + 4*k1*sqrt(x20)*sigma + 
               4*k2*x20^(3/2)*sigma + x20*rho*(4*kappa + sigma^2)))/
  (4*x20^(3/2)*(-1 + rho^2)*sigma) + 
  ((x1 - x10)*(x2 - x20)*(-36*theta*kappa*rho + 
                          24*k1*sqrt(x20)*sigma + 
                          x20*rho*(12*kappa + sigma^2)))/(48*
                                                          x20^(5/2)*(-1 + rho^2)*sigma) + 
  ((x1 - x10)^3*(x2 - x20)*sigma*(-900*theta*kappa*rho + 
                                  480*k1*sqrt(x20)*sigma + 
                                  x20*rho*(180*kappa + 17*sigma^2)))/
  (11520*x20^(7/2)*(-1 + rho^2)^2) + 
  ((x1 - 
    x10)*(x2 - x20)^2*(12*theta*kappa*rho*(-10 + 13*rho^2) + 
                       8*k1*sqrt(x20)*(8 - 11*rho^2)*sigma - 
                       3*x20*(4*kappa*rho*(-2 + 3*rho^2) + rho^3*sigma^2)))/
  (192*x20^(7/2)*(-1 + rho^2)^2*sigma) + 
  ((x2 - x20)^4*(1440*theta*kappa*(-8 + 9*rho^2) - 
                 900*k2*x20^(3/2)*rho*
                 (-1 + rho^2)*sigma - 
                 420*k1*sqrt(x20)*rho*(-13 + 15*rho^2)*sigma + 
                 x20*(-360*kappa*(-8 + 9*rho^2) + (502 - 585*rho^2)*
                      sigma^2)))/
  (11520*x20^5*(-1 + rho^2)^2*sigma^2) + 
  (1/(192*x20^3*(-1 + rho^2)^2))*((x1 - x10)^2*(x2 - x20)*
                                  (-36*theta*kappa*rho^2 + 24*k1*sqrt(x20)*rho*sigma + 
                                  x20*(12*kappa*rho^2 + (2 + rho^2)*sigma^2))) + 
  ((x2 - x20)^3*(12*theta*kappa*(16 - 17*rho^2) + 
                 24*k2*x20^(3/2)*rho*
                 (-1 + rho^2)*sigma + 
                 8*k1*sqrt(x20)*rho*(-14 + 15*rho^2)*sigma + 
                 x20*(kappa*(-64 + 68*rho^2) + (-12 + 13*rho^2)*sigma^2)))/
  (192*x20^4*(-1 + rho^2)^2*sigma^2) - (1/(5760*x20^4*(-1 + rho^2)^2))*(x1 - x10)^2*(x2 - x20)^2*(-1620*theta*kappa*rho^2 + 
                             900*k1*sqrt(x20)*rho*sigma + x20*(360*kappa*rho^2 + 
                                                               (38 + 31*rho^2)*sigma^2)) + ((x1 - 
                                                                                             x10)*(x2 - x20)^3*
                                                                                            (-1260*theta*kappa*rho*(-5 + 8*rho^2) + 720*k1*sqrt(x20)*
                                                                                            (-4 + 7*rho^2)*sigma + 
                                                                                            x20*rho*(900*kappa*(-1 + 2*rho^2) + 
                                                                                                     (47 + 150*rho^2)*sigma^2)))/(11520*
                                                                                                                                  x20^(9/2)*(-1 + rho^2)^2*
                                                                                                                                  sigma) 
  
  
  
  
  
  
  c1 <- -(((x1 - x10)*sigma*(-3*theta*kappa*rho + k1*sqrt(x20)*sigma))/(24*x20^(3/2)*(-1 + rho^2))) + (1/(48*x20^3*(-1 + rho^2)*
                                            sigma^2))*
  (x2 - x20)*(-24*theta^2*kappa^2 + 
              36*k1*sqrt(x20)*theta*kappa*rho*sigma + 
              12*k2^2*x20^3*sigma^2 + 
              3*k2*x20^(5/2)*rho*sigma*(4*kappa + sigma^2) - 
              x20^(3/2)*rho*
              sigma*(-12*k2*theta*kappa + k1*(12*kappa + sigma^2)) + 
              12*x20*((-k1^2)*sigma^2 + 
                      2*theta*kappa*(kappa + sigma^2 - rho^2*sigma^2))) + 
  (1/(96*x20^2*(-1 + rho^2)*sigma^2))*(48*theta^2*kappa^2 - 
                                       96*k1*sqrt(x20)*theta*kappa*rho*sigma + 
                                       48*k2^2*x20^3*sigma^2 + 
                                       24*k2*x20^(5/2)*rho*sigma*(4*kappa + sigma^2) + 
                                       x20^2*(48*kappa^2 + 96*k1*k2*sigma^2 - 
                                              24*kappa*(-2 + rho^2)*sigma^2 + 
                                              (13 - 10*rho^2)*sigma^4) + 24*x20^(3/2)*rho*sigma*
                                       (-4*k2*theta*kappa + k1*(4*kappa + sigma^2)) - 
                                       24*
                                       x20*(-2*k1^2*sigma^2 + 
                                            theta*kappa*(4*kappa + (4 - 3*rho^2)*sigma^2))) + 
  (1/(23040*x20^3*(-1 + rho^2)^2))*((x1 - x10)^2*
                                    (240*theta^2*kappa^2*(4 + 9*rho^2) - 
                                    4320*k1*sqrt(x20)*theta*kappa*rho*
                                    sigma - 480*k2^2*x20^3*sigma^2 - 
                                    120*k2*x20^(5/2)*rho*sigma*
                                    (4*kappa + sigma^2) + x20^2*(240*kappa^2*rho^2 + 
                                                                 120*kappa*rho^2*sigma^2 + (4 + 11*rho^2)*sigma^4) + 
                                    120*x20^(3/2)*rho*
                                    sigma*(-4*k2*theta*kappa + 3*k1*(4*kappa + sigma^2)) - 
                                    120*x20*(-12*k1^2*sigma^2 + 
                                             theta*kappa*(4*kappa*(2 + 3*rho^2) + 
                                                          (8 - 3*rho^2)*sigma^2)))) + 
  (1/(11520*x20^(7/2)*(-1 + rho^2)^2*sigma))*((x1 - x10)*(x2 - 
                                                          x20)*
                                              (-240*theta^2*kappa^2*rho*(4 + 9*rho^2) + 
                                              4320*k1*sqrt(x20)*theta*kappa*
                                              rho^2*sigma + 480*k2^2*x20^3*rho*sigma^2 + 
                                              120*k2*x20^(5/2)*rho^2*
                                              sigma*(4*kappa + sigma^2) - 
                                              x20^2*(240*kappa^2*rho^3 + 
                                                     120*kappa*rho^3*sigma^2 + 
                                                     rho*(4 + 11*rho^2)*sigma^4) + 
                                              120*x20*
                                              rho*(-12*k1^2*sigma^2 + 
                                                   theta*kappa*(4*kappa*(2 + 3*rho^2) + 
                                                                (17 - 12*rho^2)*sigma^2)) - 
                                              120*x20^(3/2)*sigma*
                                              (-4*k2*theta*kappa*rho^2 + 
                                              k1*(12*kappa*rho^2 + (2 + rho^2)*
                                                  sigma^2)))) + (1/(23040*
                                                                    x20^4*(-1 + rho^2)^2*sigma^2))*
  ((x2 - x20)^2*(240*theta^2*kappa^2*(-56 + 69*rho^2) - 
                 1440*k1*sqrt(x20)*theta*kappa*rho*(-10 + 13*rho^2)*sigma - 
                 480*k2^2*x20^3*(-2 + 3*rho^2)*sigma^2 - 
                 120*k2*x20^(5/2)*rho*
                 (-4 + 5*rho^2)*sigma*(4*kappa + sigma^2) + 
                 x20^2*(240*kappa^2*rho^2 + 120*kappa*rho^2*sigma^2 + 
                        (4 + 11*rho^2)*sigma^4) + 
                 120*x20^(3/2)*rho*sigma*
                 (4*k2*theta*kappa*(8 - 9*rho^2) + 
                 3*k1*(-8*kappa + 12*kappa*rho^2 + 
                       rho^2*sigma^2)) - 
                 120*x20*(4*k1^2*(8 - 11*rho^2)*sigma^2 + 
                          5*theta*
                          kappa*(4*
                                 kappa*(-4 + 5*rho^2) + (-16 + 35*rho^2 - 18*rho^4)*
                                 sigma^2)))) 
  
  
  c2 <- (240*theta^2*kappa^2*(20 - 9*rho^2) - 
        2880*k1*sqrt(x20)*theta*kappa*rho*sigma - 
        960*k2*x20^(3/2)*theta*kappa*rho*sigma + 480*k2^2*x20^3*sigma^2 - 
        x20^2*(240*kappa^2*rho^2 + 
        120*kappa*rho^2*sigma^2 + (4 + 11*rho^2)*sigma^4) +
        120*x20*(4*k1^2*sigma^2 + 3*theta*kappa*(4*kappa*(-2 + rho^2) + 
        (-8 + 7*rho^2)*sigma^2)))/(11520*x20^2*(-1 + rho^2))
  output <- -log(2*pi*del) - Dv + cm1/del + c0 + c1*del + c2*(del^2/2)
  return(output)
}
