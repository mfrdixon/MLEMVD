#' This file provides the R implementation of each pricing function

#' Compute the price of an option using the Heston model 
#' 
#' This Heston pricing function is provided for transparency and comparision with the CUDA implementation. This pricing function is not called by the error_function because the performance in R is too slow. The Fourier-Cosine method is used to compute the price of a European Option
#' @param inputParameter1 The underlying value of the asset \code{inputParameter1}
#' @param inputParameter2 The strike price of the option \code{inputParameter2}
#' @param inputParameter3 The maturity of the option \code{inputParameter3}
#' @param inputParameter4 The risk free annual short-rate  \code{inputParameter4}
#' @param inputParameter5 Heston parameter: sigma \code{inputParameter5}
#' @param inputParameter6 Heston parameter: kappa \code{inputParameter6}
#' @param inputParameter7 Heston Parameter: theta \code{inputParameter7}
#' @param inputParameter8 Heston Parameter: v0 \code{inputParameter8}
#' @param inputParameter9 Heston Parameters: rho \code{inputParameter9}
#' @param inputParameter10 Contract type a 'C'all or a 'P'ut \code{inputParameter10}
#' @param inputParameter11 Number of terms in Fourier cosine series \code{inputParameter11}
#'
#' @return output the option price 
#'
#' @keywords keywords
#'
#' @export HestonCOS
#' 
#' @examples
#' V<-HestonCOS(10.0,9.5,0.5,0.01,0.0,0.5,2.0,1.0,0.05,-0.8,'C') 

HestonCOS<-function(S,K,T,r,q,sigma,kappa,theta,v0,rho,otype, N=128){
  j  <- as.complex(1i)
  c1 <- r*T+(1-exp(-kappa*T))*(theta-v0)/(2.0*kappa)-0.5*theta*T
  c2 <- 1.0/(8.0*kappa**3)*(sigma*T*kappa*exp(-kappa*T)*(v0-theta)*(8.0*kappa*rho-4.0*sigma)+kappa*rho*sigma*(1-exp(-kappa*T))*(16.0*theta-8.0*v0)+2.0*theta*kappa*T*(-4.0*kappa*rho*sigma+sigma**2+4.0*kappa**2)+sigma**2*((theta-2.0*v0)*exp(-2.0*kappa*T)+theta*(6.0*exp(-kappa*T)-7.0)+2.0*v0)+8.0*kappa**2*(v0-theta)*(1-exp(-kappa*T)))
  a <- c1-12.0*sqrt(abs(c2))
  b <- c1+12.0*sqrt(abs(c2))
  x <- log(S/K)
  k <- seq(0,N-1)
  
  if (otype == "C")
    U <- 2.0/(b-a)*(xi(k,a,b,0,b) - psi(k,a,b,0,b))
  else
    U <- 2.0/(b-a)*(xi(k,a,b,0,a) - psi(k,a,b,0,a))
  
  unit <- rep(1,N)
  unit[1] <- 0.5
  ret <- as.complex(0)
  # Note that HestonCF is independent of the strike
  HCF <- HestonCF(k*pi/(b-a),T,r,q,sigma,kappa,theta,v0,rho)
  
  for (i in 1:N)
    ret <- ret + unit[i]*HCF[i]*exp(j*k[i]*pi*(x-a)/(b-a))*U[i]
  
  return (Re(K*exp(-r*T)*ret))
}

#' Compute the characteristic function of the Heston model 
#' 
#' This characteristic function is provided for transparency and comparision with the CUDA implementation. This function is not called by the error_function because the performance in R is too slow. 
#'
#' @return 
#'
#' @keywords HestonCF 
#'
#'  
#' @examples
#' 
HestonCF<-function(u,T,r,q,sigma,kappa,theta,v0,rho){
  j  <- as.complex(1i)
  a <- kappa*theta
  b <- kappa
  d <- sqrt((j*rho*sigma*u-b)**2+(u**2+j*u)*sigma**2)
  g <- (b-j*rho*sigma*u-d)/(b-j*rho*sigma*u+d)
  ret <- exp(j*u*(r-q)*T)
  ret <- ret*exp((a/sigma**2)*((b - rho*j*sigma*u - d)*T - 2.0*log((1-g*exp(-d*T))/(1-g))))
  return (ret*exp((v0/sigma**2)*(b - rho*j*sigma*u - d)*(1-exp(-d*T))/(1-g*exp(-d*T))))
}
#' Fourier Cosine functions
xi<-function(k,a,b,c,d){
  ret <- 1.0/(1+(k*pi/(b-a))**2)*(cos(k*pi*(d-a)/(b-a))*exp(d)-cos(k*pi*(c-a)/(b-a))*exp(c)+k*pi/(b-a)*sin(k*pi*(d-a)/(b-a))*exp(d)-k*pi/(b-a)*sin(k*pi*(c-a)/(b-a))*exp(c))
  return (ret)
}
#' Psi function
psi<-function(k,a,b,c,d){
  N <- length(k)
  idx <- seq(2, N)
  ret <- rep(0,N)
  ret[1] <- d-c
  ret[idx] <-(sin(k[idx]*pi*(d-a)/(b-a))-sin(k[idx]*pi*(c-a)/(b-a)))*(b-a)/(k[idx]*pi)
  return(ret)
}