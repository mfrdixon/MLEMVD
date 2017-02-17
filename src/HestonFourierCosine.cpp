/*
* Sample code implementing the utility functions 
* for Heston Calibration. The Heston model is approximated with the Heston Fourier
* Cosine method. The error function relies on compiler vectorization,
* and OpenMP parallelisation.
*
* Authors: Matthew Dixon
* 
* Date: Nov. 30 2016
*/


#include <Rcpp.h>
#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>

#define pi 3.1415926535897932384626433832795

using namespace Rcpp;
// The Heston Characteristic Function
inline std::complex<double> HestonCF(
        double u,
        double T,
        double r,
        double q,
        double sigma,
        double lmbda,
        double meanV,
        double v0,
        double rho
        ) 
{
    std::complex<double> j1(0, 1);
    double a                    = lmbda *meanV;
    double b                    = lmbda;
    double sigma2               = sigma*sigma;	
    //double e      		= pow(rho*sigma*u-b,2);
    //double d      		= sqrt(e+(u*u+u)*sigma2);		
    std::complex<double> d      = sqrt((j1*rho*sigma*u-b)*(j1*rho*sigma*u-b) + (u*u+j1*u)*sigma2);
    std::complex<double> g      = (b-j1*rho*sigma*u - d)/(b-j1*rho*sigma*u+d);
    std::complex<double> ret    = exp(j1*u*(r-q)*T);
    ret                         *= exp((a/sigma2)*((b - rho*j1*sigma*u - d)*T - 2.0*log((1.0-g*exp(-d*T))/(1.0-g))));
    return ret*exp((v0/sigma2)*(b - rho*j1*sigma*u - d)*(1.0-exp(-d*T))/(1.0-g*exp(-d*T)));
}

// The Derivative of the Heston Characteristic Function w.r.t. to the v0 
inline std::complex<double> HestonCFdu0(
        double u,
        double T,
        double r,
        double q,
        double sigma,
        double lmbda,
        double meanV,
        double v0,
        double rho
        ) 
{
    std::complex<double> j1(0, 1);
    double a                    = lmbda *meanV;
    double b                    = lmbda;
    double sigma2               = sigma*sigma;	
    //double e      		= pow(rho*sigma*u-b,2);
    //double d      		= sqrt(e+(u*u+u)*sigma2);		
    std::complex<double> d      = sqrt((j1*rho*sigma*u-b)*(j1*rho*sigma*u-b) + (u*u+j1*u)*sigma2);
    std::complex<double> g      = (b-j1*rho*sigma*u - d)/(b-j1*rho*sigma*u+d);
    std::complex<double> ret    = exp(j1*u*(r-q)*T);
    ret                         *= exp((a/sigma2)*((b - rho*j1*sigma*u - d)*T - 2.0*log((1.0-g*exp(-d*T))/(1.0-g))));
    ret                         *= exp((v0/sigma2)*(b - rho*j1*sigma*u - d)*(1.0-exp(-d*T))/(1.0-g*exp(-d*T)));
    ret                         *=((1.0-exp(-d*T))/(1.0-g*exp(-d*T)))*(b-j1*rho*sigma*u-d)/(sigma2);
    return(ret);
}
// Functions required by the Fourier-Cosine method
inline double xi(
        double k,
        double a,
	      double b,
	      double c,
        double d) 
{
    double ret = 1.0/(1.0+pow(k*pi/(b-a),2))*(cos(k*pi*(d-a)/(b-a))*exp(d)-cos(k*pi*(c-a)/(b-a))*exp(c)+k*pi/(b-a)*sin(k*pi*(d-a)/(b-a))*exp(d)-k*pi/(b-a)*sin(k*pi*(c-a)/(b-a))*exp(c));
    return ret;
}


// Functions required by the Fourier-Cosine method
inline double psi(
        double k,
        double a,
        double b,
        double c,
        double d) 
{
    double ret=0.0;
    if (k==0)
       ret = d-c;
    else
       ret = (sin(k*pi*(d-a)/(b-a))-sin(k*pi*(c-a)/(b-a)))*(b-a)/(k*pi);
    return ret;
}
// The entry point for the Heston Fourier Cosine method
// [[Rcpp::export]]
double HestonCOS(
        double S,
        double K,
        double T,
        double r,
        double q,
        double sigma,
        double lmbda,
        double meanV,
        double v0,
        double rho,
        char otype,
        int N) 
{
    std::complex<double> j1(0, 1);
    double sigma2 = sigma*sigma;
    double lmbda2 = lmbda * lmbda;
    double c1 = r*T+(1-exp(-lmbda*T))*(meanV-v0)/(2.0*lmbda)-0.5*meanV*T;
    double c2 = (1.0/(8.0*lmbda2*lmbda))*(sigma*T*lmbda*exp(-lmbda*T)*(v0-meanV)*(8.0*lmbda*rho-4.0*sigma)+lmbda*rho*sigma*(1-exp(-lmbda*T))*(16.0*meanV-8.0*v0)+2.0*meanV*lmbda*T*(-4.0*lmbda*rho*sigma+sigma2+4.0*lmbda2)+sigma2*((meanV-2.0*v0)*exp(-2.0*lmbda*T)+meanV*(6.0*exp(-lmbda*T)-7.0)+2.0*v0)+8.0*lmbda2*(v0-meanV)*(1-exp(-lmbda*T)));
//c2 <- 1.0/(8.0*kappa**3)*(sigma*T*kappa*exp(-kappa*T)*(v0-theta)*(8.0*kappa*rho-4.0*sigma)+kappa*rho*sigma*(1-exp(-kappa*T))*(16.0*theta-8.0*v0)+2.0*theta*kappa*T*(-4.0*kappa*rho*sigma+sigma**2+4.0*kappa**2)+sigma**2*((theta-2.0*v0)*exp(-2.0*kappa*T)+theta*(6.0*exp(-kappa*T)-7.0)+2.0*v0)+8.0*kappa**2*(v0-theta)*(1-exp(-kappa*T)))
 
    double a = c1-12.0*sqrt(fabs(c2));
    double b = c1+12.0*sqrt(fabs(c2));	
    double x = log(S/K);
    double U, unit = 0.5;
    std::complex<double> ret (0.0,0.0);
 // Note that this for loop is independent of the strike
    for (int k=0; k < N; k++)
    {     
      U = 2.0/(b-a)*(xi(k,a,b,0,b) - psi(k,a,b,0,b));             
      std::complex<double> HCF = HestonCF(k*pi/(b-a),T,r,q,sigma,lmbda,meanV,v0,rho);
      ret += unit*HCF*exp(j1*double(k)*pi*(x-a)/(b-a))*U;
      unit = 1.0;
    }
    return K*exp(-r*T)*ret.real();
}
// The entry point for the Heston Fourier Cosine vega estimate 
// [[Rcpp::export]]
double HestonVega(
        double S,
        double K,
        double T,
        double r,
        double q,
        double sigma,
        double lmbda,
        double meanV,
        double v0,
        double rho,
        char otype,
        int N) 
{
    std::complex<double> j1(0, 1);
    double sigma2 = sigma*sigma;
    double lmbda2 = lmbda * lmbda;
    double c1 = r*T+(1-exp(-lmbda*T))*(meanV-v0)/(2.0*lmbda)-0.5*meanV*T;
    double c2 = 1.0/(8.0*lmbda2*lmbda)*(sigma*T*lmbda*exp(-lmbda*T)*(v0-meanV)*(8.0*lmbda*rho-4.0*sigma)+lmbda*rho*sigma*(1-exp(-lmbda*T))*(16.0*meanV-8.0*v0)+2.0*meanV*lmbda*T*(-4.0*lmbda*rho*sigma+sigma2+4.0*lmbda2)+sigma2*((meanV-2.0*v0)*exp(-2.0*lmbda*T)+meanV*(6.0*exp(-lmbda*T)-7.0)+2.0*v0)+8.0*lmbda2*(v0-meanV)*(1-exp(-lmbda*T)));
    double a = c1-12.0*sqrt(fabs(c2));
    double b = c1+12.0*sqrt(fabs(c2));	
    double x = log(S/K);
    double U, unit = 0.5;
    std::complex<double> ret (0.0,0.0);
 // Note that this for loop is independent of the strike
    for (int k=0; k < N; k++)
    {     
      U = 2.0/(b-a)*(xi(k,a,b,0,b) - psi(k,a,b,0,b));             
      std::complex<double> HCF = HestonCFdu0(k*pi/(b-a),T,r,q,sigma,lmbda,meanV,v0,rho);
      ret += unit*HCF*exp(j1*double(k)*pi*(x-a)/(b-a))*U;
      unit = 1.0;
    }
    return K*exp(-r*T)*ret.real();
}

// [[Rcpp::export]]
double getImpliedVolatility(double S,
        double V,
        double K,
        double T,
        double r,
        double q,
        double sigma,
        double lmbda,
        double meanV,
        double v0,
        double rho,
        char otype,
        double x_abs_err,
        int max_iter,
        int N)
{
  double v = v0;
  double v_n = v0;
  double error = x_abs_err;
  int iter  = 0;

  while(fabs(error)>=x_abs_err && iter < max_iter){
    double f = HestonCOS(S,K,T,r,q,sigma,lmbda,meanV,v,rho,otype,N)-V; 
    double g = HestonVega(S,K,T,r,q,sigma,lmbda,meanV,v,rho,otype,N); 
    v_n -= f/g; 
    error = v_n - v;
    v = v_n;
    iter++;
  }
  return(v);
}

int main(){
  double rate = 0.1;
  double q    = 0.0;
  double T_0  = 1.0;
  double sigma_0 = 0.00001; 
  double kappa_0 = 0.00001; 
  double v_0     = 0.1; 
  double theta_0 = v_0; 
  double rho_0   = 0.0;
  double HC = HestonCOS(45.22024623130299,45.22024623130299, T_0,rate,q,sigma_0,kappa_0,theta_0,v_0,rho_0,'C',128);
  std::cout<<HC<<std::endl;
  
  return(0);
}
