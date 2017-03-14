#' Compute the maximum likelihood estimate of the Heston model
#'
#' Heston model:
#' dln(S_t) = (\mu - \sigma^2/2) dt + \sqrt{V_t}dW_t^{1}
#' dV_t = \kappa(\theta - V_t)dt  + \sigma \sqrt{V_t}dW_t^{2}
#' 4 parameters to be estimated: (rho, kappa, theta, sigma)
#'
#' @param x: Observation of the state variable at time t
#' @param x0: Observation of the state variable at time t-1
#' @param del: The time step between the current and previous observation
#' @param param: The parameter 4-vector (rho,kappa>0,theta>0,sigma>0)
#' @export output: a list with a llk variable storing the result of the log likelihood calculation
#' @examples
#' ModelHeston(0.4,0.3,0.1,c(0.1,0.3,0.2,0.1))
#'

ModelHeston<- function(x,x0,del,param,args){

  output <- list()
  rho   <- param[1]
  kappa <- param[2]
  theta <- param[3]
  sigma <- param[4]


  if (args$mode=='implied')
  {
    S<-exp(x[1])
    K<-S
    objfun<-function(v){
      objfun<- HestonCOS(S,K,T_0,rate,q,param[4],param[2],param[3],v,param[1],args$callput, args$N)
    }

    v_0 <- getImpliedVolatility(S,x[3],K,T_0,rate,q,param[4],param[2],param[3],v_0,param[1],args$callput, 0.01,100, args$N) # the implied volatility
    #print(paste('implied vol estimate: ',v_0))
    # calculate the vega to obtain the Jacobian
    output$v <- v_0
    dVdv0 <- HestonVega(S,K,T_0,rate,q,sigma,kappa,theta,v_0,rho,args$callput,args$N)

    #print(paste('Jacobian: ',dVdv0))
    #print(paste('Numerical Jacobian', grad(objfun,x[2])))
    J <- dVdv0
    if (is.nan(log(J))){
      print(paste('Error: NAN occured in ModelHeston.R: ',J,sigma, kappa, theta,rho, v_0,T_0))
      output$llk <-0
    }
    else
      output$llk <- -log(J)
  }

  # m,a,b,s,r
  param_prime <- c(rate, kappa*theta, kappa, sigma, rho)
  model <- ModelB6(x,x0,del,param_prime)$llk
  if(is.nan(model)){
    print(param_prime)
    model <- 0
  }
  output$llk <- output$llk + model

  return(output)
}

