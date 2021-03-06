% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModelB3.R
\name{ModelB3}
\alias{ModelB3}
\title{Compute the maximum likelihood estimate of Model B3}
\usage{
ModelB3(x, x0, del, param)
}
\arguments{
\item{x}{Observation of the state variable at time t}

\item{x0}{Observation of the state variable at time t-1}

\item{del}{The time step between the current and previous observation}

\item{param}{The parameter 6-vector (mu,alpha,beta,sigma,rho,gamma)}
}
\description{
Compute the maximum likelihood estimate of Model B3
}
\examples{
ModelB3(0.4,0.3,0.1,c(0.1,0.2,0.3,0.4,0.5,0.6))

}
