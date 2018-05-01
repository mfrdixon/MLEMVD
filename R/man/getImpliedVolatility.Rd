\name{getImpliedVolatility}
\alias{getImpliedVolatility}
\docType{package}
\title{
  Approximate the implied volatility of a derivative using the Heston model with the cosine method
}
\description{
  Simple function using Rcpp
}
\usage{
  getImpliedVolatility(S, V, K, T, r, q, sigma, lmbda, meanV, v0, rho, otype, x_abs_err, max_iter, N)
}

\arguments{
  \item{S}{Underlying price}

  \item{S}{Option price}

  \item{K}{Strike price of the option}

  \item{T}{The time to maturity of the option}

  \item{r}{The risk-free rate}

  \item{q}{The annualized divident yield}

  \item{sigma}{The volatility of the underlying}

  \item{lmbda}{The mean reversion rate of the volatility}

  \item{meanV}{Equilibrium volatility level}

  \item{v0}{The volatility of the volalility}

  \item{rho}{The correlation between the Brownian motions}

  \item{otype}{The option type: (P)ut or (C)all}

  \item{x_abs_err}{The target tolerance on the absolute error}

  \item{max_iter}{Max iterations for the Newton-Rhapson method}

  \item{N}{The number of Fourier-Cosine expansion terms}

}

\examples{
  \dontrun{
    getImpliedVolatility(S, V, K, T, r, q, sigma, lmbda, meanV, v0, rho, otype, x_abs_err, max_iter, N)
  }
}
