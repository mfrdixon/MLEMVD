\name{HestonVega}
\alias{HestonVega}
\docType{package}
\title{
  Approximate the vega of a derivative using the Heston model with the cosine method
}
\description{
  Simple function using Rcpp
}
\usage{
  HestonVega(S, K, T, r, q, sigma, lmbda, meanV, v0, rho, otype, N)
}

\arguments{
  \item{S}{Underlying price}

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

  \item{N}{The number of Fourier-Cosine expansion terms}

}

\examples{
  \dontrun{
    HestonVega(S, K, T, r, q, sigma, lmbda, meanV, v0, rho, otype, N)
  }
}
