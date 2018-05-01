#' 20507 Quotes of the ZNGA Option Chain.
#'
#' A dataset containing the quotes and other attributes of ZNGA options over multiple
#' 30 second snapshots.
#'
#' @format A data frame with 20507 rows and 8 variables:
#' \describe{
#'   \item{symbol}{underlying symbol name.}
#'   \item{timestamp}{POSIX time stamp of quote.}
#'   \item{type}{Option type: (P)ut or (C)all.}
#'   \item{maturity}{POSIX time stamp of option maturity.}
#'   \item{strike}{Strike price.}
#'   \item{underlying}{underlying price.}
#'   \item{bid}{bid price.}
#'   \item{ask}{ask price.}
#'   ...
#' }
#' @source undisclosed.
"ZNGA"
